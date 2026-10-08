"""
    GmshLibm

Transcendental functions resolved against the platform C math library — the
same `libm` OpenCASCADE and Gmsh binaries link against. Julia's `sin`, `cos`,
`atan`, `asin`, `acos`, `tan`, `atan2`, `exp`, `log`, `pow`, `sinh`, `cosh`,
`tanh` use openlibm-derived kernels that can differ from the system library by
one ulp on a substantial fraction of inputs (measured on this platform:
`atan2` ~34%, `tan` ~40%, `acos` ~18%, `sin`/`cos` ~4%). Every Gmsh- or
OCCT-emulating formula calls these shims so results stay bit-identical to a
`libgmsh` build on the same platform.

IEEE-exact operations (`sqrt`, `abs`, `fma`, `floor`, `ceil`, `trunc`, `rem`,
integer casts) are implementation-independent and are not shimmed. When no
system math library can be loaded the wrappers fall back to Julia's
implementations — still correctly rounded-grade math, only not guaranteed
bit-identical to a co-platform Gmsh binary.
"""
module GmshLibm

using Libdl

export _gm_sin, _gm_cos, _gm_tan, _gm_asin, _gm_acos, _gm_atan, _gm_atan2,
       _gm_sinh, _gm_cosh, _gm_tanh, _gm_exp, _gm_log, _gm_log10, _gm_pow,
       _gm_sincos,
       _gm87_sin, _gm87_cos, _gm87_exp, _gm87_log, _gm87_pow, _gm87_atan2,
       _gm87_sincos

function _gm_mathlib_handle()
    candidates = Sys.isapple()    ? ("libSystem.B.dylib", "libm.dylib") :
                 Sys.iswindows() ? ("msvcrt.dll", "ucrtbase.dll") :
                 ("libm.so.6", "libm.so", "libc.so.6", "libSystem.B.dylib")
    for name in candidates
        handle = Libdl.dlopen(name; throw_error=false)
        handle === nothing || return handle
    end
    return nothing
end

const _GM_LIBM_FUNS = (:sin, :cos, :tan, :asin, :acos, :atan, :atan2,
                       :sinh, :cosh, :tanh, :exp, :log, :log10, :pow)

# Entry-point addresses are process state, not constants: a precompile image
# that serialized them would call whatever occupied the recorded address after
# the next boot's ASLR shuffle (dlopen bases are per-boot on Windows). They
# live in `Ref` cells repopulated by `__init__` on every module load; the
# top-level population covers the precompiling process itself.
const _GM_LIBM = NamedTuple{_GM_LIBM_FUNS}(
    ntuple(_ -> Ref{Ptr{Cvoid}}(Ptr{Cvoid}(0)), length(_GM_LIBM_FUNS)))

# `sin`+`cos` of one operand fused by the C++ toolchain: clang emits
# `__sincos_stret` (macOS) and gcc emits `sincos` (glibc) when a translation
# unit calls both on the same argument. The fused sin can differ one ulp from
# standalone `sin`, so paired evaluations need the combined entry point.
const _GM_SINCOS_STRET = Ref{Ptr{Cvoid}}(Ptr{Cvoid}(0))
const _GM_SINCOS = Ref{Ptr{Cvoid}}(Ptr{Cvoid}(0))

function _gm_resolve_libm!()
    handle = _gm_mathlib_handle()
    for name in _GM_LIBM_FUNS
        _GM_LIBM[name][] = handle === nothing ? Ptr{Cvoid}(0) :
            something(Libdl.dlsym(handle, String(name);
                                throw_error=false), Ptr{Cvoid}(0))
    end
    _GM_SINCOS_STRET[] = (!Sys.isapple() || handle === nothing) ?
        Ptr{Cvoid}(0) :
        something(Libdl.dlsym(handle, "__sincos_stret"; throw_error=false),
                  Ptr{Cvoid}(0))
    _GM_SINCOS[] = (Sys.isapple() || Sys.iswindows() || handle === nothing) ?
        Ptr{Cvoid}(0) :
        something(Libdl.dlsym(handle, "sincos"; throw_error=false),
                  Ptr{Cvoid}(0))
    return nothing
end

_gm_resolve_libm!()
__init__() = (_gm_resolve_libm!(); _x87_init_consts!())

# Whether the platform C math library resolved (bit-parity available).
_gm_libm_available() = _GM_LIBM.sin[] != C_NULL

# ---------------------------------------------------------------------------
# Windows: the shipped `gmsh.exe` is pure MinGW — its only CRT import is
# `msvcrt.dll`. The import table shows `acos, asin, atan, cosh, log10, sinh,
# tan, tanh` resolved from `msvcrt.dll`, while `sin, cos, exp, log, pow, atan2`
# are statically linked: mingwex's x87-FPU implementations. Those compute in
# 80-bit extended registers — a 64-bit significand — and trigonometric
# argument reduction uses the FPU's 66-bit π constant (0x3.243F6A8885A308D3),
# not msvcrt's algorithms. `pow` uses binary squaring for integer exponents
# and the extended-precision `y·log2(x)` chain otherwise.
#
# The functions below emulate that hardware semantics: each elementary step
# is rounded to the 64-bit extended significand through `BigFloat` (MPFR
# precision 64), then stored back to `Float64`. Verified bit-for-bit against
# the 4.15.2 Windows binary on fuzzed differential sweeps.
# ---------------------------------------------------------------------------

const _X87_LIM = 9.223372036854776e18    # 2^63 — FPU trig reduction range limit
const _X87_PI_HI = 3.141592653589793     # Float64(π_hw), top 53 bits
const _X87_PI_LO = 1.2246063538223773e-16 # π_hw − _X87_PI_HI (bits 54–66)

const _X87_PI_HW = Ref{BigFloat}()       # exact 66-bit hardware π
const _X87_LN2 = Ref{BigFloat}()         # FLDLN2 constant: round64(ln 2)
const _X87_L2E = Ref{BigFloat}()         # FLDL2E constant: round64(log₂ e)
const _X87_2PI_MANT = Ref{BigInt}()      # `fldpi`+`fadd` modulus, odd mantissa
const _X87_2PI_EXPO = Ref{Int}()         # ... times 2^this

# v = m·2^e with m an odd BigInt — exact decomposition of a dyadic value.
function _big_mant_exp(v::BigFloat)
    v == 0 && return (BigInt(0), 0)
    p = precision(v)
    e = exponent(v)                      # v ∈ [2^(e-1), 2^e)
    m = BigInt(ldexp(v, p - e))
    sh = trailing_zeros(m)
    return (m >> sh, e - p + sh)
end
function _f64_mant_exp(v::Float64)       # |v| = m·2^e, m odd BigInt
    f, e = frexp(v)                      # v = f·2^e, f ∈ [0.5,1)
    m = BigInt(ldexp(abs(f), 53))
    sh = trailing_zeros(m)
    return (m >> sh, e - 53 + sh)
end

function _x87_init_consts!()
    _X87_PI_HW[] = BigFloat(_X87_PI_HI; precision=256) +
                   BigFloat(_X87_PI_LO; precision=256)
    ln2 = log(BigFloat(2; precision=256))
    _X87_LN2[] = BigFloat(ln2; precision=64)
    _X87_L2E[] = BigFloat(BigFloat(1; precision=256) / ln2; precision=64)
    # `fldpi` loads π rounded to the 64-bit extended significand and the
    # fallback's `fadd %st(0)` doubles it exactly → M = 2·round64(π_hw).
    _X87_2PI_MANT[], _X87_2PI_EXPO[] =
        _big_mant_exp(2 * BigFloat(_X87_PI_HW[]; precision=64))
    return nothing
end

_x87_init_consts!()

@inline _x87b(x::Float64) = BigFloat(x; precision=256)
@inline _x87e(x::BigFloat) = BigFloat(x; precision=64)  # round to extended

# IEEE remainder `x − q·M` with q = round-to-nearest-even(x/M), computed
# exactly in BigInt arithmetic (x and M are both dyadic rationals) and
# rounded once to the 64-bit extended significand — the `fprem1` loop in
# mingw's `sinl_internal.S`/`cosl_internal.S` fallback for |x| ≥ 2^63.
function _x87_fprem1_2pi(x::Float64)
    mx, ex = _f64_mant_exp(x)
    mm, em = _X87_2PI_MANT[], _X87_2PI_EXPO[]
    if ex >= em
        num = mx << (ex - em); den = mm
    else
        num = mx; den = mm << (em - ex)
    end
    q0, rem = divrem(abs(num), den)
    c = cmp(2 * rem, den)
    q = (c > 0 || (c == 0 && isodd(q0))) ? q0 + 1 : q0
    rint = num - q * den                 # residue = rint · 2^min(ex,em)
    rexp = min(ex, em)
    prec = max(128, ndigits(abs(rint), base=2) + 66)
    r = ldexp(BigFloat(rint; precision=prec), rexp)
    return signbit(x) ? -_x87e(r) : _x87e(r)
end

# `fsin`/`fcos` on a 64-bit-extended argument: reduce mod the 66-bit π_hw,
# kernel-evaluate the residue at extended precision, apply the parity flip.
function _x87_trig(r::BigFloat)
    k = round(BigInt, _x87e(r / _X87_PI_HW[]))
    rr = _x87e(r - k * _X87_PI_HW[])
    rb = BigFloat(rr; precision=256)
    s = Float64(_x87e(sin(rb)))
    c = Float64(_x87e(cos(rb)))
    return isodd(k) ? (-s, -c) : (s, c)
end

function _win_sincos(x::Float64)
    ax = abs(x)
    if ax < _X87_LIM
        ax <= _X87_PI_HI/2 && return (Float64(_x87e(sin(_x87b(x)))),
                                      Float64(_x87e(cos(_x87b(x)))))
        return _x87_trig(_x87b(x))
    end
    isnan(x) && return (x, x)
    isinf(x) && return ((x - x) / (x - x), (x - x) / (x - x))
    # |x| ≥ 2^63: the FPU range check fails and the assembly fallback
    # reduces by `fprem1` modulo `2·fldpi` before `fsin`/`fcos`.
    return _x87_trig(_x87_fprem1_2pi(x))
end

function _win_exp(x::Float64)
    (isnan(x) || x == Inf) && return x
    x == -Inf && return 0.0
    v = _x87e(_x87b(x) * _X87_L2E[])
    vf = Float64(v)
    abs(vf) >= _X87_LIM && return vf > 0 ? Inf : 0.0
    k = round(Int64, vf)
    f = _x87e(v - k)
    s = _x87e(exp2(BigFloat(f; precision=256)))
    return Float64(s * exp2(BigFloat(k; precision=256)))
end

function _win_log(x::Float64)
    isnan(x) && return x
    x < 0.0 && return NaN
    x == 0.0 && return -Inf
    x == Inf && return x
    l2 = _x87e(log2(_x87b(x)))
    return Float64(_x87e(l2 * _X87_LN2[]))
end

# Binary squaring matching the binary's integer-exponent path. Mantissas are
# plain `Float64` products — bit-identical to naive squaring — but excess
# exponent is tracked separately so `x^n` past the overflow boundary (e.g.
# `2^-1074`) still scales correctly instead of collapsing through `Inf`.
@inline function _powi_extract(v::Float64)
    isfinite(v) || return (v, 0)
    s = exponent(v)
    return ldexp(v, -s), s
end
function _win_powi(x::Float64, n::Int64)
    neg = n < 0
    e = neg ? -n : n
    # A finite Float64 base contributes at most 1023 binary exponent bits
    # per factor. Int64 integer powers can therefore require more than 64
    # exponent bits; Int128 keeps every squaring and sum exact without
    # allocating. Wrapping an Int accumulator reverses overflow/underflow.
    r, re = 1.0, Int128(0)
    b, be = x, Int128(0)
    # Normalize an already-large base before its first square. Otherwise
    # x*x can overflow before extraction, losing a representable reciprocal
    # power such as (2^512)^(-2) = 2^(-1024).
    if abs(b) >= 7.237005577332262e75
        b, s = _powi_extract(b); be += s
    end
    while e > 0
        if isodd(e)
            r *= b; re += be
            if abs(r) >= 7.237005577332262e75      # 2^256 — keep products finite
                r, s = _powi_extract(r); re += s
            end
        end
        b *= b; be += be
        if abs(b) >= 7.237005577332262e75
            b, s = _powi_extract(b); be += s
        end
        e >>= 1
    end
    m = neg ? 1.0 / r : r
    final_exponent = neg ? -re : re
    # ldexp takes a machine Int. Exponents outside that range are already
    # far beyond Float64's range; retain the sign of the actual mantissa.
    final_exponent > typemax(Int) && return copysign(Inf, m)
    final_exponent < typemin(Int) && return copysign(0.0, m)
    # Julia's ldexp can prematurely return zero for m<1 at exponent -1074.
    # Scale to units of the minimum subnormal first (still an exact normal
    # operation), then let one hardware multiplication round the final value.
    # The bounded mantissa makes the intermediate safe from overflow.
    final_exponent < -1022 &&
        return ldexp(m, Int(final_exponent)+1074) * 5.0e-324
    return ldexp(m, Int(final_exponent))
end

function _win_pow(x::Float64, y::Float64)
    y == 0.0 && return 1.0                          # NaN^0 = 1, like fdlibm
    (isnan(x) || isnan(y)) && return x + y
    if isinf(y)
        absx = abs(x)
        absx == 1.0 && return 1.0                   # pow(±1, ±Inf) = 1
        return (absx > 1) == (y > 0) ? Inf : 0.0
    end
    if isinf(x)
        # (-Inf)^y: an odd-integer exponent flips the sign; every other
        # exponent (non-integer included) behaves like +Inf^y. Doubles
        # ≥ 2^53 are all even, so `isodd` settles parity wherever it applies.
        if x < 0.0 && isinteger(y) && abs(y) < _X87_LIM &&
           isodd(Int64(y))
            return y > 0.0 ? -Inf : -0.0
        end
        return y > 0.0 ? Inf : 0.0
    end
    if isinteger(y) && abs(y) < _X87_LIM
        return _win_powi(x, Int64(y))
    end
    if x < 0.0
        # Integer exponents past the Int64 range are all even — the binary
        # squares |x|. Non-integer exponents on a negative base are the
        # x87 indefinite NaN.
        isinteger(y) || return (x - x) / (x - x)
        x = -x
    end
    x == 0.0 && return y < 0 ? Inf : 0.0
    l2 = _x87e(log2(_x87b(x)))
    v = _x87e(_x87b(y) * l2)
    vf = Float64(v)
    vf >= _X87_LIM && return Inf
    vf <= -_X87_LIM && return 0.0
    k = round(Int64, vf)
    f = _x87e(v - k)
    s = _x87e(BigFloat(2; precision=256) ^ BigFloat(f; precision=256))
    return Float64(s * exp2(BigFloat(k; precision=256)))
end

function _win_atan2(y::Float64, x::Float64)
    (isnan(y) || isnan(x)) && return y + x
    return Float64(_x87e(atan(_x87b(y), _x87b(x))))
end

@inline _gm_sin(x::Float64) = _GM_LIBM.sin[] == C_NULL ? sin(x) :
    ccall(_GM_LIBM.sin[], Float64, (Float64,), x)
@inline _gm_cos(x::Float64) = _GM_LIBM.cos[] == C_NULL ? cos(x) :
    ccall(_GM_LIBM.cos[], Float64, (Float64,), x)
@inline _gm_tan(x::Float64) = _GM_LIBM.tan[] == C_NULL ? tan(x) :
    ccall(_GM_LIBM.tan[], Float64, (Float64,), x)
@inline _gm_asin(x::Float64) = _GM_LIBM.asin[] == C_NULL ? asin(x) :
    ccall(_GM_LIBM.asin[], Float64, (Float64,), x)
@inline _gm_acos(x::Float64) = _GM_LIBM.acos[] == C_NULL ? acos(x) :
    ccall(_GM_LIBM.acos[], Float64, (Float64,), x)
@inline _gm_atan(x::Float64) = _GM_LIBM.atan[] == C_NULL ? atan(x) :
    ccall(_GM_LIBM.atan[], Float64, (Float64,), x)
@inline _gm_sinh(x::Float64) = _GM_LIBM.sinh[] == C_NULL ? sinh(x) :
    ccall(_GM_LIBM.sinh[], Float64, (Float64,), x)
@inline _gm_cosh(x::Float64) = _GM_LIBM.cosh[] == C_NULL ? cosh(x) :
    ccall(_GM_LIBM.cosh[], Float64, (Float64,), x)
@inline _gm_tanh(x::Float64) = _GM_LIBM.tanh[] == C_NULL ? tanh(x) :
    ccall(_GM_LIBM.tanh[], Float64, (Float64,), x)
@inline _gm_exp(x::Float64) = _GM_LIBM.exp[] == C_NULL ? exp(x) :
    ccall(_GM_LIBM.exp[], Float64, (Float64,), x)
@inline _gm_log(x::Float64) = _GM_LIBM.log[] == C_NULL ? log(x) :
    ccall(_GM_LIBM.log[], Float64, (Float64,), x)
@inline _gm_log10(x::Float64) = _GM_LIBM.log10[] == C_NULL ? log10(x) :
    ccall(_GM_LIBM.log10[], Float64, (Float64,), x)

@inline _gm_atan2(y::Float64, x::Float64) = _GM_LIBM.atan2[] == C_NULL ?
    atan(y, x) :
    ccall(_GM_LIBM.atan2[], Float64, (Float64, Float64), y, x)
@inline _gm_pow(x::Float64, y::Float64) = _GM_LIBM.pow[] == C_NULL ? x^y :
    ccall(_GM_LIBM.pow[], Float64, (Float64, Float64), x, y)

# Returns `(sin(x), cos(x))` through the toolchain's fused entry point —
# `__sincos_stret` on macOS (two-double register struct return), `sincos` on
# Linux (out-param form). Falls back to the separate shims, which in turn fall
# back to Julia's builtins when no system library resolves.
@inline function _gm_sincos(x::Float64)
    if _GM_SINCOS_STRET[] != C_NULL
        return ccall(_GM_SINCOS_STRET[], NTuple{2,Float64}, (Float64,), x)
    elseif _GM_SINCOS[] != C_NULL
        s = Ref{Float64}(); c = Ref{Float64}()
        ccall(_GM_SINCOS[], Cvoid, (Float64, Ptr{Float64}, Ptr{Float64}), x, s, c)
        return (s[], c[])
    end
    return (_gm_sin(x), _gm_cos(x))
end

# ---------------------------------------------------------------------------
# `_gm87_*` — the statically-linked mingwex x87 implementations. `.geo` FExpr
# evaluation mirrors code compiled inside `libgmsh` itself: on Windows that
# binary is pure MinGW, so expression-visible `Sin`/`Cos`/`Exp`/`Log`/`^`/
# `Atan2` values take the emulated `_win_*` semantics verified bit-for-bit
# against gmsh.exe (the evaluator already allocates for `BigFloat`). Off
# Windows upstream links the platform `libm` dynamically, so the family
# aliases `_gm_*`. Everywhere else stays on `_gm_*`: transform matrices and
# geometry-evaluation internals run in zero-allocation hot paths and their
# captured upstream outputs were established under the system-libm
# resolution.
# ---------------------------------------------------------------------------
@inline _gm87_sin(x::Float64) = Sys.iswindows() ? _win_sincos(x)[1] : _gm_sin(x)
@inline _gm87_cos(x::Float64) = Sys.iswindows() ? _win_sincos(x)[2] : _gm_cos(x)
@inline _gm87_exp(x::Float64) = Sys.iswindows() ? _win_exp(x) : _gm_exp(x)
@inline _gm87_log(x::Float64) = Sys.iswindows() ? _win_log(x) : _gm_log(x)
@inline _gm87_pow(x::Float64, y::Float64) =
    Sys.iswindows() ? _win_pow(x, y) : _gm_pow(x, y)
@inline _gm87_atan2(y::Float64, x::Float64) =
    Sys.iswindows() ? _win_atan2(y, x) : _gm_atan2(y, x)
@inline _gm87_sincos(x::Float64) =
    Sys.iswindows() ? _win_sincos(x) : _gm_sincos(x)

end # module GmshLibm
