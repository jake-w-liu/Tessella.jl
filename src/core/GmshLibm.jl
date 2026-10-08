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
# not msvcrt's algorithms. `pow` uses binary squaring for signed Int32
# integer exponents and the extended-register `y·log2(x)` chain otherwise.
#
# The existing trigonometric/exponential shims retain their MPFR emulation.
# Power executes the native x87 instructions on Windows x86_64: transcendental
# instructions retain extended values while arithmetic observes the caller's
# control word. Hosted Gmsh DLLs inherit that word; the standalone executable
# starts with 64-bit arithmetic precision. A task-local CLI context selects
# that precision only during the non-yielding native power kernel, so file I/O
# and task scheduling cannot leak hardware state into another Julia task.
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

# MinGW's pow.def.h admits the complete signed Int32 range, inclusively.
# Larger integral exponents use log2l/exp2l, whose rounding is observably
# different from repeated Float64 squaring near one.
const _WIN_X87_NATIVE = Sys.iswindows() && Sys.ARCH === :x86_64
const _WIN_CLI_PRECISION_KEY = gensym(:TessellaCLIPrecision)

@static if Sys.iswindows() && Sys.ARCH === :x86_64
    @inline function _win_x87_control_word()
        Base.llvmcall(raw"""
            %word = alloca i16, align 2
            call void asm sideeffect "fnstcw ($0)", "r,~{memory}"(ptr %word)
            %result = load i16, ptr %word, align 2
            ret i16 %result
            """,UInt16,Tuple{})
    end

    @inline function _win_x87_set_control_word!(word::UInt16)
        Base.llvmcall(raw"""
            %slot = alloca i16, align 2
            store i16 %0, ptr %slot, align 2
            call void asm sideeffect "fldcw ($0)", "r,~{memory},~{fpsr}"(ptr %slot)
            ret void
            """,Cvoid,Tuple{UInt16},word)
    end

    @inline function _win_x87_sqrt(x::Float64)
        Base.llvmcall(raw"""
            %input = alloca double, align 8
            %out = alloca double, align 8
            store double %0, ptr %input, align 8
            call void asm sideeffect "fldl ($1); fsqrt; fstpl ($0)", "r,r,~{st},~{memory},~{fpsr}"(ptr %out,ptr %input)
            %result = load double, ptr %out, align 8
            ret double %result
            """,Float64,Tuple{Float64},x)
    end

    # Independent native implementation of the log2l/exp2l instruction path.
    # Positive finite x and finite y are admitted by _win_pow. The logarithm
    # chooses FYL2XP1 for |x-1|<=0.29 in the caller's arithmetic precision;
    # otherwise it uses FYL2X. No extended intermediate is stored as Float64.
    # Only exponent truncation changes the control word, which is restored
    # before F2XM1/FSCALE. All x87 stack entries are popped before returning.
    @inline function _win_x87_powlog(x::Float64,y::Float64)
        Base.llvmcall(raw"""
            %xs = alloca double, align 8
            %ys = alloca double, align 8
            %out = alloca double, align 8
            %limit = alloca double, align 8
            %cw = alloca [2 x i16], align 2
            %cwnew = getelementptr [2 x i16], ptr %cw, i64 0, i64 1
            store double %0, ptr %xs, align 8
            store double %1, ptr %ys, align 8
            store double 0x3FD28F5C28F5C28F, ptr %limit, align 8
            call void asm sideeffect "fld1; fldl ($1); fld %st(0); fsub %st(2), %st(0); fld %st(0); fabs; fldl ($5); fcomip %st(1), %st(0); fstp %st(0); jae 0f; fstp %st(0); fyl2x; jmp 1f; 0: fstp %st(1); fyl2xp1; 1: fldl ($2); fmulp; fld %st(0); fnstcw ($3); movzwl ($3), %eax; orb $$12, %ah; movw %ax, ($4); fldcw ($4); frndint; fldcw ($3); fsubr %st(0), %st(1); fxch; f2xm1; fld1; faddp; fscale; fstp %st(1); fstpl ($0)", "r,r,r,r,r,r,~{rax},~{st},~{st(1)},~{st(2)},~{st(3)},~{st(4)},~{memory},~{dirflag},~{fpsr},~{flags}"(ptr %out,ptr %xs,ptr %ys,ptr %cw,ptr %cwnew,ptr %limit)
            %result = load double, ptr %out, align 8
            ret double %result
            """,Float64,Tuple{Float64,Float64},x,y)
    end
end

@inline function _with_win_cli_precision(callback::F) where F
    @static if Sys.iswindows() && Sys.ARCH === :x86_64
        # Base restores the previous task-local value in finally, including
        # nested calls and exceptions. Hardware state is untouched here: the
        # callback may yield or migrate to a different OS thread during I/O.
        return task_local_storage(callback,_WIN_CLI_PRECISION_KEY,true)
    else
        return callback()
    end
end

@inline function _win_cli_precision()
    # Reading an absent key must not initialize a task's storage on a hot API
    # path. The task field is nothing until a task-local value has been set.
    storage=current_task().storage
    return storage!==nothing && get(storage,_WIN_CLI_PRECISION_KEY,false)===true
end

function _win_powi_product(base::Float64,exponent::UInt32)
    result=isodd(exponent) ? base : 1.0
    exponent>>=1
    while true
        base*=base
        isodd(exponent) && (result*=base)
        exponent>>=1
        exponent==0 && return result
    end
end

function _win_powi(x::Float64,n::Int32)
    base=abs(x)
    # Widen before negation: typemin(Int32) belongs to the primary fast path.
    exponent=UInt32(n<0 ? -Int64(n) : n)
    result=exponent==0 ? 1.0 : exponent==1 ? base :
           _win_powi_product(base,exponent)
    if n<0
        # MinGW retries an overflowing reciprocal power with 1/base first.
        # This retains representable subnormals with the primary's rounding.
        result=isinf(result) && base>1.0 ?
            _win_powi_product(1.0/base,exponent) : 1.0/result
    end
    return signbit(x) && isodd(n) ? -result : result
end

@inline _win_pow_odd(y::Float64)=isinteger(y) &&
    abs(y)<9007199254740992.0 && isodd(Int64(y))

function _win_pow(x::Float64,y::Float64)
    _WIN_X87_NATIVE || return _gm_pow(x,y)
    if _win_cli_precision()
        saved=_win_x87_control_word()
        _win_x87_set_control_word!((saved&0xfcff)|UInt16(0x0300))
        try
            # This kernel contains no Julia allocation, callbacks, or yielding
            # operations. Restore this OS thread before task scheduling resumes.
            return _win_pow_native(x,y)
        finally
            _win_x87_set_control_word!(saved)
        end
    end
    return _win_pow_native(x,y)
end

function _win_pow_native(x::Float64,y::Float64)
    (y==0.0 || x==1.0) && return 1.0
    isnan(x) && return x
    isnan(y) && return y
    if x==0.0
        magnitude=y<0.0 ? Inf : 0.0
        return signbit(x) && _win_pow_odd(y) ? -magnitude : magnitude
    end
    if isinf(y)
        abs(x)==1.0 && return 1.0
        return (abs(x)>1.0)==(y>0.0) ? Inf : 0.0
    end
    if isinf(x)
        magnitude=y>0.0 ? Inf : 0.0
        return signbit(x) && _win_pow_odd(y) ? -magnitude : magnitude
    end
    isinteger(y) && typemin(Int32)<=y<=typemax(Int32) &&
        return _win_powi(x,Int32(y))
    x<0.0 && !isinteger(y) && return -NaN
    y==0.5 && return _win_x87_sqrt(x)
    magnitude=_win_x87_powlog(abs(x),y)
    return x<0.0 && _win_pow_odd(y) ? -magnitude : magnitude
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
