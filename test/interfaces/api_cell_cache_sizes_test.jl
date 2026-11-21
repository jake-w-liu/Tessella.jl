module CellCacheSizesTests
using Test, Tessella, SHA
const P = Tessella.API

# The blocked solves behind an element's interpolation matrix mirror Eigen, as
# Gmsh uses it: block sizes follow the processor's L1/L2/L3 cache capacities,
# and past 120 rows the blocking decides the order of the additions, so the
# last bits of the result follow the machine. Gmsh's digits on the machine it
# runs on are the reference, so that stays.
#
# What must not happen is a compiled and an interpreted run on one machine
# seeing different capacities. The capacities come from CPUID, which needs the
# compiler (`llvmcall`); an interpreted run (`--compile=min`, as this suite
# runs) used to take Eigen's defaults silently and now asks the operating
# system for the capacities of the processor it is on.
#
# On a processor with two kinds of cores the capacities depend on the core
# the thread is on when they are read (48 KiB and 32 KiB L1 data caches on a
# hybrid Intel part), compiled or not. The tests below therefore compare
# against the operating system's list for all processors.
const COMPILED = Base.JLOptions().compile_enabled in (1, 2)
const X86 = Sys.ARCH in (:x86_64, :i686)
const SMALL = (32768, 262144, 2097152)       # Eigen's defaults
const LARGE = (49152, 3145728, 31457280)     # a 48 KiB L1 data cache
matrix(n) = [((131i + 71j + (i * j) % 17) % 97) / 97 + (i == j ? n : 0) for i in 1:n, j in 1:n]
digest(a) = bytes2hex(sha256(reinterpret(UInt8, vec(a))))

@testset "Cache capacities agree with the operating system's" begin
    detected = P._native_cell_cache_sizes()
    @test all(>(0), detected)
    @test detected == P._native_cell_cache_sizes()          # read once per process
    reported = P._native_cell_os_cache_sizes()
    # Windows always answers; a Linux guest may not expose the cache directories.
    X86 && Sys.iswindows() && Sys.WORD_SIZE == 64 && @test reported !== nothing
    if X86 && reported !== nothing
        capacities = [P._native_cell_cache_defaults(last(entry)...) for entry in reported]
        @test issorted(first.(reported))
        # Compiled, `detected` is what CPUID said; interpreted, it is the
        # operating system's answer. Either way it is one of these.
        COMPILED && @test P._native_cell_cpuid(UInt32(0), UInt32(0))[1] > 0
        @test detected in capacities
    end
    @test P._native_cell_cache_defaults(0, 0, 0) == SMALL
    @test P._native_cell_cache_defaults(49152, 0, 31457280) == (49152, 262144, 31457280)
end

# With the capacities given, the arithmetic is the same compiled and
# interpreted: these digests were recorded in a compiled run.
@testset "Blocked triangular solve is fixed by the capacities" begin
    tri = matrix(125) ./ 64
    right_hand_sides() = [((17i + 5j) % 23) / 23 for i in 1:125, j in 1:3]
    # 125 rows exceed the 120-row block of a 32 KiB L1 cache and fit the
    # 184-row block of a 48 KiB one, so the two results differ in the last bits.
    @test P._native_cell_blocking(125, 125, 3, 4, SMALL)[1] == 64
    @test P._native_cell_blocking(125, 125, 3, 4, LARGE)[1] == 125
    small = P._native_cell_lower_solve!(tri, right_hand_sides(), SMALL)
    large = P._native_cell_lower_solve!(tri, right_hand_sides(), LARGE)
    @test digest(small) == "be0055f98707334c96fd8e4d5d5609cef4e3694de68d329a9fddf25ad2f7683f"
    @test digest(large) == "2453de85ba18ea38d657f303efb8587bcee3fbfa05dcc99a5bc0ccfeefea0f02"
    @test small[125, 3] == -0.2345338301449582
    @test large[125, 3] == -0.23453383014495807
    @test maximum(abs.(small .- large)) <= 1e-15
end

# The inverse of a 125-node element's matrix (a fourth-order hexahedron has
# 125 nodes) takes minutes interpreted, so the suite checks it compiled only.
# The same two digests were obtained once under `--compile=min -O0`.
@testset "Inverse of a 125-row matrix is fixed by the capacities" begin
    if COMPILED
        a = matrix(125)
        @test digest(P._native_cell_inverse(a, SMALL)) ==
              "54f6898573cf60f00124472cc58b5ed94851695d6f279414d03baf6783936138"
        @test digest(P._native_cell_inverse(a, LARGE)) ==
              "d14fb0f19eaa2d90f227b7d2eff0f4af083b5115736e47ea00ba7050d1008389"
    end
end
end
