using Test
using SHA: sha256
using Tessella

if !isdefined(Tessella, :TransfiniteCurve)
    Base.include(Tessella, joinpath(@__DIR__, "..", "src", "TransfiniteCurve.jl"))
end
using Tessella.TransfiniteCurve: transfinite_curve_parameters

function _parameter_crc(groups)
    stream = IOBuffer()
    for values in groups
        write(stream, htol(UInt64(length(values))))
        for value in values
            write(stream, htol(reinterpret(UInt64, value)))
        end
    end
    return bytes2hex(sha256(take!(stream)))
end

function _big_reference(mesh_type, coefficient, count)
    setprecision(BigFloat, 256) do
        segments = count - 1
        magnitude = abs(BigFloat(coefficient))
        reverse_orientation = coefficient < 0
        result = Vector{Float64}(undef, count)
        for k in 0:segments
            u = BigFloat(k) / segments
            value = if magnitude == 1 || (mesh_type === :beta && magnitude < 1)
                u
            elseif mesh_type === :progression
                ratio = reverse_orientation ? inv(magnitude) : magnitude
                (ratio^k - 1) / (ratio^segments - 1)
            elseif mesh_type === :bump
                if magnitude > 1
                    q = sqrt(magnitude - 1)
                    1//2 + tan((2u - 1) * atan(q)) / (2q)
                else
                    q = sqrt(1 - magnitude)
                    1//2 + tanh((2u - 1) * atanh(q)) / (2q)
                end
            else
                amplitude = log((magnitude + 1) / (magnitude - 1)) / 2
                reverse_orientation ?
                    magnitude * tanh(u * amplitude) :
                    1 - magnitude * tanh((1 - u) * amplitude)
            end
            result[k + 1] = Float64(value)
        end
        result[1] = 0.0
        result[end] = 1.0
        return result
    end
end

@noinline function _curve_allocated(count, mesh_type, coefficient)
    GC.gc()
    return @allocated transfinite_curve_parameters(
        count; mesh_type=mesh_type, coefficient=coefficient)
end

@noinline function _rejected_allocated()
    GC.gc()
    return @allocated try
        transfinite_curve_parameters(100_000; max_nodes=10)
    catch err
        err isa ArgumentError || rethrow()
    end
end

@testset "Gmsh straight transfinite curve laws" begin
    @testset "values, orientation, aliases, and deterministic CRC" begin
        progression = transfinite_curve_parameters(
            6; mesh_type=:progression, coefficient=2.0)
        expected_progression = Float64[
            0.0, 0.03225806451612903, 0.09677419354838712,
            0.22580645161290322, 0.4838709677419355, 1.0]
        @test progression == expected_progression
        @test transfinite_curve_parameters(
            6; mesh_type=:power, coefficient=2.0) == progression

        progression_reversed = transfinite_curve_parameters(
            6; mesh_type=:progression, coefficient=-2.0)
        @test progression_reversed ≈ 1.0 .- reverse(progression) atol=eps(Float64)
        @test transfinite_curve_parameters(
            6; mesh_type=:progression, coefficient=0.5) ==
              progression_reversed
        @test transfinite_curve_parameters(
            6; mesh_type=:progression, coefficient=-0.5) == progression

        bump = transfinite_curve_parameters(6; mesh_type=:bump, coefficient=2.0)
        @test bump == Float64[
            0.0, 0.2452372752527856, 0.4208077798377319,
            0.5791922201622681, 0.7547627247472144, 1.0]
        @test transfinite_curve_parameters(
            6; mesh_type=:bump, coefficient=-2.0) == bump
        @test bump == 1.0 .- reverse(bump)

        beta = transfinite_curve_parameters(6; mesh_type=:beta, coefficient=2.0)
        @test beta == Float64[
            0.0, 0.173631544092455, 0.36370669761585006,
            0.5674929709256276, 0.7811572746487717, 1.0]
        beta_reversed = transfinite_curve_parameters(
            6; mesh_type=:beta, coefficient=-2.0)
        @test beta_reversed == 1.0 .- reverse(beta)

        uniform = collect(range(0.0, 1.0; length=6))
        for mesh_type in (:progression, :bump, :beta), coefficient in (1.0, -1.0)
            @test transfinite_curve_parameters(
                6; mesh_type=mesh_type, coefficient=coefficient) == uniform
        end
        for coefficient in (0.25, -0.25, prevfloat(1.0), -prevfloat(1.0))
            @test transfinite_curve_parameters(
                6; mesh_type=:beta, coefficient=coefficient) == uniform
        end

        groups = (progression, progression_reversed, bump, beta, beta_reversed,
                  transfinite_curve_parameters(9; mesh_type=:bump,
                                                coefficient=0.05),
                  transfinite_curve_parameters(9; mesh_type=:beta,
                                                coefficient=1.01))
        @test _parameter_crc(groups) ==
              "b525628945d55f49bc5151ad313d697f9c32f7b3325015e5d0080a69260ffd0e"
        @test groups == (
            transfinite_curve_parameters(6; mesh_type=:progression, coefficient=2.0),
            transfinite_curve_parameters(6; mesh_type=:progression, coefficient=-2.0),
            transfinite_curve_parameters(6; mesh_type=:bump, coefficient=2.0),
            transfinite_curve_parameters(6; mesh_type=:beta, coefficient=2.0),
            transfinite_curve_parameters(6; mesh_type=:beta, coefficient=-2.0),
            transfinite_curve_parameters(9; mesh_type=:bump, coefficient=0.05),
            transfinite_curve_parameters(9; mesh_type=:beta, coefficient=1.01))
    end

    @testset "independent high-precision primitive inverses" begin
        cases = ((:progression, 0.7), (:progression, -0.7),
                 (:progression, 1.3), (:progression, -1.3),
                 (:bump, 0.05), (:bump, 0.5), (:bump, 2.0),
                 (:bump, 20.0), (:bump, -2.0),
                 (:beta, 0.5), (:beta, -0.5), (:beta, 1.01),
                 (:beta, 1.2), (:beta, 2.0), (:beta, -2.0),
                 (:beta, 20.0))
        for count in (3, 6, 17), (mesh_type, coefficient) in cases
            actual = transfinite_curve_parameters(
                count; mesh_type=mesh_type, coefficient=coefficient)
            reference = _big_reference(mesh_type, coefficient, count)
            @test actual ≈ reference atol=16eps(Float64) rtol=16eps(Float64)
            @test actual[1] === 0.0
            @test actual[end] === 1.0
            @test all(isfinite, actual)
            @test all(>(0.0), diff(actual))
        end

        ratio_case = transfinite_curve_parameters(
            17; mesh_type=:progression, coefficient=1.3)
        lengths = diff(ratio_case)
        @test all(isapprox(lengths[index + 1] / lengths[index], 1.3;
                           atol=64eps(Float64), rtol=64eps(Float64))
                  for index in 1:length(lengths)-1)

        # These coefficients exercise both cancellation-sensitive sides of the
        # Bump law while still having distinct Float64 reference parameters.
        for (coefficient, count) in ((1e-15, 257), (1e30, 17))
            actual = transfinite_curve_parameters(
                count; mesh_type=:bump, coefficient=coefficient)
            reference = _big_reference(:bump, coefficient, count)
            @test actual ≈ reference atol=64eps(Float64) rtol=64eps(Float64)
            @test all(>(0.0), diff(actual))
        end
    end

    @testset "validation, resource limits, and representability blockers" begin
        @test transfinite_curve_parameters(2) == [0.0, 1.0]
        @test transfinite_curve_parameters(
            2; mesh_type=:bump, coefficient=floatmax(Float64), max_nodes=2) ==
              [0.0, 1.0]

        for count in (true, 1, 0, -1, 2.0, big(typemax(Int32)) + 1)
            @test_throws ArgumentError transfinite_curve_parameters(count)
        end
        @test_throws ArgumentError transfinite_curve_parameters(3; max_nodes=2)
        @test_throws ArgumentError transfinite_curve_parameters(3; max_nodes=-1)
        @test_throws ArgumentError transfinite_curve_parameters(3; max_nodes=true)
        @test_throws ArgumentError transfinite_curve_parameters(3; max_nodes=3.0)
        @test_throws ArgumentError transfinite_curve_parameters(
            3; max_nodes=big(typemax(Int32)) + 1)

        for mesh_type in (:Progression, :unknown, "Progression", nothing)
            @test_throws ArgumentError transfinite_curve_parameters(
                3; mesh_type=mesh_type)
        end
        for coefficient in (true, 0.0, -0.0, NaN, Inf, -Inf, 1 + 2im,
                            "2", BigFloat(2)^20_000)
            @test_throws ArgumentError transfinite_curve_parameters(
                3; coefficient=coefficient)
        end

        for coefficient in (floatmin(Float64), floatmax(Float64))
            @test_throws ArgumentError transfinite_curve_parameters(
                100; mesh_type=:progression, coefficient=coefficient)
            @test_throws ArgumentError transfinite_curve_parameters(
                100; mesh_type=:bump, coefficient=coefficient)
        end
        @test transfinite_curve_parameters(
            100; mesh_type=:beta, coefficient=floatmin(Float64)) ==
              [k / 99 for k in 0:99]
        @test all(>(0.0), diff(transfinite_curve_parameters(
            100; mesh_type=:beta, coefficient=floatmax(Float64))))

        transfinite_curve_parameters(100_000; max_nodes=100_000)
        rejected = _rejected_allocated()
        @test rejected < 64_000
        @test isempty(Test.detect_ambiguities(
            Tessella.TransfiniteCurve; recursive=true))
    end

    @testset "allocation growth is linear in the returned vector" begin
        for (mesh_type, coefficient) in
            ((:progression, 1.00001), (:bump, 0.5), (:beta, 2.0))
            transfinite_curve_parameters(
                128; mesh_type=mesh_type, coefficient=coefficient)
            small = _curve_allocated(20_000, mesh_type, coefficient)
            large = _curve_allocated(40_000, mesh_type, coefficient)
            @test small > 0
            @test large > small
            @test large <= 2.15small + 65_536
            @info "transfinite curve allocation ratchet" mesh_type small large
        end
    end
end
