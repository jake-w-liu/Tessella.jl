using Test
using SHA: sha256
using Tessella

if !isdefined(Tessella, :TransfiniteCurve)
    Base.include(Tessella, joinpath(@__DIR__, "..", "src", "TransfiniteCurve.jl"))
end
using Tessella.TransfiniteCurve: transfinite_curve_parameters,
                                 transfinite_curve_hwall

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

function _big_hwall_reference(count, wall_height, curve_length;
                              orientation=:start)
    setprecision(BigFloat, 256) do
        segments = count - 1
        fraction = BigFloat(wall_height) / BigFloat(curve_length)
        target = inv(fraction)
        ratio = if fraction == inv(BigFloat(segments))
            BigFloat(1)
        else
            geometric_sum(r) = sum(r^k for k in 0:segments-1)
            lo = BigFloat(1)
            hi = BigFloat(2)
            while geometric_sum(hi) < target
                hi *= 2
            end
            for _ in 1:512
                mid = (lo + hi) / 2
                if geometric_sum(mid) < target
                    lo = mid
                else
                    hi = mid
                end
            end
            (lo + hi) / 2
        end

        result = Vector{Float64}(undef, count)
        result[1] = 0.0
        accumulated = BigFloat(0)
        segment = fraction
        for index in 2:count-1
            accumulated += segment
            result[index] = Float64(accumulated)
            segment *= ratio
        end
        result[end] = 1.0
        if orientation === :end
            result = 1.0 .- reverse(result)
            result[1] = 0.0
            result[end] = 1.0
        end
        return result
    end
end

function _big_bump_hwall_reference(count, wall_height, curve_length)
    setprecision(BigFloat, 256) do
        wall_fraction = BigFloat(wall_height) / BigFloat(curve_length)
        primitive(x) = if iszero(x)
            count * wall_fraction
        else
            coefficient = exp(-x)
            q = sqrt(1 - coefficient)
            amplitude = atanh(q)
            delta = atanh(2q * wall_fraction /
                           (coefficient + 2q^2 * wall_fraction))
            count * delta / (2amplitude)
        end
        lo = BigFloat(0)
        hi = BigFloat(1)
        while primitive(hi) < 1
            hi *= 2
        end
        for _ in 1:512
            mid = (lo + hi) / 2
            if primitive(mid) < 1
                lo = mid
            else
                hi = mid
            end
        end
        coefficient = exp(-(lo + hi) / 2)
        q = sqrt(1 - coefficient)
        amplitude = atanh(q)
        segments = count - 1
        result = Vector{Float64}(undef, count)
        for k in 0:segments
            u = BigFloat(k) / segments
            result[k + 1] = Float64(
                1//2 + tanh((2u - 1) * amplitude) / (2q))
        end
        result[1] = 0.0
        result[end] = 1.0
        return result
    end
end

function _big_beta_hwall_reference(count, wall_height, curve_length;
                                   orientation=:start)
    setprecision(BigFloat, 256) do
        wall_fraction = BigFloat(wall_height) / BigFloat(curve_length)
        count_ratio = BigFloat(count - 1) / count
        equation(amplitude) = log1p(-wall_fraction) + log(tanh(amplitude)) -
                              log(tanh(count_ratio * amplitude))
        lo = BigFloat(0)
        hi = BigFloat(1)
        while equation(hi) > 0
            hi *= 2
        end
        for _ in 1:512
            mid = (lo + hi) / 2
            if equation(mid) > 0
                lo = mid
            else
                hi = mid
            end
        end
        amplitude = (lo + hi) / 2
        segments = count - 1
        result = Vector{Float64}(undef, count)
        for k in 0:segments
            u = BigFloat(k) / segments
            result[k + 1] = Float64(
                sinh(u * amplitude) /
                (sinh(amplitude) * cosh((1 - u) * amplitude)))
        end
        result[1] = 0.0
        result[end] = 1.0
        if orientation === :end
            result = 1.0 .- reverse(result)
            result[1] = 0.0
            result[end] = 1.0
        end
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

@noinline function _hwall_allocated(count, mesh_type)
    GC.gc()
    return @allocated transfinite_curve_hwall(
        count; mesh_type=mesh_type, wall_height=1 / (2count), curve_length=1.0,
        max_nodes=count)
end

@noinline function _hwall_rejected_allocated()
    GC.gc()
    return @allocated try
        transfinite_curve_hwall(
            100_000; wall_height=1e-6, curve_length=1.0, max_nodes=10)
    catch err
        err isa ArgumentError || rethrow()
    end
end


@noinline function _hwall_nominal_rejected_allocated(mesh_type)
    GC.gc()
    return @allocated try
        transfinite_curve_hwall(
            1_000_000; mesh_type=mesh_type, wall_height=0.5,
            curve_length=1.0, max_nodes=1_000_000)
    catch err
        err isa ArgumentError || rethrow()
    end
end

@noinline function _hwall_infeasible_allocated()
    GC.gc()
    return @allocated try
        transfinite_curve_hwall(
            1_000_000; wall_height=0.5, curve_length=1.0,
            max_nodes=1_000_000)
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

    @testset "Progression HWall law and independent geometric-sum oracle" begin
        exact = transfinite_curve_hwall(
            6; wall_height=1.0, curve_length=31.0)
        @test exact ≈ [0.0, 1/31, 3/31, 7/31, 15/31, 1.0]
        @test all(isapprox(diff(exact)[index + 1] / diff(exact)[index], 2.0;
                           rtol=64eps(Float64), atol=0.0)
                  for index in 1:4)

        mirrored = transfinite_curve_hwall(
            6; wall_height=1.0, curve_length=31.0, orientation=:end)
        @test mirrored ≈ 1.0 .- reverse(exact) atol=8eps(Float64)
        @test mirrored[end] - mirrored[end-1] ≈ 1/31 rtol=64eps(Float64)
        @test transfinite_curve_hwall(
            6; wall_height=2.0, curve_length=10.0) ==
              collect(range(0.0, 1.0; length=6))
        @test transfinite_curve_hwall(
            2; wall_height=3.5, curve_length=3.5) == [0.0, 1.0]
        @test transfinite_curve_hwall(
            3; wall_height=0.01, curve_length=1.0) ≈ [0.0, 0.01, 1.0]

        cases = ((3, 0.125, 1.0), (4, 1.0, 13.0),
                 (6, 0.1, 1.0), (17, 0.01, 1.0),
                 (65, prevfloat(1 / 64), 1.0))
        groups = Vector{Vector{Float64}}()
        for (count, wall_height, curve_length) in cases,
            orientation in (:start, :end)
            actual = transfinite_curve_hwall(
                count; wall_height=wall_height, curve_length=curve_length,
                orientation=orientation, max_nodes=count)
            reference = _big_hwall_reference(
                count, wall_height, curve_length; orientation=orientation)
            @test actual ≈ reference atol=32eps(Float64) rtol=128eps(Float64)
            @test actual[1] === 0.0
            @test actual[end] === 1.0
            @test all(isfinite, actual)
            @test all(>(0.0), diff(actual))
            wall_parameter = orientation === :start ? actual[2] : 1 - actual[end-1]
            @test wall_parameter ≈ wall_height / curve_length rtol=256eps(Float64)
            push!(groups, actual)
        end

        # The normalized law is invariant under a common physical-length scale.
        @test transfinite_curve_hwall(
            6; wall_height=1e299 / 31, curve_length=1e299) ≈ exact
        @test transfinite_curve_hwall(
            6; wall_height=1e-299 / 31, curve_length=1e-299) ≈ exact
        @test Tessella.transfinite_curve_hwall === transfinite_curve_hwall
        @test _parameter_crc(groups) ==
              "802ae6dd95259c50b087d03e7b7567b555f6040f8e62b1a2afc3aed6bca22379"
    end

    @testset "Bump/Beta HWall laws and high-precision primitive oracles" begin
        bump = transfinite_curve_hwall(
            6; mesh_type=:bump, wall_height=0.05, curve_length=1.0)
        @test bump ≈ [0.0, 0.06918509550105913, 0.3047609650010744,
                       0.6952390349989256, 0.9308149044989409, 1.0]
        @test bump ≈ 1.0 .- reverse(bump) atol=eps(Float64)
        @test transfinite_curve_hwall(
            6; mesh_type=:bump, wall_height=0.05, curve_length=1.0,
            orientation=:end) == bump

        beta = transfinite_curve_hwall(
            6; mesh_type=:beta, wall_height=0.05, curve_length=1.0)
        @test beta ≈ [0.0, 0.06343392236633141, 0.17768169506999182,
                       0.3685407841246565, 0.6505836059034631, 1.0]
        beta_end = transfinite_curve_hwall(
            6; mesh_type=:beta, wall_height=0.05, curve_length=1.0,
            orientation=:end)
        @test beta_end ≈ 1.0 .- reverse(beta) atol=8eps(Float64)
        for mesh_type in (:bump, :beta)
            normalized = transfinite_curve_hwall(
                6; mesh_type=mesh_type, wall_height=0.05,
                curve_length=1.0)
            @test transfinite_curve_hwall(
                6; mesh_type=mesh_type, wall_height=5e297,
                curve_length=1e299) ≈ normalized
            @test transfinite_curve_hwall(
                6; mesh_type=mesh_type, wall_height=5e-301,
                curve_length=1e-299) ≈ normalized
        end

        cases = ((6, 0.01, 1.0), (6, 0.05, 1.0), (6, 0.1, 1.0),
                 (17, 0.01, 1.0), (17, 0.05, 1.0),
                 (65, 0.005, 1.0))
        groups = Vector{Vector{Float64}}()
        for (count, wall_height, curve_length) in cases
            actual_bump = transfinite_curve_hwall(
                count; mesh_type=:bump, wall_height=wall_height,
                curve_length=curve_length, max_nodes=count)
            reference_bump = _big_bump_hwall_reference(
                count, wall_height, curve_length)
            @test actual_bump ≈ reference_bump atol=32eps(Float64) rtol=128eps(Float64)
            @test actual_bump ≈ 1.0 .- reverse(actual_bump) atol=eps(Float64)
            @test all(>(0.0), diff(actual_bump))
            push!(groups, actual_bump)

            for orientation in (:start, :end)
                actual_beta = transfinite_curve_hwall(
                    count; mesh_type=:beta, wall_height=wall_height,
                    curve_length=curve_length, orientation=orientation,
                    max_nodes=count)
                reference_beta = _big_beta_hwall_reference(
                    count, wall_height, curve_length; orientation=orientation)
                @test actual_beta ≈ reference_beta atol=32eps(Float64) rtol=128eps(Float64)
                @test all(>(0.0), diff(actual_beta))
                push!(groups, actual_beta)
            end
        end

        @test transfinite_curve_hwall(
            6; mesh_type=:power, wall_height=0.05, curve_length=1.0) ==
              transfinite_curve_hwall(
                  6; mesh_type=:progression, wall_height=0.05,
                  curve_length=1.0)
        for mesh_type in (:bump, :beta), count in (3, 6, 17)
            near_limit = transfinite_curve_hwall(
                count; mesh_type=mesh_type,
                wall_height=prevfloat(1.0 / count), curve_length=1.0)
            @test isapprox(
                near_limit, collect(range(0.0, 1.0; length=count));
                atol=16eps(Float64), rtol=16eps(Float64))
            @test all(>(0.0), diff(near_limit))
        end
        @test _parameter_crc(groups) ==
              "04169f75cdcfabf540e88477ce54b55d0a6eabcfd5ca48f19b332afaac0fb59a"
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

    @testset "HWall validation, resource, and representability blockers" begin
        @test_throws ArgumentError transfinite_curve_hwall(
            6; curve_length=1.0)
        @test_throws ArgumentError transfinite_curve_hwall(
            6; wall_height=0.1)
        for value in (true, 0.0, -0.0, -1.0, NaN, Inf, -Inf, 1 + 2im,
                      "0.1", BigFloat(2)^20_000)
            @test_throws ArgumentError transfinite_curve_hwall(
                6; wall_height=value, curve_length=1.0)
            @test_throws ArgumentError transfinite_curve_hwall(
                6; wall_height=0.1, curve_length=value)
        end
        for orientation in (:middle, "start", nothing, true)
            @test_throws ArgumentError transfinite_curve_hwall(
                6; wall_height=0.1, curve_length=1.0,
                orientation=orientation)
        end
        for mesh_type in (:Progression_HWall, :unknown, "bump", nothing, true)
            @test_throws ArgumentError transfinite_curve_hwall(
                6; mesh_type=mesh_type, wall_height=0.05, curve_length=1.0)
        end
        for count in (true, 1, 0, -1, 2.0, big(typemax(Int32)) + 1)
            @test_throws ArgumentError transfinite_curve_hwall(
                count; wall_height=0.1, curve_length=1.0)
        end
        @test_throws ArgumentError transfinite_curve_hwall(
            6; wall_height=nextfloat(0.2), curve_length=1.0)
        @test_throws ArgumentError transfinite_curve_hwall(
            2; wall_height=0.5, curve_length=1.0)
        for mesh_type in (:bump, :beta)
            @test_throws ArgumentError transfinite_curve_hwall(
                2; mesh_type=mesh_type, wall_height=0.1, curve_length=1.0)
            @test_throws ArgumentError transfinite_curve_hwall(
                6; mesh_type=mesh_type, wall_height=1 / 6,
                curve_length=1.0)
            @test_throws ArgumentError transfinite_curve_hwall(
                6; mesh_type=mesh_type, wall_height=0.2, curve_length=1.0)
            @test _hwall_nominal_rejected_allocated(mesh_type) < 64_000
        end
        @test_throws ArgumentError transfinite_curve_hwall(
            6; wall_height=floatmin(Float64), curve_length=floatmax(Float64))
        @test_throws ArgumentError transfinite_curve_hwall(
            3; wall_height=nextfloat(0.0), curve_length=nextfloat(0.0))
        @test_throws ArgumentError transfinite_curve_hwall(
            6; wall_height=1e-20, curve_length=1.0, orientation=:end)
        @test_throws ArgumentError transfinite_curve_hwall(
            6; mesh_type=:bump, wall_height=1e-20, curve_length=1.0)
        @test_throws ArgumentError transfinite_curve_hwall(
            6; mesh_type=:beta, wall_height=1e-20, curve_length=1.0,
            orientation=:end)
        beta_extreme = transfinite_curve_hwall(
            6; mesh_type=:beta, wall_height=1e-20, curve_length=1.0)
        @test all(>(0.0), diff(beta_extreme))
        @test_throws ArgumentError transfinite_curve_hwall(
            3; wall_height=0.1, curve_length=1.0, max_nodes=2)
        @test _hwall_rejected_allocated() < 64_000
        @test _hwall_infeasible_allocated() < 64_000
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

        for mesh_type in (:progression, :bump, :beta)
            transfinite_curve_hwall(
                128; mesh_type=mesh_type, wall_height=1 / 256,
                curve_length=1.0, max_nodes=128)
            hwall_small = _hwall_allocated(20_000, mesh_type)
            hwall_large = _hwall_allocated(40_000, mesh_type)
            @test hwall_small > 0
            @test hwall_large > hwall_small
            @test hwall_large <= 2.15hwall_small + 65_536
            @info "transfinite curve HWall allocation ratchet" mesh_type hwall_small hwall_large
        end
    end
end
