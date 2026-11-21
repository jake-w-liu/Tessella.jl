# The suite includes every test file into Main; this file keeps its aliases to itself.
module CurveFaceIndexTests
using Test
using Tessella
import Tessella.Model as M

# The faces of a Curve are read from an index built once per state of the
# model's Surface and Loop stores. It must give what a scan of every Surface
# gives, in the same order, after every kind of edit to those stores.

function agrees(m)
    curves = sort!(collect(keys(m.curves)))
    push!(curves, maximum(curves; init=0) + 7)      # a tag that bounds nothing
    return all(c -> M._model_curve_faces(m, c) == M._model_curve_faces_scan(m, c), curves)
end

# An n by n grid of unit squares; returns the tag of square (i,j).
function grid!(m, n)
    point(i, j) = i * (n + 1) + j + 1
    for i in 0:n, j in 0:n
        M.add_point!(m, Float64(i), Float64(j), 0.0; tag=point(i, j))
    end
    hline(i, j) = i * n + j + 1
    vline(i, j) = (n + 1) * n + i * (n + 1) + j + 1
    for i in 0:n, j in 0:n-1
        M.add_line!(m, point(i, j), point(i, j + 1); tag=hline(i, j))
    end
    for i in 0:n-1, j in 0:n
        M.add_line!(m, point(i, j), point(i + 1, j); tag=vline(i, j))
    end
    square = 0
    for i in 0:n-1, j in 0:n-1
        square += 1
        M.add_curve_loop!(m, [hline(i, j), vline(i, j + 1), -hline(i + 1, j), -vline(i, j)]; tag=square)
        M.add_plane_surface!(m, [square]; tag=square)
    end
    return hline, vline
end

@testset "Curve faces come from an index that follows the model" begin
    M._MODEL_DICT_COUNTS_EDITS || @warn "Dict has no edit counter on this Julia: the scan is used"
    m = M.GeoModel()
    @test isempty(M._model_curve_faces(m, 1))
    hline, vline = grid!(m, 4)
    @test agrees(m)
    # An interior line bounds two squares, a border line one.
    @test length(M._model_curve_faces(m, hline(2, 1))) == 2
    @test length(M._model_curve_faces(m, hline(0, 1))) == 1
    @test length(M._model_curve_faces(m, vline(1, 0))) == 1
    # Repeated queries reuse one index.
    index = M._model_curve_face_index(m)
    @test M._model_curve_face_index(m) === index

    # A Surface added afterwards, with a hole made of another square's loop.
    outer = M.add_curve_loop!(m, [hline(0, 0), hline(0, 1), hline(0, 2), hline(0, 3),
                                  vline(0, 4), vline(1, 4), vline(2, 4), vline(3, 4),
                                  -hline(4, 3), -hline(4, 2), -hline(4, 1), -hline(4, 0),
                                  -vline(3, 0), -vline(2, 0), -vline(1, 0), -vline(0, 0)])
    big = M.add_plane_surface!(m, [outer, 6])
    @test agrees(m)
    @test big in M._model_curve_faces(m, hline(0, 0))
    @test big in M._model_curve_faces(m, hline(1, 1))     # through the hole loop
    @test M._model_curve_face_index(m) !== index

    # A loop replaced in place of the old one, and a Surface given other loops.
    m.loops[3] = [hline(0, 0), vline(0, 1), -hline(1, 0), -vline(0, 0)]
    @test agrees(m)
    m.surfaces[5] = [5, 7]
    @test agrees(m)

    # Surfaces and loops removed.
    delete!(m.surfaces, big)
    @test agrees(m)
    @test !(big in M._model_curve_faces(m, hline(0, 0)))
    delete!(m.loops, 9)
    @test agrees(m)

    # A second model has its own index; a shallow copy sharing the stores agrees too.
    other = M.GeoModel()
    grid!(other, 2)
    @test agrees(other) && agrees(m)
    shared = M.GeoModel(ntuple(i -> getfield(m, i), Val(fieldcount(M.GeoModel)))...)
    @test agrees(shared)
    shared.surfaces[99] = [10]
    @test agrees(shared) && agrees(m)

    # The store itself replaced.
    m.surfaces = copy(m.surfaces)
    delete!(m.surfaces, 1)
    @test agrees(m)
    empty!(m.surfaces)
    @test agrees(m)
    @test isempty(M._model_curve_faces(m, hline(2, 1)))
end

# The straight-Line prepass shortcut answers only where the integral cannot be zero.
@testset "A straight Line is known to have length without integrating it" begin
    m = M.GeoModel()
    for (tag, xyz) in enumerate(((0.0, 0.0, 0.0), (1.0, 2.0, 3.0), (1e3, 0.0, 0.0), (1e3 + 0.5, 0.0, 0.0),
                                 (1e-3, 0.0, 0.0), (1e120, 0.0, 0.0), (5.0, 5.0, 5.0), (5.0, 5.0, 5.0),
                                 (1e6, 0.0, 0.0), (1e6 + 1e-3, 0.0, 0.0), (1e-130, 0.0, 0.0)))
        M.add_point!(m, xyz...; tag)
    end
    line(a, b, tag) = (m.curves[tag] = (a, b); tag)       # distinct tags at one position are legal
    sure(tag) = M._model_line_surely_has_length(m, tag, m.curves[tag]..., 0.0, 1.0)
    @test sure(line(1, 2, 1))
    @test sure(line(3, 4, 2))                   # chord 0.5 at 1e3
    @test sure(line(1, 5, 3))                   # a short Line near the origin
    @test !sure(line(7, 8, 4))                  # coincident ends: zero length
    @test !sure(line(9, 10, 5))                 # 1e-3 at 1e6: too short for the stencil to be trusted
    @test !sure(line(1, 6, 6))                  # extreme coordinates
    @test !sure(line(1, 11, 7))                 # near the underflow range
    @test !M._model_line_surely_has_length(m, 1, 1, 2, 0.0, 0.5)   # not the Line's own range
    # Where it answers, the integral it replaces is indeed positive.
    options = M._model_default_mesh1d_options(m, "test")
    for tag in (1, 2, 3)
        @test M._model_curve_prepass_length(m, tag, 0.0, 1.0, options, "test") > 0
    end
    # An OpenCASCADE line: its own range, a positive direction, comparable magnitudes.
    @test M._model_occ_line_surely_has_length((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 0.0, 2.0, 0.0, 2.0)
    @test !M._model_occ_line_surely_has_length((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 0.0, 2.0, 0.0, 1.0)
    @test !M._model_occ_line_surely_has_length((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), 0.0, 2.0, 0.0, 2.0)
    @test !M._model_occ_line_surely_has_length((1e9, 0.0, 0.0), (0.0, 0.0, 1.0), 0.0, 1.0, 0.0, 1.0)
end
end
