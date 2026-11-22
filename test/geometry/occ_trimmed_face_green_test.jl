# The area of a trimmed tensor face: Green's theorem on the edges' own 3-D
# points taken back onto the patch, not a triangulation of the sampled
# pcurve polygon. A flat quadratic patch trimmed by a circle of radius 0.4
# measured 0.50206 (the projected pcurve's samples) where pi * 0.16 is
# 0.502654824574367; the faces of examples/boolean/component8.step moved to
# within 1e-7 of Gmsh 4.15.2's (face 6: 531.12126653614 vs 531.121317265642).
module OccTrimmedFaceGreenTest
using Test, Tessella
const M=Tessella.Model

@testset "trimmed tensor face area by Green's theorem" begin
    m=M.GeoModel()
    pts=Int[]
    for y in (-1.0,0.5,2.0), x in (-1.0,0.5,2.0)
        push!(pts,M._model_occ_add_point!(m,x,y,0.0))
    end
    c=M._model_occ_add_circle!(m,0.5,0.5,0.0,0.4)
    w=M._model_occ_add_wire!(m,[c])
    s=M._model_occ_add_bspline_surface!(m,pts,3;wire_tags=[w])
    @test M._model_occ_get_mass(m,2,s)≈π*0.16 rtol=1e-13
end
end
