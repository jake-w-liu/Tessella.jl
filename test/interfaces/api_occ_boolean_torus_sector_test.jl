module OCCBooleanTorusSectorChecks
using Test,Tessella
const A=Tessella.API
const M=Tessella.Model

# A plane through a torus' axis leaves a genuine half-turn sector whose
# carrier record still says `angle=2π`. The face area and the face box used
# the carrier's sweep instead of the wire's window, so the half torus
# reported the whole torus' area and a box reaching across the cutting plane.
# The face is built here from a half-turn primitive whose record is given the
# full-turn carrier angle — the state `occ.intersect(box, torus)` leaves —
# so the check needs no Boolean.
@testset "OCC torus sector narrower than its carrier's sweep" begin
    R=1.0;r=0.3
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        occ.add_torus(0,0,0,R,r;angle=π)
        kernel=A.CURRENT[].occ_kernel
        id=only(id for ((d,_),id) in kernel.bindings if d==2 &&
            get(kernel.shapes[id].geometry.geometry,:occ,nothing)===:torus)
        face=kernel.shapes[id]
        carrier=merge(face.geometry.geometry,(angle=2π,))
        kernel.shapes[id]=M._ModelOCCShape(2,
            merge(face.geometry,(geometry=carrier,)),face.children)
        window=M._occ_torus_sector_window(kernel,id,carrier)
        @test collect(window)≈[0,π] atol=1e-12
        area,centroid=M._model_occ_face_measure(kernel,id)
        @test area≈2π^2*R*r rtol=1e-12
        @test centroid[2]>0.5
        @test collect(M._model_occ_shape_bbox(kernel,id))≈
            [-R-r,0,-r,R+r,R+r,r] atol=1e-12
    finally
        A.finalize()
    end
end
end
