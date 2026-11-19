module OCCBooleanTorusCavityChecks
using Test,Tessella
const A=Tessella.API
const M=Tessella.Model

# A solid bounded by one torus face and plane faces was always measured as a
# torus sector (Pappus), so a box with a toroidal cavity — what
# `occ.cut(box, torus inside it)` returns — reported the cavity's volume
# (0.10) instead of its own (7.90). The solid is assembled here from the two
# primitives' own shells, the torus shell reversed, so the check needs no
# Boolean (a box×torus Boolean traces six sections).
@testset "OCC volume of a box with a toroidal cavity" begin
    torus_volume=2π^2*0.26*0.14^2
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        occ.add_box(-1,-1,-1,2,2,2)
        occ.add_torus(0.26,0,0,0.26,0.14)
        # sector solids keep the closed form
        @test occ.get_mass(3,2)≈torus_volume rtol=1e-12
        kernel=A.CURRENT[].occ_kernel
        box=kernel.shapes[kernel.bindings[(3,1)]]
        torus=kernel.shapes[kernel.bindings[(3,2)]]
        shell=kernel.shapes[only(torus.children)]
        cavity=kernel.next_shape;solid=cavity+1
        kernel.shapes[cavity]=M._ModelOCCShape(-2,
            (uses=[(face,!reversed) for (face,reversed) in shell.geometry.uses],),
            copy(shell.children))
        kernel.shapes[solid]=M._ModelOCCShape(3,(kind=:solid,),
            vcat(box.children,[cavity]))
        # volume only: the torus wall then reduces to boundary integrals
        volume,_=M._model_occ_solid_measure(kernel,solid;want_centroid=false)
        @test volume≈8-torus_volume rtol=1e-6
    finally
        A.finalize()
    end
end
end
