# `model.occ.import_shapes` — ISO-10303-21 STEP solids bound as native
# OCC records. Planar fixtures only: shared vertices/edges dedupe so the
# shell's two-surfaces-per-curve rule holds by construction.
using Test
using Tessella
const A=Tessella.API
const FIX=joinpath(@__DIR__,"..","fixtures")

function _errtext(f)
    try
        f()
    catch err
        return sprint(showerror,err)
    end
    return ""
end

@testset "occ.import_shapes" begin
    A.initialize(); occ=A.model.occ
    try
        # ── OCC-authored unit cube STEP ─────────────────────────────
        out=occ.import_shapes(joinpath(FIX,"brep_occ_box.step"))
        @test length(out)==1
        d,t=out[1]
        @test d==3
        @test occ.get_mass(3,t)>0.0
        bb=collect(occ.get_bounding_box(3,t))
        @test all(isfinite,bb)

        # ── planar through-hole (BREP_WITH_VOIDS-style fixture) ─────
        out2=occ.import_shapes(joinpath(FIX,"brep_occ_planar_through_hole.step"))
        @test !isempty(out2)
        for (d2,t2) in out2
            @test occ.get_mass(d2,t2)>0.0
        end

        # ── imported volume is a real OCC solid: it booleans ────────
        d1,t1=out[1]
        clip=occ.add_box(bb[1]+(bb[4]-bb[1])/2,bb[2],bb[3],
                         bb[4]-bb[1],bb[5]-bb[2],bb[6]-bb[3])
        cut_out,_=occ.cut([(3,t1)],[(3,clip)])
        @test !isempty(cut_out)

        # ── errors ──────────────────────────────────────────────────
        @test occursin("missing file",
                       _errtext(()->occ.import_shapes("no_such_file_xyz.step")))
        @test occursin("ISO-10303-21",
                       _errtext(()->occ.import_shapes(joinpath(FIX,"..","..","Project.toml"))))
        # a curved STEP (cylinder) rejects loudly, not silently
        @test occursin("planar",
                       _errtext(()->occ.import_shapes(joinpath(FIX,"brep_occ_cylinder.step"))))
    finally
        A.finalize()
    end
end
