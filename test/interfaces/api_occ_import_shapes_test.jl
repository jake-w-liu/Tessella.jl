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
        # a curved STEP (cylinder) imports through the structural path —
        # native cylinder face + rim circles + projected pcurves
        out3=occ.import_shapes(joinpath(FIX,"brep_occ_cylinder.step"))
        @test length(out3)==1
        d3,t3=out3[1]
        @test d3==3
        @test occ.get_mass(3,t3)>0.0
        @test all(isfinite,collect(occ.get_bounding_box(3,t3)))
    finally
        A.finalize()
    end
end

@testset "occ.import_shapes — STEP surface models bind dim-2 faces" begin
    A.initialize(); occ=A.model.occ
    try
        # `brep_occ_surface_sheet.step`: one ADVANCED_FACE under an
        # OPEN_SHELL/SHELL_BASED_SURFACE_MODEL — a 0.5×0.5 plane sheet at
        # z=0.5 over x,y∈[0.25,0.75]. Upstream `_multiBind` reaches the
        # surface tier when no transferable solid is present: each member
        # face binds as its own dim-2 entity.
        out=occ.import_shapes(joinpath(FIX,"brep_occ_surface_sheet.step"))
        @test length(out)==1 && out[1][1]==2
        stag=out[1][2]
        @test occ.get_mass(2,stag)≈0.25
        bb=collect(occ.get_bounding_box(2,stag))
        @test bb[1]≈0.25 && bb[4]≈0.75 && bb[3]≈bb[6]≈0.5

        # the imported sheet is a real dim-2 Boolean operand — the same
        # mixed path a native `add_rectangle` sheet takes. A sheet wholly
        # inside the box never partitions it (Extent 0 → preserved), and
        # the sheet stays bound at its import tag — identical to the
        # native-sheet disposition pinned in `occ.fragment — interior
        # and outside sheets` upstream.
        btag=occ.add_box(0,0,0,1,1,1)
        occ.synchronize()
        fout,fmap=occ.fragment([(3,btag)],[(2,stag)])
        @test fout==[(3,btag),(2,stag)]
        @test fmap==[[(3,btag)],[(2,stag)]]
        @test occ.get_mass(3,btag)≈1.0
        @test occ.get_mass(2,stag)≈0.25
        # intersect volume∩surface: the whole sheet lies inside the box —
        # the object produces no image and the sheet maps onto its own
        # rebound tag, like the native-sheet direction.
        cout,cmap=occ.intersect([(3,btag)],[(2,stag)])
        @test cmap==[Tuple{Int,Int}[],[(2,stag)]]
        @test cout==[(2,stag)]
        @test occ.get_mass(2,stag)≈0.25
    finally
        A.finalize()
    end
end
