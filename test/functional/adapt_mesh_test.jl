# Native adaptive remeshing — `AdaptMesh` (`.geo`) and CLI batch mode `-4`.
# Bounded Coupez-style loop: metric edge splits via `refine_to_size`, interior
# edge collapse, embedded-boundary-protected quality swaps — exercising the
# actual authority carrier, entity ownership and physical-tag transport.

const _ADAPT_BOX_GEO = """
lc = 1.3;
Point(1) = {0,0,0,lc}; Point(2) = {1,0,0,lc}; Point(3) = {1,1,0,lc}; Point(4) = {0,1,0,lc};
Point(5) = {0,0,1,lc}; Point(6) = {1,0,1,lc}; Point(7) = {1,1,1,lc}; Point(8) = {0,1,1,lc};
Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
Line(5)={5,6};Line(6)={6,7};Line(7)={7,8};Line(8)={8,5};
Line(9)={1,5};Line(10)={2,6};Line(11)={3,7};Line(12)={4,8};
Curve Loop(1)={1,2,3,4}; Plane Surface(1)={1};
Curve Loop(2)={-7,-6,-5,-8}; Plane Surface(2)={2};
Curve Loop(3)={9,5,-10,-1}; Plane Surface(3)={3};
Curve Loop(4)={10,6,-11,-2}; Plane Surface(4)={4};
Curve Loop(5)={11,7,-12,-3}; Plane Surface(5)={5};
Curve Loop(6)={12,8,-9,-4}; Plane Surface(6)={6};
Surface Loop(1)={1,2,3,4,5,6}; Volume(1)={1};
Physical Surface(10)={1}; Physical Volume(20)={1};
Mesh 3;
Field[1] = MathEval;
Field[1].F = "0.35 + 0.3*x";
Background Field = 1;
"""

function _adapt_write_geo(body::AbstractString)
    path=joinpath(tempdir(),"tessella_adapt_"*string(hash(body),base=16)*".geo")
    open(path,"w") do io
        print(io,_ADAPT_BOX_GEO)
        print(io,body)
    end
    return path
end

# Physical bindings on the actual (authority) mesh product — per-entity tags.
_adapt_physical_tags(am)=Dict((d,t)=>copy(e.physical_tags)
    for ((d,t),e) in am.entity_data.entities)

@testset "adaptive remeshing" begin
    # ── mesh-level operators: split + collapse + swap preserve validity and tags ──
    r0=Tessella.execute_geo(_adapt_write_geo(""))
    carrier=r0.mesh
    field=Tessella.FunctionSize((x,y,z)->0.35+0.3*x)
    adapted,stats=Tessella.Mesh3D.adapt_mesh(carrier,field;niter=4)
    @test stats.splits>0
    @test stats.collapses+stats.swaps+stats.splits>0
    @test stats.unresolved_edges==0
    @test size(adapted.tets,2)>size(carrier.tets,2)
    @test Tessella.validate(adapted).ok
    # metric postcondition: no edge exceeds the field target
    @test all(eachcol(adapted.tets)) do v
        all(Tessella.metric_edge_length(
            field,(adapted.coords[1,a],adapted.coords[2,a],adapted.coords[3,a]),
            (adapted.coords[1,b],adapted.coords[2,b],adapted.coords[3,b]))<=1+1e-9
            for (a,b) in ((v[1],v[2]),(v[1],v[3]),(v[1],v[4]),(v[2],v[3]),(v[2],v[4]),(v[3],v[4])))
    end

    # ── .geo AdaptMesh: entity ownership + physical tags survive rebind ──
    r=Tessella.execute_geo(_adapt_write_geo("AdaptMesh {1} {} {{}} {2, 0};\n"))
    m=r.mesh
    @test m isa Tessella.Mesh
    @test size(m.tets,2)>0
    @test Tessella.validate(m).ok
    # owner tags bound per family: 12 curves, 6 surfaces, 1 volume
    @test sort!(unique(m.seg_tag))==Int32.(1:12)
    @test sort!(unique(m.tri_tag))==Int32.(1:6)
    @test sort!(unique(m.tet_tag))==Int32[1]
    auth=r.model.meshing.mesh_authority
    @test auth!==nothing && auth.mesh isa Tessella.MixedMesh
    @test auth.carrier===m
    @test length(keys(auth.mesh.entity_data.entities))==27
    phys=_adapt_physical_tags(auth.mesh)
    @test phys[(3,1)]==Int32[20]          # Physical Volume binding survives
    @test phys[(2,1)]==Int32[10]          # Physical Surface binding survives
    # graded output: tetrahedra grew where the field is fine
    @test size(m.tets,2)>size(carrier.tets,2)

    # ── loud bounded rejections: recorded `Error` diagnostics, mesh untouched ──
    for bad in ("AdaptMesh {2} {} {{}} {2, 0};\n",      # unsupported technique
                "AdaptMesh {1} {} {{}} {2};\n",         # malformed `{niter, meshAll}`
                "AdaptMesh {1} {99} {{}} {2, 0};\n")    # undefined Field tag
        badr=Tessella.execute_geo(_adapt_write_geo(bad))
        @test !isempty(badr.errors)
        @test any(contains(e,"AdaptMesh") for e in badr.errors)
        # failure atomicity: the generated mesh is the pre-AdaptMesh product
        @test badr.mesh isa Tessella.Mesh
        @test size(badr.mesh.tets,2)==size(carrier.tets,2)
        @test badr.mesh.coords==carrier.coords
    end

    # ── CLI -4: session path adapts the loaded/generated mesh ──
    src=_adapt_write_geo("")
    out=joinpath(tempdir(),"tessella_adapt4_"*string(hash(src),base=16)*".msh")
    isfile(out) && rm(out)
    dest=Tessella.CLI.main([src,"-4","-o",out])
    @test dest==out && isfile(out)
    written=Tessella.read_mixed_msh(out)
    @test Tessella.validate(written).ok
    @test any(b->b.msh==4,written.blocks)
end
