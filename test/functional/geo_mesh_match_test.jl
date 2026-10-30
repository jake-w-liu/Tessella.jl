# `Geometry.MatchGeomAndMesh` — `GeomMeshMatcher::match` parity: an opened
# `.msh` grafts onto the current model's entities by vertex position,
# matched-endpoint curves, boundary-set faces and complete boundary regions.
# Physical records follow upstream's direction rule (faces/regions take the
# file's memberships, points/curves keep the model's) and a failed match
# leaves the model and mesh untouched.

function _match_box_geo(lc=1.5)
    return """
lc = $lc;
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
"""
end

_shifted_geo = """
lc = 1.5;
Point(1) = {10,0,0,lc}; Point(2) = {11,0,0,lc}; Point(3) = {11,1,0,lc}; Point(4) = {10,1,0,lc};
Point(5) = {10,0,1,lc}; Point(6) = {11,0,1,lc}; Point(7) = {11,1,1,lc}; Point(8) = {10,1,1,lc};
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
"""

function _match_generate_msh(dir,name,geo,extra="")
    path=joinpath(dir,name)
    open(path,"w") do io
        print(io,geo);println(io,"Mesh 3;");println(io,"Save \"$(joinpath(dir,name*".msh"))\";")
        print(io,extra)
    end
    r=Tessella.execute_geo(path)
    @test isempty(r.errors)
    return joinpath(dir,name*".msh")
end

_dim_tets(mesh)=sum(size(b.nodes,2) for b in mesh.blocks if
    Tessella.Elements.msh_family(b.msh)===:tet;init=0)

@testset "geometry/mesh matching" begin
    mktempdir() do dir
        box_msh=_match_generate_msh(dir,"box",_match_box_geo(1.5))

        # ── graft onto an unmeshed model — every entity re-keys to its
        # model tag, physicals follow upstream's direction rule ──
        geo=joinpath(dir,"match.geo")
        open(geo,"w") do io
            print(io,_match_box_geo(1.5))
            println(io,"Geometry.MatchGeomAndMesh = 1;")
            println(io,"Merge \"$box_msh\";")
        end
        r=Tessella.execute_geo(geo)
        @test isempty(r.errors)
        @test any(w->occursin("matched mesh to CAD",w),r.warnings)
        auth=r.model.meshing.mesh_authority
        @test auth!==nothing && auth.mesh isa Tessella.MixedMesh
        data=auth.mesh.entity_data
        owners=keys(data.entities)
        for key in ((0,i) for i in 1:8);@test key in owners;end
        for key in ((1,i) for i in 1:12);@test key in owners;end
        for key in ((2,i) for i in 1:6);@test key in owners;end
        @test (3,1) in owners
        # Direction rule: the surface/volume take the file's physicals.
        @test data.entities[(2,1)].physical_tags==Int32[10]
        @test data.entities[(3,1)].physical_tags==Int32[20]
        @test isempty(data.entities[(1,1)].physical_tags)
        # Model-side records adopted the file's memberships on dim>=2.
        @test get(r.model.entity_physicals,(3,1),nothing)==[20]
        @test get(r.model.physical,(2,10),nothing)==[1]
        @test get(r.model.physical,(3,20),nothing)==[1]
        @test Tessella.validate(auth.mesh).ok
        @test auth.carrier isa Tessella.Mesh
        @test _dim_tets(auth.mesh)>0
        # Ownership is per-entity: every tet lives on Volume 1, every seg on
        # a Line, and point cells classify to the eight Points.
        for (bi,block) in enumerate(auth.mesh.blocks)
            fam=Tessella.Elements.msh_family(block.msh)
            fam===:tet && @test data.block_entities[bi]==fill(Int32(1),
                size(block.nodes,2))
        end

        # ── replace path — a meshed model adopts the new mesh; the previous
        # cells on matched entities are dropped (upstream deleteMesh) ──
        fine_msh=_match_generate_msh(dir,"fine",_match_box_geo(0.7))
        geo2=joinpath(dir,"match_replace.geo")
        open(geo2,"w") do io
            print(io,_match_box_geo(1.5))
            println(io,"Mesh 3;")
            println(io,"Geometry.MatchGeomAndMesh = 1;")
            println(io,"Merge \"$fine_msh\";")
        end
        r2=Tessella.execute_geo(geo2)
        @test isempty(r2.errors)
        auth2=r2.model.meshing.mesh_authority
        # The adopted mesh is the fine file's, not an accumulation of both.
        @test _dim_tets(auth2.mesh)==
            _dim_tets(Tessella.Elements.read_mixed_msh(fine_msh))
        @test Tessella.validate(auth2.mesh).ok

        # ── failure atomicity — a translated mesh matches nothing; the
        # error is recorded and the previously generated mesh survives ──
        shifted_msh=_match_generate_msh(dir,"shifted",_shifted_geo)
        geo3=joinpath(dir,"match_bad.geo")
        open(geo3,"w") do io
            print(io,_match_box_geo(1.5))
            println(io,"Mesh 3;")
            println(io,"Geometry.MatchGeomAndMesh = 1;")
            println(io,"Merge \"$shifted_msh\";")
        end
        before=Tessella.execute_geo(joinpath(dir,"baseline.geo") |> p ->
            (open(p,"w") do io
                print(io,_match_box_geo(1.5));println(io,"Mesh 3;");p
            end))
        r3=Tessella.execute_geo(geo3)
        @test !isempty(r3.errors)
        @test any(e->occursin("match",e),r3.errors)
        auth3=r3.model.meshing.mesh_authority
        @test auth3!==nothing
        @test size(auth3.mesh.coords,2)==
            size(before.model.meshing.mesh_authority.mesh.coords,2)

        # ── API surface — `Geometry.MatchGeomAndMesh` is a stored session
        # option; `API.merge` grafts through the same path ──
        Tessella.API.initialize()
        try
            plain_geo=joinpath(dir,"plain.geo");write(plain_geo,_match_box_geo(1.5))
            Tessella.API.open(plain_geo)
            Tessella.API.option("Geometry.MatchGeomAndMesh",1)
            Tessella.API.merge(box_msh)
            api_auth=lock(Tessella.API.STATE_LOCK) do
                Tessella.API._model_locked().meshing.mesh_authority
            end
            @test api_auth!==nothing && api_auth.mesh isa Tessella.MixedMesh
            @test (3,1) in keys(api_auth.mesh.entity_data.entities)
            @test Tessella.validate(api_auth.mesh).ok
            # A failing API merge leaves the published model untouched.
            Tessella.API.model.add("again")
            Tessella.API.open(plain_geo)
            Tessella.API.mesh.generate(3)
            @test_throws Exception Tessella.API.merge(shifted_msh)
            still=lock(Tessella.API.STATE_LOCK) do
                Tessella.API._model_locked().meshing.mesh_authority
            end
            @test still!==nothing && _dim_tets(still.mesh)>0
        finally
            Tessella.API.finalize()
        end
    end
end
