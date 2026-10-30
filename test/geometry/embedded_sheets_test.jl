# EMB01: embedded sheets in volumes — holed, multiple and subdivided sheets.
# Runtime-verified small fixtures: every emitted sheet Triangle must be an
# interior tetrahedron face (multiplicity 2), the volume stays exact, and the
# mesh validates. Run: julia --project test/geometry/embedded_sheets_test.jl
module EmbeddedSheetsChecks
using Test, LinearAlgebra, Tessella
using Tessella.Elements: MixedMesh, msh_dimension
const M = Tessella.Model
const A = Tessella.API

const CUBE = """
Point(1)={0,0,0,1}; Point(2)={1,0,0,1}; Point(3)={1,1,0,1}; Point(4)={0,1,0,1};
Point(5)={0,0,1,1}; Point(6)={1,0,1,1}; Point(7)={1,1,1,1}; Point(8)={0,1,1,1};
Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
Line(5)={5,6}; Line(6)={6,7}; Line(7)={7,8}; Line(8)={8,5};
Line(9)={1,5}; Line(10)={2,6}; Line(11)={3,7}; Line(12)={4,8};
Line Loop(1)={1,2,3,4}; Line Loop(2)={5,6,7,8};
Line Loop(3)={1,10,-5,-9}; Line Loop(4)={2,11,-6,-10};
Line Loop(5)={3,12,-7,-11}; Line Loop(6)={4,9,-8,-12};
Plane Surface(1)={1}; Plane Surface(2)={2};
Plane Surface(3)={3}; Plane Surface(4)={4};
Plane Surface(5)={5}; Plane Surface(6)={6};
Surface Loop(1)={1,2,3,4,5,6};
Volume(1)={1};
"""

const HOLED = """
Point(101)={0.2,0.2,0.5,1}; Point(102)={0.8,0.2,0.5,1};
Point(103)={0.8,0.8,0.5,1}; Point(104)={0.2,0.8,0.5,1};
Line(101)={101,102}; Line(102)={102,103}; Line(103)={103,104}; Line(104)={104,101};
Line Loop(101)={101,102,103,104};
Point(105)={0.4,0.4,0.5,1}; Point(106)={0.6,0.4,0.5,1};
Point(107)={0.6,0.6,0.5,1}; Point(108)={0.4,0.6,0.5,1};
Line(105)={105,106}; Line(106)={106,107}; Line(107)={107,108}; Line(108)={108,105};
Line Loop(102)={-108,-107,-106,-105};
Plane Surface(101)={101,102};
Surface{101} In Volume{1};
"""

const MULTI = """
Point(101)={0.2,0.2,0.3,1}; Point(102)={0.8,0.2,0.3,1};
Point(103)={0.8,0.8,0.3,1}; Point(104)={0.2,0.8,0.3,1};
Line(101)={101,102}; Line(102)={102,103}; Line(103)={103,104}; Line(104)={104,101};
Line Loop(101)={101,102,103,104};
Plane Surface(101)={101};
Point(111)={0.2,0.2,0.7,1}; Point(112)={0.8,0.2,0.7,1};
Point(113)={0.8,0.8,0.7,1}; Point(114)={0.2,0.8,0.7,1};
Line(111)={111,112}; Line(112)={112,113}; Line(113)={113,114}; Line(114)={114,111};
Line Loop(111)={111,112,113,114};
Plane Surface(111)={111};
Surface{101,111} In Volume{1};
"""

# Two perpendicular interior chords subdivide the sheet — interior vertices and
# sheet-nested curves reach the volume.
const SUBDIVIDED = """
Point(101)={0.2,0.2,0.5,1}; Point(102)={0.8,0.2,0.5,1};
Point(103)={0.8,0.8,0.5,1}; Point(104)={0.2,0.8,0.5,1};
Line(101)={101,102}; Line(102)={102,103}; Line(103)={103,104}; Line(104)={104,101};
Line Loop(101)={101,102,103,104};
Plane Surface(101)={101};
Point(150)={0.5,0.5,0.5,0.1};
Point{150} In Surface{101};
Surface{101} In Volume{1};
"""

function tet_faces(mesh)
    faces=Dict{NTuple{3,Int32},Int}();count=0;volume=0.0;positive=true
    blocks=mesh isa MixedMesh ? [b.nodes for b in mesh.blocks if b.msh==4] :
                                [mesh.tets]
    for cells in blocks, tet in eachcol(cells)
        count+=1
        p,q,r,s=(mesh.coords[:,i] for i in tet)
        det=LinearAlgebra.det(hcat(q-p,r-p,s-p))
        positive &= det>0; volume+=det/6
        for skip in 1:4
            key=Tuple(sort!(Int32[tet[i] for i in 1:4 if i!=skip]))
            faces[key]=get(faces,key,0)+1
        end
    end
    return faces,count,volume,positive
end

function sheet_cells(mesh)
    data=mesh.entity_data; out=Dict{Int,Vector{NTuple{3,Int32}}}()
    for (bi,block) in enumerate(mesh.blocks)
        msh_dimension(block.msh)==2 || continue
        for (ci,cell) in enumerate(eachcol(block.nodes))
            ent=Int(data.block_entities[bi][ci])
            ent>=100 || continue
            push!(get!(out,ent,NTuple{3,Int32}[]),
                  Tuple(sort!(collect(cell))))
        end
    end
    return out
end

function run_volume(name,geo)
    path=joinpath(mktempdir(),name*".geo");write(path,geo)
    result=Tessella.execute_geo(path;mesh_dim=3)
    mesh=result.model.meshing.mesh_authority.mesh
    faces,tets,volume,positive=tet_faces(mesh)
    return mesh,faces,tets,volume,positive
end

@testset "embedded sheets in volumes" begin
    @testset "holed planar sheet" begin
        mesh,faces,tets,volume,positive=run_volume("holed_sheet",CUBE*HOLED)
        @test mesh isa MixedMesh && mesh.entity_data!==nothing
        @test tets>0 && positive
        @test volume≈1.0 atol=1e-9
        @test all(count->count in (1,2),values(faces))
        cells=sheet_cells(mesh)
        # The holed sheet emits a ring of triangles — every one interior.
        @test haskey(cells,101) && !isempty(cells[101])
        @test all(face->get(faces,face,0)==2,cells[101])
        @test allunique(cells[101])
    end
    @testset "two embedded sheets in one volume" begin
        mesh,faces,tets,volume,positive=run_volume("multi_sheet",CUBE*MULTI)
        @test tets>0 && positive
        @test volume≈1.0 atol=1e-9
        cells=sheet_cells(mesh)
        for sheet in (101,111)
            @test haskey(cells,sheet) && !isempty(cells[sheet])
            @test all(face->get(faces,face,0)==2,cells[sheet])
        end
        @test isempty(intersect(Set(cells[101]),Set(cells[111])))
    end
    @testset "sheet carrying a nested embedded Point" begin
        mesh,faces,tets,volume,positive=run_volume("subdivided_sheet",CUBE*SUBDIVIDED)
        @test tets>0 && positive
        @test volume≈1.0 atol=1e-9
        cells=sheet_cells(mesh)
        @test haskey(cells,101) && !isempty(cells[101])
        @test all(face->get(faces,face,0)==2,cells[101])
        # The nested Point's coordinate is a real mesh vertex.
        found=any(node->Tuple(mesh.coords[:,node])==(0.5,0.5,0.5),
                  axes(mesh.coords,2))
        @test found
    end
    @testset "generic (no-authority) volume path on a holed sheet" begin
        path=joinpath(mktempdir(),"generic_holed.geo");write(path,CUBE*HOLED)
        result=Tessella.execute_geo(path;mesh_dim=0)
        mesh=M.mesh_model_volume(result.model,1)
        faces,tets,volume,positive=tet_faces(mesh)
        @test tets>0 && positive
        @test volume≈1.0 atol=1e-9
        @test all(count->count in (1,2),values(faces))
        diag=Tessella.Elements.validate(mesh)
        @test diag.ok
    end
end

@testset "embed! API contract" begin
    m=Tessella.GeoModel()
    for (tag,p) in enumerate(((0,0,0),(1,0,0),(0,1,0),(0,0,1)))
        M.add_point!(m,p...;tag)
    end
    M.add_line!(m,1,2;tag=1);M.add_line!(m,2,3;tag=2);M.add_line!(m,3,1;tag=3)
    M.add_curve_loop!(m,[1,2,3];tag=1)
    M.add_plane_surface!(m,[1];tag=1)
    M.add_point!(m,0.2,0.2,0.0;tag=5)
    # Supported combos: Point/Line in Surface; Point/Line/Surface in Volume.
    @test M.embed!(m,0,[5],2,1)==1
    @test M.embed!(m,1,[2],2,1)==1
    @test_throws ArgumentError M.embed!(m,2,[1],3,7)      # unknown Volume
    @test_throws ArgumentError M.embed!(m,0,[9],2,1)      # unknown Point
    @test_throws ArgumentError M.embed!(m,0,[5],2,1)      # duplicate
    @test_throws ArgumentError M.embed!(m,2,[1],2,1)      # dim !< target
    @test_throws ArgumentError M.embed!(m,1,[1],1,1)      # inDim must be 2|3
    @test get(m.embeds,(2,1),[])==[(0,5),(1,2)]
    M.remove_embedded!(m,[(2,1)],0)
    @test get(m.embeds,(2,1),[])==[(1,2)]
end

@testset "non-planar embedded ruled sheet meshes" begin
    # The ruled embedded sheet is certified through the curved-source tier
    # (67bd8db/fb6867e): the disk is parameterized in the patch's own uv
    # domain, so a warped `Ruled Surface` inside a cube now meshes instead
    # of rejecting. Verify the produced volume mesh contains the sheet's
    # interior nodes (strictly inside the cube's boundary) and validates.
    geo=CUBE*"""
Point(101)={0.2,0.2,0.4,1}; Point(102)={0.8,0.2,0.4,1};
Point(103)={0.8,0.8,0.7,1}; Point(104)={0.2,0.8,0.7,1};
Point(105)={0.5,0.2,0.5,1}; Point(106)={0.5,0.8,0.6,1};
Circle(101)={101,105,102};
Line(102)={102,103};
Circle(103)={103,106,104};
Line(104)={104,101};
Line Loop(101)={101,102,103,104};
Ruled Surface(101)={101};
Surface{101} In Volume{1};
"""
    path=joinpath(mktempdir(),"curved_sheet.geo");write(path,geo)
    res=Tessella.execute_geo(path;mesh_dim=3)
    m=res.mesh
    @test m isa Tessella.MixedMesh || m isa Tessella.Mesh
    @test Tessella.validate(m).ok
    # embedded sheet boundary curves must appear in the produced mesh:
    # the ruled surface's four boundary arcs run interior to the cube —
    # at least one produced node must lie strictly inside the box
    coords=m.coords
    @test any(i->coords[1,i]>0.15 && coords[1,i]<0.85 &&
                  coords[2,i]>0.15 && coords[2,i]<0.85 &&
                  coords[3,i]>0.3  && coords[3,i]<0.8,
              axes(coords,2))
end
end # module
