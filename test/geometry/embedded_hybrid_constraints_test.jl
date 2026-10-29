# EMB01: embedded constraints in the native hybrid (quad-bounded) route —
# warped ruled sheets, nonlinear-boundary disks, and nonconvex shells. Every
# emitted sheet facet must be an interior cell face (multiplicity 2), the
# embedded Point/Curve coordinates must be real mesh vertices whose segments
# are cell edges, and an entity outside a nonconvex shell rejects loudly.
# Run: julia --project test/geometry/embedded_hybrid_constraints_test.jl
module EmbeddedHybridConstraintChecks
using Test, LinearAlgebra, Tessella
using Tessella.Elements: MixedMesh, msh_dimension

const CUBE = """
Point(1)={0,0,0,0.4}; Point(2)={1,0,0,0.4}; Point(3)={1,1,0,0.4}; Point(4)={0,1,0,0.4};
Point(5)={0,0,1,0.4}; Point(6)={1,0,1,0.4}; Point(7)={1,1,1,0.4}; Point(8)={0,1,1,0.4};
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
Recombine Surface{1,2,3,4,5,6};
"""

# A genuinely warped ruled patch: the z=0.4/z=0.6 ends lift the sheet off any
# plane, so the embedded disk is certified in its declared (u,v) domain.
const WARPED = """
Point(101)={0.2,0.2,0.4,0.3}; Point(102)={0.8,0.2,0.4,0.3};
Point(103)={0.8,0.8,0.6,0.3}; Point(104)={0.2,0.8,0.6,0.3};
Line(101)={101,102}; Line(102)={102,103}; Line(103)={103,104}; Line(104)={104,101};
Line Loop(101)={101,102,103,104};
Surface(101)={101};
Surface{101} In Volume{1};
"""

# A planar disk bounded by three nonlinear Circles — the sheet is flat but its
# declared boundary is not a Line chain.
const ARCDISK = """
Point(101)={0.75,0.5,0.5,0.25}; Point(102)={0.375,0.7165063509461097,0.5,0.25};
Point(103)={0.375,0.28349364905389027,0.5,0.25}; Point(104)={0.5,0.5,0.5,0.25};
Circle(101)={101,104,102}; Circle(102)={102,104,103}; Circle(103)={103,104,101};
Line Loop(101)={101,102,103};
Plane Surface(101)={101};
Surface{101} In Volume{1};
"""

# L-shaped (nonconvex) all-quad shell.
const LSHAPE = """
Point(1)={0,0,0,0.5}; Point(2)={2,0,0,0.5}; Point(3)={2,1,0,0.5};
Point(4)={1,1,0,0.5}; Point(5)={1,2,0,0.5}; Point(6)={0,2,0,0.5};
Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,5}; Line(5)={5,6}; Line(6)={6,1};
Point(7)={0,0,1,0.5}; Point(8)={2,0,1,0.5}; Point(9)={2,1,1,0.5};
Point(10)={1,1,1,0.5}; Point(11)={1,2,1,0.5}; Point(12)={0,2,1,0.5};
Line(7)={7,8}; Line(8)={8,9}; Line(9)={9,10}; Line(10)={10,11}; Line(11)={11,12}; Line(12)={12,7};
Line(13)={1,7}; Line(14)={2,8}; Line(15)={3,9}; Line(16)={4,10}; Line(17)={5,11}; Line(18)={6,12};
Line Loop(1)={1,2,3,4,5,6}; Line Loop(2)={7,8,9,10,11,12};
Line Loop(3)={1,14,-7,-13}; Line Loop(4)={2,15,-8,-14};
Line Loop(5)={3,16,-9,-15}; Line Loop(6)={4,17,-10,-16};
Line Loop(7)={5,18,-11,-17}; Line Loop(8)={6,13,-12,-18};
Plane Surface(1)={1}; Plane Surface(2)={2};
Plane Surface(3)={3}; Plane Surface(4)={4}; Plane Surface(5)={5};
Plane Surface(6)={6}; Plane Surface(7)={7}; Plane Surface(8)={8};
Surface Loop(1)={1,2,3,4,5,6,7,8};
Volume(1)={1};
Recombine Surface{1,2,3,4,5,6,7,8};
"""

const LPOINT = """
Point(50)={0.5,0.5,0.5,0.3};
Point{50} In Volume{1};
"""

const LCURVE = """
Point(50)={0.5,0.5,0.5,0.3}; Point(51)={1.5,0.5,0.5,0.3};
Line(50)={50,51};
Curve{50} In Volume{1};
"""

# (1.5,1.5,0.5) sits in the missing leg of the L — strictly outside the shell
# but inside its bounding box, so only the exact parity test rejects it.
const LOUTSIDE = """
Point(50)={1.5,1.5,0.5,0.3};
Point{50} In Volume{1};
"""

function tet_faces(mesh)
    faces=Dict{NTuple{3,Int32},Int}();count=0;volume=0.0;positive=true
    blocks=mesh isa MixedMesh ? [b.nodes for b in mesh.blocks if b.msh==4] :
                                [mesh.tets]
    pyramids=mesh isa MixedMesh ? [b.nodes for b in mesh.blocks if b.msh==7] :
                                Matrix{Int32}[]
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
    # Pyramid side faces meet Tet faces; the base splits on the same kept
    # diagonal the kernel fans use, so boundary quads count once per fan.
    for cells in pyramids, cell in eachcol(cells)
        count+=1
        b1,b2,b3,b4,a=(mesh.coords[:,cell[i]] for i in 1:5)
        for (x,y,z) in ((b1,b2,b3),(b1,b3,b4))
            det=LinearAlgebra.det(hcat(y-x,z-x,a-x))
            positive &= det>0; volume+=det/6
        end
        for side in ((1,2,5),(2,3,5),(3,4,5),(4,1,5),(1,2,3),(1,3,4))
            key=Tuple(sort!(Int32[cell[i] for i in side]))
            faces[key]=get(faces,key,0)+1
        end
    end
    return faces,count,volume,positive
end

function tet_edges(mesh)
    edges=Set{NTuple{2,Int32}}()
    blocks=mesh isa MixedMesh ? [b.nodes for b in mesh.blocks if b.msh==4] :
                                [mesh.tets]
    pyramids=mesh isa MixedMesh ? [b.nodes for b in mesh.blocks if b.msh==7] :
                                Matrix{Int32}[]
    for cells in blocks, tet in eachcol(cells)
        for (i,j) in ((1,2),(1,3),(1,4),(2,3),(2,4),(3,4))
            push!(edges,minmax(tet[i],tet[j]))
        end
    end
    for cells in pyramids, cell in eachcol(cells)
        for (i,j) in ((1,2),(2,3),(3,4),(4,1),(1,5),(2,5),(3,5),(4,5))
            push!(edges,minmax(cell[i],cell[j]))
        end
    end
    return edges
end

function sheet_cells(mesh)
    data=mesh.entity_data; out=Dict{Int,Vector{NTuple{3,Int32}}}()
    for (bi,block) in enumerate(mesh.blocks)
        msh_dimension(block.msh)==2 || continue
        for (ci,cell) in enumerate(eachcol(block.nodes))
            ent=Int(data.block_entities[bi][ci])
            ent>=100 || continue
            push!(get!(out,ent,NTuple{3,Int32}[]),Tuple(sort!(collect(cell))))
        end
    end
    return out
end

entity_nodes(mesh,dim,tag)=
    Int[row for (row,owner) in enumerate(mesh.entity_data.node_entities)
        if owner==(dim,Int32(tag))]

function run_volume(name,geo)
    path=joinpath(mktempdir(),name*".geo");write(path,geo)
    result=Tessella.execute_geo(path;mesh_dim=3)
    return result.model.meshing.mesh_authority.mesh
end

@testset "hybrid embedded surfaces and nonconvex shells" begin
    @testset "warped ruled sheet in a hybrid Region" begin
        mesh=run_volume("hybrid_warped",CUBE*WARPED)
        @test mesh isa MixedMesh && mesh.entity_data!==nothing
        faces,tets,volume,positive=tet_faces(mesh)
        @test tets>0 && positive
        @test volume≈1.0 atol=1e-9
        @test all(count->count in (1,2),values(faces))
        cells=sheet_cells(mesh)
        # The literal warped facets emit and every one is interior.
        @test haskey(cells,101) && !isempty(cells[101])
        @test all(face->get(faces,face,0)==2,cells[101])
        @test allunique(cells[101])
    end
    @testset "planar disk with nonlinear Circle boundary in a hybrid Region" begin
        mesh=run_volume("hybrid_arcdisk",CUBE*ARCDISK)
        @test mesh isa MixedMesh && mesh.entity_data!==nothing
        faces,tets,volume,positive=tet_faces(mesh)
        @test tets>0 && positive
        @test volume≈1.0 atol=1e-9
        cells=sheet_cells(mesh)
        @test haskey(cells,101) && !isempty(cells[101])
        @test all(face->get(faces,face,0)==2,cells[101])
    end
    @testset "embedded Point in a nonconvex hybrid shell" begin
        mesh=run_volume("hybrid_nc_point",LSHAPE*LPOINT)
        @test mesh isa MixedMesh && mesh.entity_data!==nothing
        faces,tets,volume,positive=tet_faces(mesh)
        @test tets>0 && positive
        @test volume≈3.0 atol=1e-9
        @test all(count->count in (1,2),values(faces))
        rows=entity_nodes(mesh,0,50)
        @test length(rows)==1
        @test Tuple(mesh.coords[:,only(rows)])==(0.5,0.5,0.5)
    end
    @testset "embedded Curve in a nonconvex hybrid shell" begin
        mesh=run_volume("hybrid_nc_curve",LSHAPE*LCURVE)
        @test mesh isa MixedMesh && mesh.entity_data!==nothing
        faces,tets,volume,positive=tet_faces(mesh)
        @test tets>0 && positive
        @test volume≈3.0 atol=1e-9
        @test all(count->count in (1,2),values(faces))
        # The embedded chain runs from the first declared endpoint to the
        # last; every consecutive actual segment remains a real cell edge.
        rows=sort!(collect(Set(entity_nodes(mesh,0,50) ∪ entity_nodes(mesh,0,51) ∪
                               entity_nodes(mesh,1,50)));by=row->mesh.coords[1,row])
        @test length(rows)>=2
        @test Tuple(mesh.coords[:,first(rows)])==(0.5,0.5,0.5)
        @test Tuple(mesh.coords[:,last(rows)])==(1.5,0.5,0.5)
        edges=tet_edges(mesh)
        @test all(index->minmax(rows[index],rows[index+1]) in edges,
                  1:length(rows)-1)
    end
    @testset "entity in the nonconvex notch rejects loudly" begin
        path=joinpath(mktempdir(),"hybrid_nc_outside.geo");write(path,LSHAPE*LOUTSIDE)
        @test_throws ArgumentError Tessella.execute_geo(path;mesh_dim=3)
    end
end
end # module
