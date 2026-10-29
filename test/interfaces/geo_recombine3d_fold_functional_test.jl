module GeoRecombine3DFoldFunctionalChecks
using Test,LinearAlgebra,Tessella
const E=Tessella.Elements

const _GEO_CUBE=raw"""
Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
Point(3)={1,1,0,1}; Point(4)={0,1,0,1};
Point(5)={0,0,1,1}; Point(6)={1,0,1,1};
Point(7)={1,1,1,1}; Point(8)={0,1,1,1};
Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
Line(5)={5,6}; Line(6)={6,7}; Line(7)={7,8}; Line(8)={8,5};
Line(9)={1,5}; Line(10)={2,6}; Line(11)={3,7}; Line(12)={4,8};
Curve Loop(1)={1,2,3,4}; Curve Loop(2)={-8,-7,-6,-5};
Curve Loop(3)={9,5,-10,-1}; Curve Loop(4)={10,6,-11,-2};
Curve Loop(5)={11,7,-12,-3}; Curve Loop(6)={12,8,-9,-4};
Plane Surface(1)={1}; Plane Surface(2)={2}; Plane Surface(3)={3};
Plane Surface(4)={4}; Plane Surface(5)={5}; Plane Surface(6)={6};
Surface Loop(1)={1,2,3,4,5,6}; Volume(1)={1};
Transfinite Curve{:}=3;
Transfinite Surface{:}; Recombine Surface{:};
"""

function _execute_recombine3d(source::AbstractString;mesh_dim=3)
    return mktemp() do path,io
        write(io,source)
        close(io)
        Tessella.execute_geo(path;mesh_dim=mesh_dim)
    end
end

function _cell_inventory(mesh)
    counts=Dict{Int,Int}()
    if mesh isa Tessella.MixedMesh
        for block in mesh.blocks
            block isa E.ElementBlock || continue
            counts[block.msh]=get(counts,block.msh,0)+size(block.nodes,2)
        end
    else
        size(mesh.segs,2)>0 && (counts[1]=size(mesh.segs,2))
        size(mesh.tris,2)>0 && (counts[2]=size(mesh.tris,2))
        size(mesh.tets,2)>0 && (counts[4]=size(mesh.tets,2))
    end
    return counts
end

function _surface_part(execution,tag)
    part=Tessella.geo_entity_mesh(execution,2,tag)
    part===nothing && return Dict{Int,Int}()
    # Restrict to facet families — part products also carry node/edge blocks.
    return Dict(k=>v for (k,v) in _cell_inventory(
        Tessella.Model.model_to_mixed(execution.model,part,2,tag)) if k in (2,3))
end

# Signed tet volumes sum to the unit cube — an independent certificate of the
# boundary fold (the fill sees only triangles, like upstream's
# `quadsToTriangles` before `MeshDelaunayVolume`).
function _check_tet_volume(mesh)
    tets=NTuple{4,Int32}[]
    if mesh isa Tessella.MixedMesh
        for block in mesh.blocks
            block isa E.ElementBlock && block.msh==4 || continue
            for column in axes(block.nodes,2)
                push!(tets,Tuple(block.nodes[:,column]))
            end
        end
    else
        for column in axes(mesh.tets,2)
            push!(tets,Tuple(mesh.tets[:,column]))
        end
    end
    @test !isempty(tets)
    total=0.0
    for (a,b,c,d) in tets
        pa,pb,pc,pd=(mesh.coords[:,row] for row in (a,b,c,d))
        total+=det(hcat(pb-pa,pc-pa,pd-pa))/6
    end
    @test total≈1.0 atol=1e-12
    return nothing
end

@testset "Recombine Volume folds boundary quads to triangles" begin
    folded=_execute_recombine3d(_GEO_CUBE * raw"""
        Recombine Volume{1};
        Mesh 3;
        """)
    mesh=folded.mesh
    counts=_cell_inventory(mesh)
    @test !haskey(counts,3) && !haskey(counts,7)
    @test counts[2]==48 && get(counts,4,0)>0
    for surface in 1:6
        part=_surface_part(folded,surface)
        @test part==Dict(2=>8)
    end
    @test Tessella.validate(mesh).ok
    _check_tet_volume(mesh)
end

@testset "Mesh.Recombine3DAll selects the same fold" begin
    folded=_execute_recombine3d(_GEO_CUBE * raw"""
        Mesh.Recombine3DAll=1;
        Mesh 3;
        """)
    mesh=folded.mesh
    counts=_cell_inventory(mesh)
    @test !haskey(counts,3) && !haskey(counts,7)
    @test counts[2]==48
    for surface in 1:6
        @test _surface_part(folded,surface)==Dict(2=>8)
    end
    @test Tessella.validate(mesh).ok
    _check_tet_volume(mesh)
end

@testset "retained Mesh 2 products fold at Mesh 3 fill time" begin
    # Upstream `Mesh3D` splits the STORED boundary quadrangles of a
    # `recombine3D` Region before `MeshDelaunayVolume` — faces meshed before
    # `Recombine Volume` took effect still feed the fill an all-triangle PLC
    # and stay split in the published products and constructor receipts.
    retained=_execute_recombine3d(_GEO_CUBE * raw"""
        Mesh 2;
        Recombine Volume{1};
        Mesh 3;
        """)
    mesh=retained.mesh
    counts=_cell_inventory(mesh)
    @test !haskey(counts,3) && !haskey(counts,7)
    @test counts[2]==48 && get(counts,4,0)>0
    for surface in 1:6
        @test _surface_part(retained,surface)==Dict(2=>8)
    end
    @test Tessella.validate(mesh).ok
    _check_tet_volume(mesh)
    # The same fold honors `Mesh.Recombine3DAll` toggled between the two
    # meshing passes.
    late_global=_execute_recombine3d(_GEO_CUBE * raw"""
        Mesh 2;
        Mesh.Recombine3DAll=1;
        Mesh 3;
        """)
    counts=_cell_inventory(late_global.mesh)
    @test !haskey(counts,3) && !haskey(counts,7)
    @test counts[2]==48 && get(counts,4,0)>0
    for surface in 1:6
        @test _surface_part(late_global,surface)==Dict(2=>8)
    end
    @test Tessella.validate(late_global.mesh).ok
    _check_tet_volume(late_global.mesh)
end

@testset "unflagged volumes keep quads and Pyramid transitions" begin
    plain=_execute_recombine3d(_GEO_CUBE * "Mesh 3;\n")
    mesh=plain.mesh
    counts=_cell_inventory(mesh)
    # Quad boundary faces admit Pyramid tops feeding the tetrahedral fill —
    # the hybrid admission the recombine3D fold deliberately bypasses.
    @test counts[3]==24 && get(counts,7,0)>0 && get(counts,4,0)>0
    for surface in 1:6
        @test _surface_part(plain,surface)==Dict(3=>4)
    end
    @test Tessella.validate(mesh).ok
    mesh isa Tessella.MixedMesh &&
        @test Tessella.validate(mesh;certify_jacobians=true).ok
end

@testset "Transfinite and Layers volumes are not Delaunay-admitted" begin
    # `Transfinite Volume` never enters upstream's delaunay list, so its faces
    # keep their quadrangles even under the global recombine3D flag.
    structured=_execute_recombine3d(_GEO_CUBE * raw"""
        Transfinite Volume{1};
        Mesh.Recombine3DAll=1;
        Mesh 3;
        """)
    mesh=structured.mesh
    counts=_cell_inventory(mesh)
    @test get(counts,5,0)==8 && !haskey(counts,4)
    @test counts[3]==24
    for surface in 1:6
        @test _surface_part(structured,surface)==Dict(3=>4)
    end
    @test Tessella.validate(mesh).ok
    # A Layers-extruded volume is likewise outside the delaunay set — its
    # boundary quads are not split by `Recombine Volume` on other entities.
    extruded=_execute_recombine3d(raw"""
        Mesh.Recombine3DAll=1;
        Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
        Point(3)={1,1,0,1}; Point(4)={0,1,0,1};
        Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};
        Plane Surface(1)={1};
        Transfinite Curve{:}=3;
        Transfinite Surface{1}; Recombine Surface{1};
        out[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;};
        """)
    counts=_cell_inventory(extruded.mesh)
    @test !haskey(counts,4) && get(counts,5,0)+get(counts,6,0)>0
end

@testset "session API exposes the stored Recombine3D controls" begin
    A=Tessella.API
    A.initialize()
    try
        @test A.option("Mesh.Recombine3DAll")==0.0
        @test A.option("Mesh.Recombine3DLevel")==0.0
        @test A.option("Mesh.Recombine3DConformity")==0.0
        @test A.option("Mesh.Recombine3DAll",1.0)==1.0
        @test A.option("Mesh.Recombine3DAll")==1.0
        @test A.option("Mesh.Recombine3DAll",0.0)==0.0
        @test A.option("Mesh.Recombine3DLevel",2.9)==2.0
        @test A.option("Mesh.Recombine3DConformity",1.0)==1.0
        # `API.mesh.set_recombine(3, tag)` records `GRegion::recombine3D`; the
        # delaunay-admission predicate folds that Region's boundary quads.
        box=A.model.add_box(0,0,0,1,1,1)
        m=A._model_locked()
        @test !Tessella.Model._model_volume_recombine3d(m,box)
        A.mesh.set_recombine(3,box)
        @test m.meshing.recombine[(3,box)]==45.0
        @test Tessella.Model._model_volume_recombine3d(m,box)
        # `Mesh.Recombine3DAll` selects every Delaunay-admitted Region —
        # transport is exercised end-to-end by the `.geo` testsets above.
        m.meshing.recombine_3d_all=true
        empty!(m.meshing.recombine)
        @test Tessella.Model._model_volume_recombine3d(m,box)
        m.meshing.recombine_3d_all=false
    finally
        A.finalize()
    end
end
end
