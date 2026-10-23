# CLI `-barycentric_refine` — verified against the Gmsh 4.15.2 oracle:
# only linear simplices on entities at the model dimension subdivide —
# triangles split (v3,v1,c),(v1,v2,c),(v2,v3,c) and tetrahedra split
# (v1,v2,v3,c),(v2,v3,v4,c),(v3,v4,v1,c),(v4,v1,v2,c) with one fresh
# barycenter node per parent. Boundary/lower-dimensional cells, quads and
# higher-order cells pass through; children inherit the parent physical tag
# and element tags renumber sequentially in emission order.
module CliBarycentricRefineChecks
using Test, Tessella
using Tessella.CLI: main
using Tessella.Elements: read_mixed_msh, validate, mixed_crc, barycentric_refine,
    ElementBlock, MixedMesh
using Tessella.MeshTypes: Mesh, ntets, ntris

const _SQUARE_GEO="""
Point(1) = {0, 0, 0, 0.3};
Point(2) = {1, 0, 0, 0.3};
Point(3) = {1, 1, 0, 0.3};
Point(4) = {0, 1, 0, 0.3};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Physical Surface("domain", 7) = {1};
"""

const _BOX_GEO="""
Point(1) = {0,0,0,1}; Point(2) = {1,0,0,1}; Point(3) = {1,1,0,1}; Point(4) = {0,1,0,1};
Point(5) = {0,0,1,1}; Point(6) = {1,0,1,1}; Point(7) = {1,1,1,1}; Point(8) = {0,1,1,1};
Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
Line(5) = {5,6}; Line(6) = {6,7}; Line(7) = {7,8}; Line(8) = {8,5};
Line(9) = {1,5}; Line(10) = {2,6}; Line(11) = {3,7}; Line(12) = {4,8};
Curve Loop(1) = {1,2,3,4}; Curve Loop(2) = {5,6,7,8};
Curve Loop(3) = {1,10,-5,-9}; Curve Loop(4) = {2,11,-6,-10};
Curve Loop(5) = {3,12,-7,-11}; Curve Loop(6) = {4,9,-8,-12};
Plane Surface(1) = {1}; Plane Surface(2) = {2}; Plane Surface(3) = {3};
Plane Surface(4) = {4}; Plane Surface(5) = {5}; Plane Surface(6) = {6};
Surface Loop(1) = {1,2,3,4,5,6};
Volume(1) = {1};
Physical Volume("cube", 3) = {1};
"""

function _cells(mesh::MixedMesh,msh::Int)
    out=Tuple{Int,Vector{Int32}}[]
    for block in mesh.blocks
        block isa ElementBlock && block.msh==msh || continue
        for j in axes(block.nodes,2)
            push!(out,(Int(block.tags[j]),collect(block.nodes[:,j])))
        end
    end
    return out
end
_ncells(mesh::MixedMesh,msh::Int)=
    sum(block isa ElementBlock && block.msh==msh ? size(block.nodes,2) : 0
        for block in mesh.blocks)

@testset "CLI -barycentric_refine: 2-D simplex subdivision" begin
    mktempdir() do dir
        geo=joinpath(dir,"sq.geo");write(geo,_SQUARE_GEO)
        base=joinpath(dir,"base.msh");out=joinpath(dir,"ref.msh")
        main([geo,"-2","-o",base])
        before=read_mixed_msh(base)
        ntri=_ncells(before,2);ntri>0 ||
            throw(ErrorException("fixture produced no triangles"))
        nseg=_ncells(before,1);npnt=_ncells(before,15)
        # `-save_all` keeps non-physical boundary cells in the output —
        # upstream's default `Mesh.SaveAll=0` emits physical-entity cells
        # only (verified against Gmsh 4.15.2 `Save`).
        main([base,"-save_all","-barycentric_refine","-o",out])
        refined=read_mixed_msh(out)
        @test validate(refined).ok
        # Every parent triangle emits three centroid children; lower cells
        # pass through unchanged (oracle: lines and points are inert).
        @test _ncells(refined,2)==3*ntri
        @test _ncells(refined,1)==nseg
        @test _ncells(refined,15)==npnt
        @test size(refined.coords,2)==size(before.coords,2)+ntri
        # Child connectivity follows the oracle's (v3,v1),(v1,v2),(v2,v3)
        # cyclic order around the appended centroid.
        before_cells=_cells(before,2)
        refined_cells=_cells(refined,2)
        for (i,(phys,v)) in enumerate(before_cells)
            c=size(before.coords,2)+i
            expected=((v[3],v[1],c),(v[1],v[2],c),(v[2],v[3],c))
            for (k,conn) in enumerate(expected)
                _,got=refined_cells[(i-1)*3+k]
                @test collect(conn)==got
            end
            # The appended centroid is the barycenter of the parent nodes.
            cx=(before.coords[1,v[1]]+before.coords[1,v[2]]+before.coords[1,v[3]])/3
            cy=(before.coords[2,v[1]]+before.coords[2,v[2]]+before.coords[2,v[3]])/3
            @test refined.coords[1,c]≈cx
            @test refined.coords[2,c]≈cy
        end
        # Physical tags are inherited; every refined cell keeps the parent's tag.
        for (i,(phys,v)) in enumerate(before_cells)
            for k in 1:3
                @test refined_cells[(i-1)*3+k][1]==phys
            end
        end
    end
end

@testset "CLI -barycentric_refine: volume tets, inert boundary" begin
    mktempdir() do dir
        geo=joinpath(dir,"box.geo");write(geo,_BOX_GEO)
        base=joinpath(dir,"base.msh");out=joinpath(dir,"ref.msh")
        main([geo,"-3","-o",base])
        before=read_mixed_msh(base)
        ntet=_ncells(before,4);ntet>0 ||
            throw(ErrorException("fixture produced no tetrahedra"))
        ntri=_ncells(before,2);nseg=_ncells(before,1);npnt=_ncells(before,15)
        main([base,"-save_all","-barycentric_refine","-o",out])
        refined=read_mixed_msh(out)
        @test validate(refined).ok
        # Oracle: tets subdivide 1→4; the surface triangles stay untouched.
        @test _ncells(refined,4)==4*ntet
        @test _ncells(refined,2)==ntri
        @test _ncells(refined,1)==nseg
        @test _ncells(refined,15)==npnt
        @test size(refined.coords,2)==size(before.coords,2)+ntet
        before_cells=_cells(before,4)
        refined_cells=_cells(refined,4)
        for (i,(phys,v)) in enumerate(before_cells)
            c=size(before.coords,2)+i
            expected=((v[1],v[2],v[3],c),(v[2],v[3],v[4],c),
                      (v[3],v[4],v[1],c),(v[4],v[1],v[2],c))
            for (k,conn) in enumerate(expected)
                @test collect(conn)==refined_cells[(i-1)*4+k][2]
            end
            for k in 1:4
                @test refined_cells[(i-1)*4+k][1]==phys
            end
            cx=(before.coords[1,v[1]]+before.coords[1,v[2]]+
                before.coords[1,v[3]]+before.coords[1,v[4]])/4
            @test refined.coords[1,c]≈cx
        end
    end
end

@testset "CLI -barycentric_refine: batch-action contract" begin
    mktempdir() do dir
        geo=joinpath(dir,"sq.geo");write(geo,_SQUARE_GEO)
        base=joinpath(dir,"base.msh")
        main([geo,"-2","-o",base])
        err=try
            main([geo,"-2","-barycentric_refine","-o",joinpath(dir,"x.msh")]);""
        catch e;sprint(showerror,e);end
        @test occursin("conflicting batch action",err)
        err=try
            main([base,"-barycentric_refine","-barycentric_refine"]); ""
        catch e;sprint(showerror,e);end
        @test occursin("duplicate -barycentric_refine",err)
    end
end

@testset "barycentric_refine: flat simplex mesh" begin
    coords=reshape([0.0 0 0; 1 0 0; 0 1 0; 0 0 1]',3,4)
    mesh=Mesh(coords;tets=Int32[1 2 3 4]',tet_tag=Int32[5])
    refined=barycentric_refine(mesh)
    @test ntets(refined)==4
    @test size(refined.coords,2)==5
    @test all(refined.tet_tag.==5)
    # Boundary tris are inert in a volume mesh — a tri+quad only mesh is
    # treated as a 2-D model and its tris split while quads pass through.
    flat=Mesh(coords[:,1:3];tris=Int32[1 2 3]',tri_tag=Int32[9])
    flat_refined=barycentric_refine(flat)
    @test ntris(flat_refined)==3
    @test size(flat_refined.coords,2)==4
end
end
