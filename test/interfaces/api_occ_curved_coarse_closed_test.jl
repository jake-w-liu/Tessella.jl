module APIOCCCurvedCoarseClosedChecks
using Test, Tessella
const A=Tessella.API

# Sizes larger than the model on closed curved faces. The sized fill queued a
# cell by its Cartesian circumradius, but a parameter cell that runs more than
# half a turn around a sphere has a small flat world image: once a few
# vertices existed next to one seam column, the cells covering the whole far
# side counted as resolved. The result was a lens (sphere h=1: 12 nodes,
# enclosed volume 0.16 of 4.19) or, after the weld dropped duplicate cells,
# an open surface (h=2: 4 triangles, `generate(3)` failed with "triangle edge
# (1, 2) has incidence 1"). Gmsh 4.15.2 returns closed, roughly spherical
# meshes: 27, 17 and 16 nodes with volume 3.23, 2.85, 2.75 at h=1, 2, 10.

# Node count, triangle count, enclosed volume and closedness of the dim-2 mesh.
function shell(build,h;dimension=2)
    A.initialize()
    try
        build();A.model.occ.synchronize()
        A.option("Mesh.MeshSizeMin",0.0)
        A.option("Mesh.MeshSizeMax",h);A.option("Mesh.MeshSizeMin",h)
        A.mesh.generate(dimension)
        tags,xyz,_=A.mesh.get_nodes()
        row=Dict(tag=>index for (index,tag) in enumerate(tags))
        point(tag)=(xyz[3row[tag]-2],xyz[3row[tag]-1],xyz[3row[tag]])
        directed=Dict{Tuple{UInt64,UInt64},Int}()
        used=Set{UInt64}()
        volume=0.0;cells=0
        for (_,surface) in A.model.get_entities(2)
            types,_,nodes=A.mesh.get_elements(2,surface)
            for (type,connectivity) in zip(types,nodes)
                type==2 || continue
                for k in 1:3:length(connectivity)
                    a,b,c=connectivity[k],connectivity[k+1],connectivity[k+2]
                    pa,pb,pc=point(a),point(b),point(c)
                    volume+=(pa[1]*(pb[2]*pc[3]-pb[3]*pc[2])-
                             pa[2]*(pb[1]*pc[3]-pb[3]*pc[1])+
                             pa[3]*(pb[1]*pc[2]-pb[2]*pc[1]))/6
                    for edge in ((a,b),(b,c),(c,a))
                        directed[edge]=get(directed,edge,0)+1
                    end
                    push!(used,a,b,c)
                    cells+=1
                end
            end
        end
        closed=all(==(1),values(directed)) &&
            all(edge->haskey(directed,(edge[2],edge[1])),keys(directed))
        tets=0
        if dimension==3
            _,tet_tags,_=A.mesh.get_elements(3)
            tets=sum(length,tet_tags;init=0)
        end
        return (nodes=length(used),cells=cells,volume=volume,closed=closed,
                tets=tets)
    finally
        A.finalize()
    end
end

@testset "coarse sizes keep closed curved faces closed and round" begin
    sphere=()->A.model.occ.add_sphere(0,0,0,1)
    @testset "sphere h=1 (was a lens)" begin
        s=shell(sphere,1.0)
        @test s.closed
        @test s.nodes>=12                 # Gmsh: 27
        @test s.cells==2s.nodes-4         # Euler: a sphere
        @test 2.2<s.volume<4π/3           # Gmsh: 3.23; was 0.16
    end
    @testset "sphere h=2 (was open) feeds generate(3)" begin
        s=shell(sphere,2.0;dimension=3)
        @test s.closed
        @test s.cells==2s.nodes-4
        @test 2.2<s.volume<4π/3           # Gmsh: 2.85
        @test s.tets>0
    end
    @testset "torus h=2" begin
        s=shell(()->A.model.occ.add_torus(0,0,0,2,0.5),2.0)
        @test s.closed
        @test s.cells==2s.nodes           # Euler: a torus
        @test 2.5<s.volume<2π^2*2*0.25
    end
    @testset "cone with apex h=2" begin
        s=shell(()->A.model.occ.add_cone(0,0,0,0,0,1,1,0),2.0)
        @test s.closed
        @test s.cells==2s.nodes-4
        @test 0.8<s.volume<π/3            # Gmsh: 0.912 (7-gon pyramid)
    end
end
end
