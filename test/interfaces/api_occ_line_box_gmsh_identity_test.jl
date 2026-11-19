module APIOCCLineBoxGmshIdentityChecks
using Test, Tessella
const A=Tessella.API

# Planar faces bounded by straight OpenCASCADE edges are meshed by the port of
# Gmsh's own 1-D and 2-D algorithms, so the result can be Gmsh's result to the
# last bit. Two things kept it from being so:
#
# * a straight OCC edge was graded by the generic curve kernel (numerical
#   tangent, convex-combination inversion): its nodes sat up to 2e-12 off
#   Gmsh's, which is enough to change the frontal mesh of both faces it bounds;
# * the faces of a `Box` (the GEO OpenCASCADE `Box` and `model.add_box`) took a
#   separate refinement arm (1616 triangles on the unit box at size 0.1 where
#   Gmsh has 1470), and their boundary was listed one edge later than the wire
#   OpenCASCADE builds, which changes the initial triangulation.
#
# Every number below is Gmsh 4.15.2 (`gmsh.model.occ.addRectangle` /
# `addBox`, `Mesh.MeshSizeMin = Mesh.MeshSizeMax`, `generate(2)`); with these
# inputs the complete node and element arrays are identical, the test pins a
# sample of them.
@testset "straight OCC edges and Box faces reproduce Gmsh 4.15.2" begin
    function generate(build,h)
        A.initialize()
        try
            build()
            A.option("Mesh.MeshSizeMin",0.0)
            A.option("Mesh.MeshSizeMax",h);A.option("Mesh.MeshSizeMin",h)
            A.mesh.generate(2)
        catch
            A.finalize()
            rethrow()
        end
    end
    triangles(surface)=begin
        types,tags,nodes=A.mesh.get_elements(2,surface)
        @test types==[2]
        (tags[1],nodes[1])
    end

    @testset "unit OCC rectangle, size 0.1" begin
        generate(()->(A.model.occ.add_rectangle(0,0,0,1,1);
                      A.model.occ.synchronize()),0.1)
        try
            tags,xyz,_=A.mesh.get_nodes()
            @test length(tags)==144
            # Edge 1 leaves the origin along +x: exact tenths. Edge 4 comes
            # down the y axis from (0,1,0): 1 - t in floating point.
            _,first_edge,_=A.mesh.get_nodes(1,1)
            @test first_edge[1:3:end]==[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
            _,last_edge,_=A.mesh.get_nodes(1,4)
            @test last_edge[2:3:end]==[0.9,0.8,0.7,0.6,0.5,0.4,
                0.30000000000000004,0.19999999999999996,0.09999999999999998]
            face_tags,face_xyz,_=A.mesh.get_nodes(2,1)
            @test face_tags[1:3]==[41,42,43]
            @test face_xyz[1:9]==[0.45033512242488166,0.08552481249556187,0.0,
                                  0.09500008811816107,0.5360897063089831,0.0,
                                  0.9154700538379252,0.5499999999999999,0.0]
            @test length(triangles(1)[1])==246
        finally
            A.finalize()
        end
    end

    @testset "unit Box, size 0.25" begin
        generate(()->A.model.add_box(0,0,0,1,1,1),0.25)
        try
            tags,xyz,_=A.mesh.get_nodes()
            @test length(tags)==134
            point(tag)=begin
                row=findfirst(==(tag),tags)
                (xyz[3row-2],xyz[3row-1],xyz[3row])
            end
            # First face-owned node of faces 1, 2 and 3.
            @test point(45)==(0.0,0.625,0.21650635094610965)
            @test point(60)==(1.0,0.7890759238962054,0.61533121635129673)
            @test point(75)==(0.21169301126395795,0.0,0.62394065767978091)
            for surface in 1:6
                @test length(triangles(surface)[1])==44
            end
            elements,nodes=triangles(1)
            @test elements[1:2]==[57,58]
            @test nodes[1:6]==[45,49,51,49,45,50]
            elements,nodes=triangles(2)
            @test elements[1:2]==[101,102]
            @test nodes[1:6]==[25,65,26,65,71,26]
        finally
            A.finalize()
        end
    end
end
end
