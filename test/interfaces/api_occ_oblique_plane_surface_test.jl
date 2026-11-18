module APIOCCObliquePlaneSurfaceChecks
using Test, Tessella
const A=Tessella.API
const M=Tessella.Model

# The planar-surface gate certified coplanarity with the exact orient3
# predicate. Every plane that is not axis aligned — a rotated Box, the caps of
# a cylinder or cone with an oblique axis, the sector faces of a partial
# sphere — carries a few ulps of rounding on its vertices, so `generate(2)`
# failed with "Point[2] is not coplanar with Surface[3]" on all of them
# (Gmsh 4.15.2 meshes each one). On a cap bounded by a single circle the
# plane basis was additionally built from the edge vertex and an evaluated
# rim sample of the same point, i.e. from pure rounding noise.

function shell(build;h=0.3)
    A.initialize()
    try
        build();A.model.occ.synchronize()
        A.option("Mesh.MeshSizeMin",h);A.option("Mesh.MeshSizeMax",h)
        A.mesh.generate(2)
        tags,xyz,_=A.mesh.get_nodes()
        row=Dict(tag=>index for (index,tag) in enumerate(tags))
        point(tag)=(xyz[3row[tag]-2],xyz[3row[tag]-1],xyz[3row[tag]])
        undirected=Dict{Tuple{UInt64,UInt64},Int}()
        area=0.0
        for (_,surface) in A.model.get_entities(2)
            types,_,nodes=A.mesh.get_elements(2,surface)
            for (type,connectivity) in zip(types,nodes)
                type==2 || continue
                for k in 1:3:length(connectivity)
                    a,b,c=connectivity[k],connectivity[k+1],connectivity[k+2]
                    u=point(b).-point(a);v=point(c).-point(a)
                    area+=0.5*sqrt((u[2]*v[3]-u[3]*v[2])^2+
                                   (u[3]*v[1]-u[1]*v[3])^2+
                                   (u[1]*v[2]-u[2]*v[1])^2)
                    for edge in ((a,b),(b,c),(c,a))
                        key=minmax(edge...)
                        undirected[key]=get(undirected,key,0)+1
                    end
                end
            end
        end
        return (area=area,closed=all(==(2),values(undirected)),
                surfaces=length(A.model.get_entities(2)))
    finally
        A.finalize()
    end
end

@testset "planes carrying rounding noise pass the planar gate" begin
    @testset "rotated Box" begin
        s=shell() do
            box=A.model.occ.add_box(0,0,0,1,1,1)
            A.model.occ.rotate([(3,box)],0,0,0,1,1,1,0.3)
        end
        @test s.surfaces==6
        @test s.closed
        @test s.area≈6.0 atol=1e-9
    end
    @testset "cylinder with an oblique axis" begin
        s=shell(()->A.model.occ.add_cylinder(0,0,0,1,2,3,0.5))
        @test s.surfaces==3
        @test s.closed
        # Inscribed polyhedron: slightly below 2πr·|axis| + 2πr².
        exact=2pi*0.5*sqrt(14.0)+2pi*0.25
        @test 0.95exact<s.area<exact
    end
    @testset "sphere sector with plane side faces" begin
        s=shell(()->A.model.occ.add_sphere(0,0,0,1;angle3=pi/2))
        @test s.closed
        exact=pi+pi           # quarter sphere + two half disks
        @test 0.95exact<s.area<exact
    end
    @testset "noise is admitted, warp is still rejected" begin
        # The oblique cylinder's cap: its circle vertex, an evaluated rim
        # sample of the same point two ulps away, and two more rim samples.
        cap=NTuple{3,Float64}[
            (1.237053211105853e-16,0.4160251471689219,-0.2773500981126146),
            (-1.7225868909353747e-18,0.41602514716892186,-0.2773500981126146),
            (-0.4172614801981402,-0.1438184997078239,0.23496615987126265),
            (0.4172614801981401,-0.2722066474610981,0.04238393824135202)]
        anchor,second,third,_=M._model_surface_projection(cap,[2,2,2,2],3,"test")
        # The basis skips the near-duplicate of the anchor.
        @test maximum(abs.(second.-anchor))>0.1
        @test maximum(abs.(third.-anchor))>0.1
        # A corner lifted 1e-6 out of a unit square is genuine warp.
        warped=NTuple{3,Float64}[(0.0,0.0,0.0),(1.0,0.0,0.0),(1.0,1.0,1e-6),
                                 (0.0,1.0,0.0)]
        @test_throws ArgumentError M._model_surface_projection(
            warped,[1,2,3,4],1,"test")
        # ... and so is one lifted 1e-9 out of it.
        warped[3]=(1.0,1.0,1e-9)
        @test_throws ArgumentError M._model_surface_projection(
            warped,[1,2,3,4],1,"test")
        # The lifted corner may be the one that spans the widest triangle
        # with the anchor (here with the diagonal corner): it lies on that
        # triangle's plane by construction, so every other corner has to be
        # measured against it too. This wire is a ruled patch, not a plane.
        lifted=NTuple{3,Float64}[(0.0,0.0,0.0),(1.0,0.0,0.0),(1.0,1.0,0.0),
                                 (0.0,1.0,0.5)]
        @test_throws ArgumentError M._model_surface_projection(
            lifted,[1,2,3,4],1,"test")
        # Exactly coplanar input keeps its first distinct point as basis.
        square=NTuple{3,Float64}[(0.0,0.0,0.5),(1.0,0.0,0.5),(1.0,1.0,0.5),
                                 (0.0,1.0,0.5)]
        a,b,c,axes=M._model_surface_projection(square,[1,2,3,4],1,"test")
        @test (a,b)==(square[1],square[2]) && axes==(1,2)
    end
end
end
