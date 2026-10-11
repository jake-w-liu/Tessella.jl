module RegionAttributeDispatchFunctionalTests
using Test, LinearAlgebra, Tessella
const A=Tessella.API

function session(callback)
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.option("Mesh.Optimize",0)
        A.option("Mesh.Smoothing",0)
        callback()
    finally
        A.finalize()
    end
end

function shell(points,edges,loops;signed_faces=collect(eachindex(loops)),recombine=false)
    for (tag,p) in enumerate(points)
        A.model.add_point(p...;tag,meshSize=3.0)
    end
    for (tag,ends) in enumerate(edges)
        A.model.add_line(ends...;tag)
        A.mesh.set_transfinite_curve(tag,2)
    end
    for (tag,lines) in enumerate(loops)
        A.model.add_curve_loop(collect(lines);tag)
        A.model.add_plane_surface([tag];tag)
        A.mesh.set_transfinite_surface(tag)
        recombine && A.mesh.set_recombine(2,tag)
    end
    A.model.add_surface_loop(signed_faces;tag=1)
    A.model.add_volume([1];tag=1)
end

function tetrahedron()
    shell(((0.,0.,0.),(1.,0.,0.),(0.,1.,0.),(0.,0.,1.)),
          ((1,2),(2,3),(3,1),(1,4),(2,4),(3,4)),
          ((1,2,3),(1,5,-4),(2,6,-5),(3,4,-6));signed_faces=[-1,2,3,4])
end

function cube()
    shell(((0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.),
           (0.,0.,1.),(1.,0.,1.),(1.,1.,1.),(0.,1.,1.)),
          ((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),(7,8),(8,5),
           (1,5),(2,6),(3,7),(4,8)),
          ((1,2,3,4),(5,6,7,8),(1,10,-5,-9),
           (2,11,-6,-10),(3,12,-7,-11),(4,9,-8,-12));recombine=true)
end

function vertices()
    tags,coordinates,_=A.mesh.get_nodes()
    return Dict(tag=>coordinates[3i-2:3i] for (i,tag) in enumerate(tags))
end

function signed_tet(points)
    a,b,c,d=points
    return det(hcat(b-a,c-a,d-a))/6
end

function volume_cells()
    points=vertices()
    families,labels,connectivity=A.mesh.get_elements(3,1)
    volumes=Float64[]
    for (family,ids,nodes) in zip(families,labels,connectivity)
        @test family in (4,7)
        width=family==4 ? 4 : 5
        for cell in eachindex(ids)
            p=[points[uid] for uid in nodes[width*(cell-1)+1:width*cell]]
            volume=family==4 ? signed_tet(p) :
                signed_tet(p[[1,2,3,5]])+signed_tet(p[[1,3,4,5]])
            push!(volumes,volume)
        end
    end
    return families,volumes
end

@testset "native attribute dimension dispatch" begin
    session() do
        before=(A.model.get_entities(),A.mesh.get_nodes(),A.mesh.get_elements())
        for dim in (-1,0,1,3,4)
            @test A.mesh.set_smoothing(dim,999,-7)===nothing
        end
        for dim in (-1,0,3,4)
            @test A.mesh.set_reverse(dim,999,true)===nothing
        end
        @test (A.model.get_entities(),A.mesh.get_nodes(),A.mesh.get_elements())==before
        @test_throws ArgumentError A.mesh.set_smoothing(3.0,999,5)
        @test_throws ArgumentError A.mesh.set_reverse(true,999,true)
        @test_throws ArgumentError A.mesh.set_smoothing(big(typemax(Int))+1,999,5)
    end
end

@testset "Volume attributes are public native no-ops on Tet and hybrid meshes" begin
    for (build,expected,families) in ((tetrahedron,1/6,Set([4])),(cube,1.,Set([4,7])))
        session() do
            build()
            @test A.mesh.set_reverse(3,1,true)===nothing
            @test A.mesh.set_smoothing(3,1,5)===nothing
            A.mesh.generate(3)
            actual,volumes=volume_cells()
            @test Set(actual)==families
            @test !isempty(volumes) && all(>(0),volumes)
            @test sum(volumes)≈expected
            @test Tessella.validate(A.mesh.get()).ok
        end
    end
end

@testset "Surface smoothing dispatch and reverse remain supported" begin
    session() do
        for (tag,p) in enumerate(((0.,0.,0.),(1.,0.,0.),(0.,1.,0.)))
            A.model.add_point(p...;tag,meshSize=3.)
        end
        for (tag,ends) in enumerate(((1,2),(2,3),(3,1)))
            A.model.add_line(ends...;tag)
            A.mesh.set_transfinite_curve(tag,2)
        end
        A.model.add_curve_loop([1,2,3];tag=1)
        A.model.add_plane_surface([1];tag=1)
        A.mesh.set_transfinite_surface(1)
        @test A.mesh.set_smoothing(2,1,2)===nothing
        @test_throws ArgumentError A.mesh.set_smoothing(2,999,2)
        @test A.mesh.set_reverse(2,1,true)===nothing
        @test_throws ArgumentError A.mesh.set_reverse(2,999,true)
        A.mesh.generate(2)
        points=vertices()
        families,labels,connectivity=A.mesh.get_elements(2,1)
        @test families==[2] && length(only(labels))==1
        a,b,c=[points[uid] for uid in only(connectivity)]
        @test ((b[1]-a[1])*(c[2]-a[2])-(b[2]-a[2])*(c[1]-a[1]))/2≈-0.5
    end
end
end
