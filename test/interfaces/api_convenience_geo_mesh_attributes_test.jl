module ConvenienceGeoMeshAttributesControlsV1
using Tessella, Test
const A=Tessella.API
function curve_laws()
    return lock(A.STATE_LOCK) do
        deepcopy(A.CURRENT[].meshing.transfinite_curves)
    end
end
@testset "Convenience CAD factories preserve live Curve mesh laws" begin
    A.initialize()
    try
        for (tag,(x,y)) in enumerate(((0.,0.),(1.,0.),(0.,1.)))
            A.model.add_point(x,y,0;tag)
        end
        for (tag,(a,b)) in enumerate(((1,2),(2,3),(3,1)))
            prior=curve_laws()
            A.model.add_line(a,b;tag)
            @test curve_laws()==prior
            A.mesh.set_transfinite_curve(tag,5)
        end
        before=curve_laws()
        A.model.add_curve_loop([1,2,3];tag=1)
        @test curve_laws()==before
        A.model.add_plane_surface([1];tag=1)
        @test curve_laws()==before
        A.mesh.set_transfinite_surface(1)
        @test size(A.mesh.generate(2).coords,2)==21
        A.option("Mesh.TransfiniteTri",1)
        @test size(A.mesh.generate(2).coords,2)==15
        A.model.synchronize()
        @test isempty(curve_laws())
        @test Set(A.model.get_entities(1))==Set(((1,1),(1,2),(1,3)))
    finally
        A.finalize()
    end
end
@testset "Explicit GEO synchronization restores internal Curve declaration" begin
    A.initialize()
    try
        A.model.add_point(0.,0.,0.;tag=1)
        A.model.add_point(1.,0.,0.;tag=2)
        A.model.add_line(1,2;tag=11)
        A.model.geo.mesh.set_transfinite_curve(11,7)
        A.mesh.set_transfinite_curve(11,9)
        A.model.add_point(2.,1.,0.;tag=3)
        @test curve_laws()[11].num_nodes==9
        A.model.geo.synchronize()
        @test curve_laws()[11].num_nodes==7
        A.model.set_tag(1,11,91)
        A.model.add_point(3.,1.,0.;tag=4)
        @test Set(A.model.get_entities(1))==Set(((1,91),))
        @test curve_laws()[91].num_nodes==7
        A.model.geo.synchronize()
        @test Set(A.model.get_entities(1))==Set(((1,11),))
        @test curve_laws()[11].num_nodes==7
    finally
        A.finalize()
    end
end
end

@testset "immediate factories preserve live removal and pending GEO boundaries" begin
    A=Tessella.API
    A.finalize();A.initialize()
    try
        for (tag,(x,y)) in enumerate(((0.,0.),(1.,0.),(1.,1.),(0.,1.)))
            A.model.add_point(x,y,0.;tag)
        end
        for tag in 1:4;A.model.add_line(tag,mod1(tag+1,4);tag);end
        A.model.add_curve_loop([1,2,3,4];tag=1)
        A.model.add_plane_surface([1];tag=1)
        A.model.remove_entities([(2,1)],true)
        @test isempty(A.model.get_entities())
        @test A.model.add_point(0.,0.,0.;tag=1)==1
        @test A.model.get_entities()==[(0,1)]
        @test A.model.add_point(1.,0.,0.;tag=2)==2
        @test A.model.add_point(1.,1.,0.;tag=3)==3
        @test A.model.add_point(0.,1.,0.;tag=4)==4
        for tag in 1:4;@test A.model.add_line(tag,mod1(tag+1,4);tag)==tag;end
        @test A.model.add_curve_loop([1,2,3,4];tag=1)==1
        @test A.model.add_plane_surface([1];tag=1)==1
        @test sort(A.model.get_entities())==vcat([(0,i) for i in 1:4],[(1,i) for i in 1:4],[(2,1)])
    finally
        A.finalize()
    end
    A.initialize()
    try
        @test A.model.geo.add_point(0.,0.,0.;tag=11)==11
        @test A.model.geo.add_point(1.,0.,0.;tag=12)==12
        @test isempty(A.model.get_entities())
        @test A.model.add_point(.5,1.,0.;tag=1)==1
        @test A.model.get_entities()==[(0,1)]
        @test_throws ArgumentError A.model.add_line(11,12;tag=11)
        @test A.model.get_entities()==[(0,1)]
        A.model.geo.synchronize()
        @test sort(A.model.get_entities())==[(0,1),(0,11),(0,12)]
        @test A.model.add_line(11,12;tag=11)==11
        A.model.geo.mesh.set_transfinite_curve(11,5)
        A.mesh.set_transfinite_curve(11,7)
        @test A.CURRENT[].meshing.transfinite_curves[11].num_nodes==7
        @test A.model.add_point(2.,0.,0.;tag=13)==13
        @test A.CURRENT[].meshing.transfinite_curves[11].num_nodes==7
        A.model.synchronize()
        @test A.CURRENT[].meshing.transfinite_curves[11].num_nodes==5
    finally
        A.finalize()
    end
end
