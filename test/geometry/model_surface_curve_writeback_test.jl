module ModelSurfaceCurveWritebackTests

using Test, Tessella, Random
using Tessella.Model: mesh_model_surface, model_to_mixed
const Model=Tessella.Model

function fixture(axes,b,law,winding)
    corners=((0.,0.),(1.,0.),(1.,1.),(0.,1.))
    points=join(("Point($i)={"*join(ntuple(d->d==axes[1] ? p[1] :
        d==axes[2] ? p[2] : 0.,3),",")*",1};" for (i,p) in enumerate(corners)))
    loop=winding==1 ? "1,2,3,4" : "-4,-3,-2,-1"
    return "Geometry.AutoCoherence=0;"*points*
        "Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};"*
        "Curve Loop(1)={$loop};Plane Surface(1)={1};"*
        "Transfinite Curve{1,3}=3 $law;Transfinite Curve{2,4}=$(b+1) $law;"*
        "Transfinite Surface{1};Recombine Surface{1};"
end

function parse_model(text)
    path=tempname()*".geo"
    write(path,text)
    try
        return execute_geo(path;mesh_dim=0).model
    finally
        rm(path;force=true)
    end
end

function match_all(existing,queries,tolerance)
    ordered=issorted(existing) && all(isfinite,existing)
    return [Model._model_surface_curve_parameter_match(
        existing,parameter,tolerance,ordered) for parameter in queries]
end

@testset "Stored curve parameter matching keeps the first admissible value" begin
    # A tolerance interval can contain several stored values: identity is the
    # first stored index, rather than the nearest value or a sorted replacement.
    for existing in ([0.,0.25,0.5,0.75,1.], [0.75,0.25,0.5,0.,1.],
                     [0.,0.,nextfloat(0.),1.], [0.,NaN,0.5,Inf,1.],
                     [-floatmax(),0.,floatmax()], Float64[]),
        tolerance in (0.,eps(),0.3,Inf,-1.,NaN),
        parameter in (-Inf,-0.,0.,nextfloat(0.),0.4,1.,Inf,NaN)
        expected=findall(value->abs(value-parameter)<=tolerance,existing)
        wanted=isempty(expected) ? 0 : first(expected)
        @test only(match_all(existing,[parameter],tolerance))==wanted
    end
    rng=MersenneTwister(418)
    for _ in 1:100
        existing=sort!(randn(rng,30))
        queries=vcat(existing,randn(rng,20))
        for values in (existing,reverse(existing)),tolerance in (0.,1e-12,0.1)
            actual=match_all(values,queries,tolerance)
            expected=[something(findfirst(value->abs(value-parameter)<=tolerance,
                                          values),0) for parameter in queries]
            @test actual==expected
        end
    end
end

@testset "Native surface writeback retains actual sampled chains and owners" begin
    for axes in ((1,2,3),(2,3,1),(3,1,2)),winding in (-1,1),
        law in ("", "Using Progression 2", "Using Bump 2", "Using Beta 1.2")
        model=parse_model(fixture(axes,4,law,winding))
        first_mesh=mesh_model_surface(model,1)
        coordinates=copy(first_mesh.coords)
        cells=copy(only(first_mesh.blocks).nodes)
        parameters=deepcopy(model.curve_params)
        for curve in 1:4
            values=parameters[curve]
            @test issorted(values) && all(diff(values).>0)
            @test first(values)==0. && last(values)==1.
            @test length(values)==(isodd(curve) ? 3 : 5)
        end
        repeated=mesh_model_surface(model,1)
        @test reinterpret(UInt64,vec(repeated.coords))==reinterpret(UInt64,vec(coordinates))
        @test only(repeated.blocks).nodes==cells
        @test model.curve_params==parameters
        projected=model_to_mixed(model,repeated,1)
        @test reinterpret(UInt64,vec(projected.coords))==reinterpret(UInt64,vec(coordinates))
        @test [count(owner->owner[1]==dimension,
                     projected.entity_data.node_entities) for dimension in 0:2]==[4,8,3]
        @test model.curve_params==parameters
    end
end

end # module
