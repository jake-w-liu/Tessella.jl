using Test
using Tessella
using SHA

const _NLS_MODEL=Tessella.Model

function _nls_execute(source;dim=0)
    mktemp() do path,io
        write(io,source)
        dim==1 && write(io,"\nMesh 1;\n")
        close(io)
        return Tessella.execute_geo(path;mesh_dim=dim==1 ? 0 : dim)
    end
end

# Literal Gmsh 4.15.2 recipes from the independently retained M7/N3
# captures. Keep CRLF bytes for the original recipe digest.
function _nls_skew_source(direction)
    text="""
    General.NumThreads=1;Mesh.ElementOrder=1;Mesh.HighOrderOptimize=0;Mesh.Renumber=1;
    Geometry.AutoCoherence=1;Geometry.OldRuledSurface=0;
    Point(1)={0.0,0.0,0.0,1};
    Point(2)={0.0,2.0,0.0,1};
    Point(3)={0.0,2.5,1.5,1};
    Point(4)={0.0,0.5,1.5,1};
    Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
    Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
    Transfinite Curve{1,3}=8 Using Progression 4.0;Transfinite Curve{2,4}=2;
    Transfinite Surface{1}={1,2,3,4};Recombine Surface{1};
    sweep[]=Extrude{$(direction*0.75),0.0,0.0}{Surface{1};Layers{{1,2},{0.2,1.0}};Recombine;QuadTriNoNewVerts RecombLaterals;};
    """
    return replace(text,"\n"=>"\r\n")
end

@testset "Native GEO Line density placement uses primary derivative samples" begin
    # These are actual stored Gmsh parameters, rather than analytic
    # progression fractions. Offset Curve3 exposes the source-coordinate
    # error hidden by a parameter-only tolerance on a length-two Line.
    primary_params=(
        Float64[0,0.00018312569201412332,0.0009156220746376075,
            0.0038455820633995197,0.015565319851520962,
            0.062443862336257205,0.24995639760397706,1],
        Float64[0,0.00018312569201355813,0.0009156220746330919,
            0.0038455820633724355,0.015565319851374652,
            0.06244386233556305,0.24995639760103924,1])
    primary_y=(
        Float64[0,0.00036625138402824665,0.001831244149275215,
            0.0076911641267990395,0.031130639703041925,
            0.12488772467251441,0.4999127952079541,2],
        Float64[2.5,2.499633748615973,2.4981687558507337,
            2.492308835873255,2.468869360297251,
            2.375112275328874,2.0000872047979215,0.5])
    numerical_lengths=(1.9999999999961713,2.0000000000020033)
    recipe_hashes=(
        "0b2f5ed0c62b1f658d42d12d3013e4578c00811c12e8e056d0ec2bcf722c2794",
        "4c7f25480372e3d3f18dc41139d73be3e3da42cc9822237a6bc6baa4b9a62876")
    for (direction,recipe_hash) in zip((-1,1),recipe_hashes)
        source=_nls_skew_source(direction)
        @test bytes2hex(sha256(source))==recipe_hash
        model=_nls_execute(source).model
        generated=Vector{Float64}[]
        for (index,curve) in enumerate((1,3))
            control=model.meshing.transfinite_curves[curve]
            params=_NLS_MODEL._model_curve_transfinite_native_params(
                model,curve,0.,1.,control,"native Line sample regression")
            push!(generated,params)
            @test length(params)==control.num_nodes==8
            @test first(params)==0. && last(params)==1.
            @test all(diff(params).>0)
            @test params≈primary_params[index] atol=2e-11 rtol=0
            @test _NLS_MODEL._model_native_line_sampling_length(
                model,curve,"native Line sample regression")≈numerical_lengths[index] atol=2e-13 rtol=0
            for node in eachindex(params)
                actual=_NLS_MODEL._model_curve_point(
                    model,curve,params[node],"native Line sample regression")
                expected=(0.,primary_y[index][node],index==1 ? 0. : 1.5)
                @test maximum(abs.(actual.-expected))<=2e-11
                copied=(direction*0.75,actual[2],actual[3])
                expected_top=(direction*0.75,expected[2],expected[3])
                @test maximum(abs.(copied.-expected_top))<=2e-11
            end
        end
        source_mesh=_NLS_MODEL.mesh_model_surface(model,1)
        @test size(source_mesh.coords,2)==16
        @test length(source_mesh.blocks)==1
        @test only(source_mesh.blocks).msh==3
        @test size(only(source_mesh.blocks).nodes)==(4,7)
        for (index,curve) in enumerate((1,3))
            @test model.curve_params[curve]==generated[index]
            for node in eachindex(primary_y[index])
                expected=(0.,primary_y[index][node],index==1 ? 0. : 1.5)
                @test minimum(maximum(abs.(source_mesh.coords[:,slot].-expected))
                    for slot in axes(source_mesh.coords,2))<=2e-11
            end
        end
    end
end

@testset "Older native Line sampling controls retain their primary gates" begin
    controls=(
        ("Progression",2.,Float64[0,.03225806585200275,.09677419010543135,
            .2258064432066724,.4838709572937826,1],2e-11),
        ("Bump",2.,Float64[0,.2113246604743475,.3660253944175726,
            .4999999999989506,.6339746055809568,.788675339524792,1],2e-11),
        ("Beta",1.2,Float64[0,.0865292929320693,.2036267615379391,
            .3559899524133273,.544421476246288,.7633519289437464,1],2e-11),
        ("Progression",2.,Float64[0,.06666666781522945,
            .19999999863373244,.46666666554782726,1],3e-12))
    for (law,coefficient,expected,tolerance) in controls
        execution=_nls_execute("""
            Point(1)={0,0,0,1};Point(2)={1,0,0,1};Line(1)={1,2};
            Transfinite Curve{1}=$(length(expected)) Using $law $coefficient;
            """;dim=1)
        mesh=Tessella.geo_entity_mesh(execution,1,1)
        params=execution.model.curve_params[1]
        @test size(mesh.coords,2)==length(expected)
        @test size(mesh.segs,2)==length(expected)-1
        @test params≈expected atol=tolerance rtol=0
        @test mesh.coords[1,:]≈expected atol=tolerance rtol=0
        @test mesh.coords[2:3,:]==zeros(2,length(expected))
        @test params==collect(mesh.coords[1,:])
    end
    # Nonpositive, unit and direct/ungraded arms bypass native density
    # placement, retaining their established uniform counts and endpoints.
    for (law,coefficient) in (("Progression",-2.),("Progression",1.),("Beta",.5))
        execution=_nls_execute("""
            Point(1)={0,0,0,1};Point(2)={1,0,0,1};Line(1)={1,2};
            Transfinite Curve{1}=6 Using $law $coefficient;
            """;dim=1)
        @test execution.model.curve_params[1]≈collect(range(0.,1.;length=6)) atol=2e-15 rtol=0
        @test size(Tessella.geo_entity_mesh(execution,1,1).segs,2)==5
    end
end

@testset "Native Line sampling keeps robust cold-scale evaluation" begin
    # The primary squared norm over/underflows at these finite scales.
    # Placement keeps the previous robust evaluator; count-only primary
    # threshold arithmetic remains a separate unchanged operation.
    for scale in (1e-200,1e200)
        first=(0.,0.,0.);last=(scale,0.,0.)
        gamma=t->(scale*t,0.,0.)
        for t in (0.,.125,.5,.875,1.)
            primary=_NLS_MODEL._transfinite_native_line_speed(first,last,t)
            @test primary==0. || !isfinite(primary)
            previous=_NLS_MODEL._length_point(gamma,nothing,t,0.,1.,"cold-scale regression").xp
            speed=_NLS_MODEL._model_native_line_sampling_speed(first,last,t,gamma,"cold-scale regression")
            @test speed==previous
            @test isfinite(speed) && speed>0
            @test speed≈scale rtol=2e-10
        end
    end
    first=(2.0^50+2.,0.,0.);last=(2.0^50+6.,0.,0.)
    gamma=t->(first[1]+4t,0.,0.)
    @test _NLS_MODEL._transfinite_native_line_speed(first,last,.5)==0.
    @test _NLS_MODEL._model_native_line_sampling_speed(first,last,.5,gamma,"offset regression")==4.
end
