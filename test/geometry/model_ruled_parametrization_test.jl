using Test,Tessella
using Tessella.Model: model_value,model_parametrization,add_ruled_surface!

const _RUP_API=Tessella.API

function _rup_model(corners)
    model=GeoModel()
    for (tag,p) in pairs(corners);add_point!(model,p...;tag);end
    for tag in eachindex(corners)
        add_line!(model,tag,mod1(tag+1,length(corners));tag)
    end
    add_curve_loop!(model,collect(eachindex(corners));tag=1)
    add_ruled_surface!(model,[1];tag=1)
    return model
end

function _rup_with_api(source,action)
    api=_RUP_API;api.initialize()
    try
        mktempdir() do directory
            path=joinpath(directory,"ruled.geo");write(path,source)
            geometry=api.open_geo!(path;mesh_dim=0)
            action(api,geometry)
        end
    finally
        api.finalize()
    end
end

function _rup_has_box(value)
    value===Core.Box && return true
    value isa GlobalRef && return value.mod===Core && value.name===:Box
    value isa QuoteNode && return _rup_has_box(value.value)
    value isa Expr && return any(_rup_has_box,value.args)
    value isa Core.CodeInfo && return any(_rup_has_box,value.code)
    return false
end

@testset "native ruled inverse and computed mesh parameters" begin
    @testset "ruled geometry round trip and native branches" begin
        # Gmsh4.15.2 getValue/getParametrization on this three-generatrix
        # patch: new [.4,.2] -> [.4,.4,.2], old -> [.64,.16,.08].
        triangle=_rup_model(((0.,0.,0.),(2.,0.,0.),(0.,2.,1.)))
        quad=_rup_model(((0.,0.,0.),(2.,0.,0.),(2.,3.,1.),(0.,3.,0.)))
        for old in (false,true),model in (triangle,quad)
            parameters=[.4,.2,.6,.2,.7,.6]
            xyz=model_value(model,2,1,parameters;old_ruled_surface=old)
            recovered=model_parametrization(model,2,1,xyz;old_ruled_surface=old)
            @test recovered≈parameters atol=2e-8 rtol=0
            @test model_value(model,2,1,recovered;old_ruled_surface=old)≈xyz atol=2e-8 rtol=0
            @test model_parametrization(model,2,1,Float64[];old_ruled_surface=old)==Float64[]
            @test_throws ArgumentError model_parametrization(model,2,1,[0.,0.];old_ruled_surface=old)
            @test_throws ArgumentError model_parametrization(model,2,1,[NaN,0.,0.];old_ruled_surface=old)
        end
        @test model_value(triangle,2,1,[.4,.2])≈[.4,.4,.2] atol=4eps(Float64)
        @test model_value(triangle,2,1,[.4,.2];old_ruled_surface=true)≈[.64,.16,.08] atol=4eps(Float64)
        @test model_value(quad,2,1,[.4,.2])≈[.8,.6,.08] atol=4eps(Float64)
        plane=GeoModel()
        for (tag,p) in pairs(((0.,0.,0.),(2.,0.,0.),(0.,2.,0.)))
            add_point!(plane,p...;tag)
        end
        for tag in 1:3;add_line!(plane,tag,mod1(tag+1,3);tag);end
        add_curve_loop!(plane,[1,2,3];tag=1);add_plane_surface!(plane,[1];tag=1)
        for dimension in (1,2)
            point=[.5,.3,1.]
            @test model_parametrization(plane,dimension,1,point;old_ruled_surface=true)==
                model_parametrization(plane,dimension,1,point)
        end
    end

    @testset "global ruled option reaches supported API evaluators" begin
        source="""
        Point(1)={0,0,0};Point(2)={2,0,0};Point(3)={0,2,1};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,1};
        Curve Loop(1)={1,2,3};Surface(1)={1};
        """
        _rup_with_api(source,(api,_)->begin
            @test api.option("Geometry.OldRuledSurface")==0.
            for old in (0,1)
                @test api.option("Geometry.OldRuledSurface",old)==old
                uv=[.4,.2]
                xyz=api.model.get_value(2,1,uv)
                @test xyz≈(old==0 ? [.4,.4,.2] : [.64,.16,.08]) atol=4eps(Float64)
                @test api.model.get_parametrization(2,1,xyz)≈uv atol=2e-8 rtol=0
                model=api.CURRENT[]
                @test api.model.get_derivative(2,1,uv)==Tessella.Model.model_derivative(model,2,1,uv;old_ruled_surface=old!=0)
                @test api.model.get_second_derivative(2,1,uv)==Tessella.Model.model_second_derivative(model,2,1,uv;old_ruled_surface=old!=0)
                @test api.model.get_normal(1,uv)==Tessella.Model.model_normal(model,1,uv;old_ruled_surface=old!=0)
                @test api.model.reparametrize_on_surface(1,3,[.2],1)≈(old==0 ? [.8,.8] : [.8,1.]) atol=4eps(Float64)
                @test api.model.get_parametrization_bounds(2,1)==([0.,0.],[1.,1.])
            end
            @test api.option("Geometry.OldRuledSurface",1.8)==1.
            for invalid in (NaN,Inf,true)
                @test_throws ArgumentError api.option("Geometry.OldRuledSurface",invalid)
                @test api.option("Geometry.OldRuledSurface")==1.
            end
            @test api.LAST_MESH[]===nothing
        end)
    end

    @testset "OCC and discrete inverses retain their own parameters" begin
        _rup_with_api("SetFactory(\"OpenCASCADE\");Box(1)={1,2,3,4,5,6};",(api,_)->begin
            for old in (0,1)
                api.option("Geometry.OldRuledSurface",old)
                parameters=[.25,.75]
                point=api.model.get_value(2,1,parameters)
                @test api.model.get_parametrization(2,1,point)≈parameters atol=2e-12 rtol=0
            end
        end)
        api=_RUP_API;api.initialize()
        try
            api.model.add_discrete_entity(2,1)
            api.mesh.add_nodes(2,1,[11,22,33],[0.,0.,0.,2.,0.,0.,0.,2.,0.],
                [0.,0.,1.,0.,0.,1.])
            api.mesh.add_elements_by_type(1,2,[55],[11,22,33])
            for old in (0,1)
                api.option("Geometry.OldRuledSurface",old)
                @test api.model.get_parametrization(2,1,[.5,.5,0.])==[.25,.25]
                @test api.model.get_value(2,1,[.25,.25])==[.5,.5,0.]
            end
        finally
            api.finalize()
        end
    end

    @testset "Prism18 lateral-node queries use the actual ruled inverse" begin
        # Pinned Gmsh4.15.2 stores [.5,1/3,.5,2/3] for Surface11's
        # two preexisting horizontal-edge nodes, and empty parameters on
        # its three new face-center nodes. Legacy native Tessella queries
        # deliberately compute parameters for every requested native node;
        # this regression checks that contract, not stored-UV provenance.
        for winding in ("1,2,3","-3,-2,-1")
            source="""
            Mesh.TransfiniteTri=1;
            Point(1)={0,0,0,1};Point(2)={1,0,0,1};Point(3)={0,1,0,1};
            Line(1)={1,2};Line(2)={2,3};Line(3)={3,1};
            Curve Loop(1)={$winding};Plane Surface(1)={1};
            Transfinite Curve{:}=2;Transfinite Surface{1};
            s[]=Extrude{0,0,1}{Surface{1};Layers{3};Recombine;QuadTriNoNewVerts RecombLaterals;};
            """
            _rup_with_api(source,(api,geometry)->begin
                api.mesh.generate(3);api.mesh.set_order(2)
                cache=api.LAST_MESH[];class=api.LAST_MESH_CLASS[]
                before=(deepcopy(api.CURRENT[].curve_params),api.mesh.get_elements(),
                    api.mesh.get_max_node_tag(),api.mesh.get_max_element_tag())
                for old in (0,1)
                    api.option("Geometry.OldRuledSurface",old)
                    for surface in Int.(geometry.lists["s"][3:end])
                        own_tags,own_xyz,own_uv=api.mesh.get_nodes(2,surface,false,true)
                        @test length(own_tags)==5 && length(own_uv)==10
                        @test api.model.get_value(2,surface,own_uv)≈own_xyz atol=2e-8 rtol=0
                        @test api.model.get_parametrization(2,surface,own_xyz)≈own_uv atol=2e-8 rtol=0
                        for (position,node) in pairs(own_tags)
                            xyz,uv,dimension,entity=api.mesh.get_node(node)
                            @test (dimension,entity)==(2,surface)
                            @test xyz==own_xyz[3position-2:3position]
                            @test uv==own_uv[2position-1:2position]
                            @test uv≈[.5,xyz[3]] atol=2e-8 rtol=0
                        end
                        for boundary in (false,true),parametric in (false,true)
                            tags,xyz,uv=api.mesh.get_nodes(2,surface,boundary,parametric)
                            @test length(tags)==(boundary ? 21 : 5)
                            @test length(uv)==(parametric ? 2length(tags) : 0)
                            parametric && @test api.model.get_value(2,surface,uv)≈xyz atol=2e-8 rtol=0
                        end
                    end
                end
                @test api.LAST_MESH[]===cache && api.LAST_MESH_CLASS[]===class
                @test (api.CURRENT[].curve_params,api.mesh.get_elements(),
                    api.mesh.get_max_node_tag(),api.mesh.get_max_element_tag())==before
            end)
        end
    end

    @testset "parametrization dispatch has no boxed locals" begin
        @test _rup_has_box(GlobalRef(Core,:Box))
        @test !_rup_has_box(:(sin(1)))
        for function_value in (model_parametrization,_RUP_API._api_model_parametrization,
                _RUP_API._mesh_entity_parameters,_RUP_API._append_p2_midnodes!)
            for method in methods(function_value)
                @test !_rup_has_box(Base.uncompressed_ast(method))
            end
        end
    end
end
