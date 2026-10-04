using Test,Tessella

function _mcp_point_model(points)
    model=Tessella.Model.GeoModel()
    for (tag,p) in enumerate(points)
        Tessella.Model.add_point!(model,p...;tag)
    end
    return model
end

function _mcp_execute(source;dim=3)
    mktempdir() do directory
        path=joinpath(directory,"coherence.geo")
        write(path,source)
        return execute_geo(path;mesh_dim=dim)
    end
end

@testset "coherence uses a local GEO characteristic length" begin
    model_module=Tessella.Model
    # Pinned Gmsh 4.15.2 GEO Coherence uses the unpadded diagonal: the
    # synchronized API's separate padded characteristic length is not used.
    for plane in (:xy,:xz,:x),tol in (1e-8,1e-3),factor in (0.5,1.2,1.5,2.5)
        corners=plane===:xy ? [(0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.)] :
            plane===:xz ? [(0.,0.,0.),(1.,0.,0.),(1.,0.,1.),(0.,0.,1.)] :
                         [(0.,0.,0.),(1.,0.,0.)]
        count=length(corners)
        push!(corners,(tol*factor,0.,0.))
        for shift in (0.,2.0^50)
            # Shift the unused plane coordinate so the test separation remains
            # representable and independent of world-coordinate magnitude.
            points=[plane===:xy ? (p[1],p[2],p[3]+shift) :
                                 (p[1],p[2]+shift,p[3]) for p in corners]
            model=_mcp_point_model(points)
            scale=plane===:x ? 1. : sqrt(2.)
            @test model_module._coherence_eps(model,tol)≈tol*scale rtol=2eps(Float64)
            merged=plane===:x ? factor<1. : factor<sqrt(2.)
            @test model_module.coherence!(model;tol)==merged
            @test length(model.points)==count+1-Int(merged)
            @test haskey(model.points,count+1)==!merged
        end
    end
    for shift in (0.,2.0^50),tol in (0.,1e-8)
        model=_mcp_point_model([(0.,-0.,shift),(-0.,0.,shift)])
        @test model_module._coherence_eps(model,tol)==tol
        @test model_module.coherence!(model;tol)
        @test collect(keys(model.points))==[1]
    end
    # Tiny positive tolerances and exact-zero underflow retain distinct Points
    # without either overflowing Int bins or introducing approximate merges.
    for tol in (0.,nextfloat(0.),1e-300)
        model=_mcp_point_model([(0.,0.,0.),(1.,0.,0.),(1.,0.,0.),
                                (2.,0.,0.),(2.,0.,0.)])
        @test model_module.coherence!(model;tol)
        @test sort!(collect(keys(model.points)))==[1,2,4]
    end
    model=_mcp_point_model([(0.,0.,0.),(1e-200,0.,0.),(1e-200,0.,0.)])
    @test model_module._coherence_eps(model,nextfloat(0.))==0.
    @test model_module.coherence!(model;tol=nextfloat(0.))
    @test sort!(collect(keys(model.points)))==[1,2]
    # Rounded origin subtraction previously placed these close Points two
    # bins apart at a quotient near 2^53, despite satisfying actual tolerance.
    model=_mcp_point_model([(-1.,0.,0.),(nextfloat(1.),0.,0.),
        (nextfloat(nextfloat(1.)),0.,0.)])
    tolerance=2.0^-53
    @test model_module._points_close(model.points[2],model.points[3],
                                    model_module._coherence_eps(model,tolerance))
    @test model_module.coherence!(model;tol=tolerance)
    @test sort!(collect(keys(model.points)))==[1,2]
    origin=(-floatmax(Float64),0.,0.)
    positive=model_module._coherence_spatial_cell((floatmax(Float64),0.,0.),origin,1e-300)
    negative=model_module._coherence_spatial_cell(origin,(floatmax(Float64),0.,0.),1e-300)
    @test positive==(model_module._COHERENCE_BIN_LIMIT,0,0)
    @test negative==(-model_module._COHERENCE_BIN_LIMIT,0,0)
    @test positive[1]+1>positive[1] && negative[1]-1<negative[1]
    model=_mcp_point_model([origin,(floatmax(Float64),0.,0.),origin])
    @test isfinite(model_module._coherence_eps(model))
    @test model_module.coherence!(model)
    @test sort!(collect(keys(model.points)))==[1,2]
    extreme=floatmax(Float64)
    model=_mcp_point_model([(-extreme,-extreme,-extreme),
                           (extreme,extreme,extreme),(-extreme,-extreme,-extreme)])
    @test isfinite(model_module._coherence_eps(model))
    @test model_module._coherence_eps(model)≈Float64(
        2BigFloat(extreme)*sqrt(BigFloat(3))*BigFloat(1e-8)) rtol=4eps(Float64)
    @test model_module.coherence!(model)
    @test sort!(collect(keys(model.points)))==[1,2]
end

@testset "coherence keeps closed OCC-circle extents" begin
    model_module=Tessella.Model
    model=model_module.GeoModel()
    p=model_module.add_point!(model,1e8,0,2.0^50)
    model_module._add_occ_circle!(model,p,p,(0.,0.,2.0^50),(0.,0.,1.),
        (1.,0.,0.),(0.,1.,0.),1e8,0.,2π)
    copy=model_module.add_point!(model,1e8,-2e-8,2.0^50)
    @test model_module._coherence_eps(model)≈2sqrt(2.) rtol=2eps(Float64)
    @test model_module.coherence!(model)
    @test haskey(model.points,p) && !haskey(model.points,copy)
    # Seam lookup shares the relative bins and still checks actual distance.
    point=(1.,2.,2.0^50)
    tolerance=1e-8
    cell=model_module._coherence_spatial_cell(point,point,tolerance)
    index=(Dict(cell=>[point]),tolerance,point)
    @test model_module._extrude_quadtri_snap_seam(index,point,"precision test")==point
    @test_throws r"does not match" model_module._extrude_quadtri_snap_seam(
        index,(2.,2.,2.0^50),"precision test")
end

@testset "direct GEO extrusion preserves shifted source identities" begin
    for laterals in (false,true)
        source="""
        B=2^50;
        Point(1)={0,0,B,1};Point(2)={1,0,B,1};
        Point(3)={1,1,B,1};Point(4)={0,1,B,1};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
        Transfinite Curve{1,2,3,4}=2;Transfinite Surface{1} Left;
        Extrude{0,0,1}{Surface{1};Layers{4};Recombine;
          QuadTriNoNewVerts $(laterals ? "RecombLaterals" : "");}
        """
        result=_mcp_execute(source)
        @test length(result.model.points)==8
        @test all(result.model.curves[c][1]!=result.model.curves[c][2] for c in 1:4)
        volume=geo_entity_mesh(result,3,1)
        @test Tessella.MeshTypes.nnodes(volume)==20
        @test Tessella.MeshTypes.validate(volume).ok
        @test sort!(unique(volume.coords[3,:])).-2.0^50==[0.,0.25,0.5,0.75,1.]
        if laterals
            @test only(volume.blocks).msh==6
            @test size(only(volume.blocks).nodes,2)==8
        else
            @test size(volume.tets,2)==24
        end
    end
end
