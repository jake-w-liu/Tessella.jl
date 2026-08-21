# Gmsh-compatible scalar/size fields and field-driven tetrahedral refinement.

using Test
using LinearAlgebra
using Tessella
using Tessella.Mesh1D: mesh_segment
using Tessella.MeshTypes: Mesh, nnodes, ntets, unique_edges, validate, node, tet_volume
using Tessella.Geometry: box_surface
using Tessella.IO: read_geo_params, GeoParams, GeoFieldSpec

function _sizefield_intersection_checksum(a::Metric3,b::Metric3,n::Int)
    checksum=0.0
    @inbounds for _ in 1:n
        checksum+=Tessella.SizeField.intersection_alauzet(a,b).m11
    end
    return checksum
end

function _nearest_point_checksum(field,point,n::Int)
    distance=0.0;identity=0
    @inbounds for _ in 1:n
        d,i=Tessella.SizeField._nearest_point(field,point)
        distance+=d;identity+=i
    end
    return distance+identity
end

function _field_value_checksum(field,point,n::Int)
    value=0.0
    @inbounds for _ in 1:n
        value+=field_value(field,point)
    end
    return value
end

function _size_at_checksum(field,point,n::Int)
    value=0.0
    @inbounds for _ in 1:n
        value+=size_at(field,point)
    end
    return value
end

function _metric_eigenvalue_checksum(metric,n::Int)
    value=0.0
    @inbounds for _ in 1:n
        eigenvalues=Tessella.SizeField.metric_eigenvalues(metric)
        value+=eigenvalues[1]+eigenvalues[2]+eigenvalues[3]
    end
    return value
end

function _attractor_metric_checksum(field,n::Int)
    value=0.0
    @inbounds for _ in 1:n
        value+=metric_at(field,0.2,0.3,0.4).m11
    end
    return value
end

function _postview_checksum(field,point,n::Int)
    value=0.0
    @inbounds for _ in 1:n
        value+=field_value(field,point...)
    end
    return value
end

@testset "Gmsh-compatible size fields" begin
    @testset "distance to discrete entities" begin
        fp=DistanceField(points=[(0.0,0.0,0.0),(10.0,0.0,0.0)])
        @test field_value(fp,3.0,4.0,0.0)==5.0
        @test field_value(fp,9.0,0.0,0.0)==1.0

        fs=DistanceField(segments=[((0.0,0.0,0.0),(2.0,0.0,0.0))])
        @test field_value(fs,1.0,3.0,4.0)==5.0
        @test field_value(fs,-1.0,0.0,0.0)==1.0

        ft=DistanceField(triangles=[((0.0,0.0,0.0),(1.0,0.0,0.0),(0.0,1.0,0.0))])
        @test field_value(ft,0.25,0.25,2.0)==2.0
        @test field_value(ft,2.0,0.0,0.0)==1.0

        tm=Mesh(Float64[0 1 0;0 0 1;0 0 0];tris=reshape(Int32[1,2,3],3,1))
        @test field_value(DistanceField(tm),0.25,0.25,0.5)==0.5
        @test_throws ArgumentError DistanceField()
        @test_throws ArgumentError field_value(fp,NaN,0.0,0.0)
        @test_throws ArgumentError field_value(fp,(0.0,))
        @test_throws ArgumentError size_at(ConstantSize(1.0),(0.0,))

        # The acceleration hierarchy must return the identical minimum as an
        # exhaustive scan of the same primitive kernels, including across leaves.
        points=[(sin(i),cos(2i),0.1i) for i in 1:24]
        segments=[((0.2i,-1.0,0.1i),(-0.1i,1.0,-0.05i)) for i in 1:18]
        triangles=[((0.1i,0.0,-0.2),(0.1i+0.3,0.1,0.0),
                    (0.1i-0.1,0.4,0.2)) for i in 1:18]
        accelerated=DistanceField(;points=points,segments=segments,triangles=triangles)
        @test length(accelerated.bvh_left)>1
        for k in 1:30
            q=(sin(0.7k)*2,cos(0.3k)*1.5,sin(0.2k))
            expected=minimum(vcat(
                [hypot(q[1]-p[1],q[2]-p[2],q[3]-p[3]) for p in points],
                [Tessella.SizeField._distance_segment(q,s...) for s in segments],
                [Tessella.SizeField._distance_triangle(q,t...) for t in triangles]))
            @test field_value(accelerated,q)==expected
        end
        _field_value_checksum(accelerated,(0.1,0.2,0.3),1)
        @test (@allocated _field_value_checksum(
            accelerated,(0.1,0.2,0.3),10_000))<=64
    end

    @testset "Threshold, Box, Min/Max, and final bounds" begin
        distance=DistanceField(points=[(0.0,0.0,0.0)])
        linear=ThresholdField(distance;dist_min=1.0,dist_max=3.0,
                              size_min=0.1,size_max=0.5)
        @test size_at(linear,0.0,0.0,0.0)==0.1
        @test size_at(linear,2.0,0.0,0.0)≈0.3
        @test size_at(linear,4.0,0.0,0.0)==0.5

        sigmoid=ThresholdField(distance;dist_min=1.0,dist_max=3.0,
                               size_min=0.1,size_max=0.5,sigmoid=true)
        logistic(r)=exp(12r-6)/(1+exp(12r-6))
        @test size_at(sigmoid,0.0,0.0,0.0)≈0.1*(1-logistic(0.0))+0.5logistic(0.0)
        @test size_at(sigmoid,2.0,0.0,0.0)≈0.3
        stopped=ThresholdField(distance;dist_min=1.0,dist_max=3.0,
                               size_min=0.1,size_max=0.5,stop_at_dist_max=true)
        @test size_at(stopped,3.0,0.0,0.0)==Tessella.SizeField.GMSH_MAX_SIZE

        box=BoxField(0.0,1.0,0.0,1.0,0.0,1.0;vin=0.2,vout=1.0,thickness=0.5)
        @test size_at(box,0.5,0.5,0.5)==0.2
        @test size_at(box,1.0,0.5,0.5)==0.2
        @test size_at(box,1.25,0.5,0.5)≈0.6
        @test size_at(box,2.0,0.5,0.5)==1.0
        corner_distance=hypot(0.25,0.25)
        @test size_at(box,1.25,1.25,0.5)≈0.2+corner_distance/0.5*0.8

        mn=MinSize([linear,box]); mx=MaxSize([linear,box])
        @test size_at(mn,0.5,0.5,0.5)==min(size_at(linear,0.5,0.5,0.5),0.2)
        @test size_at(mx,0.5,0.5,0.5)==max(size_at(linear,0.5,0.5,0.5),0.2)
        bounded=BoundedSize(distance;size_min=0.05,size_max=0.4,factor=2.0)
        @test size_at(bounded,0.0,0.0,0.0)==0.1
        @test size_at(bounded,10.0,0.0,0.0)==0.8

        @test_throws ArgumentError ThresholdField(distance;dist_min=1.0,dist_max=1.0)
        reversed=BoxField(1,0,0,1,0,1;vin=0.2,vout=1.0)
        @test field_value(reversed,0.5,0.5,0.5)==1.0
        reversed_layer=BoxField(1,0,0,1,0,1;vin=0.2,vout=1.0,thickness=0.5)
        @test field_value(reversed_layer,0.5,0.5,0.5)==0.2
        @test field_value(BoxField(0,1,0,1,0,1;vin=0.2,vout=1.0,
                                  thickness=-1),1.1,0.5,0.5)==1.0
        @test field_value(MinSize(AbstractField[]),0,0,0)==
              Tessella.SizeField.GMSH_MAX_SIZE
        @test field_value(MaxSize(AbstractField[]),0,0,0)==
              -Tessella.SizeField.GMSH_MAX_SIZE
        @test_throws ArgumentError BoundedSize(distance;size_min=1.0,size_max=0.5)
    end

    @testset "Ball, finite Cylinder, and Frustum" begin
        ball=BallField((1.0,2.0,3.0),2.0;vin=0.1,vout=0.9,thickness=1.0)
        @test size_at(ball,1.0,2.0,3.0)==0.1
        @test size_at(ball,3.0,2.0,3.0)==0.1 # transition starts at the open sphere
        @test size_at(ball,3.5,2.0,3.0)≈0.5
        @test size_at(ball,4.0,2.0,3.0)==0.9
        sharp=BallField((0.0,0.0,0.0),1.0;vin=0.2,vout=1.0)
        @test size_at(sharp,1.0,0.0,0.0)==1.0

        cylinder=CylinderField((0.0,0.0,0.0),(0.0,0.0,2.0),1.0;
                               vin=0.2,vout=1.0)
        @test size_at(cylinder,0.5,0.0,1.0)==0.2
        @test size_at(cylinder,1.0,0.0,0.0)==1.0
        @test size_at(cylinder,0.0,0.0,2.0)==1.0
        @test size_at(cylinder,0.0,0.0,-1.999)==0.2
        @test field_value(BallField((0,0,0),-1;vin=0.2,vout=1),0,0,0)==1
        @test field_value(BallField((0,0,0),1;vin=0.2,vout=1,
                                   thickness=-1),1,0,0)==1
        @test_throws ArgumentError BallField(0,1;vin=1)
        @test field_value(CylinderField((0,0,0),(0,0,0),1;
                                       vin=0.2,vout=1),0,0,0)==1
        @test field_value(CylinderField((0,0,0),(0,0,1),-1;
                                       vin=0.2,vout=1),0,0,0)==0.2

        frustum=FrustumField((0.0,0.0,0.0),(0.0,0.0,2.0);
            inner_r1=0.0,outer_r1=1.0,inner_r2=0.0,outer_r2=1.0,
            inner_v1=0.1,outer_v1=1.0,inner_v2=0.2,outer_v2=0.8)
        @test size_at(frustum,0.0,0.0,1.0)≈0.15
        @test size_at(frustum,0.5,0.0,1.0)≈0.525
        @test size_at(frustum,1.0,0.0,1.0)≈0.9
        @test size_at(frustum,1.01,0.0,1.0)==Tessella.SizeField.GMSH_MAX_SIZE
        @test size_at(frustum,0.0,0.0,-0.01)==Tessella.SizeField.GMSH_MAX_SIZE
        @test field_value(FrustumField((0,0,0),(0,0,0)),0,0,0)==
              Tessella.SizeField.GMSH_MAX_SIZE
        reversed_frustum=FrustumField((0,0,0),(0,0,1);
            inner_r1=1,outer_r1=0,inner_r2=1,outer_r2=0,
            inner_v1=0.1,outer_v1=0.9,inner_v2=0.1,outer_v2=0.9)
        @test field_value(reversed_frustum,0.25,0,0.5)≈0.7
        @test field_value(FrustumField((0,0,0),(0,0,1);
            inner_r1=1,outer_r1=1,inner_r2=1,outer_r2=1),1,0,0.5)==
              Tessella.SizeField.GMSH_MAX_SIZE
    end

    @testset "graded 1-D and local 3-D refinement" begin
        distance=DistanceField(points=[(0.0,0.0,0.0)])
        field=ThresholdField(distance;dist_min=0.0,dist_max=1.0,
                             size_min=0.1,size_max=0.5)
        points,_=mesh_segment((-2.0,0.0,0.0),(2.0,0.0,0.0),field)
        gaps=[abs(points[i+1][1]-points[i][1]) for i in 1:length(points)-1]
        mids=[(points[i+1][1]+points[i][1])/2 for i in 1:length(points)-1]
        near=[gaps[i] for i in eachindex(gaps) if abs(mids[i])<0.25]
        far=[gaps[i] for i in eachindex(gaps) if abs(mids[i])>1.25]
        @test !isempty(near) && !isempty(far)
        @test maximum(near)<minimum(far)

        coarse=Tessella.Mesh3D.mesh_box(0.0,1.0,0.0,1.0,0.0,1.0;hmax=1.0)
        localfield=BoxField(0.0,0.5,0.0,1.0,0.0,1.0;vin=0.4,vout=2.0)
        refined=refine_to_size(coarse,localfield)
        @test validate(refined).ok
        @test nnodes(refined)>nnodes(coarse) && ntets(refined)>ntets(coarse)
        edges=unique_edges(refined.tris,refined.tets)
        left=Float64[];right=Float64[]
        for (a,b) in edges
            p=Tuple(refined.coords[:,a]);q=Tuple(refined.coords[:,b])
            mid=((p[1]+q[1])/2,(p[2]+q[2])/2,(p[3]+q[3])/2)
            len=hypot(p[1]-q[1],p[2]-q[2],p[3]-q[3])
            target=min(size_at(localfield,p),size_at(localfield,q),size_at(localfield,mid))
            @test len<=target
            mid[1]<0.4 && push!(left,len)
            mid[1]>0.75 && push!(right,len)
        end
        @test !isempty(left) && !isempty(right)
        @test maximum(left)<=0.4
        @test maximum(right)>0.4

        surfaced=mesh_sized(box_surface(0.0,1.0,0.0,1.0,0.0,1.0);field=localfield)
        @test validate(surfaced).ok
        @test_throws ArgumentError mesh_sized(box_surface(0,1,0,1,0,1))
        @test_throws ArgumentError mesh_sized(box_surface(0,1,0,1,0,1);
                                              hmax=1.0,field=localfield)

        # A constrained boundary edge cannot be split by moving the new point a few
        # ULPs off the edge: doing so changes the piecewise-linear domain. This
        # enclosure-derived sliver has no distinct, exactly collinear Float64 split
        # point at the attempted midpoint and must block instead of changing volume.
        a=(0.171184375,0.025476562500000004,0.14882812499999998)
        b=(0.17144062499999999,0.024763020833333337,0.151171875)
        c=(0.17131249999999998,0.02511979166666667,0.15)
        d=(0.171184375,0.025476562500000004,0.151171875)
        midpoint=ntuple(i -> a[i]+(b[i]-a[i])/2,3)
        @test midpoint==c
        sliver=Mesh(hcat(collect(a),collect(b),collect(c),collect(d));
                    tets=reshape(Int32[1,2,4,3],4,1))
        err=try
            refine_to_size(sliver,0.0024)
            nothing
        catch caught
            caught
        end
        @test err isa ErrorException
        @test occursin("constrained edge geometry cannot be moved",sprint(showerror,err))

        # 0.2 is not the exact rational midpoint of Float64(0.1) and Float64(0.3),
        # but it is exactly collinear on this axis-aligned boundary edge. Preserve
        # that valid representable split instead of over-rejecting it.
        ea=(0.1,0.0,0.0);eb=(0.3,0.0,0.0)
        ec=(0.2,0.01,0.0);ed=(0.2,0.0,0.01)
        ep=ntuple(i -> ea[i]+(eb[i]-ea[i])/2,3)
        R=Rational{BigInt}
        @test R(ep[1])!=(R(ea[1])+R(eb[1]))/2
        boundary_sliver=Mesh(hcat(collect(ea),collect(eb),collect(ec),collect(ed));
                             tets=reshape(Int32[1,2,3,4],4,1))
        boundary_split=refine_to_size(boundary_sliver,0.15)
        volume(mesh)=sum(tet_volume(node(mesh,mesh.tets[1,t]),node(mesh,mesh.tets[2,t]),
                                    node(mesh,mesh.tets[3,t]),node(mesh,mesh.tets[4,t]))
                         for t in 1:ntets(mesh))
        @test validate(boundary_split).ok
        @test ntets(boundary_split)==2 && volume(boundary_split)≈volume(boundary_sliver)

        # For a general 3-D edge, independently rounding the three canonical
        # midpoint coordinates can put the represented point infinitesimally off the
        # exact represented line. A unique canonical midpoint with valid exact child
        # orientations remains a usable boundary split; only shifted collision
        # fallbacks are forbidden above.
        ca=(0.7071067811865476,0.7071067811865475,0.0)
        cb=(0.3826834323650899,0.9238795325112867,2.0)
        cp=ntuple(i -> ca[i]+(cb[i]-ca[i])/2,3)
        crossxy=(R(cb[1])-R(ca[1]))*(R(cp[2])-R(ca[2]))-
                (R(cb[2])-R(ca[2]))*(R(cp[1])-R(ca[1]))
        @test !iszero(crossxy)
        cu=(cp[1]+0.05,cp[2],cp[3]);cv=(cp[1],cp[2]+0.05,cp[3])
        canonical_tet=Mesh(hcat(collect(ca),collect(cb),collect(cu),collect(cv));
                           tets=reshape(Int32[1,2,3,4],4,1))
        canonical_split=refine_to_size(canonical_tet,1.5)
        @test validate(canonical_split).ok
        @test nnodes(canonical_split)==5 && ntets(canonical_split)==2
        @test volume(canonical_split)≈volume(canonical_tet)
    end

    @testset ".geo background-field graph" begin
        gp=read_geo_params(joinpath(@__DIR__,"fixtures","enclosure_coax_junction.geo"))
        surface(x)=Mesh(Float64[x x+0.001 x; 0 0.001 0.001; 0 0 0];
                        tris=reshape(Int32[1,2,3],3,1))
        line=Mesh(Float64[0.1708 0.1716;0.1605 0.1605;0.15 0.15];
                  segs=reshape(Int32[1,2],2,1))
        entities=Dict{Tuple{Int,String},Mesh}(
            (2,"sm_coax_pin_pec")=>surface(0.17),
            (2,"sm_p1_surface_t")=>surface(0.17),
            (1,"sm_p1_line_t")=>line)
        field=build_geo_size_field(gp,entities)
        @test field isa AbstractSizeField
        @test size_at(field,0.1708,0.1605,0.15)≈gp.mesh_size_min
        @test size_at(field,0.0,0.0,0.0)==0.009
        @test size_at(field,1.0,1.0,1.0)==0.012
        @test_throws ArgumentError build_geo_size_field(gp,Dict())

        badkind=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),
            Dict(1=>GeoFieldSpec(1,"Unsupported",Dict{String,String}())),1)
        @test_throws ArgumentError build_geo_size_field(badkind,Dict())
        cycle=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),Dict(
            1=>GeoFieldSpec(1,"Min",Dict("FieldsList"=>"{2}")),
            2=>GeoFieldSpec(2,"Max",Dict("FieldsList"=>"{1}"))),1)
        @test_throws ArgumentError build_geo_size_field(cycle,Dict())
        rawdistance=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),
            Dict(1=>GeoFieldSpec(1,"Distance",Dict("PointsList"=>"{p}"))),1)
        point=Mesh(reshape(Float64[0,0,0],3,1))
        point_entities=Dict((0,"p")=>point)
        rawleaf=Tessella.SizeField._build_geo_field(rawdistance,point_entities,1)
        @test field_value(rawleaf,0.0,0.0,0.0)==0.0
        rawfield=build_geo_size_field(rawdistance,point_entities)
        @test_throws ArgumentError size_at(rawfield,0.0,0.0,0.0)
        @test size_at(rawfield,1.0,0.0,0.0)==1.0
        legacy_distance=GeoParams(NaN,NaN,1.0,0,
            Dict{Tuple{Int,Int},String}(),Dict(1=>GeoFieldSpec(1,"Distance",
                Dict("FieldX"=>"7","FieldY"=>"8","FieldZ"=>"9"))),1)
        @test size_at(build_geo_size_field(legacy_distance,Dict()),0,0,0)==1.0
        bad_legacy_distance=GeoParams(NaN,NaN,1.0,0,
            Dict{Tuple{Int,Int},String}(),Dict(1=>GeoFieldSpec(1,"Distance",
                Dict("FieldX"=>"not-an-integer"))),1)
        @test_throws ArgumentError build_geo_size_field(bad_legacy_distance,Dict())
        mktempdir() do dir
            aliases=joinpath(dir,"distance_aliases.geo")
            write(aliases,"""
                Field[1] = Distance;
                Field[1].PointsList = {p1};
                Field[1].NodesList = {p2};
                Background Field = 1;
                """)
            p1=Mesh(reshape(Float64[0,0,0],3,1))
            p2=Mesh(reshape(Float64[10,0,0],3,1))
            aliased=build_geo_size_field(read_geo_params(aliases),
                Dict((0,"p1")=>p1,(0,"p2")=>p2))
            @test size(aliased.input.points,2)==1
            @test aliased.input.points[:,1]==[10.0,0.0,0.0]
        end
        curve_mesh=Mesh(Float64[0 1;0 0;0 0];segs=reshape(Int32[1,2],2,1))
        for (sampling,expected) in ((2,1.0),(3,0.5),(20,1/19))
            sampled=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),Dict(
                1=>GeoFieldSpec(1,"Distance",Dict("CurvesList"=>"{c}",
                                                   "Sampling"=>string(sampling)))),1)
            @test size_at(build_geo_size_field(sampled,Dict((1,"c")=>curve_mesh)),
                          0,0,0)≈expected
        end
        oversampled=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),Dict(
            1=>GeoFieldSpec(1,"Distance",Dict("CurvesList"=>"{c}",
                                               "Sampling"=>"20"))),1)
        @test_throws ArgumentError build_geo_size_field(oversampled,
            Dict((1,"c")=>curve_mesh);max_distance_samples=17)
        surface_mesh=Mesh(Float64[0 1 0;0 0 1;0 0 0];
                          tris=reshape(Int32[1,2,3],3,1))
        sampled_surface=GeoParams(NaN,NaN,1.0,0,
            Dict{Tuple{Int,Int},String}(),Dict(1=>GeoFieldSpec(1,"Distance",
                Dict("SurfacesList"=>"{s}","Sampling"=>"2"))),1)
        @test size_at(build_geo_size_field(sampled_surface,
            Dict((2,"s")=>surface_mesh)),1,0,0)≈0.5

        analytic=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),Dict(
            1=>GeoFieldSpec(1,"Ball",Dict("XCenter"=>"0","YCenter"=>"0",
                "ZCenter"=>"0","Radius"=>"1","Thickness"=>"1",
                "VIn"=>"0.1","VOut"=>"0.8")),
            2=>GeoFieldSpec(2,"Cylinder",Dict("XCenter"=>"0","YCenter"=>"0",
                "ZCenter"=>"0","XAxis"=>"2","YAxis"=>"0","ZAxis"=>"0",
                "Radius"=>"0.5","VIn"=>"0.2","VOut"=>"0.9")),
            3=>GeoFieldSpec(3,"Min",Dict("FieldsList"=>"{1,2}"))),3)
        analytic_field=build_geo_size_field(analytic,Dict())
        @test size_at(analytic_field,0.0,0.0,0.0)==0.1
        @test size_at(analytic_field,1.5,0.0,0.0)==0.2
        @test size_at(analytic_field,0.0,1.5,0.0)≈0.45
        @test size_at(analytic_field,3.0,0.0,0.0)==0.8

        frustum_graph=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),Dict(
            1=>GeoFieldSpec(1,"Frustum",Dict("Z1"=>"0","Z2"=>"2",
                "R1_inner"=>"0","R1_outer"=>"1","R2_inner"=>"0",
                "R2_outer"=>"1","V1_inner"=>"0.1","V1_outer"=>"1",
                "V2_inner"=>"0.2","V2_outer"=>"0.8"))),1)
        parsed_frustum=build_geo_size_field(frustum_graph,Dict())
        @test size_at(parsed_frustum,0.5,0.0,1.0)≈0.525
        aliased_frustum=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),Dict(
            1=>GeoFieldSpec(1,"Frustum",Dict("InnerR1"=>"0.75","R1_inner"=>"0.25"),
                            ["R1_inner","InnerR1"])),1)
        @test build_geo_size_field(aliased_frustum,Dict()).input.inner_r1==0.75
    end

    @testset "remaining Gmsh 4.15.2 field catalog" begin
        # MathEval vs independent Julia evaluation of the same closed form.
        me=MathEvalField("cos(4*3.14*x)*sin(4*3.14*y)/10+0.101")
        oracle_me(x,y)=cos(4*3.14*x)*sin(4*3.14*y)/10+0.101
        @test size_at(me,0.1,0.2,0.0)≈oracle_me(0.1,0.2)
        @test size_at(me,0.0,0.0,0.0)≈oracle_me(0.0,0.0)
        @test_throws ArgumentError MathEvalField("")
        @test_throws ArgumentError field_value(MathEvalField("log(-1.0)"),1,0,0)
        @test field_value(MathEvalField("Pi+x"),1.0,0.0,0.0)≈π+1
        @test field_value(MathEvalField("Sin(Pi/2)"),0.0,0.0,0.0)≈1.0
        @test field_value(MathEvalField("-2^2"),0.0,0.0,0.0)==4.0
        @test field_value(MathEvalField("5%2 + sign(0) + step(0)"),0.0,0.0,0.0)==3.0
        @test field_value(MathEvalField("sum(1,2,3)+max(2,5,4)+min(2,5,4)+med(2,4)+e"),
                          0.0,0.0,0.0)≈16+ℯ
        @test field_value(MathEvalField("fac(4.2)+round(-1.7)"),0.0,0.0,0.0)==23.0
        @test_throws ArgumentError MathEvalField("2^3^2")
        @test_throws ArgumentError MathEvalField("2^-2")
        @test_throws ArgumentError MathEvalField("0.1+0.01*-2^2")
        @test_throws ArgumentError MathEvalField("--2")
        @test_throws ArgumentError MathEvalField("rand()")
        @test_throws ArgumentError MathEvalField("atan2(1,2)")
        for incompatible in ("π+x","X","PI","E","SIN(x)")
            @test_throws ArgumentError MathEvalField(incompatible)
        end
        _size_at_checksum(me,(0.1,0.2,0.0),1)
        @test (@allocated _size_at_checksum(me,(0.1,0.2,0.0),10_000))<=64

        # F-tag composition: F1 is Distance to origin, MathEval F1+0.05.
        dist=DistanceField(points=[(0.0,0.0,0.0)])
        composed=MathEvalField("F1 + 0.05"; fields=Dict(1=>dist))
        @test size_at(composed,3.0,4.0,0.0)==5.05
        @test field_value(MathEvalField("F999+1"),0,0,0)==
              Tessella.SizeField.GMSH_MAX_SIZE+1
        _size_at_checksum(composed,(3.0,4.0,0.0),1)
        @test (@allocated _size_at_checksum(composed,(3.0,4.0,0.0),10_000))<=64

        r2=MathEvalField("x^2 + y^2 + z^2")
        # Independent derivative/laplacian/mean of r².
        @test field_value(GradientField(r2;kind=0,delta=1e-4),1.0,0.0,0.0)≈2.0 atol=1e-6
        @test field_value(GradientField(r2;kind=3,delta=1e-4),0.0,3.0,4.0)≈10.0 atol=1e-5
        @test field_value(LaplacianField(r2;delta=1e-3),0.5,-0.2,0.3)≈6.0 atol=1e-6
        @test field_value(MeanField(ConstantSize(0.4);delta=0.1),0,0,0)≈0.4
        @test_throws ArgumentError GradientField(r2;kind=4)
        huge=ConstantSize(1e308)
        @test_throws ArgumentError field_value(MeanField(huge),0.0,0.0,0.0)
        @test_throws ArgumentError field_value(LaplacianField(huge),0.0,0.0,0.0)
        @test_throws ArgumentError field_value(MaxEigenHessianField(huge),0.0,0.0,0.0)

        # div(grad r / |grad r|) = 2/r for r=|x|.
        radial=MathEvalField("sqrt(x^2+y^2+z^2)")
        @test field_value(CurvatureField(radial;delta=1e-4),3.0,0.0,0.0)≈2/3 atol=2e-3
        @test field_value(MaxEigenHessianField(MathEvalField("x^2");delta=1e-3),0,0,0)≈2 atol=1e-4

        lon=LonLatField(MathEvalField("0.2 + 0.1*x"))  # x here is longitude
        @test field_value(lon,1.0,0.0,0.0)≈0.2 + 0.1*atan(0.0,1.0)
        @test field_value(LonLatField(ConstantSize(0.2);radius=1.0),0.0,0.0,2.0)==
              Tessella.SizeField.GMSH_MAX_SIZE
        @test field_value(CurvatureField(ConstantSize(0.2)),0.0,0.0,0.0)==
              Tessella.SizeField.GMSH_MAX_SIZE
        par=ParametricField(MathEvalField("x+y+z"); fx="2*x", fy="y", fz="0")
        @test field_value(par,1.0,3.0,9.0)==2.0+3.0+0.0

        origin=(0.0,0.0,0.0); delta=(1.0,1.0,1.0)
        data=Array{Float64}(undef,2,2,2)
        # v = x + 2y + 3z on the unit cube corners.
        for i in 0:1, j in 0:1, k in 0:1
            data[i+1,j+1,k+1]=i+2j+3k
        end
        st=StructuredField(origin,delta,data)
        @test field_value(st,0.0,0.0,0.0)==0.0
        @test field_value(st,1.0,1.0,1.0)==6.0
        @test field_value(st,0.5,0.5,0.5)≈0.5+1.0+1.5
        stout=StructuredField(origin,delta,data;outside=7.0)
        @test field_value(stout,2.0,0.0,0.0)==7.0
        singleton=StructuredField(origin,(0.0,0.0,0.0),reshape([0.2],1,1,1))
        @test field_value(singleton,9.0,-2.0,4.0)==0.2
        @test_throws ArgumentError StructuredField(origin,(0.0,1.0,1.0),
                                                    zeros(2,1,1))
        @test_throws ArgumentError StructuredField(origin,delta,fill(NaN,1,1,1))
        mktempdir() do dir
            ascii=joinpath(dir,"grid.txt")
            write(ascii,"0 0 0\n1 1 1\n2 2 2\n0 3 2 5 1 4 3 6\n")
            loaded=StructuredField(ascii;text=true,max_samples=8)
            @test field_value(loaded,0.5,0.5,0.5)≈3.0
            @test_throws ArgumentError StructuredField(ascii;text=true,max_samples=7)
            binary=joinpath(dir,"grid.bin")
            open(binary,"w") do io
                write(io,Float64[0,0,0,1,1,1])
                write(io,Int32[2,2,2])
                write(io,Float64[0,3,2,5,1,4,3,6])
            end
            loadedbinary=StructuredField(binary;text=false,max_samples=8)
            @test field_value(loadedbinary,0.5,0.5,0.5)≈3.0
            open(binary,"a") do io; write(io,UInt8(0)); end
            @test_throws ArgumentError StructuredField(binary;text=false,max_samples=8)
        end

        box=BoxField(0,1,0,1,0,1;vin=0.2,vout=1.0)
        rst=RestrictField(box; volumes=[1])
        @test size_at(rst,0.5,0.5,0.5)==0.2                 # no entity → unrestricted
        @test size_at(rst,0.5,0.5,0.5,(3,1))==0.2
        @test size_at(rst,0.5,0.5,0.5,(3,2))==Tessella.SizeField.GMSH_MAX_SIZE
        cf=ConstantField(;vin=0.3,vout=0.9,volumes=[5])
        @test field_value(cf,0,0,0)==Tessella.SizeField.GMSH_MAX_SIZE
        @test size_at(cf,0,0,0,(3,5))==0.3
        @test size_at(cf,0,0,0,(3,6))==0.9
        topology=Dict((3,5)=>[(2,7)],(2,7)=>[(1,9)],(3,6)=>[(1,11)])
        embedded=Dict((3,5)=>[(1,12)])
        cftopo=ConstantField(;vin=0.3,vout=0.9,volumes=[5],
            entity_boundaries=topology,entity_embedded=embedded)
        @test size_at(cftopo,0,0,0,(2,7))==0.3
        @test size_at(cftopo,0,0,0,(1,9))==0.3
        @test size_at(cftopo,0,0,0,(1,12))==0.3
        @test size_at(cftopo,0,0,0,(1,11))==0.9
        embedded_boundaries=Dict((2,7)=>[(1,8)],(1,8)=>[(0,3)],
                                 (1,9)=>[(0,4),(0,5)])
        embedded_children=Dict((2,7)=>[(1,9)])
        embedded_constant=ConstantField(;vin=0.25,vout=1.0,surfaces=[7],
            include_boundary=false,include_embedded=true,
            entity_boundaries=embedded_boundaries,entity_embedded=embedded_children)
        @test size_at(embedded_constant,0,0,0,(1,8))==1.0
        @test size_at(embedded_constant,0,0,0,(0,3))==1.0
        @test size_at(embedded_constant,0,0,0,(1,9))==0.25
        @test size_at(embedded_constant,0,0,0,(0,4))==0.25
        embedded_restrict=RestrictField(ConstantSize(0.25);surfaces=[7],
            include_boundary=false,include_embedded=true,
            entity_boundaries=embedded_boundaries,entity_embedded=embedded_children)
        @test size_at(embedded_restrict,0,0,0,(0,4))==0.25
        @test size_at(embedded_restrict,0,0,0,(0,3))==
              Tessella.SizeField.GMSH_MAX_SIZE
        @test_throws ArgumentError ConstantField(;surfaces=[7],
            entity_embedded=Dict((2,7)=>[(2,9)]))
        nested_constant=ConstantField(;vin=0.25,vout=0.75,volumes=[7])
        nested_restrict=RestrictField(nested_constant;volumes=[7])
        nested_min=BoundedSize(MinSize((nested_restrict,ConstantSize(0.5))))
        @test size_at(nested_min,0,0,0,(3,7))==0.25
        @test size_at(nested_min,0,0,0,(3,8))==0.5
        nested_math=MathEvalField("F1+0.1";fields=Dict(1=>nested_restrict))
        @test size_at(nested_math,0,0,0,(3,7))==0.35
        @test size_at(nested_math,0,0,0,(3,8))>1e21

        ext=ExtendField([(0.0,0.0,0.0)],[0.1]; dist_max=1.0, size_max=0.5, power=1)
        @test size_at(ext,0.0,0.0,0.0)==Tessella.SizeField.GMSH_MAX_SIZE
        @test size_at(ext,0.0,0.0,0.0,(1,1))==Tessella.SizeField.GMSH_MAX_SIZE
        @test size_at(ext,0.0,0.0,0.0,(2,1))==0.1
        @test size_at(ext,0.5,0.0,0.0,(2,1))≈0.5*0.1+0.5*0.5
        @test size_at(ext,1.0,0.0,0.0,(3,1))==0.5
        scaled_ext=ExtendField(;curve_seeds=[(0.0,0.0,0.0)],curve_sizes=[0.2],
            dist_max=1.0,size_max=0.5,global_factor=2.0,
            entity_factors=Dict((2,7)=>2.0))
        @test size_at(scaled_ext,0,0,0,(2,7))==0.05
        @test size_at(scaled_ext,0,0,0,(3,7))==Tessella.SizeField.GMSH_MAX_SIZE

        oct=OctreeField(box,-0.5,1.5,-0.5,1.5,-0.5,1.5; max_level=3)
        @test field_value(oct,0.5,0.5,0.5)==0.2
        @test field_value(oct,2.0,2.0,2.0)==1.0
        lineoct=OctreeField(ConstantSize(0.2),0.0,1.0,0.0,0.0,0.0,0.0;
                            max_level=1)
        @test field_value(lineoct,0.5,0.0,0.0)==0.2
        tinyvarying=FunctionSize((x,y,z)->1e-20*(1+100*abs(x-0.5)))
        tinyoct=OctreeField(tinyvarying,0.0,1.0,0.0,0.0,0.0,0.0;max_level=4)
        @test field_value(tinyoct,0.5,0.0,0.0)≈1.1953125e-20
        @test_throws ArgumentError OctreeField(box,0.0,1.0,0.0,1.0,0.0,1.0;
                                               max_level=33)
        @test_throws ArgumentError OctreeField(box,0.0,1.0,0.0,1.0,0.0,1.0;
                                               max_level=1,max_cells=1)

        pv=PostViewField([0.0 1.0; 0.0 0.0; 0.0 0.0],[0.2,0.8])
        @test pv.values==[0.2,0.8]
        @test field_value(pv,0.0,0.0,0.0)==0.2
        @test field_value(pv,0.9,0.0,0.0)==0.8
        pvexact=PostViewField([0.0 1.0;0.0 0.0;0.0 0.0],[0.2,0.8];
                              use_closest=false)
        @test field_value(pvexact,0.0,0.0,0.0)==0.2
        @test field_value(pvexact,eps(),0.0,0.0)==Tessella.SizeField.GMSH_MAX_SIZE

        # Gmsh 4.15.2 scalar-list oracles (`gmsh.view.probe`) for SL/ST/SS:
        # linear interpolation in containing first-order simplices, followed by
        # closest-node fallback only when requested.
        pvline=PostViewField([0.0 2.0;0.0 0.0;0.0 0.0],[0.2,0.8];
                             lines=[(1,2)],use_closest=false)
        @test field_value(pvline,0.5,0.0,0.0)≈0.35
        @test field_value(pvline,2.0+1e-7,0.0,0.0)≈0.80000003
        @test field_value(pvline,3.0,0.0,0.0)==Tessella.SizeField.GMSH_MAX_SIZE
        @test field_value(pvline,0.5,0.005,0.0)==Tessella.SizeField.GMSH_MAX_SIZE
        pvline_closest=PostViewField([0.0 2.0;0.0 0.0;0.0 0.0],[0.2,0.8];
                                     lines=reshape(Int32[1,2],2,1))
        @test field_value(pvline_closest,3.0,0.0,0.0)==0.8
        @test field_value(pvline_closest,0.5,0.005,0.0)==0.2

        pvtriangle=PostViewField([0.0 1.0 0.0;0.0 0.0 1.0;0.0 0.0 0.0],
            [0.2,0.5,0.8];triangles=[(1,2,3)],use_closest=false)
        @test field_value(pvtriangle,0.25,0.25,0.0)≈0.425
        @test field_value(pvtriangle,0.5,0.5000001,0.0)≈0.65000006
        # List-triangle lookup uses Gmsh's 1%-thickened element box and its
        # dominant-plane reference map, so this nearby off-plane query projects.
        @test field_value(pvtriangle,0.2,0.2,0.005)≈0.38
        @test field_value(pvtriangle,0.2,0.2,0.02)==
              Tessella.SizeField.GMSH_MAX_SIZE
        pvtriangle_closest=PostViewField(
            [0.0 1.0 0.0;0.0 0.0 1.0;0.0 0.0 0.0],[0.2,0.5,0.8];
            triangles=reshape(Int32[1,2,3],3,1))
        @test field_value(pvtriangle_closest,2.0,2.0,0.0)==0.5

        pvtetra=PostViewField(
            [0.0 1.0 0.0 0.0;0.0 0.0 1.0 0.0;0.0 0.0 0.0 1.0],
            [0.2,0.4,0.6,0.8];tetrahedra=[(1,2,3,4)],use_closest=false)
        @test field_value(pvtetra,0.1,0.2,0.3)≈0.48
        @test field_value(pvtetra,1/3,1/3,1/3)≈0.6
        @test field_value(pvtetra,1/3,1/3,1/3+1e-7)≈0.60000006
        @test field_value(pvtetra,2.0,2.0,2.0)==Tessella.SizeField.GMSH_MAX_SIZE
        pvtetra_closest=PostViewField(
            [0.0 1.0 0.0 0.0;0.0 0.0 1.0 0.0;0.0 0.0 0.0 1.0],
            [0.2,0.4,0.6,0.8];tetrahedra=reshape(Int32[1,2,3,4],4,1))
        @test field_value(pvtetra_closest,2.0,2.0,2.0)==0.4

        # Public Gmsh 4.15.2 `view.probe` oracles for first-order SQ/SH/SI/SY.
        # The physical queries below were generated from the listed reference
        # coordinates with Gmsh's own shape functions; the returned values were
        # 0.516, 0.603, 0.695 and 0.645, respectively.
        quadrangle_coords=[0.0 2.0 2.4 -0.2;
                           0.0 0.0 1.5  1.0;
                           0.0 0.0 0.0  0.0]
        pvquadrangle=PostViewField(quadrangle_coords,[0.2,0.4,0.8,1.0];
            quadrangles=[(1,2,3,4)],use_closest=false)
        @test field_value(pvquadrangle,1.256,0.455,0.0)≈0.516
        @test field_value(pvquadrangle,1.256,0.455,1e-7)≈0.516
        @test field_value(pvquadrangle,1.256,0.455,1e-6)==
              Tessella.SizeField.GMSH_MAX_SIZE
        @test field_value(pvquadrangle,2.240000118,0.900000015,0.0)≈
              0.639999998
        affine_origin=[0.1,-0.3,0.7]
        affine_a=[sqrt(2.0),pi/7,-0.4];affine_b=[-0.2,sqrt(3.0),0.6]
        affine_quadrangle=hcat(affine_origin,affine_origin+affine_a,
            affine_origin+affine_a+affine_b,affine_origin+affine_b)
        affine_query=affine_quadrangle*[0.26,0.39,0.21,0.14]
        pvaffine_quadrangle=PostViewField(affine_quadrangle,[0.2,0.4,0.8,1.0];
            quadrangles=[(1,2,3,4)],use_closest=false)
        @test field_value(pvaffine_quadrangle,affine_query...)≈0.516

        hexahedron_coords=[0.0 2.0 2.2 -0.1  0.1 2.1 2.3 0.0;
                           0.0 0.0 1.2  1.0 -0.1 0.1 1.1 1.2;
                           0.0 0.0 0.1  0.0  1.1 1.0 1.4 1.2]
        pvhexahedron=PostViewField(hexahedron_coords,
            [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9];
            hexahedra=[ntuple(identity,8)],use_closest=false)
        @test field_value(pvhexahedron,1.298,0.406,0.8029)≈0.603

        prism_coords=[0.0 2.0 0.0  0.2 2.3 -0.1;
                      0.0 0.0 1.0 -0.1 0.2  1.2;
                      0.0 0.0 0.0  1.2 1.1  1.3]
        pvprism=PostViewField(prism_coords,[0.2,0.3,0.5,0.7,0.9,1.1];
            prisms=[ntuple(identity,6)],use_closest=false)
        @test field_value(pvprism,0.491,0.335,0.847)≈0.695

        pyramid_coords=[-1.0  1.0 1.0 -1.0  0.2;
                        -1.0 -1.0 1.0  1.0 -0.1;
                         0.0  0.0 0.2  0.1  1.4]
        pvpyramid=PostViewField(pyramid_coords,[0.2,0.3,0.5,0.7,1.1];
            pyramids=[ntuple(identity,5)],use_closest=false)
        @test field_value(pvpyramid,0.28,-0.34,0.585)≈0.645
        pvpyramid_closest=PostViewField(pyramid_coords,[0.2,0.3,0.5,0.7,1.1];
            pyramids=[ntuple(identity,5)])
        @test field_value(pvpyramid_closest,10.0,10.0,10.0)==0.5

        # Gmsh interpolates vector components before taking their norm. Its
        # PostView scalar operator returns MAX_LC for tensor views.
        vector_values=[3.0 0.0 0.0 1.0;
                       0.0 4.0 0.0 2.0;
                       0.0 0.0 5.0 2.0]
        pvvector=PostViewField(quadrangle_coords,vector_values;
            quadrangles=[(1,2,3,4)],use_closest=false)
        @test field_value(pvvector,0.584,0.92,0.0)≈2.862236887471056
        @test field_value(pvvector,0.584,0.92,0.0,(2,7))≈2.862236887471056
        @test field_value(pvvector,10.0,10.0,10.0)==
              Tessella.SizeField.GMSH_MAX_SIZE
        vector_line=PostViewField([0.0 2.0;zeros(2,2)],
            [3.0 0.0;0.0 4.0;0.0 0.0];lines=[(1,2)],use_closest=false)
        @test field_value(vector_line,1.0,0.0,0.0)==2.5
        zero_vector=PostViewField(reshape(Float64[0,0,0],3,1),zeros(3,1))
        @test field_value(zero_vector,0.0,0.0,0.0)==
              Tessella.SizeField.GMSH_MAX_SIZE
        @test field_value(PostViewField(reshape(Float64[0,0,0],3,1),zeros(3,1);
                          crop_negative=false),0.0,0.0,0.0)==0.0
        tensor_values=reshape(
            [10(k-1)+(component-1)+0.25 for component in 1:9,k in 1:5],9,5)
        pvtensor=PostViewField(pyramid_coords,tensor_values;
            pyramids=[ntuple(identity,5)],use_closest=false)
        @test field_value(pvtensor,0.0,0.0,0.0)==Tessella.SizeField.GMSH_MAX_SIZE
        @test field_value(pvtensor,100.0,100.0,100.0)==
              Tessella.SizeField.GMSH_MAX_SIZE

        # SS/SH/SI/SY, then ST/SQ, SL and SP are searched in Gmsh's exact
        # list-class precedence order even when their coordinates overlap.
        pvmixed=PostViewField(
            [0.0 1.0 0.0  0.0 1.0;0.0 0.0 1.0  0.0 0.0;zeros(1,5)],
            [0.2,0.5,0.8, 1.2,1.8];lines=[(4,5)],triangles=[(1,2,3)])
        @test field_value(pvmixed,0.5,0.0,0.0)==0.35
        @test field_value(pvmixed,2.0,0.0,0.0)==0.5

        pvcrop=PostViewField([0.0 1.0;0.0 0.0;0.0 0.0],[-1.0,1.0];
                             lines=[(1,2)])
        @test field_value(pvcrop,0.5,0.0,0.0)==Tessella.SizeField.GMSH_MAX_SIZE
        @test field_value(PostViewField([0.0 1.0;0.0 0.0;0.0 0.0],[-1.0,1.0];
                          lines=[(1,2)],crop_negative=false),0.5,0.0,0.0)==0.0
        pvactive=PostViewField([0.0 2.0 100.0;0.0 0.0 0.0;0.0 0.0 0.0],
                               [0.2,0.8,0.01];lines=[(1,2)])
        @test field_value(pvactive,99.0,0.0,0.0)==0.8
        forest_coords=zeros(3,24);forest_values=zeros(24)
        forest_lines=Matrix{Int32}(undef,2,12)
        for j in 1:12
            forest_coords[:,2j-1].=(0.0,3j,0.0)
            forest_coords[:,2j].=(2.0,3j,0.0)
            forest_values[2j-1]=j/100
            forest_values[2j]=j/100+0.2
            forest_lines[:,j].=(2j-1,2j)
        end
        pvforest=PostViewField(forest_coords,forest_values;lines=forest_lines)
        @test field_value(pvforest,1.0,30.0,0.0)==0.2

        pvneg=PostViewField(reshape(Float64[0,0,0],3,1),[-1.0])
        @test field_value(pvneg,0,0,0)==Tessella.SizeField.GMSH_MAX_SIZE
        @test_throws ArgumentError PostViewField(reshape(Float64[0,0,0],3,1),[NaN])
        @test_throws ArgumentError PostViewField(reshape(Float64[0,0,0],3,1),[Inf])
        @test_throws ArgumentError PostViewField([0.0 1.0;0.0 0.0;0.0 0.0],
                                                  [0.2,0.8];lines=[(1,3)])
        @test_throws ArgumentError PostViewField([0.0 0.0;0.0 0.0;0.0 0.0],
                                                  [0.2,0.8];lines=[(1,2)])
        @test_throws ArgumentError PostViewField(
            [0.0 1.0 2.0;0.0 0.0 0.0;0.0 0.0 0.0],[0.2,0.5,0.8];
            triangles=[(1,2,3)])
        @test_throws ArgumentError PostViewField(
            [0.0 1.0 0.0 1.0;0.0 0.0 1.0 1.0;0.0 0.0 0.0 0.0],
            [0.2,0.4,0.6,0.8];tetrahedra=[(1,2,3,4)])
        @test_throws ArgumentError PostViewField(
            [0.0 1.0 1.0 0.0;0.0 0.0 1.0 1.0;0.0 0.0 0.2 0.0],
            [0.2,0.4,0.6,0.8];quadrangles=[(1,2,3,4)])
        @test_throws ArgumentError PostViewField(
            [0.0 1.0 1.0 0.0;0.0 0.0 1.0 1.0;zeros(1,4)],
            [0.2,0.4,0.6,0.8];quadrangles=[(1,2,3,1)])
        @test_throws ArgumentError PostViewField(zeros(3,8),collect(1.0:8.0);
            hexahedra=[ntuple(identity,8)])
        @test_throws ArgumentError PostViewField(zeros(3,6),collect(1.0:6.0);
            prisms=[ntuple(identity,6)])
        @test_throws ArgumentError PostViewField(zeros(3,5),collect(1.0:5.0);
            pyramids=[ntuple(identity,5)])
        @test_throws ArgumentError PostViewField(zeros(3,2),zeros(2,2))
        @test_throws ArgumentError PostViewField(zeros(3,2),zeros(3,3))
        @test_throws ArgumentError PostViewField(zeros(3,2),"unsupported")
        @test_throws ArgumentError PostViewField(reshape(Float64[0,0,0],3,1),
            fill(NaN,3,1))
        @test_throws ArgumentError PostViewField(quadrangle_coords,
            [0.2,0.4,0.8,1.0];points=[1],quadrangles=[(1,2,3,4)],
            max_elements=1)
        @test_throws ArgumentError PostViewField(reshape(Float64[0,0,0],3,1),
            [0.2];points=Int[])
        @test_throws ArgumentError PostViewField(reshape(Float64[0,0,0],3,1),
            [0.2];max_nodes=0)
        @test_throws ArgumentError PostViewField([0.0 1.0;0.0 0.0;0.0 0.0],
            [0.2,0.8];max_elements=1)
        @test_throws ArgumentError PostViewField(reshape(Float64[0,0,0],3,1),
            [0.2];reference_tolerance=-1)
        @test_throws ArgumentError PostViewField(reshape(Float64[0,0,0],3,1),
            [0.2];reference_tolerance=NaN)
        @test_throws ArgumentError PostViewField([0.0 Inf;0.0 0.0;0.0 0.0],
            [0.2,0.8];lines=[(1,2)])
        @test_throws ArgumentError PostViewField([0.0 1.0;0.0 0.0;0.0 0.0],
            [0.2,0.8];lines=[(1,true)])
        @test_throws ArgumentError PostViewField([0.0 1.0;0.0 0.0;0.0 0.0],
            [0.2,0.8];lines=[(1,2)],max_nodes=1)
        farpv=PostViewField([0.0 5e22; 0.0 0.0; 0.0 0.0],[0.2,0.8])
        @test field_value(farpv,1e23,0.0,0.0)==0.8
        overflowpv=PostViewField([-1e308 -prevfloat(1e308); 0.0 0.0; 0.0 0.0],[0.2,0.8])
        @test field_value(overflowpv,1e308,0.0,0.0)==0.8
        wideline=PostViewField([-1e308 1e308;0.0 0.0;0.0 0.0],[0.2,0.8];
                               lines=[(1,2)],use_closest=false)
        @test field_value(wideline,0.0,0.0,0.0)≈0.5
        widetriangle=PostViewField(
            [-1e308 1e308 -1e308;-1e308 -1e308 1e308;0.0 0.0 0.0],
            [0.2,0.5,0.8];triangles=[(1,2,3)],use_closest=false)
        @test field_value(widetriangle,-5e307,-5e307,0.0)≈0.425
        largest=floatmax(Float64)
        widequadrangle=PostViewField(
            [-largest largest largest -largest;-largest -largest largest largest;
             0.0 0.0 0.0 0.0],[0.2,0.4,0.8,1.0];
            quadrangles=[(1,2,3,4)],use_closest=false)
        @test field_value(widequadrangle,0.0,0.0,0.0)≈0.6
        widehexahedron=PostViewField(
            [-largest largest largest -largest -largest largest largest -largest;
             -largest -largest largest largest -largest -largest largest largest;
             -largest -largest -largest -largest largest largest largest largest],
            [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9];
            hexahedra=[ntuple(identity,8)],use_closest=false)
        @test field_value(widehexahedron,0.0,0.0,0.0)≈0.55
        tiny=Float64(2)^-1000
        tinyquadrangle=PostViewField([0.0 tiny tiny 0.0;0.0 0.0 tiny tiny;
                                      0.0 0.0 0.0 0.0],
            [0.2,0.4,0.8,1.0];quadrangles=[(1,2,3,4)],use_closest=false)
        @test field_value(tinyquadrangle,tiny/2,tiny/2,0.0)≈0.6
        tinytetra=PostViewField(
            [0.0 tiny 0.0 0.0;0.0 0.0 tiny 0.0;0.0 0.0 0.0 tiny],
            [0.2,0.4,0.6,0.8];tetrahedra=[(1,2,3,4)],use_closest=false)
        @test field_value(tinytetra,tiny/10,tiny/5,3tiny/10)≈0.48
        robust_vector=PostViewField(reshape(Float64[0,0,0],3,1),
            reshape(fill(floatmax(Float64)/2,3),3,1))
        @test field_value(robust_vector,0,0,0)≈
              hypot(floatmax(Float64)/2,floatmax(Float64)/2,floatmax(Float64)/2)
        overflowing_vector=PostViewField(reshape(Float64[0,0,0],3,1),
            reshape(fill(floatmax(Float64),3),3,1))
        @test_throws ArgumentError field_value(overflowing_vector,0,0,0)
        _nearest_point_checksum(pv.tree,(0.2,0.3,0.4),1)
        # Julia 1.11 can box the single checksum returned through `@allocated`;
        # the bound remains constant over 10,000 nearest-point queries.
        @test (@allocated _nearest_point_checksum(
            pv.tree,(0.2,0.3,0.4),10_000))<=256
        for (postview,point) in ((pv,(0.0,0.0,0.0)),
                                 (pvline,(0.5,0.0,0.0)),
                                 (pvtriangle,(0.25,0.25,0.0)),
                                 (pvtetra,(0.1,0.2,0.3)),
                                 (pvquadrangle,(1.256,0.455,0.0)),
                                 (pvhexahedron,(1.298,0.406,0.8029)),
                                 (pvprism,(0.491,0.335,0.847)),
                                 (pvpyramid,(0.28,-0.34,0.585)),
                                 (pvvector,(0.584,0.92,0.0)),
                                 (pvforest,(1.0,30.0,0.0)))
            _postview_checksum(postview,point,1)
            @test (@allocated _postview_checksum(postview,point,10_000))<=64
        end

        iso1=ConstantSize(0.2); iso2=ConstantSize(0.5)
        @test metric_size(metric_at(MinAnisoField((iso1,iso2)),0,0,0))≈0.2
        @test metric_size(metric_at(IntersectAnisoField((iso1,iso2)),0,0,0))≈0.2
        negative_scalar=MathEvalField("-2")
        for composed in (MinAnisoField((negative_scalar,)),
                         IntersectAnisoField((negative_scalar,)))
            @test metric_eigenvalues(metric_at(composed,0,0,0))==(0.25,0.25,0.25)
        end
        @test_throws ArgumentError metric_at(
            MinAnisoField((MathEvalField("0"),)),0,0,0)
        @test metric_size(metric_at(MinAnisoField(()),0,0,0))≈1e11
        @test metric_eigenvalues(metric_at(IntersectAnisoField(()),0,0,0))==
              (1.0,1.0,1.0)
        minaniso_a=MathEvalAnisoField(m11="1",m22="4",m33="9")
        minaniso_b=MathEvalAnisoField(m11="4",m22="1",m33="16",
                                      m12="1.2",m13="0.5",m23="-0.7")
        minaniso_nonproportional=MinAnisoField((minaniso_a,minaniso_b))
        @test field_value(minaniso_nonproportional,0,0,0)≈
              0.2499015799854048 rtol=2e-14
        @test field_value(MathEvalField("F1";
            fields=Dict(1=>minaniso_nonproportional)),0,0,0)≈
              0.2499015799854048 rtol=2e-14
        @test metric_size(metric_at(minaniso_nonproportional,0,0,0))≈
              0.24105168550141742 rtol=2e-14
        nearly_proportional=MinAnisoField((
            MathEvalAnisoField(m11="1.009",m22="1.005",m33="1"),
            MathEvalAnisoField(m11="1",m22="1",m33="1")))
        @test field_value(nearly_proportional,0,0,0)==1.0
        # Packed isotropic metric M=I/h² has eigenvalues 1/h².
        m=isotropic_metric(0.25)
        @test all(isapprox.(metric_eigenvalues(m),(16.0,16.0,16.0); atol=1e-10))
        @test metric_size(m)==0.25
        anisotropic_max_input=MathEvalAnisoField(m11="25",m22="4",m33="1")
        @test field_value(anisotropic_max_input,0,0,0)==25.0
        @test field_value(MathEvalField("F1";fields=Dict(1=>anisotropic_max_input)),
                          0,0,0)==25.0
        @test field_value(MinSize((anisotropic_max_input,ConstantSize(0.3))),0,0,0)==0.2
        @test field_value(MaxSize((anisotropic_max_input,ConstantSize(0.3))),0,0,0)==1.0

        negative_constant=ConstantField(;vin=-0.2,vout=-0.5,surfaces=[7])
        absolute_constant=MathEvalField("Abs(F1)";fields=Dict(1=>negative_constant))
        @test field_value(negative_constant,0,0,0,(2,7))==-0.2
        @test size_at(absolute_constant,0,0,0,(2,7))==0.2
        @test_throws ArgumentError size_at(negative_constant,0,0,0,(2,7))
        negative_fields=(
            ThresholdField(ConstantSize(0.5);dist_min=0,dist_max=1,
                           size_min=-1,size_max=-2),
            BoxField(0,1,0,1,0,1;vin=-1,vout=-2),
            BallField((0,0,0),1;vin=-1,vout=-2),
            CylinderField((0,0,0),(0,0,1),1;vin=-1,vout=-2),
            FrustumField((0,0,0),(0,0,1);inner_v1=-1,outer_v1=-2,
                         inner_v2=-3,outer_v2=-4))
        @test all(field_value(f,0,0,0)<0 for f in negative_fields)
        @test_throws ArgumentError Metric3(1.0,1.0,-1.0,0.0,0.0,0.0)
        @test_throws ArgumentError Metric3(0.031210066175308172,0.04548846499171803,
            1.0,0.03767888006038275,0.0,0.0)
        exact_spd=Metric3(464.144700333979,0.3315716601763125,1.0,
                          12.405532187366026,0.0,0.0)
        @test minimum(metric_eigenvalues(exact_spd))>0
        cancellation_spd=Metric3(3.152064392711487e12,3.2511546865103975e14,
                                  1.0,3.2012261592312562e13,0.0,0.0)
        cancellation_norm=hypot(cancellation_spd.m12,cancellation_spd.m11)
        cancellation_direction=(cancellation_spd.m12/cancellation_norm,
                                -cancellation_spd.m11/cancellation_norm,0.0)
        @test Tessella.SizeField._metric_displacement(cancellation_spd,
            cancellation_direction...,"test")≈0.00021189863351576526 rtol=1e-12
        tiny=Metric3(2e-20,2e-20,1e-20,1e-20,0.0,0.0)
        @test all(isapprox.(metric_eigenvalues(tiny),(1e-20,1e-20,3e-20);rtol=2e-14))
        @test metric_size(tiny)≈inv(sqrt(3e-20)) rtol=2e-14
        @test metric_size(isotropic_metric(1000.0))≈1000.0
        @test_throws ArgumentError isotropic_metric(1e308)
        extreme_diag=Metric3(1e308,1e-308,1.0,0.0,0.0,0.0)
        @test metric_eigenvalues(extreme_diag)==(1e-308,1.0,1e308)
        fine=isotropic_metric(1e-100)
        @test Tessella.SizeField.intersection_alauzet(isotropic_metric(1e100),fine)==fine
        @test metric_size(metric_at(MinAnisoField((ConstantSize(1e20),)),0,0,0))≈1e11
        @test metric_size(metric_at(IntersectAnisoField((ConstantSize(0.8),
            ConstantSize(0.5),ConstantSize(0.2))),0,0,0))≈0.2

        # Independent LinearAlgebra eigensolver certifies Loewner dominance.
        ma=Metric3(4.0,1.5,2.0,0.4,0.2,0.1)
        mb=Metric3(1.2,5.0,3.0,-0.3,0.1,0.4)
        mc=Tessella.SizeField.intersection_alauzet(ma,mb)
        matrix(m)=Float64[m.m11 m.m12 m.m13; m.m12 m.m22 m.m23;
                            m.m13 m.m23 m.m33]
        scale=max(opnorm(matrix(ma),Inf),opnorm(matrix(mb),Inf))
        @test eigmin(Symmetric(matrix(mc)-matrix(ma)))>=-2e-13*scale
        @test eigmin(Symmetric(matrix(mc)-matrix(mb)))>=-2e-13*scale
        strictbase=Metric3(2.0,2.0,1.0,0.5,0.0,0.0)
        strictscale=1-1e-12
        strictsmall=Metric3(strictscale*strictbase.m11,strictscale*strictbase.m22,
            strictscale*strictbase.m33,strictscale*strictbase.m12,0.0,0.0)
        strictintersection=Tessella.SizeField.intersection_alauzet(strictsmall,strictbase)
        @test eigmin(Symmetric(matrix(strictintersection)-matrix(strictbase)))>=0
        F=prevfloat(floatmax(Float64))
        boundary_a=Metric3(0.625F,0.625F,0.1F,0.375F,0.0,0.0)
        boundary_b=Metric3(0.625F,0.625F,0.1F,-0.375F,0.0,0.0)
        for boundary in (Tessella.SizeField.intersection_alauzet(boundary_a,boundary_b),
                         Tessella.SizeField.intersection_alauzet(boundary_b,boundary_a))
            @test boundary==Metric3(F,F,0.1F,0.0,0.0,0.0)
            @test Tessella.SizeField._metric_exact_dominates(boundary,boundary_a)
            @test Tessella.SizeField._metric_exact_dominates(boundary,boundary_b)
        end
        false_positive_b=Metric3(1.4921023375593005e308,6.573400430653564e307,
            2.354581448184908e307,5.90323140770533e307,0.0,0.0)
        false_positive_c=Metric3(1.7976931348623127e308,1.7976931348623125e308,
            2.354581448184912e307,-3.855613692892678e292,0.0,0.0)
        @test !Tessella.SizeField._metric_dominates(false_positive_c,false_positive_b)
        repaired=Tessella.SizeField.intersection_alauzet(false_positive_c,false_positive_b)
        @test Tessella.SizeField._metric_exact_dominates(repaired,false_positive_c)
        @test Tessella.SizeField._metric_exact_dominates(repaired,false_positive_b)
        guarded_b=Metric3(9.169084656150875e307,3.375470769655791e307,
            3.0467013359016653e307,5.563270372841562e307,0.0,0.0)
        guarded_c=Metric3(9.169215246617614e307,3.3755188447166945e307,
            3.0467447284718497e307,5.563349607585468e307,0.0,0.0)
        @test !Tessella.SizeField._metric_dominates(guarded_c,guarded_b)
        guarded_upper=Tessella.SizeField.intersection_alauzet(guarded_c,guarded_b)
        @test Tessella.SizeField._metric_exact_dominates(guarded_upper,guarded_c)
        @test Tessella.SizeField._metric_exact_dominates(guarded_upper,guarded_b)
        conserve_a=Metric3(2.360718490485931e228,6.146417278634859e230,
            7.8598868218991e231,3.795146938844832e229,-1.359315533190524e230,
            -2.1974630417556248e231)
        conserve_b=Metric3(1.30162270189748e232,5.014110953537154e230,
            4.549433197104851e230,-2.554360706831325e231,
            -2.433144076915384e231,4.775039284740944e230)
        for (first,second) in ((conserve_a,conserve_b),(conserve_b,conserve_a))
            upper=Tessella.SizeField.intersection_conserve_mostaniso(first,second)
            @test Tessella.SizeField._metric_exact_dominates(upper,conserve_a)
            @test Tessella.SizeField._metric_exact_dominates(upper,conserve_b)
        end
        hard_a=Metric3((125/128)*F,(1/2)*F,(1/2)*F,(7/32)*F,0.0,0.0)
        hard_b=Metric3((125/128)*F,(1/2)*F,(1/2)*F,(9/32)*F,0.0,0.0)
        hard_upper=Tessella.SizeField.intersection_alauzet(hard_a,hard_b)
        @test Tessella.SizeField._metric_exact_dominates(hard_upper,hard_a)
        @test Tessella.SizeField._metric_exact_dominates(hard_upper,hard_b)
        impossible_a=Metric3(1.3909615777407798e308,8.478924518090991e307,
            2.0585825161099052e307,-6.215415599727288e307,0.0,0.0)
        impossible_b=Metric3(8.478924518090991e307,1.3909615777407798e308,
            2.0585825161099052e307,6.215415599727288e307,0.0,0.0)
        @test_throws ErrorException Tessella.SizeField.intersection_alauzet(impossible_a,
                                                                            impossible_b)
        _metric_eigenvalue_checksum(mc,1)
        Tessella.SizeField.intersection_alauzet(ma,mb)
        @test (@allocated _metric_eigenvalue_checksum(mc,10_000))<=64
        expected_checksum=10_000*Tessella.SizeField.intersection_alauzet(ma,mb).m11
        @test _sizefield_intersection_checksum(ma,mb,10_000)≈expected_checksum
        _sizefield_intersection_checksum(ma,mb,1)
        @test (@allocated _sizefield_intersection_checksum(ma,mb,10_000))<=64
        mext1=Metric3(1e300,1.0,1.0,0.0,0.0,0.0)
        ea=(1e280+1e270)/2; eb=(1e280-1e270)/2
        mext2=Metric3(1e270,ea,ea,0.0,0.0,eb)
        mext=Tessella.SizeField.intersection_conserve_mostaniso(mext1,mext2)
        quad(m,v)=sum(v .* (matrix(m)*v))
        for v in ([1.0,0.0,0.0],[0.0,inv(sqrt(2.0)),inv(sqrt(2.0))],
                  [0.0,inv(sqrt(2.0)),-inv(sqrt(2.0))])
            @test quad(mext,v)>=quad(mext1,v)*(1-2e-12)
            @test quad(mext,v)>=quad(mext2,v)*(1-2e-12)
        end

        att=AttractorAnisoCurveField([(0.0,0.0,0.0),(1.0,0.0,0.0)],
                                     [(1.0,0.0,0.0),(1.0,0.0,0.0)];
                                     dist_min=0.0,dist_max=1.0,
                                     size_min_tangent=0.2,size_max_tangent=0.8,
                                     size_min_normal=0.1,size_max_normal=0.4)
        # Metric consumers see the normal length; Gmsh's distinct scalar
        # operator (used by F-references) returns max(curve distance, 0.05).
        @test metric_size(metric_at(att,0.0,0.0,0.0))≈0.1
        @test size_at(att,0.0,0.0,0.0)==0.05
        @test field_value(att,0.0,0.2,0.0)==0.2
        @test field_value(MathEvalField("F1";fields=Dict(1=>att)),
                          0.0,0.2,0.0)==0.2
        oblique=AttractorAnisoCurveField([(0.0,0.0,0.0)],[(1.0,1.0,0.0)];
            dist_min=0.0,dist_max=1.0,size_min_tangent=0.5,size_max_tangent=0.5,
            size_min_normal=0.1,size_max_normal=0.1)
        @test all(isapprox.(metric_eigenvalues(metric_at(oblique,0,0,0)),
                            (4.0,50.0,50.0)))
        _attractor_metric_checksum(att,1)
        @test (@allocated _attractor_metric_checksum(att,10_000))<=64

        bl=BoundaryLayerField(dist; hwall=0.05, ratio=1.2, thickness=0.2, hfar=0.5)
        @test size_at(bl,0.0,0.0,0.0)==0.05
        @test size_at(bl,0.1,0.0,0.0)≈0.1*(1.2-1)+0.05
        @test size_at(bl,1.0,0.0,0.0)==Tessella.SizeField.GMSH_MAX_SIZE
        unit_ratio=BoundaryLayerField(dist;hwall=0.1,ratio=1.0,
                                      thickness=1.0,hfar=1.0)
        @test field_value(unit_ratio,0.5,0.0,0.0)==0.1
        zero_ratio=BoundaryLayerField(dist;hwall=0.1,ratio=0.0,
                                      thickness=1.0,hfar=1.0)
        @test field_value(zero_ratio,0.0,0.0,0.0)==0.1
        @test field_value(zero_ratio,0.1,0.0,0.0)==
              Tessella.SizeField.GMSH_MAX_SIZE
        disabled_layer=BoundaryLayerField(dist;hwall=0.1,ratio=1.1,
                                          thickness=-1.0,hfar=1.0)
        @test field_value(disabled_layer,0.0,0.0,0.0)==
              Tessella.SizeField.GMSH_MAX_SIZE
        raw_layer=BoundaryLayerField(dist;hwall=0.1,ratio=0.5,
                                     thickness=1.0,hfar=1.0)
        @test field_value(raw_layer,0.2,0.0,0.0)==0.0
        @test_throws ArgumentError size_at(raw_layer,0.2,0.0,0.0)
        negative_layer=BoundaryLayerField(dist;hwall=-0.2,ratio=1.0,
                                          thickness=1.0,hfar=-0.1)
        @test field_value(negative_layer,0.2,0.0,0.0)==-0.2
        overflow_support=BoundaryLayerField(dist;hwall=0.1,ratio=1e308,
                                            thickness=1e308,hfar=1.0)
        @test field_value(overflow_support,0.0,0.0,0.0)==0.1

        # Regular octahedron inscribed in the unit sphere: curvature ~ 1.
        C=Float64[1 -1 0 0 0 0; 0 0 1 -1 0 0; 0 0 0 0 1 -1]
        F=Int32[1 1 1 1 2 2 2 2; 3 5 4 6 3 6 4 5; 5 4 6 3 6 4 5 3]
        octsurf=Mesh(C; tris=F)
        auto=AutomaticMeshSizeField(octsurf; n_nodes_per_circle=40, hmin=0.01, hmax=2.0)
        h0=size_at(auto,1.0,0.0,0.0)
        @test 2π/40*0.5 < h0 < 2π/40*2.5   # sphere-fit radius near 1
        for radius in (1.0,1e-6,1e-9)
            scaled=Mesh(radius.*C;tris=F)
            autoscaled=AutomaticMeshSizeField(scaled;n_nodes_per_circle=40,
                                               hmin=1e-15,hmax=2.0)
            @test size_at(autoscaled,radius,0.0,0.0)/radius≈2π/40 rtol=1e-12
        end
        flat_auto_mesh=Mesh(Float64[0 1 1 0;0 0 1 1;0 0 0 0];
            tris=Int32[1 1;2 3;3 4])
        flat_auto=AutomaticMeshSizeField(flat_auto_mesh)
        @test all(==(0.05),flat_auto.h)

        helper=joinpath(@__DIR__,"tmp","extproc_size.jl")
        mkpath(dirname(helper))
        write(helper,"""
            while true
                eof(stdin) && break
                x=read(stdin,Float64); y=read(stdin,Float64); z=read(stdin,Float64)
                any(isnan,(x,y,z)) && break
                write(stdout, hypot(x,y,z)+0.05); flush(stdout)
            end
            """)
        extp=ExternalProcessField("julia --startup-file=no $helper")
        @test size_at(extp,3.0,4.0,0.0)≈5.05
        @test isopen(extp)
        close(extp)
        @test !isopen(extp)
        @test_throws ArgumentError size_at(extp,0.0,0.0,0.0)
        # Explicit close waits for the protocol peer, so a clean peer cannot
        # lose its final stdout flush and a non-zero peer exit is observable.
        for _ in 1:3
            clean=ExternalProcessField("julia --startup-file=no $helper")
            @test size_at(clean,1.0,2.0,2.0)≈3.05
            @test isnothing(close(clean))
        end
        failed=ExternalProcessField("julia --startup-file=no $helper; exit 3")
        @test size_at(failed,3.0,4.0,0.0)≈5.05
        @test_throws ArgumentError close(failed)
        nested_leaf=ExternalProcessField("julia --startup-file=no $helper")
        nested=BoundedSize(MathEvalField("F1+0.01";fields=Dict(1=>nested_leaf));
                           size_min=0.001,size_max=10.0)
        @test size_at(nested,0.0,0.0,0.0)≈0.06
        @test isopen(nested)
        @test isnothing(close(nested))
        @test !isopen(nested_leaf) && !isopen(nested)
        @test_throws ArgumentError ExternalProcessField("")
        scalar_helper=joinpath(@__DIR__,"tmp","extproc_scalar.jl")
        write(scalar_helper,"""
            while true
                eof(stdin) && break
                x=read(stdin,Float64); y=read(stdin,Float64); z=read(stdin,Float64)
                any(isnan,(x,y,z)) && break
                write(stdout,x); flush(stdout)
            end
            """)
        scalar_leaf=ExternalProcessField("julia --startup-file=no $scalar_helper")
        @test field_value(scalar_leaf,-0.2,0.0,0.0)==-0.2
        @test_throws ArgumentError size_at(scalar_leaf,-0.2,0.0,0.0)
        scalar_abs=MathEvalField("Abs(F1)";fields=Dict(1=>scalar_leaf))
        @test size_at(scalar_abs,-0.2,0.0,0.0)==0.2
        @test isnothing(close(scalar_leaf))
        mktempdir() do dir
            marker=joinpath(dir,"started")
            slow=joinpath(dir,"slow_extproc.jl")
            write(slow,"""
                marker=ARGS[1]
                x=read(stdin,Float64); y=read(stdin,Float64); z=read(stdin,Float64)
                write(marker,"started")
                sleep(0.5)
                write(stdout,x+y+z+1); flush(stdout)
                while !eof(stdin)
                    x=read(stdin,Float64); y=read(stdin,Float64); z=read(stdin,Float64)
                    any(isnan,(x,y,z)) && break
                    write(stdout,x+y+z+1); flush(stdout)
                end
                """)
            interrupted=ExternalProcessField(
                "julia --startup-file=no $slow $marker")
            task=@async try
                size_at(interrupted,1.0,0.0,0.0)
            catch err
                err
            end
            for _ in 1:500
                isfile(marker) && break
                sleep(0.01)
            end
            @test isfile(marker)
            schedule(task,InterruptException();error=true)
            @test fetch(task) isa InterruptException
            @test !isopen(interrupted)
            @test_throws ArgumentError size_at(interrupted,2.0,0.0,0.0)
        end

        aniso=MathEvalAnisoField(; m11="1/(0.2*0.2)", m22="1/(0.4*0.4)", m33="1/(0.4*0.4)")
        @test metric_size(metric_at(aniso,0,0,0))≈0.2
        @test size_at(aniso,0,0,0)≈25.0
        @test directional_size(aniso,(0,0,0),(1,0,0))≈0.2
        @test directional_size(aniso,(0,0,0),(0,1,0))≈0.4
        @test metric_edge_length(aniso,(0,0,0),(0.2,0,0))≈1.0
        @test metric_edge_length(aniso,(0,0,0),(0,0.4,0))≈1.0

        gp=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),Dict(
            1=>GeoFieldSpec(1,"MathEval",Dict("F"=>"0.1+0.2*x"))),1)
        gf=build_geo_size_field(gp,Dict())
        @test size_at(gf,1.0,0.0,0.0)≈0.3
        bgm_scalar=GeoParams(0.1,10.0,2.0,0,Dict{Tuple{Int,Int},String}(),
            Dict(1=>GeoFieldSpec(1,"MathEval",Dict("F"=>"5"))),1)
        @test size_at(build_geo_size_field(bgm_scalar,Dict();
            model_characteristic_length=1.0),0,0,0)==2.0
        bgm_aniso=GeoParams(0.1,10.0,2.0,0,Dict{Tuple{Int,Int},String}(),
            Dict(1=>GeoFieldSpec(1,"MathEvalAniso",Dict(
                "M11"=>"10000","M22"=>"4","M33"=>"1",
                "M12"=>"0","M13"=>"0","M23"=>"0"))),1)
        bgm_metric=metric_at(build_geo_size_field(bgm_aniso,Dict();
            model_characteristic_length=1.0),0,0,0)
        @test metric_eigenvalues(bgm_metric)==(0.25,1.0,2500.0)
        bgm_nearly_proportional=GeoParams(0.0,10.0,1.0,0,
            Dict{Tuple{Int,Int},String}(),Dict(1=>GeoFieldSpec(1,"MathEvalAniso",
                Dict("M11"=>"1.009","M22"=>"1.005","M33"=>"1",
                     "M12"=>"0","M13"=>"0","M23"=>"0"))),1)
        @test metric_eigenvalues(metric_at(build_geo_size_field(
            bgm_nearly_proportional,Dict();model_characteristic_length=1.0),
            0,0,0))==(1.0,1.0,1.0)
        for (kind,native,point,expected) in (
            ("PostView",pv,(0.0,0.0,0.0),0.2),
            ("AttractorAnisoCurve",att,(0.0,0.0,0.0),0.1),
            ("AutomaticMeshSizeField",auto,(1.0,0.0,0.0),h0))
            contextgp=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),
                Dict(1=>GeoFieldSpec(1,kind,Dict{String,String}())),1)
            @test size_at(build_geo_size_field(contextgp,Dict();
                context_fields=Dict(1=>native),
                model_characteristic_length=Tessella.SizeField.GMSH_MAX_SIZE),
                point...)≈expected
            @test_throws ArgumentError build_geo_size_field(contextgp,Dict())
        end
        extendgp=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),
            Dict(1=>GeoFieldSpec(1,"Extend",Dict("CurvesList"=>"{1}"))),1)
        @test_throws ArgumentError build_geo_size_field(extendgp,Dict())
        built_extend=build_geo_size_field(extendgp,Dict();
            context_fields=Dict(1=>ext),
            model_characteristic_length=Tessella.SizeField.GMSH_MAX_SIZE)
        @test size_at(built_extend,0,0,0,(2,1))≈0.1
        disabled_extend=GeoParams(NaN,NaN,1.0,0,
            Dict{Tuple{Int,Int},String}(),Dict(
                1=>GeoFieldSpec(1,"Extend",Dict(
                    "CurvesList"=>"{1}","DistMax"=>"0"))),1)
        @test size_at(build_geo_size_field(disabled_extend,Dict()),0,0,0)==1.0
        empty_extend_graph=GeoParams(NaN,NaN,1.0,0,
            Dict{Tuple{Int,Int},String}(),Dict(
                1=>GeoFieldSpec(1,"Extend",Dict{String,String}()),
                2=>GeoFieldSpec(2,"MathEval",Dict(
                    "F"=>"0.01+Abs((F1-1e22)/1e22)"))),2)
        @test size_at(build_geo_size_field(empty_extend_graph,Dict()),0,0,0)≈0.01
        invalid_postview=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),
            Dict(1=>GeoFieldSpec(1,"PostView",Dict(
                "ViewIndex"=>"definitely_not_an_int",
                "CropNegativeValues"=>"maybe"))),1)
        @test_throws ArgumentError build_geo_size_field(invalid_postview,Dict();
            context_fields=Dict(1=>pv))
        captured_config=Ref{Any}()
        configured_postview=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),
            Dict(1=>GeoFieldSpec(1,"PostView",Dict("ViewIndex"=>"2",
                "ViewTag"=>"9","CropNegativeValues"=>"0","UseClosest"=>"1"))),1)
        build_geo_size_field(configured_postview,Dict();context_fields=
            (spec,config,entities,params)->(captured_config[]=config;pv))
        @test captured_config[]==(view_index=2,view_tag=9,
            crop_negative_values=false,use_closest=true)
        coerced_postview=GeoParams(NaN,NaN,1.0,0,
            Dict{Tuple{Int,Int},String}(),Dict(
                1=>GeoFieldSpec(1,"PostView",Dict("ViewIndex"=>"2.7",
                    "ViewTag"=>"-1.9","CropNegativeValues"=>"2",
                    "UseClosest"=>"-3"))),1)
        build_geo_size_field(coerced_postview,Dict();context_fields=
            (spec,config,entities,params)->(captured_config[]=config;pv))
        @test captured_config[]==(view_index=2,view_tag=-1,
            crop_negative_values=true,use_closest=true)
        automatic_config=Ref{Any}()
        automatic_defaults=GeoParams(NaN,NaN,1.0,0,
            Dict{Tuple{Int,Int},String}(),
            Dict(1=>GeoFieldSpec(1,"AutomaticMeshSizeField",Dict{String,String}())),1)
        build_geo_size_field(automatic_defaults,Dict();context_fields=
            (spec,config,entities,params)->(automatic_config[]=config;auto),
            model_characteristic_length=Tessella.SizeField.GMSH_MAX_SIZE)
        @test automatic_config[].points_per_gap==0
        @test automatic_config[].points_per_circle==20
        automatic_snapshot=GeoParams(NaN,NaN,1.0,0,
            Dict{Tuple{Int,Int},String}(),Dict(
                1=>GeoFieldSpec(1,"AutomaticMeshSizeField",
                    Dict{String,String}(),String[],37)),1)
        build_geo_size_field(automatic_snapshot,Dict();context_fields=
            (spec,config,entities,params)->(automatic_config[]=config;auto),
            model_characteristic_length=Tessella.SizeField.GMSH_MAX_SIZE)
        @test automatic_config[].points_per_circle==37
        blparams=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),
            Dict(1=>GeoFieldSpec(1,"BoundaryLayer",Dict{String,String}())),0,[1])
        selected_layers=build_geo_boundary_layer_fields(blparams,Dict())
        @test length(selected_layers)==1
        @test field_value(selected_layers[1],0,0,0)==Tessella.SizeField.GMSH_MAX_SIZE
        active_blparams=GeoParams(NaN,NaN,1.0,0,
            Dict{Tuple{Int,Int},String}(),Dict(
                1=>GeoFieldSpec(1,"BoundaryLayer",Dict("CurvesList"=>"{1}"))),0,[1])
        @test build_geo_boundary_layer_fields(active_blparams,Dict();
            context_fields=Dict(1=>bl))==(bl,)
        bl_config=Ref{Any}()
        custom_fan_default=GeoParams(NaN,NaN,1.0,0,
            Dict{Tuple{Int,Int},String}(),active_blparams.fields,0,[1],NaN,30)
        build_geo_boundary_layer_fields(custom_fan_default,Dict();
            context_fields=(spec,config,entities,params)->
                (bl_config[]=config;bl))
        @test bl_config[].fan_elements==30
        raw_blparams=GeoParams(NaN,NaN,1.0,0,
            Dict{Tuple{Int,Int},String}(),Dict(
                1=>GeoFieldSpec(1,"BoundaryLayer",Dict("CurvesList"=>"{1}",
                    "Size"=>"-0.2","Ratio"=>"0.5","SizeFar"=>"-0.1",
                    "Thickness"=>"-1"))),0,[1])
        build_geo_boundary_layer_fields(raw_blparams,Dict();
            context_fields=(spec,config,entities,params)->
                (bl_config[]=config;bl))
        @test (bl_config[].size,bl_config[].ratio,bl_config[].size_far,
               bl_config[].thickness)==(-0.2,0.5,-0.1,-1.0)
        fan_only_bl=GeoParams(NaN,NaN,1.0,0,
            Dict{Tuple{Int,Int},String}(),Dict(
                1=>GeoFieldSpec(1,"BoundaryLayer",Dict(
                    "FanPointsList"=>"{1}"))),1)
        @test size_at(build_geo_size_field(fan_only_bl,Dict()),0,0,0)==1.0
        @test_throws ArgumentError build_geo_size_field(blparams,Dict();
            context_fields=Dict(1=>bl))
        wrong_bl=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),
            Dict(1=>GeoFieldSpec(1,"MathEval",Dict("F"=>"1"))),0,[1])
        @test_throws ArgumentError build_geo_boundary_layer_fields(wrong_bl,Dict())
        paramgp=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),Dict(
            1=>GeoFieldSpec(1,"MathEval",Dict("F"=>"x")),
            2=>GeoFieldSpec(2,"Param",Dict("InField"=>"1","FX"=>"F3",
                                            "FY"=>"0","FZ"=>"0")),
            3=>GeoFieldSpec(3,"MathEval",Dict("F"=>"2"))),2)
        @test size_at(build_geo_size_field(paramgp,Dict();
            model_characteristic_length=10.0),0.0,0.0,0.0)==2.0
        badparam=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),Dict(
            1=>GeoFieldSpec(1,"MathEval",Dict("F"=>"0.2")),
            2=>GeoFieldSpec(2,"Param",Dict("InField"=>"1"))),2)
        @test size_at(build_geo_size_field(badparam,Dict()),0,0,0)==0.2
        defaultmath=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),Dict(
            1=>GeoFieldSpec(1,"MathEval",Dict{String,String}())),1)
        @test size_at(build_geo_size_field(defaultmath,Dict()),0,0,0)==1.0
        defaultaniso=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),Dict(
            1=>GeoFieldSpec(1,"MathEvalAniso",Dict{String,String}())),1)
        @test size_at(build_geo_size_field(defaultaniso,Dict()),0,0,0)==1.0
        defaultthreshold=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),Dict(
            1=>GeoFieldSpec(1,"Threshold",Dict{String,String}())),1)
        @test size_at(build_geo_size_field(defaultthreshold,Dict()),0,0,0)==1.0
        unresolved_threshold=GeoParams(NaN,NaN,1.0,0,
            Dict{Tuple{Int,Int},String}(),Dict(
                1=>GeoFieldSpec(1,"Threshold",Dict{String,String}()),
                2=>GeoFieldSpec(2,"MathEval",Dict(
                    "F"=>"0.01+Abs(F1-1e22)/1e22"))),2)
        @test size_at(build_geo_size_field(unresolved_threshold,Dict()),0,0,0)==0.01
        negative_background=GeoParams(NaN,NaN,1.0,0,
            Dict{Tuple{Int,Int},String}(),
            Dict(1=>GeoFieldSpec(1,"MathEval",Dict("F"=>"-0.2"))),1)
        @test_throws ArgumentError size_at(
            build_geo_size_field(negative_background,Dict()),0,0,0)
        bounded_negative=GeoParams(0.1,NaN,1.0,0,
            Dict{Tuple{Int,Int},String}(),negative_background.fields,1)
        @test size_at(build_geo_size_field(bounded_negative,Dict()),0,0,0)==0.1
        explicit_unresolved_threshold=GeoParams(NaN,NaN,1.0,0,
            Dict{Tuple{Int,Int},String}(),
            Dict(1=>GeoFieldSpec(1,"Threshold",Dict("InField"=>"9"))),1)
        @test size_at(build_geo_size_field(explicit_unresolved_threshold,Dict()),
                      0,0,0)==1.0
        for kind in ("Min","MinAniso","IntersectAniso")
            empty_graph=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),
                Dict(1=>GeoFieldSpec(1,kind,Dict{String,String}())),1)
            empty_field=build_geo_size_field(empty_graph,Dict();
                model_characteristic_length=Tessella.SizeField.GMSH_MAX_SIZE)
            expected=kind=="IntersectAniso" ? 1.0 :
                     kind=="MinAniso" ? 1e11 : Tessella.SizeField.GMSH_MAX_SIZE
            @test size_at(empty_field,0,0,0)≈expected
        end
        empty_max=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),Dict(
            1=>GeoFieldSpec(1,"Max",Dict{String,String}())),1)
        @test_throws ArgumentError size_at(build_geo_size_field(empty_max,Dict()),0,0,0)
        composed_empty_max=GeoParams(NaN,NaN,1.0,0,
            Dict{Tuple{Int,Int},String}(),Dict(
                1=>GeoFieldSpec(1,"Max",Dict{String,String}()),
                2=>GeoFieldSpec(2,"MathEval",Dict(
                    "F"=>"0.1+Abs((F1+1e22)/1e22)"))),2)
        @test size_at(build_geo_size_field(composed_empty_max,Dict()),0,0,0)≈0.1
        for kind in ("Distance","Attractor")
            empty_graph=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),Dict(
                1=>GeoFieldSpec(1,kind,Dict{String,String}()),
                2=>GeoFieldSpec(2,"MathEval",Dict(
                    "F"=>"0.01+Abs((F1-1e22)/1e22)"))),2)
            @test size_at(build_geo_size_field(empty_graph,Dict()),0,0,0)≈0.01
        end
        for options in (Dict{String,String}(),Dict("CommandLine"=>"\"\""))
            empty_external=GeoParams(NaN,NaN,1.0,0,
                Dict{Tuple{Int,Int},String}(),Dict(
                    1=>GeoFieldSpec(1,"ExternalProcess",options),
                    2=>GeoFieldSpec(2,"MathEval",Dict(
                        "F"=>"0.01+Abs((F1-1e22)/1e22)"))),2)
            @test size_at(build_geo_size_field(empty_external,Dict()),0,0,0)≈0.01
        end
        scalar_child=GeoFieldSpec(2,"Box",Dict("VIn"=>"0.05","VOut"=>"0.05"))
        for kind in ("Min","Max","MinAniso","IntersectAniso"),
            fields_list in ("{1,2}","{99,2}")
            filtered=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),Dict(
                1=>GeoFieldSpec(1,kind,Dict("FieldsList"=>fields_list)),
                2=>scalar_child),1)
            @test size_at(build_geo_size_field(filtered,Dict()),0,0,0)≈0.05
        end
        for kind in ("MinAniso","IntersectAniso"),
            fields_list in ("{1,2}","{2,99}","{99,2}")
            scalar_placeholder=GeoParams(NaN,NaN,1.0,0,
                Dict{Tuple{Int,Int},String}(),Dict(
                    1=>GeoFieldSpec(1,kind,Dict("FieldsList"=>fields_list)),
                    2=>GeoFieldSpec(2,"Box",Dict("VIn"=>"2","VOut"=>"2")),
                    3=>GeoFieldSpec(3,"MathEval",Dict("F"=>"F1"))),3)
            @test size_at(build_geo_size_field(scalar_placeholder,Dict();
                model_characteristic_length=10.0),0,0,0)==1.0
        end
        for kind in ("Gradient","Laplacian","Mean","Curvature",
                     "MaxEigenHessian","LonLat","Param","Restrict")
            defaults=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),Dict(
                1=>GeoFieldSpec(1,"MathEval",Dict("F"=>"0.2")),
                2=>GeoFieldSpec(2,kind,Dict{String,String}())),2)
            @test build_geo_size_field(defaults,Dict()) isa AbstractSizeField
            self_default=GeoParams(NaN,NaN,1.0,0,
                Dict{Tuple{Int,Int},String}(),
                Dict(1=>GeoFieldSpec(1,kind,Dict{String,String}())),1)
            @test size_at(build_geo_size_field(self_default,Dict()),0,0,0)==1.0
            unresolved=GeoParams(NaN,NaN,1.0,0,
                Dict{Tuple{Int,Int},String}(),
                Dict(1=>GeoFieldSpec(1,kind,Dict("InField"=>"9"))),1)
            @test size_at(build_geo_size_field(unresolved,Dict()),0,0,0)==1.0
        end
        for (flag,expected) in ((0,1.0),(1,1+asin(0.6)),(2,1.0),(-1,1.0))
            lonlatgp=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),Dict(
                1=>GeoFieldSpec(1,"MathEval",Dict("F"=>"1+y")),
                2=>GeoFieldSpec(2,"LonLat",Dict("InField"=>"1",
                    "FromStereo"=>string(flag),"RadiusStereo"=>"1"))),2)
            @test size_at(build_geo_size_field(lonlatgp,Dict();
                model_characteristic_length=10.0),1,0,0)≈expected
        end
        @test_throws ArgumentError Tessella.SizeField.parse_matheval("f1")
        windows_path=raw"C:\new\tab\grid.bin"
        @test Tessella.SizeField._geo_string_value("\""*windows_path*"\"","test")==
              windows_path
        @test Tessella.SizeField._geo_string_value("'grid.dat'","test")=="grid.dat"
        @test_throws ArgumentError Tessella.SizeField._geo_string_value("'grid.dat\"","test")
        gradgp=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),Dict(
            1=>GeoFieldSpec(1,"MathEval",Dict("F"=>"x")),
            2=>GeoFieldSpec(2,"Gradient",Dict("InField"=>"1","Kind"=>"0"))),2)
        gradbuilt=build_geo_size_field(gradgp,Dict();model_bbox=(0,10,0,0,0,0))
        @test gradbuilt.input.delta≈sqrt(500)/1e4
        pointbuilt=build_geo_size_field(gradgp,Dict();model_bbox=(0,0,0,0,0,0))
        @test pointbuilt.input.delta≈sqrt(8)/1e4
        @test size_at(gradbuilt,2.0,0.0,0.0)≈1.0
        mktempdir() do dir
            path=joinpath(dir,"quoted_matheval.geo")
            write(path,"""
                Field[1] = MathEval;
                Field[1].F = "0.1 + 0.2*x";
                Background Field = 1;
                """)
            parsed=build_geo_size_field(read_geo_params(path),Dict())
            @test size_at(parsed,1.0,0.0,0.0)≈0.3
        end
    end

    @testset "directional anisotropic 1-D/2-D/3-D sizing" begin
        aniso=MathEvalAnisoField(;m11="25",m22="4",m33="4") # hx=.2, hy=hz=.5
        px,_=mesh_segment((0.0,0.0,0.0),(1.0,0.0,0.0),aniso;
                          nsample=100,anisotropic_metric=true)
        py,_=mesh_segment((0.0,0.0,0.0),(0.0,1.0,0.0),aniso;
                          nsample=100,anisotropic_metric=true)
        @test length(px)==6
        @test length(py)==3

        xs=[0.0,1.0,1.0,0.0]; ys=[0.0,0.0,1.0,1.0]
        segs=[(1,2),(2,3),(3,4),(4,1)]
        planar=mesh_planar(xs,ys,segs;field=aniso,min_angle_deg=15.0)
        @test validate(planar).ok
        for (a,b) in unique_edges(planar.tris,planar.tets)
            @test metric_edge_length(aniso,Tuple(planar.coords[:,a]),
                                      Tuple(planar.coords[:,b]))<=1+2e-14
        end

        seed=Mesh(Float64[0 1 0 0; 0 0 1 0; 0 0 0 1];
                  tets=reshape(Int32[1,2,3,4],4,1))
        volume=refine_to_size(seed,aniso)
        @test validate(volume).ok
        for (a,b) in unique_edges(volume.tris,volume.tets)
            @test metric_edge_length(aniso,Tuple(volume.coords[:,a]),
                                      Tuple(volume.coords[:,b]))<=1+2e-14
        end
    end

    @testset "entity-aware 1-D/2-D/3-D sizing" begin
        line_field=ConstantField(;vin=0.25,vout=0.5,volumes=[7])
        fine_line,_=mesh_segment((0.0,0.0,0.0),(1.0,0.0,0.0),line_field;
                                 nsample=100,entity=(3,7))
        coarse_line,_=mesh_segment((0.0,0.0,0.0),(1.0,0.0,0.0),line_field;
                                   nsample=100,entity=(3,8))
        @test length(fine_line)>length(coarse_line)

        xs=[0.0,1.0,1.0,0.0]; ys=[0.0,0.0,1.0,1.0]
        segs=[(1,2),(2,3),(3,4),(4,1)]
        surface_field=ConstantField(;vin=0.35,vout=2.0,surfaces=[7])
        fine_surface=mesh_planar(xs,ys,segs;min_angle_deg=0.0,
                                 field=surface_field,entity=(2,7))
        coarse_surface=mesh_planar(xs,ys,segs;min_angle_deg=0.0,
                                   field=surface_field,entity=(2,8))
        @test validate(fine_surface).ok && validate(coarse_surface).ok
        @test size(fine_surface.tris,2)>size(coarse_surface.tris,2)

        loops=[[(0.0,0.0,0.0),(1.0,0.0,0.0),(1.0,1.0,0.0),(0.0,1.0,0.0)]]
        curve_field=ConstantField(;vin=0.25,vout=10.0,curves=[1],
                                  include_boundary=false)
        classified_boundary=Tessella.MeshSurface.mesh_planar_face(loops,curve_field;
            min_angle_deg=0.0,entity=(2,7),boundary_entities=(1,1))
        unclassified_boundary=Tessella.MeshSurface.mesh_planar_face(loops,curve_field;
            min_angle_deg=0.0,entity=(2,7))
        @test validate(classified_boundary).ok && validate(unclassified_boundary).ok
        @test nnodes(classified_boundary)==16
        @test nnodes(unclassified_boundary)==4
        @test_throws ArgumentError Tessella.MeshSurface.mesh_planar_face(loops,
            curve_field;min_angle_deg=0.0,entity=(2,7),boundary_entities=[[]])
        classified_pslg=mesh_planar(xs,ys,segs;min_angle_deg=0.0,
            field=curve_field,entity=(2,7),
            boundary_entities=[(1,1),(1,2),(1,3),(1,4)])
        @test validate(classified_pslg).ok
        @test count(i->classified_pslg.coords[2,i]==0.0,
                    axes(classified_pslg.coords,2))==5
        @test_throws ArgumentError mesh_planar(xs,ys,segs;field=curve_field,
            entity=(2,7),boundary_entities=Dict(2=>(1,2)))

        cylinder_unclassified=Tessella.MeshSurface.mesh_cylinder_face(
            (0.,0.,0.),(0.,0.,1.),1.0,1.0,curve_field;
            min_angle_deg=0.0,entity=(2,7))
        cylinder_classified=Tessella.MeshSurface.mesh_cylinder_face(
            (0.,0.,0.),(0.,0.,1.),1.0,1.0,curve_field;
            min_angle_deg=0.0,entity=(2,7),boundary_entities=[(1,1),(1,2)])
        @test validate(cylinder_unclassified).ok && validate(cylinder_classified).ok
        @test nnodes(cylinder_classified)>nnodes(cylinder_unclassified)
        @test_throws ArgumentError Tessella.MeshSurface.mesh_cylinder_face(
            (0.,0.,0.),(0.,0.,1.),1.0,1.0,curve_field;
            min_angle_deg=0.0,entity=(2,7),boundary_entities=[(1,1)])
        localized=BallField((1.0,0.0,0.5),0.25;vin=0.1,vout=1.0)
        localized_cylinder=Tessella.MeshSurface.mesh_cylinder_face(
            (0.,0.,0.),(0.,0.,1.),1.0,1.0,localized;min_angle_deg=0.0)
        @test validate(localized_cylinder).ok
        for (a,b) in unique_edges(localized_cylinder.tris,localized_cylinder.tets)
            @test metric_edge_length(localized,Tuple(localized_cylinder.coords[:,a]),
                Tuple(localized_cylinder.coords[:,b]))<=1+256eps(Float64)
        end

        volume_field=ConstantField(;vin=0.7,vout=2.0,volumes=[7])
        seed=Tessella.mesh_box(0,1,0,1,0,1;hmax=1.0)
        fine_volume=refine_to_size(seed,volume_field;entity=(3,7))
        coarse_volume=refine_to_size(seed,volume_field;entity=(3,8))
        @test validate(fine_volume).ok && validate(coarse_volume).ok
        @test size(fine_volume.tets,2)>size(coarse_volume.tets,2)

        coords=Float64[0 1 0 0 3 4 3 3;
                       0 0 1 0 0 0 1 0;
                       0 0 0 1 0 0 0 1]
        cells=Int32[1 5;2 6;3 7;4 8]
        regions=Mesh(coords;tets=cells,tet_tag=Int32[7,8])
        selective=refine_to_size(regions,
            ConstantField(;vin=0.45,vout=10.0,volumes=[7],include_boundary=false);
            entity_resolver=Dict(7=>(3,7),8=>(3,8)))
        @test validate(selective).ok
        @test count(==(Int32(7)),selective.tet_tag)==82
        @test count(==(Int32(8)),selective.tet_tag)==1
        @test_throws ArgumentError refine_to_size(regions,volume_field;
            entity=(3,7),entity_resolver=Dict(7=>(3,7),8=>(3,8)))
        @test_throws ArgumentError refine_to_size(regions,volume_field;
            entity_resolver=Dict(7=>(3,7)))

        surface_seed=Mesh(Float64[0 1 0 0;0 0 1 0;0 0 0 1];
            tris=Int32[1 1 1 2;3 2 4 3;2 4 3 4],
            tets=reshape(Int32[1,2,3,4],4,1),
            tri_tag=Int32[1,2,3,4],tet_tag=Int32[7])
        surface_only=ConstantField(;vin=0.5,vout=10.0,surfaces=[1,2,3,4],
                                   include_boundary=false)
        surface_ignored=refine_to_size(surface_seed,surface_only;entity=(3,7))
        resolver=Dict{Any,Any}((3,7)=>(3,7),
            ((2,tag)=>(2,tag) for tag in 1:4)...)
        surface_applied=refine_to_size(surface_seed,surface_only;
                                       entity_resolver=resolver)
        @test validate(surface_ignored).ok && validate(surface_applied).ok
        @test size(surface_ignored.tets,2)==1
        @test size(surface_applied.tets,2)>1
        for f in axes(surface_applied.tris,2)
            context=(2,Int(surface_applied.tri_tag[f]))
            for (i,j) in ((1,2),(2,3),(3,1))
                a=surface_applied.tris[i,f];b=surface_applied.tris[j,f]
                @test metric_edge_length(surface_only,
                    Tuple(surface_applied.coords[:,a]),Tuple(surface_applied.coords[:,b]);
                    entity=context)<=1+256eps(Float64)
            end
        end
        boundary_only=Mesh(copy(surface_seed.coords);tris=copy(surface_seed.tris),
                           tri_tag=copy(surface_seed.tri_tag))
        filled_resolver=Dict{Any,Any}((3,0)=>(3,7),
            ((2,tag)=>(2,tag) for tag in 1:4)...)
        sized_boundary=mesh_sized(boundary_only;field=surface_only,
                                  entity_resolver=filled_resolver)
        @test validate(sized_boundary).ok
        @test size(sized_boundary.tris,2)>0 && size(sized_boundary.tets,2)>1

        point_only=ConstantField(;vin=0.45,vout=10.0,points=[1],
                                 include_boundary=false)
        point_ignored=refine_to_size(surface_seed,point_only;entity=(3,7))
        point_applied=refine_to_size(surface_seed,point_only;entity=(3,7),
                                     vertex_entities=Dict(1=>(0,1)))
        @test validate(point_ignored).ok && validate(point_applied).ok
        @test size(point_ignored.tets,2)==1
        @test size(point_applied.tets,2)>1
        point_sized=mesh_sized(boundary_only;field=point_only,entity=(3,7),
                               vertex_entities=Dict(1=>(0,1)))
        @test validate(point_sized).ok
        @test size(point_sized.tris,2)>0 && size(point_sized.tets,2)>1
        @test_throws ArgumentError refine_to_size(surface_seed,point_only;
            entity=(3,7),vertex_entities=[(0,1)])
        @test_throws ArgumentError refine_to_size(surface_seed,point_only;
            entity=(3,7),vertex_entities=Dict(1=>(1,1)))
        @test_throws ArgumentError refine_to_size(surface_seed,point_only;
            entity=(3,7),vertex_entities=Dict(1=>(0,0)))
        @test_throws ArgumentError refine_to_size(surface_seed,point_only;
            entity=(3,7),vertex_entities=Dict(1=>(0,big(typemax(Int))+1)))
    end

    @testset "field-driven 1-D/2-D/3-D sizing" begin
        field=MathEvalField("0.2 + 0.3*abs(x)")
        pts,_=mesh_segment((-1.0,0.0,0.0),(1.0,0.0,0.0),field)
        @test length(pts)>=4
        for i in 1:length(pts)-1
            a=pts[i]; b=pts[i+1]
            mid=((a[1]+b[1])/2,(a[2]+b[2])/2,(a[3]+b[3])/2)
            len=hypot(a[1]-b[1],a[2]-b[2],a[3]-b[3])
            target=min(size_at(field,a),size_at(field,b),size_at(field,mid))
            @test len<=target*1.25
        end
        gaps=[abs(pts[i+1][1]-pts[i][1]) for i in 1:length(pts)-1]
        mids=[(pts[i+1][1]+pts[i][1])/2 for i in 1:length(pts)-1]
        near=[gaps[i] for i in eachindex(gaps) if abs(mids[i])<0.25]
        far=[gaps[i] for i in eachindex(gaps) if abs(mids[i])>0.7]
        @test !isempty(near) && !isempty(far)
        @test maximum(near)<minimum(far)

        xs=[0.0,1.0,1.0,0.0]; ys=[0.0,0.0,1.0,1.0]
        segs=[(1,2),(2,3),(3,4),(4,1)]
        planar=mesh_planar(xs,ys,segs; field=BoxField(0,0.5,0,1,0,1;vin=0.15,vout=0.5),
                           min_angle_deg=20.0)
        @test validate(planar).ok
        @test size(planar.tris,2)>0
        localf=BoxField(0,0.5,0,1,0,1;vin=0.15,vout=0.5)
        # Mesh2D's size contract: every interior triangle's longest edge is ≤ the
        # field at the centroid (the kernel sample location).
        for t in 1:size(planar.tris,2)
            i,j,k=Int(planar.tris[1,t]),Int(planar.tris[2,t]),Int(planar.tris[3,t])
            a=Tuple(planar.coords[:,i]); b=Tuple(planar.coords[:,j]); c=Tuple(planar.coords[:,k])
            emax=max(hypot(a[1]-b[1],a[2]-b[2]),hypot(b[1]-c[1],b[2]-c[2]),
                     hypot(c[1]-a[1],c[2]-a[2]))
            cx=(a[1]+b[1]+c[1])/3; cy=(a[2]+b[2]+c[2])/3
            @test emax<=size_at(localf,cx,cy,0.0)
        end

        surf=Tessella.Geometry.box_surface(0.0,1.0,0.0,1.0,0.0,1.0)
        vol=mesh_sized(surf; field=MathEvalField("0.45"))
        @test validate(vol).ok
        @test size(vol.tets,2)>0
        counts=Tessella.tets_per_region(vol)
        @test !isempty(counts) && all(>(0), values(counts))
        for (a,b) in unique_edges(vol.tris,vol.tets)
            p=Tuple(vol.coords[:,a]); q=Tuple(vol.coords[:,b])
            mid=((p[1]+q[1])/2,(p[2]+q[2])/2,(p[3]+q[3])/2)
            len=hypot(p[1]-q[1],p[2]-q[2],p[3]-q[3])
            target=min(size_at(MathEvalField("0.45"),p),
                       size_at(MathEvalField("0.45"),q),
                       size_at(MathEvalField("0.45"),mid))
            @test len<=target
        end
    end
end
