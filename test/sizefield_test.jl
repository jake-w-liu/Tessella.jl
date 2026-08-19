# Gmsh-compatible scalar/size fields and field-driven tetrahedral refinement.

using Test
using Tessella
using Tessella.Mesh1D: mesh_segment
using Tessella.MeshTypes: Mesh, nnodes, ntets, unique_edges, validate, node, tet_volume
using Tessella.Geometry: box_surface
using Tessella.IO: read_geo_params, GeoParams, GeoFieldSpec

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
        field_value(accelerated,(0.1,0.2,0.3))
        @test (@allocated field_value(accelerated,(0.1,0.2,0.3)))==0
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
        @test_throws ArgumentError BoxField(1,0,0,1,0,1;vin=1.0)
        @test_throws ArgumentError MinSize(AbstractField[])
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
        @test_throws ArgumentError BallField((0,0,0),-1;vin=1)
        @test_throws ArgumentError BallField(0,1;vin=1)
        @test_throws ArgumentError CylinderField((0,0,0),(0,0,0),1;vin=1)

        frustum=FrustumField((0.0,0.0,0.0),(0.0,0.0,2.0);
            inner_r1=0.0,outer_r1=1.0,inner_r2=0.0,outer_r2=1.0,
            inner_v1=0.1,outer_v1=1.0,inner_v2=0.2,outer_v2=0.8)
        @test size_at(frustum,0.0,0.0,1.0)≈0.15
        @test size_at(frustum,0.5,0.0,1.0)≈0.525
        @test size_at(frustum,1.0,0.0,1.0)≈0.9
        @test size_at(frustum,1.01,0.0,1.0)==Tessella.SizeField.GMSH_MAX_SIZE
        @test size_at(frustum,0.0,0.0,-0.01)==Tessella.SizeField.GMSH_MAX_SIZE
        @test_throws ArgumentError FrustumField((0,0,0),(0,0,0))
        @test_throws ArgumentError FrustumField((0,0,0),(0,0,1);outer_r1=0)
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
        @test field isa BoundedSize
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
        @test_throws ArgumentError build_geo_size_field(rawdistance,Dict((0,"p")=>point))

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
        conflicting=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),Dict(
            1=>GeoFieldSpec(1,"Frustum",Dict("InnerR1"=>"0","R1_inner"=>"0"))),1)
        @test_throws ArgumentError build_geo_size_field(conflicting,Dict())
    end
end
