module ConstantMetricBoxChecks
using Test, Tessella, LinearAlgebra

function audit(mesh,matrix,first,last)
    volumes=Float64[];faces=Dict{NTuple{3,Int},Int}();maximum_metric=0.0
    for vertices in eachcol(mesh.tets)
        a,b,c,d=(mesh.coords[:,vertex] for vertex in vertices)
        push!(volumes,det(hcat(b-a,c-a,d-a))/6)
        for i in 1:3,j in i+1:4
            difference=mesh.coords[:,vertices[j]]-mesh.coords[:,vertices[i]]
            maximum_metric=max(maximum_metric,sqrt(dot(difference,matrix*difference)))
        end
        for selected in ((1,2,3),(1,2,4),(1,3,4),(2,3,4))
            face=Tuple(sort(Int[vertices[index] for index in selected]))
            faces[face]=get(faces,face,0)+1
        end
    end
    boundary=[face for (face,count) in faces if count==1]
    edges=Dict{Tuple{Int,Int},Int}()
    for face in boundary,(i,j) in ((1,2),(1,3),(2,3))
        edge=(face[i],face[j]);edges[edge]=get(edges,edge,0)+1
    end
    boundary_nodes=Set(vertex for face in boundary for vertex in face)
    planar=all(face->any(axis->all(mesh.coords[axis,vertex]==first[axis] for vertex in face) ||
        all(mesh.coords[axis,vertex]==last[axis] for vertex in face),1:3),boundary)
    counts=ntuple(axis->length(unique(mesh.coords[axis,:]))-1,3)
    return (;volumes,maximum_metric,counts,planar,boundary,
        boundary_closed=all(==(2),values(edges)) && all(count->count in (1,2),values(faces)),
        euler=length(boundary_nodes)-length(edges)+length(boundary))
end

@testset "constant diagonal metrics change actual box grid directions" begin
    for (metric,matrix,expected) in (
        (Metric3(9,1,1,0,0,0),diagm([9.,1.,1.]),(6,2,2)),
        (Metric3(1,9,1,0,0,0),diagm([1.,9.,1.]),(2,6,2)))
        before=(metric.m11,metric.m22,metric.m33,metric.m12,metric.m13,metric.m23)
        mesh=mesh_box(0,1,0,1,0,1,metric)
        checked=audit(mesh,matrix,(0.,0.,0.),(1.,1.,1.))
        @test checked.counts==expected && size(mesh.tets,2)==144 && size(mesh.coords,2)==63
        @test all(>(0),checked.volumes) && sum(checked.volumes)≈1.0
        @test checked.maximum_metric<=1.0
        @test checked.boundary_closed && checked.euler==2
        @test length(checked.boundary)==4*(expected[1]*expected[2]+expected[1]*expected[3]+expected[2]*expected[3])
        @test checked.planar && all(0<=coordinate<=1 for coordinate in mesh.coords)
        @test validate(mesh).ok
        @test (metric.m11,metric.m22,metric.m33,metric.m12,metric.m13,metric.m23)==before
    end
end

@testset "constant SPD cross terms participate in every rounded grid edge bound" begin
    metric=Metric3(4,1,1,-1,1,.25)
    matrix=[4. -1. 1.;-1. 1. .25;1. .25 1.]
    mesh=mesh_box(2,3,-1,0,4,5,metric)
    checked=audit(mesh,matrix,(2.,-1.,4.),(3.,0.,5.))
    @test checked.counts==(5,3,3) && size(mesh.tets,2)==270
    @test checked.maximum_metric<=1.0
    @test all(>(0),checked.volumes) && sum(checked.volumes)≈1.0
    @test checked.boundary_closed && checked.euler==2 && checked.planar
    @test all(2<=mesh.coords[1,node]<=3 && -1<=mesh.coords[2,node]<=0 && 4<=mesh.coords[3,node]<=5 for node in axes(mesh.coords,2))
    @test validate(mesh).ok
    scaled=mesh_box(0,1,0,1,0,1,Metric3(9,1,1,0,0,0);max_metric_edge=2)
    checked=audit(scaled,diagm([9.,1.,1.]),(0.,0.,0.),(1.,1.,1.))
    @test checked.counts==(3,1,1) && size(scaled.tets,2)==18
    @test checked.maximum_metric<=2.0 && all(>(0),checked.volumes) && sum(checked.volumes)≈1.0
end

@testset "constant metric box limits fail explicitly and physical hmax keeps its contract" begin
    metric=Metric3(9,1,1,0,0,0)
    for target in (0.,-1.,NaN,Inf,true)
        @test_throws ArgumentError mesh_box(0,1,0,1,0,1,metric;max_metric_edge=target)
    end
    @test_throws ArgumentError mesh_box(0,1,0,1,0,1,metric;max_nodes=62)
    @test_throws ArgumentError mesh_box(0,1,0,1,0,1,metric;max_tets=143)
    @test_throws ArgumentError mesh_box(0,1,0,1,0,1,metric;max_nodes=true)
    @test_throws ArgumentError mesh_box(0,1,0,1,0,1,metric;max_tets=1.5)
    @test_throws ArgumentError mesh_box(0,1,0,1,0,1,metric;max_metric_edge=1e-100)
    @test_throws ArgumentError mesh_box(1,0,0,1,0,1,metric)
    @test_throws ArgumentError mesh_box(NaN,1,0,1,0,1,metric)
    @test_throws ArgumentError mesh_box(false,1,0,1,0,1,metric)
    @test_throws ArgumentError mesh_box(1e16,1e16+2,0,1,0,1,isotropic_metric(1))
    physical=mesh_box(0,1,0,1,0,1;hmax=1)
    checked=audit(physical,Matrix{Float64}(I,3,3),(0.,0.,0.),(1.,1.,1.))
    @test checked.counts==(2,2,2) && size(physical.tets,2)==48
    @test checked.maximum_metric<=1.0 && sum(checked.volumes)≈1.0 && checked.boundary_closed
end
end
