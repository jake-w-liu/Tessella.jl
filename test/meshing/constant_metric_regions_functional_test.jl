module ConstantMetricRegionsChecks
using Test, Tessella, LinearAlgebra

function audit(mesh,matrix)
    volumes=Float64[];bytag=Dict{Int32,Float64}()
    faces=Dict{NTuple{3,Int},Vector{Int}}();maximum_metric=0.0
    for (cell,vertices) in enumerate(eachcol(mesh.tets))
        a,b,c,d=(mesh.coords[:,vertex] for vertex in vertices)
        volume=det(hcat(b-a,c-a,d-a))/6
        push!(volumes,volume)
        tag=mesh.tet_tag[cell];bytag[tag]=get(bytag,tag,0.0)+volume
        for i in 1:3,j in i+1:4
            difference=mesh.coords[:,vertices[j]]-mesh.coords[:,vertices[i]]
            maximum_metric=max(maximum_metric,sqrt(dot(difference,matrix*difference)))
        end
        for selected in ((1,2,3),(1,2,4),(1,3,4),(2,3,4))
            face=Tuple(sort(Int[vertices[index] for index in selected]))
            push!(get!(faces,face,Int[]),cell)
        end
    end
    boundary=[face for (face,cells) in faces if length(cells)==1]
    edges=Dict{Tuple{Int,Int},Int}()
    for face in boundary,(i,j) in ((1,2),(1,3),(2,3))
        edge=(face[i],face[j]);edges[edge]=get(edges,edge,0)+1
    end
    boundary_nodes=Set(vertex for face in boundary for vertex in face)
    return (;volumes,bytag,maximum_metric,faces,boundary,
        counts=ntuple(axis->length(unique(mesh.coords[axis,:]))-1,3),
        closed=all(==(2),values(edges)) && all(cells->length(cells) in (1,2),values(faces)),
        euler=length(boundary_nodes)-length(edges)+length(boundary))
end

@testset "constant metric shared grid preserves adjacent material interfaces" begin
    regions=[BoxRegion(0,1,0,1,0,1,3),BoxRegion(1,2,0,1,0,1,7)]
    before=[(r.x0,r.x1,r.y0,r.y1,r.z0,r.z1,r.tag) for r in regions]
    metric=Metric3(9,1,1,0,0,0)
    mesh=mesh_box_regions(regions,metric)
    checked=audit(mesh,diagm([9.,1.,1.]))
    interface=[(face,cells) for (face,cells) in checked.faces
               if all(mesh.coords[1,vertex]==1.0 for vertex in face)]
    @test checked.counts==(12,2,2) && size(mesh.tets,2)==288 && size(mesh.coords,2)==117
    @test Set(mesh.tet_tag)==Set(Int32[3,7]) && all(>(0),checked.volumes)
    @test checked.bytag[3]≈1.0 && checked.bytag[7]≈1.0
    @test checked.maximum_metric<=1.0
    @test length(interface)==8 && all(length(cells)==2 &&
        Set(mesh.tet_tag[cells])==Set(Int32[3,7]) for (_,cells) in interface)
    @test checked.closed && checked.euler==2 && validate(mesh).ok
    @test all(0<=mesh.coords[1,node]<=2 && 0<=mesh.coords[2,node]<=1 &&
              0<=mesh.coords[3,node]<=1 for node in axes(mesh.coords,2))
    @test [(r.x0,r.x1,r.y0,r.y1,r.z0,r.z1,r.tag) for r in regions]==before
    swapped=mesh_box_regions(regions,Metric3(1,9,1,0,0,0))
    checked_swap=audit(swapped,diagm([1.,9.,1.]))
    @test checked_swap.counts==(4,6,2) && size(swapped.tets,2)==288
    @test checked_swap.maximum_metric<=1 && all(>(0),checked_swap.volumes) &&
        checked_swap.bytag[3]≈1 && checked_swap.bytag[7]≈1
end

@testset "constant metric Region priority carves a closed literal-plane void" begin
    regions=[BoxRegion(0,1,0,1,0,1,3),BoxRegion(.25,.75,.25,.75,.25,.75,0)]
    mesh=mesh_box_regions(regions,Metric3(9,1,1,0,0,0))
    checked=audit(mesh,diagm([9.,1.,1.]))
    cavity=[face for face in checked.boundary if all(
        .25<=mesh.coords[axis,vertex]<=.75 for axis in 1:3,vertex in face)]
    @test Set(mesh.tet_tag)==Set(Int32[3]) && all(>(0),checked.volumes)
    @test checked.bytag[3]≈.875 && sum(checked.volumes)≈.875
    @test checked.maximum_metric<=1.0
    @test checked.closed && checked.euler==4 && validate(mesh).ok
    @test length(cavity)==28 && all(face->any(axis->all(
        mesh.coords[axis,vertex]==.25 for vertex in face) || all(
        mesh.coords[axis,vertex]==.75 for vertex in face),1:3),cavity)
    @test all(value->value in mesh.coords,.25:.5:.75) &&
        all(vertices->!all(.25<sum(mesh.coords[axis,vertices])/4<.75 for axis in 1:3),eachcol(mesh.tets))
    filled=mesh_box_regions(reverse(regions),Metric3(9,1,1,0,0,0))
    @test sum(audit(filled,diagm([9.,1.,1.])).volumes)≈1.0
end

@testset "constant SPD cross terms certify the actual shared Region grid" begin
    regions=[BoxRegion(2,2.5,-1,-.5,4,4.5,4),BoxRegion(2.5,3,-1,-.5,4,4.5,9)]
    metric=Metric3(4,1,1,-1,1,.25)
    matrix=[4. -1. 1.;-1. 1. .25;1. .25 1.]
    mesh=mesh_box_regions(regions,metric)
    checked=audit(mesh,matrix)
    @test checked.maximum_metric<=1.0
    @test all(>(0),checked.volumes) && checked.bytag[4]≈.125 && checked.bytag[9]≈.125
    @test checked.closed && checked.euler==2 && validate(mesh).ok
    @test all(2<=mesh.coords[1,node]<=3 && -1<=mesh.coords[2,node]<=-.5 &&
              4<=mesh.coords[3,node]<=4.5 for node in axes(mesh.coords,2))
    @test any(face->length(checked.faces[face])==2 &&
        all(mesh.coords[1,vertex]==2.5 for vertex in face),keys(checked.faces))
end

@testset "constant metric Region limits and physical hmax semantics" begin
    regions=[BoxRegion(0,1,0,1,0,1,3),BoxRegion(1,2,0,1,0,1,7)]
    metric=Metric3(9,1,1,0,0,0)
    for target in (0.,-1.,NaN,Inf,true)
        @test_throws ArgumentError mesh_box_regions(regions,metric;max_metric_edge=target)
    end
    @test_throws ArgumentError mesh_box_regions(regions,metric;max_nodes=116)
    @test_throws ArgumentError mesh_box_regions(regions,metric;max_tets=287)
    @test_throws ArgumentError mesh_box_regions(regions,metric;max_nodes=true)
    @test_throws ArgumentError mesh_box_regions(regions,metric;max_tets=1.5)
    @test_throws ArgumentError mesh_box_regions(regions,metric;max_metric_edge=1e-100)
    @test_throws ArgumentError mesh_box_regions(BoxRegion[],metric)
    @test_throws ArgumentError mesh_box_regions([BoxRegion(0,1,0,1,0,1,0)],metric)
    @test_throws ArgumentError mesh_box_regions(
        [BoxRegion(1e16,1e16+2,0,1,0,1,3)],isotropic_metric(1))
    @test_throws ArgumentError mesh_box_regions(
        [BoxRegion(0,1,0,1,0,1,3),BoxRegion(100,101,100,101,100,101,7)],
        isotropic_metric(1);max_nodes=100)
    signed_zero=mesh_box_regions([BoxRegion(-0.,1,0,1,0,1,3),
        BoxRegion(0.,.5,0,1,0,1,7)],isotropic_metric(1))
    zero_checked=audit(signed_zero,Matrix{Float64}(I,3,3))
    @test zero_checked.bytag[3]≈.5 && zero_checked.bytag[7]≈.5 && validate(signed_zero).ok
    physical=mesh_box_regions(regions;hmax=1)
    checked=audit(physical,Matrix{Float64}(I,3,3))
    @test checked.maximum_metric<=1.0 && all(>(0),checked.volumes) &&
        checked.bytag[3]≈1 && checked.bytag[7]≈1 && checked.closed
end
end
