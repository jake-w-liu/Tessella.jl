using Test
using Tessella
using Tessella.MeshTypes: Mesh, boundary_edges, mesh_crc, nnodes, node, nsegs,
                          ntris, triangle_area, validate

if !isdefined(Tessella,:Transfinite)
    Base.include(Tessella,joinpath(@__DIR__,"..","src","Transfinite.jl"))
end
using Tessella.Transfinite: mesh_transfinite_patch

function _transfinite_rectangle(L::Int,H::Int;
                               xmin=0.0,xmax=Float64(L),
                               ymin=0.0,ymax=Float64(H))
    L>0&&H>0 || throw(ArgumentError("positive cell counts required"))
    bottom=[(xmin+(xmax-xmin)*i/L,ymin,0.0) for i in 0:L]
    right=[(xmax,ymin+(ymax-ymin)*j/H,0.0) for j in 0:H]
    top=[(xmax-(xmax-xmin)*i/L,ymax,0.0) for i in 0:L]
    left=[(xmin,ymax-(ymax-ymin)*j/H,0.0) for j in 0:H]
    return bottom,right,top,left
end

function _surface_area(mesh)
    sum(triangle_area(node(mesh,mesh.tris[1,t]),node(mesh,mesh.tris[2,t]),
                      node(mesh,mesh.tris[3,t])) for t in 1:ntris(mesh);init=0.0)
end

function _polygon_area(sides)
    ring=NTuple{3,Float64}[]
    for side in sides
        append!(ring,@view side[1:end-1])
    end
    area=0.0
    for i in eachindex(ring)
        p=ring[i];q=ring[mod1(i+1,length(ring))]
        area+=p[1]*q[2]-p[2]*q[1]
    end
    return abs(area)/2
end

function _canonical_triangles(mesh)
    result=NTuple{3,Int32}[]
    for t in axes(mesh.tris,2)
        values=sort(mesh.tris[:,t])
        push!(result,(values[1],values[2],values[3]))
    end
    sort!(result)
end

function _expected_triangles(L,H,arrangement)
    id(i,j)=Int32(i+1+j*(L+1))
    result=NTuple{3,Int32}[]
    for i in 0:L-1,j in 0:H-1
        v1=id(i,j);v2=id(i+1,j);v3=id(i+1,j+1);v4=id(i,j+1)
        right=arrangement===:right ||
              (arrangement===:alternate_right&&isodd(i+j)) ||
              (arrangement===:alternate_left&&iseven(i+j))
        first,second=right ? ((v1,v2,v3),(v3,v4,v1)) :
                             ((v1,v2,v4),(v4,v2,v3))
        for triangle in (first,second)
            values=sort(collect(triangle))
            push!(result,(values[1],values[2],values[3]))
        end
    end
    sort!(result)
end

function _edge_set(matrix)
    result=Set{NTuple{2,Int32}}()
    for i in axes(matrix,2)
        a=matrix[1,i];b=matrix[2,i]
        push!(result,a<b ? (a,b) : (b,a))
    end
    result
end

struct _CountOnlySide <: AbstractVector{NTuple{3,Float64}}
    count::Int
end
Base.size(side::_CountOnlySide)=(side.count,)
Base.getindex(::_CountOnlySide,::Int)=error("resource counts were not checked first")

@noinline function _transfinite_allocated(sides)
    GC.gc()
    return @allocated mesh_transfinite_patch(sides...;arrangement=:alternate_left)
end

@testset "four-sided planar transfinite patches" begin
    @testset "Gmsh arrangements, deterministic CRCs, and physical tags" begin
        sides=_transfinite_rectangle(3,2)
        expected_crc=Dict(
            :left=>"bc5c8cf1ba3bcbd22315c3095990658c9b95efe33ba8f7dcee99b9cc6dadf9e4",
            :right=>"b968747720c166a0a8f25a787528277a5265f8c187db70f7c73903354f05e151",
            :alternate_left=>"4265a3d7e815fb0510d77e6b6bda1b96f776c18bfa8b75fc03fa0a6a0df13f2d",
            :alternate_right=>"8c7c0070649022cbc5dcd81ae8abf0898d091e352ca8fb86b7e8a20e262dc8de")
        for arrangement in (:left,:right,:alternate_left,:alternate_right)
            mesh=mesh_transfinite_patch(sides...;arrangement=arrangement,
                                        face_tag=21,side_tags=(11,12,13,14))
            @test validate(mesh).ok
            @test (nnodes(mesh),nsegs(mesh),ntris(mesh))==(12,10,12)
            @test mesh_crc(mesh).bbox==((0.0,0.0,0.0),(3.0,2.0,0.0))
            @test mesh_crc(mesh).sha==expected_crc[arrangement]
            @test _canonical_triangles(mesh)==_expected_triangles(3,2,arrangement)
            @test mesh.tri_tag==fill(Int32(21),12)
            @test mesh.seg_tag==Int32[11,11,11,12,12,13,13,13,14,14]
            boundary,maxincidence=boundary_edges(mesh.tris)
            @test maxincidence==2
            @test Set(boundary)==_edge_set(mesh.segs)
            @test _surface_area(mesh)==6.0
            @test mesh_crc(mesh)==mesh_crc(mesh_transfinite_patch(
                sides...;arrangement=arrangement,face_tag=21,
                side_tags=(11,12,13,14)))
        end
    end

    @testset "average-chord Coons interpolation and boundary conservation" begin
        bottom=[(0.,0.,0.),(0.7,-0.25,0.),(2.1,-0.55,0.),(3.2,-0.2,0.),(4.,0.,0.)]
        right=[(4.,0.,0.),(4.35,0.8,0.),(4.2,2.1,0.),(4.,3.,0.)]
        top=[(4.,3.,0.),(3.1,3.55,0.),(2.0,3.35,0.),(0.8,3.15,0.),(0.,3.,0.)]
        left=[(0.,3.,0.),(-0.3,2.15,0.),(-0.5,0.9,0.),(0.,0.,0.)]
        sides=(bottom,right,top,left)
        mesh=mesh_transfinite_patch(sides...;arrangement=:alternate_left)
        @test (nnodes(mesh),nsegs(mesh),ntris(mesh))==(20,14,24)
        @test validate(mesh).ok
        @test _surface_area(mesh)≈_polygon_area(sides) atol=64eps(Float64)

        # Boundary nodes are copied bit-for-bit, in deterministic row-major order.
        @test [node(mesh,i) for i in 1:5]==bottom
        @test [node(mesh,5+5j) for j in 0:3]==right
        @test [node(mesh,i+15) for i in 1:5]==reverse(top)
        @test [node(mesh,1+5j) for j in 0:3]==reverse(left)

        # Independent direct evaluation of the documented average-chord formula.
        top_grid=reverse(top);left_grid=reverse(left)
        chord(a,b)=hypot(a[1]-b[1],a[2]-b[2],a[3]-b[3])
        ustep=[0.5(chord(bottom[i+1],bottom[i])+chord(top_grid[i+1],top_grid[i]))
               for i in 1:4]
        vstep=[0.5(chord(right[j+1],right[j])+chord(left_grid[j+1],left_grid[j]))
               for j in 1:3]
        u=sum(ustep[1:2])/sum(ustep);v=vstep[1]/sum(vstep)
        c1=bottom[1];c2=bottom[end];c3=top_grid[end];c4=top_grid[1]
        expected=ntuple(3) do d
            (1-u)*left_grid[2][d]+u*right[2][d]+
            (1-v)*bottom[3][d]+v*top_grid[3][d]-
            ((1-u)*(1-v)*c1[d]+u*(1-v)*c2[d]+u*v*c3[d]+(1-u)*v*c4[d])
        end
        @test all(isapprox(node(mesh,8)[d],expected[d];
                           atol=16eps(Float64),rtol=16eps(Float64)) for d in 1:3)
    end

    @testset "tilted and clockwise patches" begin
        base=_transfinite_rectangle(4,3;xmax=2.0,ymax=1.5)
        ex=(inv(sqrt(2.0)),inv(sqrt(2.0)),0.0)
        ey=(-inv(sqrt(6.0)),inv(sqrt(6.0)),2inv(sqrt(6.0)))
        origin=(1.0,-2.0,3.0)
        transform(p)=(origin[1]+p[1]*ex[1]+p[2]*ey[1],
                      origin[2]+p[1]*ex[2]+p[2]*ey[2],
                      origin[3]+p[1]*ex[3]+p[2]*ey[3])
        tilted=map(side->transform.(side),base)
        mesh=mesh_transfinite_patch(tilted...;arrangement=:right)
        @test validate(mesh).ok
        @test _surface_area(mesh)≈3.0 atol=256eps(Float64)

        clockwise=([(0.,0.,0.),(0.,1.,0.),(0.,2.,0.)],
                   [(0.,2.,0.),(1.,2.,0.),(2.,2.,0.),(3.,2.,0.)],
                   [(3.,2.,0.),(3.,1.,0.),(3.,0.,0.)],
                   [(3.,0.,0.),(2.,0.,0.),(1.,0.,0.),(0.,0.,0.)])
        reversed_mesh=mesh_transfinite_patch(clockwise...;
                                              arrangement=:alternate_right)
        @test validate(reversed_mesh).ok
        @test _surface_area(reversed_mesh)==6.0

        long_patch=mesh_transfinite_patch(_transfinite_rectangle(2048,1)...)
        @test (nnodes(long_patch),ntris(long_patch))==(4098,4096)
        @test validate(long_patch).ok
    end

    @testset "validated blockers and resource limits" begin
        sides=_transfinite_rectangle(3,2)
        @test_throws ArgumentError mesh_transfinite_patch(
            [(0.,0.,0.)],sides[2],sides[3],sides[4])
        @test_throws ArgumentError mesh_transfinite_patch(
            sides[1],sides[2],[(3.,2.,0.),(0.,2.,0.)],sides[4])
        mismatched=copy(sides[2]);insert!(mismatched,2,(3.,0.5,0.))
        @test_throws ArgumentError mesh_transfinite_patch(
            sides[1],mismatched,sides[3],sides[4])
        broken=copy(sides[2]);broken[1]=(3.,nextfloat(0.0),0.)
        @test_throws ArgumentError mesh_transfinite_patch(
            sides[1],broken,sides[3],sides[4])
        nonfinite=copy(sides[1]);nonfinite[2]=(NaN,0.,0.)
        @test_throws ArgumentError mesh_transfinite_patch(
            nonfinite,sides[2],sides[3],sides[4])
        extra=Any[(point...,point==(1.,0.,0.) ? NaN : 0.0)
                  for point in sides[1]]
        @test_throws ArgumentError mesh_transfinite_patch(
            extra,sides[2],sides[3],sides[4])
        duplicate=copy(sides[1]);duplicate[2]=duplicate[1]
        @test_throws ArgumentError mesh_transfinite_patch(
            duplicate,sides[2],sides[3],sides[4])
        nonplanar=copy(sides[1]);nonplanar[2]=(1.,0.,1e-6)
        @test_throws ArgumentError mesh_transfinite_patch(
            nonplanar,sides[2],sides[3],sides[4])

        bowtie=([(0.,0.,0.),(1.,1.,0.)],
                [(1.,1.,0.),(0.,1.,0.)],
                [(0.,1.,0.),(1.,0.,0.)],
                [(1.,0.,0.),(0.,0.,0.)])
        @test_throws ArgumentError mesh_transfinite_patch(bowtie...)
        backtrack=([(0.,0.,0.),(1.,0.,0.),(0.5,0.,0.),(2.,0.,0.)],
                   [(2.,0.,0.),(2.,1.,0.)],
                   [(2.,1.,0.),(1.5,1.,0.),(1.,1.,0.),(0.,1.,0.)],
                   [(0.,1.,0.),(0.,0.,0.)])
        @test_throws ArgumentError mesh_transfinite_patch(backtrack...)
        folded=([(0.,0.,0.),(1.138676204796158,-0.8998793245683743,0.),
                 (1.,0.,0.)],
                [(1.,0.,0.),(0.9753676924206318,-0.6302141619964883,0.),
                 (1.,1.,0.)],
                [(1.,1.,0.),(-0.14986440109197496,1.9513598146977102,0.),
                 (0.,1.,0.)],
                [(0.,1.,0.),(-0.9877714621697844,-0.56016161339637,0.),
                 (0.,0.,0.)])
        @test_throws ArgumentError mesh_transfinite_patch(folded...)

        @test_throws ArgumentError mesh_transfinite_patch(sides...;arrangement=:alternate)
        @test_throws ArgumentError mesh_transfinite_patch(sides...;arrangement="Left")
        @test_throws ArgumentError mesh_transfinite_patch(sides...;face_tag=true)
        @test_throws ArgumentError mesh_transfinite_patch(sides...;face_tag=-1)
        @test_throws ArgumentError mesh_transfinite_patch(
            sides...;face_tag=big(typemax(Int32))+1)
        @test_throws ArgumentError mesh_transfinite_patch(sides...;side_tags=(1,2,3))
        @test_throws ArgumentError mesh_transfinite_patch(sides...;side_tags=(1,2,true,4))
        @test_throws ArgumentError mesh_transfinite_patch(sides...;side_tags=(1,2,-1,4))
        @test_throws ArgumentError mesh_transfinite_patch(sides...;max_nodes=true)
        @test_throws ArgumentError mesh_transfinite_patch(sides...;max_triangles=false)
        @test_throws ArgumentError mesh_transfinite_patch(sides...;max_nodes=-1)
        @test_throws ArgumentError mesh_transfinite_patch(
            sides...;max_nodes=big(typemax(Int32))+1)
        @test_throws TypeError mesh_transfinite_patch(sides...;max_nodes=12.0)
        @test_throws ArgumentError mesh_transfinite_patch(sides...;max_nodes=11)
        @test_throws ArgumentError mesh_transfinite_patch(sides...;max_triangles=11)
        bounded=mesh_transfinite_patch(sides...;max_nodes=12,max_triangles=12)
        @test (nnodes(bounded),ntris(bounded))==(12,12)

        huge=_CountOnlySide(typemax(Int))
        @test_throws ArgumentError mesh_transfinite_patch(huge,huge,huge,huge)
        wide=_CountOnlySide(100_000)
        @test_throws ArgumentError mesh_transfinite_patch(wide,wide,wide,wide)

        maximum=Float64(floatmax(Float64))
        overflowing=([(maximum,0.,0.),(-maximum,0.,0.)],
                     [(-maximum,0.,0.),(-maximum,1.,0.)],
                     [(-maximum,1.,0.),(maximum,1.,0.)],
                     [(maximum,1.,0.),(maximum,0.,0.)])
        @test_throws ArgumentError mesh_transfinite_patch(overflowing...)
    end

    @testset "allocation growth remains linear in output size" begin
        small_sides=_transfinite_rectangle(64,64)
        large_sides=_transfinite_rectangle(128,64)
        mesh_transfinite_patch(small_sides...;arrangement=:alternate_left)
        mesh_transfinite_patch(large_sides...;arrangement=:alternate_left)
        small=_transfinite_allocated(small_sides)
        large=_transfinite_allocated(large_sides)
        @test small>0
        @test large>small
        @test large<=2.25small+262_144
        @info "transfinite allocation ratchet" small_bytes=small large_bytes=large
    end
end
