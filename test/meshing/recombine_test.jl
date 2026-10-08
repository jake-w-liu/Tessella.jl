using Test
using Tessella

module RecombineQualityRangeRegressionTests
using Test,Tessella

function geometry_oracle(coords)
    xy=[(Rational{BigInt}(coords[1,i]),Rational{BigInt}(coords[2,i])) for i in 1:4]
    edges=[(xy[mod1(i+1,4)][1]-xy[i][1],xy[mod1(i+1,4)][2]-xy[i][2]) for i in 1:4]
    squared=[edge[1]^2+edge[2]^2 for edge in edges]
    corners=[edges[i][1]*edges[mod1(i+1,4)][2]-edges[i][2]*edges[mod1(i+1,4)][1] for i in 1:4]
    twice_area=sum(xy[i][1]*xy[mod1(i+1,4)][2]-xy[i][2]*xy[mod1(i+1,4)][1] for i in 1:4)
    @test twice_area>0 && all(>(0),corners) && all(>(0),squared)
    @test all(==(coords[3,1]),coords[3,:])
    quality=setprecision(BigFloat,256) do
        lengths=sqrt.(BigFloat.(squared))
        sines=[BigFloat(corners[i])/(lengths[i]*lengths[mod1(i+1,4)]) for i in 1:4]
        # Exact constant-z convex geometry has positive parallel triangle
        # normals. Edge ratios and corner sines use exact represented inputs.
        Float64(min(minimum(lengths)/maximum(lengths),minimum(sines),BigFloat(1)))
    end
    return twice_area/2,quality
end

function exact_area(mesh)
    coords=Rational{BigInt}.(mesh.coords)
    cells=mesh isa Mesh ? (mesh.tris,) :
        (block.nodes for block in mesh.blocks if block.msh in (2,3))
    area=zero(Rational{BigInt})
    for nodes in cells,cell in axes(nodes,2),i in axes(nodes,1)
        a=nodes[i,cell];b=nodes[mod1(i+1,size(nodes,1)),cell]
        area+=(coords[1,a]*coords[2,b]-coords[1,b]*coords[2,a])/2
    end
    return area
end

function fixtures()
    h=1e-200;s=nextfloat(0.0)
    values=[
        Float64[0 1 1 0;0 0 h h;0 0 0 0],
        Float64[0 1 1 0;0 0 h h;1e100 1e100 1e100 1e100],
        Float64[0 1 1 0;0 0 h h;1e200 1e200 1e200 1e200],
        Float64[0 1 1 0;0 0 1 1;1e200 1e200 1e200 1e200],
        Float64[0 h h 0;0 0 h 1;0 0 0 0],
        Float64[0 h h 0;0 0 h h;1 1 1 1],
        Float64[1 2 2 1;h h 2h 2h;0 0 0 0],
        Float64[h 2h 2h h;-h -h 0 0;1 1 1 1],
        Float64[0 1 1 0;0 0 s s;0 0 0 0],
        Float64[1e12 1e12+1 1e12+1 1e12;1e12 1e12 1e12+2 1e12+2;0 0 0 0]]
    for steps in (1,2,4,16,100)
        right=1.0+steps*eps(1.0)
        push!(values,Float64[1 right right 1;0 0 1 1;3 3 3 3])
    end
    return values
end

@noinline function quality_sum(coords,repetitions)
    nodes=(Int32(1),Int32(2),Int32(3),Int32(4))
    result=0.0
    for _ in 1:repetitions
        result+=Tessella.Recombine._quad_quality(coords,nodes)
    end
    return result
end

@testset "Represented geometry survives quality underflow and cancellation" begin
    segs=Int32[1 2 3 4;2 3 4 1]
    tris=Int32[1 1;2 3;3 4]
    nodes=(Int32(1),Int32(2),Int32(3),Int32(4))
    for coords in fixtures()
        area,quality=geometry_oracle(coords)
        @test quality>0
        @test Tessella.Recombine._quad_quality(coords,nodes)==quality
        mesh=Mesh(coords;segs=segs,tris=tris,seg_tag=Int32[2,2,3,3],tri_tag=Int32[7,7])
        snapshot=(copy(mesh.coords),copy(mesh.segs),copy(mesh.tris),
            copy(mesh.seg_tag),copy(mesh.tri_tag))
        before=mesh_crc(mesh)
        for algorithm in (:greedy,:blossom)
            result=recombine_triangles(mesh;algorithm=algorithm,min_quality=quality,
                full_quad=algorithm===:blossom)
            @test Tessella.Elements.validate(result).ok
            @test [block.msh for block in result.blocks]==[1,3]
            @test result.blocks[2].nodes==reshape(Int32[1,2,3,4],4,1)
            @test result.blocks[2].tags==Int32[7]
            @test result.blocks[1].nodes==segs && result.blocks[1].tags==mesh.seg_tag
            @test exact_area(result)==exact_area(mesh)==area
            @test result.coords!==mesh.coords && result.blocks[1].nodes!==mesh.segs
            @test result.blocks[1].tags!==mesh.seg_tag
            @test mesh_crc(mesh)==before
            result.coords[1,1]=17.0
            result.blocks[1].tags[1]=Int32(19)
            @test (mesh.coords,mesh.segs,mesh.tris,mesh.seg_tag,mesh.tri_tag)==snapshot
        end
        if quality<1
            rejected=recombine_triangles(mesh;min_quality=nextfloat(quality))
            @test [block.msh for block in rejected.blocks]==[1,2]
            @test rejected.blocks[2].nodes==tris
            @test_throws ArgumentError recombine_triangles(mesh;algorithm=:blossom,
                full_quad=true,min_quality=nextfloat(quality))
        end
        protected=recombine_triangles(mesh;protected_edges=Set([(Int32(1),Int32(3))]))
        @test [block.msh for block in protected.blocks]==[1,2]
        @test exact_area(protected)==area
        quality_sum(coords,1000)
        @test minimum(@allocated(quality_sum(coords,1000)) for _ in 1:3)==0
    end
    # This represented rectangle has width eps(1). Divide-before-subtract at
    # constant z=3 formerly inflated its score by 50%, accepting this threshold.
    cancellation=Mesh(Float64[1 nextfloat(1.) nextfloat(1.) 1;0 0 1 1;3 3 3 3];tris=tris)
    @test only(recombine_triangles(cancellation;min_quality=2.8e-16,
        preserve_segments=false).blocks).msh==2
    h=1e-200
    folded=Float64[0 1 1 0;0 0 h h;0 2h 0 2h]
    @test Tessella.Recombine._quad_quality(folded,nodes)==0
end
end # module RecombineQualityRangeRegressionTests

module RecombineZeroRoundedQualityRegressionTests
using Test,Tessella
using ..RecombineQualityRangeRegressionTests: geometry_oracle,exact_area,quality_sum

function fixtures()
    values=Matrix{Float64}[]
    for (width,height) in ((1e200,1e-200),(1e300,1e-300),
            (floatmax(Float64),nextfloat(0.0))),z in (0.0,1.0,1e200)
        push!(values,Float64[0 width width 0;0 0 height height;z z z z])
    end
    for width in (1e150,1e200,1e300)
        height=1/width
        push!(values,Float64[0 width 2width width;0 0 height height;0 0 0 0])
    end
    values
end

@testset "Positive represented geometry with a zero-rounded quality score" begin
    nodes=(Int32(1),Int32(2),Int32(3),Int32(4))
    segs=Int32[1 2 3 4;2 3 4 1];tris=Int32[1 1;2 3;3 4]
    for coords in fixtures()
        area,quality=geometry_oracle(coords)
        @test area>0 && quality>=0
        @test Tessella.Recombine._quad_quality(coords,nodes)==quality
        mesh=Mesh(coords;segs,tris,seg_tag=Int32[2,2,3,3],tri_tag=Int32[7,7])
        @test validate(mesh).ok
        snapshot=(copy(mesh.coords),copy(mesh.segs),copy(mesh.tris),
            copy(mesh.seg_tag),copy(mesh.tri_tag))
        before=mesh_crc(mesh)
        for algorithm in (:greedy,:blossom)
            result=recombine_triangles(mesh;algorithm,min_quality=quality,
                full_quad=algorithm===:blossom)
            @test Tessella.Elements.validate(result).ok
            @test [block.msh for block in result.blocks]==[1,3]
            @test result.blocks[2].nodes==reshape(Int32[1,2,3,4],4,1)
            @test result.blocks[2].tags==Int32[7]
            @test result.blocks[1].nodes==segs && result.blocks[1].tags==mesh.seg_tag
            @test exact_area(result)==exact_area(mesh)==area
            @test result.coords!==mesh.coords && result.blocks[1].nodes!==mesh.segs
            @test result.blocks[1].tags!==mesh.seg_tag
            @test mesh_crc(mesh)==before
            result.coords[1,1]=17.0;result.blocks[1].tags[1]=Int32(19)
            @test (mesh.coords,mesh.segs,mesh.tris,mesh.seg_tag,mesh.tri_tag)==snapshot
            rejected=recombine_triangles(mesh;algorithm,min_quality=nextfloat(quality))
            @test [block.msh for block in rejected.blocks]==[1,2]
            @test rejected.blocks[2].nodes==tris
        end
        @test_throws ArgumentError recombine_triangles(mesh;algorithm=:blossom,
            full_quad=true,min_quality=nextfloat(quality))
        protected=recombine_triangles(mesh;protected_edges=Set([(Int32(1),Int32(3))]))
        @test [block.msh for block in protected.blocks]==[1,2]
        @test exact_area(protected)==area
        # Axis-aligned extreme cells need no exact-arithmetic allocation.
        if coords[1,4]==coords[1,1]
            quality_sum(coords,1000)
            @test minimum(@allocated(quality_sum(coords,1000)) for _ in 1:3)==0
        end
    end
    for width in (1e200,1e300)
        height=1/width
        folded=Float64[0 width width 0;0 0 height height;0 2height 0 2height]
        @test Tessella.Recombine._quad_quality(folded,nodes)==0
        @test only(recombine_triangles(Mesh(folded;tris);
            preserve_segments=false).blocks).msh==2
        collinear=Float64[0 width 2width 3width;0 0 0 0;0 0 0 0]
        @test Tessella.Recombine._quad_quality(collinear,nodes)==0
    end
end
end # module RecombineZeroRoundedQualityRegressionTests

function _recombine_square(;width=1.0,tags=Int32[7,7],reverse_second=false)
    coords=Float64[0 width width 0;
                   0 0     1     1;
                   0 0     0     0]
    tris=reverse_second ? Int32[1 1;2 4;3 3] : Int32[1 1;2 3;3 4]
    segs=Int32[1 2 3 4;2 3 4 1]
    return Mesh(coords;segs=segs,tris=tris,seg_tag=Int32[2,2,3,3],tri_tag=tags)
end

module RecombineAngleRegressionTests
using Test,Tessella

function pentagon()
    Mesh(Float64[0 1 1 0 -1;0 0 1 1 .5;0 0 0 0 0];
        tris=Int32[1 1 1;2 3 4;3 4 5],tri_tag=fill(Int32(7),3),
        segs=Int32[1 2 3 4 5;2 3 4 5 1],seg_tag=fill(Int32(9),5))
end

function exact_area(mesh)
    coords=Rational{BigInt}.(mesh.coords)
    area=zero(Rational{BigInt})
    cells=mesh isa Mesh ? (mesh.tris,) :
        (block.nodes for block in mesh.blocks if block.msh in (2,3))
    for nodes in cells,cell in axes(nodes,2)
        for i in axes(nodes,1)
            a=nodes[i,cell];b=nodes[mod1(i+1,size(nodes,1)),cell]
            area+=(coords[1,a]*coords[2,b]-coords[1,b]*coords[2,a])/2
        end
    end
    return area
end

@testset "Pinned odd-count Blossom greedy fallback" begin
    # Gmsh 4.15.2 warns that Blossom cannot run on three triangles. The
    # angle .5 greedy pass selects the nonsquare pair (1,3,4,5); above 1
    # it selects the square (1,2,3,4). Pair identities come from primary replay.
    mesh=pentagon();snapshot=(copy(mesh.coords),copy(mesh.tris),copy(mesh.tri_tag))
    for (angle,triangle,quad) in ((.5,Int32[1,2,3],Int32[1,3,4,5]),
                                 (1.,Int32[1,2,3],Int32[1,3,4,5]),
                                 (1.01,Int32[1,4,5],Int32[1,2,3,4]))
        result=recombine_triangles(mesh;algorithm=:blossom,recombine_angle=angle)
        @test [b.msh for b in result.blocks]==[1,2,3]
        @test result.blocks[2].nodes==reshape(triangle,3,1)
        @test result.blocks[3].nodes==reshape(quad,4,1)
        @test result.blocks[2].tags==result.blocks[3].tags==Int32[7]
        @test result.blocks[1].nodes==mesh.segs && result.blocks[1].tags==mesh.seg_tag
        @test exact_area(result)==exact_area(mesh)==3//2
        @test Tessella.Elements.validate(result).ok
        @test (mesh.coords,mesh.tris,mesh.tri_tag)==snapshot
        @test_throws ArgumentError recombine_triangles(mesh;algorithm=:blossom,
            recombine_angle=angle,full_quad=true)
    end
    protected=Set([(Int32(1),Int32(4))])
    retained=recombine_triangles(mesh;algorithm=:blossom,recombine_angle=.5,
        protected_edges=protected)
    @test [b.msh for b in retained.blocks]==[1,2]
    @test retained.blocks[2].nodes==mesh.tris
    @test exact_area(retained)==3//2
    # The no-angle public core continues to provide maximum-cardinality
    # matching independently of Gmsh's explicit-angle fallback behavior.
    matched=recombine_triangles(mesh;algorithm=:blossom)
    @test matched.blocks[3].nodes==reshape(Int32[1,2,3,4],4,1)
    square=Mesh(mesh.coords[:,1:4];tris=mesh.tris[:,1:2])
    for angle in (0.,.5,1.)
        result=recombine_triangles(square;algorithm=:blossom,
            recombine_angle=angle,full_quad=true,preserve_segments=false)
        @test only(result.blocks).msh==3
        @test only(result.blocks).nodes==reshape(Int32[1,2,3,4],4,1)
    end
    for angle in (-1.,90.01,NaN,Inf,true,".5")
        @test_throws ArgumentError recombine_triangles(mesh;algorithm=:blossom,
            recombine_angle=angle)
    end
end

@testset "Public recombine forwards the odd-count angle and edge exclusions" begin
    api=Tessella.API
    for (protect,expected_types,triangle,quad) in
            ((false,Int32[2,3],UInt64[1,2,3],UInt64[1,3,4,5]),
             (true,Int32[2],UInt64[1,2,3,1,3,4,1,4,5],UInt64[]))
        try
            mesh=pentagon();api.initialize();api.option("Mesh.Renumber",0)
            api.option("Mesh.RecombinationAlgorithm",1)
            protect && api.model.add_discrete_entity(1,11)
            api.model.add_discrete_entity(2,1,protect ? [11] : Int[])
            api.mesh.add_nodes(2,1,UInt64.(1:5),vec(mesh.coords))
            api.mesh.add_elements_by_type(1,2,[101,102,103],UInt64.(vec(mesh.tris)))
            protect && api.mesh.add_elements_by_type(11,1,[104],[1,4])
            before=api.mesh.get_nodes()
            api.mesh.set_recombine(2,1,.5);api.mesh.recombine()
            types,ids,nodes=api.mesh.get_elements(2,1)
            @test types==expected_types
            @test nodes[1]==triangle
            @test api.mesh.get_nodes()==before
            if protect
                @test ids==[UInt64[101,102,103]]
                @test api.mesh.get_elements(1,11)[3]==[UInt64[1,4]]
            else
                @test ids[1]==UInt64[101]
                @test Set(nodes[2])==Set(quad)
                points=[api.mesh.get_node(node)[1] for node in nodes[2]]
                @test sum(points[i][1]*points[mod1(i+1,4)][2]-
                    points[mod1(i+1,4)][1]*points[i][2] for i in 1:4)/2==1
            end
        finally
            api.finalize()
        end
    end
end

@noinline function measure_allocated(coords,nodes,count)
    f=Tessella.Recombine._gmsh_recombine_pair_measure
    total=0.0
    for _ in 1:count
        total+=f(coords,nodes)
    end
    return total
end

@testset "Fixed corner-angle measure uses bounded workspace" begin
    nodes=(Int32(1),Int32(2),Int32(3),Int32(4))
    # Installed Gmsh 4.15.2 getElementQualities(...,"gamma") exposes eta.
    # Its raw cross-product norm underflows/overflows at the extreme scales;
    # this is a compatibility measure, distinct from geometric quad quality.
    for (scale,square_eta,trapezoid_eta) in
            ((1e-300,0.,0.),(1e-150,0.,0.),(1.,1.,.5),
             (1e150,1.,1.),(1e300,1.,.5))
        square=Float64[0 1 1 0;0 0 1 1;0 0 0 0].*scale
        trapezoid=Float64[0 2 1 0;0 0 1 1;0 0 0 0].*scale
        @test Tessella.Recombine._gmsh_recombine_pair_measure(square,nodes)==square_eta
        @test Tessella.Recombine._gmsh_recombine_pair_measure(trapezoid,nodes)==trapezoid_eta
        measure_allocated(square,nodes,1000)
        bytes=minimum(@allocated(measure_allocated(square,nodes,1000)) for _ in 1:5)
        @test bytes<=64
    end
end

@testset "Raw corner eta and stable geometric quality agree at admission" begin
    nodes=(Int32(1),Int32(2),Int32(3),Int32(4))
    triangles=Int32[1 1;2 3;3 4]
    # Independent installed-primary recombination at all three thresholds.
    for (coords,eta,quality) in
            ((Float64[0 1 1 0;0 0 1 1;1e200 1e200 1e200 1e200],1.,1.),
             (Float64[0 1 1 0;0 0 1e-200 1e-200;0 0 0 0],0.,1e-200),
             (Float64[0 1e-200 1e-200 0;0 0 1e-200 1e-200;1 1 1 1],0.,1.))
        @test Tessella.Recombine._gmsh_recombine_pair_measure(coords,nodes)==eta
        @test Tessella.Recombine._quad_quality(coords,nodes)≈quality rtol=8eps()
        mesh=Mesh(coords;tris=triangles,tri_tag=Int32[7,7])
        snapshot=(copy(mesh.coords),copy(mesh.tris),copy(mesh.tri_tag))
        for angle in (.5,1.,1.01)
            result=recombine_triangles(mesh;recombine_angle=angle,preserve_segments=false)
            block=only(result.blocks)
            @test block.msh==(eta<angle ? 3 : 2)
            @test block.nodes==(eta<angle ? reshape(Int32[1,2,3,4],4,1) : triangles)
            @test Tessella.Elements.validate(result).ok
            @test result.coords==coords
        end
        @test (mesh.coords,mesh.tris,mesh.tri_tag)==snapshot
    end
end
end # module RecombineAngleRegressionTests

@static if Sys.iswindows() && Sys.ARCH===:x86_64
    @testset "Pinned raw eta bits at strict native angle boundaries" begin
        # Literal public getElementQualities(...,"gamma") observations from
        # Gmsh 4.15.2. Each score is an admission threshold, not a tolerance.
        fixtures=(
            (Float64[0 2e-80 1e-80 0;0 0 1e-80 1e-80;0 0 0 0],0x3fdffff8917cd1c2),
            (Float64[0 1.1028256440530324 .954072810820371 0;
                     0 0 1.1686172353169644 1.1686172353169644;0 0 0 0],0x3fed6bb62243421b),
            (Float64[0 3.4520631498890326 1.2609112869906363 0;
                     0 0 .49343752331195523 .49343752331195523;0 0 0 0],0x3fc20cab0ea2b51c))
        nodes=(Int32(1),Int32(2),Int32(3),Int32(4))
        for (coords,bits) in fixtures
            threshold=reinterpret(Float64,bits)
            @test reinterpret(UInt64,Tessella.Recombine._gmsh_recombine_pair_measure(coords,nodes))==bits
            mesh=Mesh(coords;tris=Int32[1 1;2 3;3 4])
            for angle in (prevfloat(threshold),threshold,nextfloat(threshold))
                result=recombine_triangles(mesh;recombine_angle=angle,preserve_segments=false)
                @test only(result.blocks).msh==(angle>threshold ? 3 : 2)
            end
        end
    end
end

function _recombine_grid(n::Int)
    n>0 || throw(ArgumentError("grid extent must be positive"))
    side=n+1;coords=Matrix{Float64}(undef,3,side^2)
    node(i,j)=Int32(i+1+j*side)
    @inbounds for j in 0:n,i in 0:n
        id=Int(node(i,j));coords[1,id]=i;coords[2,id]=j;coords[3,id]=0
    end
    tris=Matrix{Int32}(undef,3,2n^2);tags=Vector{Int32}(undef,2n^2)
    cell=0
    @inbounds for j in 0:n-1,i in 0:n-1
        a=node(i,j);b=node(i+1,j);c=node(i+1,j+1);d=node(i,j+1)
        cell+=1
        tris[1,2cell-1]=a;tris[2,2cell-1]=b;tris[3,2cell-1]=c
        tris[1,2cell]=a;tris[2,2cell]=c;tris[3,2cell]=d
        tags[2cell-1]=Int32(1);tags[2cell]=Int32(1)
    end
    return Mesh(coords;tris=tris,tri_tag=tags)
end

function _brute_matching_size(adjacency)
    n=length(adjacency);memo=Dict{UInt,Int}()
    function search(mask::UInt)
        haskey(memo,mask) && return memo[mask]
        vertex=findfirst(i->iszero(mask&(UInt(1)<<(i-1))),1:n)
        vertex===nothing && return 0
        occupied=mask|(UInt(1)<<(vertex-1))
        best=search(occupied)
        for neighbour in adjacency[vertex]
            iszero(mask&(UInt(1)<<(neighbour-1))) || continue
            best=max(best,1+search(occupied|(UInt(1)<<(neighbour-1))))
        end
        memo[mask]=best
        return best
    end
    return search(UInt(0))
end

function _exhaustive_matching_mismatch()
    cases=0
    for n in 0:6
        edges=[(i,j) for i in 1:n-1 for j in i+1:n]
        for bits in UInt(0):(UInt(1)<<length(edges))-1
            adjacency=[Int[] for _ in 1:n]
            for (bit,(i,j)) in pairs(edges)
                iszero(bits&(UInt(1)<<(bit-1))) && continue
                push!(adjacency[i],j);push!(adjacency[j],i)
            end
            mate=Tessella.Recombine._edmonds_matching(n,adjacency)
            got=count(i->mate[i]>i,1:n)
            expected=_brute_matching_size(adjacency)
            cases+=1
            got==expected || return (;cases,n,bits,got,expected,adjacency,mate)
        end
    end
    return nothing
end

@noinline function _recombine_allocated(n::Int)
    mesh=_recombine_grid(n)
    recombine_triangles(mesh)
    GC.gc()
    return @allocated recombine_triangles(mesh)
end

@noinline function _recombine_angle_grid_allocated(mesh,angle)
    recombine_triangles(mesh;recombine_angle=angle)
    return minimum(@allocated(recombine_triangles(mesh;recombine_angle=angle)) for _ in 1:5)
end

@testset "deterministic triangle-to-quadrangle recombination" begin
    @testset "square oracle and metadata preservation" begin
        mesh=_recombine_square()
        names=Dict((1,2)=>"horizontal",(1,3)=>"vertical",(2,7)=>"face")
        result=recombine_triangles(mesh;physical_names=names)
        @test Tessella.Elements.validate(result).ok
        @test [block.msh for block in result.blocks]==[1,3]
        @test result.blocks[1].nodes==mesh.segs
        @test result.blocks[1].tags==mesh.seg_tag
        @test result.blocks[2].nodes==reshape(Int32[1,2,3,4],4,1)
        @test result.blocks[2].tags==Int32[7]
        @test result.physical_names==
              Dict((1,2)=>"horizontal",(1,3)=>"vertical",(2,7)=>"face")
        @test result.physical_names!==names
        @test result.coords!==mesh.coords
        @test result.blocks[1].nodes!==mesh.segs
        @test result.blocks[1].tags!==mesh.seg_tag
        @test Tessella.Elements.mixed_crc(result)==
              Tessella.Elements.mixed_crc(recombine_triangles(mesh;
                  physical_names=result.physical_names))
        names[(2,7)]="mutated"
        @test result.physical_names[(2,7)]=="face"

        without_segments=recombine_triangles(mesh;preserve_segments=false)
        @test [block.msh for block in without_segments.blocks]==[3]
        @test without_segments.blocks[1].nodes==reshape(Int32[1,2,3,4],4,1)
    end

    @testset "selection contracts" begin
        different_tags=recombine_triangles(
            _recombine_square(tags=Int32[7,8]);preserve_segments=false)
        @test [block.msh for block in different_tags.blocks]==[2]
        @test different_tags.blocks[1].nodes==_recombine_square().tris
        @test different_tags.blocks[1].tags==Int32[7,8]

        inconsistent=recombine_triangles(
            _recombine_square(reverse_second=true);preserve_segments=false)
        @test [block.msh for block in inconsistent.blocks]==[2]
        @test size(inconsistent.blocks[1].nodes,2)==2

        elongated=_recombine_square(width=10.0)
        accepted=recombine_triangles(elongated;min_quality=0.09,
                                     preserve_segments=false)
        rejected=recombine_triangles(elongated;min_quality=0.11,
                                     preserve_segments=false)
        @test [block.msh for block in accepted.blocks]==[3]
        @test [block.msh for block in rejected.blocks]==[2]

        concave=Mesh(Float64[0 1 .2 0;0 0 .2 1;0 0 0 0];
                     tris=Int32[1 1;2 3;3 4],tri_tag=Int32[4,4])
        concave_result=recombine_triangles(concave;preserve_segments=false)
        @test [block.msh for block in concave_result.blocks]==[2]

        vertical=Mesh(Float64[0 0 0 0;0 1 1 0;0 0 1 1];
                      tris=Int32[1 1;2 3;3 4],tri_tag=Int32[6,6])
        vertical_result=recombine_triangles(vertical;preserve_segments=false)
        @test vertical_result.blocks[1].msh==3
        @test vertical_result.blocks[1].nodes==reshape(Int32[1,2,3,4],4,1)
    end

    @testset "empty, curve-only, and error paths" begin
        empty_mesh=Mesh(zeros(3,0))
        empty_result=recombine_triangles(empty_mesh)
        @test size(empty_result.coords)==(3,0)
        @test isempty(empty_result.blocks)
        @test Tessella.Elements.validate(empty_result).ok

        curve=Mesh(Float64[0 1;0 0;0 0];segs=reshape(Int32[1,2],2,1),
                   seg_tag=Int32[5])
        curve_result=recombine_triangles(curve)
        @test curve_result.blocks[1].msh==1
        @test curve_result.blocks[1].nodes==curve.segs
        @test curve_result.blocks[1].tags==curve.seg_tag

        tetra=Mesh(Float64[0 1 0 0;0 0 1 0;0 0 0 1];
                   tets=reshape(Int32[1,2,3,4],4,1))
        @test_throws ArgumentError recombine_triangles(tetra)
        @test_throws ArgumentError recombine_triangles(_recombine_square();min_quality=-eps())
        @test_throws ArgumentError recombine_triangles(_recombine_square();min_quality=1.1)
        @test_throws ArgumentError recombine_triangles(_recombine_square();min_quality=NaN)
        @test_throws ArgumentError recombine_triangles(_recombine_square();min_quality=true)
        @test_throws ArgumentError recombine_triangles(
            _recombine_square();preserve_segments=1)
        @test_throws ArgumentError recombine_triangles(
            _recombine_square();full_quad=1)
        @test_throws ArgumentError recombine_triangles(
            _recombine_square();algorithm="greedy")
        @test_throws ArgumentError recombine_triangles(
            _recombine_square();physical_names=Dict((2,0)=>"bad"))

        duplicate=Mesh(_recombine_square().coords;
                       tris=Int32[1 1;2 2;3 3],tri_tag=Int32[1,1])
        @test !validate(duplicate).ok
        @test_throws ArgumentError recombine_triangles(duplicate)
    end

    @testset "structured grid and allocation growth" begin
        grid=_recombine_grid(12)
        result=recombine_triangles(grid;preserve_segments=false)
        @test Tessella.Elements.validate(result).ok
        @test length(result.blocks)==1
        @test result.blocks[1].msh==3
        @test size(result.blocks[1].nodes)==(4,144)
        @test all(==(Int32(1)),result.blocks[1].tags)
        @test Tessella.Elements.mixed_crc(result).sha==
              "dbb1bf17965d4e011e7f51a452c6a03e4018a628ffc7e7d9b33d9fc6b922439f"

        small=_recombine_allocated(20)
        large=_recombine_allocated(40)
        @test small>0
        @test large>small
        @test large<=5.25small+262_144
        # Warmed complete public calls, including validation/output ownership.
        # The measured former three-vector corner workspace exceeds this budget
        # at all sizes; leave room for allocator size classes and fixed setup.
        for n in (16,32,64),angle in (.5,1.01)
            mesh=_recombine_grid(n)
            bytes=_recombine_angle_grid_allocated(mesh,angle)
            @test bytes<=1800n^2+16_384
        end
    end

    @testset "scale and translation invariance" begin
        for scale in (1e-300,1e-150,1.0,1e150)
            coords=Float64[0 scale scale 0;0 0 scale scale;0 0 0 0]
            mesh=Mesh(coords;tris=Int32[1 1;2 3;3 4],tri_tag=Int32[7,7])
            result=recombine_triangles(mesh;min_quality=0.99,
                                       preserve_segments=false)
            @test length(result.blocks)==1
            @test result.blocks[1].msh==3
            @test result.blocks[1].nodes==reshape(Int32[1,2,3,4],4,1)
        end
        for offset in (1e-100,1e100,1e150)
            width=16eps(offset)
            coords=Float64[offset offset+width offset+width offset;
                           offset offset offset+width offset+width;
                           0 0 0 0]
            mesh=Mesh(coords;tris=Int32[1 1;2 3;3 4],tri_tag=Int32[7,7])
            result=recombine_triangles(mesh;min_quality=0.99,
                                       preserve_segments=false)
            @test result.blocks[1].msh==3
        end
    end

    @testset "blossom matching and full-quad" begin
        matching_size(mate)=count(i->mate[i]>i, eachindex(mate))
        consistent(mate)=all(i->mate[i]==0 || mate[mate[i]]==i, eachindex(mate))

        path=Tessella.Recombine._edmonds_matching(4,[Int[2],Int[1,3],Int[2,4],Int[3]])
        @test matching_size(path)==2
        @test consistent(path)

        triangle=Tessella.Recombine._edmonds_matching(3,[Int[2,3],Int[1,3],Int[1,2]])
        @test matching_size(triangle)==1
        @test consistent(triangle)

        flower=Tessella.Recombine._edmonds_matching(4,[Int[2,3,4],Int[1,3],Int[1,2],Int[1]])
        @test matching_size(flower)==2
        @test consistent(flower)
        @test flower[4]!=0

        isolated=Tessella.Recombine._edmonds_matching(1,[Int[]])
        @test isolated==[0]
        @test Tessella.Recombine._edmonds_matching(0,Vector{Int}[])==Int[]

        # Regression for an alternating-tree/blossom collision that previously
        # indexed parent vertex 0. Exhaustive simple graphs through six vertices
        # are compared with an independent subset-search oracle below.
        regression=[Int[3,5],Int[3,5],Int[1,2,4,5],Int[3],Int[1,2,3]]
        regression_mate=Tessella.Recombine._edmonds_matching(5,regression)
        @test matching_size(regression_mate)==2
        @test consistent(regression_mate)
        @test _exhaustive_matching_mismatch()===nothing

        grid=_recombine_grid(12)
        blossom=recombine_triangles(grid; algorithm=:blossom, full_quad=true,
                                    preserve_segments=false)
        greedy=recombine_triangles(grid; preserve_segments=false)
        @test Tessella.Elements.validate(blossom).ok
        @test length(blossom.blocks)==1
        @test blossom.blocks[1].msh==3
        @test size(blossom.blocks[1].nodes)==(4,144)
        @test all(==(Int32(1)), blossom.blocks[1].tags)
        @test Tessella.Elements.mixed_crc(blossom)==
              Tessella.Elements.mixed_crc(recombine_triangles(grid;
                  algorithm=:blossom, full_quad=true, preserve_segments=false))
        @test size(greedy.blocks[1].nodes)==(4,144)

        @test_throws ArgumentError recombine_triangles(grid; full_quad=true)
        @test_throws ArgumentError recombine_triangles(grid; algorithm=:nope)
        odd=Mesh(Float64[0 1 0;0 0 1;0 0 0];
                 tris=reshape(Int32[1,2,3],3,1), tri_tag=Int32[1])
        @test_throws ArgumentError recombine_triangles(odd; algorithm=:blossom,
                                                       full_quad=true)
        split=recombine_triangles(_recombine_square(tags=Int32[7,8]);
                                  algorithm=:blossom, preserve_segments=false)
        @test [block.msh for block in split.blocks]==[2]
        @test size(split.blocks[1].nodes,2)==2
    end

    @testset "public documentation" begin
        @test isempty(Base.Docs.undocumented_names(Tessella.Recombine;private=false))
        @test isempty(Test.detect_ambiguities(Tessella.Recombine;recursive=true))
    end
end

module RecombineNativePriorityRegressionTests
using Test,Tessella
const R=Tessella.Recombine
const api=Tessella.API
include(joinpath(@__DIR__,"..","artifacts","recombine_priority_primary.txt"))

canonical(nodes)=minimum(ntuple(i->nodes[mod1(i+shift,length(nodes))],length(nodes))
    for shift in 0:length(nodes)-1)

function records(elements)
    result=Dict{Int,Vector{Tuple}}()
    for (kind,labels,nodes) in zip(elements...)
        width=kind==2 ? 3 : kind==3 ? 4 : error("Unexpected surface element $kind")
        result[Int(kind)]=sort!([(Int(labels[i]),canonical(nodes[width*(i-1)+1:width*i]))
            for i in eachindex(labels)])
    end
    result
end

function public_recombine(coords,triangles,labels,algorithm)
    input=(copy(coords),copy(triangles),copy(labels))
    try
        api.initialize();api.option("Mesh.Renumber",0)
        api.option("Mesh.RecombinationAlgorithm",algorithm)
        api.model.add_discrete_entity(2,1)
        api.mesh.add_nodes(2,1,labels,vec(coords))
        api.mesh.add_elements_by_type(1,2,collect(30001:30000+size(triangles,2)),
            labels[vec(triangles)])
        # Deleting a point element preserves the native factory high-water mark.
        api.model.add_discrete_entity(0,99)
        api.mesh.add_nodes(0,99,[70001],[10.,10.,10.])
        api.mesh.add_elements_by_type(99,15,[90001],[70001])
        api.mesh.remove_elements(0,99,[90001])
        before=api.mesh.get_nodes()
        api.mesh.set_recombine(2,1,45.);api.mesh.recombine()
        @test api.mesh.get_nodes()==before
        @test (coords,triangles,labels)==input
        return records(api.mesh.get_elements(2,1)),api.mesh.get_max_element_tag()
    finally
        api.finalize()
    end
end

function patches(count,kind)
    width=kind===:healthy ? 1. : kind===:nan ? 1e160 : 1e200
    height=kind===:healthy ? 1. : kind===:nan ? 1e-160 : 1e-200
    nodes_per=kind===:healthy ? 4 : 5
    cells_per=kind===:healthy ? 2 : 3
    points=Matrix{Float64}(undef,3,nodes_per*count)
    cells=Matrix{Int32}(undef,3,cells_per*count)
    for patch in 1:count
        offset=nodes_per*(patch-1)
        if kind===:healthy
            points[:,offset+1:offset+4]=Float64[0 1 1 0;0 0 1 1;0 0 0 0].+
                [3.0*(patch-1),0.,0.]
            cells[:,2patch-1:2patch]=Int32[1 1;2 3;3 4].+Int32(offset)
        else
            points[:,offset+1:offset+5]=Float64[0 width 2width 0 2width;
                0 width 2width 0 2width;0 0 2height height 3height].+
                [6width*(patch-1),6width*(patch-1),5height*(patch-1)]
            cells[:,3patch-2:3patch]=Int32[1 1 3;2 3 5;3 4 4].+Int32(offset)
        end
    end
    points,cells
end

function projected_area(coords,labels,entries,axes)
    positions=Dict(Int(labels[i])=>i for i in eachindex(labels))
    points=Rational{BigInt}.(coords[collect(axes),:])
    area=zero(Rational{BigInt})
    for block in values(entries),(_,cell) in block
        nodes=[positions[Int(label)] for label in cell]
        area+=sum(points[1,nodes[i]]*points[2,nodes[mod1(i+1,length(nodes))]]-
            points[1,nodes[mod1(i+1,length(nodes))]]*points[2,nodes[i]]
            for i in eachindex(nodes))/2
    end
    area
end

@noinline function priority_sum(coords,triangles,candidate,repetitions)
    result=0.
    for _ in 1:repetitions
        result+=R._gmsh_recombine_pair_priority(coords,triangles,candidate)
    end
    result
end

@noinline function sort_workspace!(workspace,input;heap=false)
    copyto!(workspace,input)
    if heap
        R._gmsh_pair_introsort!(workspace,1,length(workspace)+1,0)
    else
        R._gmsh_pair_sort!(workspace)
    end
    workspace
end

@testset "Pinned native directed-edge priority and public labels" begin
    base=Int32[1 1 3;2 3 5;3 4 4]
    orders=((1,2,3),(1,3,2),(2,1,3),(2,3,1),(3,1,2),(3,2,1))
    tagsets=(Int32[10003+7i for i in 0:4],Int32[10003+7i for i in 4:-1:0],
        Int32[10024,10003,10031,10010,10017])
    for (wind,reverse_winding) in enumerate((false,true)),(tag_index,labels) in enumerate(tagsets),
            (order_index,order) in enumerate(orders),algorithm in (0,1)
        coords=Float64[0 1e160 2e160 0 2e160;0 1e160 2e160 0 2e160;
            0 0 2e-160 1e-160 3e-160]
        triangles=base[:,collect(order)]
        reverse_winding && (triangles=triangles[[1,3,2],:])
        leftover=PRIMARY_PRIORITY_SINGLE_LEFTOVER[wind][tag_index,order_index]
        quad=leftover==3 ? (1,2,3,4) : (1,3,5,4)
        reverse_winding && (quad=(quad[1],quad[4],quad[3],quad[2]))
        remaining=findfirst(==(leftover),order)
        expected=Dict(2=>Tuple[(30000+remaining,canonical(labels[triangles[:,remaining]]))],
            3=>Tuple[(90004,canonical(labels[collect(quad)]))])
        actual,maximum=public_recombine(coords,triangles,labels,algorithm)
        @test actual==expected
        @test maximum==90004
        input=Dict(2=>Tuple[(30000+i,Tuple(labels[triangles[:,i]])) for i in axes(triangles,2)])
        @test projected_area(coords,labels,actual,(1,3))==projected_area(coords,labels,input,(1,3))
    end

    # Eighteen through sixty-six constructors cross native introsort's
    # sixteen-entry partition boundary. These are captured primary identities,
    # including creation labels, rather than expectations from the sorter.
    for kind in (:healthy,:nan,:zero),count in (9,17,33),algorithm in (0,1)
        coords,triangles=patches(count,kind)
        labels=Int32[10003+7i for i in 0:size(coords,2)-1]
        nodes_per=kind===:healthy ? 4 : 5
        constructors=kind===:healthy ? count : 2count
        creation=PRIMARY_PRIORITY_LARGE_CREATION[(kind,count,algorithm)]
        expected=Dict(3=>Tuple[(90001+constructors+i,
            canonical(labels[collect(local_nodes).+nodes_per*(patch-1)]))
            for (i,(patch,local_nodes)) in enumerate(creation)])
        if kind!==:healthy
            expected[2]=sort!(Tuple[(30000+3patch-(local_nodes==(1,2,3,4) ? 0 : 2),
                canonical(labels[triangles[:,3patch-(local_nodes==(1,2,3,4) ? 0 : 2)]]))
                for (patch,local_nodes) in creation])
        end
        actual,maximum=public_recombine(coords,triangles,labels,algorithm)
        @test actual==expected
        @test maximum==90001+constructors+count
        input=Dict(2=>Tuple[(30000+i,Tuple(labels[triangles[:,i]])) for i in axes(triangles,2)])
        projection_axes=kind===:healthy ? (1,2) : (1,3)
        @test projected_area(coords,labels,actual,projection_axes)==
            projected_area(coords,labels,input,projection_axes)
    end
    # Native publication of successful Blossom pairs is visited in ascending
    # high triangle endpoint order, including permuted original triangle tags.
    for count in (17,33),order in (:reverse,:split),tag_order in (:forward,:reverse,:shift),
            algorithm in (0,1)
        coords,triangles=patches(count,:healthy)
        triangles=order===:reverse ? reverse(triangles;dims=2) :
            triangles[:,vcat(collect(2:2:2count),collect(2count-1:-2:1))]
        labels=Int32[10003+7i for i in 0:size(coords,2)-1]
        tag_order===:reverse && reverse!(labels)
        tag_order===:shift && (labels=circshift(labels,5))
        creation=PRIMARY_PRIORITY_PERMUTED_CREATION[(count,order,tag_order,algorithm)]
        expected=Dict(3=>Tuple[(90001+count+i,
            canonical(labels[collect(local_nodes).+4*(patch-1)]))
            for (i,(patch,local_nodes)) in enumerate(creation)])
        actual,maximum=public_recombine(coords,triangles,labels,algorithm)
        @test actual==expected
        @test maximum==90001+2count
        input=Dict(2=>Tuple[(30000+i,Tuple(labels[triangles[:,i]])) for i in axes(triangles,2)])
        @test projected_area(coords,labels,actual,(1,2))==projected_area(coords,labels,input,(1,2))
    end
end

@testset "Priority and preallocated native sort resource bounds" begin
    for width in (1.,1e160,1e200)
        coords=Float64[0 width 2width 0;0 width 2width 0;0 0 2/width 1/width]
        triangles=Int32[1 1;2 3;3 4]
        candidate=R._QuadCandidate(Int32(1),Int32(2),(Int32(1),Int32(2),Int32(3),Int32(4)),
            .5,(Int32(1),Int32(3)))
        priority_sum(coords,triangles,candidate,1000)
        @test minimum(@allocated(priority_sum(coords,triangles,candidate,1000)) for _ in 1:3)==0
    end
    for count in (0,1,16,17,33,129,1024)
        input=[(priority=Float64(mod(37i,count+1)),identity=i) for i in 1:count]
        # The selected sizes are coprime to 37, hence all finite keys are unique.
        expected=sort(input;by=pair->pair.priority)
        workspace=copy(input)
        @test sort_workspace!(workspace,input)==expected
        @test minimum(@allocated(sort_workspace!(workspace,input)) for _ in 1:3)==0
        if count>16
            @test sort_workspace!(workspace,input;heap=true)==expected
            @test minimum(@allocated(sort_workspace!(workspace,input;heap=true)) for _ in 1:3)==0
        end
        @test workspace!==input
    end
end
end # module RecombineNativePriorityRegressionTests
