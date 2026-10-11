module NativeMatchingAdapterTestsV1
using Tessella, Test
const M = Tessella.Model

# Pinned SDK 4.15.2: API PC53, angle 45, Blossom, topology/repositioning off,
# minimum quality zero, sparse public labels, and actual embedded Line cells.
# The input graph and output labels are independent oracle observations. The
# null entry maps to triangle one; selected closure edges emit no Quad.
const PRIMARY_DRIVER_SHA256 = "9e572b13e3fb082909ab063b8d9bf7ab746aaad4017391b102edf522ae36ed5b"
const PRIMARY_RECEIPT_SHA256 = "57cd6de4e933c823dcf3b66a83ca90dbce882287c77e5831fbec7e3b947f0983"
const TRIANGLES = NTuple{3,Int32}[(1,2,5),(1,5,4),(2,3,6),(2,6,5)]
const NODE_TAGS = UInt64[10003,10010,10017,10024,10031,10038]
function state()
    world = Float64[0 1 2 0 1 2; 0 0 0 1 1 1; 0 0 0 0 0 0]
    return M._NativeSurfaceRecombine(world,copy(world[1:2,:]),
        fill((2,Int32(1)),6),copy(NODE_TAGS),copy(TRIANGLES),
        UInt64[30001,30002,30003,30004],NTuple{4,Int32}[],UInt64[],
        trues(6),UInt64(90001),nothing,"primary periodic closure control")
end
function oriented_cycle(nodes)
    first = argmin(nodes)
    return ntuple(i -> nodes[mod1(first+i-1,length(nodes))], length(nodes))
end
function exact_area(world,cells)
    total = Rational{BigInt}(0)
    for nodes in cells
        for i in eachindex(nodes)
            a = nodes[i]; b = nodes[mod1(i+1,length(nodes))]
            total += Rational{BigInt}(world[1,a])*Rational{BigInt}(world[2,b]) -
                Rational{BigInt}(world[2,a])*Rational{BigInt}(world[1,b])
        end
    end
    return total/2
end
@testset "Native matching closure: strict primary birth and geometry controls" begin
    for (name,edges,pair_count,periodic,quad_tags,tri_tags) in (
        ("plain",Tuple{Int32,Int32}[],3,[(1,2),(1,3),(2,4),(3,4)],
            UInt64[90005,90006],UInt64[]),
        ("outer only",Tuple{Int32,Int32}[(1,2)],3,[(2,0),(3,0),(2,4),(3,4)],
            UInt64[90005,90006],UInt64[]),
        ("diagonal only",Tuple{Int32,Int32}[(1,5)],2,[(1,2),(1,3),(2,4),(3,4)],
            UInt64[90004],UInt64[30001,30002]),
        ("outer and diagonal",Tuple{Int32,Int32}[(1,2),(1,5)],2,
            [(2,0),(3,0),(2,4),(3,4)],UInt64[90004],UInt64[30001,30002]))
        @testset "$name" begin
            protected = Set(edges)
            graph_state = state()
            pairs,closure = M._native_recombine_pairs!(graph_state,
                graph_state.node_tags,protected)
            @test length(pairs) == pair_count
            @test closure == periodic
            @test graph_state.high_water == UInt64(90001+pair_count)
            if pair_count == 2
                @test [(p.first,p.second,p.cost) for p in pairs] ==
                    [(3,4,367),(1,4,606)]
            end
            if name == "outer and diagonal"
                result = M._native_recombine_matching_result(4,pairs,closure,name)
                @test result !== nothing
                @test result.edges == [(3,4,367),(1,4,606),(2,1,100000),
                    (3,1,100000),(2,4,100000),(3,4,100000)]
                @test result.selected == [3,1]
                @test sum(result.edges[i][3] for i in result.selected) == 100367
                @test result.physical_count == 2
            end
            actual = state()
            world = copy(actual.world); uv = copy(actual.uv)
            owners = copy(actual.owners); labels = copy(actual.node_tags)
            before_area = exact_area(world,actual.triangles)
            M._native_recombine_initial!(actual,actual.node_tags,protected,true,45.0)
            expected_quads = length(quad_tags) == 2 ?
                [(Int32(1),Int32(2),Int32(5),Int32(4)),
                 (Int32(2),Int32(3),Int32(6),Int32(5))] :
                [(Int32(2),Int32(3),Int32(6),Int32(5))]
            @test oriented_cycle.(actual.quads) == expected_quads
            @test actual.quad_tags == quad_tags
            @test actual.triangle_tags == tri_tags
            @test actual.triangles == (isempty(tri_tags) ? NTuple{3,Int32}[] : TRIANGLES[1:2])
            @test actual.high_water == last(quad_tags)
            @test reinterpret(UInt64,vec(actual.world)) == reinterpret(UInt64,vec(world))
            @test reinterpret(UInt64,vec(actual.uv)) == reinterpret(UInt64,vec(uv))
            @test actual.owners == owners
            @test actual.node_tags == labels == NODE_TAGS
            @test all(actual.alive)
            @test exact_area(actual.world,actual.triangles) +
                exact_area(actual.world,actual.quads) == before_area == 2
            # Mutating the detached fixture must leave the retained inputs alone.
            actual.world[1,1] = 99
            actual.node_tags[1] = 99
            @test world[1,1] == 0 && NODE_TAGS[1] == 10003
            @test TRIANGLES == NTuple{3,Int32}[(1,2,5),(1,5,4),(2,3,6),(2,6,5)]
        end
    end
end

@testset "Directed closure identity, safe null self-loop, and validation" begin
    for physical_cost in (100000,200001)
        pairs = [(first=1,second=2,cost=physical_cost)]
        closure = [(1,2)]
        before = deepcopy((pairs,closure))
        result = M._native_recombine_matching_result(2,pairs,closure,"parallel closure")
        @test result.edges == [(1,2,physical_cost),(1,2,100000)]
        @test result.selected == [2]
        @test result.physical_count == 1
        @test all(index > result.physical_count for index in result.selected)
        @test (pairs,closure) == before
    end
    # The primary shifted control hung and is not a passing oracle. This is an
    # independent safe-loop contract: a loop covers no vertex, the other edge
    # retains its original direction and closes the two-vertex graph.
    result = M._native_recombine_matching_result(2,NamedTuple[],[(1,0),(2,0)],"safe null loop")
    @test result.edges == [(1,1,100000),(2,1,100000)]
    @test result.selected == [2]
    @test result.physical_count == 0
    @test M._native_recombine_matching(2,NamedTuple[],[(1,0),(2,0)],"safe null loop") == [(2,1)]
    @test M._native_recombine_matching(0,NamedTuple[],Tuple{Int,Int}[],"empty") == Tuple{Int,Int}[]
    @test M._native_recombine_matching(3,NamedTuple[],Tuple{Int,Int}[],"odd") === nothing
    @test M._native_recombine_matching(4,[(first=1,second=2,cost=0)],
        Tuple{Int,Int}[],"uncovered") === nothing
    @test_throws ArgumentError M._native_recombine_matching(2,
        [(first=1,second=3,cost=0)],Tuple{Int,Int}[],"invalid endpoint")
    @test_throws ArgumentError M._native_recombine_matching(2,
        [(first=1,second=2,cost=typemax(Int))],Tuple{Int,Int}[],"invalid cost")
end

function matching_measure(n)
    pairs = [(first=i,second=i+1,cost=-i) for i in 1:2:n]
    expected = [(i,i+1) for i in 1:2:n]
    for _ in 1:2
        @test M._native_recombine_matching(n,pairs,Tuple{Int,Int}[],"resource control") == expected
    end
    samples = Int[]
    for _ in 1:2
        push!(samples,@allocated M._native_recombine_matching(n,pairs,
            Tuple{Int,Int}[],"resource control"))
    end
    return maximum(samples)
end
@testset "Matching adapter allocation scales with independent components" begin
    matching_measure(32)
    previous = 0
    for n in (1024,2048,4096)
        bytes = matching_measure(n)
        @test previous == 0 || bytes <= 2.15previous + 65536
        println("MATCHING_ADAPTER_RESOURCE n=",n," bytes=",bytes)
        previous = bytes
    end
end
end
