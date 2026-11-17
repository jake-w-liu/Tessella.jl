module OrderedMatchingNativeStateTests
using Test, Tessella
const O=Tessella.OrderedMatching

# Literal per-vertex states are (label, edge, dual, flag), with 0-based edges.
# Label 2 marks an integral match, 1 a directed odd half-cycle, and 0 a free
# vertex. The rounded state supplies the seed for the Blossom continuation.
const CASES=[
    # No retained cycles after the fractional search.
    (name="weighted_39",count=6,ne=6,nm=3,
        edges=[(1,2,3),(1,3,5),(1,4,7),(1,6,0),(2,5,6),(3,4,10)],
        selected=[6,5,4],cost=16,cycles=0,
        init=[(2,3,0,0),(2,4,6,0),(0,-1,10,0),(0,-1,10,0),(2,4,6,0),(2,3,0,0)],
        fractional=[(2,3,0,0),(2,4,6,0),(2,5,10,0),(2,5,10,0),(2,4,6,0),(2,3,0,0)],
        rounded=[(2,3,0,0),(2,4,6,0),(2,5,10,0),(2,5,10,0),(2,4,6,0),(2,3,0,0)]),
    # Fractional duals legitimately go negative on this case.
    (name="weighted_43",count=6,ne=7,nm=3,
        edges=[(1,2,5),(1,3,7),(1,5,0),(1,6,2),(2,4,6),(2,5,8),(3,4,1)],
        selected=[7,6,4],cost=11,cycles=0,
        init=[(2,2,0,0),(0,-1,10,0),(2,6,1,0),(2,6,1,0),(2,2,0,0),(0,-1,4,0)],
        fractional=[(2,3,-3,0),(2,5,13,0),(2,6,3,0),(2,6,-1,0),(2,5,3,0),(2,3,7,0)],
        rounded=[(2,3,-3,0),(2,5,13,0),(2,6,3,0),(2,6,-1,0),(2,5,3,0),(2,3,7,0)]),
    # Odd cycles close and reopen inside the fractional search itself; the
    # published fractional state is again fully integral.
    (name="weighted_59",count=6,ne=8,nm=3,
        edges=[(1,2,2),(1,3,4),(1,5,8),(2,3,1),(2,4,3),(2,5,5),(2,6,7),(3,4,9)],
        selected=[8,3,7],cost=24,cycles=0,
        init=[(0,-1,3,0),(2,3,1,0),(2,3,1,0),(0,-1,5,0),(0,-1,9,0),(0,-1,13,0)],
        fractional=[(2,2,1,0),(2,6,-5,0),(2,7,7,0),(2,7,11,0),(2,2,15,0),(2,6,19,0)],
        rounded=[(2,2,1,0),(2,6,-5,0),(2,7,7,0),(2,7,11,0),(2,2,15,0),(2,6,19,0)]),
    # Two retained half-rings cover every vertex. The rounded record keeps two
    # ring edges and releases one vertex per ring; the Blossom phase then
    # joins the two free vertices through edge (2,4).
    (name="weighted_353",count=6,ne=8,nm=3,
        edges=[(1,2,6),(1,6,3),(2,3,5),(2,4,7),(2,6,0),(3,4,2),(3,5,4),(4,5,8)],
        selected=[4,7,2],cost=14,cycles=2,
        init=[(0,-1,6,0),(2,4,0,0),(2,5,2,0),(2,5,2,0),(0,-1,6,0),(2,4,0,0)],
        fractional=[(1,1,9,0),(1,0,3,0),(1,6,-2,0),(1,5,6,0),(1,7,10,0),(1,4,-3,0)],
        rounded=[(2,1,9,0),(0,-1,3,0),(2,6,-2,0),(0,-1,6,0),(2,6,10,0),(2,1,-3,0)]),
    (name="weighted_355",count=6,ne=8,nm=3,
        edges=[(1,2,7),(1,3,9),(1,6,4),(2,5,10),(2,6,1),(3,4,3),(3,5,5),(4,5,9)],
        selected=[6,4,3],cost=17,cycles=2,
        init=[(0,-1,7,0),(2,4,1,0),(2,5,3,0),(2,5,3,0),(0,-1,7,0),(2,4,1,0)],
        fractional=[(1,2,10,0),(1,0,4,0),(1,6,-1,0),(1,5,7,0),(1,7,11,0),(1,4,-2,0)],
        rounded=[(2,2,10,0),(0,-1,4,0),(2,6,-1,0),(0,-1,7,0),(2,6,11,0),(2,2,-2,0)]),
    (name="weighted_557",count=6,ne=9,nm=3,
        edges=[(1,2,9),(1,4,2),(1,5,4),(1,6,6),(2,3,8),(2,4,10),(2,5,1),
            (3,5,7),(4,6,2)],
        selected=[5,3,9],cost=14,cycles=2,
        init=[(2,1,2,0),(2,6,1,0),(0,-1,13,0),(2,1,2,0),(2,6,1,0),(0,-1,2,0)],
        fractional=[(1,1,6,0),(1,6,2,0),(1,4,14,0),(1,8,-2,0),(1,7,0,0),(1,3,6,0)],
        rounded=[(2,1,6,0),(2,6,2,0),(0,-1,14,0),(2,1,-2,0),(2,6,0,0),(0,-1,6,0)]),
    (name="weighted_837",count=6,ne=9,nm=3,
        edges=[(1,2,6),(1,4,10),(1,6,3),(2,5,9),(2,6,0),(3,4,2),(3,5,4),
            (4,5,8),(4,6,10)],
        selected=[6,4,3],cost=14,cycles=2,
        init=[(0,-1,6,0),(2,4,0,0),(2,5,2,0),(2,5,2,0),(0,-1,6,0),(2,4,0,0)],
        fractional=[(1,2,9,0),(1,0,3,0),(1,6,-2,0),(1,5,6,0),(1,7,10,0),(1,4,-3,0)],
        rounded=[(2,2,9,0),(0,-1,3,0),(2,6,-2,0),(0,-1,6,0),(2,6,10,0),(2,2,-3,0)]),
    # A wider no-cycle case exercising checked dual arithmetic.
    (name="5/0/0",count=16,ne=28,nm=8,
        edges=[(9,12,367),(10,11,367),(5,7,367),(6,8,367),(6,15,999),
            (1,5,1000),(3,9,1000),(7,16,999),(2,8,1000),(10,13,999),
            (4,11,1000),(12,14,999),(7,9,1000),(8,10,999),(3,16,999),
            (2,13,10000),(11,12,1000),(4,14,999),(5,6,1000),(1,15,10000),
            (13,2,100000),(4,14,100000),(3,16,100000),(1,15,100000),
            (13,4,100000),(14,3,100000),(16,1,100000),(2,15,100000)],
        selected=[6,9,13,17,10,18,5,15],cost=7996,cycles=0,
        init=[(2,5,1633,0),(0,-1,1633,0),(2,14,999,0),(2,17,999,0),
            (2,5,367,0),(2,3,367,0),(0,-1,367,0),(2,3,367,0),
            (2,0,367,0),(2,1,367,0),(2,1,367,0),(2,0,367,0),
            (0,-1,1631,0),(2,17,999,0),(0,-1,1631,0),(2,14,999,0)],
        fractional=[(2,5,4165,0),(2,8,1633,0),(2,14,2899,0),(2,17,1633,0),
            (2,5,-2165,0),(2,4,367,0),(2,12,2899,0),(2,8,367,0),
            (2,12,-899,0),(2,9,367,0),(2,16,367,0),(2,16,1633,0),
            (2,9,1631,0),(2,17,365,0),(2,4,1631,0),(2,14,-901,0)],
        rounded=[(2,5,4165,0),(2,8,1633,0),(2,14,2899,0),(2,17,1633,0),
            (2,5,-2165,0),(2,4,367,0),(2,12,2899,0),(2,8,367,0),
            (2,12,-899,0),(2,9,367,0),(2,16,367,0),(2,16,1633,0),
            (2,9,1631,0),(2,17,365,0),(2,4,1631,0),(2,14,-901,0)])]

# Active retained rings only: `cycles` keeps reopened rings in storage.
active_cycles(state)=[cycle for (id,cycle) in enumerate(state.cycles)
    if state.cycle_owner[first(cycle[1])]==id]

@testset "Literal greedy, fractional and rounded matching states" begin
    for fixture in CASES
        count=fixture.count
        before=copy(fixture.edges);graph=O._graph(count,fixture.edges)

        # Production greedy init reproduces `S init` exactly: matched flags,
        # 0-based matched edge identity, and the raised vertex duals.
        init_dual=Vector{Int64}(undef,count);init_mate=zeros(Int,count)
        O._greedy_singleton_init!(graph,init_dual,init_mate)
        for vertex in 1:count
            label,edge,dual,flag=fixture.init[vertex]
            @test flag==0
            @test label==(init_mate[vertex]!=0 ? 2 : 0)
            @test edge==init_mate[vertex]-1
            @test dual==init_dual[vertex]
        end

        # The fractional state reproduces `S fractional`: identical duals,
        # identical matched edge indices, and identical directed ring edges
        # for every vertex on a retained half-cycle.
        state=O._fractional_singleton_state(graph)
        @test state!==nothing
        @test length(active_cycles(state))==fixture.cycles
        for vertex in 1:count
            label,edge,dual,flag=fixture.fractional[vertex]
            @test flag==0
            @test dual==state.dual[vertex]
            if label==1
                id=state.cycle_owner[vertex]
                @test id!=0
                vertices,edges=state.cycles[id]
                at=only(findall(==(vertex),vertices))
                @test edge==edges[at]-1
            else
                @test edge==state.mate[vertex]-1
                @test state.cycle_owner[vertex]==0
            end
        end

        # The Blossom seed reproduces `S rounded`: unchanged duals, the same
        # rounded ring mates, and the same released base vertex per ring.
        seed=O._fractional_singleton_seed(graph,state)
        for vertex in 1:count
            label,edge,dual,flag=fixture.rounded[vertex]
            @test flag==0
            @test dual==state.dual[vertex]
            if label==0 || edge<0
                @test seed.mate[vertex]==0
            else
                @test seed.mate[vertex]!=0
                a,b,_=fixture.edges[edge+1]
                @test minmax(a,b)==minmax(vertex,seed.mate[vertex])
            end
        end

        # Directed output identity and certified optimal cost.
        result=O.ordered_perfect_matching(count,fixture.edges)
        @test result!==nothing && result.selected==fixture.selected
        @test fixture.edges==before
        @test sort([vertex for i in result.selected for vertex in fixture.edges[i][1:2]])==collect(1:count)
        @test sum(fixture.edges[i][3] for i in result.selected)==fixture.cost
        # Weak-duality bound: the certified singleton duals cover the integral
        # optimum from below (retained rings can leave it strictly below).
        @test sum(Int128,state.dual)<=2Int128(fixture.cost)
        if O._bipartite(graph)
            # Bipartite inputs take the ordered alternating-tree path; its
            # greedy init shares `_greedy_singleton_init!` and the resulting
            # cover must still reach the same optimum.
            bip=O._bipartite_perfect(graph)
            @test bip!==nothing
            if bip!==nothing
                @test sort!([minmax(v,O._other(graph,bip[v],v))
                        for v in 1:count if bip[v]!=0 &&
                        O._other(graph,bip[v],v)>v])==
                    sort!([minmax(fixture.edges[i][1],fixture.edges[i][2])
                        for i in fixture.selected])
            end
        end
    end
end

end # module OrderedMatchingNativeStateTests
