module OrderedMatchingSeededTests
using Test, Random
const Tessella = getfield(parentmodule(@__MODULE__), :Tessella)
const O = Tessella.OrderedMatching
const W = Tessella.WeightedMatching

# Enumerate independent vertex pairings, allowing unmatched vertices. The
# objective is cardinality first, then exact integer cost; no production
# traversal, duals, or blossom state enter this oracle.
function exact_matching(count, edges)
    costs = Dict{Tuple{Int,Int},Int128}()
    for (a, b, cost) in edges
        a == b && continue
        pair = minmax(a, b)
        costs[pair] = min(get(costs, pair, typemax(Int128)), Int128(cost))
    end
    function search(remaining)
        isempty(remaining) && return (0, Int128(0))
        first = remaining[1]
        tail = remaining[2:end]
        best = search(tail)
        for (index, second) in enumerate(tail)
            cost = get(costs, minmax(first, second), nothing)
            cost === nothing && continue
            rest = [tail[i] for i in eachindex(tail) if i != index]
            pairs, subtotal = search(rest)
            trial = (pairs + 1, subtotal + cost)
            if trial[1] > best[1] || (trial[1] == best[1] && trial[2] < best[2])
                best = trial
            end
        end
        return best
    end
    return search(collect(1:count))
end

function matching_objective(mate, edges)
    costs = Dict{Tuple{Int,Int},Int128}()
    for (a, b, cost) in edges
        a == b && continue
        pair = minmax(a, b)
        costs[pair] = min(get(costs, pair, typemax(Int128)), Int128(cost))
    end
    pairs = 0
    total = Int128(0)
    for vertex in eachindex(mate)
        other = mate[vertex]
        other == 0 && continue
        1 <= other <= length(mate) && other != vertex && mate[other] == vertex ||
            error("matching does not pair distinct reciprocal vertices")
        vertex < other || continue
        pairs += 1
        total += costs[minmax(vertex, other)]
    end
    return pairs, total
end

@testset "Seeded odd cycles preserve finite optimum and native directed choices" begin
    # Literal native cases that keep fractional rings and therefore exercise
    # the continuation, rather than the integral singleton fast path.
    fixtures = (
        (edges=[(1,2,6),(1,6,3),(2,3,5),(2,4,7),(2,6,0),
                (3,4,2),(3,5,4),(4,5,8)], selected=[4,7,2]),
        (edges=[(1,3,5),(1,4,7),(1,6,0),(2,4,4),(2,5,6),
                (3,4,10),(3,6,3),(4,5,5),(4,6,7)], selected=[2,5,7]),
        (edges=[(1,5,8),(1,6,10),(2,3,1),(2,4,3),(2,6,7),
                (3,4,9),(5,6,8)], selected=[6,1,5]),
        (edges=[(1,2,2),(1,3,4),(1,4,6),(1,6,10),(2,3,1),
                (2,4,3),(2,5,5),(3,5,0),(3,6,2),(5,6,8)], selected=[6,8,4]),
    )
    for fixture in fixtures
        before = copy(fixture.edges)
        state = O._fractional_singleton_state(O._graph(6, fixture.edges))
        @test state !== nothing && any(!iszero, state.cycle_owner)
        result = O.ordered_perfect_matching(6, fixture.edges)
        @test result !== nothing && result.selected == fixture.selected
        @test matching_objective(result.mate, fixture.edges) == exact_matching(6, fixture.edges)
        @test fixture.edges == before
    end

    triangles = [(1,2,-2),(2,3,-2),(3,1,-2),(4,5,-2),(5,6,-2),(6,4,-2)]
    @test O.ordered_perfect_matching(6, triangles) === nothing
    state = O._fractional_singleton_state(O._graph(6, triangles))
    @test state !== nothing && any(!iszero, state.cycle_owner)
    seed = O._fractional_singleton_seed(O._graph(6, triangles), state)
    @test matching_objective(W.minimum_weight_matching(6, triangles, seed), triangles) ==
        exact_matching(6, triangles)

    linked = vcat(triangles, [(3,4,3),(4,3,3),(1,1,-100)])
    @test matching_objective(O.ordered_perfect_matching(6, linked).mate, linked) ==
        exact_matching(6, linked)
    @test O.ordered_perfect_matching(0, Tuple{Int,Int,Int}[]) == (mate=Int[], selected=Int[])
    @test O.ordered_perfect_matching(3, [(1,2,0),(2,3,0),(3,1,0)]) === nothing
    @test_throws ArgumentError O.ordered_perfect_matching(2, [(0,1,0)])
    @test_throws ArgumentError O.ordered_perfect_matching(2, [(1,2,true)])
end

# The retained-cycle certificate publishes each ring anchored at its base and
# traversing the closing edge against its stored direction — the orientation
# the opaque native state records show. The weighted continuation consumes the
# anchored reversal of the same ring. Both conventions plus the rounded mates,
# certified duals, and deterministic selection are pinned on this literal
# weighted_353 small graph.
@testset "Retained odd-cycle rings publish observed orientation; seed consumes reversal" begin
    edges=[(1,2,6),(1,6,3),(2,3,5),(2,4,7),(2,6,0),(3,4,2),(3,5,4),(4,5,8)]
    graph=O._graph(6,edges)
    state=O._fractional_singleton_state(graph)
    @test state!==nothing
    live=[(vertices,ring_edges) for (id,(vertices,ring_edges)) in
        enumerate(state.cycles) if state.cycle_owner[first(vertices)]==id]
    @test live==[([1,6,2],[2,5,1]),([5,4,3],[8,6,7])]
    for (vertices,ring_edges) in live
        for index in eachindex(vertices)
            a,b,_=graph.edges[ring_edges[index]]
            @test Set((a,b))==Set((vertices[index],
                vertices[mod1(index+1,length(vertices))]))
        end
    end
    seed=O._fractional_singleton_seed(graph,state)
    @test seed.mate==[6,0,5,0,3,1]
    @test seed.dual==[9,3,-2,6,10,-3]
    @test seed.cycles==[[1,2,6],[5,3,4]]
    @test seed.cycles==[vcat(vs[1],reverse(vs[2:end])) for (vs,_) in live]
    result=O.ordered_perfect_matching(6,edges)
    @test result!==nothing && result.selected==[4,7,2] && result.mate==[6,4,5,2,3,1]
    @test matching_objective(result.mate,edges)==exact_matching(6,edges)
    # The same literal graph through a fresh container and a rebuilt state
    # reproduces the identical certificate and selection deterministically.
    clone=O.ordered_perfect_matching(6,Tuple(edges))
    @test clone!==nothing && clone.mate==result.mate && clone.selected==result.selected
    state2=O._fractional_singleton_state(O._graph(6,edges))
    @test [(vs,es) for (id,(vs,es)) in enumerate(state2.cycles) if
        state2.cycle_owner[first(vs)]==id]==live
end

# Named fixtures for each remaining production arm of the fractional
# singleton search, so a regression fails on the arm it broke rather than
# only through aggregate corpus counts.
@testset "Fractional singleton search covers every production arm" begin
    # Dual descent between two outer vertices: the (3,1) edge is non-tight
    # until the outer-outer half-slack update tightens it, then it closes
    # the ring. The inner vertex drops below zero through `checked_sub`.
    state=O._fractional_singleton_state(O._graph(3,[(1,2,0),(2,3,0),(3,1,1)]))
    @test state!==nothing && state.dual==[1,-1,1]
    @test state.cycles==[([3,2,1],[2,1,3])]

    # A ring whose base still holds its stem mate: closing releases the base
    # and re-augments the stem (the `prefix != 0` arm), and the stored
    # joining-edge direction publishes the flipped traversal.
    tree=[(1,2,10),(2,3,0),(3,4,0),(3,5,0),(4,6,0),(5,7,0),(6,7,0)]
    state=O._fractional_singleton_state(O._graph(7,tree))
    @test state!==nothing && state.cycles==[([3,5,7,6,4],[4,6,7,5,3])]
    @test state.mate==[1,1,0,0,0,0,0]
    seed=O._fractional_singleton_seed(O._graph(7,tree),state)
    @test seed.mate==[2,1,5,0,3,7,6] && seed.cycles==[[3,4,6,7,5]]

    # Fractional infeasibility exits through the exhausted-delta arm: the
    # alternating tree over the path component owns only inner partners, so
    # neither a tight nor a slack edge remains to grow it.
    disjoint=[(1,2,0),(2,3,0),(3,1,0),(4,5,0),(5,6,0)]
    @test O._fractional_singleton_state(O._graph(6,disjoint))===nothing
    @test O.ordered_perfect_matching(6,disjoint)===nothing
    # An isolated vertex fails before the search; the cold kernel confirms
    # no perfect matching exists either.
    @test O._fractional_singleton_state(O._graph(4,[(1,2,0)]))===nothing
    @test O.ordered_perfect_matching(4,[(1,2,0)])===nothing
end

# The oracle-generated 3x3 structured triangulation is a seeded two-cycle
# graph whose optimum admits multiple members. The oracle deterministically
# selects one member; production deterministically selects another at the
# identical certified cost. The oracle member is not a distinct convention:
# a single input-edge permutation of the same literal multiset reaches it,
# so both outputs sit inside the same ordering-sensitive optimal family —
# the corpus-wide `GENERAL_NATIVE_TIE_IDENTITY_REMAINS_OPEN` class rather
# than a recoverable chronology divergence.
@testset "Seeded-cycle equal-cost ties stay inside the ordering family" begin
    edges=[(1,2,411),(4,17,513),(4,18,513),(3,15,513),(3,16,513),
        (6,24,536),(5,26,536),(5,23,536),(6,25,536),(2,22,580),
        (1,21,580),(1,19,580),(2,20,580),(7,9,581),(8,10,581),
        (11,13,581),(12,14,581),(13,24,614),(8,23,614),(7,25,614),
        (14,26,614),(9,15,679),(10,18,679),(12,16,679),(11,17,679),
        (19,26,710),(20,24,710),(22,25,710),(21,23,710),(16,19,785),
        (17,20,785),(18,21,785),(15,22,785),(13,11,100000),(7,9,100000),
        (12,14,100000),(8,10,100000),(13,6,100000),(7,6,100000),
        (9,3,100000),(12,3,100000),(14,5,100000),(8,5,100000),
        (10,4,100000),(11,4,100000)]
    oracle_pairs=Set([(1,21),(2,22),(3,16),(4,18),(5,23),(6,24),(7,25),
        (8,10),(9,15),(11,13),(12,14),(17,20),(19,26)])
    selected_pairs(list,result)=Set(Tuple(sort!([list[i][1],list[i][2]]))
        for i in result.selected)

    # The fractional state retains two half-cycles and hands an integral
    # seed to the continuation; the deterministic output is certified
    # optimal at the same cost as the recorded oracle member.
    graph=O._graph(26,edges)
    @test !O._bipartite(graph)
    state=O._fractional_singleton_state(graph)
    @test state!==nothing &&
        count(id->state.cycle_owner[first(state.cycles[id][1])]==id,
            eachindex(state.cycles))==2
    result=O.ordered_perfect_matching(26,edges)
    @test result!==nothing
    @test sum(Int128(edges[i][3]) for i in result.selected)==7789
    @test matching_objective(result.mate,edges)[2]==7789
    # The recorded oracle member is itself a valid selection of the same
    # literal edge multiset at the identical optimum.
    @test all(any(minmax(edges[i][1],edges[i][2])==pair for i in
        eachindex(edges)) for pair in oracle_pairs)
    @test sum(minimum(edges[i][3] for i in eachindex(edges) if
        minmax(edges[i][1],edges[i][2])==pair) for pair in oracle_pairs)==7789

    # Permuting only the input edge order — the same literal multiset —
    # reaches the oracle member, so the divergence is ordering chronology
    # inside one optimal family, not an algorithmic gap.
    rng=Random.MersenneTwister(12345)
    @test any(1:20) do _
        permuted=edges[randperm(rng,length(edges))]
        selected_pairs(permuted,
            O.ordered_perfect_matching(26,permuted))==oracle_pairs
    end
    # Determinism is preserved: the unpermuted run and its clone agree.
    @test O.ordered_perfect_matching(26,Tuple(edges)).selected==result.selected
end
end
