module OrderedMatchingSeededTests
using Test
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
end
