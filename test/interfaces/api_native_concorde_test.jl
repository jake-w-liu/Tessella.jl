module NativeMatchingContractTestsV1
using Tessella, Test, TOML
const M = Tessella.Model

# These are the original pinned graph inputs and native directed outputs. The
# independent engines must cover every vertex with the recorded integer cost.
# Private fractional duals and contraction records describe the former solver;
# they are parsed and retained in the immutable corpus, not a public contract.
struct GraphProof
    name::String
    count::Int
    edges::Vector{Tuple{Int,Int,Int}}
    matching::Vector{NTuple{4,Int}}
    states::Int
    expands::Int
    cost::Int
end
function read_graph_proofs(path)
    rows = GraphProof[]
    open(path) do io
        while !eof(io)
            header = split(readline(io))
            header[1] == "CASE" || error("invalid pinned graph case")
            n, ne, nm, _, _, expands, ns, cost = parse.(Int, header[3:end])
            edges = Tuple{Int,Int,Int}[]
            for _ in 1:ne
                fields = split(readline(io))
                fields[1] == "E" || error("invalid pinned graph edge")
                push!(edges, Tuple(parse.(Int, fields[2:4])))
            end
            matching = NTuple{4,Int}[]
            for _ in 1:nm
                fields = split(readline(io))
                fields[1] == "M" || error("invalid pinned graph matching")
                push!(matching, Tuple(parse.(Int, fields[2:5])))
            end
            for _ in 1:ns
                fields = split(readline(io))
                fields[1] == "S" || error("invalid preserved solver state")
                for cell in fields[3:end]
                    length(parse.(Int, split(cell, ','))) == 4 ||
                        error("invalid preserved solver state tuple")
                end
            end
            push!(rows, GraphProof(header[2], n, edges, matching, ns, expands, cost))
        end
    end
    return rows
end

# Independently color the undirected graph. Loops cannot cover two vertices.
function is_bipartite(n, edges)
    neighbors = [Int[] for _ in 1:n]
    for (a, b, _) in edges
        a == b && continue
        push!(neighbors[a], b); push!(neighbors[b], a)
    end
    colors = zeros(Int8, n)
    for root in 1:n
        colors[root] == 0 || continue
        colors[root] = 1
        queue = Int[root]
        for v in queue
            for other in neighbors[v]
                if colors[other] == 0
                    colors[other] = -colors[v]; push!(queue, other)
                elseif colors[other] == colors[v]
                    return false
                end
            end
        end
    end
    return true
end

# A subset recurrence supplies an independent exact minimum and counts optimal
# endpoint matchings (capped at two). Equal-minimum parallel identities are
# checked separately before an optimum is called a unique directed result.
function independent_optimum(n, edges)
    costs = Dict{Tuple{Int,Int},Int128}()
    occurrences = Dict{Tuple{Int,Int},Int}()
    for (a, b, input_cost) in edges
        a == b && continue
        key = minmax(a, b)
        previous = get(costs, key, nothing)
        if previous === nothing || input_cost < previous
            costs[key] = Int128(input_cost); occurrences[key] = 1
        elseif input_cost == previous
            occurrences[key] += 1
        end
    end
    memo = Dict{UInt64,Tuple{Int128,Int}}()
    function visit(mask::UInt64)
        mask == 0 && return (Int128(0), 1)
        haskey(memo, mask) && return memo[mask]
        first = trailing_zeros(mask) + 1
        rest = mask & ~(UInt64(1) << (first - 1))
        best = typemax(Int128); optimal_ways = 0
        for other in first+1:n
            bit = UInt64(1) << (other - 1)
            rest & bit == 0 && continue
            edge_cost = get(costs, (first, other), nothing)
            edge_cost === nothing && continue
            sub_value, sub_ways = visit(rest & ~bit)
            sub_ways == 0 && continue
            candidate_cost = sub_value + edge_cost
            if candidate_cost < best
                best = candidate_cost; optimal_ways = sub_ways
            elseif candidate_cost == best
                optimal_ways = min(2, optimal_ways + sub_ways)
            end
        end
        return memo[mask] = (best, optimal_ways)
    end
    total_value, total_ways = visit((UInt64(1) << n) - UInt64(1))
    return total_value, total_ways == 1 && all(==(1), values(occurrences))
end

@testset "Independent matching adapter: preserved pinned graph corpus" begin
    rows = read_graph_proofs(joinpath(@__DIR__, "..", "..", "validation",
        "native_concorde", "native_concorde_integral_literal_v1.txt"))
    @test length(rows) == 5178
    @test count(row -> row.states > 0, rows) == 1550
    @test count(row -> row.expands > 0, rows) == 533
    bipartite_count = 0; independently_checked = 0; unique_count = 0
    observed_general = Dict{String,Any}[]
    for row in rows
        @testset "$(row.name)" begin
            pairs = [(first=a, second=b, cost=cost) for (a,b,cost) in row.edges]
            before = copy(row.edges)
            result = M._native_recombine_matching_result(row.count, pairs,
                Tuple{Int,Int}[], "pinned matching corpus")
            @test result !== nothing
            result === nothing && continue
            @test result.edges == row.edges
            @test row.edges == before
            @test result.physical_count == length(row.edges)
            @test length(result.selected) == row.count ÷ 2
            @test all(1 <= index <= length(row.edges) for index in result.selected)
            @test length(unique(result.selected)) == length(result.selected)
            selected = row.edges[result.selected]
            @test all(a != b for (a,b,_) in selected)
            @test sort!([vertex for (a,b,_) in selected for vertex in (a,b)]) ==
                collect(1:row.count)
            @test sum(Int128(cost) for (_,_,cost) in selected; init=Int128(0)) == row.cost
            @test issorted([b for (_,b,_) in selected])
            expected = [(a,b) for (a,b,_,_) in row.matching]
            actual = [(a,b) for (a,b,_) in selected]
            # Validate every original primary matching record against its input.
            @test all(row.edges[index] == (a,b,cost)
                for (a,b,cost,index) in row.matching)
            @test sum(Int128(cost) for (_,_,cost,_) in row.matching;
                init=Int128(0)) == row.cost
            bipartite = is_bipartite(row.count, row.edges)
            bipartite_count += bipartite
            unique_directed = false
            if row.count <= 14
                value, unique_directed = independent_optimum(row.count, row.edges)
                @test value == row.cost
                independently_checked += 1
                unique_count += unique_directed
            end
            if bipartite || unique_directed
                @test actual == expected
            else
                # Preserve an explicit observation of every unresolved general
                # tie identity. Minimum-cost correctness is a separate claim.
                push!(observed_general, Dict("case" => row.name,
                    "same_directed_output" => actual == expected,
                    "primary_edge_indices" => [entry[4] for entry in row.matching],
                    "actual_edge_indices" => result.selected))
            end
        end
    end
    @test bipartite_count == 248
    @test independently_checked == 3934
    @test unique_count == 3324
    println("GENERAL_NATIVE_TIE_IDENTITY_REMAINS_OPEN")
    TOML.print(stdout, Dict("general_identity_observations" => observed_general))
end

# The existing public infeasible controls retain every original Quad cycle,
# constructor label, node payload, and failed-perfect-match retry assertion.
include(joinpath(@__DIR__, "..", "..", "validation", "native_concorde",
    "native_recombine_infeasible_focus_v1.jl"))
end
