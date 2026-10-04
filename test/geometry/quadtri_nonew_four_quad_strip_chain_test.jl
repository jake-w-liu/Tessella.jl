module QuadTriNoNewFourQuadStripChainTests

using Test,Tessella
const M=Tessella.Model

# The immutable 315-record domain has independent full-hex geometry tests.
# Treat its faces as input data here: enumerate real compatible pairs and join
# their middle interface, without a production selector or dynamic program.
decode(state)=ntuple(number->((state-1)÷3^(number-1))%3,4)
encode(states)=1+states[1]+3states[2]+9states[3]+27states[4]
reverse_state(state)=state==0 ? 0 : 3-state

function incidence(order)
    # Match first visitation in the original source-cell/side order. Each pair
    # itself is sorted by the original stored cell indices, as in the catalog.
    edge_order=Int[]
    for original in order,side in 1:4
        edge=side==2 && original<4 ? original : side==4 && original>1 ? original-1 : 0
        edge>0 && !(edge in edge_order) && push!(edge_order,edge)
    end
    shared_cells=ntuple(3) do index
        edge=edge_order[index]
        a=findfirst(==(edge),order);b=findfirst(==(edge+1),order)
        (UInt8(min(a,b)),UInt8(max(a,b)))
    end
    shared_sides=ntuple(3) do index
        edge=edge_order[index]
        ntuple(2) do position
            original=order[Int(shared_cells[index][position])]
            UInt8(original==edge ? 2 : 4)
        end
    end
    exterior=ntuple(4) do number
        [side for side in 1:4 if !any(shared_cells[edge][position]==number &&
            shared_sides[edge][position]==side for edge in 1:3 for position in 1:2)]
    end
    return shared_cells,shared_sides,exterior
end

function canonical_candidates(preferred)
    exterior=((1,3,4),(1,3),(1,3),(1,2,3))
    templates=M._extrude_nonew_templates()
    return ntuple(4) do number
        groups=[Tuple{NTuple{6,UInt8},Int}[] for _ in 1:9]
        for template in templates
            all(side->template.faces[side]!=0,exterior[number]) || continue
            cost=sum(side->Int(template.faces[side]!=preferred[number][side]),exterior[number])
            key=1+Int(template.faces[5])+3Int(template.faces[6])
            push!(groups[key],(template.faces,cost))
        end
        groups
    end
end

function pair_costs(first,second;middle_on_second)
    costs=fill(typemax(Int),9,9,3)
    for cap_a in 1:9,cap_b in 1:9,a in first[cap_a],b in second[cap_b]
        a[1][2]==reverse_state(b[1][4]) || continue
        middle=middle_on_second ? b[1][2] : a[1][4]
        costs[cap_a,cap_b,Int(middle)+1]=min(costs[cap_a,cap_b,Int(middle)+1],a[2]+b[2])
    end
    return costs
end

function canonical_costs(preferred)
    candidates=canonical_candidates(preferred)
    left=pair_costs(candidates[1],candidates[2];middle_on_second=true)
    right=pair_costs(candidates[3],candidates[4];middle_on_second=false)
    costs=fill(typemax(Int),81,81)
    for upper in 1:81,lower in 1:81
        lo=decode(lower);hi=decode(upper)
        keys=ntuple(number->1+lo[number]+3hi[number],4)
        for state in 0:2
            first=left[keys[1],keys[2],state+1]
            second=right[keys[3],keys[4],reverse_state(state)+1]
            first==typemax(Int) || second==typemax(Int) ||
                (costs[lower,upper]=min(costs[lower,upper],first+second))
        end
    end
    return costs
end

function permuted_costs(canonical,order,top)
    # A cap is indexed in original source-cell order, not along the dual path.
    inverse=ntuple(number->findfirst(==(number),order),4)
    indices=[encode(ntuple(number->decode(state)[inverse[number]],4)) for state in 1:81]
    costs=canonical[indices,indices]
    alignment=[sum(Int(actual!=wanted) for (actual,wanted) in zip(decode(upper),top))
               for upper in 1:81]
    return costs,alignment,indices
end

function exhaustive_path_score(costs,alignment,intervals,final)
    best=(typemax(Int),typemax(Int))
    for intermediate in Iterators.product(ntuple(_->1:81,intervals-1)...)
        lower=1;cost=0;align=0
        for upper in (intermediate...,final)
            cost+=costs[lower,upper];align+=alignment[upper];lower=upper
        end
        best=min(best,(cost,align))
    end
    return best
end

@testset "Four-quad cap chain matches independent paired candidates and short paths" begin
    path_preferences=((0x01,0x00,0x02,0x01),(0x02,0x00,0x01,0x00),
                      (0x01,0x00,0x02,0x00),(0x02,0x01,0x01,0x00))
    path_top=(0x01,0x02,0x01,0x02)
    factories=M._extrude_nonew_templates()
    @test length(factories)==315
    canonical=canonical_costs(path_preferences)
    @test size(canonical)==(81,81)
    @test all(!=(typemax(Int)),canonical)
    orders=[Tuple(order) for order in Iterators.product(ntuple(_->1:4,4)...)
            if allunique(order)]
    @test length(orders)==24
    for order in orders
        shared_cells,shared_sides,exterior=incidence(order)
        preferred=ntuple(number->path_preferences[order[number]],4)
        top=ntuple(number->path_top[order[number]],4)
        costs,alignment,indices=permuted_costs(canonical,order,top)
        @test allunique(indices) && indices[1]==1
        @test indices[encode(top)]==encode(path_top)
        for intervals in 1:3
            picks=M._extrude_nonew_four_quad_strip_chain(intervals,shared_cells,shared_sides,
                preferred,top,"independent paired-template exhaustive paths")
            @test size(picks)==(4,intervals)
            @test picks==M._extrude_nonew_four_quad_strip_chain(intervals,shared_cells,shared_sides,
                preferred,top,"deterministic repeat")
            @test all(index->1<=index<=length(factories),picks)
            lower=(0,0,0,0);actual_cost=0;actual_alignment=0
            for interval in 1:intervals
                templates=ntuple(number->factories[Int(picks[number,interval])],4)
                @test ntuple(number->Int(templates[number].faces[5]),4)==lower
                for edge in 1:3
                    a=templates[Int(shared_cells[edge][1])].faces[Int(shared_sides[edge][1])]
                    b=templates[Int(shared_cells[edge][2])].faces[Int(shared_sides[edge][2])]
                    @test a==reverse_state(b)
                end
                for number in 1:4,side in exterior[number]
                    @test templates[number].faces[side]!=0
                    actual_cost+=templates[number].faces[side]!=preferred[number][side]
                end
                lower=ntuple(number->Int(templates[number].faces[6]),4)
                actual_alignment+=sum(Int(actual!=wanted) for (actual,wanted) in zip(lower,top))
            end
            @test lower==top
            @test (actual_cost,actual_alignment)==
                exhaustive_path_score(costs,alignment,intervals,encode(top))
        end
    end
    @test_throws r"positive intervals" M._extrude_nonew_four_quad_strip_chain(0,
        incidence((1,2,3,4))[1:2]...,path_preferences,path_top,"invalid interval count")
end

end
