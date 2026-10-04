module QuadTriNoNewThreeQuadStripChainTests

using Test,Tessella
const M=Tessella.Model

# The immutable factory domain is independently checked by the full-hex tests.
# Enumerate actual candidate triples here, then enumerate every short cap path;
# this oracle does not use the production macro relations or dynamic program.
decode(state)=((state-1)%3,((state-1)÷3)%3,(state-1)÷9)
encode(states)=1+states[1]+3states[2]+9states[3]
reverse_state(state)=state==0 ? 0 : 3-state

function incidence(order)
    shared_cells=ntuple(2) do edge
        a=findfirst(==(edge),order);b=findfirst(==(edge+1),order)
        (UInt8(min(a,b)),UInt8(max(a,b)))
    end
    shared_sides=ntuple(2) do edge
        ntuple(2) do incidence
            original=order[Int(shared_cells[edge][incidence])]
            UInt8(original==edge ? 2 : 4)
        end
    end
    return shared_cells,shared_sides
end

function oracle_costs(shared_cells,shared_sides,preferred,top)
    templates=M._extrude_nonew_templates()
    exterior=ntuple(3) do number
        [side for side in 1:4 if !any(shared_cells[edge][i]==number &&
            shared_sides[edge][i]==side for edge in 1:2 for i in 1:2)]
    end
    candidates=ntuple(3) do number
        groups=[Tuple{NTuple{6,UInt8},Int}[] for _ in 1:9]
        for template in templates
            all(side->template.faces[side]!=0,exterior[number]) || continue
            cost=sum(side->Int(template.faces[side]!=preferred[number][side]),exterior[number])
            push!(groups[1+Int(template.faces[5])+3Int(template.faces[6])],(template.faces,cost))
        end
        groups
    end
    changes=fill(typemax(Int),27,27)
    for upper in 1:27,lower in 1:27
        lo=decode(lower);hi=decode(upper)
        choices=ntuple(n->candidates[n][1+lo[n]+3hi[n]],3)
        for a in choices[1],b in choices[2],c in choices[3]
            triple=(a,b,c)
            all(edge->triple[Int(shared_cells[edge][1])][1][Int(shared_sides[edge][1])]==
                reverse_state(triple[Int(shared_cells[edge][2])][1][Int(shared_sides[edge][2])]),1:2) || continue
            changes[lower,upper]=min(changes[lower,upper],a[2]+b[2]+c[2])
        end
    end
    alignment=[sum(Int(state!=wanted) for (state,wanted) in zip(decode(upper),top)) for upper in 1:27]
    return changes,alignment,exterior
end

function exhaustive_path_score(changes,alignment,intervals,final)
    best=(typemax(Int),typemax(Int))
    for intermediate in Iterators.product(ntuple(_->1:27,intervals-1)...)
        lower=1;cost=0;align=0
        for upper in (intermediate...,final)
            cost+=changes[lower,upper];align+=alignment[upper];lower=upper
        end
        best=min(best,(cost,align))
    end
    return best
end

@testset "Three-quad cap chain matches exhaustive short paths" begin
    path_preferences=((0x01,0x00,0x02,0x01),(0x02,0x00,0x01,0x00),(0x01,0x02,0x01,0x00))
    path_top=(0x01,0x02,0x01)
    for order in ((1,2,3),(1,3,2),(2,1,3),(2,3,1),(3,1,2),(3,2,1))
        shared_cells,shared_sides=incidence(order)
        preferred=ntuple(n->path_preferences[order[n]],3)
        top=ntuple(n->path_top[order[n]],3)
        changes,alignment,exterior=oracle_costs(shared_cells,shared_sides,preferred,top)
        @test all(!=(typemax(Int)),changes)
        for intervals in 1:4
            picks=M._extrude_nonew_three_quad_strip_chain(intervals,shared_cells,shared_sides,
                preferred,top,"independent exhaustive path")
            @test size(picks)==(3,intervals)
            @test picks==M._extrude_nonew_three_quad_strip_chain(intervals,shared_cells,shared_sides,
                preferred,top,"deterministic repeat")
            lower=(0,0,0);actual_changes=0;actual_alignment=0
            for interval in 1:intervals
                templates=ntuple(n->M._extrude_nonew_templates()[Int(picks[n,interval])],3)
                @test ntuple(n->Int(templates[n].faces[5]),3)==lower
                for edge in 1:2
                    a=templates[Int(shared_cells[edge][1])].faces[Int(shared_sides[edge][1])]
                    b=templates[Int(shared_cells[edge][2])].faces[Int(shared_sides[edge][2])]
                    @test a==reverse_state(b)
                end
                for n in 1:3,side in exterior[n]
                    @test templates[n].faces[side]!=0
                    actual_changes+=templates[n].faces[side]!=preferred[n][side]
                end
                lower=ntuple(n->Int(templates[n].faces[6]),3)
                actual_alignment+=sum(Int(state!=wanted) for (state,wanted) in zip(lower,top))
            end
            @test lower==top
            @test (actual_changes,actual_alignment)==
                exhaustive_path_score(changes,alignment,intervals,encode(top))
        end
    end
end

end
