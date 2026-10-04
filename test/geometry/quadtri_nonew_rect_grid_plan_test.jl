module RectGridPlanIndependentTests
using Tessella, Test, Random
const M=Tessella.Model

# The immutable factories are input data. Their typed cells and all cube masks
# are independently certified by quadtri_nonew_templates_test.jl. This test
# does not invoke the production mask lookup, cap, or physical-state helpers.
const EXPORTED=[Dict{String,Any}("index"=>index,"kind"=>Int(template.kind),
    "states"=>Int.(template.faces),
    "cells"=>[Dict{String,Any}("type"=>Int(cell.msh),
        "nodes"=>Int.(cell.nodes[1:(cell.msh==4 ? 4 : cell.msh==5 ? 8 : cell.msh==6 ? 6 : 5)]))
        for cell in template.cells[1:Int(template.ncells)]])
    for (index,template) in enumerate(M._extrude_nonew_templates())]
const REF=Dict{NTuple{6,UInt8},Any}()
for record in EXPORTED
    states=Tuple(UInt8.(record["states"]))
    previous=get(REF,states,nothing)
    if previous===nothing || (record["kind"],record["index"])<(previous["kind"],previous["index"])
        REF[states]=record
    end
end

const FACEPATTERNS=Dict(4=>((1,3,2),(1,2,4),(2,3,4),(3,1,4)),
    5=>((1,4,3,2),(5,6,7,8),(1,2,6,5),(2,3,7,6),(3,4,8,7),(4,1,5,8)),
    6=>((1,3,2),(4,5,6),(1,2,5,4),(2,3,6,5),(3,1,4,6)),
    7=>((1,4,3,2),(1,2,5),(2,3,5),(3,4,5),(4,1,5)))
const MACROFACES=((1,2,6,5),(2,3,7,6),(3,4,8,7),(4,1,5,8),(1,2,3,4),(5,6,7,8))

function fixture(a,b;seed=0,reverse_nodes=false,reverse_cells=false,rotations=true)
    vertices=(a+1)*(b+1)
    mapping=seed==0 ? collect(1:vertices) : randperm(MersenneTwister(seed),vertices)
    reverse_nodes && reverse!(mapping)
    index(i,j)=Int32(mapping[1+i+(a+1)*j])
    coords=Vector{NTuple{3,Float64}}(undef,vertices)
    boundary=falses(vertices)
    for j in 0:b,i in 0:a
        node=index(i,j);coords[node]=(Float64(i)/a,Float64(j)/b,0.0)
        boundary[node]=i==0 || i==a || j==0 || j==b
    end
    chains=([index(i,0) for i in 0:a],[index(a,j) for j in 0:b],
        [index(i,b) for i in a:-1:0],[index(0,j) for j in b:-1:0])
    cycle=reduce(vcat,(chain[1:end-1] for chain in chains))
    native=Dict{NTuple{2,Int32},Tuple{UInt8,Int8}}()
    for (slot,chain) in enumerate(chains),k in 1:length(chain)-1
        u,v=chain[k],chain[k+1]
        native[minmax(u,v)]=(UInt8(slot),Int8(u<v ? 1 : -1))
    end
    cells=NTuple{4,Int32}[(index(i,j),index(i+1,j),index(i+1,j+1),index(i,j+1)) for j in 0:b-1 for i in 0:a-1]
    reverse_cells && reverse!(cells)
    if seed!=0;shuffle!(MersenneTwister(seed+1),cells);end
    if rotations
        for k in eachindex(cells)
            rotation=mod(k+seed,4);original=cells[k]
            cells[k]=ntuple(p->original[mod1(p+rotation,4)],4)
        end
    end
    edges=NTuple{2,Int32}[];edge_cells=NTuple{2,Int32}[];edge_sides=NTuple{2,UInt8}[]
    indices=Vector{NTuple{4,Int32}}(undef,length(cells));seen=Dict{NTuple{2,Int32},Int}()
    for (c,cell) in enumerate(cells)
        cell_edges=zeros(Int32,4)
        for side in 1:4
            key=minmax(cell[side],cell[mod1(side+1,4)])
            edge=get(seen,key,0)
            if edge==0
                push!(edges,key);push!(edge_cells,(Int32(c),Int32(0)));push!(edge_sides,(UInt8(side),UInt8(0)))
                edge=length(edges);seen[key]=edge
            else
                edge_cells[edge]=(edge_cells[edge][1],Int32(c));edge_sides[edge]=(edge_sides[edge][1],UInt8(side))
            end
            cell_edges[side]=Int32(edge)
        end
        indices[c]=(cell_edges[1],cell_edges[2],cell_edges[3],cell_edges[4])
    end
    mask=UInt8[sum(UInt8(boundary[cell[k]])<<(k-1) for k in 1:4) for cell in cells]
    category=UInt8[count_ones(m) for m in mask]
    slots=UInt8[get(native,e,(UInt8(0),Int8(0)))[1] for e in edges]
    directions=Int8[get(native,e,(UInt8(0),Int8(0)))[2] for e in edges]
    return M._ExtrudeNoNewRectGridSource(cells,coords,cycle,(1,2,3,4),chains,
        (a+1,b+1,a+1,b+1),edges,indices,edge_cells,edge_sides,slots,directions,
        boundary,mask,category,(a,b),Int8(3),Int8(1),Int8(1))
end

pair(a,b)=minmax(a,b)
vertex(node,level)=(Int32(node),Int32(level))
function opposite_cycle(first,second)
    length(first)==length(second) || return false
    reverse_second=reverse(second)
    for shift in 0:length(first)-1
        all(first[k]==reverse_second[mod1(k+shift,length(first))] for k in eachindex(first)) && return true
    end
    return false
end
function expected(source,N,recombined)
    counts=zeros(Int,4);face_rows=Dict{Tuple,Vector{Tuple}}();all_diagonals=Set{Tuple}()
    picks=zeros(UInt16,length(source.source_cells),N);top_states=zeros(UInt8,length(source.source_cells))
    carry=Union{Nothing,NTuple{2,Int32}}[nothing for _ in source.source_cells]
    for interval in 1:N
        physical=Set{Tuple}()
        for (edge,(u,v)) in enumerate(source.edges)
            terminal=interval==N
            if source.edge_cells[edge][2]==0
                recombined && continue
                chain=source.boundary_chains[Int(source.edge_curve[edge])]
                position=findfirst(==(u),chain)
                if position<length(chain) && chain[position+1]==v
                    push!(physical,pair(vertex(u,interval),vertex(v,interval+1)))
                else
                    @assert position>1 && chain[position-1]==v
                    push!(physical,pair(vertex(v,interval),vertex(u,interval+1)))
                end
            elseif !recombined || terminal
                ub,vb=source.boundary_vertices[u],source.boundary_vertices[v]
                if ub!=vb
                    boundary,inside=ub ? (u,v) : (v,u)
                    push!(physical,pair(vertex(boundary,interval),vertex(inside,interval+1)))
                elseif terminal
                    low,other=minmax(u,v)
                    push!(physical,pair(vertex(low,interval+1),vertex(other,interval)))
                end
            end
        end
        union!(all_diagonals,physical)
        for (c,cell) in enumerate(source.source_cells)
            boundary=Bool[source.boundary_vertices[n] for n in cell]
            eligible=findall(!,boundary)
            kind=count(boundary)
            upper=if interval==N
                anchor=eligible[argmin(cell[p] for p in eligible)]
                pair(cell[anchor],cell[mod1(anchor+2,4)])
            elseif recombined || kind==0
                nothing
            elseif kind==3
                anchor=only(eligible);pair(cell[anchor],cell[mod1(anchor+2,4)])
            elseif carry[c]!==nothing
                carry[c]
            else
                side=only(p for p in 1:4 if boundary[p] && boundary[mod1(p+1,4)])
                next=mod1(side+1,4)
                forward=pair(vertex(cell[side],interval),vertex(cell[next],interval+1)) in physical
                anchor=forward ? next : side
                pair(cell[anchor],cell[mod1(anchor+2,4)])
            end
            local_diagonals=copy(physical)
            if carry[c]!==nothing
                u,v=carry[c];push!(local_diagonals,pair(vertex(u,interval),vertex(v,interval)))
            end
            if upper!==nothing
                u,v=upper;diagonal=pair(vertex(u,interval+1),vertex(v,interval+1))
                push!(local_diagonals,diagonal);push!(all_diagonals,diagonal)
            end
            corners=(ntuple(k->vertex(cell[k],interval),4)...,ntuple(k->vertex(cell[k],interval+1),4)...)
            faces=ntuple(6) do face
                a,b,d,e=MACROFACES[face]
                ac=pair(corners[a],corners[d]) in local_diagonals
                bd=pair(corners[b],corners[e]) in local_diagonals
                @assert !(ac && bd)
                UInt8(ac ? 1 : bd ? 2 : 0)
            end
            record=REF[faces];picks[c,interval]=UInt16(record["index"])
            interval==N && (top_states[c]=faces[6])
            for emitted in record["cells"]
                family=emitted["type"];counts[family-3]+=1
                nodes=Tuple(corners[k] for k in emitted["nodes"])
                for pattern in FACEPATTERNS[family]
                    cycle=Tuple(nodes[k] for k in pattern);key=Tuple(sort!(collect(cycle)))
                    push!(get!(face_rows,key,Tuple[]),cycle)
                end
            end
            carry[c]=upper
        end
    end
    @test all(rows->length(rows) in (1,2),values(face_rows))
    @test all(rows->length(rows)==1 || opposite_cycle(rows[1],rows[2]),values(face_rows))
    exterior=count(rows->length(rows)==1,values(face_rows))
    border=sum(length(chain)-1 for chain in source.boundary_chains)
    @test exterior==(recombined ? 1 : 2)*border*N+3length(source.source_cells)
    return picks,top_states,Tuple(counts),length(face_rows),length(all_diagonals)
end

@testset "Rectangular planner independent pair-based reference" begin
    for (a,b) in ((2,2),(2,3),(3,2),(3,3),(2,7),(4,3),(5,4)),
            seed in 0:5,N in 1:4,recombined in (false,true)
        source=fixture(a,b;seed,reverse_nodes=isodd(seed),reverse_cells=seed%3==1)
        levels=[Float64(k)/N for k in 0:N]
        refs=NTuple{2,Int32}[(Int32(k<2 ? 0 : 1),Int32(k<2 ? k : k-2)) for k in 0:N-1]
        plan=M._extrude_nonew_rect_grid_plan(source,levels,refs,recombined,"independent test")
        picks,top,counts,faces,diagonals=expected(source,N,recombined)
        @test plan.source===source
        @test plan.template_indices==picks
        @test plan.top_states==top
        @test plan.cell_counts==counts
        @test plan.face_capacity==faces
        @test plan.diagonal_count==diagonals
        for edge in eachindex(source.edges),interval in 1:N
            c=Int(source.edge_cells[edge][1]);s=Int(source.edge_sides[edge][1])
            state=Tuple(UInt8.(EXPORTED[Int(picks[c,interval])]["states"]))[s]
            want=source.source_cells[c][s]==source.edges[edge][1] ? state : state==0 ? UInt8(0) : UInt8(3-state)
            @test M._extrude_nonew_rect_grid_edge_state(plan,edge,interval)==want
        end
    end
end

@testset "Planner layer limits precede detached picks" begin
    source=fixture(2,3)
    @test_throws ArgumentError M._extrude_nonew_rect_grid_plan(source,[0.,1.],NTuple{2,Int32}[],false,"test")
    @test_throws ArgumentError M._extrude_nonew_rect_grid_plan(source,[0.,0.,1.],[(Int32(0),Int32(0)),(Int32(1),Int32(0))],false,"test")
    @test_throws ArgumentError M._extrude_nonew_rect_grid_plan(source,[0.,2.],[(Int32(0),Int32(0))],false,"test")
    @test_throws ArgumentError M._extrude_nonew_rect_grid_plan(source,[0.,NaN,1.],[(Int32(0),Int32(0)),(Int32(1),Int32(0))],false,"test")
    @test_throws ArgumentError M._extrude_nonew_rect_grid_plan(source,[.1,1.],[(Int32(0),Int32(0))],false,"test")
    @test_throws ArgumentError M._extrude_nonew_rect_grid_top_state((Int32(1),Int32(2),Int32(3),Int32(4)),0x01,0x01,"test")
    @test_throws ArgumentError M._extrude_nonew_rect_grid_top_state((Int32(1),Int32(2),Int32(3),Int32(4)),0x05,0x02,"test")
end

end # module
