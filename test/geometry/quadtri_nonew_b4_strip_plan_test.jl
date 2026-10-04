module QuadTriNoNewB4StripPlanTests

using Test,Tessella
const M=Tessella.Model
const Q=Rational{BigInt}
const FACES=Dict(4=>((1,3,2),(1,2,4),(2,3,4),(3,1,4)),
    5=>((1,4,3,2),(5,6,7,8),(1,2,6,5),(2,3,7,6),(3,4,8,7),(4,1,5,8)),
    7=>((1,4,3,2),(1,2,5),(2,3,5),(3,4,5),(4,1,5)))
key(face)=Tuple(sort!(collect(face)))
cycle(face)=minimum(Tuple(circshift(collect(face),k)) for k in 0:length(face)-1)
det(a,b,c)=a[1]*(b[2]*c[3]-b[3]*c[2])-a[2]*(b[1]*c[3]-b[3]*c[1])+a[3]*(b[1]*c[2]-b[2]*c[1])

# Integer-coordinate regular disks provide certified source-data inputs for
# the private planner. Native sampling and Source admission are independently
# exercised by the public geometry suite; this fixture substitutes no native
# sample fractions in those tests.
function source(cells;transpose=false,permutation=1,rotation=0)
    count=cells;vertices=2(count+1)
    coords=[(Float64(i),Float64(j),0.) for j in 0:1 for i in 0:count]
    nodes=[(Int32(i),Int32(i+1),Int32(i+count+2),Int32(i+count+1)) for i in 1:count]
    chains=(Int32.(1:count+1),Int32[count+1,vertices],
        Int32.(vertices:-1:count+2),Int32[count+2,1])
    if transpose
        coords=[(p[2],p[1],p[3]) for p in coords]
        nodes=[Tuple(reverse(c)) for c in nodes]
        chains=(reverse(chains[4]),reverse(chains[3]),reverse(chains[2]),reverse(chains[1]))
    end
    ordinal=permutation==1 ? collect(1:vertices) : permutation==2 ? collect(vertices:-1:1) :
        permutation==3 ? circshift(collect(1:vertices),1) :
        vcat(collect(1:2:vertices),collect(2:2:vertices))
    inverse=invperm(ordinal);coords=coords[ordinal]
    nodes=[Tuple(Int32.(inverse[collect(c)])) for c in nodes]
    chains=Tuple(Int32.(inverse[chain]) for chain in chains)
    order=iseven(permutation) ? collect(count:-1:1) : circshift(collect(1:count),1)
    nodes=[Tuple(circshift(collect(nodes[order[i]]),mod(rotation+i-1,4))) for i in 1:count]
    edges=NTuple{2,Int32}[];owners=NTuple{2,Int32}[];sides=NTuple{2,UInt8}[]
    lookup=Dict{NTuple{2,Int32},Int32}();indices=NTuple{4,Int32}[]
    for (number,c) in enumerate(nodes)
        localindices=Int32[]
        for side in 1:4
            edge=minmax(c[side],c[mod1(side+1,4)])
            index=get(lookup,edge,Int32(0))
            if index==0
                push!(edges,edge);index=Int32(length(edges));lookup[edge]=index
                push!(owners,(Int32(number),Int32(0)));push!(sides,(UInt8(side),UInt8(0)))
            else
                owners[index]=(owners[index][1],Int32(number))
                sides[index]=(sides[index][1],UInt8(side))
            end
            push!(localindices,index)
        end
        push!(indices,Tuple(localindices))
    end
    edgecurve=zeros(UInt8,length(edges));direction=zeros(Int8,length(edges))
    for curve in 1:4,i in 1:length(chains[curve])-1
        a,b=chains[curve][i],chains[curve][i+1];index=lookup[minmax(a,b)]
        edgecurve[index]=UInt8(curve);direction[index]=Int8(a<b ? 1 : -1)
    end
    boundary=vcat((chain[1:end-1] for chain in chains)...)
    grid=transpose ? (1,count) : (count,1)
    M._ExtrudeNoNewRectGridSource(nodes,coords,boundary,(11,17,23,29),chains,
        Tuple(length.(chains)),edges,indices,owners,sides,edgecurve,direction,
        trues(vertices),fill(UInt8(15),count),fill(UInt8(4),count),grid,Int8(3),Int8(1),Int8(1))
end

# Build the independently specified B4 terminal fan: five actual Quad-faced
# pyramids plus two top-face tetrahedra, preceded by whole Hex macros. Verify
# its actual typed complex, not only a formula copied from the planner.
function independent_pattern(src,intervals,states)
    vertices=length(src.source_coordinates);coords=NTuple{3,Q}[]
    for layer in 0:intervals,p in src.source_coordinates
        push!(coords,(Q(p[1]),Q(p[2]),Q(layer)))
    end
    emitted=Tuple{Int,Tuple}[]
    for (parent,cell) in enumerate(src.source_cells)
        for layer in 1:intervals
            corners=Tuple(vcat(collect(Int.(cell)).+(layer-1)*vertices,
                collect(Int.(cell)).+layer*vertices))
            if layer<intervals
                push!(emitted,(5,corners));continue
            end
            center=ntuple(d->sum(coords[n][d] for n in corners)/8,3)
            push!(coords,center);apex=length(coords)
            for face in (FACES[5][1],FACES[5][3:end]...)
                base=Tuple(corners[k] for k in reverse(face))
                push!(emitted,(7,(base...,apex)))
            end
            top=corners[5:8]
            triangles=states[parent]==1 ? ((top[1],top[2],top[3]),(top[1],top[3],top[4])) :
                ((top[2],top[3],top[4]),(top[2],top[4],top[1]))
            for tri in triangles;push!(emitted,(4,(tri[2],tri[1],tri[3],apex)));end
        end
    end
    faces=Dict{Tuple,Vector{Tuple}}();families=zeros(Int,4);volume=Q(0)
    for (msh,nodes) in emitted
        families[msh==4 ? 1 : msh==5 ? 2 : msh==6 ? 3 : 4]+=1
        p=Tuple(coords[n] for n in nodes)
        if msh==4
            jac=det(p[2].-p[1],p[3].-p[1],p[4].-p[1])
            @test jac>0;volume+=jac/6
        elseif msh==5
            @test all(p[k+4].-p[k]==(Q(0),Q(0),Q(1)) for k in 1:4)
            jac=det(p[2].-p[1],p[4].-p[1],p[5].-p[1])
            @test jac>0;volume+=jac
        else
            turns=[det(p[mod1(k+1,4)].-p[k],p[mod1(k-1,4)].-p[k],p[5].-p[k]) for k in 1:4]
            @test all(>(0),turns);volume+=sum(turns)/12
        end
        for pattern in FACES[msh]
            face=Tuple(nodes[k] for k in pattern);push!(get!(faces,key(face),Tuple[]),face)
        end
    end
    for uses in values(faces)
        @test length(uses) in (1,2)
        length(uses)==2 && @test cycle(uses[1])==cycle(reverse(uses[2]))
    end
    @test volume==length(src.source_cells)*intervals
    @test length(coords)==vertices*(intervals+1)+length(src.source_cells)
    Tuple(families),length(faces)
end

@testset "B4 planner original ordinal caps and independent fan counts" begin
    for count in (5,7,9),transpose in (false,true),permutation in 1:4,rotation in 0:3
        src=source(count;transpose,permutation,rotation)
        before=deepcopy(src)
        for intervals in (1,2,3)
            levels=Float64.(collect(0:intervals)./intervals)
            refs=[(Int32(1),Int32(k)) for k in 1:intervals]
            catalog=M._extrude_nonew_b4_strip_plan(src,levels,refs,true,"B4 planner independent proof")
            @test catalog isa M._ExtrudeNoNewDynamicGridCatalog
            @test catalog.source===src && catalog.levels===levels && catalog.layer_refs===refs
            @test catalog.recombined && catalog.diagonal_count==count
            for (cell,state) in zip(src.source_cells,catalog.top_states)
                selected=state==1 ? (cell[1],cell[3]) : (cell[2],cell[4])
                @test minimum(cell) in selected
                for shift in 0:3
                    shifted=Tuple(circshift(collect(cell),shift))
                    other=M._extrude_nonew_b4_strip_top_state(shifted)
                    @test Set(selected)==Set(other==1 ? (shifted[1],shifted[3]) : (shifted[2],shifted[4]))
                end
            end
            expected,faces=independent_pattern(src,intervals,catalog.top_states)
            @test catalog.cell_counts==expected && catalog.face_capacity==faces
            repeat=M._extrude_nonew_b4_strip_plan(src,levels,refs,true,"B4 deterministic repeat")
            @test repeat.top_states==catalog.top_states && repeat.cell_counts==catalog.cell_counts
        end
        @test all(getfield(src,k)==getfield(before,k) for k in fieldnames(typeof(src)))
    end
end

@testset "B4 private source admission and checked planner failures" begin
    for count in (5,7,9),shape in ((1,count),(count,1))
        @test M._extrude_nonew_rect_grid_sizes(shape...,"B4 admission";mode=:b4_strip)==
            (2count+2,count,3count+1,2count+2)
        @test_throws r"at least two actual quadrangles" M._extrude_nonew_rect_grid_sizes(shape...,"default rectangular")
    end
    @test M._extrude_nonew_rect_grid_sizes(2,3,"unchanged default")==
        M._extrude_nonew_rect_grid_sizes(2,3,"explicit default";mode=:rectangular)==(12,6,17,10)
    for shape in ((1,4),(4,1),(2,3),(0,5),(-1,5))
        @test_throws r"exactly one direction" M._extrude_nonew_rect_grid_sizes(shape...,"wrong B4 shape";mode=:b4_strip)
    end
    @test_throws r"overflows Int" M._extrude_nonew_rect_grid_sizes(1,typemax(Int),"checked source";mode=:b4_strip)
    src=source(5);levels=[0.,.125,.25,1.];refs=[(Int32(1),Int32(1)),(Int32(1),Int32(2)),(Int32(2),Int32(1))]
    before=deepcopy(src);copy_levels=copy(levels);copy_refs=copy(refs)
    @test_throws r"free lateral faces" M._extrude_nonew_b4_strip_plan(src,levels,refs,false,"free pending")
    @test_throws r"positive intervals" M._extrude_nonew_b4_strip_counts(src,0,"zero intervals")
    @test_throws r"overflows Int" M._extrude_nonew_b4_strip_counts(src,typemax(Int),"overflow")
    @test_throws r"checked node or Int32 limit" M._extrude_nonew_b4_strip_counts(src,M._EXTRUDE_NONEW_MAX_NODES,"node budget")
    for invalid in ([0.,1.],[.1,.2,.3,1.],[0.,.25,.25,1.],[0.,NaN,.5,1.],[0.,Inf,.5,1.],[0.,.2,.3,.9])
        @test_throws r"level|normalized" M._extrude_nonew_b4_strip_plan(src,invalid,refs,true,"invalid levels")
    end
    wrong=deepcopy(src);wrong.category[1]=UInt8(3)
    @test_throws r"four actual boundary vertices" M._extrude_nonew_b4_strip_plan(wrong,levels,refs,true,"wrong category")
    wrong=deepcopy(src);wrong.boundary_masks[1]=UInt8(7)
    @test_throws r"four actual boundary vertices" M._extrude_nonew_b4_strip_plan(wrong,levels,refs,true,"wrong mask")
    wrong=deepcopy(src);pop!(wrong.source_cells)
    @test_throws r"dimensions differ" M._extrude_nonew_b4_strip_plan(wrong,levels,refs,true,"missing certified cell")
    @test all(getfield(src,k)==getfield(before,k) for k in fieldnames(typeof(src)))
    @test levels==copy_levels && refs==copy_refs
end

end
