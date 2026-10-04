if !isdefined(@__MODULE__, :QuadTriNoNewQuadStripCertificates)
    include("quadtri_nonew_quad_strip_certificates.jl")
end

module QuadTriNoNewRectGridCertificates

using Tessella
using Tessella.Elements: MixedMesh, ElementBlock, msh_dimension
using Tessella.MeshTypes: nnodes
using ..QuadTriNoNewQuadStripCertificates
const Strip=QuadTriNoNewQuadStripCertificates
const Patch=Strip.Patch
const Quad=Strip.Quad
const Q=Strip.Q
const FACES=Strip.FACES
const point=Strip.point
const position=Strip.position
const exact=Strip.exact
const key=Strip.key
const cycle=Strip.cycle
const det=Strip.det
const orient=Strip.orient
const blocks=Strip.blocks
const execute=Strip.execute
const public_volume=Strip.public_volume
const packed_digest=Strip.packed_digest
const crosses=Strip.crosses
const map_positive=Strip.map_positive
const cell_volume=Strip.cell_volume
const macro_boundary=Strip.macro_boundary
require(ok,message)=ok || throw(ArgumentError("Rectangular-grid certificate: "*message))

function fixture(name;grid=(2,3),plane=:XY,shape=:unit,winding=1,height=1.,laterals=false,
                 layers=:three,pins=:auto,curve_reverse=0,tags=:dense,
                 offset=(0.,0.,0.),scale=1.,arrangement="Left",progression=1.)
    a,b=Int.(grid);require(a>=2 && b>=2,"fixture needs at least two cells per direction")
    f=Strip.fixture(name;plane,shape,winding,height,laterals,layers,pins,curve_reverse,
                    tags,offset,scale,arrangement,progression)
    xcurves=join((f.curve_tags[1],f.curve_tags[3]),',')
    ycurves=join((f.curve_tags[2],f.curve_tags[4]),',')
    law=progression==1. ? "" : " Using Progression $progression"
    before="Transfinite Curve{$xcurves}=3$law;Transfinite Curve{$ycurves}=2;"
    after="Transfinite Curve{$xcurves}=$(a+1)$law;Transfinite Curve{$ycurves}=$(b+1);"
    require(occursin(before,f.source),"fixture grading declaration missing")
    levels=layers===:nonbinary ? [0.,.2,.2+(1.0-.2)/2,1.] : f.levels
    return merge(f,(;source=replace(f.source,before=>after),levels,grid=(a,b),
        intervals=length(levels)-1,base=(curves=f.curve_tags,surface=f.surface)))
end

function source_complex(source,f)
    a,b=f.grid;v=(a+1)*(b+1);nf=a*b;nb=2(a+b)
    require(nnodes(source)==v && all(isfinite,source.coords),"wrong finite source-node count")
    require(length(Set(point(source,i) for i in 1:v))==v,"source coordinates coincide")
    cells=Tuple(Tuple(Int.(cell)) for block in blocks(source) for cell in eachcol(block.nodes) if block.msh==3)
    require(length(cells)==nf && sum(size(block.nodes,2) for block in blocks(source))==nf,
        "source is not the complete linear Quad4 complex")
    require(all(source.coords[f.axes[3],:].==source.coords[f.axes[3],1]),"source is not exactly axis planar")
    xy=Tuple(Tuple(Q(source.coords[axis,i]) for axis in f.axes[1:2]) for i in 1:v)
    incidence=Dict{Tuple,Vector{Tuple}}();parents=Dict{Tuple,Vector{Int}}();areas=Q[];used=Set{Int}()
    for (parent,cell) in enumerate(cells)
        require(allunique(cell) && all(i->1<=i<=v,cell),"foreign or repeated source node")
        union!(used,cell)
        require(all(sign(orient(xy[cell[k]],xy[cell[mod1(k+1,4)]],xy[cell[mod1(k+2,4)]]))==f.winding for k in 1:4),
            "actual source Quad4 is not coherently strictly convex")
        push!(areas,abs(sum(xy[cell[k]][1]*xy[cell[mod1(k+1,4)]][2]-xy[cell[mod1(k+1,4)]][1]*xy[cell[k]][2] for k in 1:4))/2)
        for k in 1:4
            edge=(cell[k],cell[mod1(k+1,4)]);support=key(edge)
            push!(get!(incidence,support,Tuple[]),edge)
            push!(get!(parents,support,Int[]),parent)
        end
    end
    require(used==Set(1:v) && length(Set(key(c) for c in cells))==nf,"unused or duplicate source cell")
    require(length(incidence)==a*(b+1)+b*(a+1),"wrong rectangular source edge count")
    boundary=Set{Tuple}();shared=Tuple[];adjacency=[Int[] for _ in 1:nf]
    for (edge,uses) in incidence
        require(length(uses) in (1,2),"source edge has excess incidence")
        if length(uses)==1;push!(boundary,only(uses))
        else
            require(uses[1]==reverse(uses[2]),"shared source edge has equal direction")
            push!(shared,edge);left,right=parents[edge]
            push!(adjacency[left],right);push!(adjacency[right],left)
        end
    end
    border=Set(node for edge in boundary for node in edge)
    require(length(boundary)==nb && length(border)==nb && v-length(incidence)+nf==1,"wrong disk Euler/boundary counts")
    following=Dict(first(edge)=>last(edge) for edge in boundary)
    require(length(following)==nb && length(Set(values(following)))==nb,"source boundary branches")
    walk=Int[minimum(border)]
    for _ in 1:nb-1;push!(walk,following[last(walk)]);end
    require(allunique(walk) && following[last(walk)]==first(walk),"source boundary is not one cycle")
    visited=Set([1]);pending=[1]
    while !isempty(pending)
        for neighbor in adjacency[pop!(pending)]
            if !(neighbor in visited);push!(visited,neighbor);push!(pending,neighbor);end
        end
    end
    require(length(visited)==nf,"source cell graph is disconnected")
    # Exact actual edge-contact check includes zero-area contacts, not only
    # positive polygon overlap. Allowed contacts are precisely shared IDs.
    edges=collect(keys(incidence))
    for i in eachindex(edges),j in i+1:length(edges)
        first_edge,second_edge=edges[i],edges[j]
        common=intersect(first_edge,second_edge)
        if isempty(common)
            require(!crosses(xy[first_edge[1]],xy[first_edge[2]],xy[second_edge[1]],xy[second_edge[2]]),
                "nonincident source edges touch or cross")
        else
            shared_node=only(common)
            p=only(node for node in first_edge if node!=shared_node)
            q=only(node for node in second_edge if node!=shared_node)
            if orient(xy[shared_node],xy[p],xy[q])==0
                dp=xy[p].-xy[shared_node];dq=xy[q].-xy[shared_node]
                require(sum(dp[k]*dq[k] for k in 1:2)<0,"incident source edges overlap")
            end
        end
    end
    for i in 1:nf,j in i+1:nf
        common=intersect(cells[i],cells[j])
        require(any(all(begin
                side=f.winding*orient(xy[p],xy[q],xy[node])
                isempty(common) ? side<0 : side<=0
            end for node in other)
            for (own,other) in ((cells[i],cells[j]),(cells[j],cells[i]))
            for (p,q) in ((own[k],own[mod1(k+1,4)]) for k in 1:4)),"source cell interiors overlap")
    end
    polygon=abs(sum(xy[p][1]*xy[q][2]-xy[q][1]*xy[p][2] for (p,q) in boundary))/2
    require(sum(areas)==polygon,"source cells do not cover their actual boundary")
    counts=[count(node->node in border,cell) for cell in cells]
    require(count(==(3),counts)==4 && count(==(2),counts)==2(a+b-4) && count(==(0),counts)==(a-2)*(b-2),
        "source boundary categories are not B3/adjacent-B2/B0")
    require(all(cell->begin
        flags=[node in border for node in cell]
        sum(flags)!=2 || count(k->flags[k] && flags[mod1(k+1,4)],1:4)==1
    end,cells),"opposite-boundary B2 category")
    cad=Tuple(only(i for i in 1:v if point(source,i)==position(p)) for p in f.corners)
    neighbors=Dict(i=>Int[] for i in border)
    for (p,q) in boundary;push!(neighbors[p],q);push!(neighbors[q],p);end
    chains=Tuple(begin
        first_node,last_node=cad[k],cad[mod1(k+1,4)];paths=Tuple[]
        for neighbor in neighbors[first_node]
            path=Int[first_node,neighbor]
            while last(path)!=last_node && !(last(path) in cad)
                push!(path,only(node for node in neighbors[last(path)] if node!=path[end-1]))
                require(length(path)<=nb,"boundary chain did not terminate")
            end
            last(path)==last_node && push!(paths,Tuple(path))
        end
        require(length(paths)==1,"CAD corner pair lacks one actual chain")
        only(paths)
    end for k in 1:4)
    require(all(length(chains[k])==(isodd(k) ? a+1 : b+1) for k in 1:4),"actual chain width differs")
    return (;cells,boundary,shared=Tuple(sort(shared)),boundary_cycle=Tuple(walk),chains,
        xy,areas,area=polygon,border,interior=setdiff(Set(1:v),border),incidence,parents,counts,adjacency)
end

function state_of(group,carrier)
    length(group)==1 && return 0
    require(length(group)==2 && all(length(face)==3 for face in group),"invalid typed macro face split")
    diagonal=intersect(group[1],group[2]);require(length(diagonal)==2,"split face lacks a diagonal")
    diagonal=Set(diagonal)
    return diagonal==Set((carrier[1],carrier[3])) ? 1 :
        diagonal==Set((carrier[2],carrier[4])) ? 2 : error("foreign typed face diagonal")
end

function acyclic(nodes,relations)
    following=Dict(node=>Int[] for node in nodes);degree=Dict(node=>0 for node in nodes)
    for (a,b) in Set(relations);push!(following[a],b);degree[b]+=1;end
    pending=[node for node in nodes if degree[node]==0];seen=0
    while !isempty(pending)
        node=pop!(pending);seen+=1
        for other in following[node];degree[other]-=1;degree[other]==0 && push!(pending,other);end
    end
    return seen==length(nodes)
end

function global_rank_exists(src,top_diagonals,terminal_edge_diagonals)
    relations=Tuple{Int,Int}[];alternatives=Vector{Vector{Vector{Tuple{Int,Int}}}}()
    for (edge,winner) in terminal_edge_diagonals
        loser=only(node for node in edge if node!=winner);push!(relations,(winner,loser))
    end
    for (cell,diagonal) in zip(src.cells,top_diagonals)
        eligible=[node for node in cell if !(node in src.border)]
        selected=[node for node in eligible if node in diagonal]
        require(!isempty(selected),"top diagonal misses every eligible corner")
        push!(alternatives,[[(winner,other) for other in eligible if other!=winner] for winner in selected])
    end
    function search(k,edges)
        acyclic(1:length(src.xy),edges) || return false
        k>length(alternatives) && return true
        return any(additions->search(k+1,vcat(edges,additions)),alternatives[k])
    end
    return search(1,relations)
end

function certify(mesh,f,source;oracle=false,native_top=!oracle)
    src=source_complex(source,f);a,b=f.grid;n=length(f.levels)-1
    v=(a+1)*(b+1);nf=a*b;nb=2(a+b)
    require(nnodes(mesh)==v*(n+1) && all(isfinite,mesh.coords),"added/missing/orphan output nodes")
    intended=[source.coords[f.axes[3],1]+Float64(t)*f.height for t in f.levels]
    planes=oracle ? sort!(unique(mesh.coords[f.axes[3],:]);rev=f.height<0) : intended
    require(length(planes)==n+1 && all(sign.(diff(planes)).==sign(f.height)),"actual planes collapse or reverse")
    require(all(abs(planes[k]-intended[k])<=2e-11 for k in eachindex(planes)),"actual heights differ from declared layers")
    lookup=Dict(point(mesh,node)=>node for node in 1:nnodes(mesh))
    require(length(lookup)==nnodes(mesh),"actual output coordinates coincide")
    grid=[lookup[ntuple(k->k==f.axes[3] ? planes[row] : source.coords[k,col],3)] for row in 1:n+1,col in 1:v]
    require(length(Set(vec(grid)))==v*(n+1),"actual column bijection differs")
    location=Dict(grid[row,col]=>(row,col) for row in 1:n+1,col in 1:v)
    domains=Dict{Tuple,Vector{NamedTuple}}();faces=Dict{Tuple,Vector{Tuple}}();edges=Set{Tuple}()
    used=Set{Int}();families=Dict{Int,Int}();total=zero(Q)
    for block in blocks(mesh),cell in eachcol(block.nodes)
        msh_dimension(block.msh)==3 || continue
        msh=Int(block.msh);require(haskey(FACES,msh) && allunique(cell),"unexpected actual family/node repetition")
        columns=[location[Int(node)][2] for node in cell];rows=[location[Int(node)][1] for node in cell]
        layer=minimum(rows);require(maximum(rows)==layer+1,"cell crosses nonadjacent planes")
        candidates=[parent for parent in 1:nf if Set(columns)⊆Set(src.cells[parent])]
        require(length(candidates)==1,"cell has ambiguous or foreign macro membership")
        parent=only(candidates);p=Tuple(exact(point(mesh,node)) for node in cell)
        require(map_positive(p,msh),"nonpositive whole P1 reference map")
        volume=cell_volume(p,msh);require(volume>0,"nonpositive exact actual-map volume")
        total+=volume;union!(used,Int.(cell));families[msh]=get(families,msh,0)+1
        push!(get!(domains,(parent,layer),NamedTuple[]),(;p,msh,volume,nodes=Tuple(Int.(cell))))
        for pattern in FACES[msh]
            face=Tuple(Int(cell[k]) for k in pattern);push!(get!(faces,key(face),Tuple[]),face)
            for k in eachindex(face);push!(edges,key((face[k],face[mod1(k+1,length(face))])));end
        end
    end
    require(used==Set(1:nnodes(mesh)),"unused actual output node")
    if f.laterals
        expected=Dict(4=>4nf-2nb,7=>nf+nb);n>1 && (expected[5]=nf*(n-1))
        require(families==expected,"wrong B3/B2/B0 recombined families")
    end
    states=Dict{Tuple,NTuple{6,Int}}();macro_groups=Dict{Tuple,Vector{Vector{Tuple}}}()
    for parent in 1:nf,layer in 1:n
        cells=domains[(parent,layer)]
        require(sum(cell.volume for cell in cells)==src.areas[parent]*abs(Q(planes[layer+1])-Q(planes[layer])),
            "actual maps do not exactly partition their macro")
        vertices=Tuple(grid[row,col] for row in (layer,layer+1) for col in src.cells[parent])
        groups,_=macro_boundary(cells,vertices,mesh)
        # Independent face order matches the geometric reference cycles,
        # regardless of any production template index/selected family order.
        logical=((1,2,6,5),(2,3,7,6),(3,4,8,7),(4,1,5,8),(1,2,3,4),(5,6,7,8))
        ordered=[only(group for group in groups if !isempty(group) && Set(first(group))⊆Set(vertices[k] for k in pattern)) for pattern in logical]
        states[(parent,layer)]=ntuple(k->state_of(ordered[k],Tuple(vertices[j] for j in logical[k])),6)
        macro_groups[(parent,layer)]=ordered
    end
    boundary=Dict{Tuple,Tuple}();internal=Dict{Tuple,Vector{Tuple}}()
    for (face,uses) in faces
        require(length(uses) in (1,2),"typed face incidence exceeds two")
        if length(uses)==2
            require(cycle(uses[1])==cycle(reverse(uses[2])),"shared outward cycles agree")
            internal[face]=uses
        else;boundary[face]=only(uses)
        end
    end
    require(length(boundary)==(f.laterals ? nb*n+3nf : 2nb*n+3nf),"wrong exterior typed face count")
    ncell=sum(values(families));require(length(used)-length(edges)+length(faces)-ncell==1,"actual complex is not a ball")
    top_diagonals=Vector{Set{Int}}()
    for (parent,cell) in enumerate(src.cells)
        require(key(grid[1,node] for node in cell) in keys(boundary),"source Quad4 is missing")
        upper=states[(parent,n)][6];require(upper in (1,2),"top cap is not two triangles")
        diagonal=Set(upper==1 ? (cell[1],cell[3]) : (cell[2],cell[4]));push!(top_diagonals,diagonal)
        eligible=[node for node in cell if !(node in src.border)]
        if native_top;require(minimum(eligible) in diagonal,"native cap misses its least original eligible source ID");end
        for layer in 1:n-1
            require(states[(parent,layer)][6]==states[(parent,layer+1)][5],"intermediate cap carry differs")
        end
    end
    terminal_edge_diagonals=Dict{Tuple,Int}()
    for edge in src.shared
        left=first(src.parents[edge]);cell=src.cells[left]
        side=only(k for k in 1:4 if key((cell[k],cell[mod1(k+1,4)]))==edge)
        p,q=cell[side],cell[mod1(side+1,4)]
        for layer in 1:n
            state=states[(left,layer)][side]
            boundary_count=count(node->node in src.border,edge)
            if f.laterals && layer<n;require(state==0,"recombined lower interval is not a whole lateral")
            elseif boundary_count==1
                require(state==(p in src.border ? 1 : 2),"boundary/interior lateral points toward the wrong upper endpoint")
            elseif layer<n;require(state==0,"nonterminal interior/interior lateral split")
            else
                require(state in (1,2),"terminal interior/interior lateral is whole")
                winner=state==1 ? q : p;terminal_edge_diagonals[edge]=winner
                native_top && require(winner==minimum(edge),"terminal lateral misses the least upper source ID")
            end
        end
    end
    require(global_rank_exists(src,top_diagonals,terminal_edge_diagonals),"no single global vertex rank explains cap/shared-face choices")
    require(total==src.area*abs(Q(last(planes))-Q(first(planes))),"represented swept volume differs")
    return (;source=src,grid,heights=planes,boundary,internal,total,ncell,families,domains,
        states,macro_groups,top_diagonals,terminal_edge_diagonals)
end

include("quadtri_nonew_rect_grid_p2_certificates.jl")

end
