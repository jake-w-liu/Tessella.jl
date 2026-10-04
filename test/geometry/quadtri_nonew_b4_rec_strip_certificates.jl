# Independent B4 source disk, emitted-center shell and actual P1/P2 maps.
if !isdefined(@__MODULE__,:QuadTriNoNewRectGridCertificates)
    include("quadtri_nonew_rect_grid_certificates.jl")
end
module QuadTriNoNewB4RecStripCertificates
using Tessella,SHA
using Tessella.Elements:MixedMesh,ElementBlock,msh_dimension
using Tessella.MeshTypes:nnodes
using ..QuadTriNoNewRectGridCertificates
const Rect=QuadTriNoNewRectGridCertificates
const Strip=Rect.Strip
const Quad=Rect.Quad
const Q=Rect.Q
const FACES=Rect.FACES
const point=Rect.point
const position=Rect.position
const exact=Rect.exact
const key=Rect.key
const cycle=Rect.cycle
const det=Rect.det
const orient=Rect.orient
const blocks=Rect.blocks
const crosses=Rect.crosses
const map_positive=Rect.map_positive
const cell_volume=Rect.cell_volume
const separated=Strip.separated
const macro_boundary=Rect.macro_boundary
const execute=Rect.execute
const public_volume=Rect.public_volume
const packed_digest=Rect.packed_digest
const quadratic_map_certificate=Rect.quadratic_map_certificate
require(ok,message)=ok || throw(ArgumentError("B4 recombined-strip certificate: "*message))

function fixture(name;strip_length=5,plane=:XY,shape=:unit,winding=1,height=1.,
        laterals=true,layers=:three,pins=:auto,curve_reverse=0,tags=:dense,
        offset=(0.,0.,0.),scale=1.,arrangement="Left",law="Progression",
        coefficient=1.,long_pair=(1,3))
    require(strip_length>=5 && long_pair in ((1,3),(2,4)),"fixture needs a valid long opposite pair")
    # Start from the independently used two-strip recipe builder; replace the
    # one explicit native grading declaration, not its sampled coordinates.
    f=Strip.fixture(name;plane,shape,winding,height,laterals,layers,pins,curve_reverse,
        tags,offset,scale,arrangement,progression=1.,long_pair)
    curves=join((f.curve_tags[k] for k in long_pair),',')
    before="Transfinite Curve{$curves}=3;"
    require(occursin(before,f.source),"missing explicit long-chain control")
    require(law in ("Progression","Bump","Beta"),"unrecognized native grading law")
    grading=law=="Progression" && coefficient==1. ? "" : " Using $law $coefficient"
    text=replace(f.source,before=>"Transfinite Curve{$curves}=$(strip_length+1)$grading;")
    levels=layers===:nonbinary ? [0.,.2,.2+(1.0-.2)/2,1.] : f.levels
    grid=long_pair==(1,3) ? (strip_length,1) : (1,strip_length)
    return merge(f,(;source=text,levels,grid,strip_length,law,coefficient,
        intervals=length(levels)-1,base=(curves=f.curve_tags,surface=f.surface)))
end

# Saved top diagonals must admit one common total rank; the rank is never
# guessed from public Gmsh tags. This small test-only search is not a planner.
function global_rank_exists(cells,diagonals)
    alternatives=[[[ (winner,other) for other in cell if other!=winner ]
        for winner in cell if winner in diagonal] for (cell,diagonal) in zip(cells,diagonals)]
    nodes=unique(vcat(collect.(cells)...))
    function visit(k,relations)
        Rect.acyclic(nodes,relations) || return false
        k>length(alternatives) && return true
        return any(more->visit(k+1,vcat(relations,more)),alternatives[k])
    end
    return visit(1,Tuple{Int,Int}[])
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
    require(min(a,b)==1 && max(a,b)>=5 && all(==(4),counts),
        "source is not a complete all-boundary B4 strip")
    require(length(shared)==nf-1 && sort!(length.(adjacency))==vcat([1,1],fill(2,nf-2)),
        "source dual is not one path")
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

function certify(mesh,f,source;oracle=false,native_top=!oracle,check_separation=true)
    src=source_complex(source,f);n=length(f.levels)-1
    m=prod(f.grid);v=2(m+1);require(f.laterals,"B4 draft only certifies recombined laterals")
    expected_nodes=v*(n+1)+m
    require(nnodes(mesh)==expected_nodes && all(isfinite,mesh.coords),"wrong primary/centroid node count")
    intended=[source.coords[f.axes[3],1]+Float64(t)*f.height for t in f.levels]
    planes=if oracle
        columns=Set(Tuple(source.coords[a,v] for a in f.axes[1:2]) for v in 1:v)
        sort!(unique(mesh.coords[f.axes[3],v] for v in 1:nnodes(mesh) if
            Tuple(mesh.coords[a,v] for a in f.axes[1:2]) in columns);rev=f.height<0)
    else
        intended
    end
    require(length(planes)==n+1 && all(sign.(diff(planes)).==sign(f.height)),"actual planes collapse or reverse")
    require(all(abs(planes[i]-intended[i])<=2e-11 for i in eachindex(planes)),"actual layer heights differ from their declared levels")
    lookup=Dict(point(mesh,v)=>v for v in 1:nnodes(mesh))
    require(length(lookup)==nnodes(mesh),"actual output contains coincident node coordinates")
    grid=[lookup[ntuple(k->k==f.axes[3] ? planes[row] : source.coords[k,col],3)]
        for row in 1:n+1,col in 1:v]
    primary=Set(vec(grid));centers=sort!(collect(setdiff(Set(1:nnodes(mesh)),primary)))
    require(length(primary)==v*(n+1) && length(centers)==m,"column bijection differs")
    center_parent=Dict{Int,Int}()
    for center in centers
        p=exact(point(mesh,center));xy=Tuple(p[a] for a in f.axes[1:2])
        require(min(Q(planes[end-1]),Q(planes[end]))<p[f.axes[3]]<max(Q(planes[end-1]),Q(planes[end])),
            "actual terminal centroid is not strictly between its planes")
        parents=[i for (i,cell) in enumerate(src.cells) if all(sign(orient(src.xy[cell[k]],src.xy[cell[mod1(k+1,4)]],xy))==f.winding for k in 1:4)]
        require(length(parents)==1,"actual terminal centroid lacks one strict macro owner")
        parent=only(parents)
        mean=ntuple(k->sum(Q(source.coords[k,node]) for node in src.cells[parent])/4,3)
        expected=ntuple(k->k==f.axes[3] ? (Q(planes[end-1])+Q(planes[end]))/2 : mean[k],3)
        require(all(abs(Q(p[k])-expected[k])<=Q(2e-11) for k in 1:3),
            "actual center differs from the exact mean of its eight corners")
        center_parent[center]=parent
    end
    require(Set(values(center_parent))==Set(1:m) ,"terminal centroids do not occupy distinct macros")
    location=Dict(grid[row,col]=>(row,col) for row in 1:n+1,col in 1:v)
    domains=Dict{Tuple,Vector{NamedTuple}}();faces=Dict{Tuple,Vector{Tuple}}();edges=Set{Tuple}()
    used=Set{Int}();families=Dict{Int,Int}();total=zero(Q)
    for block in blocks(mesh),cell in eachcol(block.nodes)
        msh_dimension(block.msh)==3 || continue
        msh=Int(block.msh);require(haskey(FACES,msh) && allunique(cell),"unexpected family or repeated actual node IDs")
        body=[Int(v) for v in cell if !(v in primary)]
        column_ids=[location[Int(v)][2] for v in cell if v in primary]
        levels=[location[Int(v)][1] for v in cell if v in primary]
        layer=isempty(body) ? minimum(levels) : n
        require(isempty(body) ? maximum(levels)==layer+1 :
            all(row->row in (n,n+1),levels),"cell crosses nonadjacent planes")
        parents=[i for i in 1:m if Set(column_ids)⊆Set(src.cells[i]) && all(center_parent[v]==i for v in body)]
        require(length(parents)==1,"actual cell has ambiguous or foreign macro membership")
        parent=only(parents)
        require(isempty(body) || (layer==n && length(body)==1),"body centroid escaped its terminal macro")
        p=Tuple(exact(point(mesh,v)) for v in cell)
        require(map_positive(p,msh),"nonpositive whole P1 reference map")
        volume=cell_volume(p,msh);require(volume>0,"nonpositive exact actual-map volume")
        total+=volume;union!(used,Int.(cell));families[msh]=get(families,msh,0)+1
        push!(get!(domains,(parent,layer),NamedTuple[]),(;p,msh,volume,nodes=Tuple(Int.(cell))))
        for pattern in FACES[msh]
            face=Tuple(Int(cell[k]) for k in pattern);push!(get!(faces,key(face),Tuple[]),face)
            for k in eachindex(face);push!(edges,key((face[k],face[mod1(k+1,length(face))])));end
        end
    end
    require(used==Set(1:nnodes(mesh)),"unused actual primary/centroid node")
    if f.laterals
        expected=Dict(4=>2m,7=>5m);n>1 && (expected[5]=m*(n-1))
        require(families==expected,"wrong recombined terminal-fan families")
    end
    separated_pairs=0
    for parent in 1:m,layer in 1:n
        cells=domains[(parent,layer)]
        require(sum(c.volume for c in cells)==src.areas[parent]*abs(Q(planes[layer+1])-Q(planes[layer])),"actual reference maps do not exactly partition the macro")
        vertices=Tuple(grid[row,col] for row in (layer,layer+1) for col in src.cells[parent])
        groups,reference=macro_boundary(cells,vertices,mesh)
        if layer<n
            require(length(cells)==1 && only(cells).msh==5,"earlier macro is not one whole Hex8")
        else
            center=only(center for (center,owner) in center_parent if owner==parent)
            require(length(cells)==7 && count(cell->cell.msh==4,cells)==2 &&
                count(cell->cell.msh==7,cells)==5 && all(center in cell.nodes for cell in cells),
                "terminal macro is not its real nine-vertex/seven-cell fan")
            radial=Set(key((center,node)) for node in vertices)
            require(radial⊆edges,"terminal fan lacks an actual radial edge")
            for pattern in (reference[1],reference[3:6]...)
                patches=only(group for group in groups if !isempty(group) && Set(first(group))⊆Set(pattern))
                require(length(patches)==1 && length(only(patches))==4,"terminal lower/lateral face is split")
            end
        end
        if check_separation
            for i in 1:length(cells),j in i+1:length(cells)
                # A planar supporting face gives a second separation proof.
                # Warped Pri6 interfaces instead use the degree-one shell
                # above: identical opposite typed maps cancel, all six outer
                # patches cover the convex macro exactly, positive full maps
                # and nonnegative reference weights permit no extra covering.
                separated_pairs+=separated(cells[i],cells[j])
            end
        end
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
    require(length(boundary)==2(m+1)*n+3m,"wrong exterior typed face count")
    ncell=sum(values(families));require(length(used)-length(edges)+length(faces)-ncell==1,"actual complex is not a ball")
    remaining=Set(keys(boundary))
    function subfaces(nodes)
        a,b,c,d=nodes
        [face for pattern in ((a,b,c),(a,b,d),(a,c,d),(b,c,d),(a,b,c,d))
            for face in (key(pattern),) if face in remaining]
    end
    top_diagonals=Set{Int}[]
    for cell in src.cells
        base=key(grid[1,v] for v in cell);require(base in remaining,"source Quad4 is missing");delete!(remaining,base)
        selected=subfaces(Tuple(grid[n+1,v] for v in cell))
        require(length(selected)==2 && all(length(face)==3 for face in selected),"copied cap is not two actual Tri3")
        diagonal=intersect(Set(selected[1]),Set(selected[2]));require(length(diagonal)==2,"top cap diagonal is missing")
        push!(top_diagonals,Set(col for col in cell if grid[n+1,col] in diagonal))
        if native_top
            start=findfirst(==(minimum(cell)),cell);opposite=cell[mod1(start+2,4)]
            require(diagonal==Set((grid[n+1,minimum(cell)],grid[n+1,opposite])),"top cap does not use original source-node identity")
        end
        setdiff!(remaining,selected)
    end
    for layer in 1:n,(a,b) in src.boundary
        selected=subfaces((grid[layer,a],grid[layer,b],grid[layer+1,a],grid[layer+1,b]))
        require(length(selected)==(f.laterals ? 1 : 2),"sampled boundary edge has wrong lateral split")
        setdiff!(remaining,selected)
    end
    require(isempty(remaining),"unrelated exterior typed face")
    require(global_rank_exists(src.cells,top_diagonals),"no single total source-vertex rank explains all saved cap diagonals")
    require(total==src.area*abs(Q(last(planes))-Q(first(planes))),"represented swept volume differs")
    return (;source=src,grid,heights=planes,centers,boundary,internal,total,ncell,families,domains,separated_pairs,top_diagonals,center_parent)
end

function certify_quadratic(mesh,f,source;oracle=false,native_top=!oracle)
    widths=Dict(11=>4,12=>8,13=>6,14=>5)
    require(all(haskey(widths,Int(block.msh)) for block in blocks(mesh)),"unexpected actual P2 family")
    primary=sort!(unique(Int(node) for block in blocks(mesh) for cell in eachcol(block.nodes)
        for node in cell[1:widths[Int(block.msh)]]))
    positions=Dict(node=>position for (position,node) in enumerate(primary))
    linear_blocks=[ElementBlock(Int(block.msh)-7,Int32[positions[Int(block.nodes[row,column])]
        for row in 1:widths[Int(block.msh)],column in axes(block.nodes,2)]) for block in blocks(mesh)]
    linear=MixedMesh(mesh.coords[:,primary],linear_blocks)
    underlying=certify(linear,f,source;oracle,native_top)
    support=Quad.certify_quadratic(mesh,linear)
    maps=[quadratic_map_certificate(Int(block.msh),Tuple(point(mesh,node) for node in cell))
        for block in blocks(mesh) for cell in eachcol(block.nodes)]
    require(length(maps)==underlying.ncell,"actual P2 cell set changed")
    require(nnodes(mesh)==(12f.strip_length+6)*f.intervals+14f.strip_length+3,
        "actual P2 support count differs")
    # Rounded actual P2 placement has its own exact dyadic map integral.
    # It is not replaced with P1 or forced to unit volume on nonunit shapes.
    total=sum(map.volume for map in maps)
    require(abs(total-underlying.total)<=Q(2e-11),"actual P2 integral differs from its certified P1 shell")
    return merge(underlying,(;total,p1_total=underlying.total,maps,nsupport=support.nsupport,linear))
end
function saved_fixture(record)
    raw=record["payload"];controls=record["derived_actual_curve_controls"]
    variant=raw["variant"];m=Int(controls["grid"][1]);axes=Tuple(Int.(variant["axes"]).+1)
    long=only(row for row in controls["statements"] if row["curves"]==controls["long_pair"])
    require(long["declared_count"]==m+1 && controls["long_pair"]==[1,3],"saved actual long controls differ")
    require(bytes2hex(sha256(raw["input_geo"]))==record["input_literal_sha256"],"saved literal input SHA differs")
    require(long["law"]==variant["law"] && long["coefficient"]==variant["coefficient"],"saved display metadata differs from literal controls")
    f=fixture(raw["name"];strip_length=m,layers=raw["intervals"]==1 ? :one : :nonbinary,
        height=Float64(variant["height"]))
    corners=Tuple(ntuple(k->k==axes[1] ? Float64(p[1]) : k==axes[2] ? Float64(p[2]) : 0.,3)
        for p in variant["corners"])
    # The helper call above supplies structural descriptor defaults only.
    # Exact captured input is ALWAYS substituted verbatim before execution;
    # never reconstruct it from law/display metadata or ideal sample fractions.
    return merge(f,(;source=raw["input_geo"],axes,corners,law=long["law"],
        coefficient=Float64(long["coefficient"]),actual_controls=controls,
        source_parity_provenance=record["future_source_parity_provenance"]))
end

end
