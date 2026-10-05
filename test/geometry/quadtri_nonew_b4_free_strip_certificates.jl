# Actual whole maps and typed macro shells; no fixed free family signature.
if !isdefined(@__MODULE__, :QuadTriNoNewB4RecStripCertificates)
    include("quadtri_nonew_b4_rec_strip_certificates.jl")
end
module QuadTriNoNewB4FreeStripCertificates
using Tessella, SHA
using Tessella.Elements: MixedMesh, ElementBlock, msh_dimension
using Tessella.MeshTypes: nnodes
using ..QuadTriNoNewB4RecStripCertificates
const Rec=QuadTriNoNewB4RecStripCertificates
const Rect=Rec.Rect
const Strip=Rec.Strip
const Quad=Rec.Quad
const Q=Rec.Q
const FACES=Rec.FACES
const point=Rec.point
const position=Rec.position
const exact=Rec.exact
const key=Rec.key
const cycle=Rec.cycle
const det=Rec.det
const orient=Rec.orient
const blocks=Rec.blocks
const crosses=Rec.crosses
const map_positive=Rec.map_positive
const cell_volume=Rec.cell_volume
const separated=Rec.separated
const macro_boundary=Rec.macro_boundary
const execute=Rec.execute
const public_volume=Rec.public_volume
const packed_digest=Rec.packed_digest
const quadratic_map_certificate=Rec.quadratic_map_certificate
const source_complex=Rec.source_complex
require(ok,message)=ok || throw(ArgumentError("B4 free-strip certificate: "*message))

function fixture(name;laterals=false,kwargs...)
    require(!laterals,"requires free lateral faces")
    return Rec.fixture(name;laterals,kwargs...)
end

function saved_fixture(record)
    raw=record["payload"];controls=record["derived_actual_curve_controls"]
    variant=raw["variant"];grid=Tuple(Int.(controls["grid"]));m=prod(grid)
    axes=Tuple(Int.(variant["axes"]).+1)
    long=only(row for row in controls["statements"] if row["curves"]==controls["long_pair"])
    require(long["declared_count"]==m+1 && controls["long_pair"]==[1,3],"saved actual long controls differ")
    require(bytes2hex(sha256(raw["input_geo"]))==record["input_literal_sha256"],"saved literal input SHA differs")
    require(long["law"]==variant["law"] && long["coefficient"]==variant["coefficient"],"literal law differs from display metadata")
    f=fixture(raw["name"];strip_length=m,height=Float64(variant["height"]),layers=:one)
    corners=Tuple(ntuple(k->k==axes[1] ? Float64(p[1]) : k==axes[2] ? Float64(p[2]) : 0.,3) for p in variant["corners"])
    levels=Float64[0.];previous=0.
    for (count,height) in zip(variant["groups"],variant["heights"])
        count=Int(count);height=Float64(height)
        for step in 1:count;push!(levels,previous+(height-previous)*step/count);end
        previous=height
    end
    require(length(levels)==Int(raw["intervals"])+1,"saved layer groups differ")
    # Use the original literal, including its CRLF bytes and actual sampling
    # controls. The fixture builder supplies descriptor defaults only.
    return merge(f,(;source=raw["input_geo"],axes,corners,grid,levels,
        intervals=length(levels)-1,law=long["law"],coefficient=Float64(long["coefficient"]),
        actual_controls=controls))
end

function certify(mesh,f,source;oracle=false,native_top=!oracle,check_separation=true)
    src=source_complex(source,f);n=length(f.levels)-1;m=prod(f.grid);v=2(m+1)
    require(!f.laterals && all(isfinite,mesh.coords),"requires a finite free product")
    require(v*(n+1)<=nnodes(mesh)<=v*(n+1)+m*n,"wrong possible primary/center node count")
    intended=[source.coords[f.axes[3],1]+Float64(t)*f.height for t in f.levels]
    columns=Set(Tuple(source.coords[a,node] for a in f.axes[1:2]) for node in 1:v)
    planes=oracle ? sort!(unique(mesh.coords[f.axes[3],node] for node in 1:nnodes(mesh)
        if Tuple(mesh.coords[a,node] for a in f.axes[1:2]) in columns);rev=f.height<0) : intended
    require(length(planes)==n+1 && all(sign.(diff(planes)).==sign(f.height)),"actual planes collapse or reverse")
    require(all(abs(planes[i]-intended[i])<=2e-11 for i in eachindex(planes)),"actual heights differ from declared levels")
    lookup=Dict(point(mesh,node)=>node for node in 1:nnodes(mesh))
    require(length(lookup)==nnodes(mesh),"actual output has coincident coordinates")
    require(all(haskey(lookup,ntuple(k->k==f.axes[3] ? planes[row] : source.coords[k,col],3))
        for row in 1:n+1,col in 1:v),"a required original source column is missing")
    grid=[lookup[ntuple(k->k==f.axes[3] ? planes[row] : source.coords[k,col],3)] for row in 1:n+1,col in 1:v]
    primary=Set(vec(grid));centers=sort!(collect(setdiff(Set(1:nnodes(mesh)),primary)))
    require(length(primary)==v*(n+1),"column bijection differs")
    center_parent=Dict{Int,Tuple{Int,Int}}()
    for center in centers
        p=exact(point(mesh,center));xy=Tuple(p[a] for a in f.axes[1:2])
        parents=[i for (i,cell) in enumerate(src.cells) if all(sign(orient(src.xy[cell[k]],src.xy[cell[mod1(k+1,4)]],xy))==f.winding for k in 1:4)]
        intervals=[layer for layer in 1:n if min(Q(planes[layer]),Q(planes[layer+1]))<p[f.axes[3]]<max(Q(planes[layer]),Q(planes[layer+1]))]
        require(length(parents)==length(intervals)==1,"center lacks one strict macro/interval owner")
        parent=only(parents);layer=only(intervals)
        expected=ntuple(k->k==f.axes[3] ? (Q(planes[layer])+Q(planes[layer+1]))/2 : sum(Q(source.coords[k,node]) for node in src.cells[parent])/4,3)
        require(all(abs(p[k]-expected[k])<=Q(2e-11) for k in 1:3),"center differs from the exact mean of its eight corners")
        center_parent[center]=(parent,layer)
    end
    require(length(Set(values(center_parent)))==length(centers),"multiple centers occupy one macro")
    location=Dict(grid[row,col]=>(row,col) for row in 1:n+1,col in 1:v)
    domains=Dict{Tuple,Vector{NamedTuple}}();faces=Dict{Tuple,Vector{Tuple}}();edges=Set{Tuple}()
    used=Set{Int}();families=Dict{Int,Int}();total=zero(Q)
    for block in blocks(mesh),cell in eachcol(block.nodes)
        require(msh_dimension(block.msh)==3,"certificate expects actual volume cells only")
        msh=Int(block.msh);require(haskey(FACES,msh) && allunique(cell),"unexpected family or repeated node IDs")
        require(all(node->1<=node<=nnodes(mesh),cell),"cell references a foreign node")
        body=[Int(node) for node in cell if !(node in primary)]
        column_ids=[location[Int(node)][2] for node in cell if node in primary]
        rows=[location[Int(node)][1] for node in cell if node in primary]
        require(!isempty(rows) && length(body)<=1,"cell has multiple body nodes or no source columns")
        layer=isempty(body) ? minimum(rows) : center_parent[only(body)][2]
        require(all(row->row in (layer,layer+1),rows) && (!isempty(body) ||
            minimum(rows)==layer && maximum(rows)==layer+1),"cell crosses nonadjacent planes")
        parents=[i for i in 1:m if Set(column_ids)⊆Set(src.cells[i]) && all(center_parent[node]==(i,layer) for node in body)]
        require(length(parents)==1,"cell has ambiguous or foreign macro membership")
        parent=only(parents);p=Tuple(exact(point(mesh,node)) for node in cell)
        require(map_positive(p,msh),"nonpositive whole P1 reference map")
        volume=cell_volume(p,msh);require(volume>0,"nonpositive exact actual-map volume")
        total+=volume;union!(used,Int.(cell));families[msh]=get(families,msh,0)+1
        push!(get!(domains,(parent,layer),NamedTuple[]),(;p,msh,volume,nodes=Tuple(Int.(cell))))
        for pattern in FACES[msh]
            face=Tuple(Int(cell[k]) for k in pattern);push!(get!(faces,key(face),Tuple[]),face)
            for k in eachindex(face);push!(edges,key((face[k],face[mod1(k+1,length(face))])));end
        end
    end
    require(used==Set(1:nnodes(mesh)),"unused primary or center node")
    require(length(domains)==m*n,"empty or unrelated macro interval")
    states=Dict{Tuple,NTuple{6,Int}}();macro_groups=Dict{Tuple,Vector{Vector{Tuple}}}();separated_pairs=0
    logical=((1,2,6,5),(2,3,7,6),(3,4,8,7),(4,1,5,8),(1,2,3,4),(5,6,7,8))
    for parent in 1:m,layer in 1:n
        cells=domains[(parent,layer)]
        require(sum(cell.volume for cell in cells)==src.areas[parent]*abs(Q(planes[layer+1])-Q(planes[layer])),"actual maps do not exactly partition their macro")
        vertices=Tuple(grid[row,col] for row in (layer,layer+1) for col in src.cells[parent])
        groups,_=macro_boundary(cells,vertices,mesh)
        ordered=[only(group for group in groups if !isempty(group) && Set(first(group))⊆Set(vertices[k] for k in pattern)) for pattern in logical]
        states[(parent,layer)]=ntuple(k->Rect.state_of(ordered[k],Tuple(vertices[j] for j in logical[k])),6)
        macro_groups[(parent,layer)]=ordered
        owners=[center for (center,owner) in center_parent if owner==(parent,layer)]
        if !isempty(owners)
            center=only(owners)
            require(all(cell.msh in (4,7) && center in cell.nodes for cell in cells),"retained macro is not its actual center fan")
            bases=[key(node for node in cell.nodes if node!=center) for cell in cells]
            shell=[key(face) for group in ordered for face in group]
            require(length(bases)==length(shell) && Set(bases)==Set(shell),"retained fan cells do not bijectively cover its final face patches")
            require(all(key((center,node)) in edges for node in vertices),"retained fan lacks an actual radial edge")
        end
        if check_separation
            # Planar supporting faces give a second proof where available.
            # Otherwise positive whole maps, opposite identical typed traces,
            # the exact macro shell and degree-one integral certify covering.
            for i in 1:length(cells),j in i+1:length(cells);separated_pairs+=separated(cells[i],cells[j]);end
        end
    end
    boundary=Dict{Tuple,Tuple}();internal=Dict{Tuple,Vector{Tuple}}()
    for (face,uses) in faces
        require(length(uses) in (1,2),"typed face incidence exceeds two")
        if length(uses)==2
            require(cycle(uses[1])==cycle(reverse(uses[2])),"shared outward typed cycles agree")
            internal[face]=uses
        else;boundary[face]=only(uses)
        end
    end
    require(length(boundary)==4(m+1)*n+3m,"wrong exterior typed face count")
    ncell=sum(values(families));require(length(used)-length(edges)+length(faces)-ncell==1,"actual complex is not a ball")
    remaining=Set(keys(boundary));top_diagonals=Set{Int}[]
    function subfaces(nodes)
        a,b,c,d=nodes
        return [face for pattern in ((a,b,c),(a,b,d),(a,c,d),(b,c,d),(a,b,c,d)) for face in (key(pattern),) if face in remaining]
    end
    for (parent,cell) in enumerate(src.cells)
        base=key(grid[1,node] for node in cell);require(base in remaining,"source Quad4 is missing");delete!(remaining,base)
        selected=subfaces(Tuple(grid[n+1,node] for node in cell))
        require(length(selected)==2 && all(length(face)==3 for face in selected),"copied cap is not two actual triangles")
        diagonal=intersect(Set(selected[1]),Set(selected[2]));require(length(diagonal)==2,"top diagonal is missing")
        push!(top_diagonals,Set(node for node in cell if grid[n+1,node] in diagonal))
        if native_top
            start=findfirst(==(minimum(cell)),cell);opposite=cell[mod1(start+2,4)]
            require(diagonal==Set((grid[n+1,minimum(cell)],grid[n+1,opposite])),"top cap differs from original source-node identity")
        end
        setdiff!(remaining,selected)
        for layer in 1:n-1;require(states[(parent,layer)][6]==states[(parent,layer+1)][5],"intermediate cap carry differs");end
    end
    for layer in 1:n,(a,b) in src.boundary
        selected=subfaces((grid[layer,a],grid[layer,b],grid[layer+1,a],grid[layer+1,b]))
        require(length(selected)==2 && all(length(face)==3 for face in selected),"exterior lateral is not two actual triangles")
        setdiff!(remaining,selected)
    end
    require(isempty(remaining),"unrelated exterior typed face")
    require(Rec.global_rank_exists(src.cells,top_diagonals),"no total source rank explains all saved cap diagonals")
    require(total==src.area*abs(Q(last(planes))-Q(first(planes))),"represented swept volume differs")
    return (;source=src,grid,heights=planes,centers,center_parent,boundary,internal,total,ncell,families,domains,states,macro_groups,separated_pairs,top_diagonals)
end

function quadratic_support_certificate(mesh,linear,primary)
    widths=Dict(11=>4,12=>8,13=>6,14=>5)
    support_nodes=Dict{Tuple,Int}();node_supports=Dict{Int,Tuple}()
    used=Set{Int}();faces=Dict{Tuple,Vector{Tuple}}()
    for block in blocks(mesh),cell in eachcol(block.nodes)
        msh=Int(block.msh);width=widths[msh]
        reference=Rect._P2_REFERENCE_TWICE[msh]
        require(length(cell)==length(reference) && allunique(cell),"actual quadratic width or node identity differs")
        source=Tuple(Int(node) for node in cell[1:width]);local_supports=Tuple[]
        for (index,ref) in enumerate(reference)
            u,v,w=Tuple(Q(value)/2 for value in ref)
            weights=Tuple(Q(value) for value in Quad.primary_weights(msh-7,u,v,w))
            support=Tuple(sort!([(source[i],weights[i]) for i in 1:width if weights[i]!=0]))
            expected=ntuple(k->sum(weights[i]*Q(mesh.coords[k,source[i]]) for i in 1:width),3)
            require(all(abs(Q(mesh.coords[k,cell[index]])-expected[k])<=Q(2e-11) for k in 1:3),"actual quadratic support differs from its independent primary interpolation")
            require(get(support_nodes,support,Int(cell[index]))==cell[index],"shared actual quadratic support is split")
            require(get(node_supports,Int(cell[index]),support)==support,"unrelated actual quadratic supports are welded")
            support_nodes[support]=Int(cell[index]);node_supports[Int(cell[index])]=support
            push!(local_supports,Tuple(first(pair) for pair in support));push!(used,Int(cell[index]))
        end
        for pattern in FACES[msh-7]
            carrier=key(source[k] for k in pattern)
            trace=key(Int(cell[k]) for k in eachindex(cell) if Set(local_supports[k])⊆Set(carrier))
            require(length(trace)==(length(carrier)==3 ? 6 : 9),"actual quadratic face trace has wrong width")
            push!(get!(faces,carrier,Tuple[]),trace)
        end
    end
    require(used==Set(1:nnodes(mesh)) && length(support_nodes)==nnodes(mesh),"unused or duplicate actual quadratic supports")
    require(all(length(rows) in (1,2) && (length(rows)==1 || first(rows)==last(rows)) for rows in values(faces)),"opposite actual quadratic face traces differ")
    require(Set(primary)==Set(node for (node,support) in node_supports if length(support)==1 && last(only(support))==1),"actual quadratic primary identity differs")
    return (;nsupport=length(support_nodes),node_supports,support_nodes,faces)
end

function certify_quadratic(mesh,f,source;oracle=false,native_top=!oracle)
    widths=Dict(11=>4,12=>8,13=>6,14=>5)
    require(all(haskey(widths,Int(block.msh)) for block in blocks(mesh)),"unexpected actual P2 family")
    primary=sort!(unique(Int(node) for block in blocks(mesh) for cell in eachcol(block.nodes) for node in cell[1:widths[Int(block.msh)]]))
    positions=Dict(node=>position for (position,node) in enumerate(primary))
    linear_blocks=[ElementBlock(Int(block.msh)-7,Int32[positions[Int(block.nodes[row,column])]
        for row in 1:widths[Int(block.msh)],column in axes(block.nodes,2)]) for block in blocks(mesh)]
    linear=MixedMesh(mesh.coords[:,primary],linear_blocks)
    underlying=certify(linear,f,source;oracle,native_top)
    support=quadratic_support_certificate(mesh,linear,primary)
    maps=[quadratic_map_certificate(Int(block.msh),Tuple(point(mesh,node) for node in cell)) for block in blocks(mesh) for cell in eachcol(block.nodes)]
    require(length(maps)==underlying.ncell,"actual P2 cell set changed")
    total=sum(map.volume for map in maps)
    require(abs(total-underlying.total)<=Q(2e-11),"actual P2 integral differs from its certified P1 shell")
    return merge(underlying,(;total,p1_total=underlying.total,maps,nsupport=support.nsupport,linear,support))
end
include("quadtri_nonew_b4_free_strip_retained_products.jl")
end
