if !isdefined(@__MODULE__, :QuadTriNoNewQuadStripCertificates)
    include("quadtri_nonew_quad_strip_certificates.jl")
end

module QuadTriNoNewThreeQuadStripCertificates

using Tessella
using Tessella.Elements: MixedMesh,ElementBlock,msh_dimension,msh_spec
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
const prism_det=Strip.prism_det
const pyramid_det=Strip.pyramid_det
const map_positive=Strip.map_positive
const cell_volume=Strip.cell_volume
const separated=Strip.separated
const macro_boundary=Strip.macro_boundary
require(ok,message)=ok || throw(ArgumentError("Three-Quad strip certificate: "*message))

function fixture(name;plane=:XY,shape=:unit,winding=1,height=1.,laterals=false,
                 layers=:three,pins=:auto,curve_reverse=0,tags=:dense,
                 offset=(0.,0.,0.),scale=1.,arrangement="Left",progression=1.,long_pair=(1,3))
    f=Strip.fixture(name;plane,shape,winding,height,laterals,layers,pins,curve_reverse,
                    tags,offset,scale,arrangement,progression,long_pair)
    long_tags=[f.curve_tags[i] for i in long_pair]
    grading=progression==1. ? "" : " Using Progression $progression"
    from="Transfinite Curve{$(join(long_tags,','))}=3$grading;"
    to="Transfinite Curve{$(join(long_tags,','))}=4$grading;"
    require(occursin(from,f.source),"fixture source grading declaration missing")
    # The declared second group interpolates from represented .2 to 1.
    # Its midpoint is one ULP above the decimal literal .6.
    represented_levels=layers===:nonbinary ? [0.,.2,.2+(1.0-.2)/2,1.] : f.levels
    return merge(f,(;source=replace(f.source,from=>to),levels=represented_levels,intervals=length(f.levels)-1,
        translation=ntuple(k->k==f.axes[3] ? Float64(height) : 0.,3),
        base=(curves=f.curve_tags,surface=f.surface)))
end

function source_complex(source,f)
    require(nnodes(source)==8 && all(isfinite,source.coords),"source must have eight finite nodes")
    require(length(Set(point(source,i) for i in 1:8))==8,"source coordinates coincide")
    cells=Tuple(Tuple(Int.(cell)) for b in blocks(source) for cell in eachcol(b.nodes) if b.msh==3)
    require(length(cells)==3 && sum(size(b.nodes,2) for b in blocks(source))==3,
            "source must contain exactly three Quad4")
    require(all(source.coords[f.axes[3],:].==source.coords[f.axes[3],1]),"source is not exactly axis planar")
    xy=Tuple(Tuple(Q(source.coords[a,i]) for a in f.axes[1:2]) for i in 1:8)
    incidence=Dict{Tuple,Vector{Tuple}}();parents=Dict{Tuple,Vector{Int}}();areas=Q[];used=Set{Int}()
    for (parent,cell) in enumerate(cells)
        require(allunique(cell) && all(i->1<=i<=8,cell),"foreign or repeated source node")
        union!(used,cell)
        require(all(sign(orient(xy[cell[i]],xy[cell[mod1(i+1,4)]],xy[cell[mod1(i+2,4)]]))==f.winding for i in 1:4),
                "actual Quad4 is not strictly convex with CAD winding")
        push!(areas,abs(sum(xy[cell[i]][1]*xy[cell[mod1(i+1,4)]][2]-xy[cell[mod1(i+1,4)]][1]*xy[cell[i]][2] for i in 1:4))/2)
        for i in 1:4
            edge=(cell[i],cell[mod1(i+1,4)]);support=key(edge)
            push!(get!(incidence,support,Tuple[]),edge)
            push!(get!(parents,support,Int[]),parent)
        end
    end
    require(used==Set(1:8) && length(incidence)==10,"source is not V8/E10/F3")
    boundary=Set{Tuple}();shared=Tuple[];adjacency=[Int[] for _ in 1:3]
    for (edge,uses) in incidence
        require(length(uses) in (1,2),"source edge has excess incidence")
        if length(uses)==1;push!(boundary,only(uses))
        else
            require(uses[1]==reverse(uses[2]),"shared source edge has equal cycles")
            push!(shared,edge);left,right=parents[edge]
            push!(adjacency[left],right);push!(adjacency[right],left)
        end
    end
    require(length(boundary)==8 && length(shared)==2,"wrong source shared/exterior edge count")
    require(sort(length.(adjacency))==[1,1,2] && all(allunique,adjacency),"source dual is not a three-cell path")
    border=Set(v for edge in boundary for v in edge)
    require(border==Set(1:8),"source has a nonboundary node")
    next=Dict(a=>b for (a,b) in boundary)
    require(length(next)==8 && length(Set(values(next)))==8,"source boundary branches")
    start=minimum(keys(next));walk=Int[start]
    for _ in 1:7;push!(walk,next[last(walk)]);end
    require(allunique(walk) && next[last(walk)]==start,"source boundary is not one eight-node cycle")
    for i in 1:8,j in i+1:8
        a,b=walk[i],walk[mod1(i+1,8)];c,d=walk[j],walk[mod1(j+1,8)]
        isempty(intersect((a,b),(c,d))) || continue
        require(!crosses(xy[a],xy[b],xy[c],xy[d]),"actual sampled boundary crosses")
    end
    for edge in shared
        ids=parents[edge]
        sides=Tuple(Tuple(sign(orient(xy[edge[1]],xy[edge[2]],xy[v])) for v in cells[i] if !(v in edge)) for i in ids)
        require(all(s->all(==(s[1]),s) && s[1]!=0,sides) && sides[1][1]==-sides[2][1],
                "adjacent actual Quads do not have disjoint interiors")
    end
    # A path graph alone does not exclude overlapping nonneighbor cells.
    for i in 1:3,j in i+1:3
        require(any(all(f.winding*orient(xy[a],xy[b],xy[v])<=0 for v in other)
            for (own,other) in ((cells[i],cells[j]),(cells[j],cells[i]))
            for (a,b) in ((own[k],own[mod1(k+1,4)]) for k in 1:4)),
            "source Quad interiors overlap")
    end
    polygon=abs(sum(xy[a][1]*xy[b][2]-xy[b][1]*xy[a][2] for (a,b) in boundary))/2
    require(sum(areas)==polygon,"source cells do not cover the actual boundary")
    cad=Tuple(only(i for i in 1:8 if point(source,i)==position(p)) for p in f.corners)
    neighbors=Dict(i=>Int[] for i in 1:8)
    for (a,b) in boundary;push!(neighbors[a],b);push!(neighbors[b],a);end
    chains=Tuple(begin
        first,last=cad[i],cad[mod1(i+1,4)];paths=Tuple[]
        for neighbor in neighbors[first]
            path=Int[first,neighbor]
            while path[end]!=last && !(path[end] in cad)
                push!(path,only(v for v in neighbors[path[end]] if v!=path[end-1]))
                length(path)<=8 || error("source boundary walk did not terminate")
            end
            path[end]==last && push!(paths,Tuple(path))
        end
        require(length(paths)==1,"CAD corner pair lacks one actual boundary chain")
        only(paths)
    end for i in 1:4)
    require(sort(length.(collect(chains)))==[2,2,4,4],"native chain widths differ")
    require(all(length(chains[i])==(i in f.long_pair ? 4 : 2) for i in 1:4),"long-chain CAD identity differs")
    return (;cells,boundary,shared=Tuple(sort(shared)),boundary_cycle=Tuple(walk),chains,
            xy,areas,area=polygon,adjacency=Tuple(Tuple(row) for row in adjacency),incidence)
end

function certify(mesh,f,source;oracle=false,native_top=!oracle,check_separation=true)
    src=source_complex(source,f);n=length(f.levels)-1
    expected_nodes=8(n+1)+(f.laterals ? 3 : 0)
    require(nnodes(mesh)==expected_nodes && all(isfinite,mesh.coords),"wrong primary/centroid node count")
    intended=[source.coords[f.axes[3],1]+Float64(t)*f.height for t in f.levels]
    planes=if oracle
        columns=Set(Tuple(source.coords[a,v] for a in f.axes[1:2]) for v in 1:8)
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
        for row in 1:n+1,col in 1:8]
    primary=Set(vec(grid));centers=sort!(collect(setdiff(Set(1:nnodes(mesh)),primary)))
    require(length(primary)==8(n+1) && length(centers)==(f.laterals ? 3 : 0),"column bijection differs")
    center_parent=Dict{Int,Int}()
    for center in centers
        p=exact(point(mesh,center));xy=Tuple(p[a] for a in f.axes[1:2])
        require(min(Q(planes[end-1]),Q(planes[end]))<p[f.axes[3]]<max(Q(planes[end-1]),Q(planes[end])),
            "actual terminal centroid is not strictly between its planes")
        parents=[i for (i,cell) in enumerate(src.cells) if all(sign(orient(src.xy[cell[k]],src.xy[cell[mod1(k+1,4)]],xy))==f.winding for k in 1:4)]
        require(length(parents)==1,"actual terminal centroid lacks one strict macro owner")
        center_parent[center]=only(parents)
    end
    require(Set(values(center_parent))==Set(f.laterals ? (1,2,3) : ()) ,"terminal centroids do not occupy distinct macros")
    location=Dict(grid[row,col]=>(row,col) for row in 1:n+1,col in 1:8)
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
        parents=[i for i in 1:3 if Set(column_ids)⊆Set(src.cells[i]) && all(center_parent[v]==i for v in body)]
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
        expected=Dict(4=>6,7=>15);n>1 && (expected[5]=3(n-1))
        require(families==expected,"wrong recombined terminal-fan families")
    end
    separated_pairs=0
    for parent in 1:3,layer in 1:n
        cells=domains[(parent,layer)]
        require(sum(c.volume for c in cells)==src.areas[parent]*abs(Q(planes[layer+1])-Q(planes[layer])),"actual reference maps do not exactly partition the macro")
        vertices=Tuple(grid[row,col] for row in (layer,layer+1) for col in src.cells[parent])
        macro_boundary(cells,vertices,mesh)
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
    require(length(boundary)==(f.laterals ? 8n+9 : 16n+9),"wrong exterior typed face count")
    ncell=sum(values(families));require(length(used)-length(edges)+length(faces)-ncell==1,"actual complex is not a ball")
    remaining=Set(keys(boundary))
    function subfaces(nodes)
        a,b,c,d=nodes
        [face for pattern in ((a,b,c),(a,b,d),(a,c,d),(b,c,d),(a,b,c,d))
            for face in (key(pattern),) if face in remaining]
    end
    for cell in src.cells
        base=key(grid[1,v] for v in cell);require(base in remaining,"source Quad4 is missing");delete!(remaining,base)
        selected=subfaces(Tuple(grid[n+1,v] for v in cell))
        require(length(selected)==2 && all(length(face)==3 for face in selected),"copied cap is not two actual Tri3")
        diagonal=intersect(Set(selected[1]),Set(selected[2]));require(length(diagonal)==2,"top cap diagonal is missing")
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
    require(total==src.area*abs(Q(last(planes))-Q(first(planes))),"represented swept volume differs")
    return (;source=src,grid,heights=planes,centers,boundary,internal,total,ncell,families,domains,separated_pairs)
end


end
