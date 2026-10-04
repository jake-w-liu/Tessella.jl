if !isdefined(@__MODULE__, :QuadTriNoNewQuadPatchCertificates)
    include("quadtri_nonew_quad_patch_certificates.jl")
end

module QuadTriNoNewQuadStripCertificates

using Tessella
using Tessella.Elements: MixedMesh,ElementBlock,msh_dimension,msh_spec
using Tessella.MeshTypes: nnodes
using ..QuadTriNoNewQuadPatchCertificates
const Patch=QuadTriNoNewQuadPatchCertificates
const Quad=Patch.Quad
const Q=Patch.Q
const FACES=merge(Patch.FACES,Dict(6=>((1,3,2),(4,5,6),(1,2,5,4),(2,3,6,5),(3,1,4,6))))
const point=Patch.point
const position=Patch.position
const exact=Patch.exact
const key=Patch.key
const cycle=Patch.cycle
const det=Patch.det
const orient=Patch.orient
const blocks=Patch.blocks
const execute=Patch.execute
const public_volume=Patch.public_volume
const packed_digest=Patch.packed_digest
require(ok,message)=ok || throw(ArgumentError("Quad-strip certificate: "*message))

function crosses(a,b,c,d)
    on(p,q,r)=orient(p,q,r)==0 &&
        all(min(p[k],q[k])<=r[k]<=max(p[k],q[k]) for k in 1:2)
    ac,ad,ca,cb=orient(a,b,c),orient(a,b,d),orient(c,d,a),orient(c,d,b)
    return sign(ac)*sign(ad)<0 && sign(ca)*sign(cb)<0 ||
        on(a,b,c) || on(a,b,d) || on(c,d,a) || on(c,d,b)
end

function fixture(name;plane=:XY,shape=:unit,winding=1,height=1.,laterals=false,
                 layers=:three,pins=:auto,curve_reverse=0,tags=:dense,
                 offset=(0.,0.,0.),scale=1.,arrangement="Left",progression=1.,long_pair=(1,3))
    base_plane=plane===:XZ ? :ZX : plane
    base_offset=plane===:XZ ? (offset[3],offset[2],offset[1]) : offset
    f=Patch.fixture(name;plane=base_plane,shape,winding,height,laterals,layers,pins,
                    curve_reverse,tags,offset=base_offset,scale,arrangement)
    if plane===:XZ
        corners=Tuple((p[3],p[2],p[1]) for p in f.corners)
        source=f.source
        for (tag,old,new) in zip(f.point_tags,f.corners,corners)
            source=replace(source,"Point($tag)={$(join(old,',')),1};"=>"Point($tag)={$(join(new,',')),1};")
        end
        f=merge(f,(;source,corners,axes=(1,3,2),plane,offset))
    end
    long_tags=[f.curve_tags[i] for i in long_pair]
    short_tags=[c for c in f.curve_tags if !(c in long_tags)]
    grading=progression==1. ? "" : " Using Progression $progression"
    text=replace(f.source,"Transfinite Curve{$(join(f.curve_tags,','))}=3;"=>
        "Transfinite Curve{$(join(long_tags,','))}=3$grading;Transfinite Curve{$(join(short_tags,','))}=2;")
    return merge(f,(;source=text,progression=Float64(progression),long_pair))
end

function source_complex(source,f)
    require(nnodes(source)==6 && all(isfinite,source.coords),"source must have six finite nodes")
    require(length(Set(point(source,i) for i in 1:6))==6,"source coordinates coincide")
    cells=Tuple(Tuple(Int.(cell)) for b in blocks(source) for cell in eachcol(b.nodes)
        if b.msh==3)
    require(length(cells)==2 && sum(size(b.nodes,2) for b in blocks(source))==2,
        "source must contain exactly two Quad4")
    require(all(source.coords[f.axes[3],:].==source.coords[f.axes[3],1]),"source is not exactly axis planar")
    xy=Tuple(Tuple(Q(source.coords[a,i]) for a in f.axes[1:2]) for i in 1:6)
    incidence=Dict{Tuple,Vector{Tuple}}();areas=Q[];used=Set{Int}()
    for cell in cells
        require(allunique(cell) && all(i->1<=i<=6,cell),"foreign or repeated source node")
        union!(used,cell)
        require(all(sign(orient(xy[cell[i]],xy[cell[mod1(i+1,4)]],xy[cell[mod1(i+2,4)]]))==f.winding for i in 1:4),
            "actual Quad4 is not strictly convex with CAD winding")
        push!(areas,abs(sum(xy[cell[i]][1]*xy[cell[mod1(i+1,4)]][2]-xy[cell[mod1(i+1,4)]][1]*xy[cell[i]][2] for i in 1:4))/2)
        for i in 1:4
            edge=(cell[i],cell[mod1(i+1,4)])
            push!(get!(incidence,key(edge),Tuple[]),edge)
        end
    end
    require(used==Set(1:6) && length(incidence)==7,"source is not V6/E7/F2")
    boundary=Set{Tuple}();shared=Tuple[]
    for (edge,uses) in incidence
        require(length(uses) in (1,2),"source edge has excess incidence")
        if length(uses)==1;push!(boundary,only(uses))
        else
            require(uses[1]==reverse(uses[2]),"shared source edge has equal cycles")
            push!(shared,edge)
        end
    end
    require(length(boundary)==6 && length(shared)==1,"wrong source shared/exterior edge count")
    next=Dict(a=>b for (a,b) in boundary)
    require(length(next)==6,"source boundary branches")
    start=minimum(keys(next));walk=Int[start]
    for _ in 1:5;push!(walk,next[last(walk)]);end
    require(allunique(walk) && next[last(walk)]==start,"source boundary is not one six-node cycle")
    for i in 1:6,j in i+1:6
        a,b=walk[i],walk[mod1(i+1,6)];c,d=walk[j],walk[mod1(j+1,6)]
        isempty(intersect((a,b),(c,d))) || continue
        require(!crosses(xy[a],xy[b],xy[c],xy[d]),"actual sampled boundary crosses")
    end
    common=only(shared)
    sides=Tuple(Tuple(sign(orient(xy[common[1]],xy[common[2]],xy[v])) for v in cell if !(v in common)) for cell in cells)
    require(all(s->all(==(s[1]),s) && s[1]!=0,sides) && sides[1][1]==-sides[2][1],
        "actual Quads do not have disjoint interiors across their shared edge")
    polygon=abs(sum(xy[a][1]*xy[b][2]-xy[b][1]*xy[a][2] for (a,b) in boundary))/2
    require(sum(areas)==polygon,"source cells do not cover the actual boundary")
    cad=Tuple(only(i for i in 1:6 if point(source,i)==position(p)) for p in f.corners)
    border_edges=Set(key(edge) for edge in boundary)
    chains=Tuple(begin
        a,b=cad[i],cad[mod1(i+1,4)]
        key((a,b)) in border_edges ? (a,b) :
            (a,only(v for v in 1:6 if v!=a && v!=b && key((a,v)) in border_edges && key((v,b)) in border_edges),b)
    end for i in 1:4)
    require(sort!(length.(collect(chains)))==[2,2,3,3],"native chain widths differ")
    return (;cells,boundary,shared=common,boundary_cycle=Tuple(walk),chains,xy,areas,area=polygon)
end

# Pri6's determinant is affine on each triangular slice and quadratic in w.
# Its exact extrema therefore reduce to three independent quadratic minima.
function prism_det(p,u,v,w)
    du=ntuple(d->(1-w)*(p[2][d]-p[1][d])+w*(p[5][d]-p[4][d]),3)
    dv=ntuple(d->(1-w)*(p[3][d]-p[1][d])+w*(p[6][d]-p[4][d]),3)
    dw=ntuple(d->(1-u-v)*(p[4][d]-p[1][d])+u*(p[5][d]-p[2][d])+v*(p[6][d]-p[3][d]),3)
    return det(du,dv,dw)
end
function pyramid_det(p,u,v)
    bits=((0,0),(1,0),(1,1),(0,1))
    base=ntuple(d->sum((a==0 ? 1-u : u)*(b==0 ? 1-v : v)*p[i][d]
        for (i,(a,b)) in enumerate(bits)),3)
    du=ntuple(d->sum((a==0 ? -1 : 1)*(b==0 ? 1-v : v)*p[i][d]
        for (i,(a,b)) in enumerate(bits)),3)
    dv=ntuple(d->sum((a==0 ? 1-u : u)*(b==0 ? -1 : 1)*p[i][d]
        for (i,(a,b)) in enumerate(bits)),3)
    return det(du,dv,p[5].-base)
end
function map_positive(p,msh)
    msh==7 && return all(pyramid_det(p,u,v)>0 for (u,v) in ((0,0),(1,0),(1,1),(0,1)))
    msh!=6 && return Patch.map_positive(p,msh)
    for (u,v) in ((0,0),(1,0),(0,1))
        lo=prism_det(p,u,v,0);mid=prism_det(p,u,v,1//2);hi=prism_det(p,u,v,1)
        min(lo,hi)>0 || return false
        a=2(lo+hi-2mid);b=hi-lo-a
        if a>0
            w=-b/(2a)
            0<w<1 && lo+b*w+a*w*w<=0 && return false
        end
    end
    return true
end
function cell_volume(p,msh)
    # The collapsed Pyramid5 map contributes (1-w)^2 times a bilinear
    # determinant. Integrate it exactly; a warped Quad base is not replaced
    # by two flat triangles.
    msh==7 && return sum(pyramid_det(p,u,v) for (u,v) in ((0,0),(1,0),(1,1),(0,1)))/12
    msh!=6 && return Patch.cell_volume(p,msh)
    # Exact triangle integration averages its three vertex determinants;
    # Simpson integrates their quadratic w dependence without approximation.
    return sum(prism_det(p,u,v,0)+4prism_det(p,u,v,1//2)+prism_det(p,u,v,1)
               for (u,v) in ((0,0),(1,0),(0,1)))/36
end
function separated(left,right)
    # Only actual planar supporting faces count. A warped prism's Quad side
    # is never silently replaced by a planar face or fixed tetra partition.
    any(all(det(a,b,q.-p[face[1]])>=0 for q in other) &&
        all(det(a,b,q.-p[face[1]])<=0 for q in p)
        for (p,msh,other) in ((left.p,left.msh,right.p),(right.p,right.msh,left.p))
        for face in FACES[msh] for (a,b) in ((p[face[2]].-p[face[1]],p[face[3]].-p[face[1]]),)
        if all(det(a,b,p[k].-p[face[1]])==0 for k in face))
end

function macro_boundary(cells,vertices,coordinates)
    incidence=Dict{Tuple,Vector{Tuple}}()
    for cell in cells,pattern in FACES[cell.msh]
        face=Tuple(cell.nodes[i] for i in pattern)
        push!(get!(incidence,key(face),Tuple[]),face)
    end
    reference=[Tuple(vertices[i] for i in pattern) for pattern in FACES[5]]
    p=Tuple(exact(point(coordinates,v)) for v in vertices)
    det(p[2].-p[1],p[4].-p[1],p[5].-p[1])<0 && (reference=reverse.(reference))
    groups=[Tuple[] for _ in 1:6]
    for (support,faces) in incidence
        require(length(faces) in (1,2),"macro face has excess incidence")
        if length(faces)==2
            require(cycle(faces[1])==cycle(reverse(faces[2])),"macro internal maps have equal cycles")
        else
            carriers=findall(face->Set(support)⊆Set(face),reference)
            require(length(carriers)==1,"macro exposes an interior or foreign face")
            push!(groups[only(carriers)],only(faces))
        end
    end
    for (carrier,faces) in zip(reference,groups)
        require(length(faces) in (1,2),"macro carrier has a missing or excess face patch")
        q=Tuple(exact(point(coordinates,v)) for v in carrier)
        axes=first((a,b) for (a,b) in ((1,2),(1,3),(2,3)) if
            orient(Tuple(q[1][k] for k in (a,b)),Tuple(q[2][k] for k in (a,b)),Tuple(q[3][k] for k in (a,b)))!=0)
        xy=Tuple(Tuple(v[k] for k in axes) for v in q)
        area= sum(xy[i][1]*xy[mod1(i+1,4)][2]-xy[mod1(i+1,4)][1]*xy[i][2] for i in 1:4)
        edges=Dict{Tuple,Vector{Tuple}}();actual=zero(Q)
        for face in faces
            r=Tuple(Tuple(Q(coordinates.coords[k,v]) for k in axes) for v in face)
            require(all(sign(orient(r[i],r[mod1(i+1,length(r))],r[mod1(i+2,length(r))]))==sign(area) for i in eachindex(r)),
                "macro carrier face folds or has inconsistent outward orientation")
            actual+=sum(r[i][1]*r[mod1(i+1,length(r))][2]-r[mod1(i+1,length(r))][1]*r[i][2] for i in eachindex(r))
            for i in eachindex(face)
                edge=(face[i],face[mod1(i+1,length(face))])
                push!(get!(edges,key(edge),Tuple[]),edge)
            end
        end
        outer=Set{Tuple}()
        for (support,uses) in edges
            require(length(uses) in (1,2),"macro carrier edge has excess incidence")
            if length(uses)==2
                require(uses[1]==reverse(uses[2]),"macro carrier diagonal has equal direction")
            else
                push!(outer,only(uses))
            end
        end
        require(outer==Set((carrier[i],carrier[mod1(i+1,4)]) for i in 1:4) && actual==area,
            "macro carrier faces do not exactly cover their outward boundary")
    end
    return groups,reference
end

function certify(mesh,f,source;oracle=false,native_top=!oracle,check_separation=true)
    src=source_complex(source,f);n=length(f.levels)-1
    expected_nodes=6(n+1)+(f.laterals ? 2 : 0)
    require(nnodes(mesh)==expected_nodes && all(isfinite,mesh.coords),"wrong primary/centroid node count")
    intended=[source.coords[f.axes[3],1]+Float64(t)*f.height for t in f.levels]
    planes=if oracle
        columns=Set(Tuple(source.coords[a,v] for a in f.axes[1:2]) for v in 1:6)
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
        for row in 1:n+1,col in 1:6]
    primary=Set(vec(grid));centers=sort!(collect(setdiff(Set(1:nnodes(mesh)),primary)))
    require(length(primary)==6(n+1) && length(centers)==(f.laterals ? 2 : 0),"column bijection differs")
    center_parent=Dict{Int,Int}()
    for center in centers
        p=exact(point(mesh,center));xy=Tuple(p[a] for a in f.axes[1:2])
        require(min(Q(planes[end-1]),Q(planes[end]))<p[f.axes[3]]<max(Q(planes[end-1]),Q(planes[end])),
            "actual terminal centroid is not strictly between its planes")
        parents=[i for (i,cell) in enumerate(src.cells) if all(sign(orient(src.xy[cell[k]],src.xy[cell[mod1(k+1,4)]],xy))==f.winding for k in 1:4)]
        require(length(parents)==1,"actual terminal centroid lacks one strict macro owner")
        center_parent[center]=only(parents)
    end
    require(Set(values(center_parent))==Set(f.laterals ? (1,2) : ()) ,"terminal centroids do not occupy distinct macros")
    location=Dict(grid[row,col]=>(row,col) for row in 1:n+1,col in 1:6)
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
        parents=[i for i in 1:2 if Set(column_ids)⊆Set(src.cells[i]) && all(center_parent[v]==i for v in body)]
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
        expected=Dict(4=>4,7=>10);n>1 && (expected[5]=2(n-1))
        require(families==expected,"wrong recombined terminal-fan families")
    end
    separated_pairs=0
    for parent in 1:2,layer in 1:n
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
    require(length(boundary)==(f.laterals ? 6n+6 : 12n+6),"wrong exterior typed face count")
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
