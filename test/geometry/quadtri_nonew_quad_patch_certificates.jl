if !isdefined(@__MODULE__, :QuadTriNoNewCertificates)
    include("quadtri_nonew_certificates.jl")
end

module QuadTriNoNewQuadPatchCertificates

using Tessella, SHA
using Tessella.Elements: MixedMesh, ElementBlock, msh_dimension, msh_spec
using Tessella.MeshTypes: nnodes
using ..QuadTriNoNewCertificates
const Quad = QuadTriNoNewCertificates
const Q = Rational{BigInt}
const FRAMES = Dict(:XY=>(1,2,3), :YZ=>(2,3,1), :ZX=>(3,1,2))
const SHAPES = Dict(:unit=>((0.,0.),(1.,0.),(1.,1.),(0.,1.)),
    :skew=>((0.,0.),(2.,0.),(2.5,1.5),(.5,1.5)),
    :trapezoid=>((0.,0.),(2.,0.),(1.5,1.),(.25,1.5)),
    :rounded=>((.1,.2),(2.3,.2),(1.7,1.3),(.25,1.3)))
const FACES = Dict(4=>((1,3,2),(1,2,4),(2,3,4),(3,1,4)),
    5=>((1,4,3,2),(5,6,7,8),(1,2,6,5),(2,3,7,6),(3,4,8,7),(4,1,5,8)),
    7=>((1,4,3,2),(1,2,5),(2,3,5),(3,4,5),(4,1,5)))
require(ok, message) = ok || throw(ArgumentError("Quad-patch certificate: "*message))
position(p) = ntuple(k->iszero(p[k]) ? 0. : Float64(p[k]), 3)
point(mesh, node) = position(Tuple(mesh.coords[:,node]))
exact(p) = Tuple(Q(x) for x in p)
key(nodes) = Tuple(sort!(Int.(collect(nodes))))
cycle(nodes) = minimum(Tuple(circshift(collect(nodes), i)) for i in 0:length(nodes)-1)
det(a,b,c) = a[1]*(b[2]*c[3]-b[3]*c[2])-a[2]*(b[1]*c[3]-b[3]*c[1])+a[3]*(b[1]*c[2]-b[2]*c[1])
orient(a,b,c) = (b[1]-a[1])*(c[2]-a[2])-(b[2]-a[2])*(c[1]-a[1])
blocks(mesh) = mesh isa MixedMesh ? mesh.blocks :
    [ElementBlock(msh,nodes) for (msh,nodes) in ((1,mesh.segs),(2,mesh.tris),(4,mesh.tets)) if !isempty(nodes)]
execute(source; dim=3) = Quad.execute(source; dim)

function fixture(name; plane=:XY, shape=:unit, winding=1, height=1., laterals=false,
                 layers=:three, pins=:auto, curve_reverse=0, tags=:dense,
                 offset=(0.,0.,0.), scale=1., arrangement="Left")
    axes=FRAMES[plane]
    corners=Tuple(ntuple(k->offset[k]+scale*(k==axes[1] ? p[1] : k==axes[2] ? p[2] : 0.),3)
                  for p in SHAPES[shape])
    layer_text,levels=layers===:one ? ("Layers{1}",[0.,1.]) :
        layers===:three ? ("Layers{3}",[0.,1/3,2/3,1.]) :
        layers===:four ? ("Layers{4}",[0.,.25,.5,.75,1.]) :
        layers===:graded ? ("Layers{{2,1},{0.25,1}}",[0.,.125,.25,1.]) :
        layers===:nonbinary ? ("Layers{{1,2},{0.2,1}}",[0.,.2,.6,1.]) : error("unknown layer profile")
    pt,cv,loop,surface=tags===:dense ? ((1,2,3,4),(1,2,3,4),1,1) :
        ((17,71,203,311),(101,205,407,809),431,701)
    points=join(["Point($(pt[i]))={$(join(corners[i],',')),1};" for i in 1:4],"\n")
    lines=join([begin
        a,b=pt[i],pt[mod1(i+1,4)]
        (curve_reverse & (1<<(i-1)))!=0 && ((a,b)=(b,a))
        "Line($(cv[i]))={$a,$b};"
    end for i in 1:4],"\n")
    signed=Tuple(((curve_reverse & (1<<(i-1)))!=0 ? -1 : 1)*cv[i] for i in 1:4)
    loop_nodes=winding>0 ? signed : Tuple(-c for c in reverse(signed))
    pin_order=pins===:opposite ? (3,2,1,4) : pins
    pin_text=pins===:auto ? "" : "={$(join((pt[i] for i in pin_order),','))}"
    delta=ntuple(k->k==axes[3] ? Float64(height) : 0.,3)
    source="""
    Geometry.AutoCoherence=0;Geometry.OldRuledSurface=0;Mesh.ElementOrder=1;
    $points
    $lines
    Curve Loop($loop)={$(join(loop_nodes,','))};Plane Surface($surface)={$loop};
    Transfinite Curve{$(join(cv,','))}=3;
    Transfinite Surface{$surface}$pin_text $arrangement;Recombine Surface{$surface};
    sweep[]=Extrude{$(join(delta,','))}{Surface{$surface};$layer_text;Recombine;
      QuadTriNoNewVerts$(laterals ? " RecombLaterals" : "");};
    """
    return (;name,source,axes,corners,levels,height=Float64(height),laterals,winding,
            plane,shape,pins,curve_reverse,tags,point_tags=pt,curve_tags=cv,surface,
            layer_text,offset,scale,arrangement)
end

function source_complex(source,f)
    require(nnodes(source)==9 && all(isfinite,source.coords),"source must have nine finite nodes")
    require(length(Set(point(source,i) for i in 1:9))==9,"source coordinates coincide")
    cells=Tuple(Tuple(Int.(cell)) for b in blocks(source) for cell in eachcol(b.nodes)
                if msh_dimension(b.msh)==2 && b.msh==3)
    require(length(cells)==4 && sum(size(b.nodes,2) for b in blocks(source))==4,
            "source must contain exactly four Quad4")
    require(all(source.coords[f.axes[3],:].==source.coords[f.axes[3],1]),"source is not exactly axis planar")
    xy=Tuple(Tuple(Q(source.coords[a,i]) for a in f.axes[1:2]) for i in 1:9)
    incidence=Dict{Tuple,Vector{Tuple}}();areas=Q[];used=Set{Int}()
    for cell in cells
        require(length(unique(cell))==4 && all(i->1<=i<=9,cell),"foreign or repeated source node")
        union!(used,cell)
        turns=Tuple(orient(xy[cell[i]],xy[cell[mod1(i+1,4)]],xy[cell[mod1(i+2,4)]]) for i in 1:4)
        require(all(v->sign(v)==f.winding,turns),"actual Quad4 is not strictly convex with CAD winding")
        push!(areas,abs(sum(xy[cell[i]][1]*xy[cell[mod1(i+1,4)]][2]-xy[cell[mod1(i+1,4)]][1]*xy[cell[i]][2] for i in 1:4))/2)
        for i in 1:4
            edge=(cell[i],cell[mod1(i+1,4)])
            push!(get!(incidence,key(edge),Tuple[]),edge)
        end
    end
    require(used==Set(1:9) && length(incidence)==12,"source is not V9/E12/F4")
    boundary=Set(only(uses) for uses in values(incidence) if length(uses)==1)
    shared=Set(edge for (edge,uses) in incidence if length(uses)==2)
    require(length(boundary)==8 && length(shared)==4,"wrong actual source boundary/shared edges")
    require(all(uses->length(uses) in (1,2),values(incidence)),"source edge incidence exceeds two")
    require(all(incidence[e][1]==reverse(incidence[e][2]) for e in shared),"source shared edge directions agree")
    border=Set(v for e in boundary for v in e)
    require(length(border)==8,"source exterior is not eight vertices")
    pivot=only(setdiff(Set(1:9),border))
    require(all(pivot in c for c in cells) && all(pivot in e for e in shared),"source has no common pivot star")
    next=Dict(a=>b for (a,b) in boundary)
    require(length(next)==8 && length(Set(values(next)))==8,"source exterior branches")
    boundary_cycle=Int[minimum(border)]
    for _ in 1:7;push!(boundary_cycle,next[last(boundary_cycle)]);end
    require(length(Set(boundary_cycle))==8 && next[last(boundary_cycle)]==first(boundary_cycle),"source exterior is not one cycle")
    for i in 1:8,j in i+1:8
        a,b=boundary_cycle[i],boundary_cycle[mod1(i+1,8)]
        c,d=boundary_cycle[j],boundary_cycle[mod1(j+1,8)]
        isempty(intersect(Set((a,b)),Set((c,d)))) || continue
        require(!(orient(xy[a],xy[b],xy[c])*orient(xy[a],xy[b],xy[d])<=0 &&
                  orient(xy[c],xy[d],xy[a])*orient(xy[c],xy[d],xy[b])<=0),
                "sampled source boundary self-intersects")
    end
    fan=Q[]
    for (a,b) in boundary
        area=orient(xy[a],xy[b],xy[pivot])/2
        require(sign(area)==f.winding,"pivot is not strictly within the sampled boundary")
        push!(fan,abs(area))
    end
    polygon=abs(sum(xy[a][1]*xy[b][2]-xy[b][1]*xy[a][2] for (a,b) in boundary))/2
    require(sum(fan)==polygon==sum(areas),"pivot fan does not cover the actual source")
    cad=Tuple(only(i for i in 1:9 if point(source,i)==position(p)) for p in f.corners)
    chains=Tuple(begin
        a,b=cad[i],cad[mod1(i+1,4)]
        middle=only(v for v in border if v!=a && v!=b && haskey(incidence,key((a,v))) &&
                    haskey(incidence,key((v,b))) && length(incidence[key((a,v))])==1 && length(incidence[key((v,b))])==1)
        (a,middle,b)
    end for i in 1:4)
    return (;cells,pivot,boundary,shared,boundary_cycle,chains,area=polygon,areas)
end

# Direct reference shape gradients. These exact derivatives do not call the
# production map/certificate, template table, or element-quality routines.
function map_positive(p,msh)
    if msh==4
        return det(p[2].-p[1],p[3].-p[1],p[4].-p[1])>0
    elseif msh==5
        shift=p[5].-p[1]
        require(all(p[i+4].-p[i]==shift for i in 1:4),"Hex8 is not the actual translated Quad product")
        bits=((0,0),(1,0),(1,1),(0,1))
        for (u,v) in ((0,0),(1,0),(1,1),(0,1))
            du=ntuple(d->sum((bits[i][1]==0 ? -1 : 1)*(bits[i][2]==0 ? 1-v : v)*p[i][d] for i in 1:4),3)
            dv=ntuple(d->sum((bits[i][2]==0 ? -1 : 1)*(bits[i][1]==0 ? 1-u : u)*p[i][d] for i in 1:4),3)
            det(du,dv,shift)>0 || return false
        end
        return true
    elseif msh==7
        require(det(p[2].-p[1],p[3].-p[1],p[4].-p[1])==0,"Pyramid5 base is not planar")
        for (u,v) in ((-1,-1),(1,-1),(1,1),(-1,1))
            signs=((-1,-1),(1,-1),(1,1),(-1,1))
            gradients=Tuple((a*(1+b*v)//4,b*(1+a*u)//4,-1//4+a*b*u*v//4) for (a,b) in signs)
            jac=ntuple(axis->ntuple(d->sum(p[i][d]*gradients[i][axis] for i in 1:4)+(axis==3 ? p[5][d] : 0),3),3)
            det(jac...)>0 || return false
        end
        return true
    end
    return false
end

function cell_volume(p,msh)
    origin=p[1]
    sum(det(p[face[1]].-origin,p[face[i]].-origin,p[face[i+1]].-origin)/6
        for face in FACES[msh] for i in 2:length(face)-1)
end

function separated(left,right)
    any(all(det(a,b,q.-p[face[1]])>=0 for q in other)
        for (p,msh,other) in ((left.p,left.msh,right.p),(right.p,right.msh,left.p))
        for face in FACES[msh] for (a,b) in ((p[face[2]].-p[face[1]],p[face[3]].-p[face[1]]),))
end

function certify(mesh,f,source)
    src=source_complex(source,f);n=length(f.levels)-1
    require(nnodes(mesh)==9(n+1) && all(isfinite,mesh.coords),"added/missing/orphan column nodes")
    heights=sort!(unique(mesh.coords[f.axes[3],:]);rev=f.height<0)
    require(length(heights)==n+1 && all(sign.(diff(heights)).==sign(f.height)),"actual planes collapse or reverse")
    require(first(heights)==source.coords[f.axes[3],1],"first layer differs from source")
    for i in eachindex(heights)
        expected=source.coords[f.axes[3],1]+f.levels[i]*f.height
        require(heights[i]==expected || abs(heights[i]-expected)<=4eps(expected),"actual native layer height differs")
    end
    grid=[only(v for v in 1:nnodes(mesh) if mesh.coords[f.axes[3],v]==heights[row] &&
               all(mesh.coords[a,v]==source.coords[a,col] for a in f.axes[1:2])) for row in 1:n+1,col in 1:9]
    require(Set(vec(grid))==Set(1:nnodes(mesh)),"column bijection differs")
    location=Dict(grid[row,col]=>(row,col) for row in 1:n+1,col in 1:9)
    domains=Dict{Tuple,Vector{NamedTuple}}();faces=Dict{Tuple,Vector{Tuple}}()
    used=Set{Int}();edges=Set{Tuple}();families=Dict{Int,Int}();total=zero(Q)
    for block in blocks(mesh),cell in eachcol(block.nodes)
        msh_dimension(block.msh)==3 || continue
        msh=Int(block.msh);require(haskey(FACES,msh),"unexpected typed family")
        require(allunique(cell),"cell repeats actual node IDs")
        levels=[location[Int(v)][1] for v in cell];layer=minimum(levels)
        require(maximum(levels)==layer+1,"cell crosses nonadjacent planes")
        columns=Set(location[Int(v)][2] for v in cell)
        parent=only(i for i in 1:4 if columns⊆Set(src.cells[i]))
        p=Tuple(exact(point(mesh,v)) for v in cell)
        require(map_positive(p,msh),"nonpositive whole P1 reference map")
        volume=cell_volume(p,msh);require(volume>0,"nonpositive actual cell volume")
        total+=volume;union!(used,Int.(cell));families[msh]=get(families,msh,0)+1
        push!(get!(domains,(parent,layer),NamedTuple[]),(;p,msh,volume))
        for pattern in FACES[msh]
            face=Tuple(Int(cell[i]) for i in pattern)
            push!(get!(faces,key(face),Tuple[]),face)
            for i in eachindex(face);push!(edges,key((face[i],face[mod1(i+1,length(face))])));end
        end
    end
    expected=f.laterals ? Dict(5=>4(n-1),7=>12) : Dict(4=>24n-8,7=>4)
    filter!(p->last(p)>0,expected)
    require(families==expected && used==Set(1:nnodes(mesh)),"wrong family counts or unused primary nodes")
    for parent in 1:4,layer in 1:n
        cells=domains[(parent,layer)]
        require(sum(c.volume for c in cells)==src.areas[parent]*abs(Q(heights[layer+1])-Q(heights[layer])),"actual macro is not exactly partitioned")
        for i in 1:length(cells),j in i+1:length(cells)
            require(separated(cells[i],cells[j]),"actual macro cells overlap")
        end
    end
    boundary=Dict{Tuple,Tuple}();internal=Dict{Tuple,Vector{Tuple}}()
    for (face,uses) in faces
        require(length(uses) in (1,2),"typed face incidence exceeds two")
        if length(uses)==2
            require(cycle(uses[1])==cycle(reverse(uses[2])),"shared outward cycles agree")
            internal[face]=uses
        else
            boundary[face]=only(uses)
        end
    end
    require(length(faces)==(f.laterals ? 16n+24 : 56n),"wrong total typed face count")
    require(length(boundary)==(f.laterals ? 8n+12 : 16n+12),"wrong exterior typed face count")
    ncell=sum(values(families))
    require(length(used)-length(edges)+length(faces)-ncell==1,"actual complex is not a ball")
    remaining=Set(keys(boundary))
    for cell in src.cells
        base=Tuple(grid[1,v] for v in cell)
        require(key(base) in remaining,"source Quad4 cap is missing")
        delete!(remaining,key(base))
        pivot=findfirst(==(src.pivot),cell);opposite=cell[mod1(pivot+2,4)]
        others=[v for v in cell if v!=src.pivot && v!=opposite]
        for v in others
            cap=key((grid[n+1,src.pivot],grid[n+1,opposite],grid[n+1,v]))
            require(cap in remaining,"top split does not retain the actual source pivot")
            delete!(remaining,cap)
        end
    end
    for layer in 1:n,(a,b) in src.boundary
        carrier=Set((grid[layer,a],grid[layer,b],grid[layer+1,a],grid[layer+1,b]))
        selected=[face for face in remaining if Set(face)⊆carrier]
        require(length(selected)==(f.laterals ? 1 : 2),"exterior sampled edge has wrong typed split")
        setdiff!(remaining,selected)
    end
    require(isempty(remaining),"unrelated exterior typed face")
    require(total==src.area*abs(Q(last(heights))-Q(first(heights))),"analytic represented swept volume differs")
    return (;source=src,grid,heights,boundary,internal,total,ncell,families,domains)
end

function public_volume(api,entity)
    tags,xyz,_=api.mesh.get_nodes();coords=reshape(xyz,3,:)
    positions=Dict(tag=>i for (i,tag) in enumerate(tags))
    types,_,connectivity=api.mesh.get_elements(3,entity)
    used=sort!(unique(vcat(connectivity...)));remap=Dict(v=>Int32(i) for (i,v) in enumerate(used))
    cells=[ElementBlock(msh,reshape(Int32[remap[v] for v in nodes],msh_spec(msh).nnodes,:))
           for (msh,nodes) in zip(types,connectivity)]
    return MixedMesh(coords[:,[positions[v] for v in used]],cells)
end

function packed_digest(records)
    bytes=UInt8[]
    for record in records
        append!(bytes,record.faces);push!(bytes,record.ncells)
        for cell in record.cells;push!(bytes,cell.msh);append!(bytes,cell.nodes);end
        push!(bytes,record.kind,record.choice)
    end
    return bytes2hex(sha256(bytes))
end

end
