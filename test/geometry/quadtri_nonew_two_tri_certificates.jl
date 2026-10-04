if !isdefined(@__MODULE__,:QuadTriNoNewCertificates)
    include("quadtri_nonew_certificates.jl")
end

module QuadTriNoNewTwoTriCertificates

using Tessella, SHA
using Tessella.Elements: ElementBlock, MixedMesh, msh_dimension, msh_spec
using Tessella.MeshTypes: Mesh, nnodes
using ..QuadTriNoNewCertificates
const Quad=QuadTriNoNewCertificates
const Q=Rational{BigInt}
const SHAPES=Dict(:unit=>((0.,0.),(1.,0.),(1.,1.),(0.,1.)),
    :skew=>((0.,0.),(2.,0.),(2.5,1.5),(.5,1.5)),
    :trapezoid=>((0.,0.),(2.,0.),(1.5,1.),(.25,1.5)))
const FRAMES=Dict(:XY=>(1,2,3),:YZ=>(2,3,1),:ZX=>(3,1,2),:XZ=>(1,3,2))
const FACES=Dict(4=>((1,3,2),(1,2,4),(2,3,4),(3,1,4)),
    6=>((1,3,2),(4,5,6),(1,2,5,4),(2,3,6,5),(3,1,4,6)))
require(ok,message)=ok || throw(ArgumentError("Two-Tri NoNew certificate: "*message))
position(p)=ntuple(k->iszero(p[k]) ? 0. : Float64(p[k]),3)
point(mesh,node)=position(Tuple(mesh.coords[:,node]))
exact(p)=Tuple(Q(x) for x in p)
key(nodes)=Tuple(sort!(Int.(collect(nodes))))
cycle(nodes)=minimum(Tuple(circshift(collect(nodes),i)) for i in 0:length(nodes)-1)
orient(a,b,c)=(b[1]-a[1])*(c[2]-a[2])-(b[2]-a[2])*(c[1]-a[1])
det(a,b,c)=a[1]*(b[2]*c[3]-b[3]*c[2])-a[2]*(b[1]*c[3]-b[3]*c[1])+a[3]*(b[1]*c[2]-b[2]*c[1])
blocks(mesh)=Quad.blocks(mesh)

function fixture(name;plane=:XY,shape=:unit,winding=1,height=1.,laterals=false,
                 layers=:three,arrangement="Left",pins=:auto,offset=(0.,0.,0.))
    axes=FRAMES[plane]
    corners=Tuple(ntuple(k->offset[k]+(k==axes[1] ? p[1] : k==axes[2] ? p[2] : 0.),3) for p in SHAPES[shape])
    layer_text,levels=layers===:one ? ("Layers{1}",[0.,1.]) :
        layers===:three ? ("Layers{3}",[0.,1/3,2/3,1.]) :
        layers===:graded ? ("Layers{{2,1},{0.25,1}}",[0.,.125,.25,1.]) :
        layers===:nonbinary ? ("Layers{{1,2},{0.2,1}}",[0.,.2,.6,1.]) : error("unknown layers")
    points=join(["Point($i)={$(join(corners[i],',')),1};" for i in 1:4],"\n")
    loop=winding>0 ? "1,2,3,4" : "-4,-3,-2,-1"
    pin_text=pins===:auto ? "" : pins===:opposite ? "={3,2,1,4}" : "={$(join(pins,','))}"
    delta=ntuple(k->k==axes[3] ? Float64(height) : 0.,3)
    source="""
    Geometry.AutoCoherence=0;
    $points
    Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
    Curve Loop(1)={$loop};Plane Surface(1)={1};
    Transfinite Curve{1,2,3,4}=2;
    Transfinite Surface{1}$pin_text $arrangement;
    sweep[]=Extrude{$(join(delta,','))}{Surface{1};$layer_text;Recombine;
      QuadTriNoNewVerts$(laterals ? " RecombLaterals" : "");};
    """
    return (;name,source,axes,corners,levels,height=Float64(height),laterals,winding,plane,shape,arrangement,pins)
end

function fixtures()
    result=NamedTuple[]
    for (pi,plane) in enumerate((:XY,:YZ,:ZX)),winding in (-1,1),direction in (-1,1),laterals in (false,true),(li,layers) in enumerate((:one,:three,:graded))
        shape=(:unit,:skew,:trapezoid)[mod1(pi+li-1,3)]
        push!(result,fixture("two_tri_$(plane)_$(shape)_$(layers)_W$(winding)_D$(direction)_R$(laterals)";
            plane,shape,winding,height=direction*1.5,laterals,layers,
            arrangement=li==1 ? "Left" : li==2 ? "Right" : "AlternateLeft"))
    end
    return result
end

execute(source;dim=3)=Quad.execute(source;dim)

function source_complex(source,f)
    require(source isa Mesh && size(source.tris)==(3,2) && nnodes(source)==4,"source must have exactly two Tri3 and four nodes")
    require(all(isfinite,source.coords),"nonfinite source")
    cad_to_dense=Tuple(only(i for i in 1:4 if point(source,i)==position(p)) for p in f.corners)
    require(length(unique(cad_to_dense))==4,"source CAD corner identities differ")
    source_points=[Tuple(Q(source.coords[axis,i]) for axis in f.axes[1:2]) for i in 1:4]
    cad_cycle=f.winding>0 ? cad_to_dense : (cad_to_dense[1],cad_to_dense[4],cad_to_dense[3],cad_to_dense[2])
    turns=[orient(source_points[cad_cycle[i]],source_points[cad_cycle[mod1(i+1,4)]],source_points[cad_cycle[mod1(i+2,4)]]) for i in 1:4]
    require(all(x->sign(x)==f.winding,turns),"CAD source is not strictly convex")
    cells=Tuple(Tuple(Int.(cell)) for cell in eachcol(source.tris))
    incidence=Dict{Tuple,Vector{Tuple}}()
    signs=Q[]
    for cell in cells
        require(length(unique(cell))==3,"repeated source ID")
        push!(signs,orient((source_points[node] for node in cell)...))
        for i in 1:3
            edge=(cell[i],cell[mod1(i+1,3)])
            push!(get!(incidence,key(edge),Tuple[]),edge)
        end
    end
    require(all(x->sign(x)==f.winding,signs),"source winding differs from CAD cycle")
    require(length(incidence)==5,"source edge complex is not V4/E5/F2")
    interior=[edge for (edge,uses) in incidence if length(uses)==2]
    require(length(interior)==1,"source has no unique shared diagonal")
    diagonal=only(interior)
    require(incidence[diagonal][1]==reverse(incidence[diagonal][2]),"source shared-edge orientation differs")
    boundary=Set(only(uses) for uses in values(incidence) if length(uses)==1)
    require(boundary==Set((cad_cycle[i],cad_cycle[mod1(i+1,4)]) for i in 1:4),"source boundary differs from CAD incidence")
    opposite=[i for i in 1:4 if !(i in diagonal)]
    require(orient(source_points[diagonal[1]],source_points[diagonal[2]],source_points[opposite[1]])*
        orient(source_points[diagonal[1]],source_points[diagonal[2]],source_points[opposite[2]])<0,"source triangles overlap across diagonal")
    polygon=abs(sum(source_points[cad_cycle[i]][1]*source_points[cad_cycle[mod1(i+1,4)]][2]-
        source_points[cad_cycle[i]][2]*source_points[cad_cycle[mod1(i+1,4)]][1] for i in 1:4))/2
    require(sum(abs,signs)/2==polygon,"source triangles do not cover convex polygon")
    return (;cells,diagonal,boundary,cad_to_dense,area=polygon,triangle_areas=abs.(signs)./2)
end

function certify(mesh,f,source)
    src=source_complex(source,f)
    nlevels=length(f.levels);n=nlevels-1
    require(nnodes(mesh)==4nlevels,"added, missing or orphan volume nodes")
    require(all(isfinite,mesh.coords),"nonfinite volume coordinates")
    heights=sort!(unique(mesh.coords[f.axes[3],:]);rev=f.height<0)
    require(length(heights)==nlevels,"wrong number of actual layer planes")
    require(all(i->sign(heights[i+1]-heights[i])==sign(f.height),1:n),"actual planes collapse or reverse")
    require(heights[1]==source.coords[f.axes[3],1],"first plane differs from actual source")
    for i in 1:nlevels
        expected=source.coords[f.axes[3],1]+f.levels[i]*f.height
        require(abs(heights[i]-expected)<=8eps(Float64)*max(1.,abs(expected)),"actual layer height differs")
    end
    grid=Matrix{Int}(undef,nlevels,4)
    for row in 1:nlevels,column in 1:4
        p=point(source,column)
        matches=[node for node in 1:nnodes(mesh) if mesh.coords[f.axes[3],node]==heights[row] &&
            all(axis->mesh.coords[axis,node]==p[axis],f.axes[1:2])]
        require(length(matches)==1,"source column missing or duplicated")
        grid[row,column]=only(matches)
    end
    require(Set(vec(grid))==Set(1:nnodes(mesh)),"unrelated volume nodes")
    location=Dict(grid[row,column]=>(row,column) for row in 1:nlevels,column in 1:4)
    faces=Dict{Tuple,Vector{Tuple}}();used=Set{Int}();edges=Set{Tuple}()
    per_domain=Dict{Tuple{Int,Int},Int}();domain_volume=Dict{Tuple{Int,Int},Q}()
    total=zero(Q);ncell=0
    for block in blocks(mesh)
        msh_dimension(block.msh)==3 || continue
        require(block.msh==(f.laterals ? 6 : 4),"unexpected two-Tri sweep family")
        for cell in eachcol(block.nodes)
            ncell+=1;union!(used,Int.(cell))
            require(length(unique(cell))==length(cell),"repeated actual cell ID")
            rows=[location[Int(node)][1] for node in cell]
            layer=minimum(rows)
            require(maximum(rows)==layer+1,"cell spans unrelated planes")
            columns=Set(location[Int(node)][2] for node in cell)
            parents=[i for i in 1:2 if columns⊆Set(src.cells[i])]
            require(length(parents)==1,"cell crosses the actual source diagonal")
            parent=only(parents);domain=(parent,layer)
            per_domain[domain]=get(per_domain,domain,0)+1
            p=Tuple(exact(point(mesh,node)) for node in cell)
            if block.msh==4
                jac=det(p[2].-p[1],p[3].-p[1],p[4].-p[1])
                require(jac>0,"nonpositive full Tet4 Jacobian")
                volume=jac/6
            else
                delta=p[4].-p[1]
                require(all(i->p[i+3].-p[i]==delta,1:3),"Pri6 is not an affine logical product")
                jac=det(p[2].-p[1],p[3].-p[1],delta)
                require(jac>0,"nonpositive whole-reference Pri6 Jacobian")
                volume=jac/2
            end
            total+=volume;domain_volume[domain]=get(domain_volume,domain,zero(Q))+volume
            for pattern in FACES[Int(block.msh)]
                face=Tuple(Int(cell[i]) for i in pattern)
                push!(get!(faces,key(face),Tuple[]),face)
                for i in eachindex(face);push!(edges,key((face[i],face[mod1(i+1,length(face))])));end
            end
        end
    end
    require(ncell==(f.laterals ? 2n : 6n),"wrong actual volume cell count")
    require(used==Set(1:nnodes(mesh)),"unused actual volume corner")
    for parent in 1:2,layer in 1:n
        domain=(parent,layer)
        require(get(per_domain,domain,0)==(f.laterals ? 1 : 3),"logical prism cell count differs")
        expected=src.triangle_areas[parent]*abs(Q(heights[layer+1])-Q(heights[layer]))
        require(domain_volume[domain]==expected,"logical prism is not exactly covered")
    end
    boundary=Dict{Tuple,Tuple}();internal=Dict{Tuple,Vector{Tuple}}()
    for (face,uses) in faces
        require(length(uses) in (1,2),"actual typed face incidence exceeds two")
        if length(uses)==2
            require(cycle(uses[1])==cycle(reverse(uses[2])),"internal cyclic orientations are not opposite")
            internal[face]=uses
        else
            boundary[face]=only(uses)
        end
    end
    require(length(boundary)==(f.laterals ? 4n+4 : 8n+4),"wrong exterior typed face count")
    require(length(internal)==(f.laterals ? 3n-2 : 8n-2),"wrong internal typed face count")
    require(length(used)-length(edges)+length(faces)-ncell==1,"volume is not a combinatorial ball")
    remaining=Set(keys(boundary))
    for row in (1,nlevels),cell in src.cells
        cap=Tuple(grid[row,node] for node in cell)
        require(key(cap) in remaining,"actual source/copy cap triangulation changed")
        delete!(remaining,key(cap))
        p=Tuple(exact(point(source,node)) for node in cell)
        step=exact(point(mesh,grid[nlevels,cell[1]])).-p[1]
        sigma=sign(det(p[2].-p[1],p[3].-p[1],step))
        same=cycle(cap)==cycle(boundary[key(cap)])
        require(same==(row==1 ? sigma<0 : sigma>0),"cap outward orientation differs")
    end
    for layer in 1:n,(a,b) in src.boundary
        corners=Set((grid[layer,a],grid[layer,b],grid[layer+1,a],grid[layer+1,b]))
        selected=[face for face in remaining if Set(face)⊆corners]
        require(length(selected)==(f.laterals ? 1 : 2),"external side has wrong typed split")
        require(all(face->length(face)==(f.laterals ? 4 : 3),selected),"external side has wrong family")
        setdiff!(remaining,selected)
    end
    require(isempty(remaining),"unexpected exterior surface")
    for layer in 1:n
        a,b=src.diagonal
        corners=Set((grid[layer,a],grid[layer,b],grid[layer+1,a],grid[layer+1,b]))
        selected=[face for face in keys(internal) if Set(face)⊆corners]
        require(length(selected)==(f.laterals ? 1 : 2),"shared source diagonal swept face is cracked")
    end
    expected_total=src.area*abs(Q(last(heights))-Q(first(heights)))
    require(total==expected_total,"analytic actual-plane product volume differs")
    return (;source=src,grid,heights,boundary,internal,total,ncell)
end

function public_volume(api,entity=1)
    tags,xyz,_=api.mesh.get_nodes();coordinates=reshape(xyz,3,:)
    positions=Dict(tag=>i for (i,tag) in enumerate(tags))
    types,_,families=api.mesh.get_elements(3,entity)
    used=sort!(unique(vcat(families...)))
    remap=Dict(tag=>Int32(i) for (i,tag) in enumerate(used))
    cells=ElementBlock[ElementBlock(msh,reshape(Int32[remap[tag] for tag in nodes],msh_spec(msh).nnodes,:)) for (msh,nodes) in zip(types,families)]
    return MixedMesh(coordinates[:,[positions[tag] for tag in used]],cells)
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
