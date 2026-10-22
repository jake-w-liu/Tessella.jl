module ModelNativeInitialTet
using ...Predicates: orient3,insphere
include("ModelNativeHilbertIntervals.jl")

const FACE=((1,2,3),(2,4,3),(3,4,1),(2,1,4))
mutable struct Cell
    vertices::NTuple{4,Int32}
    neighbors::NTuple{4,Int32}
    marked::Bool
end
mutable struct Kernel
    points::Vector{NTuple{3,Float64}}
    cells::Vector{Cell}
    count::Int
    random_state::UInt32
    touched::Vector{Int32}
end
@inline face(k,t,s)=ntuple(i->k.cells[t].vertices[FACE[s][i]],3)
@inline function key(f)
    a,b,c=f
    a>b && ((a,b)=(b,a));b>c && ((b,c)=(c,b));a>b && ((a,b)=(b,a))
    return (a,b,c)
end
@inline orientation(k,f,p)=orient3(k.points[f[1]],k.points[f[2]],k.points[f[3]],p)
@inline function determinant(a,b,c,d)
    ax=a[1]-d[1];bx=b[1]-d[1];cx=c[1]-d[1]
    ay=a[2]-d[2];by=b[2]-d[2];cy=c[2]-d[2]
    az=a[3]-d[3];bz=b[3]-d[3];cz=c[3]-d[3]
    return ax*(by*cz-bz*cy)+bx*(cy*az-cz*ay)+cx*(ay*bz-az*by)
end
function value(k,f,p)
    a=k.points[f[1]];b=k.points[f[2]];c=k.points[f[3]]
    v=determinant(a,b,c,p);sign=orient3(a,b,c,p)
    sign==0 && return 0.0
    signbit(v)==(sign<0) && !iszero(v) && return v
    R=Rational{BigInt}
    exact=determinant(ntuple(i->R(a[i]),3),ntuple(i->R(b[i]),3),
                      ntuple(i->R(c[i]),3),ntuple(i->R(p[i]),3))
    return Float64(exact)
end
function inside(k,t,node)
    f=k.cells[t].vertices;p=k.points
    sign=insphere(p[f[1]],p[f[2]],p[f[3]],p[f[4]],p[node])
    sign!=0 && return sign>0
    rows=(f...,Int32(node));swaps=0
    for last in 4:-1:1
        for i in 1:last
            if rows[i]>rows[i+1]
                a=rows[i];b=rows[i+1]
                rows=Base.setindex(Base.setindex(rows,b,i),a,i+1);swaps+=1
            end
        end
    end
    sign=orient3(p[rows[2]],p[rows[3]],p[rows[4]],p[rows[5]])
    if sign==0
        sign=-orient3(p[rows[1]],p[rows[3]],p[rows[4]],p[rows[5]])
        sign==0 && throw(ArgumentError("native initial tetrahedralization: unresolved symbolic orientation"))
    end
    return (isodd(swaps) ? -sign : sign)>0
end
@inline function rand!(k)
    k.random_state=k.random_state*UInt32(214013)+UInt32(2531011)
    return (k.random_state>>16)&UInt32(32767)
end
function bounds(points,scale)
    lo=ntuple(k->minimum(p->p[k],points),3)
    hi=ntuple(k->maximum(p->p[k],points),3)
    center=ntuple(k->(lo[k]+hi[k])*.5,3)
    return ntuple(k->(lo[k]-center[k])*scale+center[k],3),
           ntuple(k->(hi[k]-center[k])*scale+center[k],3)
end
function hilbert(points,count)
    order=collect(Int32,1:count)
    lo,hi=bounds(@view(points[1:count]),1.01)
    return HilbertIntervals.order!(order,points,lo,hi)
end
function connect!(k,t,s,pending)
    f=key(face(k,t,s));index=findfirst(c->c[1]==f,pending)
    if index===nothing
        push!(pending,(f,Int32(t),Int32(s)))
    else
        _,other,os=pending[index];deleteat!(pending,index)
        k.cells[t].neighbors=Base.setindex(k.cells[t].neighbors,other,s)
        k.cells[other].neighbors=Base.setindex(k.cells[other].neighbors,Int32(t),os)
    end
    return nothing
end
function Kernel(original;seed=1,factor=1e-12)
    points=NTuple{3,Float64}[Tuple(p) for p in original];count=length(points)
    4<=count<=typemax(Int32)-8 || throw(ArgumentError("native initial tetrahedralization: invalid point count"))
    all(p->all(isfinite,p),points) || throw(ArgumentError("native initial tetrahedralization: nonfinite input"))
    k=Kernel(points,Cell[],count,UInt32(seed),Int32[]);sizehint!(k.cells,10count)
    magnitude=ntuple(i->maximum(p->abs(p[i]),points),3)
    distance=sqrt(magnitude[1]^2+magnitude[2]^2+magnitude[3]^2)
    for i in eachindex(points)
        p=points[i]
        points[i]=ntuple(axis->p[axis]+distance*factor*Float64(rand!(k))/32767.,3)
    end
    lo,hi=bounds(points,1.3)
    corners=((0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1))
    append!(points,[ntuple(axis->corner[axis]==1 ? hi[axis] : lo[axis],3) for corner in corners])
    for row in ((8,3,4,2),(8,1,2,4),(2,7,8,3),(2,1,8,5),(5,2,6,8),(8,2,6,7))
        v=ntuple(i->Int32(count+row[i]),4)
        s=orient3(points[v[1]],points[v[2]],points[v[3]],points[v[4]])
        s<0 && (v=(v[2],v[1],v[3],v[4]))
        s==0 && throw(ArgumentError("native initial tetrahedralization: degenerate enclosing cube"))
        push!(k.cells,Cell(v,(0,0,0,0),false))
    end
    pending=Tuple{NTuple{3,Int32},Int32,Int32}[]
    for s in 1:4,t in 1:6
        connect!(k,t,s,pending)
    end
    return k
end
function walk(k,t,node)
    investigated=Set{Int32}([t]);queue=Int32[];position=1
    for _ in 1:100length(k.cells)
        t==0 && return Int32(0)
        minimum=0.;neighbor=0;count=0
        for s in 1:4
            v=value(k,face(k,t,s),k.points[node]);count+=v>=0
            if v<minimum
                other=k.cells[t].neighbors[s]
                if !(other in investigated)
                    neighbor=s;minimum=v
                else
                    push!(queue,other)
                end
            end
        end
        count==4 && inside(k,t,node) && return Int32(t)
        if neighbor!=0
            t=k.cells[t].neighbors[neighbor];push!(investigated,t)
        elseif position<=length(queue)
            t=queue[position];position+=1
        else
            return Int32(0)
        end
    end
    throw(ArgumentError("native initial tetrahedralization: walk did not terminate"))
end
@inline function mark!(k,t)
    if !k.cells[t].marked
        k.cells[t].marked=true;push!(k.touched,Int32(t))
    end
    return nothing
end
function cavity(k,start,node)
    cells=Int32[];boundary=Tuple{NTuple{3,Int32},Int32,Int32}[]
    stack=NTuple{3,Int32}[(0,start,1)]
    for t in k.touched
        k.cells[t].marked=false
    end
    empty!(k.touched)
    while !isempty(stack)
        previous,t,first=pop!(stack)
        if first==1
            mark!(k,t);push!(cells,t)
        end
        for s in first:4
            other=k.cells[t].neighbors[s]
            if other==0
                push!(boundary,(face(k,t,s),Int32(s),other))
            elseif other==previous
                continue
            elseif !inside(k,other,node)
                push!(boundary,(face(k,t,s),Int32(s),other));mark!(k,other)
            elseif !k.cells[other].marked
                push!(stack,(previous,t,Int32(s+1)));push!(stack,(t,other,1));break
            end
        end
    end
    return cells,boundary
end
function repair!(k,cells,boundary,node)
    negatives=[i for i in eachindex(boundary) if orientation(k,boundary[i][1],k.points[node])<=0]
    isempty(negatives) && return true
    containing=findfirst(t->all(s->orientation(k,face(k,t,s),k.points[node])>=0,1:4),cells)
    containing===nothing && return false
    root=cells[containing]
    while !isempty(negatives)
        for i in negatives
            bad_key=key(boundary[i][1]);index=findfirst(t->any(s->key(face(k,t,s))==bad_key,1:4),cells)
            index===nothing || deleteat!(cells,index)
        end
        root in cells || return false
        connected=Int32[root];stack=Int32[root]
        while !isempty(stack)
            t=pop!(stack)
            for other in k.cells[t].neighbors
                if other!=0 && other in cells && !(other in connected)
                    push!(connected,other);push!(stack,other)
                end
            end
        end
        if length(connected)!=length(cells)
            empty!(cells);append!(cells,connected)
        end
        empty!(boundary)
        for t in cells,s in 1:4
            f=face(k,t,s);index=findfirst(c->key(c[1])==key(f),boundary)
            if index===nothing
                push!(boundary,(f,Int32(s),k.cells[t].neighbors[s]))
            else
                deleteat!(boundary,index)
            end
        end
        negatives=[i for i in eachindex(boundary) if orientation(k,boundary[i][1],k.points[node])<=0]
    end
    return true
end
function random_cell!(k)
    while true
        t=Int(rand!(k)%UInt32(length(k.cells)))+1
        k.cells[t].vertices[1]!=0 && return Int32(t)
    end
end
function run!(k)
    choice=random_cell!(k)
    pending=Tuple{NTuple{3,Int32},Int32,Int32}[]
    for node in hilbert(k.points,k.count)
        k.cells[choice].vertices[1]==0 && (choice=random_cell!(k))
        seed=walk(k,choice,node)
        while seed==0
            choice=random_cell!(k);seed=walk(k,choice,node)
        end
        cells,boundary=cavity(k,seed,node)
        repair!(k,cells,boundary,node) || throw(ArgumentError("native initial tetrahedralization: invalid cavity"))
        choice=cells[1];empty!(pending)
        for i in eachindex(boundary)
            f,_,other=boundary[i]
            t=i<=length(cells) ? cells[i] : Int32(length(k.cells)+1)
            if t>length(k.cells)
                push!(k.cells,Cell((0,0,0,0),(0,0,0,0),false))
            end
            k.cells[t].vertices=(f...,node);k.cells[t].neighbors=(other,0,0,0)
            if other!=0
                os=findfirst(s->key(face(k,other,s))==key(f),1:4)
                os===nothing && throw(ArgumentError("native initial tetrahedralization: inconsistent cavity incidence"))
                k.cells[other].neighbors=Base.setindex(k.cells[other].neighbors,t,os)
            end
            for s in 2:4
                connect!(k,t,s,pending)
            end
        end
        for i in length(boundary)+1:length(cells)
            k.cells[cells[i]].vertices=(0,0,0,0)
        end
    end
    output=NTuple{4,Int32}[c.vertices for c in k.cells if c.vertices[1]!=0 && all(i->i<=k.count,c.vertices)]
    result=Matrix{Int32}(undef,4,length(output))
    for (column,cell) in enumerate(output),axis in 1:4
        result[axis,column]=cell[axis]
    end
    return result
end
end

import ..RecoverCDT

# Receipt of the volume-boundary subdivision write-back — see
# `_model_native_volume_boundary_writeback!` below. `constraints` is the
# literal-certificate-ready subdivided triple list; `nodes`/`node_params`
# assign each recovery-born mesh row its owning entity and parameter;
# `surfaces`/`curves` carry each product's complete post-write-back cell list
# as `(msh, connectivity, label)`; `replaced` lists the dropped source cell
# connectivity per owner so the dim-01 assembly can remove them atomically.
struct _ModelVolumeBoundaryWriteback
    constraints
    nodes::Dict{Int32,Tuple{Int,Int32}}
    node_params::Dict{Int32,Vector{Float64}}
    surfaces::Dict{Int,Vector{Tuple{Int,Vector{Int32},UInt64}}}
    curves::Dict{Int,Vector{Tuple{Int,Vector{Int32},UInt64}}}
    replaced::Dict{Tuple{Int,Int32},Set{Tuple}}
    added::Dict{Tuple{Int,Int32},Vector{Tuple{Int,Vector{Int32},UInt64,Int32}}}
end

# Input rows follow the actual MVertex number ordering. Perturbed working
# coordinates stay inside the kernel; emitted vertices retain their births.
mutable struct _ModelNativeVolumeBirths
    source::Union{Nothing,Mesh}
    node_tags::Vector{UInt64}
    node_high_water::UInt64
    temporary_elements::UInt64
    element_tags::Vector{UInt64}
    element_high_water::UInt64
    cell_source::Union{Nothing,Mesh}
    boundary_writeback::Union{Nothing,_ModelVolumeBoundaryWriteback}
end
_ModelNativeVolumeBirths(node_max::UInt64=UInt64(0),element_max::UInt64=UInt64(0))=
    _ModelNativeVolumeBirths(nothing,UInt64[],node_max,0,UInt64[],element_max,nothing,nothing)

function _model_native_volume_cell_births!(births,mesh)
    births===nothing && return mesh
    count=UInt64(ntets(mesh))
    last=Base.checked_add(births.element_high_water,count)
    first=births.element_high_water+UInt64(1)
    births.element_tags=collect(first:last)
    births.element_high_water=last
    births.cell_source=mesh
    return mesh
end

# A graft/recovery pass splices a new tetrahedron complex in place of the
# initial fill — the same vertex rows stay bound but the recorded cell source
# no longer names the mesh, which would silently deactivate the refinement
# receipt and leave refinement-born nodes without UIDs. Re-bind the receipt:
# fresh node UIDs for any appended Steiner rows and fresh cell numbers for
# the rebuilt complex.
function _model_native_volume_rebind!(births,mesh)
    births===nothing && return mesh
    births.source===nothing && return mesh
    births.cell_source===mesh && return mesh
    extra=nnodes(mesh)-length(births.node_tags)
    if extra>0
        first_uid=births.node_high_water+UInt64(1)
        births.node_high_water=Base.checked_add(births.node_high_water,
            UInt64(extra))
        append!(births.node_tags,first_uid:births.node_high_water)
    end
    _model_native_volume_cell_births!(births,mesh)
    return mesh
end

struct _ModelVolumeGeometryProduction
    mesh::Union{Mesh,MixedMesh}
    boundary_node_tags::Vector{UInt64}
    node_high_water::UInt64
    temporary_elements::UInt64
    element_tags::Vector{UInt64}
    element_high_water::UInt64
    boundary_writeback::Union{Nothing,_ModelVolumeBoundaryWriteback}
end

function _model_volume_geometry_production(mesh,births,caller)
    labels=zeros(UInt64,nnodes(mesh))
    if births.source!==nothing
        source=births.source
        nnodes(mesh)>=nnodes(source) &&
            @view(mesh.coords[:,1:nnodes(source)])==source.coords || throw(ArgumentError(
                "$caller: volume recovery lost its explicit boundary row bindings"))
        length(births.node_tags) in (nnodes(source),nnodes(mesh)) || throw(ArgumentError(
            "$caller: born volume vertex bindings do not cover producer rows"))
        labels[1:length(births.node_tags)].=births.node_tags
    end
    cells=if births.cell_source!==nothing && births.cell_source.tets==mesh.tets
        length(births.element_tags)==ntets(mesh) && allunique(births.element_tags) &&
            all(tag->0<tag<=births.element_high_water,births.element_tags) || throw(ArgumentError(
            "$caller: born volume cell bindings do not cover producer cells"))
        copy(births.element_tags)
    else
        UInt64[]
    end
    return _ModelVolumeGeometryProduction(mesh,labels,births.node_high_water,
        births.temporary_elements,cells,births.element_high_water,
        births.boundary_writeback)
end

function _model_native_initial_delaunay_fill(surface::Mesh,order,options,caller;
                                             births=nothing,kernel_output=nothing)
    length(order)==nnodes(surface) && sort(collect(order))==collect(1:nnodes(surface)) ||
        throw(ArgumentError("$caller: native initial boundary order is not a complete row permutation"))
    points=NTuple{3,Float64}[_model_mesh_coordinate(surface,row) for row in order]
    state=options===nothing ? Ref(UInt32(1)) : options.random_state
    factor=options===nothing ? 1e-12 : options.random_factor3d
    kernel=ModelNativeInitialTet.Kernel(points;seed=state[],factor)
    cells=try
        ModelNativeInitialTet.run!(kernel)
    finally
        state[]=kernel.random_state
    end
    if births!==nothing
        # removeBox=false creates eight explicitly numbered temporary MVertices
        # and one MTetrahedron for each live allocator row, including the cube.
        births.node_high_water=max(births.node_high_water,UInt64(kernel.count+8))
        temporary=UInt64(count(cell->cell.vertices[1]!=0,kernel.cells))
        births.temporary_elements=Base.checked_add(births.temporary_elements,temporary)
        births.element_high_water=Base.checked_add(births.element_high_water,temporary)
    end
    kernel_output===nothing || (kernel_output[]=kernel)
    mapped=Matrix{Int32}(undef,4,size(cells,2))
    for column in axes(cells,2)
        vertices=ntuple(slot->Int32(order[cells[slot,column]]),4)
        sign=ModelNativeInitialTet.orient3(
            (_model_mesh_coordinate(surface,row) for row in vertices)...)
        sign>0 && (vertices=(vertices[2],vertices[1],vertices[3],vertices[4]))
        for slot in 1:4
            mapped[slot,column]=vertices[slot]
        end
    end
    return Mesh(surface.coords;tets=mapped)
end

# Gmsh's initial set includes actual embedded GVertex/GEdge/GFace products.
# Append by born UID, preserving the boundary prefix and the producer rows.
function _model_native_initial_embedded_input(m,volume,geometry,authority,caller)
    labels=geometry.native_node_tags
    labels===nothing && return geometry
    embeds=get(m.embeds,(3,volume),NTuple{2,Int}[])
    isempty(embeds) && return geometry
    authority isa _ModelMeshAuthority || throw(ArgumentError(
        "$caller: native embedded inputs require actual mesh authority"))
    coordinates=NTuple{3,Float64}[
        _model_mesh_coordinate(geometry.surface,row) for row in 1:nnodes(geometry.surface)]
    tags=copy(labels);rows=Dict(tag=>Int32(row) for (row,tag) in enumerate(tags))
    for (dimension,tag) in embeds
        product=_model_mesh_product(authority,dimension,tag;caller,producer_order=true)
        product isa MixedMesh && product.entity_data!==nothing || throw(ArgumentError(
            "$caller: embedded entity ($dimension, $tag) has no actual generated product"))
        for row in 1:nnodes(product)
            uid=product.entity_data.external_node_tags[row]
            uid>0 || throw(ArgumentError("$caller: embedded vertex has no positive born UID"))
            point=_model_mesh_coordinate(product,row)
            previous=get(rows,uid,Int32(0))
            if previous==0
                length(coordinates)<typemax(Int32) || throw(ArgumentError(
                    "$caller: native initial inputs exceed the Int32 node limit"))
                push!(coordinates,point);push!(tags,uid);rows[uid]=Int32(length(tags))
            else
                coordinates[previous]==point || throw(ArgumentError(
                    "$caller: actual embedded vertex $uid has conflicting coordinates"))
            end
        end
    end
    points=Matrix{Float64}(undef,3,length(coordinates))
    for (row,point) in enumerate(coordinates),axis in 1:3
        points[axis,row]=point[axis]
    end
    input=Mesh(points;tris=geometry.surface.tris)
    return merge(geometry,(surface=input,native_node_tags=tags,native_order=sortperm(tags)))
end


function _model_native_recovered_volume_fill(surface,order,kernel,caller;protected_facets=nothing)
    inverse=Vector{Int32}(undef,length(order))
    for (native,row) in enumerate(order);inverse[row]=Int32(native);end
    facets=NTuple{3,Int32}[(inverse[surface.tris[1,column]],
        inverse[surface.tris[2,column]],inverse[surface.tris[3,column]])
        for column in axes(surface.tris,2)]
    recovered=ModelNativeTetgenArena.recover_convex(kernel,facets)
    if protected_facets!==nothing
        for face in protected_facets
            bound=ntuple(i->inverse[face[i]],3)
            ModelNativeTetgenArena.recoverface!(recovered.arena,bound)
            for i in 1:3
                push!(recovered.arena.protected_segments,
                    ModelNativeTetgenArena.edgekey(bound[i],bound[mod1(i+1,3)]))
            end
        end
    end
    ModelNativeTetgenArena.improve_quality32!(recovered.arena)
    ModelNativeTetgenArena.consistency(recovered.arena) || throw(ArgumentError(
        "$caller: native recovered quality flips lost reciprocal face bonds"))
    ordered_cells,_=ModelNativeTetgenArena.real_cells(recovered.arena)
    cells=Matrix{Int32}(undef,4,length(ordered_cells))
    for (column,cell) in enumerate(ordered_cells)
        vertices=ntuple(slot->Int32(order[cell[slot]]),4)
        sign=ModelNativeInitialTet.orient3(
            (_model_mesh_coordinate(surface,row) for row in vertices)...)
        sign>0 && (vertices=(vertices[2],vertices[1],vertices[3],vertices[4]))
        for slot in 1:4;cells[slot,column]=vertices[slot];end
    end
    return Mesh(surface.coords;tets=cells)
end

# Actual lower cell labels bind their vertices through UID, independently of
# CAD proximity, public renumbering, producer order, and geometric coverage.
function _model_native_volume_input_constraints(m,volume,geometry,authority,caller)
    tags=geometry.native_node_tags
    tags!==nothing && allunique(tags) && length(tags)==nnodes(geometry.surface) ||
        throw(ArgumentError("$caller: actual Volume[$volume] lacks unique initial row bindings"))
    rows=Dict(uid=>Int32(row) for (row,uid) in enumerate(tags))
    _,curves,surfaces,_,_,boundaries=_model_volume_embedding_inventory(m,volume,caller)
    boundary_set=Set(abs.(boundaries))
    outer=Tuple{Int,UInt64,NTuple{3,Int32}}[]
    embedded=Tuple{Int,UInt64,NTuple{3,Int32}}[]
    edges=Tuple{Int,UInt64,Tuple{Int32,Int32}}[]
    function actual_row(product,row)
        uid=product.entity_data.external_node_tags[row]
        target=get(rows,uid,Int32(0))
        target>0 || throw(ArgumentError(
            "$caller: actual Volume[$volume] input omits lower vertex $uid"))
        all(axis->isequal(product.coords[axis,row],geometry.surface.coords[axis,target]),1:3) ||
            throw(ArgumentError("$caller: actual lower vertex $uid changed coordinates"))
        return target
    end
    for surface in surfaces
        product=_model_mesh_product(authority,2,surface;caller)
        product isa MixedMesh && product.entity_data!==nothing || throw(ArgumentError(
            "$caller: actual Surface[$surface] lacks its generated cell product"))
        destination=surface in boundary_set ? outer : embedded
        folded_faces=_volume_boundary_faces(product)
        cursor=1
        for (block,labels) in zip(product.blocks,product.entity_data.external_element_tags)
            block.msh in (2,3) || isempty(block.nodes) || throw(ArgumentError(
                "$caller: actual Surface[$surface] requires linear Triangle or Quad cells"))
            for (column,label) in enumerate(labels)
                # Use the same Triangle/Quad diagonal policy as the PLC builder.
                count=block.msh==2 ? 1 : 2
                for index in cursor:(cursor+count-1)
                    face=folded_faces[index]
                    push!(destination,(surface,label,ntuple(i->actual_row(product,face[i]),3)))
                end
                cursor+=count
            end
        end
    end
    for curve in curves
        product=_model_mesh_product(authority,1,curve;caller)
        product===nothing && continue
        for (block,labels) in zip(product.blocks,product.entity_data.external_element_tags)
            block.msh==1 || isempty(block.nodes) || throw(ArgumentError(
                "$caller: actual Curve[$curve] requires linear Line cells"))
            for (column,label) in enumerate(labels)
                edge=ntuple(i->actual_row(product,block.nodes[i,column]),2)
                push!(edges,(curve,label,edge))
            end
        end
    end
    return (outer_faces=outer,embedded_faces=embedded,edges=edges)
end

# Certify the exact FEM incidences that generation promises to retain. A
# positive, closed, unit-volume tet mesh can still fail every published PLC cell.
function _model_certify_native_volume_constraints(mesh,constraints,volume,caller)
    mesh isa Mesh || throw(ArgumentError("$caller: actual Volume constraints require a tetrahedral mesh"))
    faces=Dict{NTuple{3,Int32},Int}()
    for tet in eachcol(mesh.tets),skip in 1:4
        face=_sorted_face3((tet[i] for i in 1:4 if i!=skip)...)
        faces[face]=get(faces,face,0)+1
    end
    for (surface,label,face) in constraints.outer_faces
        get(faces,_sorted_face3(face...),0)==1 || throw(ArgumentError(
            "$caller: Volume[$volume] lost actual boundary Surface[$surface] cell $label"))
    end
    expected=Set(_sorted_face3(face...) for (_,_,face) in constraints.outer_faces)
    observed=Set(face for (face,count) in faces if count==1)
    expected==observed || throw(ArgumentError(
        "$caller: Volume[$volume] tetrahedron boundary differs from its actual Surface cells"))
    for (surface,label,face) in constraints.embedded_faces
        get(faces,_sorted_face3(face...),0)==2 || throw(ArgumentError(
            "$caller: Volume[$volume] lost actual embedded Surface[$surface] cell $label"))
    end
    if !isempty(constraints.edges)
        required=Set(minmax(edge...) for (_,_,edge) in constraints.edges)
        found=Set{Tuple{Int32,Int32}}()
        sizehint!(found,length(required))
        for tet in eachcol(mesh.tets),first in 1:3,last in first+1:4
            edge=minmax(tet[first],tet[last])
            edge in required && push!(found,edge)
        end
        for (curve,label,edge) in constraints.edges
            minmax(edge...) in found || throw(ArgumentError(
                "$caller: Volume[$volume] lost actual Curve[$curve] Line cell $label"))
        end
    end
    return nothing
end


# ---------------------------------------------------------------------------
# Volume-boundary subdivision write-back
#
# Constrained recovery (the float Steiner fill, the exact conforming-Delaunay
# backend, gmsh-delaunay refinement) may place Steiner vertices on the PLC: a
# declared Surface cell (a,b,c) arrives on the tetrahedron boundary as the
# children (a,b,s), (b,c,s), (c,a,s) — a conforming triangulation of the same
# coplanar region, never a moved boundary. `_model_certify_native_volume_
# constraints` keeps its literal incidence audit (the same check meshGRegion
# runs against recovered MFace objects); what must change is the product the
# audit sees. Upstream recovers the boundary and writes the recovered children
# back onto each GFace/GEdge — new MVertex owners, fresh MFace/MEdge labels,
# stored curve discretizations — which is the dimension-2 analogue of
# `_model_surface_boundary_writeback!` (ModelMesh1D.jl). This routine derives
# that write-back from the recovered mesh:
#
#   * every tetrahedron boundary face is attributed to exactly one declared
#     Surface's exact coplanar regions (`RecoverCDT.build_regions` +
#     `in_region_closed` on the rational centroid — plane membership uses a
#     Float64 tolerance only where an exact-rational Steiner preimage was
#     rounded ~1 ulp off its supporting plane);
#   * every lost Curve cell is recovered as a chain of tet edges over the
#     same segment;
#   * each recovery-born boundary vertex is assigned its owning entity —
#     (1, Curve) on a realized crease chain, else the single Surface whose
#     region contains it (the smallest tag when coplanar siblings share an
#     edge, mirroring the ambiguous-vertex convention);
#   * children inherit fresh cell labels minted through the births element
#     high-water so publication keeps UID uniqueness;
#   * the receipt carries the complete subdivided cell lists, the dropped
#     source connectivity, node owners and parameters so the dim-01 assembly
#     can replace the written-back products atomically.
#
# It returns `(constraints, receipt)`; `constraints` is the literal-
# certificate-ready triple list (kept cells plus children) that
# `_model_certify_native_volume_constraints` audits unchanged. `receipt` is
# `nothing` when the recovered boundary already matched the declared cells.
# `persist=false` derives the subdivided constraints without consuming
# allocator labels or stored curve discretizations (used by the mid-pipeline
# audit before refinement can subdivide the boundary again).
# ---------------------------------------------------------------------------

# Exact plane membership first; the Float64 distance only admits the ~1 ulp
# rounding between an exact-rational Steiner preimage and its output image.
@inline function _model_volume_wb_on_plane(plane,pt,tol::Float64)
    RecoverCDT.side(plane,pt)==0 && return true
    n=plane.N;a0=plane.A0
    nx=Float64(n[1]);ny=Float64(n[2]);nz=Float64(n[3])
    norm=sqrt(nx*nx+ny*ny+nz*nz)
    norm==0.0 && return false
    dist=abs(nx*(Float64(pt[1])-Float64(a0[1]))+
             ny*(Float64(pt[2])-Float64(a0[2]))+
             nz*(Float64(pt[3])-Float64(a0[3])))/norm
    return dist<=tol
end

# 2D barycentric weights of `p` on triangle (a,b,c) projected on the dominant
# axis. `nothing` for a degenerate triangle; weights are Float64.
function _model_volume_wb_barycentric(p,a,b,c)
    ax,ay,az=a;bx,by,bz=b;cx,cy,cz=c
    nx=(by-ay)*(cz-az)-(bz-az)*(cy-ay)
    ny=(bz-az)*(cx-ax)-(bx-ax)*(cz-az)
    nz=(bx-ax)*(cy-ay)-(by-ay)*(cx-ax)
    drop=abs(nx)>=abs(ny) ? (abs(nx)>=abs(nz) ? 1 : 3) :
                            (abs(ny)>=abs(nz) ? 2 : 3)
    proj(q)=drop==1 ? (q[2],q[3]) : (drop==2 ? (q[1],q[3]) : (q[1],q[2]))
    a2,b2,c2,p2=proj(a),proj(b),proj(c),proj(p)
    d=(b2[1]-a2[1])*(c2[2]-a2[2])-(b2[2]-a2[2])*(c2[1]-a2[1])
    isfinite(d) && d!=0.0 || return nothing
    w1=((b2[1]-p2[1])*(c2[2]-p2[2])-(b2[2]-p2[2])*(c2[1]-p2[1]))/d
    w2=((c2[1]-p2[1])*(a2[2]-p2[2])-(c2[2]-p2[2])*(a2[1]-p2[1]))/d
    return (w1,w2,1.0-w1-w2)
end

# The declared Surface/Curve cell list the constraint triples came from.
# With the mesh authority this replays the producer fold used by
# `_model_native_volume_input_constraints`, carrying the product rows (for
# parameter interpolation) and physical tags. Without it (bare-constraint
# unit tests) every entry is treated as its own linear cell.
function _model_volume_wb_constraint_cells(volume,constraints,geometry,authority,caller)
    surfaces=Dict{Int,Vector{NamedTuple}}()
    curves=Dict{Int,Vector{NamedTuple}}()
    surface_order=Int[];curve_order=Int[]
    if authority===nothing
        for (surface,label,face) in constraints.outer_faces
            haskey(surfaces,surface) || push!(surface_order,surface)
            push!(get!(surfaces,surface,NamedTuple[]),
                (label=label,msh=2,conn=Int32[face...],faces=[face],
                 rows=Int32[],face_rows=NTuple{3,Int32}[],physical=Int32(0)))
        end
        for (curve,label,edge) in constraints.edges
            haskey(curves,curve) || push!(curve_order,curve)
            push!(get!(curves,curve,NamedTuple[]),
                (label=label,msh=1,conn=Int32[edge[1],edge[2]],
                 rows=Int32[],physical=Int32(0)))
        end
        return surfaces,curves,surface_order,curve_order
    end
    tags=geometry.native_node_tags
    tags!==nothing && length(tags)==nnodes(geometry.surface) ||
        throw(ArgumentError("$caller: actual Volume[$volume] lacks unique initial row bindings"))
    rows=Dict(uid=>Int32(row) for (row,uid) in enumerate(tags))
    slabels=Dict{Int,Set{UInt64}}();clabels=Dict{Int,Set{UInt64}}()
    for (surface,label,_) in constraints.outer_faces
        push!(get!(slabels,surface,Set{UInt64}()),label)
    end
    for (curve,label,_) in constraints.edges
        push!(get!(clabels,curve,Set{UInt64}()),label)
    end
    function actual_row(product,prow)
        uid=product.entity_data.external_node_tags[prow]
        target=get(rows,uid,Int32(0))
        target>0 || throw(ArgumentError(
            "$caller: actual Volume[$volume] input omits lower vertex $uid"))
        return target
    end
    for surface in sort!(collect(keys(slabels)))
        product=_model_mesh_product(authority,2,surface;caller)
        product isa MixedMesh && product.entity_data!==nothing || throw(ArgumentError(
            "$caller: actual Surface[$surface] lacks its generated cell product"))
        folded=_volume_boundary_faces(product)
        cells=NamedTuple[];cursor=1
        for (block,labels) in zip(product.blocks,product.entity_data.external_element_tags)
            block.msh in (2,3) || isempty(block.nodes) || throw(ArgumentError(
                "$caller: actual Surface[$surface] requires linear Triangle or Quad cells"))
            for (column,label) in enumerate(labels)
                count=block.msh==2 ? 1 : 2
                conn=Int32[actual_row(product,block.nodes[i,column])
                           for i in 1:size(block.nodes,1)]
                prows=Int32[block.nodes[i,column] for i in 1:size(block.nodes,1)]
                faces=NTuple{3,Int32}[];frows=NTuple{3,Int32}[]
                for index in cursor:(cursor+count-1)
                    face=folded[index]
                    push!(faces,ntuple(i->Int32(actual_row(product,face[i])),3))
                    push!(frows,ntuple(i->Int32(face[i]),3))
                end
                cursor+=count
                push!(cells,(label=label,msh=Int(block.msh),conn=conn,faces=faces,
                             rows=prows,face_rows=frows,
                             physical=Int32(block.tags[column])))
            end
        end
        surfaces[surface]=cells;push!(surface_order,surface)
    end
    for curve in sort!(collect(keys(clabels)))
        product=_model_mesh_product(authority,1,curve;caller)
        product===nothing && continue
        cells=NamedTuple[]
        for (block,labels) in zip(product.blocks,product.entity_data.external_element_tags)
            block.msh==1 || continue
            for (column,label) in enumerate(labels)
                prows=Int32[block.nodes[1,column],block.nodes[2,column]]
                conn=Int32[actual_row(product,prows[1]),actual_row(product,prows[2])]
                push!(cells,(label=label,msh=1,conn=conn,rows=prows,
                             physical=Int32(block.tags[column])))
            end
        end
        curves[curve]=cells;push!(curve_order,curve)
    end
    return surfaces,curves,surface_order,curve_order
end

# A child face's emitted winding follows its parent region's declared
# orientation: identical signed projected orient2 in exact arithmetic.
function _model_volume_wb_orient(face,region,pts)
    drop=region.plane.drop
    sign=RecoverCDT.orient2(RecoverCDT.proj2(pts[Int(face[1])],drop),
                            RecoverCDT.proj2(pts[Int(face[2])],drop),
                            RecoverCDT.proj2(pts[Int(face[3])],drop))
    sign==0 && return face
    a,b,c=region.tris[1]
    reference=RecoverCDT.orient2(RecoverCDT.proj2(pts[a],drop),
                                 RecoverCDT.proj2(pts[b],drop),
                                 RecoverCDT.proj2(pts[c],drop))
    sign==reference && return face
    return (face[1],face[3],face[2])
end

# Physical tag a child face inherits: the first dropped cell whose folded
# face contains its centroid in 2D barycentric coordinates, else the first
# dropped cell's tag (a straddling child keeps the region's physical).
function _model_volume_wb_parent_physical(mesh,face,dropped_cells)
    isempty(dropped_cells) && return Int32(0)
    center=((mesh.coords[1,face[1]]+mesh.coords[1,face[2]]+mesh.coords[1,face[3]])/3,
            (mesh.coords[2,face[1]]+mesh.coords[2,face[2]]+mesh.coords[2,face[3]])/3,
            (mesh.coords[3,face[1]]+mesh.coords[3,face[2]]+mesh.coords[3,face[3]])/3)
    for cell in dropped_cells
        for folded in cell.faces
            w=_model_volume_wb_barycentric(center,_pt3(mesh,folded[1]),
                _pt3(mesh,folded[2]),_pt3(mesh,folded[3]))
            w===nothing && continue
            all(x->x>=-1e-9,w) && return cell.physical
        end
    end
    return dropped_cells[1].physical
end

# Parameter for a recovery-born vertex: the declared parametrization first,
# then barycentric interpolation inside a containing declared cell (Surface)
# or along the lost cell's segment (Curve) when the authority supplies
# product parameters.
function _model_volume_wb_param(m,volume,owner,row,mesh,authority,
                                surface_cells,curve_cells,caller)
    owner[1] in (1,2) || return nothing
    point=collect(_pt3(mesh,row))
    value=try
        model_parametrization(m,owner[1],Int(owner[2]),point)
    catch err
        err isa InterruptException && rethrow()
        nothing
    end
    value!==nothing && all(isfinite,value) && return Float64.(value)
    authority===nothing && return nothing
    if owner[1]==2
        surface=Int(owner[2])
        product=_model_mesh_product(authority,2,surface;caller)
        product isa MixedMesh && product.entity_data!==nothing || return nothing
        params=product.entity_data.node_parametric
        params===nothing && return nothing
        for cell in get(surface_cells,surface,NamedTuple[])
            isempty(cell.rows) && continue
            for (fi,face) in enumerate(cell.faces)
                w=_model_volume_wb_barycentric(
                    _pt3(mesh,row),_pt3(mesh,face[1]),_pt3(mesh,face[2]),
                    _pt3(mesh,face[3]))
                w===nothing && continue
                all(x->x>=-1e-9,w) || continue
                prows=cell.face_rows[fi]
                p1=params[prows[1]];p2=params[prows[2]];p3=params[prows[3]]
                (p1===nothing || p2===nothing || p3===nothing) && continue
                interpolated=[w[1]*p1[k]+w[2]*p2[k]+w[3]*p3[k]
                              for k in eachindex(p1)]
                all(isfinite,interpolated) || continue
                return interpolated
            end
        end
        return nothing
    end
    curve=Int(owner[2])
    product=_model_mesh_product(authority,1,curve;caller)
    product isa MixedMesh && product.entity_data!==nothing || return nothing
    params=product.entity_data.node_parametric
    params===nothing && return nothing
    p=_pt3(mesh,row)
    for cell in get(curve_cells,curve,NamedTuple[])
        isempty(cell.rows) && continue
        u,v=cell.conn
        pu,pv=_pt3(mesh,u),_pt3(mesh,v)
        seg=(pv[1]-pu[1],pv[2]-pu[2],pv[3]-pu[3])
        len2=seg[1]*seg[1]+seg[2]*seg[2]+seg[3]*seg[3]
        len2==0.0 && continue
        alpha=((p[1]-pu[1])*seg[1]+(p[2]-pu[2])*seg[2]+(p[3]-pu[3])*seg[3])/len2
        (alpha<-1e-9 || alpha>1+1e-9) && continue
        tu=params[cell.rows[1]];tv=params[cell.rows[2]]
        (tu===nothing || tv===nothing || isempty(tu) || isempty(tv)) && continue
        interpolated=[tu[k]+alpha*(tv[k]-tu[k]) for k in eachindex(tu)]
        all(isfinite,interpolated) || continue
        return interpolated
    end
    return nothing
end

function _model_native_volume_boundary_writeback!(m,volume,mesh::Mesh,constraints,
                                                 geometry,authority,births,
                                                 caller;persist::Bool=true)
    face_count=Dict{NTuple{3,Int32},Int}()
    for tet in eachcol(mesh.tets),skip in 1:4
        face=_sorted_face3((tet[i] for i in 1:4 if i!=skip)...)
        face_count[face]=get(face_count,face,0)+1
    end
    boundary=Set(face for (face,count) in face_count if count==1)
    declared=Set(_sorted_face3(face...) for (_,_,face) in constraints.outer_faces)
    declared==boundary && return constraints,nothing
    input_rows=nnodes(geometry.surface)
    RB=Rational{BigInt}
    pts=RecoverCDT.RBv[(RB(mesh.coords[1,row]),RB(mesh.coords[2,row]),
                        RB(mesh.coords[3,row])) for row in 1:nnodes(mesh)]
    scale=isempty(mesh.coords) ? 1.0 : max(1.0,maximum(abs,mesh.coords))
    plane_tol=1e-9*scale
    budget=RecoverCDT._RecoveryBudget(typemax(Int)÷4,typemax(Int)÷4,String(caller))
    surface_cells,curve_cells,surface_order,curve_order=
        _model_volume_wb_constraint_cells(volume,constraints,geometry,authority,caller)
    # Kept cells contribute their folded faces literally; dropped cells leave
    # their region to the recovered children.
    kept_faces=Set{NTuple{3,Int32}}()
    kept_cell=Dict{Tuple{Int,Int},Bool}()     # (surface, cell index) => kept
    for (surface,cells) in surface_cells
        for (index,cell) in enumerate(cells)
            kept=all(face->_sorted_face3(face...) in boundary,cell.faces)
            kept_cell[(surface,index)]=kept
            kept && foreach(face->push!(kept_faces,_sorted_face3(face...)),cell.faces)
        end
    end
    extra=setdiff(boundary,kept_faces)
    isempty(extra) && return constraints,nothing
    # Exact coplanar regions per declared surface (folded-facet input).
    regions=Dict{Int,Vector{RecoverCDT.Region}}()
    for surface in surface_order
        facets=NTuple{3,Int}[Int.(face) for cell in surface_cells[surface]
                             for face in cell.faces]
        regions[surface]=RecoverCDT.build_regions(pts,facets)
    end
    face_owner=Dict{NTuple{3,Int32},Int}()
    face_region=Dict{NTuple{3,Int32},RecoverCDT.Region}()
    for face in sort!(collect(extra))
        center=((pts[Int(face[1])][1]+pts[Int(face[2])][1]+pts[Int(face[3])][1])//3,
                (pts[Int(face[1])][2]+pts[Int(face[2])][2]+pts[Int(face[3])][2])//3,
                (pts[Int(face[1])][3]+pts[Int(face[2])][3]+pts[Int(face[3])][3])//3)
        found=0
        for surface in surface_order
            for region in regions[surface]
                all(v->_model_volume_wb_on_plane(region.plane,pts[Int(v)],plane_tol),
                    face) || continue
                RecoverCDT.in_region_closed(region,pts,center,budget) || continue
                if found!=0 && found!=surface
                    throw(ArgumentError("$caller: Volume[$volume] recovered " *
                        "boundary face is claimed by both Surface[$found] and " *
                        "Surface[$surface]"))
                end
                found=surface;face_region[face]=region
            end
        end
        found==0 && throw(ArgumentError(
            "$caller: Volume[$volume] tetrahedron boundary face $face is not " *
            "covered by any actual Surface region"))
        face_owner[face]=found
    end
    # Lost Curve cells recover as tet-edge chains over their declared segment.
    tet_edges=_tet_edge_set(mesh)
    chain_children=Dict{Tuple{Int,UInt64},Vector{Tuple{Int32,Int32}}}()
    chain_nodes=Dict{Int32,Int}()
    for curve in curve_order
        for cell in curve_cells[curve]
            u,v=cell.conn
            minmax(u,v) in tet_edges && continue
            p,q=_pt3(mesh,u),_pt3(mesh,v)
            seg=(q[1]-p[1],q[2]-p[2],q[3]-p[3])
            length2=seg[1]*seg[1]+seg[2]*seg[2]+seg[3]*seg[3]
            chain=Int32[u]
            for row in 1:nnodes(mesh)
                (row==u || row==v) && continue
                _on_segment3(_pt3(mesh,row),p,q;atol=1e-9) || continue
                push!(chain,Int32(row))
            end
            push!(chain,v)
            sort!(chain;by=row->(mesh.coords[1,row]-p[1])*seg[1]+
                                (mesh.coords[2,row]-p[2])*seg[2]+
                                (mesh.coords[3,row]-p[3])*seg[3])
            for i in 1:(length(chain)-1)
                minmax(chain[i],chain[i+1]) in tet_edges || throw(ArgumentError(
                    "$caller: Volume[$volume] Curve[$curve] Line cell " *
                    "$(cell.label) has no covering tetrahedron-edge chain"))
            end
            chain_children[(curve,cell.label)]=
                [(chain[i],chain[i+1]) for i in 1:(length(chain)-1)]
            for i in 2:(length(chain)-1)
                row=chain[i]
                owner=get(chain_nodes,row,0)
                (owner==0 || owner==curve) || throw(ArgumentError(
                    "$caller: Volume[$volume] chain vertex lies on both " *
                    "Curve[$owner] and Curve[$curve]"))
                chain_nodes[row]=curve
            end
        end
    end
    # Recovery-born boundary vertex ownership: crease chains first, then the
    # unique containing Surface (smallest tag on a shared coplanar edge).
    nodes=Dict{Int32,Tuple{Int,Int32}}()
    for (row,curve) in chain_nodes
        row>input_rows && (nodes[row]=(1,Int32(curve)))
    end
    for (face,surface) in face_owner
        for row in face
            (row>input_rows && !haskey(nodes,row)) || continue
            candidates=Int[]
            for other in surface_order
                for region in regions[other]
                    _model_volume_wb_on_plane(region.plane,pts[Int(row)],plane_tol) &&
                    RecoverCDT.in_region_closed(region,pts,pts[Int(row)],budget) &&
                        push!(candidates,other)
                end
            end
            nodes[row]=(2,Int32(isempty(candidates) ? surface : minimum(candidates)))
        end
    end
    node_params=Dict{Int32,Vector{Float64}}()
    m===nothing || for (row,owner) in nodes
        param=_model_volume_wb_param(m,volume,owner,row,mesh,authority,
                                     surface_cells,curve_cells,caller)
        param===nothing || (node_params[row]=param)
    end
    # Mint child labels through the births high-water (or a throwaway counter
    # for the mid-pipeline audit), then build the new constraint triples.
    next_label=Ref(births===nothing ? UInt64(0) : births.element_high_water)
    mint!()=(next_label[]+=UInt64(1);next_label[])
    children=Dict{Int,Vector{NTuple{3,Int32}}}()
    for (face,surface) in face_owner
        oriented=_model_volume_wb_orient(face,face_region[face],pts)
        push!(get!(children,surface,NTuple{3,Int32}[]),oriented)
    end
    new_outer=Tuple{Int,UInt64,NTuple{3,Int32}}[]
    surfaces_out=Dict{Int,Vector{Tuple{Int,Vector{Int32},UInt64}}}()
    replaced=Dict{Tuple{Int,Int32},Set{Tuple}}()
    added=Dict{Tuple{Int,Int32},Vector{Tuple{Int,Vector{Int32},UInt64,Int32}}}()
    for surface in surface_order
        cells=surface_cells[surface];out=Tuple{Int,Vector{Int32},UInt64}[]
        dropped_cells=NamedTuple[]
        for (index,cell) in enumerate(cells)
            if kept_cell[(surface,index)]
                push!(out,(cell.msh,cell.conn,cell.label))
                for face in cell.faces
                    push!(new_outer,(surface,cell.label,face))
                end
            else
                push!(dropped_cells,cell)
                push!(get!(replaced,(2,Int32(surface)),Set{Tuple}()),
                      Tuple(sort!(copy(cell.conn))))
            end
        end
        for face in sort!(collect(get(children,surface,NTuple{3,Int32}[])))
            label=mint!()
            physical=_model_volume_wb_parent_physical(mesh,face,dropped_cells)
            push!(out,(2,Int32[face[1],face[2],face[3]],label))
            push!(new_outer,(surface,label,face))
            push!(get!(added,(2,Int32(surface)),
                       Tuple{Int,Vector{Int32},UInt64,Int32}[]),
                  (2,Int32[face[1],face[2],face[3]],label,physical))
        end
        surfaces_out[surface]=out
    end
    new_edges=Tuple{Int,UInt64,Tuple{Int32,Int32}}[]
    curves_out=Dict{Int,Vector{Tuple{Int,Vector{Int32},UInt64}}}()
    for curve in curve_order
        out=Tuple{Int,Vector{Int32},UInt64}[]
        for cell in curve_cells[curve]
            u,v=cell.conn
            if minmax(u,v) in tet_edges
                push!(out,(1,Int32[u,v],cell.label))
                push!(new_edges,(curve,cell.label,(u,v)))
            else
                push!(get!(replaced,(1,Int32(curve)),Set{Tuple}()),
                      Tuple(sort!(Int32[u,v])))
                for (a,b) in chain_children[(curve,cell.label)]
                    label=mint!()
                    push!(out,(1,Int32[a,b],label))
                    push!(new_edges,(curve,label,(a,b)))
                    push!(get!(added,(1,Int32(curve)),
                           Tuple{Int,Vector{Int32},UInt64,Int32}[]),
                          (1,Int32[a,b],label,cell.physical))
                end
            end
        end
        curves_out[curve]=out
    end
    # A curve without an actual product keeps its declared edges verbatim —
    # the literal certificate then audits them honestly. Curves already
    # emitted (kept or chain-subdivided) are not re-declared: a split cell's
    # source label must not resurrect its unrecovered chord.
    processed=Set(curve_order)
    for (curve,label,edge) in constraints.edges
        curve in processed || push!(new_edges,(curve,label,edge))
    end
    new_constraints=(outer_faces=new_outer,embedded_faces=constraints.embedded_faces,
                     edges=new_edges)
    receipt=_ModelVolumeBoundaryWriteback(new_constraints,nodes,node_params,
                                          surfaces_out,curves_out,replaced,added)
    persist || return new_constraints,receipt
    births!==nothing && (births.element_high_water=next_label[])
    births!==nothing && (births.boundary_writeback=receipt)
    # Upstream `meshGRegion` stores the recovered chain discretization back on
    # each GEdge (`GEdge::mesh_vertices`); mirror that on `m.curve_params` so a
    # later Curve produce reuses the realized chain.
    if m!==nothing
        for curve in curve_order
            stored=get(m.curve_params,curve,nothing)
            stored===nothing && continue
            additions=Float64[]
            for cell in curve_cells[curve]
                u,v=cell.conn
                minmax(u,v) in tet_edges && continue
                haskey(chain_children,(curve,cell.label)) || continue
                for (a,b) in chain_children[(curve,cell.label)]
                    for row in (a,b)
                        row>input_rows || continue
                        param=get(node_params,row,nothing)
                        if param===nothing && m!==nothing
                            param=_model_volume_wb_param(m,volume,(1,Int32(curve)),
                                row,mesh,authority,surface_cells,curve_cells,caller)
                        end
                        param===nothing || isempty(param) || push!(additions,param[1])
                    end
                end
            end
            isempty(additions) && continue
            tolerance=1e-12*max(1.0,maximum(abs,stored;init=0.0))
            ordered=issorted(stored) && all(isfinite,stored)
            for value in additions
                isfinite(value) || continue
                _model_surface_curve_parameter_match(stored,value,tolerance,
                                                     ordered)>0 && continue
                if ordered
                    index=searchsortedfirst(stored,value)
                    insert!(stored,index,value)
                else
                    push!(stored,value)
                end
            end
        end
    end
    return new_constraints,receipt
end


# The first owned Point vertex supplies Curve endpoints. Additional Point
# records keep their separate identities without becoming endpoint aliases.
# Bind that first source by UID; coordinates verify a bound row only.
function _model_native_volume_point_bindings(authority,geometry,point_tags,caller)
    native_rows=Dict(uid=>Int32(row) for (row,uid) in enumerate(geometry.native_node_tags))
    bindings=Dict{Int,Tuple{UInt64,Int32,NTuple{3,Float64}}}()
    actual=authority.mesh;data=actual.entity_data
    for point in point_tags
        owner=(0,Int32(point))
        product=get(authority.products,owner,nothing)
        product!==nothing && !isempty(product.rows) ||
            throw(ArgumentError("$caller: actual Point[$point] lacks its owned source vertex"))
        source_rows=get(authority.producer_rows,owner,product.rows)
        isempty(source_rows) && throw(ArgumentError("$caller: actual Point[$point] lacks its first source row"))
        source=first(source_rows)
        uid=data.external_node_tags[source]
        row=get(native_rows,uid,Int32(0))
        row>0 || throw(ArgumentError("$caller: actual Point[$point] UID $uid is absent from native input"))
        coordinate=_model_mesh_coordinate(actual,source)
        isequal(coordinate,_model_mesh_coordinate(geometry.surface,row)) ||
            throw(ArgumentError("$caller: actual Point[$point] UID $uid changed input coordinates"))
        bindings[point]=(uid,row,coordinate)
    end
    return bindings
end

function _model_native_volume_point_nodes(mesh,bindings,births,caller)
    # One UID index per stage avoids scanning every mesh row for each Point.
    native_rows=births!==nothing && length(births.node_tags)==nnodes(mesh) ?
        Dict(uid=>row for (row,uid) in enumerate(births.node_tags)) : nothing
    nodes=Dict{Int,Int}()
    for (point,(uid,initial_row,coordinate)) in bindings
        row=native_rows===nothing ? Int(initial_row) : get(native_rows,uid,0)
        row>0 || throw(ArgumentError("$caller: actual Point[$point] UID $uid lost its native row"))
        # Without native storage ordering, refinement retains the input prefix.
        # Never recover an identity by geometric proximity or CAD coincidence.
        row<=nnodes(mesh) && isequal(coordinate,_model_mesh_coordinate(mesh,row)) ||
            throw(ArgumentError("$caller: actual Point[$point] UID $uid changed its bound coordinates"))
        nodes[point]=row
    end
    return nodes
end
