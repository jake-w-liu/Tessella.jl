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
end
_ModelNativeVolumeBirths(node_max::UInt64=UInt64(0),element_max::UInt64=UInt64(0))=
    _ModelNativeVolumeBirths(nothing,UInt64[],node_max,0,UInt64[],element_max,nothing)

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

struct _ModelVolumeGeometryProduction
    mesh::Union{Mesh,MixedMesh}
    boundary_node_tags::Vector{UInt64}
    node_high_water::UInt64
    temporary_elements::UInt64
    element_tags::Vector{UInt64}
    element_high_water::UInt64
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
        births.temporary_elements,cells,births.element_high_water)
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
