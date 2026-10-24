# A transfinite fill consumes already-emitted Surface vectors. Boundary rows
# retain their public identity; only the higher assembly allocates interiors.
# The native initial-Delaunay birth ledger remains a separate producer.
# A transfinite quadrilateral grid has a logical cell before either Surface
# or Volume chooses a diagonal. Certify each complete cell partition through
# its UID incidence, independently of coordinates and diagonal selection.
include("ModelTransfiniteGridAlias.jl")

function _model_transfinite_grid_partition(panels,cells,caller)
    possible=Dict{Tuple,Int}()
    supports=Vector{Vector{UInt64}}()
    for (index,panel) in enumerate(panels)
        vertices=unique(UInt64.(collect(panel)))
        length(vertices) in (3,4) || throw(ArgumentError(
            "$caller: transfinite grid cell has no triangle or quadrilateral support"))
        push!(supports,vertices)
        keys=Tuple[Tuple(sort(vertices))]
        if length(vertices)==4
            for missing in 1:4
                push!(keys,Tuple(sort(UInt64[vertices[i] for i in 1:4 if i!=missing])))
            end
        end
        for key in keys
            haskey(possible,key) && throw(ArgumentError(
                "$caller: transfinite grid cells have ambiguous UID incidence"))
            possible[key]=index
        end
    end
    grouped=[Tuple[] for _ in panels]
    mapping=Dict{Tuple,Int}()
    for cell in cells
        face=Tuple(sort(UInt64.(collect(cell))))
        length(face) in (3,4) && allunique(face) || throw(ArgumentError(
            "$caller: transfinite boundary cell has invalid UID incidence"))
        index=get(possible,face,0)
        index>0 || throw(ArgumentError(
            "$caller: transfinite boundary cell crosses its recorded grid cells"))
        haskey(mapping,face) && throw(ArgumentError(
            "$caller: transfinite boundary cell repeats UID incidence"))
        mapping[face]=index;push!(grouped[index],face)
    end
    for (vertices,faces) in zip(supports,grouped)
        whole=Tuple(sort(vertices))
        length(faces)==1 && only(faces)==whole && continue
        length(vertices)==4 && length(faces)==2 && all(face->length(face)==3,faces) ||
            throw(ArgumentError("$caller: transfinite grid cell lacks complete boundary coverage"))
        first,second=Set(faces[1]),Set(faces[2])
        shared=intersect(first,second)
        union(first,second)==Set(vertices) && length(shared)==2 &&
            (shared==Set((vertices[1],vertices[3])) || shared==Set((vertices[2],vertices[4]))) ||
            throw(ArgumentError("$caller: transfinite grid triangle pair has an invalid diagonal"))
    end
    return mapping
end

# A compact triangle producer can join atoms on adjacent logical panels.
# Retain its original triangles, rather than assuming every emitted Quad is
# contained in one panel. The ordinary panel certificate still certifies the
# retained atoms themselves and the final Volume exterior.
function _model_transfinite_grid_compact_atoms(ids,triangles,panels,caller)
    triangles isa AbstractMatrix && size(triangles,1)==3 &&
        eltype(triangles)<:Integer && !(eltype(triangles)<:Bool) ||
        throw(ArgumentError("$caller: compact Surface lacks its triangle constructor rows"))
    atoms=NTuple{3,UInt64}[]
    sizehint!(atoms,size(triangles,2))
    for column in axes(triangles,2)
        rows=triangles[:,column]
        all(row->1<=row<=length(ids),rows) || throw(ArgumentError(
            "$caller: compact triangle constructor references an absent grid row"))
        push!(atoms,ntuple(slot->UInt64(ids[rows[slot]]),3))
    end
    _model_transfinite_grid_partition(panels,atoms,caller)
    return atoms
end

function _model_transfinite_grid_edge_chain(cells)
    chain=Dict{Tuple{UInt64,UInt64},Int}()
    for cell in cells,slot in eachindex(cell)
        first=UInt64(cell[slot]);second=UInt64(cell[mod1(slot+1,length(cell))])
        edge=minmax(first,second)
        value=get(chain,edge,0)+(first<second ? 1 : -1)
        value==0 ? delete!(chain,edge) : (chain[edge]=value)
    end
    return chain
end

function _model_transfinite_grid_atomic_partition(atoms,cells,caller)
    lookup=Dict{Tuple,Int}()
    incidence=Dict{Tuple{UInt64,UInt64},Tuple{Int,Int}}()
    for (index,atom) in enumerate(atoms)
        length(atom)==3 && allunique(atom) || throw(ArgumentError(
            "$caller: compact constructor atom has invalid UID incidence"))
        key=Tuple(sort(UInt64.(collect(atom))))
        haskey(lookup,key) && throw(ArgumentError(
            "$caller: compact constructor repeats an atomic triangle"))
        lookup[key]=index
        for slot in 1:3
            first=UInt64(atom[slot]);second=UInt64(atom[mod1(slot+1,3)])
            edge=minmax(first,second)
            count,direction=get(incidence,edge,(0,0))
            count+=1;direction+=first<second ? 1 : -1
            count<=2 && (count==1 || direction==0) || throw(ArgumentError(
                "$caller: compact constructor atoms are not an oriented conforming sheet"))
            incidence[edge]=(count,direction)
        end
    end
    consumed=falses(length(atoms))
    orientation=0
    for cell in cells
        count=length(cell)
        count in (3,4) && allunique(cell) || throw(ArgumentError(
            "$caller: compact Surface cell has invalid UID incidence"))
        covered=Int[]
        if count==3
            key=Tuple(sort(UInt64.(collect(cell))))
            index=get(lookup,key,0)
            index>0 && push!(covered,index)
        else
            for missing in 1:4
                key=Tuple(sort(UInt64[cell[slot] for slot in 1:4 if slot!=missing]))
                index=get(lookup,key,0)
                index>0 && push!(covered,index)
            end
        end
        length(covered)==count-2 || throw(ArgumentError(
            "$caller: compact Surface cell crosses its recorded atomic triangles"))
        # An oriented boundary comparison requires the two triangles to
        # share the actual Quad diagonal. It rejects bow-ties and mismatched
        # cyclic corner orders without consulting mesh coordinates.
        atom_chain=_model_transfinite_grid_edge_chain(atoms[covered])
        cell_chain=_model_transfinite_grid_edge_chain((cell,))
        current=atom_chain==cell_chain ? 1 :
            (length(atom_chain)==length(cell_chain) &&
             all(get(cell_chain,edge,0)==-value for (edge,value) in atom_chain) ? -1 : 0)
        current!=0 || throw(ArgumentError(
            "$caller: compact Surface cell changes its constructor edge orientation"))
        orientation==0 && (orientation=current)
        current==orientation || throw(ArgumentError(
            "$caller: compact Surface cells have inconsistent constructor orientation"))
        for index in covered
            !consumed[index] || throw(ArgumentError(
                "$caller: compact Surface cells overlap a recorded atomic triangle"))
            consumed[index]=true
        end
    end
    all(consumed) || throw(ArgumentError(
        "$caller: compact Surface cells leave recorded atomic triangles uncovered"))
    return nothing
end

function _model_transfinite_grid_panels(ids,extents,caller)
    width,height=extents
    width>0 && height>0 && length(ids)==(width+1)*(height+1) || throw(ArgumentError(
        "$caller: transfinite boundary grid extent does not align with its UID rows"))
    panels=Vector{UInt64}[]
    stride=width+1
    for j in 0:height-1,i in 0:width-1
        first=i+1+j*stride
        panel=UInt64[ids[first],ids[first+1],ids[first+1+stride],ids[first+stride]]
        length(unique(panel))>=3 || continue # Collapsed/compact aliased slots are not cells.
        push!(panels,panel)
    end
    return panels
end

# Keep the actual input and its ordered Surface products with the logical grids.
# The grid cells precede the independent Surface and Region diagonal choices.
struct _ModelTransfiniteGridFace
    owner::Int32
    ids::Vector{UInt64}
    extents::NTuple{2,Int}
    cells::Vector{Tuple{Int,Tuple}}
    element_tags::Vector{Vector{UInt64}}
    atomic_cells::Union{Nothing,Vector{NTuple{3,UInt64}}}
end

_ModelTransfiniteGridFace(owner,ids,extents,cells,element_tags)=
    _ModelTransfiniteGridFace(owner,ids,extents,cells,element_tags,nothing)

struct _ModelTransfiniteGridReceipt
    input::_ModelMeshAuthority
    faces::Vector{_ModelTransfiniteGridFace}
    aliases::Union{Nothing,_ModelTransfiniteGridAliasReceipt}
end

_ModelTransfiniteGridReceipt(input,faces)=_ModelTransfiniteGridReceipt(input,faces,nothing)

function _model_transfinite_grid_source_cells(input,product,caller)
    data=input.mesh.entity_data
    cells=Tuple{Int,Tuple}[]
    length(product.blocks)==length(product.element_tags) || throw(ArgumentError(
        "$caller: transfinite Surface product lacks its element bindings"))
    for (block,tags) in zip(product.blocks,product.element_tags)
        isempty(block.nodes) && continue
        block.msh in (2,3) && length(tags)==size(block.nodes,2) || throw(ArgumentError(
            "$caller: transfinite grid requires bound linear Surface cells"))
        for column in axes(block.nodes,2)
            ids=Tuple(data.external_node_tags[product.rows[block.nodes[row,column]]]
                      for row in axes(block.nodes,1))
            allunique(ids) || throw(ArgumentError(
                "$caller: transfinite Surface cell repeats a source UID"))
            push!(cells,(Int(block.msh),ids))
        end
    end
    isempty(cells) && throw(ArgumentError("$caller: transfinite Surface grid has no cells"))
    return cells
end

function _model_transfinite_grid_receipt(input,face_ids,face_tags,counts,caller;
                                       mesh=nothing,compact=false,compact_triangles=nothing)
    input isa _ModelMeshAuthority && input.mesh!==nothing && input.mesh.entity_data!==nothing ||
        throw(ArgumentError("$caller: transfinite grid lacks actual input authority"))
    nu,nv,nw=counts
    extents=((nu,nw),(nv,nw),(nu,nw),(nv,nw),(nu,nv),(nu,nv))
    length(face_ids)==length(face_tags)==6 || throw(ArgumentError(
        "$caller: transfinite grid lacks its six canonical slots"))
    faces=_ModelTransfiniteGridFace[]
    owners=Set{Int32}()
    for slot in 1:6
        ids=face_ids[slot]
        ids===nothing && continue
        owner=Int32(face_tags[slot])
        owner>0 && !(owner in owners) || throw(ArgumentError(
            "$caller: transfinite grid repeats or omits a Surface owner"))
        push!(owners,owner)
        product=get(input.products,(2,owner),nothing)
        product!==nothing && product.owner==(2,owner) || throw(ArgumentError(
            "$caller: transfinite grid lacks its owned Surface product"))
        actual_ids=input.mesh.entity_data.external_node_tags[product.rows]
        Set(ids)==Set(actual_ids) && allunique(actual_ids) || throw(ArgumentError(
            "$caller: transfinite grid does not cover its actual Surface rows"))
        panels=_model_transfinite_grid_panels(ids,extents[slot],caller)
        cells=_model_transfinite_grid_source_cells(input,product,caller)
        atoms=if compact && slot in (5,6)
            if compact_triangles===nothing
                # Existing unrecombined receipts already carry their actual
                # triangles. A recombined receipt needs the producer trace.
                all(cell->first(cell)==2,cells) || throw(ArgumentError(
                    "$caller: recombined compact Surface lacks its atomic constructor trace"))
                captured=NTuple{3,UInt64}[Tuple(UInt64.(last(cell))) for cell in cells]
                _model_transfinite_grid_partition(panels,captured,caller)
                captured
            else
                length(compact_triangles)==6 || throw(ArgumentError(
                    "$caller: compact triangle constructor slots are incomplete"))
                _model_transfinite_grid_compact_atoms(ids,compact_triangles[slot],panels,caller)
            end
        else
            nothing
        end
        push!(faces,_ModelTransfiniteGridFace(owner,copy(ids),extents[slot],cells,
                                             deepcopy(product.element_tags),atoms))
    end
    length(faces) in (5,6) || throw(ArgumentError(
        "$caller: transfinite grid does not cover its Volume slots"))
    aliases=if compact
        mesh!==nothing && face_ids[4]===nothing && length(faces)==5 || throw(ArgumentError(
            "$caller: compact diagonal aliases require their recorded five-face constructor"))
        _model_transfinite_compact_alias_events(mesh,face_ids[3],counts,face_tags[3],caller)
    else
        nothing
    end
    return _ModelTransfiniteGridReceipt(input,faces,aliases)
end

function _model_transfinite_grid_boundary(receipt,input,source_tags,actual_faces,caller;mesh=nothing)
    receipt isa _ModelTransfiniteGridReceipt && receipt.input===input || throw(ArgumentError(
        "$caller: transfinite boundary lacks its current grid receipt"))
    panels=Vector{UInt64}[]
    owners=Int32[]
    surface_cells=Tuple[]
    for face in receipt.faces
        product=get(input.products,(2,face.owner),nothing)
        product!==nothing && product.owner==(2,face.owner) || throw(ArgumentError(
            "$caller: transfinite grid Surface owner changed after emission"))
        cells=_model_transfinite_grid_source_cells(input,product,caller)
        isequal(cells,face.cells) && isequal(product.element_tags,face.element_tags) &&
            Set(input.mesh.entity_data.external_node_tags[product.rows])==Set(face.ids) ||
            throw(ArgumentError("$caller: transfinite Surface product changed after grid emission"))
        local_panels=_model_transfinite_grid_panels(face.ids,face.extents,caller)
        evidence=if face.atomic_cells===nothing
            _model_transfinite_grid_partition(local_panels,last.(cells),caller)
            last.(cells)
        else
            receipt.aliases!==nothing || throw(ArgumentError(
                "$caller: atomic Surface coverage lacks its compact grid constructor"))
            _model_transfinite_grid_partition(local_panels,face.atomic_cells,caller)
            _model_transfinite_grid_atomic_partition(face.atomic_cells,last.(cells),caller)
            face.atomic_cells
        end
        append!(panels,local_panels)
        append!(owners,fill(face.owner,length(local_panels)))
        append!(surface_cells,evidence)
    end
    # Global incidence also excludes overlapping cells on distinct owners.
    _model_transfinite_grid_partition(panels,surface_cells,caller)
    aliases=_model_transfinite_compact_alias_binding(receipt.aliases,receipt.faces,mesh,source_tags,caller)
    volume_cells=Tuple[]
    row_faces=Tuple[]
    for face in actual_faces
        ids=Tuple(iszero(source_tags[row]) ? get(aliases,row,UInt64(0)) : source_tags[row] for row in face)
        all(!iszero,ids) || throw(ArgumentError(
            "$caller: transfinite exterior face contains an unbound interior event"))
        push!(volume_cells,ids);push!(row_faces,face)
    end
    mapping=_model_transfinite_grid_partition(panels,volume_cells,caller)
    result=Dict{Tuple,Tuple{Int,Int32}}()
    for (face,ids) in zip(row_faces,volume_cells)
        result[face]=(2,owners[mapping[Tuple(sort(collect(ids)))]])
    end
    return result
end

struct _ModelTransfiniteVolumeEmission
    block::Int
    nodes::Matrix{Int32}
    columns::Vector{Int}
    rows::Vector{Int32}
    grid::Union{Nothing,_ModelTransfiniteGridReceipt}
end

_ModelTransfiniteVolumeEmission(block,nodes,columns,rows)=
    _ModelTransfiniteVolumeEmission(block,nodes,columns,rows,nothing)

struct _ModelTransfiniteVolumeProduction
    mesh::Union{Mesh,MixedMesh}
    boundary_node_tags::Vector{UInt64}
    boundary_node_retained::BitVector
    geometry::Union{Nothing,Mesh,MixedMesh}
    emission::_ModelTransfiniteVolumeEmission
    node_high_water::UInt64
    temporary_elements::UInt64
    element_tags::Vector{UInt64}
    element_high_water::UInt64
end

function _model_transfinite_face_source(input,surf,count,caller)
    input===nothing && return nothing
    haskey(input.producer_rows,(2,Int32(surf))) || throw(ArgumentError(
        "$caller: transfinite Surface[$surf] lacks its emitted grid row binding"))
    source=_model_mesh_product(input,2,surf;caller,producer_order=true)
    source!==nothing && nnodes(source)==count || throw(ArgumentError(
        "$caller: transfinite Surface[$surf] actual grid differs from its declared layout"))
    all(block->msh_spec(block.msh).order==1,source.blocks) || throw(ArgumentError(
        "$caller: transfinite Surface[$surf] boundary grid must be first order"))
    return source
end

function _model_transfinite_source_point(source,row,geometry_points)
    tag=source.entity_data.external_node_tags[row]
    actual=ntuple(axis->source.coords[axis,row],3)
    return geometry_points===nothing ? actual : get(geometry_points,tag,actual)
end

function _model_transfinite_face_slot(ordered,junctions,caller)
    found=0;slot=0
    for p in 1:8,f in 1:6
        perm=_VOLUME_FACE_PERMS[p];corners=_VOLUME_FACE_SLOTS[f]
        all(k->ordered[corners[k]]==junctions[perm[k]],1:4) || continue
        found==0 || throw(ArgumentError("$caller: actual face grid has ambiguous canonical orientation"))
        found=p;slot=f
    end
    found>0 || throw(ArgumentError("$caller: actual face grid lacks a canonical volume slot"))
    return found,slot
end

# The copied Surface keeps the actual master's lattice frame. Declared slave
# loop order can differ by a rotation or reflection; Point counterparts give
# the frame relation without inspecting or matching mesh coordinates.
function _model_transfinite_quad_producer_frame(m,surf,input,caller,visited=Set{Int}())
    surf in visited && throw(ArgumentError("$caller: periodic Surface grid has a cyclic master chain"))
    push!(visited,surf)
    spec,signed_curves,curve_points,nside,_,junctions,_,_,_=
        _transfinite_surface_sides(m,surf,caller;_mesh_input=input)
    nside==4 || throw(ArgumentError("$caller: actual transfinite Surface[$surf] is not a four-sided grid"))
    constraint=get(m.periodic,(2,surf),nothing)
    if constraint===nothing
        haskey(input.producer_rows,(2,Int32(surf))) || throw(ArgumentError(
            "$caller: periodic master Surface[$surf] lacks its actual grid row binding"))
        corners,sides=isempty(spec.corners) ?
            _transfinite_chained_tri(m,signed_curves,curve_points,surf,caller) :
            (junctions,curve_points)
        return corners,length(sides[1])-1,length(sides[2])-1
    end
    master=Int(constraint.master_entity)
    corners,L,H=_model_transfinite_quad_producer_frame(m,master,input,caller,visited)
    point_map,_,_=_model_periodic_surface_edge_map(
        m,surf,master,constraint.affine,constraint.atol,caller)
    all(corner->haskey(point_map,corner),corners) || throw(ArgumentError(
        "$caller: periodic Surface[$surf] grid loses a master corner counterpart"))
    mapped=Int[point_map[corner] for corner in corners]
    length(mapped)==4 && allunique(mapped) && Set(mapped)==Set(junctions) ||
        throw(ArgumentError("$caller: periodic Surface[$surf] grid changes its four corner identities"))
    return mapped,L,H
end

function _model_transfinite_quad_face_source!(points,m,surf,ordered,spec,
    signed_curves,curve_points,junctions,input,geometry_points,slot,caller)
    input===nothing && return nothing
    producer_junctions,producer_sides=isempty(spec.corners) ?
        _transfinite_chained_tri(m,signed_curves,curve_points,surf,caller) :
        (junctions,curve_points)
    L=length(producer_sides[1])-1;H=length(producer_sides[2])-1
    if haskey(m.periodic,(2,surf))
        producer_junctions,L,H=_model_transfinite_quad_producer_frame(m,surf,input,caller)
    end
    source=_model_transfinite_face_source(input,surf,(L+1)*(H+1),caller)
    found,producer_slot=_model_transfinite_face_slot(ordered,producer_junctions,caller)
    producer_slot==slot || throw(ArgumentError("$caller: actual Surface[$surf] changes canonical slot"))
    M=iseven(found) ? H+1 : L+1
    size(points,2)==(L+1)*(H+1) || throw(ArgumentError("$caller: actual face grid extent differs"))
    ids=zeros(UInt64,size(points,2));written=falses(size(points,2))
    for row in axes(source.coords,2)
        index=row-1;i=index%(L+1);j=index÷(L+1)
        a,b=_volume_face_canonical_coords(found,i,j,L,H);target=a+1+b*M
        written[target] && throw(ArgumentError("$caller: actual face row is emitted twice"))
        point=_model_transfinite_source_point(source,row,geometry_points)
        points[:,target].=point;ids[target]=source.entity_data.external_node_tags[row]
        written[target]=true
    end
    all(written) || throw(ArgumentError("$caller: actual face grid row binding is incomplete"))
    return ids
end

function _model_transfinite_prism_face_source(input,surf,kernel,caller)
    return _model_transfinite_face_source(input,surf,nnodes(kernel),caller)
end

# Compact face expansions repeat actual lattice identities. Number them in
# the existing kernel's first-emission order; coordinates never select an ID.
function _model_transfinite_compact_ids(ids,n,caller)
    rows=Vector{Int32}(undef,length(ids));seen=Dict{UInt64,Int32}()
    for (i,id) in enumerate(ids)
        rows[i]=get!(seen,id) do
            Int32(length(seen)+1)
        end
    end
    length(seen)==(n+1)*(n+2)÷2 || throw(ArgumentError(
        "$caller: compact face identities do not cover its lattice"))
    return rows
end

struct _ModelTransfiniteCubeIndex
    npu::Int
    npv::Int
end
struct _ModelTransfinitePrismIndex
    npv::Int
    layer::Int
end
struct _ModelTransfiniteCompactIndex
    npv::Int
    nw::Int
    lattice::Int
    block::Int
    base5::Int
    lower::Vector{Int32}
    upper::Vector{Int32}
end

@inline _model_transfinite_volume_row(index::_ModelTransfiniteCubeIndex,i,j,k)=
    Int32(i+1+index.npu*(j+index.npv*k))
@inline _model_transfinite_volume_row(index::_ModelTransfinitePrismIndex,i,j,k)=
    Int32(k*index.layer+(i==0 ? 1 : 2+(i-1)*index.npv+j))
@inline function _model_transfinite_volume_row(index::_ModelTransfiniteCompactIndex,i,j,k)
    k==0 && return index.lower[i+1+j*index.npv]
    k==index.nw && return Int32(index.base5)+index.upper[i+1+j*index.npv]
    return Int32(index.lattice+(k-1)*index.block+(i==0 ? 1 : 2+(i-1)*index.npv+j))
end

function _model_transfinite_volume_bindings(mesh,face_ids,counts,caller;
    prism=false,compact=false)
    nu,nv,nw=counts
    tags=zeros(UInt64,nnodes(mesh));retained=falses(nnodes(mesh))
    index=if compact
        nu==nv || throw(ArgumentError("$caller: compact prism grid extents differ"))
        lower=_model_transfinite_compact_ids(face_ids[5],nu,caller)
        upper=_model_transfinite_compact_ids(face_ids[6],nu,caller)
        lattice=(nu+1)*(nu+2)÷2;block=1+nu*(nv+1)
        base5=lattice+(nw-1)*block
        _ModelTransfiniteCompactIndex(nv+1,nw,lattice,block,base5,lower,upper)
    elseif prism
        _ModelTransfinitePrismIndex(nv+1,1+nu*(nv+1))
    else
        _ModelTransfiniteCubeIndex(nu+1,nv+1)
    end
    function bind(slot,column,i,j,k)
        ids=face_ids[slot];ids===nothing && throw(ArgumentError(
            "$caller: actual volume face has no source identity grid"))
        row=_model_transfinite_volume_row(index,i,j,k);0<row<=nnodes(mesh) || throw(ArgumentError(
            "$caller: canonical volume event has no emitted row"))
        id=ids[column]
        !retained[row] || tags[row]==id || throw(ArgumentError(
            "$caller: volume junction receives distinct actual boundary identities"))
        tags[row]=id;retained[row]=true
        return nothing
    end
    for k in 0:nw,i in 0:nu
        bind(1,i+1+k*(nu+1),i,0,k)
        bind(3,i+1+k*(nu+1),i,nv,k)
    end
    for k in 0:nw,j in 0:nv
        bind(2,j+1+k*(nv+1),nu,j,k)
        prism || bind(4,j+1+k*(nv+1),0,j,k)
    end
    for j in 0:nv,i in 0:nu
        bind(5,i+1+j*(nu+1),i,j,0)
        bind(6,i+1+j*(nu+1),i,j,nw)
    end
    allunique(tags[retained]) || throw(ArgumentError(
        "$caller: distinct volume rows repeat an actual boundary identity"))
    return tags,retained
end

function _model_transfinite_volume_production(m,t,mesh,bindings,input,geometry_points,
    node_high_water,element_high_water,caller)
    tags,retained,layout=bindings
    positions=Dict(tag=>row for (row,tag) in enumerate(input.mesh.entity_data.external_node_tags))
    for row in eachindex(tags)
        retained[row] || continue
        source=get(positions,tags[row],nothing)
        source!==nothing && isequal(mesh.coords[:,row],input.mesh.coords[:,source]) ||
            throw(ArgumentError("$caller: actual volume boundary loses its emitted source value"))
    end
    geometry=nothing
    if geometry_points!==nothing && !isempty(geometry_points)
        # Preserve the standalone geometry kernel's arithmetic and cell order.
        geometry=_mesh_model_volume(m,t;_mesh_input=nothing,_geometry_points=nothing,
            _boundary_regrade=false)
        nnodes(geometry)==nnodes(mesh) &&
            _model_transfinite_volume_cells(geometry)==_model_transfinite_volume_cells(mesh) ||
            throw(ArgumentError("$caller: geometry volume does not share actual emitter events"))
        coordinates=copy(geometry.coords)
        for row in eachindex(tags)
            retained[row] || continue
            point=get(geometry_points,tags[row],nothing)
            point===nothing || (coordinates[:,row].=point)
        end
        geometry=geometry isa Mesh ? Mesh(coordinates;tets=geometry.tets,tet_tag=geometry.tet_tag) :
            MixedMesh(coordinates,geometry.blocks)
    end
    emission=_model_transfinite_volume_emission(mesh,layout,caller;
        reverse=get(m.meshing.reverse,(3,t),false))
    return _ModelTransfiniteVolumeProduction(mesh,tags,retained,geometry,emission,
        UInt64(node_high_water),UInt64(0),UInt64[],UInt64(element_high_water))
end

_model_transfinite_volume_cells(mesh::Mesh)=[(4,mesh.tets)]
_model_transfinite_volume_cells(mesh::MixedMesh)=[(Int(block.msh),block.nodes) for block in mesh.blocks]

# The standalone kernels keep their original geometric row/cell storage.
# This receipt translates those exact events to Gmsh's public emission order.
struct _ModelTransfiniteVolumeLayout{I}
    index::I
    counts::NTuple{3,Int}
    prism::Bool
    compact::Bool
    grid::Union{Nothing,_ModelTransfiniteGridReceipt}
end

_ModelTransfiniteVolumeLayout(index,counts,prism,compact)=
    _ModelTransfiniteVolumeLayout(index,counts,prism,compact,nothing)

function _model_transfinite_volume_index(face_ids,counts,caller;prism=false,compact=false)
    nu,nv,nw=counts
    if compact
        nu==nv || throw(ArgumentError("$caller: compact prism grid extents differ"))
        lower=_model_transfinite_compact_ids(face_ids[5],nu,caller)
        upper=_model_transfinite_compact_ids(face_ids[6],nu,caller)
        lattice=(nu+1)*(nu+2)÷2;block=1+nu*(nv+1)
        return _ModelTransfiniteCompactIndex(nv+1,nw,lattice,block,
            lattice+(nw-1)*block,lower,upper)
    elseif prism
        return _ModelTransfinitePrismIndex(nv+1,1+nu*(nv+1))
    end
    return _ModelTransfiniteCubeIndex(nu+1,nv+1)
end

function _model_transfinite_volume_boundary(mesh,face_ids,counts,caller;
    prism=false,compact=false,input=nothing,face_tags=nothing,compact_triangles=nothing)
    tags,retained=_model_transfinite_volume_bindings(mesh,face_ids,counts,caller;prism,compact)
    index=_model_transfinite_volume_index(face_ids,counts,caller;prism,compact)
    grid=input===nothing ? nothing :
        _model_transfinite_grid_receipt(input,face_ids,face_tags,counts,caller;
                                       mesh,compact,compact_triangles)
    return tags,retained,_ModelTransfiniteVolumeLayout(index,counts,prism,compact,grid)
end

function _model_transfinite_tet_emission!(nodes,columns,source,mesh,vertices,reverse,caller)
    key=Tuple(sort!(collect(vertices)))
    column=get(source,key,0)
    column>0 || throw(ArgumentError("$caller: native TF tetrahedron loses its geometric parent"))
    push!(columns,column)
    a,b,c,d=vertices
    sign=orient3(ntuple(axis->mesh.coords[axis,a],3),ntuple(axis->mesh.coords[axis,b],3),
        ntuple(axis->mesh.coords[axis,c],3),ntuple(axis->mesh.coords[axis,d],3))
    sign!=0 || throw(ArgumentError("$caller: native TF tetrahedron has zero volume"))
    # MTetrahedron::reverse swaps the first two vertices. The geometric
    # kernel's positive-orientation permutation stays in its own carrier.
    cell=(sign>0) ⊻ reverse ? (b,a,c,d) : vertices
    nodes[:,length(columns)].=cell
    return nothing
end

function _model_transfinite_six_emissions!(nodes,columns,source,mesh,index,i,j,k,reverse,caller)
    a=_model_transfinite_volume_row(index,i,j,k)
    b=_model_transfinite_volume_row(index,i+1,j,k)
    c=_model_transfinite_volume_row(index,i,j+1,k)
    d=_model_transfinite_volume_row(index,i,j,k+1)
    e=_model_transfinite_volume_row(index,i+1,j,k+1)
    f=_model_transfinite_volume_row(index,i,j+1,k+1)
    g=_model_transfinite_volume_row(index,i+1,j+1,k)
    h=_model_transfinite_volume_row(index,i+1,j+1,k+1)
    for vertices in ((a,b,c,d),(b,c,d,e),(d,e,c,f),(b,c,e,g),(c,f,e,g),(e,f,h,g))
        _model_transfinite_tet_emission!(nodes,columns,source,mesh,vertices,reverse,caller)
    end
    return nothing
end

function _model_transfinite_volume_emission(mesh,layout,caller;reverse=false)
    index=layout.index;nu,nv,nw=layout.counts
    rows=Int32[];seen=falses(nnodes(mesh))
    # meshGRegionTransfinite.cpp emits fresh interiors in i,j,k order;
    # repeated collapsed/compact boundary slots name existing logical rows.
    for i in 0:nu,j in 0:nv,k in 0:nw
        row=_model_transfinite_volume_row(index,i,j,k)
        0<row<=nnodes(mesh) || throw(ArgumentError("$caller: native TF event has no carrier row"))
        seen[row] && continue
        seen[row]=true;push!(rows,row)
    end
    if !all(seen)
        # TransfQuadTri divides cells around centroid births appended after
        # the base lattice rows; emit them after the lattice in carrier
        # (creation) order. A plain Mesh carrier is all-lattice, and an
        # uncovered row at or before the last lattice row is an actual
        # coverage defect, not a centroid.
        lattice_end=findlast(seen)
        (mesh isa MixedMesh && lattice_end!==nothing &&
            all(row->row>lattice_end,findall(!,seen))) ||
            throw(ArgumentError("$caller: native TF rows do not cover the geometric carrier"))
        for row in eachindex(seen)
            seen[row] || (seen[row]=true;push!(rows,row))
        end
    end
    # Recombined block connectivity already has a separate Higher receipt.
    # Forward node births here without claiming a new mixed-cell contract.
    mesh isa Mesh || return _ModelTransfiniteVolumeEmission(0,zeros(Int32,0,0),Int[],rows,layout.grid)
    nodes=Matrix{Int32}(undef,4,size(mesh.tets,2));columns=Int[]
    source=Dict{NTuple{4,Int32},Int}()
    for column in axes(mesh.tets,2)
        key=Tuple(sort!(collect(mesh.tets[:,column])))
        haskey(source,key) && throw(ArgumentError("$caller: TF geometric tetrahedron is repeated"))
        source[key]=column
    end
    if !layout.prism
        for i in 0:nu-1,j in 0:nv-1,k in 0:nw-1
            _model_transfinite_six_emissions!(nodes,columns,source,mesh,index,i,j,k,reverse,caller)
        end
    elseif !layout.compact
        for j in 0:nv-1,k in 0:nw-1
            a=_model_transfinite_volume_row(index,0,0,k)
            b=_model_transfinite_volume_row(index,1,j,k)
            c=_model_transfinite_volume_row(index,1,j+1,k)
            d=_model_transfinite_volume_row(index,0,0,k+1)
            e=_model_transfinite_volume_row(index,1,j,k+1)
            f=_model_transfinite_volume_row(index,1,j+1,k+1)
            for vertices in ((a,b,c,d),(b,c,d,e),(d,f,e,c))
                _model_transfinite_tet_emission!(nodes,columns,source,mesh,vertices,reverse,caller)
            end
        end
        for i in 1:nu-1,j in 0:nv-1,k in 0:nw-1
            _model_transfinite_six_emissions!(nodes,columns,source,mesh,index,i,j,k,reverse,caller)
        end
    else
        for i in 0:nu-1,k in 0:nw-1
            a=_model_transfinite_volume_row(index,i,i,k)
            b=_model_transfinite_volume_row(index,i+1,i,k)
            d=_model_transfinite_volume_row(index,i,i,k+1)
            e=_model_transfinite_volume_row(index,i+1,i,k+1)
            g=_model_transfinite_volume_row(index,i+1,i+1,k)
            h=_model_transfinite_volume_row(index,i+1,i+1,k+1)
            for vertices in ((a,b,g,d),(b,g,d,e),(d,e,g,h))
                _model_transfinite_tet_emission!(nodes,columns,source,mesh,vertices,reverse,caller)
            end
        end
        for i in 1:nu-1,j in 0:i-1,k in 0:nw-1
            a=_model_transfinite_volume_row(index,i,j,k)
            b=_model_transfinite_volume_row(index,i+1,j,k)
            c=_model_transfinite_volume_row(index,i,j+1,k)
            d=_model_transfinite_volume_row(index,i,j,k+1)
            e=_model_transfinite_volume_row(index,i+1,j,k+1)
            f=_model_transfinite_volume_row(index,i,j+1,k+1)
            g=_model_transfinite_volume_row(index,i+1,j+1,k)
            h=_model_transfinite_volume_row(index,i+1,j+1,k+1)
            for vertices in ((a,c,d,g),(c,f,d,g),(d,f,h,g),(a,b,g,d),(b,g,d,e),(d,e,g,h))
                _model_transfinite_tet_emission!(nodes,columns,source,mesh,vertices,reverse,caller)
            end
        end
    end
    length(columns)==size(mesh.tets,2) && allunique(columns) || throw(ArgumentError(
        "$caller: native TF tetrahedron events are not bijective with the carrier"))
    return _ModelTransfiniteVolumeEmission(3,nodes,columns,rows,layout.grid)
end
