# Explicit cell/row frontend shared by native API, model records and GEO.
# Geometry providers build an STL from CAD only when the actual entity has one.
function _model_box_edge_geometry(m,edge)
    origin=m.points[edge[1]]
    endpoint=m.points[edge[2]]
    delta=ntuple(index->endpoint[index]-origin[index],3)
    length=hypot(delta...)
    isfinite(length) && length>0 || throw(ArgumentError(
        "add_box!: actual Box edge length must be finite and positive"))
    direction=ntuple(index->delta[index]/length,3)
    return (occ=:line,origin=origin,dir=direction,t0=0.0,t1=length,primitive=:box)
end

function _model_box_plane_geometry(m,tag,origin,extents,face::Int)
    axis=(face+1)÷2
    reversed=isodd(face)
    location=ntuple(index->origin[index]+(index==axis && !reversed ? extents[index] : 0.0),3)
    normal=ntuple(index->index==axis ? 1.0 : 0.0,3)
    X=axis==3 ? (1.0,0.0,0.0) : (0.0,0.0,1.0)
    Y=axis==1 ? (0.0,-1.0,0.0) : axis==2 ? (1.0,0.0,0.0) : (0.0,1.0,0.0)
    lower=axis==1 ? [0.0,-extents[2]] : [0.0,0.0]
    upper=axis==1 ? [extents[3],0.0] : axis==2 ?
        [extents[3],extents[1]] : [extents[1],extents[2]]
    geometry=(occ=:plane,center=location,axis=normal,X=X,Y=Y,reversed=reversed,
        pcurves=Dict{Int,NamedTuple}(),uvb=(lower,upper),primitive=:box)
    for loop in m.surfaces[tag],signed_curve in m.loops[loop]
        curve=abs(signed_curve)
        edge,t0,t1=_brep_edge_curve(m,curve)
        pc=_brep_synth_pcurve(geometry,edge,t0,t1)
        geometry.pcurves[curve]=(fwd=pc,rev=nothing)
    end
    return geometry
end

function _model_orientation_canonical_plane(m,tag,mesh,caller)
    geometry=get(m.surface_geometry,tag,nothing)
    geometry!==nothing && hasproperty(geometry,:occ) && geometry.occ===:plane || return mesh
    cells=if mesh isa Mesh
        mesh.tris
    else
        block=findfirst(block->msh_dimension(block.msh)==2 && !isempty(block.nodes),mesh.blocks)
        block===nothing && return mesh
        mesh.blocks[block].nodes
    end
    isempty(cells) && return mesh
    normal=_occ_surface_normal(geometry,0.0,0.0)
    a=ntuple(axis->mesh.coords[axis,cells[1,1]],3)
    b=ntuple(axis->mesh.coords[axis,cells[2,1]],3)
    c=ntuple(axis->mesh.coords[axis,cells[3,1]],3)
    dot=_model_orientation_dot(_model_orientation_cross(
        _model_orientation_sub(b,a),_model_orientation_sub(c,a)),normal)
    return dot>0 ? mesh : _reversed_surface_mesh(mesh,caller,tag)
end

function _model_orientation_surface_mode(m::GeoModel,tag::Int)
    haskey(m.periodic,(2,tag)) && return :periodic
    haskey(m.discrete,(2,tag)) && return :fully_discrete
    _surface_type(m,tag)==:boundary_layer && return :boundary_layer_surface
    return :cad
end

function _model_orientation_entities(m,regions)
    surfaces=Set{Int}()
    for region in regions,(dimension,tag) in model_boundary(m,[(3,region)],false,false,false)
        dimension==2 && push!(surfaces,tag)
    end
    return surfaces,Set{Int}(regions)
end

function _model_orientation_occ_plane(m::GeoModel,tag::Int)
    geometry=get(m.surface_geometry,tag,nothing)
    if geometry!==nothing && hasproperty(geometry,:occ)
        return geometry.occ===:plane
    end
    return false
end

function _model_orientation_cad_stl(m::GeoModel,tag::Int;old_ruled_surface::Bool=false)
    caller="Model mesh orientation CAD STL"
    geometry=get(m.surface_geometry,tag,nothing)
    occ=geometry!==nothing && hasproperty(geometry,:occ)
    _model_orientation_occ_plane(m,tag) || begin
        occ && throw(ArgumentError(
            "$caller: Surface $tag requires the native trimmed $(geometry.occ) CAD STL provider"))
        return nothing
    end
    loops=m.surfaces[tag]
    length(loops)==1 || throw(ArgumentError(
        "$caller: Surface $tag requires a trimmed planar CAD STL with interior loops"))
    all(curve->_curve_type(m,abs(curve))===:line,m.loops[only(loops)]) ||
        throw(ArgumentError("$caller: Surface $tag requires a curved-boundary planar CAD STL"))
    vertices=first(_surface_loop_vertex_chain(m,only(loops),caller))
    length(vertices)>=3 || throw(ArgumentError("$caller: Surface $tag has no planar polygon"))
    values=model_normal(m,tag,[0.0,0.0];old_ruled_surface)
    normal=(values[1],values[2],values[3])
    # A convex trimmed Plane has an exact polygon triangulation. Concave and
    # multiply trimmed CAD faces require their own admitted provider.
    sign=0
    for index in eachindex(vertices)
        a,b,c=vertices[index],vertices[mod1(index+1,length(vertices))],
              vertices[mod1(index+2,length(vertices))]
        turn=_model_orientation_dot(_model_orientation_cross(
            _model_orientation_sub(b,a),_model_orientation_sub(c,b)),normal)
        iszero(turn) && continue
        current=turn>0 ? 1 : -1
        if sign==0
            sign=current
        elseif sign!=current
            throw(ArgumentError("$caller: Surface $tag requires a concave planar CAD STL provider"))
        end
    end
    sign!=0 || throw(ArgumentError("$caller: Surface $tag has degenerate CAD boundary"))
    triangles=_ModelMeshOrientationTriangle[]
    sizehint!(triangles,length(vertices)-2)
    for index in 2:length(vertices)-1
        a,b,c=vertices[1],vertices[index],vertices[index+1]
        push!(triangles,sign>0 ? (a,b,c) : (a,c,b))
    end
    return triangles
end

function _model_orientation_stl(m::GeoModel,tag,coords,cells,cache;
        old_ruled_surface::Bool=false)
    previous=get(cache,tag,nothing)
    previous!==nothing && !isempty(previous) && return previous
    cad=_model_orientation_cad_stl(m,tag;old_ruled_surface)
    triangles=cad===nothing ? _ModelMeshOrientationTriangle[] : cad
    if cad===nothing
        for cell in cells
            cell.msh in (2,9) || continue
            push!(triangles,ntuple(slot->ntuple(axis->
                coords[axis,cell.nodes[slot]],3),3))
        end
    end
    isempty(triangles) || (cache[tag]=triangles)
    return triangles
end

function _model_orientation_node_uv(m,tag,node,coords,owners,parameters,
        old_ruled_surface)
    dimension,entity=owners[node]
    parameter=parameters[node]
    dimension==2 && entity==tag && parameter!==nothing && length(parameter)>=2 &&
        return parameter[1:2]
    if dimension==0 && haskey(m.points,Int(entity))
        return model_reparametrize_on_surface(m,0,Int(entity),Float64[],tag;
            old_ruled_surface)
    elseif dimension==1 && haskey(m.curves,Int(entity))
        # A plain vertex on a native Curve has no stored U. Upstream returns
        # false here instead of projecting its edited Cartesian coordinates.
        parameter!==nothing && !isempty(parameter) || return nothing
        return model_reparametrize_on_surface(m,1,Int(entity),parameter[1:1],tag;
            old_ruled_surface)
    end
    return model_parametrization(m,2,tag,collect(@view coords[:,node]);old_ruled_surface)
end

function _model_orientation_normals(m,tag,coords,cells,owners,parameters;
        old_ruled_surface::Bool=false,node_tags=nothing,diagnostics=nothing)
    # GFace stores all Triangles before all Quadrangles, including when raw
    # cells or imported blocks were supplied in the opposite family order.
    cells=sort(cells;by=cell->cell.msh in (2,9) ? 0 : 1,alg=Base.Sort.MergeSort)
    _model_orientation_surface_mode(m,tag)==:cad || return cells
    result=_ModelMeshOrientationCell[]
    sizehint!(result,length(cells))
    known_nonbl=false
    known_bl=!any(cell->cell.boundary_layer,cells)
    for cell in cells
        vertex=nothing
        # getGFaceOrientation evaluates a normal only while that BL/non-BL
        # group's orientation is unknown. A found normal determines its sign,
        # including a zero dot product, before later cells can be evaluated.
        if !(cell.boundary_layer ? known_bl : known_nonbl)
            for node in cell.nodes
                owners[node]==(2,tag) && parameters[node]!==nothing &&
                    length(parameters[node])>=2 || continue
                values=model_normal(m,tag,parameters[node][1:2];old_ruled_surface)
                vertex=(values[1],values[2],values[3])
                if cell.boundary_layer;known_bl=true;else;known_nonbl=true;end
                break
            end
        end
        push!(result,_ModelMeshOrientationCell(cell.key,cell.msh,cell.nodes;
            vertex_normal=vertex,
            boundary_layer=cell.boundary_layer))
    end
    for index in eachindex(result)
        cell=result[index]
        (cell.boundary_layer ? known_bl : known_nonbl) && continue
        uv=[0.0,0.0]
        complete=true
        for node in cell.nodes
            parameter=_model_orientation_node_uv(m,tag,node,coords,owners,parameters,
                                                 old_ruled_surface)
            if parameter===nothing
                if diagnostics!==nothing
                    label=node_tags===nothing ? "a mesh node" : "node $(node_tags[node])"
                    push!(diagnostics,"No parametric coordinate on $label classified on curve $(owners[node][2])")
                end
                complete=false
                break
            end
            length(parameter)==2 || throw(ArgumentError(
                "Model mesh orientation: Surface $tag node parameters must have two entries"))
            uv[1]+=parameter[1];uv[2]+=parameter[2]
        end
        complete || continue
        uv./=length(cell.nodes)
        values=model_normal(m,tag,uv;old_ruled_surface)
        result[index]=_ModelMeshOrientationCell(cell.key,cell.msh,cell.nodes;
            vertex_normal=cell.vertex_normal,barycenter_normal=(values[1],values[2],values[3]),
            boundary_layer=cell.boundary_layer)
        if cell.boundary_layer;known_bl=true;else;known_nonbl=true;end
    end
    return result
end

function _model_orientation_region_scale(m,tag,coords,volume_nodes)
    if haskey(m.discrete,(3,tag))
        nodes=get(volume_nodes,tag,Int32[])
        isempty(nodes) && return Inf
        lower=fill(Inf,3);upper=fill(-Inf,3)
        for node in nodes,axis in 1:3
            lower[axis]=min(lower[axis],coords[axis,node])
            upper[axis]=max(upper[axis],coords[axis,node])
        end
        return hypot((upper.-lower)...)
    end
    bounds=model_bounding_box(m,3,tag)
    return hypot(ntuple(axis->bounds[axis+3]-bounds[axis],3)...)
end

function _model_prepare_mesh_orientation(m::GeoModel,regions,coords,surface_cells,
        node_owners,node_parameters;volume_nodes=Dict{Int,Vector{Int32}}(),
        old_ruled_surface::Bool=false,node_tags=nothing,
        rng::Random.AbstractRNG=Random.Xoshiro(0x4d657368))
    length(node_owners)==size(coords,2)==length(node_parameters) || throw(ArgumentError(
        "Model mesh orientation: node ownership and stored parameters do not align"))
    node_tags===nothing || length(node_tags)==size(coords,2) || throw(ArgumentError(
        "Model mesh orientation: public node labels do not align with coordinate rows"))
    cache=copy(m.meshing.orientation_stl)
    flags=Dict(tag=>flag for ((dim,tag),flag) in m.meshing.reverse if dim==2)
    changed=Dict{Int,Vector{Int32}}();diagnostics=String[]
    inputs=Dict{Int,Vector{_ModelMeshOrientationCell}}()
    for region in regions
        boundaries=model_boundary(m,[(3,region)],false,false,false)
        surfaces=_ModelMeshOrientationSurface[]
        for (dimension,tag) in boundaries
            dimension==2 || throw(ArgumentError(
                "Model mesh orientation: Volume $region has a non-Surface boundary"))
            cells=get!(inputs,tag) do
                _model_orientation_normals(m,tag,coords,
                    get(surface_cells,tag,_ModelMeshOrientationCell[]),node_owners,
                    node_parameters;old_ruled_surface,node_tags,diagnostics)
            end
            stl=_model_orientation_stl(m,tag,coords,cells,cache;old_ruled_surface)
            push!(surfaces,_ModelMeshOrientationSurface(tag,coords,cells,stl;
                mode=_model_orientation_surface_mode(m,tag)))
        end
        result=_model_mesh_orientation_plan(region,surfaces;
            scaling=_model_orientation_region_scale(m,region,coords,volume_nodes),
            reverse_flags=flags,rng)
        flags=result.reverse_flags
        append!(diagnostics,result.diagnostics)
        merge!(changed,result.cells)
        # Shared Surface cells can be visited by another requested region.
        # Its current orientation changes while its first retained STL does not.
        for surface in surfaces
            cells=inputs[surface.tag]
            for index in eachindex(cells)
                cell=cells[index]
                nodes=get(result.cells,cell.key,nothing)
                nodes===nothing && continue
                cells[index]=_ModelMeshOrientationCell(cell.key,cell.msh,nodes;
                    vertex_normal=cell.vertex_normal,barycenter_normal=cell.barycenter_normal,
                    boundary_layer=cell.boundary_layer)
            end
        end
    end
    reverse=copy(m.meshing.reverse)
    for (tag,flag) in flags;reverse[(2,tag)]=flag;end
    return (cells=changed,reverse=reverse,stl=cache,diagnostics=diagnostics)
end

function _model_publish_mesh_orientation!(m::GeoModel,plan)
    m.meshing.reverse=plan.reverse
    m.meshing.orientation_stl=plan.stl
    return nothing
end

function _model_orientation_mixed_input(m,mesh::MixedMesh;covered=(),regions=nothing,
        parameters=_model_mesh_runtime_parameters(m,mesh))
    data=mesh.entity_data
    data===nothing && throw(ArgumentError(
        "Model mesh orientation: actual mesh requires owned entity metadata"))
    records=[(owner,record) for (owner,record) in _discrete_mesh_records_model(m)
        if !(owner in covered)]
    count=sum(size(block.nodes,2) for block in mesh.blocks;init=0)
    wanted_surfaces,wanted_volumes=regions===nothing ? (nothing,nothing) :
        _model_orientation_entities(m,regions)
    input=_model_orientation_record_input(m;coords=mesh.coords,
        node_tags=data.external_node_tags,owners=data.node_entities,
        parameters,key_offset=count,records,
        surface_tags=wanted_surfaces,volume_tags=wanted_volumes)
    offset=0
    for (index,block) in enumerate(mesh.blocks)
        dimension=msh_dimension(block.msh)
        if dimension in (2,3)
            for column in axes(block.nodes,2)
                tag=Int(data.block_entities[index][column])
                wanted=dimension==2 ? wanted_surfaces : wanted_volumes
                wanted!==nothing && !(tag in wanted) && continue
                if dimension==2
                    push!(get!(()->_ModelMeshOrientationCell[],input.surface_cells,tag),
                        _ModelMeshOrientationCell(offset+column,block.msh,
                            @view(block.nodes[:,column])))
                else
                    append!(get!(()->Int32[],input.volume_nodes,tag),@view(block.nodes[:,column]))
                end
            end
        end
        offset+=size(block.nodes,2)
    end
    return input
end

function _model_orientation_mixed_replacement(mesh::MixedMesh,changes)
    offset=0;changed=false;blocks=ElementBlock[]
    for block in mesh.blocks
        cells=copy(block.nodes)
        for column in axes(cells,2)
            nodes=get(changes,offset+column,nothing)
            nodes===nothing && continue
            nodes==@view(cells[:,column]) && continue
            cells[:,column]=nodes;changed=true
        end
        push!(blocks,ElementBlock(block.msh,cells,block.tags))
        offset+=size(cells,2)
    end
    changed || return mesh
    return MixedMesh(mesh.coords,blocks;physical_names=mesh.physical_names,
        entity_data=mesh.entity_data,elementary_entities=mesh.elementary_entities,
        periodic_links=mesh.periodic_links,ancillary_sections=mesh.ancillary_sections,
        data_sections=mesh.data_sections,partition_data=mesh.partition_data)
end

function _model_orientation_carrier(actual::MixedMesh,replacement,carrier,rows,changes;
        cell_tags=Vector{UInt64}[])
    carrier===nothing && return nothing
    carrier===actual && return replacement
    replacement===actual && return carrier
    isempty(cell_tags) || return _model_orientation_traced_carrier(
        actual,replacement,carrier,rows,changes,cell_tags)
    carrier isa Mesh || throw(ArgumentError(
        "Model mesh orientation: a detached non-simplex carrier needs producer cell bindings"))
    length(rows)==nnodes(carrier) && allunique(rows) &&
        all(row->0<row<=nnodes(actual),rows) || throw(ArgumentError(
            "Model mesh orientation: carrier rows must bind actual node identities"))
    triangles=Dict{NTuple{3,Int32},Int}();duplicates=Set{NTuple{3,Int32}}()
    changed_triangles=Set{Int}()
    changed_quads=NTuple{4,Int32}[]
    offset=0
    for block in actual.blocks
        if block.msh in (2,9)
            for column in axes(block.nodes,2)
                key=Tuple(sort(block.nodes[1:3,column]))
                haskey(triangles,key) && push!(duplicates,key)
                triangles[key]=offset+column
                nodes=get(changes,offset+column,nothing)
                nodes!==nothing && nodes!=@view(block.nodes[:,column]) &&
                    push!(changed_triangles,offset+column)
            end
        elseif block.msh in (3,10,16)
            for column in axes(block.nodes,2)
                nodes=get(changes,offset+column,nothing)
                (nodes===nothing || nodes==@view(block.nodes[:,column])) && continue
                push!(changed_quads,Tuple(block.nodes[1:4,column]))
            end
        end
        offset+=size(block.nodes,2)
    end
    cells=copy(carrier.tris);changed=false
    for column in axes(cells,2)
        corners=ntuple(slot->Int32(rows[cells[slot,column]]),3)
        key=Tuple(sort(collect(corners)))
        key in duplicates && throw(ArgumentError(
            "Model mesh orientation: carrier Triangle has ambiguous actual cell topology; producer binding is required"))
        element=get(triangles,key,0)
        if element==0
            any(quad->all(node->node in quad,corners),changed_quads) && throw(ArgumentError(
                "Model mesh orientation: triangulated Quadrangle carrier requires its explicit parent-cell binding"))
            continue
        end
        element in changed_triangles || continue
        cells[2,column],cells[3,column]=cells[3,column],cells[2,column]
        changed=true
    end
    changed || return carrier
    return Mesh(carrier.coords;segs=carrier.segs,tris=cells,tets=carrier.tets,
        seg_tag=carrier.seg_tag,tri_tag=carrier.tri_tag,tet_tag=carrier.tet_tag)
end

function _model_orientation_authority(previous,replacement,carrier)
    rebound=_model_mesh_authority(replacement,previous.done;
        private_rank=previous.private_rank,carrier,carrier_rows=previous.carrier_rows)
    # Rows, labels and provenance are unchanged. Preserve all owned fields,
    # including producer row bindings supplied by the higher-dimensional phase.
    names=fieldnames(typeof(previous))
    values=ntuple(length(names)) do index
        name=names[index]
        name===:mesh ? replacement : name===:products ? rebound.products :
        name===:carrier ? carrier : getfield(previous,name)
    end
    return typeof(previous)(values...)
end

# Record ingestion keeps actual source labels and connectivity, including raw
# cells referencing cache-owned rows supplied by the API caller.
function _model_orientation_record_input(m::GeoModel;coords=zeros(Float64,3,0),
        node_tags=UInt64[],owners=Tuple{Int,Int}[],
        parameters=Union{Nothing,Vector{Float64}}[],key_offset::Int=0,
        records=sort!(_discrete_mesh_records_model(m);by=first),
        surface_tags=nothing,volume_tags=nothing)
    coordinates=NTuple{3,Float64}[ntuple(axis->coords[axis,node],3)
        for node in axes(coords,2)]
    owned=Tuple{Int,Int}[(Int(d),Int(t)) for (d,t) in owners]
    stored=Union{Nothing,Vector{Float64}}[p for p in parameters]
    rows=Dict(UInt64(tag)=>Int32(index) for (index,tag) in enumerate(node_tags))
    for ((dimension,tag),record) in records
        for (index,node) in enumerate(record.node_tags)
            haskey(rows,UInt64(node)) && continue
            length(coordinates)<typemax(Int32) || throw(ArgumentError(
                "Model mesh orientation: source node rows exceed Int32"))
            rows[UInt64(node)]=Int32(length(coordinates)+1)
            push!(coordinates,ntuple(axis->record.node_coords[axis,index],3))
            push!(owned,(dimension,tag))
            parameter=_model_mesh_node_parameters(record,index)
            push!(stored,parameter===nothing ? nothing : copy(parameter))
        end
    end
    surface_cells=Dict{Int,Vector{_ModelMeshOrientationCell}}()
    volume_nodes=Dict{Int,Vector{Int32}}()
    bindings=Dict{Int,Tuple{DiscreteEntity,Int}}()
    key=key_offset
    for ((dimension,tag),record) in records
        dimension in (2,3) || continue
        wanted=dimension==2 ? surface_tags : volume_tags
        wanted!==nothing && !(tag in wanted) && continue
        for index in eachindex(record.element_tags)
            nodes=Int32[]
            for node in record.element_nodes[index]
                row=get(rows,UInt64(node),Int32(0))
                row>0 || throw(ArgumentError(
                    "Model mesh orientation: element $(record.element_tags[index]) references unknown node $node"))
                push!(nodes,row)
            end
            if dimension==2
                key=Base.checked_add(key,1)
                push!(get!(()->_ModelMeshOrientationCell[],surface_cells,tag),
                    _ModelMeshOrientationCell(key,record.element_types[index],nodes))
                bindings[key]=(record,index)
            else
                append!(get!(()->Int32[],volume_nodes,tag),nodes)
            end
        end
    end
    matrix=Matrix{Float64}(undef,3,length(coordinates))
    for (index,point) in enumerate(coordinates),axis in 1:3
        matrix[axis,index]=point[axis]
    end
    public_tags=Vector{UInt64}(undef,length(coordinates))
    for (tag,row) in rows;public_tags[row]=tag;end
    return (coords=matrix,owners=owned,parameters=stored,rows=rows,node_tags=public_tags,
        surface_cells=surface_cells,volume_nodes=volume_nodes,bindings=bindings)
end

function _model_orient_record_mesh!(m::GeoModel,tag::Int)
    surfaces,volumes=_model_orientation_entities(m,[tag])
    input=_model_orientation_record_input(m;surface_tags=surfaces,volume_tags=volumes)
    plan=_model_prepare_mesh_orientation(m,[tag],input.coords,input.surface_cells,
        input.owners,input.parameters;volume_nodes=input.volume_nodes,node_tags=input.node_tags)
    inverse=Dict(row=>Int32(node) for (node,row) in input.rows)
    replacements=Tuple{DiscreteEntity,Int,Vector{Int32}}[]
    for (key,nodes) in plan.cells
        record,index=input.bindings[key]
        push!(replacements,(record,index,Int32[inverse[row] for row in nodes]))
    end
    _model_publish_mesh_orientation!(m,plan)
    for (record,index,nodes) in replacements;record.element_nodes[index]=nodes;end
    for message in plan.diagnostics;@warn message;end
    return nothing
end

function _model_orientation_carrier_blocks(mesh::Mesh)
    return ((1,mesh.segs,mesh.seg_tag),(2,mesh.tris,mesh.tri_tag),
            (4,mesh.tets,mesh.tet_tag))
end
function _model_orientation_carrier_blocks(mesh::MixedMesh)
    return [(Int(block.msh),block.nodes,block.tags) for block in mesh.blocks]
end

function _model_orientation_traced_carrier(actual,replacement,carrier,rows,changes,cell_tags)
    caller="Model mesh orientation"
    blocks=_model_orientation_carrier_blocks(carrier)
    length(cell_tags)==length(blocks) && all(length(tags)==size(block[2],2)
        for (tags,block) in zip(cell_tags,blocks)) || throw(ArgumentError(
        "$caller: carrier cell bindings do not align with geometry blocks"))
    length(rows)==nnodes(carrier) && allunique(rows) &&
        all(row->0<row<=nnodes(actual),rows) || throw(ArgumentError(
        "$caller: carrier rows must bind actual node identities"))
    catalog=Dict{UInt64,Tuple{Int,Int,Int}}()
    changed_parents=Tuple{Symbol,Vector{Int32}}[]
    offset=0
    for (index,block) in enumerate(actual.blocks)
        if msh_dimension(block.msh)==2
            for column in axes(block.nodes,2)
                tag=actual.entity_data.external_element_tags[index][column]
                haskey(catalog,tag) && throw(ArgumentError(
                    "$caller: actual element identity is ambiguous"))
                catalog[tag]=(index,column,offset+column)
                nodes=get(changes,offset+column,nothing)
                parent=@view block.nodes[:,column]
                if nodes!==nothing && nodes!=parent
                    permutation=_model_orientation_permutation(block.msh,caller)
                    nodes==parent[collect(permutation)] || throw(ArgumentError(
                        "$caller: actual parent change is not an orientation reversal"))
                    parent_family=msh_family(block.msh)
                    parent_width=parent_family===:tri ? 3 : parent_family===:qua ? 4 : 0
                    parent_width>0 || throw(ArgumentError(
                        "$caller: unsupported actual Surface cell family"))
                    push!(changed_parents,(parent_family,collect(@view parent[1:parent_width])))
                end
            end
        end
        offset+=size(block.nodes,2)
    end
    prepared=Matrix{Int32}[]
    changed=false
    for (index,(msh,cells,_)) in enumerate(blocks)
        updated=cells
        if msh_dimension(msh)==2
            for column in axes(cells,2)
                family=msh_family(msh)
                primary=family===:tri ? 3 : family===:qua ? 4 : 0
                primary>0 || throw(ArgumentError(
                    "$caller: unsupported carrier Surface cell family"))
                bound=rows[cells[:,column]]
                if cell_tags[index][column]===nothing
                    # A missing binding supplies no producer identity. Detect affected
                    # aliases only to reject the operation; never choose a parent
                    # or update a cell from its topology.
                    affected=any(changed_parents) do (parent_family,corners)
                        compatible=family===parent_family ||
                            (family===:tri && parent_family===:qua && msh==2)
                        compatible && all(node->node in corners,@view bound[1:primary])
                    end
                    affected && throw(ArgumentError(
                        "$caller: changed Surface cell touches a carrier without its producer binding"))
                    continue
                end
                binding=get(catalog,cell_tags[index][column],nothing)
                binding===nothing && throw(ArgumentError(
                    "$caller: geometry Surface cell loses its actual parent identity"))
                bi,source_column,key=binding
                block=actual.blocks[bi]
                actual_family=msh_family(block.msh)
                proxy=family===:tri && actual_family===:qua && msh==2
                family==actual_family || proxy || throw(ArgumentError(
                    "$caller: geometry cell family disagrees with its actual parent"))
                actual_primary=actual_family===:tri ? 3 : actual_family===:qua ? 4 : 0
                primary>0 && actual_primary>0 || throw(ArgumentError(
                    "$caller: unsupported carrier Surface cell family"))
                parent=@view block.nodes[:,source_column]
                all(node->node in parent,bound) &&
                    all(node->node in @view(parent[1:actual_primary]),
                        @view(bound[1:primary])) || throw(ArgumentError(
                    "$caller: geometry cell rows disagree with their actual parent"))
                nodes=get(changes,key,nothing)
                (nodes===nothing || nodes==parent) && continue
                actual_permutation=_model_orientation_permutation(block.msh,caller)
                nodes==parent[collect(actual_permutation)] || throw(ArgumentError(
                    "$caller: actual parent change is not an orientation reversal"))
                updated===cells && (updated=copy(cells))
                updated[:,column]=cells[collect(_model_orientation_permutation(msh,caller)),column]
                changed=true
            end
        end
        push!(prepared,updated)
    end
    changed || return carrier
    if carrier isa Mesh
        return Mesh(carrier.coords;segs=prepared[1],tris=prepared[2],tets=prepared[3],
            seg_tag=carrier.seg_tag,tri_tag=carrier.tri_tag,tet_tag=carrier.tet_tag)
    end
    return MixedMesh(carrier.coords,
        ElementBlock[ElementBlock(block.msh,cells,block.tags)
            for (block,cells) in zip(carrier.blocks,prepared)];
        physical_names=carrier.physical_names,entity_data=carrier.entity_data,
        elementary_entities=carrier.elementary_entities,periodic_links=carrier.periodic_links,
        ancillary_sections=carrier.ancillary_sections,data_sections=carrier.data_sections,
        partition_data=carrier.partition_data)
end
