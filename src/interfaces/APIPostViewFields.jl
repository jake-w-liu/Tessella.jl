# A model PostView field owns a snapshot of the mesh and data bound to step 0.
# Display selection and the target/current model do not select its source.
# Records retain full geometry and value bases instead of flattening high-order
# cells to their corners or round-tripping through the lossy POS text writer.
using ..SizeField: ConstantSize, GMSH_MAX_SIZE, PostViewAnisoField

function _model_postview_scheme(cell,kind,schemes)
    family=cell.family
    family===:pnt && return nothing
    # Canonical linear cells use the field's compact native connectivity.
    # Quad geometry stays explicit so warped bilinear surfaces are retained.
    msh_spec(cell.msh).order==1 && family!==:qua && return nothing
    return get!(schemes,cell.msh) do
        geometry=_native_view_writer_basis(cell.msh)
        if kind=="ElementData"
            return Post.PosScheme(ones(1,1),zeros(Int,1,3),
                geometry.coefficients,geometry.exponents)
        end
        return Post.PosScheme(geometry.coefficients,geometry.exponents,
            geometry.coefficients,geometry.exponents)
    end
end

function _model_postview_values(table,cell,nval,caller)
    nc=table.components
    if table.kind=="NodeData"
        nval<=length(cell.node_uids) || throw(ArgumentError(
            "$caller: model-view value basis exceeds element $(cell.uid) node count"))
        all(uid->haskey(table.rows,uid),@view(cell.node_uids[1:nval])) || return nothing
        values=Matrix{Float64}(undef,nc,nval)
        for node in 1:nval
            values[:,node]=table.rows[cell.node_uids[node]]
        end
        return values
    end
    row=get(table.rows,cell.uid,nothing)
    row===nothing && return nothing
    if table.kind=="ElementData"
        # Canonical linear cells broadcast the constant to their corners.
        return repeat(reshape(row,nc,1),1,nval)
    end
    length(row)==nc*nval || throw(ArgumentError(
        "$caller: ElementNodeData for element $(cell.uid) has $(length(row)÷nc) " *
        "values but its interpolation basis requires $nval"))
    return reshape(copy(row),nc,nval)
end

function _model_post_view_field(view,tag::Int;
        crop_negative::Bool=true,use_closest::Bool=true)
    caller="build_geo_size_field: Field[$tag]"
    view=_view_data_owner(view)
    _slot_capture_locked!()
    kind=view.model_kind
    kind in ("NodeData","ElementData","ElementNodeData") || throw(ArgumentError(
        "$caller: unsupported model-view data kind"))
    key=(0,kind)
    context=get(view.model_contexts,key,nothing)
    # Like list fields, the sizing operator always samples the first time
    # step. An allocated empty initial step contributes no size constraint;
    # a populated later step and View.TimeStep must not replace it.
    context===nothing && return ConstantSize(GMSH_MAX_SIZE)
    table=_native_view_text_table(view,key)
    isempty(table.rows) && return ConstantSize(GMSH_MAX_SIZE)
    table.components in (1,3,9) || throw(ArgumentError(
        "$caller: model-view data must have 1, 3, or 9 components"))
    slot=context.model_slot
    slot===nothing && throw(ArgumentError("$caller: model view has no source model"))
    mesh=_native_view_writer_mesh(slot,caller)
    cells=_native_view_text_cells(mesh;primary_only=false)
    identities=kind=="NodeData" ? Set(mesh.entity_data.external_node_tags) :
        Set(cell.uid for cell in cells)
    for uid in keys(table.rows)
        uid in identities || throw(ArgumentError(
            "$caller: model view contains unknown " *
            "$(kind=="NodeData" ? "node" : "element") UID $uid"))
    end
    schemes=Dict{Int,Post.PosScheme}()
    records=Post.PosElement[]
    for cell in cells
        # `_LIST_MSH_FAMILY` maps record names to mesh families, so resolve its
        # inverse without guessing primary arities or simplex subdivisions.
        cellkind=nothing
        for (name,family) in pairs(_LIST_MSH_FAMILY)
            family===cell.family && (cellkind=name;break)
        end
        cellkind===nothing && throw(ArgumentError(
            "$caller: model-view element family $(cell.family) is unsupported"))
        # User-set list interpolation matrices do not override model element
        # bases (confirmed through the public Gmsh 4.15.2 SDK).
        scheme=_model_postview_scheme(cell,kind,schemes)
        ngeo=scheme===nothing ? length(cell.nodes) : size(scheme.coefgeo,1)
        nval=scheme===nothing ? length(cell.nodes) : size(scheme.coefval,1)
        ngeo<=length(cell.nodes) || throw(ArgumentError(
            "$caller: model-view geometry basis exceeds element $(cell.uid) node count"))
        values=_model_postview_values(table,cell,nval,caller)
        values===nothing && continue
        coordinates=mesh.coords[:,cell.nodes[1:ngeo]]
        push!(records,Post.PosElement(cellkind,coordinates,values;scheme))
    end
    isempty(records) && return ConstantSize(GMSH_MAX_SIZE)
    # One model-data step has a uniform component count. Construct directly
    # from the owned records: POS serialization would require rebinding family
    # schemes at every mixed-order transition and is not part of this path.
    field=PostViewField(records;crop_negative,use_closest)
    return table.components==9 ? PostViewAnisoField(field) : field
end
