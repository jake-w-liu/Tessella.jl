# View queries retain the complete value and geometry bases. The shared
# safeguarded scheme inverse supplies intrinsic coordinates and certifies the
# actual physical map; derivatives are pulled back with the geometry Jacobian.
using LinearAlgebra: svd

const _VIEW_PROBE_MSH=(point=15,line=1,triangle=2,quadrangle=3,
    tetrahedron=4,hexahedron=5,prism=6,pyramid=7)

_view_probe_exponents(matrix)=NTuple{3,Float64}[
    ntuple(axis->axis<=size(matrix,2) ? Float64(matrix[row,axis]) : 0.0,3)
    for row in axes(matrix,1)]

function _view_probe_record(kind,coordinates,values,basis,token,caller;identity=nothing)
    if kind===:point
        return (;kind,dim=0,coords=coordinates,values,cell=nothing,partition=true,token,identity)
    end
    cell=SizeField._PVSchemeCell(SizeField._POSTVIEW_SCHEME_KIND[kind],
        isempty(basis.coefgeo),kind===:pyramid,basis.coefval,basis.expval,
        basis.coefgeo,basis.expgeo,coordinates,values)
    SizeField._pv_scheme_validate(cell,1,caller)
    return (;kind,dim=SizeField._pv_scheme_dim(cell.kind),coords=coordinates,
        values,cell,partition=basis.partition,token,identity)
end

function _view_probe_basis(msh,constant,bases)
    return get!(bases,(msh,constant)) do
        native=_native_view_writer_basis(msh)
        (partition=true,coefval=constant ? ones(1,1) : native.coefficients,
         expval=constant ? [(0.0,0.0,0.0)] : _view_probe_exponents(native.exponents),
         coefgeo=native.coefficients,expgeo=_view_probe_exponents(native.exponents))
    end
end

function _view_probe_model_records(view,step,meshes,bases,caller;
        components=nothing,geometry_only=false,include_missing=false)
    kind=view.model_kind
    range=_native_view_step_range(view.model_step_ranges,step)
    range===nothing && return NamedTuple[]
    context=get(view.model_contexts,(step,kind),nothing)
    slot=context===nothing ? range.model_slot : context.model_slot
    mesh=get!(meshes,slot) do
        _native_view_writer_mesh(slot,caller)
    end
    table=context===nothing ? nothing : _native_view_text_table(view,(step,kind))
    nc=components===nothing ? (table===nothing ? range.components : table.components) : components
    nc>0 || throw(ArgumentError("$caller: invalid model-view component count"))
    records=NamedTuple[]
    for cell in _native_view_text_cells(mesh;primary_only=false)
        shape=get(_PLUGIN_FAMILY_KIND,cell.family,nothing)
        shape===nothing && throw(ArgumentError(
            "$caller: element family $(cell.family) has no supported reference map"))
        constant=kind=="ElementData"
        count=constant ? 1 : length(cell.nodes)
        values=zeros(nc,count)
        present=table===nothing || kind=="NodeData"
        if table!==nothing
            width=min(nc,table.components)
            if kind=="NodeData"
                for (node,uid) in enumerate(cell.node_uids)
                    row=get(table.rows,uid,nothing)
                    row===nothing && continue
                    values[1:width,node].=row[1:width];present=true
                end
            else
                row=get(table.rows,cell.uid,nothing)
                if row!==nothing
                    length(row)%table.components==0 || throw(ArgumentError(
                        "$caller: incomplete ElementNodeData component row"))
                    rows=constant ? 1 : min(count,length(row)÷table.components)
                    for node in 1:rows
                        values[1:width,node].=row[(table.components*(node-1)+1):(table.components*(node-1)+width)]
                    end
                    present=true
                end
            end
        end
        present || geometry_only || include_missing || continue
        basis=shape===:point ? nothing : _view_probe_basis(cell.msh,constant,bases)
        token=(cell.uid,cell.family)
        identity=(element=cell.uid,nodes=copy(cell.node_uids))
        push!(records,_view_probe_record(shape,mesh.coords[:,cell.nodes],
            values,basis,token,caller;identity))
    end
    return records
end

# All-step model queries retain the first frame's material nodes and basis.
# Missing rows in a populated later frame retain the preceding material row;
# an allocated empty time step emits zeros without erasing that row history.
function _view_probe_material_values!(values,record,view,index)
    key=(index,view.model_kind)
    haskey(view.model_contexts,key) || return false
    table=_native_view_text_table(view,key)
    width=min(size(values,1),table.components)
    if view.model_kind=="NodeData"
        for (node,uid) in enumerate(record.identity.nodes)
            row=get(table.rows,uid,nothing)
            row===nothing || (values[1:width,node].=row[1:width])
        end
    else
        row=get(table.rows,record.identity.element,nothing)
        if row!==nothing
            count=min(size(values,2),length(row)÷table.components)
            for node in 1:count
                values[1:width,node].=row[table.components*(node-1)+1:table.components*(node-1)+width]
            end
        end
    end
    return true
end

function _view_probe_material_record(record,values)
    cell=record.cell
    if cell!==nothing
        cell=SizeField._PVSchemeCell(cell.kind,cell.flat,cell.pyramid,
            cell.coefval,cell.expval,cell.coefgeo,cell.expgeo,cell.coords,values)
    end
    return merge(record,(;values,cell))
end

function _view_probe_list_records(view,step,bases,caller;elements=nothing)
    records=NamedTuple[]
    elements===nothing && (elements=_native_list_view_records(view,caller))
    for (index,element) in enumerate(elements)
        step<size(element.values,3) || continue
        scheme=element.scheme
        basis=if scheme===nothing
            element.kind===:point ? nothing : _view_probe_basis(_VIEW_PROBE_MSH[element.kind],false,bases)
        else
            geometry=isempty(scheme.coefgeo) ?
                _view_probe_basis(_VIEW_PROBE_MSH[element.kind],false,bases) : nothing
            (partition=false,coefval=scheme.coefval,expval=_view_probe_exponents(scheme.expval),
             coefgeo=geometry===nothing ? scheme.coefgeo : geometry.coefgeo,
             expgeo=geometry===nothing ? _view_probe_exponents(scheme.expgeo) : geometry.expgeo)
        end
        push!(records,_view_probe_record(element.kind,element.coords,
            element.values[:,:,step+1],basis,index,caller))
    end
    return records
end

function _view_probe_component_count(view,caller)
    if view.model_kind!==nothing
        context=get(view.model_data,(0,view.model_kind),nothing)
        return context===nothing ? (isempty(view.model_step_ranges) ? 0 :
            first(view.model_step_ranges).components) : context[3]
    end
    components=Int[]
    for type in view.list_types
        get(view.list_counts,type,0)>0 || continue
        layout=Post._pos_list_record_layout(view.interp,type)
        layout===nothing && throw(ArgumentError("$caller: unsupported list record $type"))
        push!(components,layout.ncomp)
    end
    return 1 in components ? 1 : 3 in components ? 3 : 9 in components ? 9 : 0
end

function _view_probe_coordinates(record,point)
    if record.cell===nothing
        return point==Tuple(record.coords[:,1]) ? (0.0,0.0,0.0) : nothing
    end
    cell=record.cell;scale=SizeField._pv_scheme_scale(cell)
    scale>0 || return nothing
    base=Tuple(cell.coords[:,1])
    query=SizeField._postview_scaled_query(point,base,scale)
    all(isfinite,query) || return nothing
    coordinates=SizeField._pv_scheme_inverse(cell,query,base,scale)
    coordinates===nothing && return nothing
    SizeField._pv_scheme_inside(cell.kind,coordinates...,
        MeshPointLocation.STRICT_REFERENCE_TOLERANCE) || return nothing
    return coordinates
end

@noinline function _view_probe_exact_gradient(jacobian,weights,record,scale)
    R=Rational{BigInt}
    matrix=R.(jacobian)
    reference=zeros(R,size(weights,1),size(record.values,1))
    for component in axes(reference,2),axis in axes(reference,1)
        anchor=record.partition ? R(record.values[component,1]) : R(0)
        for node in axes(weights,2)
            reference[axis,component]+=R(weights[axis,node])*(R(record.values[component,node])-anchor)
        end
    end
    result=matrix*((transpose(matrix)*matrix)\R.(reference))
    return vec(Float64.(result./R(scale)))
end

@inline function _view_probe_rescale(value,numerator,denominator)
    iszero(value) && return value
    a,ae=frexp(value);b,be=frexp(numerator);c,ce=frexp(denominator)
    return ldexp((a*b)/c,ae+be-ce)
end

function _view_probe_evaluate(record,coordinates,gradient,caller;geometry=record)
    nc=size(record.values,1)
    record.cell===nothing && return gradient ? zeros(3nc) : vec(copy(record.values[:,1]))
    cell=record.cell;u,v,w=coordinates
    weights=zeros(gradient ? record.dim : 1,size(record.values,2))
    for term in eachindex(cell.expval)
        monomials=gradient ? SizeField._pv_monomial_gradient(cell.expval[term],u,v,w,cell.pyramid) :
            (SizeField._pv_monomial(cell.expval[term],u,v,w,cell.pyramid),)
        for node in axes(record.values,2),axis in axes(weights,1)
            weights[axis,node]=muladd(cell.coefval[node,term],monomials[axis],weights[axis,node])
        end
    end
    all(isfinite,weights) || throw(ArgumentError("$caller: value basis is not finite"))
    reference=zeros(size(weights,1),nc)
    exact_gradient=gradient && any(1:nc) do component
        largest=maximum(abs,@view(record.values[component,:]))
        largest>sqrt(floatmax(Float64)) || 0<largest<sqrt(floatmin(Float64))
    end
    for component in 1:nc,axis in axes(weights,1)
        anchor=record.partition ? record.values[component,1] : 0.0
        value=record.partition && !gradient ? anchor : 0.0
        for node in axes(weights,2)
            delta=record.values[component,node]-anchor
            value=muladd(weights[axis,node],delta,value)
        end
        if !isfinite(value) && !gradient
            # Rare cancellation/overflow cases retain the same owned nodal
            # weights without overflowing intermediate products or differences.
            R=Rational{BigInt}
            exact=record.partition && !gradient ? R(anchor) : R(0)
            for node in axes(weights,2)
                exact+=R(weights[axis,node])*(R(record.values[component,node])-R(anchor))
            end
            value=Float64(exact)
        end
        reference[axis,component]=value
    end
    if gradient
        !exact_gradient && all(iszero,reference) && return zeros(3nc)
        gc=geometry.cell;scale=SizeField._pv_scheme_scale(gc)
        derivatives=SizeField._pv_scheme_jacobian(gc,u,v,w,Tuple(gc.coords[:,1]),scale)
        jacobian=[derivatives[column][axis] for axis in 1:3,column in 1:geometry.dim]
        decomposition=svd(jacobian;full=false)
        # J' g = dV/dξ; the thin SVD gives the minimum-norm tangential
        # derivative for curves/surfaces and the full derivative in 3D.
        result=if !exact_gradient && all(isfinite,reference) && all(s->isfinite(s) && s>0,decomposition.S) &&
                minimum(decomposition.S)>sqrt(eps(Float64))*maximum(decomposition.S)
            normalized=decomposition.U*((decomposition.Vt*reference)./decomposition.S)
            physical=[_view_probe_rescale(normalized[axis,component],1.0,scale)
                for axis in 1:3,component in 1:nc]
            all(isfinite,physical) ? vec(physical) : _view_probe_exact_gradient(jacobian,weights,record,scale)
        else
            _view_probe_exact_gradient(jacobian,weights,record,scale)
        end
    else
        result=vec(reference)
    end
    all(isfinite,result) || throw(ArgumentError("$caller: interpolated result is not finite"))
    return result
end

function _view_probe_find(records,point,dimension,element_coords;
        reverse_order=false,exact_match=true,nearest_match=true)
    selected=collect(eachindex(records))
    if reverse_order
        sort!(selected;by=i->(-get(_PLUGIN_SEARCH_RANK,records[i].kind,0),-i))
    end
    nearest=nothing;distance=Inf
    for index in selected
        record=records[index]
        dimension>=0 && record.dim!=dimension && continue
        element_coords===nothing || record.coords==element_coords || continue
        coordinates=exact_match ? _view_probe_coordinates(record,point) : nothing
        coordinates===nothing || return (record,coordinates),0.0
        nearest_match || continue
        for node in axes(record.coords,2)
            candidate=hypot((point.-Tuple(record.coords[:,node]))...)
            closer=nearest===nothing || candidate<distance
            if !closer && isinf(candidate) && isinf(distance)
                R=Rational{BigInt}
                trial=sum((R(point[axis])-R(record.coords[axis,node]))^2 for axis in 1:3)
                previous=sum((R(point[axis])-R(nearest[2][axis]))^2 for axis in 1:3)
                closer=trial<previous
            end
            if closer
                distance=candidate;nearest=(record,Tuple(record.coords[:,node]))
            end
        end
    end
    return nearest,distance
end

function _view_probe(model_name,tag,x,y,z)
    caller="API.view.probe"
    model_name isa AbstractString || throw(ArgumentError("$caller: modelName must be a string"))
    point=ntuple(i->_mesh_query_coordinate((x,y,z)[i],caller,"coordinate"),3)
    return lock(STATE_LOCK) do
        _slot_capture_locked!();view=_session_view_data(tag,caller)
        slot=_native_view_model_slot(String(model_name))
        slot===nothing && throw(ArgumentError("$caller: Unknown model '$model_name'"))
        context=get(view.model_contexts,(0,view.model_kind),nothing)
        if context!==nothing && context.model_slot===slot
            slot.mesh===nothing && throw(ArgumentError("$caller: probe requires a mesh"))
            meshes=IdDict{Any,Any}();bases=Dict{Any,Any}()
            records=_view_probe_model_records(view,0,meshes,bases,caller)
            match,distance=_view_probe_find(records,point,-1,nothing;nearest_match=false)
            distance==0.0 && match!==nothing && return _view_probe_evaluate(match...,false,caller)
            geometry=_view_probe_model_records(view,0,meshes,bases,caller;geometry_only=true)
            _,geometric_distance=_view_probe_find(geometry,point,-1,nothing;nearest_match=false)
            geometric_distance==0.0 || throw(ArgumentError("$caller: no element contains ($x, $y, $z)"))
        end
        message="$caller: view $tag holds no probe-able model data (step 0 ElementData/ElementNodeData/NodeData)"
        NATIVE_LOOKUP_QUERY_ERRORS[]+=1;NATIVE_MSG_LAST_ERROR[]=message
        _native_diagnostic_emit_locked(:error,message)
        throw(ArgumentError(message))
    end
end

function _view_probe_native(tag,x,y,z,step=-1,num_comp=-1,gradient=false,
        distance_max=0.0,x_element=Float64[],y_element=Float64[],z_element=Float64[],dimension=-1)
    caller="API.view.probe"
    point=ntuple(i->_mesh_query_coordinate((x,y,z)[i],caller,"coordinate"),3)
    requested=Int(_mesh_query_integer(step,caller,"step"))
    wanted=Int(_mesh_query_integer(num_comp,caller,"numComp"))
    dim=Int(_mesh_query_integer(dimension,caller,"dim"))
    dim<=3 || throw(ArgumentError("$caller: dim must be negative or in 0:3"))
    grad=_mesh_query_bool(gradient,caller,"gradient")
    radius=_mesh_query_coordinate(distance_max,caller,"distanceMax")
    vectors=(Model._model_evaluation_values(x_element,1,caller,"xElemCoord"),
        Model._model_evaluation_values(y_element,1,caller,"yElemCoord"),
        Model._model_evaluation_values(z_element,1,caller,"zElemCoord"))
    all(length(v)==length(vectors[1]) for v in vectors) || throw(ArgumentError(
        "$caller: element coordinate arrays must have equal lengths"))
    element_coords=isempty(vectors[1]) ? nothing : permutedims(hcat(vectors...))
    return lock(STATE_LOCK) do
        _slot_capture_locked!();view=_session_view_data(tag,caller)
        steps=_view_data_steps(view)
        steps==0 && return Float64[],-1.0
        requested<steps || throw(ArgumentError("$caller: step $requested is outside 0:$(steps-1)"))
        model=view.model_kind!==nothing
        source_components=_view_probe_component_count(view,caller)
        wanted>0 && !(wanted in (1,3,9)) && return Float64[],-1.0
        classes=wanted>0 ? (wanted,) : (1,3,9)
        meshes=IdDict{Any,Any}();bases=Dict{Any,Any}()
        list_elements=model ? nothing : _native_list_view_records(view,caller)
        function records_at(index)
            return model ? _view_probe_model_records(view,index,meshes,bases,caller;
                components=source_components,include_missing=true) :
                _view_probe_list_records(view,index,bases,caller;elements=list_elements)
        end
        function exact_at(records,query)
            for components in classes
                selected=[record for record in records if size(record.values,1)==components]
                candidate_match,candidate_distance=_view_probe_find(selected,query,dim,element_coords;
                    reverse_order=true,nearest_match=false)
                candidate_distance==0.0 && candidate_match!==nothing && return candidate_match
            end
            return nothing
        end
        initial=requested<0 ? 0 : requested
        records=records_at(initial)
        match=exact_at(records,point)
        distance=0.0
        if match===nothing
            # Explicit class-only searches do not request a nearest distance
            # until a nonzero radius enables the common nearest-node search.
            wanted>0 && radius==0 && return Float64[],-1.0
            nearest,distance=_view_probe_find(records,point,dim,element_coords;
                reverse_order=true,exact_match=false)
            nearest===nothing && return Float64[],-1.0
            (radius<0 || (radius>0 && distance<=radius)) || return Float64[],distance
            match=exact_at(records,nearest[2])
            match===nothing && return Float64[],distance
        end
        record,coordinates=match
        nc=size(record.values,1)
        if requested>=0
            return _view_probe_evaluate(record,coordinates,grad,caller),distance
        end
        width=Base.checked_mul(nc,grad ? 3 : 1)
        output=Vector{Float64}(undef,Base.checked_mul(steps,width))
        carried=model ? copy(record.values) : nothing
        for index in 0:steps-1
            values=if index==0
                _view_probe_evaluate(record,coordinates,grad,caller)
            elseif model
                if _view_probe_material_values!(carried,record,view,index)
                    _view_probe_evaluate(_view_probe_material_record(record,carried),
                        coordinates,grad,caller;geometry=record)
                else
                    zeros(width)
                end
            else
                current_records=records_at(index)
                current=findfirst(r->r.token==record.token && size(r.values,1)==nc,current_records)
                current===nothing ? zeros(width) : _view_probe_evaluate(
                    current_records[current],coordinates,grad,caller;geometry=record)
            end
            output[index*width+1:(index+1)*width].=values
        end
        return output,distance
    end
end
