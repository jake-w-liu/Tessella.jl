# View queries retain the complete value and geometry bases. The shared
# safeguarded scheme inverse supplies intrinsic coordinates and certifies the
# actual physical map; derivatives are pulled back with the geometry Jacobian.
using LinearAlgebra: svd, SingularException

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
        basis.coefgeo,basis.expgeo,coordinates,values;native_geometry=basis.geometry_partition)
    SizeField._pv_scheme_validate(cell,1,caller)
    return (;kind,dim=SizeField._pv_scheme_dim(cell.kind),coords=coordinates,
        values,cell,partition=basis.partition,token,identity)
end

function _view_probe_basis(msh,constant,bases)
    return get!(bases,(msh,constant)) do
        native=_native_view_writer_basis(msh)
        (partition=true,geometry_partition=true,coefval=constant ? ones(1,1) : native.coefficients,
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
        cell=SizeField._pv_scheme_with_values(cell,values)
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
            (partition=false,geometry_partition=geometry!==nothing || scheme===Post._pos_order2_scheme(element.kind),coefval=scheme.coefval,expval=_view_probe_exponents(scheme.expval),
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

const _VIEW_PROBE_RECOVERY_PRECISION=16384

mutable struct _ViewProbeLocation
    reference::NTuple{3,Float64}
    point::NTuple{3,Float64}
    precise::Union{Nothing,NTuple{3,BigFloat}}
end

_view_probe_cell_base(cell)=SizeField._pv_scheme_base(cell)

# A bounded mantissa keeps extreme integer powers from constructing enormous
# rational denominators. The do-blocks use Julia's task-local precision and
# rounding scopes; they do not change the caller's arithmetic environment.
function _view_probe_recovery(f)
    return setprecision(BigFloat,_VIEW_PROBE_RECOVERY_PRECISION) do
        setrounding(BigFloat,RoundNearest) do
            f()
        end
    end
end

function _view_probe_precise_certificate(cell,coordinates,query,base,scale)
    mapped=SizeField._pv_scheme_map(cell,coordinates...,base,scale)
    all(isfinite,mapped) || return nothing
    residual=ntuple(axis->query[axis]-mapped[axis],3)
    derivatives=SizeField._pv_scheme_jacobian(cell,coordinates...,base,scale)
    all(column->all(isfinite,column),derivatives) || return nothing
    dim=SizeField._pv_scheme_dim(cell.kind)
    reference=collect(coordinates)
    if dim<=2
        for axis in (dim==1 ? (2,3) : (3,))
            normal=derivatives[axis]
            denominator=sum(abs2,normal)
            denominator>0 || return nothing
            reference[axis]=sum(residual[i]*normal[i] for i in 1:3)/denominator
            residual=ntuple(i->residual[i]-reference[axis]*normal[i],3)
        end
    end
    all(isfinite,reference) || return nothing
    tolerance=max(BigFloat(1e-10),BigFloat(256eps(Float64))*max(BigFloat(1),maximum(abs,query)))
    hypot(residual...)<=tolerance || return nothing
    result=Tuple(reference)
    SizeField._pv_scheme_inside(cell.kind,result...,
        MeshPointLocation.STRICT_REFERENCE_TOLERANCE) || return nothing
    return result
end

function _view_probe_precise_coordinates!(location,geometry)
    location.precise===nothing || return location.precise
    return _view_probe_recovery() do
        cell=geometry.cell;scale=SizeField._pv_scheme_scale(cell)
        base=map(BigFloat,_view_probe_cell_base(cell))
        query=ntuple(axis->(BigFloat(location.point[axis])-base[axis])/BigFloat(scale),3)
        dim=geometry.dim
        coordinates=ntuple(axis->axis<=dim ? BigFloat(location.reference[axis]) : BigFloat(0),3)
        for _ in 1:64
            mapped=SizeField._pv_scheme_map(cell,coordinates...,base,scale)
            all(isfinite,mapped) || return nothing
            residual=BigFloat[query[axis]-mapped[axis] for axis in 1:3]
            derivatives=SizeField._pv_scheme_jacobian(cell,coordinates...,base,scale)
            all(column->all(isfinite,column),derivatives) || return nothing
            jacobian=BigFloat[derivatives[column][axis] for axis in 1:3,column in 1:dim]
            step=try
                (transpose(jacobian)*jacobian)\(transpose(jacobian)*residual)
            catch error
                error isa SingularException || rethrow()
                return nothing
            end
            all(isfinite,step) || return nothing
            limit=16eps(BigFloat)*max(BigFloat(1),maximum(abs,coordinates))
            if maximum(abs,step)<=limit
                certified=_view_probe_precise_certificate(cell,coordinates,query,base,scale)
                certified===nothing && return nothing
                location.precise=certified
                return certified
            end
            previous=sum(abs2,residual);alpha=BigFloat(1);advanced=false
            for _ in 1:24
                candidate=ntuple(axis->axis<=dim ? coordinates[axis]+alpha*step[axis] : BigFloat(0),3)
                trial=SizeField._pv_scheme_map(cell,candidate...,base,scale)
                if all(isfinite,trial) && sum(abs2,query.-trial)<=previous
                    coordinates=candidate;advanced=true;break
                end
                alpha/=2
            end
            advanced || return nothing
        end
        return nothing
    end
end

function _view_probe_coordinates(record,point)
    if record.cell===nothing
        return point==Tuple(record.coords[:,1]) ? (0.0,0.0,0.0) : nothing
    end
    cell=record.cell;scale=SizeField._pv_scheme_scale(cell)
    scale>0 || return nothing
    base=_view_probe_cell_base(cell)
    query=SizeField._postview_scaled_query(point,base,scale)
    all(isfinite,query) || return nothing
    coordinates=SizeField._pv_scheme_inverse(cell,query,base,scale)
    coordinates===nothing && return nothing
    SizeField._pv_scheme_inside(cell.kind,coordinates...,
        MeshPointLocation.STRICT_REFERENCE_TOLERANCE) || return nothing
    location=_ViewProbeLocation(coordinates,point,nothing)
    # Reconstructing the physical query catches a small endpoint offset lost
    # by normalization about a far-away first node.
    if any(axis->muladd(query[axis],scale,base[axis])!=point[axis],1:3)
        _view_probe_precise_coordinates!(location,record)===nothing && return nothing
    end
    return location
end

function _view_probe_raw_monomials(exponents,coordinates,pyramid,gradient)
    a,b,c=map(BigFloat,exponents);u,v,w=coordinates
    value=u^a*v^b*w^c
    ee=pyramid ? -min(a,b) : BigFloat(0)
    t=pyramid ? max(BigFloat(1e-14),1-w) : BigFloat(1)
    power=t^ee
    gradient || return (value*power,)
    du=iszero(a) ? BigFloat(0) : a*u^(a-1)*v^b*w^c*power
    dv=iszero(b) ? BigFloat(0) : b*u^a*v^(b-1)*w^c*power
    dw=u^a*v^b*((iszero(c) ? BigFloat(0) : c*w^(c-1))*power+
        (iszero(ee) ? BigFloat(0) : -ee*w^c*t^(ee-1)))
    return du,dv,dw
end

@noinline function _view_probe_raw_evaluate(record,location,gradient,caller,geometry)
    return _view_probe_recovery() do
        coordinates=_view_probe_precise_coordinates!(location,geometry)
        coordinates===nothing && throw(ArgumentError("$caller: could not certify the original physical query"))
        cell=record.cell;nc=size(record.values,1)
        # Products and sums of the stored binary coefficients/data are formed
        # before powers. Duplicate monomials can cancel without evaluating an
        # otherwise overflowing or underflowing power first.
        groups=Dict{NTuple{3,Float64},Vector{BigFloat}}()
        order=NTuple{3,Float64}[]
        for term in eachindex(cell.expval)
            exponent=map(x->iszero(x) ? 0.0 : x,cell.expval[term])
            if !haskey(groups,exponent)
                groups[exponent]=zeros(BigFloat,nc);push!(order,exponent)
            end
            coefficient=groups[exponent]
            for component in 1:nc
                anchor=record.partition ? BigFloat(record.values[component,1]) : BigFloat(0)
                for node in axes(record.values,2)
                    coefficient[component]+=BigFloat(cell.coefval[node,term])*
                        (BigFloat(record.values[component,node])-anchor)
                end
            end
        end
        reference=zeros(BigFloat,gradient ? record.dim : 1,nc)
        if !gradient && record.partition
            reference[1,:].=BigFloat.(record.values[:,1])
        end
        for exponent in order
            coefficient=groups[exponent]
            all(iszero,coefficient) && continue
            monomials=_view_probe_raw_monomials(exponent,coordinates,cell.pyramid,gradient)
            for component in 1:nc,axis in axes(reference,1)
                iszero(coefficient[component]) && continue
                reference[axis,component]+=coefficient[component]*monomials[axis]
            end
        end
        result=if gradient
            gc=geometry.cell;scale=SizeField._pv_scheme_scale(gc)
            base=map(BigFloat,_view_probe_cell_base(gc))
            derivatives=SizeField._pv_scheme_jacobian(gc,coordinates...,base,scale)
            jacobian=BigFloat[derivatives[column][axis] for axis in 1:3,column in 1:geometry.dim]
            jacobian*((transpose(jacobian)*jacobian)\reference)/BigFloat(scale)
        else
            reference
        end
        values=vec(Float64.(result))
        all(isfinite,values) || throw(ArgumentError("$caller: interpolated result is not finite"))
        return values
    end
end

@inline function _view_probe_nonzero_monomial(exponents,coordinates,axis,pyramid)
    axis==0 && return all(i->coordinates[i]!=0 || exponents[i]==0,1:3)
    if axis<=2 || !pyramid
        exponents[axis]!=0 || return false
        return all(i->coordinates[i]!=0 || exponents[i]-(i==axis)==0,1:3)
    end
    # The two w-derivative terms can cancel. A conservative uncertainty flag
    # is enough; the recovery evaluates their actual signed sum.
    a,b,c=exponents;u,v,w=coordinates
    (u!=0 || a==0) && (v!=0 || b==0) || return false
    return (c!=0 && (w!=0 || c==1)) || (min(a,b)!=0 && (w!=0 || c==0))
end

@inline function _view_probe_product_lost(a,b,product)
    return !isfinite(product) || (0<abs(product)<floatmin(Float64)) ||
        (product==0 && a!=0 && b!=0)
end

# A finite remainder can lose useful digits well before its error enclosure
# crosses zero. Recover when the conservative arithmetic budget exceeds
# 1e-12 of the result. Division avoids overflowing gamma*mass or underflowing
# 1e-12*abs(value) while classifying large or tiny finite values.
@inline function _view_probe_roundoff_excess(value,mass,gamma)
    return !isfinite(value) || !isfinite(mass) ||
        (mass>0 && abs(value)/mass<=gamma/1e-12)
end

@inline function _view_probe_rescale(value,numerator,denominator)
    iszero(value) && return value
    a,ae=frexp(value);b,be=frexp(numerator);c,ce=frexp(denominator)
    return ldexp((a*b)/c,ae+be-ce)
end

function _view_probe_evaluate(record,coordinates,gradient,caller;geometry=record)
    nc=size(record.values,1)
    record.cell===nothing && return gradient ? zeros(3nc) : vec(copy(record.values[:,1]))
    if record.partition && all(component->all(==(record.values[component,1]),
            @view(record.values[component,:])),1:nc)
        return gradient ? zeros(3nc) : vec(copy(record.values[:,1]))
    end
    cell=record.cell;u,v,w=coordinates.reference
    weights=zeros(gradient ? record.dim : 1,size(record.values,2))
    masses=zeros(size(weights))
    recover=coordinates.precise!==nothing || any(record.values) do value
        magnitude=abs(value)
        magnitude>sqrt(floatmax(Float64)) || 0<magnitude<sqrt(floatmin(Float64))
    end
    operations=32.0+8.0*length(cell.expval)
    for exponent in cell.expval,degree in exponent
        operations+=16abs(degree)
    end
    mu=operations*eps(Float64)
    gamma=mu<0.5 ? mu/(1-mu) : Inf
    recover |= mu>1e-12
    for term in eachindex(cell.expval)
        monomials=gradient ? SizeField._pv_monomial_gradient(cell.expval[term],u,v,w,cell.pyramid) :
            (SizeField._pv_monomial(cell.expval[term],u,v,w,cell.pyramid),)
        for node in axes(record.values,2),axis in axes(weights,1)
            monomial=monomials[axis];coefficient=cell.coefval[node,term]
            recover |= !isfinite(monomial) || (0<abs(monomial)<floatmin(Float64)) ||
                (monomial==0 && _view_probe_nonzero_monomial(cell.expval[term],
                    coordinates.reference,gradient ? axis : 0,cell.pyramid))
            product=coefficient*monomial
            recover |= _view_probe_product_lost(coefficient,monomial,product)
            weights[axis,node]=muladd(coefficient,monomial,weights[axis,node])
            masses[axis,node]+=abs(product)
        end
    end
    for index in eachindex(weights)
        recover |= _view_probe_roundoff_excess(weights[index],masses[index],gamma)
    end
    reference=zeros(size(weights,1),nc)
    for component in 1:nc,axis in axes(weights,1)
        anchor=record.partition ? record.values[component,1] : 0.0
        value=record.partition && !gradient ? anchor : 0.0
        mass=abs(value)
        for node in axes(weights,2)
            delta=record.values[component,node]-anchor
            product=weights[axis,node]*delta
            recover |= !isfinite(delta) || _view_probe_product_lost(weights[axis,node],delta,product)
            value=muladd(weights[axis,node],delta,value)
            # Carry the basis-product mass through the data fold: a small
            # rounded weight can hide a large error amplified by its datum.
            mass+=masses[axis,node]*abs(delta)
        end
        recover |= _view_probe_roundoff_excess(value,mass,gamma)
        reference[axis,component]=value
    end
    recover && return _view_probe_raw_evaluate(record,coordinates,gradient,caller,geometry)
    if gradient
        all(iszero,reference) && return zeros(3nc)
        gc=geometry.cell;scale=SizeField._pv_scheme_scale(gc)
        derivatives=SizeField._pv_scheme_jacobian(gc,u,v,w,_view_probe_cell_base(gc),scale)
        jacobian=[derivatives[column][axis] for axis in 1:3,column in 1:geometry.dim]
        decomposition=svd(jacobian;full=false)
        # J' g = dV/dξ; the thin SVD gives the minimum-norm tangential
        # derivative for curves/surfaces and the full derivative in 3D.
        result=if all(s->isfinite(s) && s>0,decomposition.S) &&
                minimum(decomposition.S)>sqrt(eps(Float64))*maximum(decomposition.S)
            normalized=decomposition.U*((decomposition.Vt*reference)./decomposition.S)
            physical=[_view_probe_rescale(normalized[axis,component],1.0,scale)
                for axis in 1:3,component in 1:nc]
            all(isfinite,physical) ? vec(physical) : _view_probe_raw_evaluate(record,coordinates,true,caller,geometry)
        else
            _view_probe_raw_evaluate(record,coordinates,true,caller,geometry)
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
