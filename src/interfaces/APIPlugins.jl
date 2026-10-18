# ---- session plugins ----------------------------------------------------------
#
# Gmsh 4.15.2 `plugin` module support. `_PLUGINS` registers each plugin name
# with its numeric and string option tables (upstream defaults plus the
# automatic "Run" numeric option, so `set_number(name,"Run",1)` executes like
# `Plugin(Name).Run`) and the function performing real work on session views.
# `SESSION_PLUGINS` holds the per-name configured option values; it is emptied
# next to `SESSION_VIEWS` in initialize/finalize/clear so a fresh session
# restores upstream defaults. All state reads/mutations run under `STATE_LOCK`.
#
# Implemented plugins — all compute deterministic output on session view data;
# "in-place" plugins rewrite the target view like upstream, the others create
# new views allocated through `_fresh_view_tag`:
#
#   Integrate       spatial field integral into one `*_Integrate` SP view
#                   (OverTime >= 0 integrates nodal values over physical time),
#   MinMax          per-step or over-time min/max into `*_Min`/`*_Max` SP views,
#   Probe           interpolated value at (X,Y,Z) into one `*_Probe` point view,
#   Skin            boundary facets of list records into `*_Skin`,
#   ExtractElements value/dimension filtered elements into `*_ExtractElements`,
#   Transform       in-place affine map of list node coordinates,
#   Remove          in-place element-class/string removal,
#   Smooth          in-place coincident-node value averaging,
#   NewView         constant `NodeData` model view over the current mesh.
#
# Session views carry no per-element visibility flags, so the upstream `Visible`
# options are accepted and stored but currently never hide an element.
# `Skin` only covers list-based views (`FromMesh = 1`, which would publish
# discrete model boundary entities, rejects) and `Transform`/`Remove`/`Smooth`
# reject model-based views — upstream in-place writers mutate `PViewDataList`.

mutable struct _PluginSpec
    numbers::Dict{String,Float64}
    strings::Dict{String,String}
end

const SESSION_PLUGINS = Dict{String,_PluginSpec}()

function _plugin_check_name(name,caller::AbstractString)
    name isa AbstractString || throw(ArgumentError("$caller: plugin name must be a string"))
    return String(name)
end

function _plugin_check_option(option,caller::AbstractString)
    option isa AbstractString || throw(ArgumentError("$caller: option name must be a string"))
    return String(option)
end

function _plugin_entry(name::String,caller::AbstractString)
    entry=get(_PLUGINS,name,nothing)
    entry===nothing && throw(ArgumentError("$caller: unknown plugin '$name'"))
    return entry
end

# Lazily materialize the per-name option copy so unconfigured plugins run with
# upstream defaults.
function _plugin_state_locked(name::String)
    return get!(SESSION_PLUGINS,name) do
        entry=_PLUGINS[name]
        _PluginSpec(copy(entry.number_options),copy(entry.string_options))
    end
end

function _plugin_set_number(name,option,value)
    caller="API.plugin.set_number"
    pname=_plugin_check_name(name,caller)
    opt=_plugin_check_option(option,caller)
    value isa Bool && throw(ArgumentError("$caller: value must not be Bool"))
    value isa Real || throw(ArgumentError("$caller: value must be a real number"))
    numeric=Float64(value)
    isfinite(numeric) || throw(ArgumentError("$caller: value must be finite"))
    entry=_plugin_entry(pname,caller)
    haskey(entry.number_options,opt) || throw(ArgumentError(
        "$caller: unknown numeric option '$opt' for plugin '$pname'"))
    return lock(STATE_LOCK) do
        spec=_plugin_state_locked(pname)
        spec.numbers[opt]=numeric
        if opt=="Run" && !iszero(numeric)
            entry.run_fn(spec,caller)
        end
        return nothing
    end
end

function _plugin_set_string(name,option,value)
    caller="API.plugin.set_string"
    pname=_plugin_check_name(name,caller)
    opt=_plugin_check_option(option,caller)
    value isa AbstractString || throw(ArgumentError("$caller: value must be a string"))
    entry=_plugin_entry(pname,caller)
    haskey(entry.string_options,opt) || throw(ArgumentError(
        "$caller: unknown string option '$opt' for plugin '$pname'"))
    return lock(STATE_LOCK) do
        spec=_plugin_state_locked(pname)
        spec.strings[opt]=String(value)
        return nothing
    end
end

function _plugin_run(name)
    caller="API.plugin.run"
    pname=_plugin_check_name(name,caller)
    entry=_plugin_entry(pname,caller)
    return lock(STATE_LOCK) do
        spec=_plugin_state_locked(pname)
        entry.run_fn(spec,caller)
        return nothing
    end
end

# `View` < 0 selects the most recently allocated session view, matching the
# upstream "current view" slot passed to `getView(iView, v)`.
function _plugin_view(numbers::Dict{String,Float64},caller::AbstractString)
    raw=numbers["View"]
    if raw<0
        isempty(SESSION_VIEWS) && throw(ArgumentError(
            "$caller: no views exist to run the plugin on"))
        tag=maximum(keys(SESSION_VIEWS))
    else
        raw==floor(raw) || throw(ArgumentError("$caller: View tag must be an integer"))
        tag=Int(raw)
        haskey(SESSION_VIEWS,tag) || throw(ArgumentError(
            "$caller: unknown view tag $tag"))
    end
    return tag,_session_view_data(tag,caller)
end

const _PLUGIN_KIND_DIM=(point=0,line=1,triangle=2,quadrangle=2,tetrahedron=3,
                       hexahedron=3,prism=3,pyramid=3)
const _PLUGIN_KIND_SHAPE=(point="P",line="L",triangle="T",quadrangle="Q",
                         tetrahedron="S",hexahedron="H",prism="I",pyramid="Y")
const _PLUGIN_COMP_HEAD=Dict{Int,String}(1=>"S",3=>"V",9=>"T")
# `msh_spec` family → the `_POS_RECORD_TYPES` element kind.
const _PLUGIN_FAMILY_KIND=(pnt=:point,lin=:line,tri=:triangle,qua=:quadrangle,
                           tet=:tetrahedron,hex=:hexahedron,pri=:prism,pyr=:pyramid)

# Upstream `ComputeScalarRep` with the default `tensorRep = 0`: scalars pass
# through, 3-component rows reduce to the Euclidean norm and 9-component rows
# to the Von-Mises invariant; other widths contribute 0.
function _plugin_scalar_rep(vals,base::Int,ncomp::Int)
    ncomp==1 && return vals[base+1]
    if ncomp==3
        return hypot(vals[base+1],vals[base+2],vals[base+3])
    elseif ncomp==9
        tr=(vals[base+1]+vals[base+5]+vals[base+9])/3
        s=0.0
        for k in 0:8
            d=vals[base+1+k]-(k in (0,4,8) ? tr : 0.0)
            s=muladd(d,d,s)
        end
        return sqrt(1.5*s)
    end
    return 0.0
end

# One geometric element resolved for plugin scans. `geo` holds the primary-node
# coordinates used for measures and point location; `vcoords` marks the value
# sites scanned for min/max and interpolation (the corner set for linear data).
# `vals` packs `val_arity*ncomp*nsteps` values — per step, per node the `ncomp`
# components — the `_POS_RECORD_TYPES` payload convention. `etag` is the source
# entity tag for model-based data (0 for list data).
struct _PluginElement
    kind::Symbol
    dim::Int
    geo::Matrix{Float64}
    vcoords::Matrix{Float64}
    vals::Vector{Float64}
    ncomp::Int
    etag::UInt64
end

function _plugin_list_steps(view,caller::AbstractString)
    steps=_view_plain_list_steps(view.list_types,view.list_counts,view.list_data)
    steps===nothing && throw(ArgumentError(
        "$caller: view contains malformed list element payload"))
    return steps
end

# Decode every `_POS_RECORD_TYPES` block of a list-based view into independent
# element records (coordinates use the Gmsh x-block/y-block/z-block layout).
function _plugin_list_elements(view,caller::AbstractString)
    steps=_plugin_list_steps(view,caller)
    elements=_PluginElement[]
    for type in view.list_types
        layout=get(Post._POS_RECORD_TYPES,type,nothing)
        layout===nothing && throw(ArgumentError(
            "$caller: unsupported list data type $type"))
        kind,arity,ncomp=layout
        count=view.list_counts[type];data=view.list_data[type]
        width=3arity+Base.checked_mul(arity,Base.checked_mul(ncomp,steps))
        length(data)==Base.checked_mul(count,width) || throw(ArgumentError(
            "$caller: list data payload does not align with numEle"))
        for element in 0:count-1
            base=element*width
            geo=Matrix{Float64}(undef,3,arity)
            for axis in 1:3,node in 1:arity
                geo[axis,node]=data[base+(axis-1)*arity+node]
            end
            vals=Vector{Float64}(undef,arity*ncomp*steps)
            copyto!(vals,1,data,base+3arity+1,length(vals))
            push!(elements,_PluginElement(kind,_PLUGIN_KIND_DIM[kind],
                                          geo,copy(geo),vals,ncomp,UInt64(0)))
        end
    end
    return elements,steps
end

# Populated `(step, data_type)` keys in step order for a model-based view.
function _plugin_model_keys(view)
    view.model_kind===nothing && return Tuple{Int,String}[]
    kind=view.model_kind
    return sort!(Tuple{Int,String}[key for key in keys(view.model_data)
                 if key[2]==kind && haskey(view.model_contexts,key)];by=first)
end

# Public element tag (or dense tag when no public table exists) → dense index.
function _plugin_element_index(mesh,public,tag)
    if public!==nothing
        index=get(public.element_indices,UInt64(tag),0)
        index!=0 && return index
    end
    _,_,total=_mesh_element_offsets(mesh)
    1<=tag<=total || return 0
    return Int(tag)
end

function _plugin_model_mesh(slot,caller::AbstractString)
    slot===nothing && throw(ArgumentError("$caller: model view has no model"))
    mesh=slot.mesh
    mesh===nothing && throw(ArgumentError("$caller: model view has no mesh"))
    return mesh
end

# Element records for one populated `(step, data_type)` model-data block; each
# element carries its single-step values: ElementData rows broadcast the stored
# row to every scanned node like the upstream GModel data path, ElementNodeData
# rows keep their node-major payload (variable row widths come from the step
# context `counts`), and NodeData rows gather per element node tag (elements
# with incomplete node coverage are skipped). `etag` retains the source tag.
function _plugin_model_elements(view,key,caller::AbstractString)
    ctx=view.model_contexts[key];entity_tags,values,nc=view.model_data[key]
    nc>0 || throw(ArgumentError("$caller: malformed model data component count"))
    mesh=_plugin_model_mesh(ctx.model_slot,caller)
    public=_mesh_public_tags(mesh)
    kind=view.model_kind
    counts=get(ctx,:counts,Int32[])
    offsets=Vector{Int}(undef,length(entity_tags)+1);offsets[1]=1
    for row in eachindex(entity_tags)
        width=row<=length(counts) ? Int(counts[row]) :
            length(values)÷length(entity_tags)
        offsets[row+1]=Base.checked_add(offsets[row],width)
    end
    offsets[end]-1==length(values) || throw(ArgumentError(
        "$caller: malformed model data row extents"))
    elements=_PluginElement[]
    if kind=="ElementData" || kind=="ElementNodeData"
        for (row,tag) in enumerate(entity_tags)
            index=_plugin_element_index(mesh,public,tag)
            index==0 && continue
            rec=mesh_element_record(mesh,index)
            spec=msh_spec(Int(rec.element_type))
            spec.family===:trih && continue
            ga=min(Elements._msh_primary_nodes(spec.family),length(rec.node_tags))
            geo=Matrix{Float64}(undef,3,ga)
            for axis in 1:3,node in 1:ga
                geo[axis,node]=mesh.coords[axis,Int(rec.node_tags[node])]
            end
            base=offsets[row]
            if kind=="ElementData"
                vals=Vector{Float64}(undef,ga*nc)
                for node in 1:ga,comp in 1:nc
                    vals[(node-1)*nc+comp]=values[base+comp-1]
                end
                push!(elements,_PluginElement(_PLUGIN_FAMILY_KIND[spec.family],
                    Int(spec.dim),geo,geo,vals,nc,UInt64(tag)))
            else
                nn=(offsets[row+1]-base)÷nc
                va=min(nn,length(rec.node_tags))
                vcoords=Matrix{Float64}(undef,3,va)
                for axis in 1:3,node in 1:va
                    vcoords[axis,node]=mesh.coords[axis,Int(rec.node_tags[node])]
                end
                vals=collect(Float64,@view(values[base:base+va*nc-1]))
                push!(elements,_PluginElement(_PLUGIN_FAMILY_KIND[spec.family],
                    Int(spec.dim),geo,vcoords,vals,nc,UInt64(tag)))
            end
        end
    elseif kind=="NodeData"
        row_of=Dict{UInt64,Int}(tag=>row for (row,tag) in enumerate(entity_tags))
        _,_,total=_mesh_element_offsets(mesh)
        for index in 1:total
            rec=mesh_element_record(mesh,index)
            spec=msh_spec(Int(rec.element_type))
            spec.family===:trih && continue
            nn=length(rec.node_tags)
            rows=Vector{Int}(undef,nn)
            covered=true
            for node in 1:nn
                tag=_cache_node_tag(mesh,Int(rec.node_tags[node]))
                row=get(row_of,tag,0)
                row==0 && (covered=false;break)
                rows[node]=row
            end
            covered || continue
            ga=min(Elements._msh_primary_nodes(spec.family),nn)
            geo=Matrix{Float64}(undef,3,ga)
            for axis in 1:3,node in 1:ga
                geo[axis,node]=mesh.coords[axis,Int(rec.node_tags[node])]
            end
            vcoords=Matrix{Float64}(undef,3,nn)
            for axis in 1:3,node in 1:nn
                vcoords[axis,node]=mesh.coords[axis,Int(rec.node_tags[node])]
            end
            vals=Vector{Float64}(undef,nn*nc)
            for node in 1:nn,comp in 1:nc
                vals[(node-1)*nc+comp]=values[offsets[rows[node]]+comp-1]
            end
            etag=public===nothing ? UInt64(index) :
                get(public.element_tags,index,UInt64(index))
            push!(elements,_PluginElement(_PLUGIN_FAMILY_KIND[spec.family],
                Int(spec.dim),geo,vcoords,vals,nc,etag))
        end
    else
        throw(ArgumentError("$caller: unsupported model data type '$kind'"))
    end
    return elements,nc
end

# Physical time attached to step `s` (0-based), reusing the session view's own
# time bookkeeping (context times for model views, `list_time` for list data).
_plugin_step_time(view,s::Int)=_view_time_at(view,s)

# ---- element measures ---------------------------------------------------------

@inline _plugin_sub3(a,b)=(a[1]-b[1],a[2]-b[2],a[3]-b[3])
@inline _plugin_dot3(a,b)=muladd(a[1],b[1],muladd(a[2],b[2],a[3]*b[3]))
@inline _plugin_norm3(a)=hypot(a[1],a[2],a[3])
@inline function _plugin_cross3(a,b)
    return (a[2]*b[3]-a[3]*b[2],a[3]*b[1]-a[1]*b[3],a[1]*b[2]-a[2]*b[1])
end
@inline function _plugin_tri_area_vec(a,b,c)
    e1=_plugin_sub3(b,a);e2=_plugin_sub3(c,a)
    n=_plugin_cross3(e1,e2)
    return (0.5n[1],0.5n[2],0.5n[3])
end
@inline _plugin_col(coords,i)=(coords[1,i],coords[2,i],coords[3,i])

# Signed tetrahedral volume of (a,b,c,d): det/6.
@inline function _plugin_tet_volume(a,b,c,d)
    e1=_plugin_sub3(b,a);e2=_plugin_sub3(c,a);e3=_plugin_sub3(d,a)
    det=e1[1]*(e2[2]*e3[3]-e2[3]*e3[2])-
        e1[2]*(e2[1]*e3[3]-e2[3]*e3[1])+
        e1[3]*(e2[1]*e3[2]-e2[2]*e3[1])
    return det/6
end

# Face fans about the corner barycentre give an exact signed polyhedral volume
# for hexes, prisms and pyramids without family-specific tetrahedral splits.
const _PLUGIN_POLY_FACES=(hexahedron=((1,2,3,4),(5,6,7,8),(1,2,6,5),
                                     (2,3,7,6),(3,4,8,7),(4,1,5,8)),
                        prism=((1,2,3),(4,5,6),(1,2,5,4),(2,3,6,5),(1,3,6,4)),
                        pyramid=((1,2,3,4),(1,2,5),(2,3,5),(3,4,5),(4,1,5)))

function _plugin_poly_volume(coords,faces)
    n=size(coords,2)
    c=(sum(@view(coords[1,:])),sum(@view(coords[2,:])),
       sum(@view(coords[3,:])))./n
    volume=0.0
    for face in faces
        a=_plugin_col(coords,face[1])
        for k in 2:length(face)-1
            b=_plugin_col(coords,face[k]);d=_plugin_col(coords,face[k+1])
            volume+=_plugin_tet_volume(c,a,b,d)
        end
    end
    return abs(volume)
end

# Element measure for the scalar integral: the nodal mean over a linear element
# equals the interpolated field integral (exact for simplex and straight-sided
# prism/hex geometries; the standard nodal-average elsewhere).
function _plugin_element_measure(coords,kind::Symbol)
    kind===:point && return 1.0
    if kind===:line
        return _plugin_norm3(_plugin_sub3(_plugin_col(coords,2),
                                          _plugin_col(coords,1)))
    elseif kind===:triangle
        return _plugin_norm3(_plugin_tri_area_vec(_plugin_col(coords,1),
            _plugin_col(coords,2),_plugin_col(coords,3)))
    elseif kind===:quadrangle
        n1=_plugin_tri_area_vec(_plugin_col(coords,1),_plugin_col(coords,2),
                                _plugin_col(coords,3))
        n2=_plugin_tri_area_vec(_plugin_col(coords,1),_plugin_col(coords,3),
                                _plugin_col(coords,4))
        return _plugin_norm3((n1[1]+n2[1],n1[2]+n2[2],n1[3]+n2[3]))
    elseif kind===:tetrahedron
        return abs(_plugin_tet_volume(_plugin_col(coords,1),_plugin_col(coords,2),
                                      _plugin_col(coords,3),_plugin_col(coords,4)))
    end
    faces=get(_PLUGIN_POLY_FACES,kind,nothing)
    faces===nothing && return 0.0
    return _plugin_poly_volume(coords,faces)
end

# Area vector (measure times unit normal) for the flux of a 3-component field
# through a 2-D element.
function _plugin_area_vector(coords,kind::Symbol)
    kind===:triangle && return _plugin_tri_area_vec(_plugin_col(coords,1),
        _plugin_col(coords,2),_plugin_col(coords,3))
    if kind===:quadrangle
        n1=_plugin_tri_area_vec(_plugin_col(coords,1),_plugin_col(coords,2),
                                _plugin_col(coords,3))
        n2=_plugin_tri_area_vec(_plugin_col(coords,1),_plugin_col(coords,3),
                                _plugin_col(coords,4))
        return (n1[1]+n2[1],n1[2]+n2[2],n1[3]+n2[3])
    end
    return (0.0,0.0,0.0)
end

# ---- view outputs ---------------------------------------------------------------

function _plugin_publish_list_view!(name,types,counts,datas,time;caller="API.plugin")
    tag=_fresh_view_tag()
    view=_SessionView(name)
    view.list_types=types;view.list_counts=counts;view.list_data=datas
    view.list_time=copy(time)
    steps=_view_plain_list_steps(view.list_types,view.list_counts,view.list_data)
    view.numbers["NbTimeStep"]=Float64(something(steps,0))
    SESSION_VIEWS[tag]=view
    VIEW_TAG_MAX[]=max(VIEW_TAG_MAX[],tag)
    return tag
end

# Bounding-box accumulation over value-node coordinates (the upstream data
# bounding box anchors the SP records of Integrate/MinMax/Probe outputs).
function _plugin_bbox_bounds!(lo,hi,elements)
    for element in elements
        for node in axes(element.vcoords,2),axis in 1:3
            v=element.vcoords[axis,node]
            v<lo[axis] && (lo[axis]=v)
            v>hi[axis] && (hi[axis]=v)
        end
    end
    return lo,hi
end

function _plugin_bbox_center(elements)
    lo=fill(Inf,3);hi=fill(-Inf,3)
    _plugin_bbox_bounds!(lo,hi,elements)
    any(!isfinite,lo) && return (0.0,0.0,0.0)
    return Tuple((lo.+hi)./2)
end

# `head` letter for a component count and `shape` letter for an element kind —
# the `_POS_RECORD_TYPES` naming used for plugin-emitted list records.
function _plugin_record_type(kind::Symbol,ncomp::Int)
    head=get(_PLUGIN_COMP_HEAD,ncomp,nothing)
    head===nothing && return nothing
    return head*_PLUGIN_KIND_SHAPE[kind]
end

# ---- Integrate ------------------------------------------------------------------

# Per-element contribution to the step integral: scalar fields use
# measure * nodal mean (the linear interpolant integral), 1-node records a
# plain sum (upstream `simpleSum`), 3-component line elements the circulation
# `v̄ · Δx` and 3-component triangle/quadrangle elements the flux `v̄ · n dA`;
# all other component/dimension pairs contribute nothing upstream.
function _plugin_integrate_spatial(elements,step::Int,dimension::Int)
    res=0.0
    for element in elements
        (dimension>0 && element.dim!=dimension) && continue
        arity=size(element.vcoords,2);nc=element.ncomp
        if arity==1
            res+=element.vals[(step-1)*nc+1]
            continue
        end
        if nc==1
            measure=_plugin_element_measure(element.geo,element.kind)
            measure==0.0 && continue
            total=0.0
            for node in 1:arity
                total+=element.vals[(step-1)*arity+node]
            end
            res+=measure*total/arity
        elseif nc==3 && element.dim==1 && element.kind===:line
            measure=_plugin_element_measure(element.geo,element.kind)
            measure==0.0 && continue
            meanv=zeros(3)
            for node in 1:arity,comp in 1:3
                meanv[comp]+=element.vals[(step-1)*arity*3+(node-1)*3+comp]
            end
            edge=_plugin_sub3(_plugin_col(element.geo,2),_plugin_col(element.geo,1))
            res+=_plugin_dot3((meanv[1]/arity,meanv[2]/arity,meanv[3]/arity),edge)
        elseif nc==3 && element.dim==2 && element.kind in (:triangle,:quadrangle)
            meanv=zeros(3)
            for node in 1:arity,comp in 1:3
                meanv[comp]+=element.vals[(step-1)*arity*3+(node-1)*3+comp]
            end
            normal=_plugin_area_vector(element.geo,element.kind)
            res+=_plugin_dot3((meanv[1]/arity,meanv[2]/arity,meanv[3]/arity),normal)
        end
    end
    return res
end

function _plugin_exec_integrate(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    tag,view=_plugin_view(numbers,caller)
    overtime=_mesh_query_integer(Int(numbers["OverTime"]),caller,"OverTime")
    dimension=_mesh_query_integer(Int(numbers["Dimension"]),caller,"Dimension")
    -1<=dimension<=3 || throw(ArgumentError(
        "$caller: Dimension must be -1 or in 0:3"))
    name=view.name*"_Integrate"
    if view.model_kind===nothing
        elements,steps=_plugin_list_elements(view,caller)
        if overtime<0
            cx,cy,cz=_plugin_bbox_center(elements)
            data=[cx,cy,cz]
            for step in 1:steps
                res=_plugin_integrate_spatial(elements,step,dimension)
                _native_diagnostic_emit_locked(:info,
                    "Step $(step-1): integral = $res")
                push!(data,res)
            end
            times=[_plugin_step_time(view,s) for s in 0:steps-1]
            _plugin_publish_list_view!(name,["SP"],Dict("SP"=>1),
                                       Dict("SP"=>data),times;caller)
            return nothing
        end
        # Over-time path (upstream `OverTime = i > -1`): per-node ∫v dt over
        # steps overTime..last with the rectangular rule, scalar payloads only.
        start=Int(overtime)
        start<steps || (start=steps)
        for element in elements
            element.ncomp==1 || throw(ArgumentError(
                "$caller: can only integrate scalar views over time"))
        end
        by_type=Dict{String,Vector{_PluginElement}}();order=String[]
        for element in elements
            type=_plugin_record_type(element.kind,element.ncomp)
            type===nothing && continue
            haskey(by_type,type) || (by_type[type]=_PluginElement[];push!(order,type))
            push!(by_type[type],element)
        end
        types=String[];counts=Dict{String,Int}();datas=Dict{String,Vector{Float64}}()
        for type in order
            members=by_type[type];flat=Float64[]
            time_prev=overtime>0 && start>0 ? _plugin_step_time(view,start-1) : 0.0
            dts=Float64[]
            for s in start:steps-1
                t=_plugin_step_time(view,s)
                push!(dts,t-time_prev);time_prev=t
            end
            for element in members
                arity=size(element.vcoords,2)
                for axis in 1:3,node in 1:arity
                    push!(flat,element.geo[axis,node])
                end
                for node in 1:arity
                    acc=0.0
                    for (i,s) in enumerate(start:steps-1)
                        acc+=element.vals[s*arity+node]*dts[i]
                    end
                    push!(flat,acc)
                end
            end
            push!(types,type);counts[type]=length(members);datas[type]=flat
        end
        _plugin_publish_list_view!(name,types,counts,datas,Float64[];caller)
        return nothing
    end
    # Model-based view path.
    keys=_plugin_model_keys(view)
    isempty(keys) && throw(ArgumentError(
        "$caller: view $tag holds no populated model data steps"))
    overtime>=0 && throw(ArgumentError(
        "$caller: Integrate OverTime is only implemented for list-based views"))
    all_elements=_PluginElement[];res_per=Float64[]
    for key in keys
        elements,_=_plugin_model_elements(view,key,caller)
        append!(all_elements,elements)
        push!(res_per,_plugin_integrate_spatial(elements,1,dimension))
    end
    cx,cy,cz=_plugin_bbox_center(all_elements)
    data=[cx,cy,cz];append!(data,res_per)
    times=[_plugin_step_time(view,key[1]) for key in keys]
    _plugin_publish_list_view!(name,["SP"],Dict("SP"=>1),
                               Dict("SP"=>data),times;caller)
    return nothing
end

# ---- MinMax ---------------------------------------------------------------------

# Scalar-representation scan over one step: (min, max, argmin, argmax) where the
# arg sites use the value-node coordinates like the upstream `getNode` calls.
function _plugin_minmax_scan(elements,step::Int)
    vmin=Inf;vmax=-Inf
    amin=(0.0,0.0,0.0);amax=(0.0,0.0,0.0)
    for element in elements
        arity=size(element.vcoords,2);nc=element.ncomp
        for node in 1:arity
            val=_plugin_scalar_rep(element.vals,(step-1)*arity*nc+(node-1)*nc,nc)
            if val<vmin
                vmin=val
                amin=_plugin_col(element.vcoords,node)
            end
            if val>vmax
                vmax=val
                amax=_plugin_col(element.vcoords,node)
            end
        end
    end
    return vmin,vmax,amin,amax
end

function _plugin_exec_minmax(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    tag,view=_plugin_view(numbers,caller)
    overtime=!iszero(numbers["OverTime"])
    argument=!iszero(numbers["Argument"])
    name=view.name
    # Per-step element source: list views share one element set across steps,
    # model views re-resolve each populated step.
    list_elements=_PluginElement[];keys=Tuple{Int,String}[]
    if view.model_kind===nothing
        list_elements,steps=_plugin_list_elements(view,caller)
        step_nums=collect(0:steps-1)
    else
        keys=_plugin_model_keys(view)
        isempty(keys) && throw(ArgumentError(
            "$caller: view $tag holds no populated model data steps"))
        steps=length(keys)
        step_nums=[key[1] for key in keys]
    end
    function step_elements(step)
        view.model_kind===nothing && return list_elements
        return _plugin_model_elements(view,keys[step],caller)[1]
    end
    times=[_plugin_step_time(view,s) for s in step_nums]
    if !overtime
        if argument
            mins=Float64[];maxs=Float64[]
            for step in 1:steps
                vmin,vmax,amin,amax=_plugin_minmax_scan(
                    step_elements(step),view.model_kind===nothing ? step : 1)
                append!(mins,(amin[1],amin[2],amin[3],vmin))
                append!(maxs,(amax[1],amax[2],amax[3],vmax))
            end
            _plugin_publish_list_view!(name*"_Min",["SP"],Dict("SP"=>steps),
                                       Dict("SP"=>mins),Float64[];caller)
            _plugin_publish_list_view!(name*"_Max",["SP"],Dict("SP"=>steps),
                                       Dict("SP"=>maxs),Float64[];caller)
        else
            lo=fill(Inf,3);hi=fill(-Inf,3)
            for step in 1:steps
                _plugin_bbox_bounds!(lo,hi,step_elements(step))
            end
            center=any(!isfinite,lo) ? (0.0,0.0,0.0) : Tuple((lo.+hi)./2)
            mins=[center...];maxs=[center...]
            for step in 1:steps
                vmin,vmax,_,_=_plugin_minmax_scan(
                    step_elements(step),view.model_kind===nothing ? step : 1)
                push!(mins,vmin);push!(maxs,vmax)
            end
            _plugin_publish_list_view!(name*"_Min",["SP"],Dict("SP"=>1),
                                       Dict("SP"=>mins),times;caller)
            _plugin_publish_list_view!(name*"_Max",["SP"],Dict("SP"=>1),
                                       Dict("SP"=>maxs),times;caller)
        end
        return nothing
    end
    gmin=Inf;gmax=-Inf;tmin=0.0;tmax=0.0;center=(0.0,0.0,0.0)
    for step in 1:steps
        elements=step_elements(step)
        vmin,vmax,_,_=_plugin_minmax_scan(
            elements,view.model_kind===nothing ? step : 1)
        step==1 && (center=_plugin_bbox_center(elements))
        if vmin<gmin;gmin=vmin;tmin=times[step];end
        if vmax>gmax;gmax=vmax;tmax=times[step];end
    end
    _plugin_publish_list_view!(name*"_Min",["SP"],Dict("SP"=>1),
                               Dict("SP"=>[center...,gmin]),[tmin];caller)
    _plugin_publish_list_view!(name*"_Max",["SP"],Dict("SP"=>1),
                               Dict("SP"=>[center...,gmax]),[tmax];caller)
    return nothing
end

# ---- Probe ------------------------------------------------------------------------

# Barycentric node weights at `p` for a linear element, or `nothing` when `p`
# lies outside (relative tolerance) — the list-data counterpart of the
# session-mesh `_local_coordinates` inversion.
function _plugin_local_weights(coords,kind::Symbol,p)
    T=1e-9
    a=_plugin_col(coords,1)
    if kind===:line
        e=_plugin_sub3(_plugin_col(coords,2),a)
        scale=_plugin_norm3(e)
        scale>0 || return nothing
        d=_plugin_sub3(p,a)
        u=_plugin_dot3(d,e)/_plugin_dot3(e,e)
        resid=_plugin_sub3(d,(e[1]*u,e[2]*u,e[3]*u))
        _plugin_norm3(resid)<=T*scale && -T<=u<=1+T || return nothing
        return (1.0-u,u)
    elseif kind===:triangle
        e1=_plugin_sub3(_plugin_col(coords,2),a)
        e2=_plugin_sub3(_plugin_col(coords,3),a)
        g11=_plugin_dot3(e1,e1);g12=_plugin_dot3(e1,e2);g22=_plugin_dot3(e2,e2)
        det=g11*g22-g12*g12
        scale=_plugin_norm3(_plugin_cross3(e1,e2))
        det>0 && scale>0 || return nothing
        d=_plugin_sub3(p,a)
        u=(g22*_plugin_dot3(d,e1)-g12*_plugin_dot3(d,e2))/det
        v=(g11*_plugin_dot3(d,e2)-g12*_plugin_dot3(d,e1))/det
        resid=_plugin_sub3(d,(e1[1]*u+e2[1]*v,e1[2]*u+e2[2]*v,e1[3]*u+e2[3]*v))
        _plugin_norm3(resid)<=T*scale &&
            u>=-T && v>=-T && u+v<=1+T || return nothing
        return (1.0-u-v,u,v)
    elseif kind===:quadrangle
        for nodes in ((1,2,3),(1,3,4))
            sub=coords[:,collect(nodes)]
            w=_plugin_local_weights(sub,:triangle,p)
            w===nothing && continue
            N=zeros(4)
            for (i,node) in enumerate(nodes);N[node]=w[i];end
            return Tuple(N)
        end
        return nothing
    elseif kind===:tetrahedron
        e1=_plugin_sub3(_plugin_col(coords,2),a)
        e2=_plugin_sub3(_plugin_col(coords,3),a)
        e3=_plugin_sub3(_plugin_col(coords,4),a)
        det=e1[1]*(e2[2]*e3[3]-e2[3]*e3[2])-
            e2[1]*(e1[2]*e3[3]-e1[3]*e3[2])+
            e3[1]*(e1[2]*e2[3]-e1[3]*e2[2])
        scale=abs(det)
        scale>0 || return nothing
        d=_plugin_sub3(p,a)
        u=(d[1]*(e2[2]*e3[3]-e2[3]*e3[2])-
           e2[1]*(d[2]*e3[3]-d[3]*e3[2])+
           e3[1]*(d[2]*e2[3]-d[3]*e2[2]))/det
        v=(e1[1]*(d[2]*e3[3]-d[3]*e3[2])-
           d[1]*(e1[2]*e3[3]-e1[3]*e3[2])+
           e3[1]*(e1[2]*d[3]-e1[3]*d[2]))/det
        w=(e1[1]*(e2[2]*d[3]-e2[3]*d[2])-
           e2[1]*(e1[2]*d[3]-e1[3]*d[2])+
           d[1]*(e1[2]*e2[3]-e1[3]*e2[2]))/det
        l0=1.0-u-v-w
        minimum((l0,u,v,w))>=-T && maximum((l0,u,v,w))<=1+T || return nothing
        return (l0,u,v,w)
    end
    return nothing
end

# Locate the first element (in deterministic record order) containing `p` and
# interpolate every populated step: returns (per-step value vectors, steps).
function _plugin_probe_list(view,p,caller::AbstractString)
    elements,steps=_plugin_list_elements(view,caller)
    for element in elements
        element.dim<1 && continue
        N=_plugin_local_weights(element.geo,element.kind,p)
        N===nothing && continue
        arity=size(element.vcoords,2);nc=element.ncomp
        out=Vector{Float64}[]
        for step in 1:steps
            acc=zeros(nc)
            for node in 1:min(arity,length(N)),comp in 1:nc
                acc[comp]+=N[node]*
                    element.vals[(step-1)*arity*nc+(node-1)*nc+comp]
            end
            push!(out,acc)
        end
        return out,collect(0:steps-1)
    end
    return nothing,Int[]
end

# Session-mesh probe over every populated step, mirroring `_view_probe`'s
# locator + linear-hull barycentric interpolation (per-step slot resolution).
function _plugin_probe_model(view,p,caller::AbstractString)
    keys=_plugin_model_keys(view)
    isempty(keys) && return nothing,Int[]
    first_ctx=view.model_contexts[first(keys)]
    slot=first_ctx.model_slot
    mesh=_plugin_model_mesh(slot,caller)
    locator=if mesh isa Mesh
        slot.locator===nothing && (slot.locator=SimplexLocator(mesh))
        slot.locator
    elseif mesh isa MixedMesh
        _MixedMeshLocator(mesh)
    else
        throw(ArgumentError("$caller: probe requires a mesh"))
    end
    tags=_locate_elements(locator,p,-1,true,caller)
    isempty(tags) && return nothing,Int[]
    etag=Int(first(tags))
    rec=mesh_element_record(mesh,etag)
    uvw,=_local_coordinates(mesh,etag,p,caller)
    u,v,w=uvw
    N=rec.dimension==1 ? (1.0-u,u) :
      rec.dimension==2 ? (1.0-u-v,u,v) :
                         (1.0-u-v-w,u,v,w)
    public=_mesh_public_tags(mesh)
    kind=view.model_kind
    out=Vector{Float64}[];step_nums=Int[]
    for key in keys
        ctx=view.model_contexts[key]
        get(ctx,:model_slot,nothing)===slot || continue
        entity_tags,values,nc=view.model_data[key]
        counts=get(ctx,:counts,Int32[])
        offsets=Vector{Int}(undef,length(entity_tags)+1);offsets[1]=1
        for row in eachindex(entity_tags)
            width=row<=length(counts) ? Int(counts[row]) :
                length(values)÷length(entity_tags)
            offsets[row+1]=Base.checked_add(offsets[row],width)
        end
        acc=nothing
        if kind=="ElementData" || kind=="ElementNodeData"
            for (row,tag) in enumerate(entity_tags)
                _plugin_element_index(mesh,public,tag)==etag || continue
                base=offsets[row]
                if kind=="ElementData"
                    acc=collect(Float64,@view(values[base:base+nc-1]))
                else
                    nn=(offsets[row+1]-base)÷nc
                    acc=zeros(nc)
                    for node in 1:min(nn,length(N)),comp in 1:nc
                        acc[comp]+=N[node]*values[base+(node-1)*nc+comp-1]
                    end
                end
                break
            end
        elseif kind=="NodeData"
            row_of=Dict{UInt64,Int}(tag=>row for (row,tag) in enumerate(entity_tags))
            trial=zeros(nc)
            for node in 1:min(length(rec.node_tags),length(N))
                tag=_cache_node_tag(mesh,Int(rec.node_tags[node]))
                row=get(row_of,tag,0)
                row==0 && continue
                acc===nothing && (acc=trial)
                for comp in 1:nc
                    trial[comp]+=N[node]*values[offsets[row]+comp-1]
                end
            end
        end
        acc===nothing && continue
        push!(out,acc);push!(step_nums,key[1])
    end
    return out,step_nums
end

function _plugin_probe_emit(view,name,p,values,step_nums,caller)
    isempty(values) &&
        return _plugin_publish_list_view!(name,String[],Dict{String,Int}(),
            Dict{String,Vector{Float64}}(),Float64[];caller)
    nc=length(first(values))
    type=get(_PLUGIN_COMP_HEAD,nc,nothing)
    type===nothing && throw(ArgumentError(
        "$caller: probe result with $nc components cannot be emitted"))
    record_type=type*"P"
    data=[p[1],p[2],p[3]]
    for acc in values,comp in acc
        push!(data,comp)
    end
    times=[_plugin_step_time(view,s) for s in step_nums]
    _plugin_publish_list_view!(name,[record_type],Dict(record_type=>1),
                               Dict(record_type=>data),times;caller)
end

function _plugin_exec_probe(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    p=(_mesh_query_coordinate(numbers["X"],caller,"X"),
       _mesh_query_coordinate(numbers["Y"],caller,"Y"),
       _mesh_query_coordinate(numbers["Z"],caller,"Z"))
    tag,view=_plugin_view(numbers,caller)
    name=view.name*"_Probe"
    if view.model_kind===nothing
        values,step_nums=_plugin_probe_list(view,p,caller)
    else
        values,step_nums=_plugin_probe_model(view,p,caller)
    end
    _plugin_probe_emit(view,name,p,something(values,Vector{Float64}[]),
                       step_nums,caller)
    return nothing
end

# ---- Skin -------------------------------------------------------------------------

# Upstream `getBoundary` facet tables (0-based in `Skin.cpp`, 1-based here):
# 2-D elements contribute their edges, 3-D elements their faces; other kinds
# carry no boundary facets.
const _PLUGIN_SKIN_FACETS=(
    triangle=((1,2),(2,3),(3,1)),
    quadrangle=((1,2),(2,3),(3,4),(4,1)),
    tetrahedron=((1,2,4),(1,3,2),(1,4,3),(2,3,4)),
    hexahedron=((1,2,6,5),(1,4,3,2),(1,5,8,4),(2,3,7,6),(3,4,8,7),(5,6,7,8)),
    prism=((1,2,5,4),(1,4,6,3),(2,3,6,5),(1,3,2),(4,5,6)),
    pyramid=((1,4,3,2),(1,2,5),(1,5,4),(2,3,5),(3,4,5)))

function _plugin_exec_skin(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    iszero(numbers["FromMesh"]) || throw(ArgumentError(
        "$caller: Skin FromMesh=1 (discrete model boundary extraction) " *
        "is not supported; run on a list-based view"))
    tag,view=_plugin_view(numbers,caller)
    view.model_kind===nothing || throw(ArgumentError(
        "$caller: Skin requires a list-based view"))
    elements,steps=_plugin_list_elements(view,caller)
    # Symmetric-difference bookkeeping like the upstream `std::set<ElmData>`
    # insert-or-erase scan: a facet seen an even number of times cancels.
    held=Dict{Vector{NTuple{3,Float64}},
              Tuple{_PluginElement,Tuple{Vararg{Int}}}}()
    for element in elements
        for facet in get(_PLUGIN_SKIN_FACETS,element.kind,())
            key=sort!(NTuple{3,Float64}[_plugin_col(element.geo,node)
                                        for node in facet])
            if haskey(held,key)
                delete!(held,key)
            else
                held[key]=(element,facet)
            end
        end
    end
    grouped=Dict{String,Vector{Tuple{_PluginElement,Tuple{Vararg{Int}}}}}()
    order=String[]
    for key in sort!(collect(keys(held)))
        element,facet=held[key]
        kind=length(facet)==2 ? :line : length(facet)==3 ? :triangle : :quadrangle
        record_type=_plugin_record_type(kind,element.ncomp)
        record_type===nothing && continue
        if !haskey(grouped,record_type)
            grouped[record_type]=Tuple{_PluginElement,Tuple{Vararg{Int}}}[]
            push!(order,record_type)
        end
        push!(grouped[record_type],(element,facet))
    end
    types=String[];counts=Dict{String,Int}();datas=Dict{String,Vector{Float64}}()
    for type in order
        members=grouped[type];flat=Float64[]
        for (element,facet) in members
            nc=element.ncomp;va=size(element.vcoords,2)
            for axis in 1:3,node in facet
                push!(flat,element.geo[axis,node])
            end
            for step in 1:steps,node in facet,comp in 1:nc
                push!(flat,element.vals[(step-1)*va*nc+(node-1)*nc+comp])
            end
        end
        push!(types,type);counts[type]=length(members);datas[type]=flat
    end
    times=[_plugin_step_time(view,s) for s in 0:steps-1]
    _plugin_publish_list_view!(view.name*"_Skin",types,counts,datas,times;caller)
    return nothing
end

# ---- ExtractElements ------------------------------------------------------------------

# Emit one list record per kept element, restricted to `emit_step` when it is
# non-negative (0-based source index) or covering every populated step.
function _plugin_extract_emit(name,elements,emit_step,nsteps,times,caller)
    by_type=Dict{String,Vector{_PluginElement}}();order=String[]
    for element in elements
        type=_plugin_record_type(element.kind,element.ncomp)
        type===nothing && continue
        haskey(by_type,type) || (by_type[type]=_PluginElement[];push!(order,type))
        push!(by_type[type],element)
    end
    types=String[];counts=Dict{String,Int}();datas=Dict{String,Vector{Float64}}()
    sel=emit_step<0 ? collect(1:nsteps) : [emit_step+1]
    for type in order
        members=by_type[type];flat=Float64[]
        for element in members
            arity=size(element.vcoords,2);nc=element.ncomp
            for axis in 1:3,node in 1:arity
                push!(flat,element.geo[axis,node])
            end
            for step in sel,node in 1:arity,comp in 1:nc
                push!(flat,element.vals[(step-1)*arity*nc+(node-1)*nc+comp])
            end
        end
        push!(types,type);counts[type]=length(members);datas[type]=flat
    end
    out_times=emit_step<0 ? times : Float64[times[emit_step+1]]
    _plugin_publish_list_view!(name,types,counts,datas,out_times;caller)
end

function _plugin_exec_extract_elements(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    tag,view=_plugin_view(numbers,caller)
    minval=numbers["MinVal"];maxval=numbers["MaxVal"]
    check=minval!=maxval
    this_step=_mesh_query_integer(Int(numbers["TimeStep"]),caller,"TimeStep")
    dimension=_mesh_query_integer(Int(numbers["Dimension"]),caller,"Dimension")
    -1<=dimension<=3 || throw(ArgumentError(
        "$caller: Dimension must be -1 or in 0:3"))
    name=view.name*"_ExtractElements"
    if view.model_kind===nothing
        elements,steps=_plugin_list_elements(view,caller)
        eval_step=this_step<0 ? 0 : this_step
        if eval_step>steps-1
            _native_diagnostic_emit_locked(:error,
                "Invalid time step ($eval_step) in View[$tag]: using first step instead")
            eval_step=0
        end
        kept=_PluginElement[]
        for element in elements
            (dimension>0 && element.dim!=dimension) && continue
            if check
                arity=size(element.vcoords,2);nc=element.ncomp
                d=0.0
                for node in 1:arity
                    d+=_plugin_scalar_rep(element.vals,
                        eval_step*arity*nc+(node-1)*nc,nc)
                end
                d/=arity
                (d<minval || d>=maxval) && continue
            end
            push!(kept,element)
        end
        times=[_plugin_step_time(view,s) for s in 0:steps-1]
        _plugin_extract_emit(name,kept,this_step<0 ? -1 : eval_step,steps,
                             times,caller)
        return nothing
    end
    view.model_kind=="NodeData" && throw(ArgumentError(
        "$caller: ExtractElements requires element-valued model data"))
    keys=_plugin_model_keys(view)
    isempty(keys) && throw(ArgumentError(
        "$caller: view $tag holds no populated model data steps"))
    eval_index=if this_step<0
        1
    else
        found=findfirst(key->key[1]==this_step,keys)
        if found===nothing
            _native_diagnostic_emit_locked(:error,
                "Invalid time step ($this_step) in View[$tag]: using first step instead")
            1
        else
            found
        end
    end
    eval_elements,eval_nc=_plugin_model_elements(view,keys[eval_index],caller)
    kept_tags=Set{UInt64}()
    for element in eval_elements
        (dimension>0 && element.dim!=dimension) && continue
        if check
            arity=size(element.vcoords,2)
            d=0.0
            for node in 1:arity
                d+=_plugin_scalar_rep(element.vals,(node-1)*eval_nc,eval_nc)
            end
            d/=arity
            (d<minval || d>=maxval) && continue
        end
        push!(kept_tags,element.etag)
    end
    emit_keys=this_step<0 ? keys : [keys[eval_index]]
    # Merge the kept elements across the emitted steps into single records
    # (unpopulated steps contribute zeros for a tag that lacks a row there).
    ordered=UInt64[element.etag for element in eval_elements
                   if element.etag in kept_tags]
    per_step=Vector{Dict{UInt64,_PluginElement}}(undef,length(emit_keys))
    for (i,key) in enumerate(emit_keys)
        step_elements,_=_plugin_model_elements(view,key,caller)
        per_step[i]=Dict{UInt64,_PluginElement}(e.etag=>e for e in step_elements)
    end
    merged=_PluginElement[]
    for etag in ordered
        first_element=nothing
        for table in per_step
            first_element=get(table,etag,nothing)
            first_element===nothing || break
        end
        first_element===nothing && continue
        nc=first_element.ncomp;arity=size(first_element.vcoords,2)
        vals=Vector{Float64}(undef,arity*nc*length(emit_keys))
        for (i,table) in enumerate(per_step)
            element=get(table,etag,nothing)
            for node in 1:arity,comp in 1:nc
                vals[(i-1)*arity*nc+(node-1)*nc+comp]=
                    element===nothing ? 0.0 :
                    element.vals[(node-1)*nc+comp]
            end
        end
        push!(merged,_PluginElement(first_element.kind,first_element.dim,
            first_element.geo,first_element.vcoords,vals,nc,etag))
    end
    _plugin_extract_emit(name,merged,
        this_step<0 ? -1 : 0,length(emit_keys),
        [_plugin_step_time(view,key[1]) for key in emit_keys],caller)
    return nothing
end

# ---- Transform ------------------------------------------------------------------------

function _plugin_exec_transform(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    tag,view=_plugin_view(numbers,caller)
    view.model_kind===nothing || throw(ArgumentError(
        "$caller: Transform rewrites list data only; model views are unsupported"))
    a=((numbers["A11"],numbers["A12"],numbers["A13"]),
       (numbers["A21"],numbers["A22"],numbers["A23"]),
       (numbers["A31"],numbers["A32"],numbers["A33"]))
    t=(numbers["Tx"],numbers["Ty"],numbers["Tz"])
    swap=!iszero(numbers["SwapOrientation"])
    steps=_plugin_list_steps(view,caller)
    for type in view.list_types
        layout=get(Post._POS_RECORD_TYPES,type,nothing)
        layout===nothing && throw(ArgumentError(
            "$caller: unsupported list data type $type"))
        _,arity,ncomp=layout
        count=view.list_counts[type];data=view.list_data[type]
        width=3arity+arity*ncomp*steps
        for element in 0:count-1
            base=element*width
            coords=Matrix{Float64}(undef,3,arity)
            for node in 1:arity
                x=data[base+node];y=data[base+arity+node]
                z=data[base+2arity+node]
                coords[1,node]=a[1][1]*x+a[1][2]*y+a[1][3]*z+t[1]
                coords[2,node]=a[2][1]*x+a[2][2]*y+a[2][3]*z+t[2]
                coords[3,node]=a[3][1]*x+a[3][2]*y+a[3][3]*z+t[3]
            end
            order=swap ? collect(arity:-1:1) : collect(1:arity)
            for axis in 1:3,node in 1:arity
                data[base+(axis-1)*arity+node]=coords[axis,order[node]]
            end
            swap || continue
            # Upstream `reverseElement` also reverses the node index of every
            # value block inside each step.
            for step in 1:steps
                step_base=base+3arity+(step-1)*arity*ncomp
                vals=collect(Float64,
                    @view(data[step_base+1:step_base+arity*ncomp]))
                for node in 1:arity,comp in 1:ncomp
                    data[step_base+(node-1)*ncomp+comp]=
                        vals[(order[node]-1)*ncomp+comp]
                end
            end
        end
    end
    return nothing
end

# ---- Remove ------------------------------------------------------------------------------

const _PLUGIN_REMOVE_KINDS=(point="Points",line="Lines",triangle="Triangles",
    quadrangle="Quadrangles",tetrahedron="Tetrahedra",hexahedron="Hexahedra",
    prism="Prisms",pyramid="Pyramids")
const _PLUGIN_REMOVE_COMPS=Dict{Int,String}(1=>"Scalar",3=>"Vector",9=>"Tensor")

function _plugin_exec_remove(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    tag,view=_plugin_view(numbers,caller)
    view.model_kind===nothing || throw(ArgumentError(
        "$caller: Remove rewrites list data only; model views are unsupported"))
    if !iszero(numbers["Text3D"])
        empty!(view.list_string_types);empty!(view.list_strings)
        empty!(view.list_string_coords);empty!(view.list_string_counts)
    end
    # No separate 2-D string class exists in session views; Text2D is accepted
    # and stored upstream but maps to nothing here.
    comps=Set{String}(name for name in values(_PLUGIN_REMOVE_COMPS)
                      if !iszero(numbers[name]))
    isempty(comps) && return nothing
    kept=String[]
    for type in view.list_types
        layout=get(Post._POS_RECORD_TYPES,type,nothing)
        layout===nothing && (push!(kept,type);continue)
        kind,_,ncomp=layout
        shape=get(_PLUGIN_REMOVE_KINDS,kind,nothing)
        if shape!==nothing && !iszero(numbers[shape]) &&
                get(_PLUGIN_REMOVE_COMPS,ncomp,"") in comps
            delete!(view.list_counts,type)
            delete!(view.list_data,type)
            continue
        end
        push!(kept,type)
    end
    view.list_types=kept
    return nothing
end

# ---- Smooth --------------------------------------------------------------------------------

# Average the values of all coincident value-nodes (exact coordinate equality:
# list records duplicate shared vertices verbatim), matching upstream
# `PViewDataList::smooth` for all component widths and steps.
function _plugin_exec_smooth(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    tag,view=_plugin_view(numbers,caller)
    view.model_kind===nothing || throw(ArgumentError(
        "$caller: Smooth rewrites list data only; model views are unsupported"))
    steps=_plugin_list_steps(view,caller)
    sums=Dict{NTuple{3,Float64},Vector{Float64}}()
    members=Dict{NTuple{3,Float64},Int}()
    for type in view.list_types
        layout=get(Post._POS_RECORD_TYPES,type,nothing)
        layout===nothing && continue
        _,arity,ncomp=layout
        arity>1 || continue
        data=view.list_data[type];count=view.list_counts[type]
        width=3arity+arity*ncomp*steps
        for element in 0:count-1
            base=element*width
            for node in 1:arity
                key=(data[base+node],data[base+arity+node],data[base+2arity+node])
                if !haskey(sums,key)
                    sums[key]=zeros(steps*ncomp);members[key]=0
                end
                members[key]+=1
                acc=sums[key]
                for step in 1:steps,comp in 1:ncomp
                    acc[(step-1)*ncomp+comp]+=
                        data[base+3arity+(step-1)*arity*ncomp+(node-1)*ncomp+comp]
                end
            end
        end
    end
    for (key,acc) in sums
        acc./=members[key]
    end
    for type in view.list_types
        layout=get(Post._POS_RECORD_TYPES,type,nothing)
        layout===nothing && continue
        _,arity,ncomp=layout
        arity>1 || continue
        data=view.list_data[type];count=view.list_counts[type]
        width=3arity+arity*ncomp*steps
        for element in 0:count-1
            base=element*width
            for node in 1:arity
                key=(data[base+node],data[base+arity+node],data[base+2arity+node])
                acc=get(sums,key,nothing)
                acc===nothing && continue
                for step in 1:steps,comp in 1:ncomp
                    data[base+3arity+(step-1)*arity*ncomp+(node-1)*ncomp+comp]=
                        acc[(step-1)*ncomp+comp]
                end
            end
        end
    end
    return nothing
end

# ---- NewView ---------------------------------------------------------------------------------

function _plugin_mesh_msh_types(mesh)
    mesh isa MixedMesh &&
        return Tuple(Int(block.msh) for block in mesh.blocks if !isempty(block.nodes))
    mesh isa Mesh &&
        return Tuple(m for (m,cells) in ((1,mesh.segs),(2,mesh.tris),(4,mesh.tets))
                     if !isempty(cells))
    return ()
end

function _plugin_exec_new_view(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers;strings=spec.strings
    ncomp=_mesh_query_integer(Int(numbers["NumComp"]),caller,"NumComp")
    ncomp>=1 || throw(ArgumentError(
        "$caller: bad number of components for Plugin(NewView)"))
    value=numbers["Value"]
    type=strings["Type"]
    type=="NodeData" || throw(ArgumentError(
        "$caller: unknown data type '$type' for Plugin(NewView)"))
    _slot_capture_locked!()
    index=_current_slot_index_locked()
    index==0 && throw(ArgumentError("$caller: no current model"))
    slot=MODEL_SLOTS[index]
    mesh=slot.mesh
    mesh===nothing && throw(ArgumentError(
        "$caller: no mesh available to create the view: please mesh your model"))
    public=_mesh_public_tags(mesh)
    tags=public===nothing ?
        UInt64[UInt64(i) for i in 1:nnodes(mesh)] : copy(public.node_tags)
    flat=fill(Float64(value),length(tags)*ncomp)
    counts=fill(Int32(ncomp),length(tags))
    extent=_native_model_view_extent(tags,_native_view_model_count(slot,"NodeData"))
    raw_tag=numbers["ViewTag"]
    raw_tag==floor(raw_tag) || throw(ArgumentError("$caller: ViewTag must be an integer"))
    tag=raw_tag>0 ? Int(raw_tag) : _fresh_view_tag()
    tag<=typemax(Int32) || throw(ArgumentError("$caller: ViewTag is out of range"))
    if haskey(SESSION_VIEWS,tag)
        delete!(SESSION_VIEWS,tag)
        _native_diagnostic_emit_locked(:info,"Removing existing View (tag = $tag)")
    end
    view=_SessionView("New view")
    view.model_kind="NodeData"
    view.model_step_ranges=[_NativeViewStepRange(0,0,slot,ncomp)]
    view.model_data[(0,"NodeData")]=(tags,flat,ncomp)
    view.model_contexts[(0,"NodeData")]=(
        model_slot=slot,time=0.0,counts=counts,extent=extent,
        partitions=Set{Int}([0]),
        interpolation_msh=_plugin_mesh_msh_types(mesh))
    SESSION_VIEWS[tag]=view
    VIEW_TAG_MAX[]=max(VIEW_TAG_MAX[],tag)
    return nothing
end

# ---- registry --------------------------------------------------------------------------------

const _PLUGINS=Dict{String,NamedTuple}(
    "Integrate"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"View"=>-1.0,"OverTime"=>-1.0,"Dimension"=>-1.0,"Visible"=>1.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_integrate),
    "MinMax"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"View"=>-1.0,"OverTime"=>0.0,"Argument"=>0.0,"Visible"=>1.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_minmax),
    "Probe"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"X"=>0.0,"Y"=>0.0,"Z"=>0.0,"View"=>-1.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_probe),
    "Skin"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"Visible"=>1.0,"FromMesh"=>0.0,"View"=>-1.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_skin),
    "ExtractElements"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"MinVal"=>0.0,"MaxVal"=>0.0,"TimeStep"=>0.0,"Visible"=>1.0,
        "Dimension"=>-1.0,"View"=>-1.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_extract_elements),
    "Transform"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"A11"=>1.0,"A12"=>0.0,"A13"=>0.0,
        "A21"=>0.0,"A22"=>1.0,"A23"=>0.0,
        "A31"=>0.0,"A32"=>0.0,"A33"=>1.0,
        "Tx"=>0.0,"Ty"=>0.0,"Tz"=>0.0,"SwapOrientation"=>0.0,"View"=>-1.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_transform),
    "Remove"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"Text2D"=>1.0,"Text3D"=>1.0,"Points"=>0.0,"Lines"=>0.0,
        "Triangles"=>0.0,"Quadrangles"=>0.0,"Tetrahedra"=>0.0,"Hexahedra"=>0.0,
        "Prisms"=>0.0,"Pyramids"=>0.0,"Scalar"=>1.0,"Vector"=>1.0,"Tensor"=>1.0,
        "View"=>-1.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_remove),
    "Smooth"=>(number_options=Dict{String,Float64}("Run"=>0.0,"View"=>-1.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_smooth),
    "NewView"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"NumComp"=>1.0,"Value"=>0.0,"ViewTag"=>-1.0),
        string_options=Dict{String,String}("Type"=>"NodeData"),
        run_fn=_plugin_exec_new_view),
)
