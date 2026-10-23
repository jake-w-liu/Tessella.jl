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
#   NewView         constant `NodeData` model view over the current mesh,
#   Gradient        per-node spatial gradient of scalar data into `*_Gradient`
#                   vector records,
#   CutPlane        `A x + B y + C z + D = 0` levelset sections into
#                   `*_Levelset` point/line/triangle/quadrangle records
#                   (ExtractVolume keeps the whole elements on the selected
#                   side),
#   CutSphere       `(X-Xc)^2 + (Y-Yc)^2 + (Z-Zc)^2 = R^2` levelset sections
#                   into `*_Levelset` records through the same linear-edge
#                   levelset interpolation (sections oriented by the
#                   levelset-gradient rule like upstream),
#   Isosurface      `scalar value == Value` level sets of the view's own
#                   data into `*_Levelset_<step>` point/line/triangle/
#                   quadrangle
#                   records — one list view per populated step since the
#                   section geometry is value-dependent; non-simplex
#                   elements cut through upstream's simplex decomposition
#                   (2 tris / 6 tets / 3 tets / 2 tets) and `ExtractVolume`
#                   keeps the whole elements on the selected side,
#   Warp            in-place coordinate displacement: `Factor` times the
#                   `TimeStep` vector field of `OtherView`, or (when both
#                   `View` and `OtherView` are < 0) smoothed element normals
#                   times the node values, like upstream `Plugin(Warp)`,
#   MathEval        per-node expression evaluation (`Expression0`..`Expression8`
#                   over x, y, z, `v0`..`v8` view components and `w0`..`w8`
#                   `OtherView` components) into a `*_MathEval` list view
#                   (scalar/vector/tensor records by expression count),
#   Scal2Vec        packs component 0 of `ViewX`/`ViewY`/`ViewZ` at the
#                   reference view's (element, node, step) slots into a
#                   `NameNewView` vector list view,
#   Divergence      per-node divergence ∇·v of 3-component data into
#                   `*_Divergence` scalar records,
#   Curl            per-node curl ∇×v of 3-component data into `*_Curl`
#                   vector records,
#   Lambda2         per-node `Eigenvalue`-th (descending order) eigenvalue of
#                   S·S + Ω·Ω built from the 3-component velocity gradient —
#                   or read straight from 9-component tensor records — into
#                   `*_Lambda2` scalar records,
#   SphericalRaise  in-place radial raise: each list-record node moves along
#                   the unit vector from (`Xc`,`Yc`,`Zc`) by `Offset + Raise`
#                   times its `TimeStep` scalar representation,
#   Annotate        appends the `Text` string as a `T2`/`T3` string record
#                   (in-place on list views; a new `*_Annotate` list view
#                   carrying the string and the source time table for model
#                   views).
#   Eigenvalues     per-node real eigenvalue triple of each 9-component
#                   tensor record into `*_MinEigenvalues`,
#                   `*_MidEigenvalues` and `*_MaxEigenvalues` scalar views,
#   ModulusPhase    in-place rewrite of the `RealPart`/`ImaginaryPart` time
#                   steps into per-node modulus and phase values,
#   HarmonicToTime  `RealPart`/`ImaginaryPart` steps resampled to `NumSteps`
#                   harmonic samples in a `*_HarmonicToTime` list view,
#   Density         per-node element-incidence count into a single-step
#                   `*_Density` scalar view (the upstream plugin predates
#                   Gmsh 4.15.2; `View` is its only surviving option).
#
# Divergence, Curl and Lambda2 differentiate through the corner-node linear
# map like upstream `getNumCornerNodes`/`element::interpolate*`, so
# high-order layouts emit linear records on their corners; records without
# the required component count (3, or 3/9 for Lambda2) are skipped, and
# point records emit the upstream zero gradients (identity Jacobian with
# vanishing shape gradients). Divergence and Curl reject multi-mesh model
# views like upstream `hasMultipleMeshes`; Lambda2 shares the check because
# the per-tag step merge cannot disambiguate same-numbered elements across
# meshes, and `Eigenvalues`, `ModulusPhase`, `HarmonicToTime` and `Density`
# run the same guard (`Eigenvalues`/`ModulusPhase`/`HarmonicToTime` have the
# identical upstream check). `Eigenvalues`, `HarmonicToTime` and `Density`
# emit linear corner records like the other per-element plugins.
# `SphericalRaise` rewrites list data only (upstream mutates the
# owning mesh through model views, the same restriction the other in-place
# plugins take). `Annotate` accepts the upstream `Font`/`FontSize`/`Align`
# options but session string records carry no per-string style slot, so
# they are stored without effect (the `Visible` convention below).
#
# Gradient and CutPlane read list views directly and model views through a
# per-element multi-step merge (`_plugin_elements_all_steps`, shared by
# Eigenvalues, HarmonicToTime and Density), so both emit list views covering
# every populated step. `ModulusPhase` is the only in-place plugin that also
# rewrites model payloads — its stored rows are per tag, so the upstream
# NodeData double-application guard is unnecessary. Gradient differentiates
# only linear nodal data — records whose value nodes outnumber the primary
# geometric map (high-order `ElementNodeData`), non-scalar components and
# point elements are skipped like the upstream scalar-record scan. CutPlane
# interpolates section values from the primary corner values; elements fully
# in the plane are re-emitted verbatim when they fit a list record type.
# Upstream's `RecurLevel`/`TargetError` options exist for the generic
# levelset subdivision of curved elements; linear sections are exact without
# them, so they are not registered.
#
# Session views carry no per-element visibility flags, so the upstream `Visible`
# options are accepted and stored but currently never hide an element.
# `Skin` only covers list-based views (`FromMesh = 1`, which would publish
# discrete model boundary entities, rejects) and `Transform`/`Remove`/`Smooth`/
# `Warp`/`SphericalRaise` reject model-based views — upstream in-place writers
# mutate `PViewDataList`. CutSphere and Isosurface register `RecurLevel`/
# `TargetError` like upstream;
# they steer the adaptive subdivision of *curved* elements, which our linear
# record catalog never produces, so they are stored without effect (the
# upstream linear-element path is reproduced exactly: corner levelset values
# interpolated along primary edges).

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
# `_PluginElement` kinds → the `_shape_functions` family used by the Newton
# inversion below (`:point` elements match by distance instead).
const _PLUGIN_WEIGHT_FAMILY=(quadrangle=:qua,hexahedron=:hex,prism=:pri,
                             pyramid=:pyr)

# Invert the bilinear/trilinear reference map through Newton iterations, the
# same `insideElement` + nodal shape evaluation the upstream `OctreePost`
# buckets perform for quadrangles, hexahedra, prisms and pyramids.
function _plugin_newton_weights(coords,kind::Symbol,p)
    family=_PLUGIN_WEIGHT_FAMILY[kind]
    ndim=kind===:quadrangle ? 2 : 3
    ξ=(0.0,0.0,0.0)
    scale=max(1.0,maximum(abs,coords))
    T=1e-9
    N=nothing
    for _ in 1:50
        N=collect(Float64,MeshPointLocation._shape_functions(
            family,ξ[1],ξ[2],ξ[3],"plugin"))
        resid=Float64[p[a]-LinearAlgebra.dot(N,view(coords,a,:)) for a in 1:3]
        LinearAlgebra.norm(resid,Inf)<=T*scale && break
        J=Matrix{Float64}(undef,3,ndim)
        for k in 1:ndim
            ξh=collect(ξ);ξh[k]+=1e-7
            Nh=collect(Float64,MeshPointLocation._shape_functions(
                family,ξh[1],ξh[2],ξh[3],"plugin"))
            for a in 1:3
                J[a,k]=LinearAlgebra.dot(Nh.-N,view(coords,a,:))/1e-7
            end
        end
        Δ=try
            (J'*J)\(J'*resid)
        catch
            return nothing
        end
        ξ=ξ.+ntuple(k->k<=ndim ? Δ[k] : 0.0,3)
        LinearAlgebra.norm(Δ,Inf)<=T && break
    end
    u,v,w=ξ
    inside=kind===:quadrangle ? (abs(u)<=1+T && abs(v)<=1+T) :
           kind===:hexahedron ? (abs(u)<=1+T && abs(v)<=1+T && abs(w)<=1+T) :
           kind===:prism ? (u>=-T && v>=-T && u+v<=1+T && abs(w)<=1+T) :
           (w>=-T && w<=1+T && abs(u)<=1-w+T && abs(v)<=1-w+T)
    inside || return nothing
    N=collect(Float64,MeshPointLocation._shape_functions(family,u,v,w,"plugin"))
    resid=Float64[p[a]-LinearAlgebra.dot(N,view(coords,a,:)) for a in 1:3]
    LinearAlgebra.norm(resid,Inf)<=T*scale || return nothing
    return Tuple(N)
end

function _plugin_local_weights(coords,kind::Symbol,p)
    T=1e-9
    a=_plugin_col(coords,1)
    if kind===:point
        scale=max(1.0,_plugin_norm3(a))
        _plugin_norm3(_plugin_sub3(p,a))<=T*scale || return nothing
        return (1.0,)
    elseif haskey(_PLUGIN_WEIGHT_FAMILY,kind)
        return _plugin_newton_weights(coords,kind,p)
    elseif kind===:line
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

# ---- multi-step element source ---------------------------------------------------

# Element records covering every populated step at once: list views already
# pack all steps into each record, while model views merge their per-(step,
# data type) blocks by element tag — an element absent from a step contributes
# zeros, matching the ExtractElements merge convention. Returns
# `(elements, nsteps, times)` with the per-step physical times.
function _plugin_elements_all_steps(view,caller::AbstractString)
    if view.model_kind===nothing
        elements,steps=_plugin_list_elements(view,caller)
        return elements,steps,[_plugin_step_time(view,s) for s in 0:steps-1]
    end
    keys=_plugin_model_keys(view)
    isempty(keys) && throw(ArgumentError(
        "$caller: view holds no populated model data steps"))
    per_step=Vector{Dict{UInt64,_PluginElement}}(undef,length(keys))
    firsts=Dict{UInt64,_PluginElement}();ordered=UInt64[]
    nc=nothing
    for (i,key) in enumerate(keys)
        elements,key_nc=_plugin_model_elements(view,key,caller)
        if nc===nothing
            nc=key_nc
        else
            key_nc==nc || throw(ArgumentError(
                "$caller: model data component count changes across steps"))
        end
        table=Dict{UInt64,_PluginElement}()
        for element in elements
            table[element.etag]=element
            haskey(firsts,element.etag) && continue
            firsts[element.etag]=element
            push!(ordered,element.etag)
        end
        per_step[i]=table
    end
    merged=_PluginElement[]
    for etag in ordered
        proto=firsts[etag]
        arity=size(proto.vcoords,2)
        vals=Vector{Float64}(undef,arity*nc*length(keys))
        for i in 1:length(keys)
            element=get(per_step[i],etag,nothing)
            covered=element!==nothing && size(element.vcoords,2)==arity
            for node in 1:arity,comp in 1:nc
                vals[(i-1)*arity*nc+(node-1)*nc+comp]=
                    covered ? element.vals[(node-1)*nc+comp] : 0.0
            end
        end
        push!(merged,_PluginElement(proto.kind,proto.dim,proto.geo,
                                    proto.vcoords,vals,nc,etag))
    end
    return merged,length(keys),[_plugin_step_time(view,key[1]) for key in keys]
end

# ---- Gradient ------------------------------------------------------------------

# Reference coordinates of each linear family's primary nodes — the nodal
# basis convention of `MeshFunctionSpaces._first_order_*` in MSH node order.
# Points carry no reference frame (no definable gradient) and are skipped.
const _PLUGIN_REFERENCE_NODES=(
    line=((-1.0,0.0,0.0),(1.0,0.0,0.0)),
    triangle=((0.0,0.0,0.0),(1.0,0.0,0.0),(0.0,1.0,0.0)),
    quadrangle=((-1.0,-1.0,0.0),(1.0,-1.0,0.0),(1.0,1.0,0.0),(-1.0,1.0,0.0)),
    tetrahedron=((0.0,0.0,0.0),(1.0,0.0,0.0),(0.0,1.0,0.0),(0.0,0.0,1.0)),
    hexahedron=((-1.0,-1.0,-1.0),(1.0,-1.0,-1.0),(1.0,1.0,-1.0),
                (-1.0,1.0,-1.0),(-1.0,-1.0,1.0),(1.0,-1.0,1.0),
                (1.0,1.0,1.0),(-1.0,1.0,1.0)),
    prism=((0.0,0.0,-1.0),(1.0,0.0,-1.0),(0.0,1.0,-1.0),
           (0.0,0.0,1.0),(1.0,0.0,1.0),(0.0,1.0,1.0)),
    pyramid=((-1.0,-1.0,0.0),(1.0,-1.0,0.0),(1.0,1.0,0.0),(-1.0,1.0,0.0),
             (0.0,0.0,1.0)))

# `_PluginElement.kind` → the `MeshFunctionSpaces` family symbol.
const _PLUGIN_KIND_FAMILY=(point=:pnt,line=:lin,triangle=:tri,
                           quadrangle=:qua,tetrahedron=:tet,hexahedron=:hex,
                           prism=:pri,pyramid=:pyr)

# Solve `(JᵀJ) h = r` for the assembled Jacobian columns `jcol` (`jcol[a][k]`
# holds ∂x_a/∂ξ_k). Returns the dim-tuple `h`, or `nothing` when the linear
# map is degenerate — the Gram determinant vanishes relative to its trace,
# which is dimensionally consistent for every dim.
function _plugin_metric_solve(jcol::NTuple{3,Vector{Float64}},
                              r::Vector{Float64},dim::Int)
    if dim==1
        g11=0.0
        for a in 1:3;g11=muladd(jcol[a][1],jcol[a][1],g11);end
        g11>0 || return nothing
        return (r[1]/g11,)
    elseif dim==2
        g11=g12=g22=0.0
        for a in 1:3
            g11=muladd(jcol[a][1],jcol[a][1],g11)
            g12=muladd(jcol[a][1],jcol[a][2],g12)
            g22=muladd(jcol[a][2],jcol[a][2],g22)
        end
        det=g11*g22-g12*g12;scale=g11+g22
        det>eps(Float64)*scale*scale || return nothing
        return (muladd(g22,r[1],-g12*r[2])/det,
                muladd(g11,r[2],-g12*r[1])/det)
    end
    g11=g12=g13=g22=g23=g33=0.0
    for a in 1:3
        ja=jcol[a]
        g11=muladd(ja[1],ja[1],g11);g12=muladd(ja[1],ja[2],g12)
        g13=muladd(ja[1],ja[3],g13);g22=muladd(ja[2],ja[2],g22)
        g23=muladd(ja[2],ja[3],g23);g33=muladd(ja[3],ja[3],g33)
    end
    c11=g22*g33-g23*g23;c12=g13*g23-g12*g33;c13=g12*g23-g13*g22
    det=muladd(g11,c11,muladd(g12,c12,g13*c13))
    tr=g11+g22+g33
    det>eps(Float64)*tr*tr*tr || return nothing
    c22=g11*g33-g13*g13;c23=g12*g13-g11*g23;c33=g11*g22-g12*g12
    return (muladd(c11,r[1],muladd(c12,r[2],c13*r[3]))/det,
            muladd(c12,r[1],muladd(c22,r[2],c23*r[3]))/det,
            muladd(c13,r[1],muladd(c23,r[2],c33*r[3]))/det)
end

# Spatial gradient of one scalar step at every value node,
# `∇u = J (JᵀJ)⁻¹ ∇ξu` evaluated at the node's reference location: the exact
# tangential gradient for dim<3 elements and the full gradient for volume
# elements. `nothing` when the element has no linear gradient — non-scalar
# data, value nodes outnumbering the primary map (high-order data), or a
# degenerate map.
function _plugin_element_gradient(element::_PluginElement,step::Int)
    nodes=get(_PLUGIN_REFERENCE_NODES,element.kind,nothing)
    nodes===nothing && return nothing
    element.ncomp==1 || return nothing
    arity=size(element.geo,2)
    size(element.vcoords,2)==arity || return nothing
    family=_PLUGIN_KIND_FAMILY[element.kind];dim=element.dim
    out=Vector{Float64}(undef,3arity)
    for node in 1:arity
        u,v,w=nodes[node]
        dN=MeshFunctionSpaces._first_order_gradients(
            Val(family),u,v,w,"API.plugin Gradient",node)
        jcols=ntuple(_->zeros(dim),3)
        r=zeros(dim)
        for i in 1:arity
            gi=dN[i];vi=element.vals[(step-1)*arity+i]
            for k in 1:dim
                r[k]=muladd(vi,gi[k],r[k])
                for a in 1:3
                    jcols[a][k]=muladd(element.geo[a,i],gi[k],jcols[a][k])
                end
            end
        end
        h=_plugin_metric_solve(jcols,r,dim)
        h===nothing && return nothing
        for a in 1:3
            g=0.0
            for k in 1:dim;g=muladd(jcols[a][k],h[k],g);end
            out[(node-1)*3+a]=g
        end
    end
    return out
end

function _plugin_exec_gradient(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    tag,view=_plugin_view(numbers,caller)
    elements,steps,times=_plugin_elements_all_steps(view,caller)
    out=_PluginElement[]
    for element in elements
        arity=size(element.geo,2)
        values=Vector{Float64}(undef,3arity*steps)
        ok=true
        for step in 1:steps
            part=_plugin_element_gradient(element,step)
            part===nothing && (ok=false;break)
            copyto!(values,(step-1)*3arity+1,part,1,3arity)
        end
        ok || continue
        push!(out,_PluginElement(element.kind,element.dim,element.geo,
                                 element.geo,values,3,element.etag))
    end
    _plugin_extract_emit(view.name*"_Gradient",out,-1,steps,times,caller)
    return nothing
end

# ---- CutPlane ------------------------------------------------------------------

# Primary-node edges per element kind in MSH vertex order (the linear pyramid
# fans its apex edges to node 5) — the levelset section walks them for zero
# crossings.
const _PLUGIN_CUT_EDGES=(
    point=NTuple{2,Int}[],
    line=((1,2),),
    triangle=((1,2),(2,3),(3,1)),
    quadrangle=((1,2),(2,3),(3,4),(4,1)),
    tetrahedron=((1,2),(1,3),(1,4),(2,3),(2,4),(3,4)),
    hexahedron=((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),(7,8),(8,5),
                (1,5),(2,6),(3,7),(4,8)),
    prism=((1,2),(2,3),(3,1),(4,5),(5,6),(6,4),(1,4),(2,5),(3,6)),
    pyramid=((1,2),(2,3),(3,4),(4,1),(1,5),(2,5),(3,5),(4,5)))

function _plugin_exec_cutplane(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    pa=_mesh_query_coordinate(numbers["A"],caller,"A")
    pb=_mesh_query_coordinate(numbers["B"],caller,"B")
    pc=_mesh_query_coordinate(numbers["C"],caller,"C")
    pd=_mesh_query_coordinate(numbers["D"],caller,"D")
    normal=(pa,pb,pc)
    _plugin_norm3(normal)>0 || throw(ArgumentError(
        "$caller: the plane normal (A,B,C) must be nonzero"))
    volume=_mesh_query_integer(Int(numbers["ExtractVolume"]),
                               caller,"ExtractVolume")
    tag,view=_plugin_view(numbers,caller)
    elements,steps,times=_plugin_elements_all_steps(view,caller)
    out=_PluginElement[]
    for element in elements
        ga=size(element.geo,2);va=size(element.vcoords,2);nc=element.ncomp
        # Plain `x*A + y*B + z*C + D` order — the fused-muladd chain rounds
        # differently (observable as `-0.0` on on-plane vertices).
        ls=ntuple(i->element.geo[1,i]*pa + element.geo[2,i]*pb +
                     element.geo[3,i]*pc + pd,ga)
        if !iszero(volume)
            # Whole-element extraction keeps elements whose corner levelset
            # values all lie on the selected side; upstream re-emits them
            # through the iso-simplex decomposition (a kept quadrangle is two
            # `ST` records, a kept hexahedron six `SS`) — only linear records
            # round-trip.
            if (volume>0 ? all(l->l>=0.0,ls) : all(l->l<=0.0,ls)) && va==ga
                _plugin_push_subsimplices!(out,element,steps)
            end
            continue
        end
        # The levelset section machinery (sub-simplex edge scan, `1e-12`
        # dedup, dimension gates, butterfly reorder, MAP orientation) is
        # shared with `Isosurface` — upstream emits the same records for an
        # all-on-plane element (a line gives `SP` at its first node, a
        # triangle `SL` on its first edge, a tet `ST` on its first face),
        # not the element verbatim.
        for (kind,points,weights) in _plugin_levelset_sections(element,ls)
            _plugin_record_type(kind,nc)===nothing && continue
            arity=length(points)
            geo=Matrix{Float64}(undef,3,arity)
            for i in 1:arity,axis in 1:3
                geo[axis,i]=points[i][axis]
            end
            vals=Vector{Float64}(undef,arity*nc*steps)
            for step in 1:steps,i in 1:arity,comp in 1:nc
                # Upstream value interpolation: v0 + c*(v1-v0).
                a=weights[i][1][1];b=weights[i][2][1];c=weights[i][2][2]
                v0=element.vals[(step-1)*va*nc+(a-1)*nc+comp]
                vals[(step-1)*arity*nc+(i-1)*nc+comp]=
                    v0+c*(element.vals[(step-1)*va*nc+(b-1)*nc+comp]-v0)
            end
            push!(out,_PluginElement(kind,_PLUGIN_KIND_DIM[kind],geo,geo,
                                     vals,nc,element.etag))
        end
    end
    _plugin_extract_emit(view.name*"_Levelset",out,-1,steps,times,caller)
    return nothing
end

# ---- shared view-tag option resolution ----------------------------------------------

# A secondary view-valued numeric option (`OtherView`, `ViewX`, ...): a value
# < 0 names the most recent session view like `View < 0`; a nonnegative value
# must be an existing integer tag. (Upstream's `getView(i, v)` contract.)
function _plugin_view_tag(numbers::Dict{String,Float64},opt::AbstractString,
                          caller::AbstractString)
    raw=numbers[opt]
    if raw<0
        isempty(SESSION_VIEWS) && throw(ArgumentError(
            "$caller: no views exist to run the plugin on"))
        tag=maximum(keys(SESSION_VIEWS))
    else
        raw==floor(raw) || throw(ArgumentError(
            "$caller: $opt tag must be an integer"))
        tag=Int(raw)
        haskey(SESSION_VIEWS,tag) || throw(ArgumentError(
            "$caller: unknown view tag $tag"))
    end
    return tag,_session_view_data(tag,caller)
end

# ---- CutSphere ---------------------------------------------------------------------------

function _plugin_exec_cutsphere(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    xc=_mesh_query_coordinate(numbers["Xc"],caller,"Xc")
    yc=_mesh_query_coordinate(numbers["Yc"],caller,"Yc")
    zc=_mesh_query_coordinate(numbers["Zc"],caller,"Zc")
    r=_mesh_query_coordinate(numbers["R"],caller,"R")
    volume=_mesh_query_integer(Int(numbers["ExtractVolume"]),
                               caller,"ExtractVolume")
    tag,view=_plugin_view(numbers,caller)
    elements,steps,times=_plugin_elements_all_steps(view,caller)
    center=(xc,yc,zc)
    out=_PluginElement[]
    for element in elements
        ga=size(element.geo,2);va=size(element.vcoords,2);nc=element.ncomp
        # Plain left-associative arithmetic like the CutPlane levelset —
        # fused muladd rounds differently on exact-sphere vertices.
        ls=ntuple(ga) do i
            dx=element.geo[1,i]-xc;dy=element.geo[2,i]-yc
            dz=element.geo[3,i]-zc
            return dx*dx+dy*dy+dz*dz-r*r
        end
        if !iszero(volume)
            # Whole-element extraction: `ExtractVolume > 0` keeps elements
            # fully outside the sphere (all ls >= 0), `< 0` fully inside —
            # same convention as CutPlane. Kept elements re-emit through
            # the iso-simplex decomposition like upstream.
            if (volume>0 ? all(l->l>=0.0,ls) : all(l->l<=0.0,ls)) && va==ga
                _plugin_push_subsimplices!(out,element,steps)
            end
            continue
        end
        # Shared levelset section machinery (see `_plugin_levelset_sections`)
        # with the sphere's analytic field gradient — the emitted facet's
        # geometric normal follows the outward radial direction
        # `centroid - centre` (the gradient of `|x-c|²-R²` at the section),
        # not the per-simplex affine fit.
        radial(points)=begin
            n=length(points)
            (sum(p->p[1],points)/n-xc,
             sum(p->p[2],points)/n-yc,
             sum(p->p[3],points)/n-zc)
        end
        for (kind,spoints,weights) in _plugin_levelset_sections(element,ls;
                                                              orient=radial)
            _plugin_record_type(kind,nc)===nothing && continue
            arity=length(spoints)
            geo=Matrix{Float64}(undef,3,arity)
            for i in 1:arity,axis in 1:3
                geo[axis,i]=spoints[i][axis]
            end
            vals=Vector{Float64}(undef,arity*nc*steps)
            for step in 1:steps,i in 1:arity,comp in 1:nc
                a=weights[i][1][1];b=weights[i][2][1];c=weights[i][2][2]
                v0=element.vals[(step-1)*va*nc+(a-1)*nc+comp]
                vals[(step-1)*arity*nc+(i-1)*nc+comp]=
                    v0+c*(element.vals[(step-1)*va*nc+(b-1)*nc+comp]-v0)
            end
            push!(out,_PluginElement(kind,_PLUGIN_KIND_DIM[kind],geo,geo,
                                     vals,nc,element.etag))
        end
    end
    _plugin_extract_emit(view.name*"_Levelset",out,-1,steps,times,caller)
    return nothing
end

# ---- Isosurface -----------------------------------------------------------------

# `Plugin(Isosurface)` extracts the `scalar value == Value` level set of the
# input view's own data: the levelset at each primary node is the upstream
# `ComputeScalarRep` of that node's components minus `Value` (scalars pass
# through, vectors their norm, tensors the Von-Mises invariant). The section
# geometry is value-dependent — upstream's "moving" levelset — so like
# upstream the plugin emits one list-based view per populated step, always
# named `<name>_Levelset_<step>` (verified on the oracle — a single-step view
# produces `<name>_Levelset_0`).
#
# Non-simplex elements pass through upstream's `getSimplexDec` decomposition
# (quadrangle → two triangles, hexahedron → six tetrahedra, prism → three,
# pyramid → two); every sub-simplex scans its edges for
# `levels[a]*levels[b] <= 0` crossings interpolated linearly along the edge
# and deduplicated on first-seen coordinates. The surviving crossings emit a
# point/line/triangle/quadrangle record — four crossings are reordered into
# a simple quad like upstream `reorderQuad` — with values linearly
# interpolated from the corner values, and section normals oppose the
# element's scalar gradient (upstream `MAP` orientation). `ExtractVolume`
# keeps the whole elements on the selected side verbatim — the
# CutPlane/CutSphere session convention — and elements lying entirely on the
# level set re-emit verbatim when they fit a list record type.
# `RecurLevel`/`TargetError` steer upstream's adaptive subdivision of curved
# elements; linear sections are exact without them, so they are stored
# without effect. Upstream's `OtherView`/`OtherTimeStep` second-source
# options are not registered: session views do not align element-wise
# across views.

# Upstream `getSimplexDec` sub-simplex node lists (1-based indices into the
# element's primary nodes): the point/line/triangle/tetrahedron kinds are
# already simplices.
const _PLUGIN_ISO_SIMPLICES=(
    point=((1,),),
    line=((1,2),),
    triangle=((1,2,3),),
    quadrangle=((1,2,3),(1,3,4)),
    tetrahedron=((1,2,3,4),),
    hexahedron=((1,2,4,8),(1,5,2,8),(2,5,6,8),
                (2,3,4,8),(2,7,3,8),(2,6,7,8)),
    prism=((1,2,3,5),(1,3,5,6),(1,4,5,6)),
    pyramid=((1,2,4,5),(2,3,4,5)))

# Upstream `exn` edge tables for the sub-simplex shapes, indexed by node
# count: a point scans no edges, then the line, triangle and tetrahedron
# edge orders.
const _PLUGIN_ISO_SUBEDGES=(NTuple{2,Int}[],
                          NTuple{2,Int}[(1,2)],
                          NTuple{2,Int}[(1,2),(1,3),(2,3)],
                          NTuple{2,Int}[(1,2),(1,3),(1,4),
                                        (2,3),(2,4),(3,4)])

# Upstream `gradSimplex`: the exact gradient of the linear field defined by
# the first four element nodes (Cramer solve; the zero vector on a
# degenerate configuration).
function _plugin_grad_simplex(coords,vals)
    j11=coords[1,2]-coords[1,1];j12=coords[2,2]-coords[2,1];j13=coords[3,2]-coords[3,1]
    j21=coords[1,3]-coords[1,1];j22=coords[2,3]-coords[2,1];j23=coords[3,3]-coords[3,1]
    j31=coords[1,4]-coords[1,1];j32=coords[2,4]-coords[2,1];j33=coords[3,4]-coords[3,1]
    det=j11*(j22*j33-j23*j32)-j12*(j21*j33-j23*j31)+j13*(j21*j32-j22*j31)
    det==0.0 && return (0.0,0.0,0.0)
    b1=vals[2]-vals[1];b2=vals[3]-vals[1];b3=vals[4]-vals[1]
    g1=(b1*(j22*j33-j23*j32)-j12*(b2*j33-j23*b3)+j13*(b2*j32-j22*b3))/det
    g2=(j11*(b2*j33-j23*b3)-b1*(j21*j33-j23*j31)+j13*(j21*b3-b2*j31))/det
    g3=(j11*(j22*b3-b2*j32)-j12*(j21*b3-b2*j31)+b1*(j21*j32-j22*j31))/det
    return (g1,g2,g3)
end

# `ExtractVolume`-style verbatim emission: upstream re-emits a kept element
# through the same `_PLUGIN_ISO_SIMPLICES` decomposition it cuts with, so a
# quadrangle publishes two `ST` records and a hexahedron six `SS` records —
# geometry and values are gathered, never interpolated.
const _PLUGIN_SIMPLEX_KIND=(:point,:line,:triangle,:tetrahedron)
function _plugin_push_subsimplices!(out,element,steps)
    va=size(element.vcoords,2);nc=element.ncomp
    for map in _PLUGIN_ISO_SIMPLICES[element.kind]
        kind=_PLUGIN_SIMPLEX_KIND[length(map)]
        _plugin_record_type(kind,nc)===nothing && continue
        arity=length(map)
        geo=Matrix{Float64}(undef,3,arity)
        for i in 1:arity,axis in 1:3
            geo[axis,i]=element.geo[axis,map[i]]
        end
        vals=Vector{Float64}(undef,arity*nc*steps)
        for step in 1:steps,i in 1:arity,comp in 1:nc
            vals[(step-1)*arity*nc+(i-1)*nc+comp]=
                element.vals[(step-1)*va*nc+(map[i]-1)*nc+comp]
        end
        push!(out,_PluginElement(kind,_PLUGIN_KIND_DIM[kind],geo,geo,
                                 vals,nc,element.etag))
    end
    return out
end

# One element's levelset section records as `(kind, points, weights)`:
# every sub-simplex (`_PLUGIN_ISO_SIMPLICES`) scans its edges for
# `levels[a]*levels[b] <= 0` crossings — the `la==lb` degenerate arm emits
# the edge's first vertex — and upstream `removeIdenticalNodes` keeps the
# first crossing within a per-coordinate `1e-12` box. `levels` holds the
# element-node levelset values (scalar representation minus the level);
# the constant level shift leaves the orientation gradient unchanged.
function _plugin_levelset_sections(element,levels;orient=nothing)
    kind=element.kind;geo=element.geo
    numedges=length(_PLUGIN_CUT_EDGES[kind])
    out=Tuple{Symbol,Vector{NTuple{3,Float64}},
              Vector{Vector{Tuple{Int,Float64}}}}[]
    for simplex in _PLUGIN_ISO_SIMPLICES[kind]
        points=NTuple{3,Float64}[];weights=Vector{Tuple{Int,Float64}}[]
        for (ea,eb) in _PLUGIN_ISO_SUBEDGES[length(simplex)]
            a=simplex[ea];b=simplex[eb]
            la=levels[a];lb=levels[b]
            la*lb<=0.0 || continue
            c=la==lb ? 0.0 : la/(la-lb)
            # Upstream `InterpolateIso`: coef*(b-a)+a (no fused ops).
            point=(c*(geo[1,b]-geo[1,a])+geo[1,a],
                   c*(geo[2,b]-geo[2,a])+geo[2,a],
                   c*(geo[3,b]-geo[3,a])+geo[3,a])
            # Upstream `removeIdenticalNodes`: the first-seen crossing wins.
            found=false
            for seen in points
                if abs(seen[1]-point[1])<1e-12 &&
                    abs(seen[2]-point[2])<1e-12 &&
                    abs(seen[3]-point[3])<1e-12
                    found=true;break
                end
            end
            found && continue
            push!(points,point);push!(weights,[(a,1.0-c),(b,c)])
        end
        np=length(points)
        # Upstream validity gates: 3-D elements emit only surface sections,
        # 2-D only lines, and more than four crossings cannot order into a
        # record (linear cases never exceed four after deduplication).
        numedges>4 && np<3 && continue
        numedges>1 && np<2 && continue
        (np<1 || np>4) && continue
        if np==4
            # Upstream `reorderQuad`: the edge-scan order reads a
            # butterfly; swapping the last two vertices orders the quad.
            points[3],points[4]=points[4],points[3]
            weights[3],weights[4]=weights[4],weights[3]
        end
        if np==3 || np==4
            # MAP orientation: the emitted section's normal opposes the
            # scalar gradient — upstream flips the polygon when the
            # (p3-p1)x(p2-p1) normal already points along it. The gradient
            # is the linear solve on THIS sub-simplex's four nodes (only a
            # tetrahedron simplex can reach np >= 3), not the parent's
            # first four — for a decomposed hex those can be coplanar. An
            # `orient` callback overrides the gradient for fields whose
            # orientation direction is not the affine fit (CutSphere's
            # radial normal).
            normal=_plugin_cross3(_plugin_sub3(points[3],points[1]),
                                  _plugin_sub3(points[2],points[1]))
            grad=orient!==nothing ? orient(points) :
                length(simplex)==4 ?
                _plugin_grad_simplex(
                    [geo[axis,simplex[i]] for axis in 1:3,i in 1:4],
                    (levels[simplex[1]],levels[simplex[2]],
                     levels[simplex[3]],levels[simplex[4]])) :
                (0.0,0.0,0.0)
            _plugin_dot3(grad,normal)>0.0 &&
                (reverse!(points);reverse!(weights))
        end
        push!(out,(np==1 ? :point : np==2 ? :line :
                   np==3 ? :triangle : :quadrangle,points,weights))
    end
    return out
end

# Append one element's single-step levelset section records to `out`.
function _plugin_iso_sections!(out,element,levels,step,nsteps)
    nc=element.ncomp;va=size(element.vcoords,2)
    for (outkind,points,weights) in _plugin_levelset_sections(element,levels)
        _plugin_record_type(outkind,nc)===nothing && continue
        np=length(points)
        section=Matrix{Float64}(undef,3,np)
        for i in 1:np,axis in 1:3
            section[axis,i]=points[i][axis]
        end
        vals=zeros(Float64,np*nc*nsteps)
        for i in 1:np,comp in 1:nc
            # Upstream value interpolation: v0 + c*(v1-v0) with the edge
            # coefficient — weights[i] is [(a,1-c),(b,c)] by construction.
            a=weights[i][1][1];b=weights[i][2][1];c=weights[i][2][2]
            v0=element.vals[step*va*nc+(a-1)*nc+comp]
            vals[step*np*nc+(i-1)*nc+comp]=
                v0+c*(element.vals[step*va*nc+(b-1)*nc+comp]-v0)
        end
        push!(out,_PluginElement(outkind,_PLUGIN_KIND_DIM[outkind],section,
                                 section,vals,nc,element.etag))
    end
    return out
end

function _plugin_exec_isosurface(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    value=_mesh_query_coordinate(numbers["Value"],caller,"Value")
    volume=_mesh_query_integer(Int(numbers["ExtractVolume"]),
                               caller,"ExtractVolume")
    tag,view=_plugin_view(numbers,caller)
    elements,steps,times=_plugin_elements_all_steps(view,caller)
    for step in 0:steps-1
        out=_PluginElement[]
        for element in elements
            ga=size(element.geo,2);va=size(element.vcoords,2);nc=element.ncomp
            levels=ntuple(ga) do node
                _plugin_scalar_rep(element.vals,step*va*nc+(node-1)*nc,nc)-
                    value
            end
            verbatim=va==ga && _plugin_record_type(element.kind,nc)!==nothing
            if !iszero(volume)
                # Whole-element extraction like CutPlane: the element is
                # kept verbatim when every corner levelset value lies on
                # the selected side.
                (volume>0 ? all(l->l>=0.0,levels) :
                            all(l->l<=0.0,levels)) && verbatim &&
                    push!(out,element)
                continue
            end
            all(l->l==0.0,levels) &&
                (verbatim && push!(out,element);continue)
            _plugin_iso_sections!(out,element,levels,step,steps)
        end
        name=view.name*"_Levelset_"*string(step)
        _plugin_extract_emit(name,out,step,steps,times,caller)
    end
    return nothing
end

# ---- Warp --------------------------------------------------------------------------------

# Smoothed facet normals keyed on node position — the upstream
# `smooth_normals` container: `add` folds a unit normal into the first cluster
# whose angle to it stays under the tolerance (running mean capped at 100
# contributions, order-dependent like upstream); `get` returns the first
# cluster within tolerance of the caller-supplied seed (the element's own
# facet normal, so the "same side" cluster wins at sharp creases), or the seed
# itself when nothing is close enough.
struct _PluginNormalCluster
    n::NTuple{3,Float64}
    nb::Int
end

# Angle between two vectors in degrees, `atan2(|a×b|, a·b)` like upstream.
function _plugin_angle_deg(a,b)
    na=_plugin_norm3(a);nb=_plugin_norm3(b)
    (na==0.0 || nb==0.0) && return 0.0
    return atand(_plugin_norm3(_plugin_cross3(a,b))/(na*nb),
                 _plugin_dot3(a,b)/(na*nb))
end

function _plugin_normals_add!(clusters,p,n,tol)
    bucket=get!(clusters,p,_PluginNormalCluster[])
    for (i,cluster) in enumerate(bucket)
        if tol>=180.0 || abs(_plugin_angle_deg(cluster.n,n))<tol
            cluster.nb<100 || return nothing
            c1=cluster.nb/(cluster.nb+1);c2=1/(cluster.nb+1)
            bucket[i]=_PluginNormalCluster(
                (c1*cluster.n[1]+c2*n[1],
                 c1*cluster.n[2]+c2*n[2],
                 c1*cluster.n[3]+c2*n[3]),cluster.nb+1)
            return nothing
        end
    end
    push!(bucket,_PluginNormalCluster(n,0))
    return nothing
end

function _plugin_normals_get(clusters,p,seed,tol)
    bucket=get(clusters,p,nothing)
    bucket===nothing && return seed
    for cluster in bucket
        abs(_plugin_angle_deg(cluster.n,seed))<tol && return cluster.n
    end
    return seed
end

# Unit facet normal `normal3points` of the element's first three nodes read at
# position `base` in a list record; `nothing` on a degenerate element.
function _plugin_facet_normal(data,base::Int,arity::Int)
    a=(data[base+1],data[base+arity+1],data[base+2arity+1])
    b=(data[base+2],data[base+arity+2],data[base+2arity+2])
    c=(data[base+3],data[base+arity+3],data[base+2arity+3])
    n=_plugin_cross3(_plugin_sub3(b,a),_plugin_sub3(c,a))
    s=_plugin_norm3(n)
    s>0 || return nothing
    return (n[1]/s,n[2]/s,n[3]/s)
end

function _plugin_exec_warp(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    tag,view=_plugin_view(numbers,caller)
    view.model_kind===nothing || throw(ArgumentError(
        "$caller: Warp rewrites list data only; model views are unsupported"))
    factor=_mesh_query_coordinate(numbers["Factor"],caller,"Factor")
    timestep=_mesh_query_integer(Int(numbers["TimeStep"]),caller,"TimeStep")
    tol=_mesh_query_coordinate(numbers["SmoothingAngle"],caller,"SmoothingAngle")
    # Upstream (`Warp.cpp`): `otherView < 0` first takes the `View` value, so
    # the smoothed-normal field is built only when the resolved tag is still
    # negative — i.e. both `View` and `OtherView` were left < 0 (current view
    # warped by its own surface normals times its scalar values).
    raw_other=numbers["OtherView"]
    raw_other<0 && (raw_other=numbers["View"])
    if raw_other<0
        normals_mode=true;otag=tag;oview=view
    else
        raw_other==floor(raw_other) || throw(ArgumentError(
            "$caller: OtherView tag must be an integer"))
        otag=Int(raw_other)
        haskey(SESSION_VIEWS,otag) || throw(ArgumentError(
            "$caller: unknown view tag $otag"))
        oview=_session_view_data(otag,caller)
        normals_mode=false
        oview.model_kind===nothing || throw(ArgumentError(
            "$caller: Warp requires a list-based OtherView"))
    end
    steps=_plugin_list_steps(view,caller)
    osteps=_plugin_list_steps(oview,caller)
    # Upstream sanity check: equal entity/element counts. List views carry no
    # entity records, so the check is the flattened element count; a per-node
    # arity or component mismatch below then just leaves `val` at zero like
    # upstream's `getNumNodes`/`getNumComponents` guards.
    positions=Tuple{String,Int}[]
    for type in view.list_types
        for element in 0:view.list_counts[type]-1
            push!(positions,(type,element))
        end
    end
    oelements=normals_mode ? _PluginElement[] :
        _plugin_list_elements(oview,caller)[1]
    normals_mode || length(oelements)==length(positions) || begin
        _native_diagnostic_emit_locked(:error,"Incompatible views")
        return nothing
    end
    0<=timestep<=osteps-1 || begin
        _native_diagnostic_emit_locked(:error,
            "Invalid TimeStep ($timestep) in View[$otag]")
        return nothing
    end
    # Precompute the smoothed normal field over the triangle/quadrangle
    # records (upstream `numEdges == 3 || 4`) before any coordinate moves.
    clusters=Dict{NTuple{3,Float64},Vector{_PluginNormalCluster}}()
    if normals_mode
        for type in view.list_types
            layout=get(Post._POS_RECORD_TYPES,type,nothing)
            layout===nothing && throw(ArgumentError(
                "$caller: unsupported list data type $type"))
            kind,arity,ncomp=layout
            kind in (:triangle,:quadrangle) || continue
            count=view.list_counts[type];data=view.list_data[type]
            width=3arity+arity*ncomp*steps
            for element in 0:count-1
                base=element*width
                n=_plugin_facet_normal(data,base,arity)
                n===nothing && continue
                for node in 1:arity
                    _plugin_normals_add!(clusters,
                        (data[base+node],data[base+arity+node],
                         data[base+2arity+node]),n,tol)
                end
            end
        end
    end
    for (pos,(type,element)) in enumerate(positions)
        layout=get(Post._POS_RECORD_TYPES,type,nothing)
        layout===nothing && throw(ArgumentError(
            "$caller: unsupported list data type $type"))
        kind,arity,ncomp=layout
        arity>=2 || continue # upstream skips 1-node records
        data=view.list_data[type]
        width=3arity+arity*ncomp*steps
        base=element*width
        oe=normals_mode ? nothing : oelements[pos]
        for step in 1:steps
            # The facet-normal seed reads possibly-already-warped coordinates
            # like the upstream per-step element scan.
            nseed=normals_mode ? _plugin_facet_normal(data,base,arity) :
                nothing
            for node in 1:arity
                mult=1.0
                if normals_mode
                    kind in (:triangle,:quadrangle) || continue
                    nseed===nothing && continue
                    p=(data[base+node],data[base+arity+node],
                       data[base+2arity+node])
                    val=_plugin_normals_get(clusters,p,nseed,tol)
                    mult=_plugin_scalar_rep(data,
                        base+3arity+(step-1)*arity*ncomp+(node-1)*ncomp,
                        ncomp)
                else
                    # The same flattened (element, node) slot of `OtherView`
                    # at `TimeStep`; only matching-arity 3-component records
                    # displace (upstream's `getNumComponents`/`getNumNodes`
                    # guards leave `val` at zero otherwise).
                    oa=size(oe.vcoords,2);onc=oe.ncomp
                    (oa==arity && onc==3) || continue
                    val=(oe.vals[timestep*oa*onc+(node-1)*onc+1],
                         oe.vals[timestep*oa*onc+(node-1)*onc+2],
                         oe.vals[timestep*oa*onc+(node-1)*onc+3])
                end
                for axis in 1:3
                    data[base+(axis-1)*arity+node]+=factor*mult*val[axis]
                end
            end
        end
    end
    return nothing
end

# ---- MathEval ------------------------------------------------------------------------------

import ..SizeField
import ..Mesh2D
import ..Mesh3D
import LinearAlgebra

# Per-node variable buffers for `v0`..`v8` (the view being processed) and
# `w0`..`w8` (`OtherView` at `OtherTimeStep` or the same step); bound into the
# `SizeField` MathEval evaluator through `F<tag>` references — `v<d>` parses
# as tag `d`, `w<d>` as tag `d + 16`.
mutable struct _PluginEvalVars
    v::Vector{Float64}
    w::Vector{Float64}
end

struct _PluginEvalField <: SizeField.AbstractField
    comp::Int
    other::Bool
end

# Mutable carrier for upstream `Time`/`TimeStep` expression variables in
# ModifyComponents — bound to private field tags F30/F31 and refreshed once
# per visited (step) so every bound expression sees the current values.
mutable struct _PluginEvalScalar <: SizeField.AbstractField
    value::Float64
end

@inline function SizeField.field_value(f::_PluginEvalField,x,y,z,
                                       vars::_PluginEvalVars)
    return f.other ? vars.w[f.comp] : vars.v[f.comp]
end

@inline function SizeField.field_value(f::_PluginEvalScalar,x,y,z,
                                       vars::_PluginEvalVars)
    return f.value
end

@inline _plugin_ident_char(c::Char)=
    ('a'<=c<='z')||('A'<=c<='Z')||('0'<=c<='9')||c=='_'

# Rewrite `v<d>`/`w<d>` variable names into the `F<tag>` references the
# `SizeField` expression compiler resolves (`v<d>` → `F<d>`, `w<d>` →
# `F<d + 16>`). `F<digits>` itself is not an upstream MathEval variable and is
# rejected so users cannot shadow the mapped tags; every other identifier
# falls through to the parser's own validation (function names, x, y, z,
# Pi, e). With `time_vars`, ModifyComponents' extra `Time`/`TimeStep`
# variables map to the private tags F30/F31 bound to `_PluginEvalScalar`
# carriers.
function _plugin_matheval_rewrite(expr::AbstractString,caller::AbstractString;
                                  time_vars::Bool=false)
    out=IOBuffer();i=firstindex(expr);n=lastindex(expr)
    while i<=n
        c=expr[i]
        if _plugin_ident_char(c) && !(c in ('+','-'))
            j=i
            while j<=n && _plugin_ident_char(expr[j]);j=nextind(expr,j);end
            ident=expr[i:prevind(expr,j)]
            if time_vars && ident=="Time"
                write(out,"F30");i=j;continue
            elseif time_vars && ident=="TimeStep"
                write(out,"F31");i=j;continue
            end
            if ident=="v" || ident=="w" || (length(ident)==2 && ident[1] in ('v','w') &&
                                          _plugin_ident_char(ident[2]))
                if length(ident)==2 && '0'<=ident[2]<='8'
                    d=ident[2]-'0'
                    write(out,'F');write(out,ident[1]=='v' ? string(d) : string(16+d))
                    i=j;continue
                elseif length(ident)>1 && all(ch->'0'<=ch<='9',ident[2:end])
                    throw(ArgumentError("$caller: unknown identifier $ident"))
                end
            elseif length(ident)>=2 && ident[1]=='F' &&
                all(ch->'0'<=ch<='9',ident[2:end])
                throw(ArgumentError("$caller: unknown identifier $ident"))
            end
            write(out,ident);i=j;continue
        end
        write(out,c);i=nextind(expr,i)
    end
    return String(take!(out))
end

# Component values of `view` at `stepnum` (0-based step index) evaluated at
# point `p` through the containing element's barycentric weights — the
# `OctreePost::search*` interpolation upstream uses for off-grid `w` values.
# `dimension > 0` restricts the searched element dimension like the upstream
# `search*(..., dimension)` argument. `nothing` when `p` lies outside or the
# step is unpopulated.
function _plugin_probe_at(view,p,stepnum::Int,caller::AbstractString;
                          dimension::Int=-1)
    if view.model_kind===nothing
        elements,osteps=_plugin_list_elements(view,caller)
        0<=stepnum<osteps || return nothing
        pos=stepnum+1
        for element in elements
            element.dim<1 && continue
            dimension>0 && element.dim!=dimension && continue
            N=_plugin_local_weights(element.geo,element.kind,p)
            N===nothing && continue
            oa=size(element.vcoords,2);onc=element.ncomp
            acc=zeros(onc)
            for node in 1:min(oa,length(N)),comp in 1:onc
                acc[comp]+=N[node]*
                    element.vals[(pos-1)*oa*onc+(node-1)*onc+comp]
            end
            return acc
        end
        return nothing
    end
    values,stepnums=_plugin_probe_model(view,p,caller)
    pos=findfirst(==(stepnum),stepnums)
    pos===nothing && return nothing
    return values[pos]
end

# Owning entity `(dim, tag)` of dense element `index` from the mesh
# classification, for `PhysicalGroup` filtering on model views — the
# `ge->physicals` lookup of the upstream entity-bucket scan. `nothing` when
# the classification is absent or has no owners for the family.
function _plugin_element_entity(mesh,class,index::Int)
    class===nothing && return nothing
    if mesh isa Mesh
        tri_off,tet_off,_=_mesh_element_offsets(mesh)
        if index<=tri_off
            index<=length(class.seg_entities) || return nothing
            return (1,Int(class.seg_entities[index]))
        elseif index<=tet_off
            col=index-tri_off
            col<=length(class.tri_entities) || return nothing
            return (2,Int(class.tri_entities[col]))
        end
        col=index-tet_off
        col<=length(class.tet_entities) || return nothing
        return (3,Int(class.tet_entities[col]))
    elseif mesh isa MixedMesh
        offset=0;family=Dict{Int,Int}()
        for (msh,cells,_) in _cache_native_blocks(mesh)
            rank=get(family,msh,0)
            if offset<index<=offset+size(cells,2)
                owners=get(class.cell_entities,msh,Int32[])
                pos=rank+(index-offset)
                pos>length(owners) && return nothing
                return (Int(msh_dimension(msh)),Int(owners[pos]))
            end
            offset+=size(cells,2);family[msh]=rank+size(cells,2)
        end
    end
    return nothing
end

function _plugin_exec_matheval(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers;strings=spec.strings
    tag,view=_plugin_view(numbers,caller)
    # Upstream component rule: any of Expression3..8 → tensor view (9), else
    # any of Expression1..2 → vector view (3), else scalar (1); empty slots
    # inside the active range evaluate as "0".
    raw_exprs=[strings["Expression$i"] for i in 0:8]
    if any(!isempty,raw_exprs[4:9])
        ncomp2=9
        for i in 1:9;isempty(raw_exprs[i]) && (raw_exprs[i]="0");end
    elseif any(!isempty,raw_exprs[2:3])
        ncomp2=3
        for i in 1:3;isempty(raw_exprs[i]) && (raw_exprs[i]="0");end
    else
        ncomp2=1
    end
    exprs=raw_exprs[1:ncomp2]
    progs=[SizeField.parse_matheval(_plugin_matheval_rewrite(e,caller))
           for e in exprs]
    # `OtherView` defaults to the same view like upstream (`otherView = v1`).
    otag,oview=numbers["OtherView"]<0 ? (tag,view) :
        _plugin_view_tag(numbers,"OtherView",caller)
    timestep=_mesh_query_integer(Int(numbers["TimeStep"]),caller,"TimeStep")
    ots=_mesh_query_integer(Int(numbers["OtherTimeStep"]),
                            caller,"OtherTimeStep")
    forcerep=!iszero(numbers["ForceInterpolation"])
    pgroup=_mesh_query_integer(Int(numbers["PhysicalGroup"]),
                               caller,"PhysicalGroup")
    dimension=_mesh_query_integer(Int(numbers["Dimension"]),caller,"Dimension")
    -1<=dimension<=3 || throw(ArgumentError(
        "$caller: Dimension must be -1 or in 0:3"))
    elements,steps,times=_plugin_elements_all_steps(view,caller)
    stepnums=view.model_kind===nothing ? collect(0:steps-1) :
        [key[1] for key in _plugin_model_keys(view)]
    # Upstream rejects multi-mesh views outright.
    for (mview,mkeys,mname) in ((view,
            view.model_kind===nothing ? Tuple{Int,String}[] :
                _plugin_model_keys(view),"view"),
        (oview===view ? nothing : oview,
            (oview===view || oview.model_kind===nothing) ?
                Tuple{Int,String}[] : _plugin_model_keys(oview),"other view"))
        mview===nothing && continue
        slots=Set{Any}()
        for key in mkeys
            push!(slots,get(mview.model_contexts[key],:model_slot,nothing))
        end
        length(slots)>1 && begin
            _native_diagnostic_emit_locked(:error,
                "MathEval plugin cannot be applied to multi-mesh views")
            return nothing
        end
    end
    ostepnums=if oview===view
        stepnums
    elseif oview.model_kind===nothing
        collect(0:_plugin_list_steps(oview,caller)-1)
    else
        [key[1] for key in _plugin_model_keys(oview)]
    end
    if ots<0 && length(ostepnums)!=length(stepnums)
        _native_diagnostic_emit_locked(:error,
            "Number of time steps don't match: using step 0")
        ots=0
    elseif ots>0 && !(ots in ostepnums)
        _native_diagnostic_emit_locked(:error,
            "Invalid time step ($ots) in View[$otag]: using step 0 instead")
        ots=0
    end
    # `TimeStep` selects a stored step index: only an index past the last
    # populated step is "invalid" upstream (then all steps run); a legal but
    # unpopulated slot emits nothing.
    emit_positions=if timestep<0
        collect(1:steps)
    elseif !isempty(stepnums) && timestep>last(stepnums)
        _native_diagnostic_emit_locked(:error,
            "Invalid time step ($timestep) in View[$tag]: " *
            "using all steps instead")
        collect(1:steps)
    else
        pos=findfirst(==(timestep),stepnums)
        pos===nothing ? Int[] : [pos]
    end
    isempty(emit_positions) && begin
        _plugin_extract_emit(view.name*"_MathEval",_PluginElement[],-1,0,
                             Float64[],caller)
        return nothing
    end
    # The upstream check (`getNumEntities`/`getNumElements` equality) decides
    # whether `OtherView` reads at the same (entity, element, node) slots or
    # through the interpolation octree; `ForceInterpolation` always probes.
    same_grid=false
    oelements=_PluginElement[]
    if oview===view
        same_grid=true;oelements=elements
    elseif view.model_kind===nothing && oview.model_kind===nothing
        same_grid=view.list_types==oview.list_types &&
                  view.list_counts==oview.list_counts
        if same_grid
            oelements,_,_=_plugin_elements_all_steps(oview,caller)
        end
    elseif view.model_kind!==nothing && oview.model_kind!==nothing
        merged_o,_,_=_plugin_elements_all_steps(oview,caller)
        same_grid=length(merged_o)==length(elements) &&
            all(i->merged_o[i].etag==elements[i].etag,eachindex(elements))
        same_grid && (oelements=merged_o)
    end
    interpolate=!same_grid || forcerep
    interpolate && _native_diagnostic_emit_locked(:info,
        "Other view based on different grid: interpolating...")
    # `PhysicalGroup > 0` filters by the owning entity's physical tags like
    # the upstream `ge->physicals` scan; list views carry no entity records
    # (the upstream `ge == nullptr` branch), so the filter drops everything.
    entity_lookup=nothing
    if pgroup>0 && view.model_kind!==nothing
        keys=_plugin_model_keys(view)
        slot=view.model_contexts[first(keys)].model_slot
        mesh=slot.mesh;class=slot.class
        public=mesh===nothing ? nothing : _mesh_public_tags(mesh)
        entity_lookup=(mesh,class,public,slot.model)
    end
    drop_all=pgroup>0 && entity_lookup===nothing
    vars=_PluginEvalVars(zeros(9),zeros(9))
    # Union of every referenced `F<tag>` across all component expressions;
    # each tag binds the matching `v`/`w` component buffer.
    tagset=sort!(collect(Set{Int}(
        tag for prog in progs for tag in prog.field_tags)))
    bindings=Tuple(SizeField._MEFieldBinding(
        tag,_PluginEvalField(tag<16 ? tag+1 : tag-15,tag>=16))
        for tag in tagset)
    out=_PluginElement[]
    for (eindex,element) in enumerate(elements)
        va=size(element.vcoords,2);ga=size(element.geo,2);nc=element.ncomp
        # The emitted record repeats the element's own shape; only primary
        # node layouts fit our list record catalog.
        va==ga || continue
        type=_plugin_record_type(element.kind,ncomp2)
        type===nothing && continue
        drop_all && continue
        if entity_lookup!==nothing
            mesh,class,public,model=entity_lookup
            index=_plugin_element_index(mesh,public,element.etag)
            entity=index==0 ? nothing :
                _plugin_element_entity(mesh,class,index)
            phys=Int[]
            entity!==nothing && (phys=try
                model_physical_groups_for_entity(model,entity[1],entity[2])
            catch err
                err isa InterruptException && rethrow()
                Int[]
            end)
            pgroup in phys || continue
        end
        vals=Vector{Float64}(undef,va*ncomp2*length(emit_positions))
        for (ei,pos) in enumerate(emit_positions)
            stepnum=stepnums[pos]
            step2num=ots<0 ? stepnum : ots
            step2pos=findfirst(==(step2num),ostepnums)
            for node in 1:va
                fill!(vars.v,0.0);fill!(vars.w,0.0)
                for comp in 1:nc
                    vars.v[comp]=element.vals[
                        (pos-1)*va*nc+(node-1)*nc+comp]
                end
                x=element.vcoords[1,node];y=element.vcoords[2,node]
                z=element.vcoords[3,node]
                if interpolate
                    acc=_plugin_probe_at(oview,(x,y,z),step2num,caller;
                                         dimension=dimension)
                    acc===nothing || copyto!(vars.w,1,acc,1,
                                             min(9,length(acc)))
                elseif step2pos!==nothing
                    oe=oelements[eindex]
                    oa=size(oe.vcoords,2);onc=oe.ncomp
                    for comp in 1:onc
                        vars.w[comp]=oe.vals[
                            (step2pos-1)*oa*onc+(node-1)*onc+comp]
                    end
                end
                for (c,prog) in enumerate(progs)
                    vals[(ei-1)*va*ncomp2+(node-1)*ncomp2+c]=try
                        SizeField._eval_matheval_entity(
                            prog,x,y,z,bindings,vars)
                    catch err
                        err isa InterruptException && rethrow()
                        throw(ArgumentError(
                            "$caller: cannot evaluate $(exprs[c]) at " *
                            "($x,$y,$z): $err"))
                    end
                end
            end
        end
        push!(out,_PluginElement(element.kind,element.dim,element.vcoords,
                                 element.vcoords,vals,ncomp2,element.etag))
    end
    emit_times=[times[pos] for pos in emit_positions]
    _plugin_extract_emit(view.name*"_MathEval",out,-1,
                         length(emit_positions),emit_times,caller)
    return nothing
end

# ---- Scal2Vec -----------------------------------------------------------------------------

function _plugin_exec_scal2vec(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers;strings=spec.strings
    name=strings["NameNewView"]
    cviews=Tuple{Int,Any}[]
    ref=nothing
    for (comp,opt) in enumerate(("ViewX","ViewY","ViewZ"))
        raw=numbers[opt]
        if raw<0
            push!(cviews,(0,nothing))
            continue
        end
        ctag,cview=_plugin_view_tag(numbers,opt,caller)
        push!(cviews,(ctag,cview))
        ref===nothing && (ref=(ctag,cview))
    end
    ref===nothing && throw(ArgumentError(
        "$caller: Scal2Vec plugin could not find any view"))
    _,rview=ref
    elements,steps,times=_plugin_elements_all_steps(rview,caller)
    merged=Vector{Vector{_PluginElement}}(undef,3)
    for comp in 1:3
        merged[comp]=cviews[comp][2]===nothing ? _PluginElement[] :
            _plugin_elements_all_steps(cviews[comp][2],caller)[1]
    end
    out=_PluginElement[]
    for (i,element) in enumerate(elements)
        va=size(element.vcoords,2);ga=size(element.geo,2)
        va==ga || continue
        _plugin_record_type(element.kind,3)===nothing && continue
        # Unset `ViewX/Y/Z` components stay zero like the upstream `val = 0.`
        # default.
        vals=zeros(va*3*steps)
        for comp in 1:3
            cviews[comp][2]===nothing && continue
            oe=merged[comp]
            i<=length(oe) || throw(ArgumentError(
                "$caller: view tag $(cviews[comp][1]) has incompatible " *
                "element structure (expected same-grid data)"))
            oelement=oe[i];oa=size(oelement.vcoords,2);onc=oelement.ncomp
            onc>=1 && va<=oa || throw(ArgumentError(
                "$caller: view tag $(cviews[comp][1]) has incompatible " *
                "element structure (expected same-grid data)"))
            osteps=length(oelement.vals)÷(oa*onc)
            for step in 1:steps,node in 1:va
                step<=osteps || throw(ArgumentError(
                    "$caller: view tag $(cviews[comp][1]) has no time step " *
                    "$(step-1)"))
                vals[(step-1)*va*3+(node-1)*3+comp]=
                    oelement.vals[(step-1)*oa*onc+(node-1)*onc+1]
            end
        end
        push!(out,_PluginElement(element.kind,element.dim,element.vcoords,
                                 element.vcoords,vals,3,element.etag))
    end
    _plugin_extract_emit(name,out,-1,steps,times,caller)
    return nothing
end

# ---- corner component gradients ------------------------------------------------

# Spatial gradient of scalar component `comp` (1-based) at every *corner*
# node of `element`, evaluated at that node's reference location — the
# upstream `element::interpolateGrad` backing the Divergence/Curl/Lambda2
# plugins. Like upstream `getNumCornerNodes`/`getValues`, only the first
# `ga` value nodes (the MSH-order corners) are read, so high-order layouts
# differentiate their corner map. Returns a node-major `3*ga` vector, or
# `nothing` for a point element (no reference frame), a component out of
# range, missing corner value coverage, or a degenerate geometric map.
function _plugin_element_component_gradient(element::_PluginElement,step::Int,
                                            comp::Int)
    nodes=get(_PLUGIN_REFERENCE_NODES,element.kind,nothing)
    nodes===nothing && return nothing
    nc=element.ncomp
    1<=comp<=nc || return nothing
    ga=size(element.geo,2);va=size(element.vcoords,2)
    va>=ga || return nothing
    length(element.vals)>=step*va*nc || return nothing
    family=_PLUGIN_KIND_FAMILY[element.kind];dim=element.dim
    out=Vector{Float64}(undef,3ga)
    for node in 1:ga
        u,v,w=nodes[node]
        dN=MeshFunctionSpaces._first_order_gradients(
            Val(family),u,v,w,"API.plugin",node)
        jcols=ntuple(_->zeros(dim),3)
        r=zeros(dim)
        for i in 1:ga
            gi=dN[i];vi=element.vals[(step-1)*va*nc+(i-1)*nc+comp]
            for k in 1:dim
                r[k]=muladd(vi,gi[k],r[k])
                for a in 1:3
                    jcols[a][k]=muladd(element.geo[a,i],gi[k],jcols[a][k])
                end
            end
        end
        h=_plugin_metric_solve(jcols,r,dim)
        h===nothing && return nothing
        for a in 1:3
            g=0.0
            for k in 1:dim;g=muladd(jcols[a][k],h[k],g);end
            out[(node-1)*3+a]=g
        end
    end
    return out
end

# All three component gradients of a 3-component element field for one step:
# `(∇v_x, ∇v_y, ∇v_z)` as node-major 3*ga vectors, `nothing` when the
# element is not vector-valued or a corner gradient fails. Point records
# yield the upstream zero gradient (a point's `getJacobian` is the identity
# and its shape gradients vanish), matching `interpolateGrad` on `element`.
function _plugin_vector_gradients(element::_PluginElement,step::Int)
    element.ncomp==3 || return nothing
    element.kind===:point && return (zeros(3),zeros(3),zeros(3))
    gx=_plugin_element_component_gradient(element,step,1)
    gx===nothing && return nothing
    gy=_plugin_element_component_gradient(element,step,2)
    gy===nothing && return nothing
    gz=_plugin_element_component_gradient(element,step,3)
    gz===nothing && return nothing
    return (gx,gy,gz)
end

# Upstream `hasMultipleMeshes` for a model view: every populated step's
# context must live in a single mesh slot. Returns false after emitting the
# upstream-style error diagnostic when the steps span multiple meshes.
function _plugin_check_single_mesh(view,name::AbstractString)
    view.model_kind===nothing && return true
    slots=Set{Any}()
    for key in _plugin_model_keys(view)
        push!(slots,get(view.model_contexts[key],:model_slot,nothing))
    end
    length(slots)>1 || return true
    _native_diagnostic_emit_locked(:error,
        "$name plugin cannot be run on multi-mesh views")
    return false
end

# ---- Divergence ----------------------------------------------------------------

function _plugin_exec_divergence(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    tag,view=_plugin_view(numbers,caller)
    _plugin_check_single_mesh(view,"Divergence") || return nothing
    elements,steps,times=_plugin_elements_all_steps(view,caller)
    out=_PluginElement[]
    for element in elements
        # Upstream `numComp` gate: only 3-component elements are processed.
        element.ncomp==3 || continue
        _plugin_record_type(element.kind,1)===nothing && continue
        ga=size(element.geo,2)
        size(element.vcoords,2)>=ga || continue
        vals=Vector{Float64}(undef,ga*steps)
        ok=true
        for step in 1:steps
            grads=_plugin_vector_gradients(element,step)
            if grads===nothing;ok=false;break;end
            gx,gy,gz=grads
            for node in 1:ga
                vals[(step-1)*ga+node]=
                    gx[(node-1)*3+1]+gy[(node-1)*3+2]+gz[(node-1)*3+3]
            end
        end
        ok || continue
        push!(out,_PluginElement(element.kind,element.dim,element.geo,
                                 element.geo,vals,1,element.etag))
    end
    _plugin_extract_emit(view.name*"_Divergence",out,-1,steps,times,caller)
    return nothing
end

# ---- Curl ----------------------------------------------------------------------

# Shared body of upstream `Plugin(Curl)` and the Tessella `Vorticity`
# extension below: the elementwise corner-gradient curl of a 3-component
# field, emitted under `name_suffix`.
function _plugin_curl_scan(view,caller,plugin_name::AbstractString)
    _plugin_check_single_mesh(view,plugin_name) || return nothing
    elements,steps,times=_plugin_elements_all_steps(view,caller)
    out=_PluginElement[]
    for element in elements
        element.ncomp==3 || continue
        _plugin_record_type(element.kind,3)===nothing && continue
        ga=size(element.geo,2)
        size(element.vcoords,2)>=ga || continue
        vals=Vector{Float64}(undef,ga*3*steps)
        ok=true
        for step in 1:steps
            grads=_plugin_vector_gradients(element,step)
            if grads===nothing;ok=false;break;end
            gx,gy,gz=grads
            for node in 1:ga
                # curl = (∂v_z/∂y - ∂v_y/∂z, ∂v_x/∂z - ∂v_z/∂x,
                #         ∂v_y/∂x - ∂v_x/∂y)
                vals[(step-1)*ga*3+(node-1)*3+1]=
                    gz[(node-1)*3+2]-gy[(node-1)*3+3]
                vals[(step-1)*ga*3+(node-1)*3+2]=
                    gx[(node-1)*3+3]-gz[(node-1)*3+1]
                vals[(step-1)*ga*3+(node-1)*3+3]=
                    gy[(node-1)*3+1]-gx[(node-1)*3+2]
            end
        end
        ok || continue
        push!(out,_PluginElement(element.kind,element.dim,element.geo,
                                 element.geo,vals,3,element.etag))
    end
    return out,steps,times
end

function _plugin_exec_curl(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    tag,view=_plugin_view(numbers,caller)
    scanned=_plugin_curl_scan(view,caller,"Curl")
    scanned===nothing && return nothing
    out,steps,times=scanned
    _plugin_extract_emit(view.name*"_Curl",out,-1,steps,times,caller)
    return nothing
end

# ---- Lambda2 --------------------------------------------------------------------

# Eigenvalues of a 3×3 matrix through its characteristic cubic, mirroring
# upstream `eigenvalue`: `FindCubicRoots` on `v³ - I₁ v² + I₂ v - I₃`
# (I₂ = ½(tr² - tr(T²))) followed by the descending `eigsort` selection.
# Only real parts are kept, so a nearly-symmetric matrix with a numerically
# complex pair degrades exactly like upstream.
function _plugin_eigenvalues3(g::Matrix{Float64})
    tr=g[1,1]+g[2,2]+g[3,3]
    tr2=0.0
    for i in 1:3,j in 1:3
        tr2=muladd(g[i,j],g[j,i],tr2)
    end
    det=g[1,1]*(g[2,2]*g[3,3]-g[2,3]*g[3,2])-
        g[1,2]*(g[2,1]*g[3,3]-g[2,3]*g[3,1])+
        g[1,3]*(g[2,1]*g[3,2]-g[2,2]*g[3,1])
    # The monic cubic's coefficients after dividing by a = 1.
    b=-tr;c=0.5*(tr*tr-tr2);d=-det
    q=(3c-b*b)/9
    r=(-27d+b*(9c-2*b*b))/54
    discrim=q*q*q+r*r
    term1=b/3
    d1,d2,d3=0.0,0.0,0.0
    if discrim>0
        # One real root; the complex pair contributes its shared real part
        # like upstream's real[] slots.
        s=cbrt(r+sqrt(discrim));t=cbrt(r-sqrt(discrim))
        d1=-term1+s+t;d2=d3=-(term1+(s+t)/2)
    elseif discrim==0.0
        r13=cbrt(r)
        d1=-term1+2r13;d2=d3=-(r13+term1)
    else
        q=-q
        dum1=acos(clamp(r/sqrt(q*q*q),-1.0,1.0))
        r13=2*sqrt(q)
        d1=-term1+r13*cos(dum1/3)
        d2=-term1+r13*cos((dum1+2.0*pi)/3)
        d3=-term1+r13*cos((dum1+4.0*pi)/3)
    end
    roots=[d1,d2,d3]
    for i in 1:3
        k=i;p=roots[i]
        for j in i+1:3
            if roots[j]>=p;p=roots[j];k=j;end
        end
        if k!=i;roots[k]=roots[i];roots[i]=p;end
    end
    return (roots[1],roots[2],roots[3])
end

function _plugin_exec_lambda2(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    ev=_mesh_query_integer(Int(numbers["Eigenvalue"]),caller,"Eigenvalue")
    if !(1<=ev<=3)
        _native_diagnostic_emit_locked(:error,
            "Eigenvalue should be 1, 2 or 3")
        return nothing
    end
    tag,view=_plugin_view(numbers,caller)
    # Upstream has no explicit multi-mesh check, but our per-tag step merge
    # cannot identify same-numbered elements across meshes.
    _plugin_check_single_mesh(view,"Lambda2") || return nothing
    elements,steps,times=_plugin_elements_all_steps(view,caller)
    out=_PluginElement[]
    for element in elements
        nc=element.ncomp
        nc in (3,9) || continue
        _plugin_record_type(element.kind,1)===nothing && continue
        ga=size(element.geo,2);va=size(element.vcoords,2)
        va>=ga || continue
        vals=Vector{Float64}(undef,ga*steps)
        ok=true
        for step in 1:steps
            grads=nc==3 ? _plugin_vector_gradients(element,step) : nothing
            nc==3 && grads===nothing && (ok=false;break)
            for node in 1:ga
                # The velocity gradient g[i][j] = ∂v_i/∂x_j, read straight
                # from the stored tensor for 9-component records.
                g=Matrix{Float64}(undef,3,3)
                if nc==9
                    for i in 1:3,j in 1:3
                        g[i,j]=element.vals[
                            (step-1)*va*nc+(node-1)*nc+3(i-1)+j]
                    end
                else
                    for i in 1:3,j in 1:3
                        g[i,j]=grads[i][(node-1)*3+j]
                    end
                end
                # A[m][n] = Σ_l sym[m][l]·sym[l][n] + asym[m][l]·asym[l][n]
                # over the symmetric/antisymmetric halves of g.
                a=zeros(3,3)
                for m in 1:3,n in 1:3
                    acc=0.0
                    for l in 1:3
                        sml=0.5*(g[m,l]+g[l,m]);sln=0.5*(g[l,n]+g[n,l])
                        oml=0.5*(g[m,l]-g[l,m]);oln=0.5*(g[l,n]-g[n,l])
                        acc=muladd(sml,sln,muladd(oml,oln,acc))
                    end
                    a[m,n]=acc
                end
                vals[(step-1)*ga+node]=_plugin_eigenvalues3(a)[ev]
            end
        end
        ok || continue
        push!(out,_PluginElement(element.kind,element.dim,element.geo,
                                 element.geo,vals,1,element.etag))
    end
    _plugin_extract_emit(view.name*"_Lambda2",out,-1,steps,times,caller)
    return nothing
end

# ---- SphericalRaise --------------------------------------------------------------

function _plugin_exec_spherical_raise(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    center=(_mesh_query_coordinate(numbers["Xc"],caller,"Xc"),
            _mesh_query_coordinate(numbers["Yc"],caller,"Yc"),
            _mesh_query_coordinate(numbers["Zc"],caller,"Zc"))
    raise=_mesh_query_coordinate(numbers["Raise"],caller,"Raise")
    offset=_mesh_query_coordinate(numbers["Offset"],caller,"Offset")
    timestep=_mesh_query_integer(Int(numbers["TimeStep"]),caller,"TimeStep")
    tag,view=_plugin_view(numbers,caller)
    view.model_kind===nothing || throw(ArgumentError(
        "$caller: SphericalRaise rewrites list data only; model views are unsupported"))
    steps=_plugin_list_steps(view,caller)
    if timestep<0 || timestep>steps-1
        _native_diagnostic_emit_locked(:error,
            "Invalid TimeStep ($timestep) in view")
        return nothing
    end
    for type in view.list_types
        layout=get(Post._POS_RECORD_TYPES,type,nothing)
        layout===nothing && throw(ArgumentError(
            "$caller: unsupported list data type $type"))
        _,arity,ncomp=layout
        count=view.list_counts[type];data=view.list_data[type]
        width=3arity+Base.checked_mul(arity,Base.checked_mul(ncomp,steps))
        length(data)==Base.checked_mul(count,width) || throw(ArgumentError(
            "$caller: list data payload does not align with numEle"))
        for element in 0:count-1
            base=element*width
            for node in 1:arity
                x=data[base+node];y=data[base+arity+node]
                z=data[base+2arity+node]
                rx=x-center[1];ry=y-center[2];rz=z-center[3]
                s=hypot(rx,ry,rz)
                # Upstream `norme` would turn a centred node into NaNs;
                # leaving it in place is the useful behaviour.
                s>0 || continue
                val=_plugin_scalar_rep(data,
                    base+3arity+timestep*arity*ncomp+(node-1)*ncomp,ncomp)
                coef=(offset+raise*val)/s
                data[base+node]+=coef*rx
                data[base+arity+node]+=coef*ry
                data[base+2arity+node]+=coef*rz
            end
        end
    end
    return nothing
end

# ---- Annotate ----------------------------------------------------------------------

function _plugin_exec_annotate(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers;strings=spec.strings
    x=_mesh_query_coordinate(numbers["X"],caller,"X")
    y=_mesh_query_coordinate(numbers["Y"],caller,"Y")
    z=_mesh_query_coordinate(numbers["Z"],caller,"Z")
    dim3=!iszero(numbers["ThreeD"])
    tag,view=_plugin_view(numbers,caller)
    text=strings["Text"]
    # Session string records hold one entry per element: the upstream 2-D
    # (`T2`, screen coordinates) or 3-D (`T3`, model coordinates) anchor
    # plus the text. There is no per-string style slot, so the FontSize/
    # Font/Align options are accepted but not retained. Stored coords are
    # `width` per record (3 = xyz, 4 = xyz+style); new records pad with the
    # table's existing width so mixed appends keep the flat layout aligned.
    key=dim3 ? "T3" : "T2"
    coord=dim3 ? Float64[x,y,z] : Float64[x,y,0.0]
    if view.model_kind===nothing
        if !(key in view.list_string_types)
            push!(view.list_string_types,key)
            view.list_strings[key]=String[]
            view.list_string_coords[key]=Float64[]
            view.list_string_counts[key]=0
        end
        prior=view.list_string_coords[key]
        width=view.list_string_counts[key]>0 ?
            div(length(prior),view.list_string_counts[key]) : 4
        push!(view.list_strings[key],text)
        append!(view.list_string_coords[key],
                width==4 ? vcat(coord,0.0) : coord)
        view.list_string_counts[key]+=1
        return nothing
    end
    # Model-based input: upstream creates a new list-based view carrying the
    # string record plus a copy of the source time table.
    out=_SessionView(view.name*"_Annotate")
    push!(out.list_string_types,key)
    out.list_strings[key]=String[text]
    out.list_string_coords[key]=vcat(coord,0.0)
    out.list_string_counts[key]=1
    # Upstream copies the whole `getNumTimeSteps` table — the step-range
    # extent, including unpopulated gaps.
    keys=_plugin_model_keys(view)
    steps=max(_view_data_steps(view),isempty(keys) ? 0 : last(keys)[1]+1)
    out.list_time=[_plugin_step_time(view,step) for step in 0:steps-1]
    out.numbers["NbTimeStep"]=Float64(steps)
    newtag=_fresh_view_tag()
    SESSION_VIEWS[newtag]=out
    VIEW_TAG_MAX[]=max(VIEW_TAG_MAX[],newtag)
    return nothing
end

# ---- Eigenvalues ------------------------------------------------------------------

# Upstream `Plugin(Eigenvalues)`: each 9-component element's tensor (the nine
# components filling `A[3][3]` row-major) is reduced to its three real
# eigenvalues through the shared characteristic-cubic solver
# (`_plugin_eigenvalues3` — the same FindCubicRoots + descending eigsort as
# `Numeric.cpp::eigenvalue`, so `w[1]` is the largest and `w[3]` the smallest
# root). Every populated step re-emits on the element's linear corner set
# into the three scalar list views `*_MinEigenvalues`, `*_MidEigenvalues` and
# `*_MaxEigenvalues`; the three views are created even when no element
# qualifies, and elements without 9 components are skipped like upstream's
# `numComp != 9` gate.
function _plugin_exec_eigenvalues(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    tag,view=_plugin_view(numbers,caller)
    _plugin_check_single_mesh(view,"Eigenvalues") || return nothing
    elements,steps,times=_plugin_elements_all_steps(view,caller)
    out=(_PluginElement[],_PluginElement[],_PluginElement[])
    for element in elements
        element.ncomp==9 || continue
        _plugin_record_type(element.kind,1)===nothing && continue
        ga=size(element.geo,2);va=size(element.vcoords,2)
        va>=ga || continue
        vals=(Vector{Float64}(undef,ga*steps),Vector{Float64}(undef,ga*steps),
              Vector{Float64}(undef,ga*steps))
        for step in 1:steps,node in 1:ga
            a=Matrix{Float64}(undef,3,3)
            base=(step-1)*va*9+(node-1)*9
            for i in 1:3,j in 1:3
                a[i,j]=element.vals[base+3(i-1)+j]
            end
            w=_plugin_eigenvalues3(a)
            slot=(step-1)*ga+node
            vals[1][slot]=w[3];vals[2][slot]=w[2];vals[3][slot]=w[1]
        end
        for k in 1:3
            push!(out[k],_PluginElement(element.kind,element.dim,element.geo,
                                        element.geo,vals[k],1,element.etag))
        end
    end
    _plugin_extract_emit(view.name*"_MinEigenvalues",out[1],-1,steps,times,caller)
    _plugin_extract_emit(view.name*"_MidEigenvalues",out[2],-1,steps,times,caller)
    _plugin_extract_emit(view.name*"_MaxEigenvalues",out[3],-1,steps,times,caller)
    return nothing
end

# ---- ModulusPhase -----------------------------------------------------------------

# Row offsets of one populated model-data block — the same `counts`-aware
# split `_plugin_model_elements` performs for ElementNodeData's variable row
# widths.
function _plugin_model_row_offsets(view,key,tags,values,caller::AbstractString)
    ctx=view.model_contexts[key]
    counts=get(ctx,:counts,Int32[])
    offsets=Vector{Int}(undef,length(tags)+1);offsets[1]=1
    for row in eachindex(tags)
        width=row<=length(counts) ? Int(counts[row]) :
            length(values)÷length(tags)
        offsets[row+1]=Base.checked_add(offsets[row],width)
    end
    offsets[end]-1==length(values) || throw(ArgumentError(
        "$caller: malformed model data row extents"))
    return offsets
end

# Upstream `Plugin(ModulusPhase)` runs in place: steps `RealPart` and
# `ImaginaryPart` of `View` hold a complex field; both steps are rewritten to
# the per-node modulus `sqrt(re^2+im^2)` and phase `atan2(im,re)` of every
# component. List payloads mutate step slots directly; model payloads rewrite
# the two steps' stored rows, pairing each `RealPart` row with the
# `ImaginaryPart` row of the same tag — upstream reaches each node once
# through its element scan plus NodeData tagging, which per-tag row storage
# already guarantees here. Out-of-range indices emit the upstream diagnostic
# and leave the view unchanged.
function _plugin_exec_modulus_phase(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    tag,view=_plugin_view(numbers,caller)
    _plugin_check_single_mesh(view,"ModulusPhase") || return nothing
    rindex=_mesh_query_integer(Int(numbers["RealPart"]),caller,"RealPart")
    iindex=_mesh_query_integer(Int(numbers["ImaginaryPart"]),caller,"ImaginaryPart")
    # Upstream bounds-checks the indices against `getNumTimeSteps`.
    steps=_view_data_steps(view)
    if !(0<=rindex<steps) || !(0<=iindex<steps)
        _native_diagnostic_emit_locked(:error,
            "Wrong real or imaginary part index")
        return nothing
    end
    if view.model_kind===nothing
        for type in view.list_types
            layout=get(Post._POS_RECORD_TYPES,type,nothing)
            layout===nothing && throw(ArgumentError(
                "$caller: unsupported list data type $type"))
            _,arity,ncomp=layout
            count=view.list_counts[type]
            count==0 && continue
            data=view.list_data[type]
            # Per-type step count (payload width); the global check above
            # already guarantees rindex/iindex fit every populated type.
            tsteps=(length(data)÷count-3arity)÷(arity*ncomp)
            (rindex<tsteps && iindex<tsteps) || continue
            width=3arity+arity*ncomp*tsteps
            for element in 0:count-1
                base=element*width+3arity
                for node in 1:arity,comp in 1:ncomp
                    ri=base+rindex*arity*ncomp+(node-1)*ncomp+comp
                    ii=base+iindex*arity*ncomp+(node-1)*ncomp+comp
                    vr=data[ri];vi=data[ii]
                    data[ri]=sqrt(vr*vr+vi*vi)
                    data[ii]=atan(vi,vr)
                end
            end
        end
    else
        kind=view.model_kind
        rblock=get(view.model_data,(rindex,kind),nothing)
        iblock=get(view.model_data,(iindex,kind),nothing)
        if rblock!==nothing && iblock!==nothing && rblock[3]==iblock[3]
            rtags,rvalues,_=rblock;itags,ivalues,_=iblock
            roff=_plugin_model_row_offsets(view,(rindex,kind),rtags,rvalues,caller)
            ioff=_plugin_model_row_offsets(view,(iindex,kind),itags,ivalues,caller)
            irow=Dict{UInt64,Int}(t=>i for (i,t) in enumerate(itags))
            for (r,t) in enumerate(rtags)
                i=get(irow,t,0);i==0 && continue
                width=roff[r+1]-roff[r]
                width==ioff[i+1]-ioff[i] || continue
                for k in 0:width-1
                    vr=rvalues[roff[r]+k];vi=ivalues[ioff[i]+k]
                    rvalues[roff[r]+k]=sqrt(vr*vr+vi*vi)
                    ivalues[ioff[i]+k]=atan(vi,vr)
                end
            end
        end
    end
    view.name=view.name*"_ModulusPhase"
    return nothing
end

# ---- HarmonicToTime -----------------------------------------------------------------

# Upstream `Plugin(HarmonicToTime)`: the `RealPart`/`ImaginaryPart` steps of
# `View` define a harmonic field; the new `*_HarmonicToTime` list view holds
# `NumSteps` samples `re*cos(p) + tsign*im*sin(p)` per node and component,
# with `p = 2*pi*NumPeriods*k/NumSteps` (`p = 0` when `Frequency` is 0) and
# physical sample times `2*pi*NumPeriods*k/(Frequency*NumSteps)` (0 for zero
# frequency). `tsign` is +1 for `TimeSign` > 0 and -1 otherwise.
function _plugin_exec_harmonic_to_time(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    tag,view=_plugin_view(numbers,caller)
    _plugin_check_single_mesh(view,"HarmonicToTime") || return nothing
    rindex=_mesh_query_integer(Int(numbers["RealPart"]),caller,"RealPart")
    iindex=_mesh_query_integer(Int(numbers["ImaginaryPart"]),caller,"ImaginaryPart")
    nsteps=_mesh_query_integer(Int(numbers["NumSteps"]),caller,"NumSteps")
    tsign=numbers["TimeSign"]>0 ? 1.0 : -1.0
    frequency=numbers["Frequency"]
    nperiods=_mesh_query_integer(Int(numbers["NumPeriods"]),caller,"NumPeriods")
    elements,steps,_=_plugin_elements_all_steps(view,caller)
    if !(0<=rindex<steps) || !(0<=iindex<steps)
        _native_diagnostic_emit_locked(:error,"Wrong real or imaginary part index")
        return nothing
    end
    if nsteps<=0
        _native_diagnostic_emit_locked(:error,"nSteps should be > 0")
        return nothing
    end
    out=_PluginElement[]
    for element in elements
        nc=element.ncomp
        _plugin_record_type(element.kind,nc)===nothing && continue
        ga=size(element.geo,2);va=size(element.vcoords,2)
        va>=ga || continue
        vals=Vector{Float64}(undef,ga*nc*nsteps)
        for k in 0:nsteps-1
            p=iszero(frequency) ? 0.0 : 2.0*pi*nperiods*k/nsteps
            cp=cos(p);sp=tsign*sin(p)
            for node in 1:ga,comp in 1:nc
                slot=k*ga*nc+(node-1)*nc+comp
                base=(node-1)*nc+comp
                vals[slot]=muladd(element.vals[rindex*va*nc+base],cp,
                                  element.vals[iindex*va*nc+base]*sp)
            end
        end
        push!(out,_PluginElement(element.kind,element.dim,element.geo,
                                 element.geo,vals,nc,element.etag))
    end
    times=Float64[iszero(frequency) ? 0.0 :
                  2.0*pi*nperiods*k/(frequency*nsteps) for k in 0:nsteps-1]
    _plugin_extract_emit(view.name*"_HarmonicToTime",out,-1,nsteps,times,caller)
    return nothing
end

# ---- Density -----------------------------------------------------------------------

# `Plugin(Density)` was dropped upstream long before Gmsh 4.15.2; the
# documented transform — the per-node mesh density, the number of elements
# incident on each node — is kept here as a single-step scalar `*_Density`
# list view. Incident elements are counted once each (the multi-step merge
# already deduplicates element tags) and nodes are identified by exact
# coordinates, the list-record convention `Smooth` also uses. `View` is the
# only surviving upstream option.
function _plugin_exec_density(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    tag,view=_plugin_view(numbers,caller)
    _plugin_check_single_mesh(view,"Density") || return nothing
    elements,_,_=_plugin_elements_all_steps(view,caller)
    incidence=Dict{NTuple{3,Float64},Int}()
    for element in elements
        for node in axes(element.geo,2)
            key=(element.geo[1,node],element.geo[2,node],element.geo[3,node])
            incidence[key]=get(incidence,key,0)+1
        end
    end
    out=_PluginElement[]
    for element in elements
        _plugin_record_type(element.kind,1)===nothing && continue
        ga=size(element.geo,2)
        vals=Vector{Float64}(undef,ga)
        for node in 1:ga
            vals[node]=Float64(incidence[(element.geo[1,node],
                element.geo[2,node],element.geo[3,node])])
        end
        push!(out,_PluginElement(element.kind,element.dim,element.geo,
                                 element.geo,vals,1,element.etag))
    end
    _plugin_extract_emit(view.name*"_Density",out,-1,1,Float64[0.0],caller)
    return nothing
end

# ---- MakeSimplex ----------------------------------------------------------------------

# Source element kinds → the simplex kind they decompose into; the node maps
# themselves live in `_PLUGIN_ISO_SIMPLICES` (the same 1-based tables as
# upstream `MakeSimplex::decompose`: quadrangle → 2 triangles, hexahedron → 6
# tetrahedra, prism → 3 tetrahedra, pyramid → 2 tetrahedra).
const _PLUGIN_SIMPLEX_TARGET=(quadrangle=:triangle,hexahedron=:tetrahedron,
                              prism=:tetrahedron,pyramid=:tetrahedron)

# Upstream `Plugin(MakeSimplex)` rewrites a list view in place (`decomposeList`
# clears the source record vectors and appends the simplex records onto the
# target lists, after any records that already live there). Scalar, vector
# and tensor records decompose identically and every populated step is
# carried through under the reordered node indices; point/line/triangle/
# tetrahedron records pass through. Model views are rejected like upstream's
# `getDataList` null check.
function _plugin_exec_makesimplex(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    tag,view=_plugin_view(numbers,caller)
    view.model_kind===nothing || throw(ArgumentError(
        "$caller: MakeSimplex rewrites list data only; model views are unsupported"))
    steps=_plugin_list_steps(view,caller)
    types=String[];counts=Dict{String,Int}();datas=Dict{String,Vector{Float64}}()
    for type in view.list_types
        layout=get(Post._POS_RECORD_TYPES,type,nothing)
        layout===nothing && throw(ArgumentError(
            "$caller: unsupported list data type $type"))
        kind,arity,ncomp=layout
        count=view.list_counts[type];data=view.list_data[type]
        target=get(_PLUGIN_SIMPLEX_TARGET,kind,nothing)
        if target===nothing
            if haskey(datas,type)
                # An earlier decomposed block already staged this list:
                # upstream appends the new records after the pre-existing
                # ones, so the stored payload precedes the staged output.
                counts[type]+=count
                datas[type]=vcat(data,datas[type])
            else
                push!(types,type);counts[type]=count;datas[type]=data
            end
            continue
        end
        ttype=_PLUGIN_COMP_HEAD[ncomp]*_PLUGIN_KIND_SHAPE[target]
        nsub=length(first(_PLUGIN_ISO_SIMPLICES[kind]))
        swidth=3arity+Base.checked_mul(arity,Base.checked_mul(ncomp,steps))
        length(data)==Base.checked_mul(count,swidth) || throw(ArgumentError(
            "$caller: list data payload does not align with numEle"))
        if !haskey(datas,ttype)
            push!(types,ttype);datas[ttype]=Float64[];counts[ttype]=0
        end
        out=datas[ttype]
        for element in 0:count-1
            base=element*swidth
            for map in _PLUGIN_ISO_SIMPLICES[kind]
                for axis in 1:3,node in 1:nsub
                    push!(out,data[base+(axis-1)*arity+map[node]])
                end
                for step in 1:steps,node in 1:nsub,comp in 1:ncomp
                    push!(out,data[base+3arity+(step-1)*arity*ncomp+
                                   (map[node]-1)*ncomp+comp])
                end
                counts[ttype]+=1
            end
        end
    end
    view.list_types=types;view.list_counts=counts;view.list_data=datas
    view.numbers["NbTimeStep"]=Float64(steps)
    return nothing
end

# ---- ModifyComponents -----------------------------------------------------------------

# Evaluate the nine ModifyComponents expressions at one node; upstream
# `f.eval` runs the whole vector in one call and drops every write on
# failure, so a single error cancels the node's update entirely.
function _plugin_mc_eval(progs,bindings,vars,x,y,z)
    res=Vector{Float64}(undef,9)
    for i in 1:9
        res[i]=try
            SizeField._eval_matheval_entity(progs[i],x,y,z,bindings,vars)
        catch err
            err isa InterruptException && rethrow()
            return nothing
        end
    end
    return res
end

# Upstream `Plugin(ModifyComponents)` rewrites `View` in place: `Expression<i>`
# replaces component i at every visited node while empty expressions leave
# that component untouched (they evaluate as `0.`). Variables x, y, z, `Time`
# (the step's physical time), `TimeStep` (its index), v0..v8 (the view's own
# components at `TimeStep`) and w0..w8 (`OtherView` at `OtherTimeStep`, or
# the same step when `OtherTimeStep < 0`) feed the shared math evaluator.
# `w` reads the same (element, node) slots when both views share their grid,
# otherwise `_plugin_probe_at` interpolates like the upstream OctreePost
# search (`ForceInterpolation` forces that path). Model `NodeData` rows are
# visited once per node — the upstream `tagNode` dedup — while ElementData
# rows compose per scanned node like upstream's per-node read/write slot.
# Invalid `TimeStep`/`OtherTimeStep`/unknown `OtherView` emit the upstream
# diagnostics and fall back to step 0 / self instead of raising.
function _plugin_exec_modifycomponents(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers;strings=spec.strings
    tag,view=_plugin_view(numbers,caller)
    raw_exprs=[strings["Expression$i"] for i in 0:8]
    progs=[SizeField.parse_matheval(_plugin_matheval_rewrite(
               isempty(e) ? "0." : e,caller;time_vars=true))
           for e in raw_exprs]
    active=[!isempty(e) for e in raw_exprs]
    timestep=_mesh_query_integer(Int(numbers["TimeStep"]),caller,"TimeStep")
    ots=_mesh_query_integer(Int(numbers["OtherTimeStep"]),caller,"OtherTimeStep")
    forcerep=!iszero(numbers["ForceInterpolation"])
    nsteps1=_view_data_steps(view)
    if timestep>nsteps1-1
        _native_diagnostic_emit_locked(:error,
            "Invalid time step ($timestep) in View[$tag]: using step 0 instead")
        timestep=0
    end
    raw_other=numbers["OtherView"]
    otag,oview=tag,view
    if raw_other>=0
        raw_other==floor(raw_other) || throw(ArgumentError(
            "$caller: OtherView tag must be an integer"))
        oi=Int(raw_other)
        if haskey(SESSION_VIEWS,oi)
            otag,oview=oi,_session_view_data(oi,caller)
        else
            _native_diagnostic_emit_locked(:error,
                "View[$oi] does not exist: using self")
        end
    end
    nsteps2=oview===view ? nsteps1 : _view_data_steps(oview)
    if ots<0 && nsteps2!=nsteps1
        _native_diagnostic_emit_locked(:error,
            "Number of time steps don't match: using step 0")
        ots=0
    elseif ots>nsteps2-1
        _native_diagnostic_emit_locked(:error,
            "Invalid time step ($ots) in View[$otag]: using step 0 instead")
        ots=0
    end
    # Upstream builds the octree when `ForceInterpolation` is set or the
    # entity/element counts differ; identical storage reads `OtherView` at the
    # same (entity, element, node) slot — for tag-keyed model blocks that
    # means an equal tag vector in the same mesh slot, and for list views the
    # same record layout.
    same_grid=false
    if oview===view
        same_grid=true
    elseif view.model_kind===nothing && oview.model_kind===nothing
        same_grid=view.list_types==oview.list_types &&
                  view.list_counts==oview.list_counts
    elseif view.model_kind!==nothing && oview.model_kind!==nothing &&
           view.model_kind==oview.model_kind
        same_grid=true
        for key in _plugin_model_keys(view)
            (timestep<0 || key[1]==timestep) || continue
            step2=ots<0 ? key[1] : ots
            ob=get(oview.model_data,(step2,view.model_kind),nothing)
            octx=get(oview.model_contexts,(step2,view.model_kind),nothing)
            if ob===nothing || ob[1]!=view.model_data[key][1] ||
               octx===nothing ||
               octx.model_slot!==view.model_contexts[key].model_slot
                same_grid=false;break
            end
        end
    end
    interpolate=!same_grid || forcerep
    interpolate && _native_diagnostic_emit_locked(:info,
        "Other view based on different grid: interpolating...")
    vars=_PluginEvalVars(zeros(9),zeros(9))
    timefield=_PluginEvalScalar(0.0);stepfield=_PluginEvalScalar(0.0)
    tagset=sort!(collect(Set{Int}(
        tag for prog in progs for tag in prog.field_tags)))
    bindings=Tuple(
        tag==30 ? SizeField._MEFieldBinding(tag,timefield) :
        tag==31 ? SizeField._MEFieldBinding(tag,stepfield) :
        SizeField._MEFieldBinding(tag,_PluginEvalField(
            tag<16 ? tag+1 : tag-15,tag>=16))
        for tag in tagset)
    if view.model_kind===nothing
        steps=_plugin_list_steps(view,caller)
        osteps=oview===view ? steps :
            (oview.model_kind===nothing ? _plugin_list_steps(oview,caller) : 0)
        for stepnum in 0:steps-1
            (timestep<0 || stepnum==timestep) || continue
            step2=ots<0 ? stepnum : ots
            timefield.value=_plugin_step_time(view,stepnum)
            stepfield.value=Float64(stepnum)
            for type in view.list_types
                layout=get(Post._POS_RECORD_TYPES,type,nothing)
                layout===nothing && throw(ArgumentError(
                    "$caller: unsupported list data type $type"))
                _,arity,nc=layout
                count=view.list_counts[type];data=view.list_data[type]
                width=3arity+arity*nc*steps
                length(data)==count*width || throw(ArgumentError(
                    "$caller: list data payload does not align with numEle"))
                odata=same_grid && oview.model_kind===nothing ?
                    get(oview.list_data,type,nothing) : nothing
                owidth=3arity+arity*nc*osteps
                if odata!==nothing && length(odata)!=count*owidth
                    odata=nothing
                end
                for element in 0:count-1
                    base=element*width
                    for node in 1:arity
                        x=data[base+node];y=data[base+arity+node]
                        z=data[base+2arity+node]
                        vbase=base+3arity+stepnum*arity*nc
                        fill!(vars.v,0.0);fill!(vars.w,0.0)
                        for comp in 1:nc
                            vars.v[comp]=data[vbase+(node-1)*nc+comp]
                        end
                        if interpolate
                            acc=_plugin_probe_at(oview,(x,y,z),step2,caller)
                            acc===nothing || copyto!(vars.w,1,acc,1,
                                                     min(9,length(acc)))
                        elseif odata!==nothing && 0<=step2<osteps
                            for comp in 1:min(9,nc)
                                vars.w[comp]=odata[
                                    element*owidth+3arity+step2*arity*nc+
                                    (node-1)*nc+comp]
                            end
                        end
                        res=_plugin_mc_eval(progs,bindings,vars,x,y,z)
                        res===nothing && continue
                        for comp in 1:nc
                            active[comp] &&
                                (data[vbase+(node-1)*nc+comp]=res[comp])
                        end
                    end
                end
            end
        end
        return nothing
    end
    kind=view.model_kind
    for key in _plugin_model_keys(view)
        stepnum=key[1]
        (timestep<0 || stepnum==timestep) || continue
        step2=ots<0 ? stepnum : ots
        timefield.value=_plugin_step_time(view,stepnum)
        stepfield.value=Float64(stepnum)
        ctx=view.model_contexts[key]
        entity_tags,values,nc=view.model_data[key]
        offsets=_plugin_model_row_offsets(view,key,entity_tags,values,caller)
        mesh=_plugin_model_mesh(ctx.model_slot,caller)
        public=_mesh_public_tags(mesh)
        # Aligned `OtherView` block for this step, resolved by the shared tag
        # vector (same storage ⇒ same row positions).
        oblock=nothing;ooffsets=Int[];onc=0
        if same_grid
            ob=get(oview.model_data,(step2,kind),nothing)
            octx=get(oview.model_contexts,(step2,kind),nothing)
            if ob!==nothing && ob[1]==entity_tags && octx!==nothing
                ooffsets=_plugin_model_row_offsets(oview,(step2,kind),
                                                 ob[1],ob[2],caller)
                oblock=ob;onc=ob[3]
            end
        end
        if kind=="NodeData"
            row_of=Dict{UInt64,Int}(tag=>row
                                    for (row,tag) in enumerate(entity_tags))
            seen=falses(length(entity_tags))
            _,_,total=_mesh_element_offsets(mesh)
            for index in 1:total
                rec=mesh_element_record(mesh,index)
                for nod in eachindex(rec.node_tags)
                    idx=Int(rec.node_tags[nod])
                    ntag=_cache_node_tag(mesh,idx)
                    row=get(row_of,ntag,0)
                    (row==0 || seen[row]) && continue
                    seen[row]=true
                    base=offsets[row]
                    fill!(vars.v,0.0);fill!(vars.w,0.0)
                    for comp in 1:nc
                        vars.v[comp]=values[base+comp-1]
                    end
                    x=mesh.coords[1,idx];y=mesh.coords[2,idx]
                    z=mesh.coords[3,idx]
                    if interpolate
                        acc=_plugin_probe_at(oview,(x,y,z),step2,caller)
                        acc===nothing || copyto!(vars.w,1,acc,1,
                                                 min(9,length(acc)))
                    elseif oblock!==nothing
                        for comp in 1:min(9,onc)
                            vars.w[comp]=oblock[2][ooffsets[row]+comp-1]
                        end
                    end
                    res=_plugin_mc_eval(progs,bindings,vars,x,y,z)
                    res===nothing && continue
                    for comp in 1:nc
                        active[comp] && (values[base+comp-1]=res[comp])
                    end
                end
            end
        else
            for (row,etag) in enumerate(entity_tags)
                index=_plugin_element_index(mesh,public,etag)
                index==0 && continue
                rec=mesh_element_record(mesh,index)
                base=offsets[row];rwidth=offsets[row+1]-base
                for nod in eachindex(rec.node_tags)
                    idx=Int(rec.node_tags[nod])
                    x=mesh.coords[1,idx];y=mesh.coords[2,idx]
                    z=mesh.coords[3,idx]
                    vslot=kind=="ElementData" ? base+0 :
                        base+(nod-1)*nc
                    fill!(vars.v,0.0);fill!(vars.w,0.0)
                    for comp in 1:nc
                        pos=vslot+comp-1
                        pos<base+rwidth || break
                        vars.v[comp]=values[pos]
                    end
                    if interpolate
                        acc=_plugin_probe_at(oview,(x,y,z),step2,caller)
                        acc===nothing || copyto!(vars.w,1,acc,1,
                                                 min(9,length(acc)))
                    elseif oblock!==nothing
                        obase=ooffsets[row];owidth=ooffsets[row+1]-obase
                        for comp in 1:min(9,onc)
                            pos=kind=="ElementData" ? obase+comp-1 :
                                obase+(nod-1)*onc+comp-1
                            pos<obase+owidth || break
                            vars.w[comp]=oblock[2][pos]
                        end
                    end
                    res=_plugin_mc_eval(progs,bindings,vars,x,y,z)
                    res===nothing && continue
                    for comp in 1:nc
                        active[comp] || continue
                        pos=vslot+comp-1
                        pos<base+rwidth || break
                        values[pos]=res[comp]
                    end
                end
            end
        end
    end
    return nothing
end

# ---- FieldFromAmplitudePhase ---------------------------------------------------------

# Scalar-only octree probe: upstream `searchScalar` zero-fills the slot on a
# miss and only reads scalar (1-component) element lists — non-scalar probes
# therefore contribute 0 here as well.
function _plugin_scalar_probe(view,p,stepnum::Int,caller::AbstractString)
    acc=_plugin_probe_at(view,p,stepnum,caller)
    (acc===nothing || length(acc)!=1) && return 0.0
    return acc[1]
end

# Upstream `Plugin(FieldFromAmplitudePhase)` builds the complex field
# `u = a * exp(j*k*phi)` from the real/imaginary amplitude steps of
# `AmplitudeView` (exactly two time steps) and the scalar phase of
# `PhaseView`, interpolated onto the nodes of the mesh read from `MeshFile`
# (with `MeshFile == ""` falling back to `<current model name>.msh`). The
# result is a two-step `NodeData` view named `FieldFromAPhi` (real then
# imaginary) bound to a fresh model slot that owns the imported mesh; the
# current model is hidden first like upstream `setVisibility(0)` and stays
# current (upstream `GModel::current()` keeps its index). Source view options
# resolve as session view tags like `View`/`OtherView`.
function _plugin_exec_field_from_amplitude_phase(spec::_PluginSpec,
                                                 caller::AbstractString)
    numbers=spec.numbers;strings=spec.strings
    k=numbers["Wavenumber"]
    meshfile=strings["MeshFile"]
    _slot_capture_locked!()
    cindex=_current_slot_index_locked()
    cindex==0 && throw(ArgumentError("$caller: no current model"))
    cslot=MODEL_SLOTS[cindex]
    name_model=""
    if isempty(meshfile)
        _native_diagnostic_emit_locked(:info,
            "Could not find mesh file for interpolating U=A*exp(j*k*phi). " *
            "Trying to use current model mesh, instead.")
        name_model=cslot.name
        meshfile=name_model*".msh"
    end
    atag,aview=_plugin_view_tag(numbers,"AmplitudeView",caller)
    if _view_data_steps(aview)!=2
        _native_diagnostic_emit_locked(:error,
            "Invalid number of time steps for AView, it must be complex!")
        return nothing
    end
    ptag,pview=_plugin_view_tag(numbers,"PhaseView",caller)
    _plugin_check_single_mesh(aview,"FieldFromAmplitudePhase") || return nothing
    _plugin_check_single_mesh(pview,"FieldFromAmplitudePhase") || return nothing
    set_entity_visibility!(cslot.model,model_entities(cslot.model),0)
    mixed=try
        Elements.read_mixed_msh(meshfile)
    catch err
        err isa InterruptException && rethrow()
        _native_diagnostic_emit_locked(:error,
            "Unable to open file '$meshfile'")
        return nothing
    end
    uslot=_new_slot(isempty(name_model) ? "fine" : name_model)
    uslot.mesh=mixed
    push!(MODEL_SLOTS,uslot)
    # Unique mesh nodes reached through element connectivity, in first-seen
    # order (upstream collects `MVertex*` through every entity's elements).
    seen=falses(nnodes(mixed));order=Int[]
    for (msh,cells,_) in _cache_native_blocks(mixed)
        for cell in axes(cells,2),node in axes(cells,1)
            idx=Int(cells[node,cell])
            seen[idx] && continue
            seen[idx]=true;push!(order,idx)
        end
    end
    dataR=Vector{Float64}(undef,length(order))
    dataI=Vector{Float64}(undef,length(order))
    for (i,idx) in enumerate(order)
        x=mixed.coords[1,idx];y=mixed.coords[2,idx];z=mixed.coords[3,idx]
        phi=_plugin_scalar_probe(pview,(x,y,z),0,caller)
        ar=_plugin_scalar_probe(aview,(x,y,z),0,caller)
        ai=_plugin_scalar_probe(aview,(x,y,z),1,caller)
        dataR[i]=ar*cos(k*phi)-ai*sin(k*phi)
        dataI[i]=ar*sin(k*phi)+ai*cos(k*phi)
    end
    tags=UInt64[_cache_node_tag(mixed,idx) for idx in order]
    extent=_native_model_view_extent(tags,
        _native_view_model_count(uslot,"NodeData"))
    counts=fill(Int32(1),length(tags))
    vtag=_fresh_view_tag()
    view=_SessionView("FieldFromAPhi")
    view.model_kind="NodeData"
    view.model_step_ranges=[_NativeViewStepRange(0,1,uslot,1)]
    view.model_data[(0,"NodeData")]=(tags,dataR,1)
    view.model_data[(1,"NodeData")]=(tags,dataI,1)
    for (s,t) in ((0,0.0),(1,1.0))
        view.model_contexts[(s,"NodeData")]=(
            model_slot=uslot,time=t,counts=counts,extent=extent,
            partitions=Set{Int}([0]),
            interpolation_msh=_plugin_mesh_msh_types(mixed))
    end
    SESSION_VIEWS[vtag]=view
    VIEW_TAG_MAX[]=max(VIEW_TAG_MAX[],vtag)
    return nothing
end

# ---- ExtractEdges ----------------------------------------------------------------------

# Upstream `Plugin(ExtractEdges)` scans the current model's surface-triangle
# mesh (every face's `triangles`), builds edge→triangle adjacency
# (`buildEdgeToTriangle` keeps the first and last element per edge) and emits
# scalar `SL` records valued {1,1}: first the two-triangle edges whose normal
# dihedral `atan2(|n1×n2|, n1·n2)` exceeds `Angle` degrees, sorted by
# descending angle, then every one-triangle (boundary) edge when
# `IncludeBoundary` is set. Degenerate normals contribute angle 0.
function _plugin_exec_extractedges(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    _slot_capture_locked!()
    index=_current_slot_index_locked()
    index==0 && throw(ArgumentError("$caller: no current model"))
    slot=MODEL_SLOTS[index]
    mesh=slot.mesh
    tris=NTuple{3,Int}[]
    if mesh!==nothing
        for (msh,cells,_) in _cache_native_blocks(mesh)
            msh_spec(msh).family===:tri || continue
            size(cells,1)>=3 || continue
            for cell in axes(cells,2)
                push!(tris,(Int(cells[1,cell]),Int(cells[2,cell]),
                            Int(cells[3,cell])))
            end
        end
    end
    isempty(tris) && throw(ArgumentError(
        "$caller: No triangles in mesh to extract edges from"))
    normals=Vector{NTuple{3,Float64}}(undef,length(tris))
    for (i,t) in enumerate(tris)
        ax=mesh.coords[1,t[2]]-mesh.coords[1,t[1]]
        ay=mesh.coords[2,t[2]]-mesh.coords[2,t[1]]
        az=mesh.coords[3,t[2]]-mesh.coords[3,t[1]]
        bx=mesh.coords[1,t[3]]-mesh.coords[1,t[1]]
        by=mesh.coords[2,t[3]]-mesh.coords[2,t[1]]
        bz=mesh.coords[3,t[3]]-mesh.coords[3,t[1]]
        nx=ay*bz-az*by;ny=az*bx-ax*bz;nz=ax*by-ay*bx
        nrm=sqrt(nx*nx+ny*ny+nz*nz)
        normals[i]=nrm>0 ? (nx/nrm,ny/nrm,nz/nrm) : (0.0,0.0,0.0)
    end
    adj=Dict{NTuple{2,Int},Tuple{Int,Int}}()
    for (ti,t) in enumerate(tris)
        for (a,b) in ((1,2),(2,3),(3,1))
            key=t[a]<t[b] ? (t[a],t[b]) : (t[b],t[a])
            hit=get(adj,key,nothing)
            adj[key]=hit===nothing ? (ti,0) : (hit[1],ti)
        end
    end
    detected=Tuple{Float64,NTuple{2,Int}}[];lonly=NTuple{2,Int}[]
    for (key,(t1,t2)) in adj
        if t2==0
            push!(lonly,key);continue
        end
        n1=normals[t1];n2=normals[t2]
        cx=n1[2]*n2[3]-n1[3]*n2[2]
        cy=n1[3]*n2[1]-n1[1]*n2[3]
        cz=n1[1]*n2[2]-n1[2]*n2[1]
        angle=atan(sqrt(cx*cx+cy*cy+cz*cz),n1[1]*n2[1]+n1[2]*n2[2]+n1[3]*n2[3])
        push!(detected,(angle,key))
    end
    sort!(detected;by=first,rev=true)
    sort!(lonly)
    threshold=numbers["Angle"]/180*pi
    flat=Float64[];ne=0
    for (angle,(v1,v2)) in detected
        angle<=threshold && break
        for axis in 1:3,v in (v1,v2)
            push!(flat,mesh.coords[axis,v])
        end
        push!(flat,1.0);push!(flat,1.0);ne+=1
    end
    if !iszero(numbers["IncludeBoundary"])
        for (v1,v2) in lonly
            for axis in 1:3,v in (v1,v2)
                push!(flat,mesh.coords[axis,v])
            end
            push!(flat,1.0);push!(flat,1.0);ne+=1
        end
    end
    _plugin_publish_list_view!("ExtractEdges",["SL"],Dict("SL"=>ne),
        Dict("SL"=>flat),[0.0];caller)
    return nothing
end

# ---- Vorticity (Tessella extension) ----------------------------------------------------

# Gmsh 4.15.2 registers `Curl`/`Divergence` but no `Vorticity` plugin; this
# registered alias keeps the requested name as a documented Tessella
# extension: the same corner-gradient curl as `Curl`, emitted under the
# `*_Vorticity` name.
function _plugin_exec_vorticity(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    tag,view=_plugin_view(numbers,caller)
    scanned=_plugin_curl_scan(view,caller,"Vorticity")
    scanned===nothing && return nothing
    out,steps,times=scanned
    _plugin_extract_emit(view.name*"_Vorticity",out,-1,steps,times,caller)
    return nothing
end

# ---- batch-3 shared machinery --------------------------------------------------
#
# The upstream `OctreePost` keeps one octree per record shape *and* per
# component class (scalar/vector/tensor), searching tets, hexes, prisms,
# pyramids, triangles, quadrangles, lines then points — the highest-dimension
# match wins and only elements whose component count matches the probed class
# can supply values. `_plugin_probe_class` reproduces that contract over the
# session view element list: descending search rank with stable record order.

# Upstream `OctreePost::search{Scalar,Vector,Tensor}` probes one octree per
# record family in fixed order — tets, hexes, prisms, pyramids, triangles,
# quadrangles, lines, then points — so the match rank is element TYPE, not
# merely dimension.
const _PLUGIN_SEARCH_RANK=(tetrahedron=8,hexahedron=7,prism=6,pyramid=5,
                           triangle=4,quadrangle=3,line=2,point=1)

function _plugin_sorted_elements(elements)
    return sort(elements;by=e->_PLUGIN_SEARCH_RANK[e.kind],rev=true,
                alg=Base.Sort.DEFAULT_STABLE)
end

# Sample `ncomp`-component data at `p` for `steps` steps. Missing coverage and
# wrong-class hits contribute zeros like the upstream searchScalar/Vector/
# Tensor zero-filled arrays.
function _plugin_class_sample(view,elements,ncomp::Int,p,steps::Int,caller)
    vals=zeros(ncomp*steps)
    if view.model_kind===nothing
        for element in elements
            element.ncomp==ncomp || continue
            N=_plugin_local_weights(element.geo,element.kind,p)
            N===nothing && continue
            va=size(element.vcoords,2)
            for s in 1:steps,node in 1:min(va,length(N)),comp in 1:ncomp
                vals[(s-1)*ncomp+comp]=muladd(N[node],
                    element.vals[(s-1)*va*ncomp+(node-1)*ncomp+comp],
                    vals[(s-1)*ncomp+comp])
            end
            return vals
        end
        return vals
    end
    # Model views carry one component width; only the matching class probes.
    isempty(elements) && return vals
    first(elements).ncomp==ncomp || return vals
    accs,_=_plugin_probe_model(view,p,caller)
    for (i,acc) in enumerate(accs)
        i>steps && break
        for comp in 1:min(ncomp,length(acc))
            vals[(i-1)*ncomp+comp]=acc[comp]
        end
    end
    return vals
end

# Component classes carried by the view (the upstream `getNumScalars`/
# `getNumVectors`/`getNumTensors` gates): only 1/3/9-component records feed
# the matching probe class.
function _plugin_view_classes(elements)
    classes=Set{Int}()
    for element in elements
        element.ncomp in (1,3,9) && push!(classes,element.ncomp)
    end
    return sort!(collect(classes))
end

# Push one `head`-class record: x-block/y-block/z-block then per-step per-node
# values, the `_POS_RECORD_TYPES` payload convention.
function _plugin_push_record!(flat,pnts,nodes,vals,ncomp::Int,steps::Int)
    for axis in 1:3,idx in nodes
        push!(flat,pnts[idx isa Tuple ? CartesianIndex(idx) : idx][axis])
    end
    for s in 1:steps,idx in nodes,comp in 1:ncomp
        push!(flat,
              vals[idx isa Tuple ? CartesianIndex(idx) : idx][(s-1)*ncomp+comp])
    end
    return flat
end

# `NamedTuple` of ordered output blocks shared by the grid plugins: `flat`
# accumulates per record-type payloads keyed by `head*shape`.
function _plugin_grid_blocks()
    return (blocks=Dict{String,Vector{Float64}}(),
            counts=Dict{String,Int}(),order=String[])
end
function _plugin_grid_push!(acc,kind::Symbol,ncomp::Int,pnts,nodes,vals,steps)
    type=_plugin_record_type(kind,ncomp)
    type===nothing && return acc
    haskey(acc.blocks,type) ||
        (acc.blocks[type]=Float64[];acc.counts[type]=0;push!(acc.order,type))
    _plugin_push_record!(acc.blocks[type],pnts,nodes,vals,ncomp,steps)
    acc.counts[type]+=1
    return acc
end

# ---- CutBox ---------------------------------------------------------------------
#
# Upstream `Plugin(CutBox)` samples `View` on the Nu×Nv×Nw grid of the
# parametric box through P0,P1,P2,P3 (`p = P0 + u(P1-P0) + v(P2-P0) +
# w(P3-P0)`, `u=i/(Nu-1)` — 0 when Nu==1 — and likewise v,w). `ConnectPoints`
# selects connected output: lines along the single active axis when two counts
# equal 1, quads on the degenerate plane when one count equals 1, hexahedra for
# a full 3-D grid, and plain points otherwise. `Boundary` keeps only the
# boundary records (points → six faces, connected degenerate grids → their
# boundary polylines/points, 3-D grids → six outward-facing quadrangle skins).
# Component classes are probed independently like the upstream scalar/vector/
# tensor octrees. Upstream `addInView` reads `vals[0][i][getNbW()]` past the
# last slot for the i-lines of the NbU==1 boundary case — an out-of-bounds bug
# reproduced here with the intended NbW-1 slot.

function _plugin_cutbox_addinview!(acc,connect::Bool,boundary::Bool,nu,nv,nw,
                                   pnts::Array{NTuple{3,Float64},3},
                                   vals::Array{Vector{Float64},3},
                                   ncomp::Int,steps::Int)
    pushpnt(idx)=_plugin_grid_push!(acc,:point,ncomp,pnts,(idx,),vals,steps)
    pushln(a,b)=_plugin_grid_push!(acc,:line,ncomp,pnts,(a,b),vals,steps)
    pushqd(a,b,c,d)=_plugin_grid_push!(acc,:quadrangle,ncomp,pnts,(a,b,c,d),
                                       vals,steps)
    pushhex(a,b,c,d,e,f,g,h)=_plugin_grid_push!(acc,:hexahedron,ncomp,pnts,
                                              (a,b,c,d,e,f,g,h),vals,steps)
    if !connect || (nu==1 && nv==1 && nw==1)
        if !boundary
            for i in 1:nu,j in 1:nv,k in 1:nw
                pushpnt((i,j,k))
            end
        else
            for i in 1:nu,j in 1:nv
                pushpnt((i,j,1));pushpnt((i,j,nw))
            end
            for i in 1:nu,j in 1:nw
                pushpnt((i,1,j));pushpnt((i,nv,j))
            end
            for i in 1:nv,j in 1:nw
                pushpnt((1,i,j));pushpnt((nu,i,j))
            end
        end
        return acc
    end
    if nu==1 && nv==1
        if !boundary
            for k in 1:nw-1;pushln((1,1,k),(1,1,k+1));end
        else
            pushpnt((1,1,1));pushpnt((1,1,nw))
        end
    elseif nu==1 && nw==1
        if !boundary
            for j in 1:nv-1;pushln((1,j,1),(1,j+1,1));end
        else
            pushpnt((1,1,1));pushpnt((1,nv,1))
        end
    elseif nv==1 && nw==1
        if !boundary
            for i in 1:nu-1;pushln((i,1,1),(i+1,1,1));end
        else
            pushpnt((1,1,1));pushpnt((nu,1,1))
        end
    elseif nu==1
        if !boundary
            for j in 1:nv-1,k in 1:nw-1
                pushqd((1,j,k),(1,j+1,k),(1,j+1,k+1),(1,j,k+1))
            end
        else
            for j in 1:nv-1
                pushln((1,j,1),(1,j+1,1));pushln((1,j,nw),(1,j+1,nw))
            end
            for k in 1:nw-1
                pushln((1,1,k),(1,1,k+1));pushln((1,nv,k),(1,nv,k+1))
            end
        end
    elseif nv==1
        if !boundary
            for i in 1:nu-1,k in 1:nw-1
                pushqd((i,1,k),(i+1,1,k),(i+1,1,k+1),(i,1,k+1))
            end
        else
            for i in 1:nu-1
                pushln((i,1,1),(i+1,1,1));pushln((i,1,nw),(i+1,1,nw))
            end
            for k in 1:nw-1
                pushln((1,1,k),(1,1,k+1));pushln((nu,1,k),(nu,1,k+1))
            end
        end
    elseif nw==1
        if !boundary
            for i in 1:nu-1,j in 1:nv-1
                pushqd((i,j,1),(i+1,j,1),(i+1,j+1,1),(i,j+1,1))
            end
        else
            for i in 1:nu-1
                pushln((i,1,1),(i+1,1,1));pushln((i,nv,1),(i+1,nv,1))
            end
            for j in 1:nv-1
                pushln((1,j,1),(1,j+1,1));pushln((nu,j,1),(nu,j+1,1))
            end
        end
    else
        if !boundary
            for i in 1:nu-1,j in 1:nv-1,k in 1:nw-1
                pushhex((i,j,k),(i+1,j,k),(i+1,j+1,k),(i,j+1,k),
                        (i,j,k+1),(i+1,j,k+1),(i+1,j+1,k+1),(i,j+1,k+1))
            end
        else
            # Six faces with forced exterior normals, UV then UW then VW.
            for i in 1:nu-1,j in 1:nv-1
                pushqd((i,j,1),(i,j+1,1),(i+1,j+1,1),(i+1,j,1))
                pushqd((i,j,nw),(i+1,j,nw),(i+1,j+1,nw),(i,j+1,nw))
            end
            for i in 1:nu-1,k in 1:nw-1
                pushqd((i,1,k),(i+1,1,k),(i+1,1,k+1),(i,1,k+1))
                pushqd((i,nv,k),(i,nv,k+1),(i+1,nv,k+1),(i+1,nv,k))
            end
            for j in 1:nv-1,k in 1:nw-1
                pushqd((1,j,k),(1,j,k+1),(1,j+1,k+1),(1,j+1,k))
                pushqd((nu,j,k),(nu,j+1,k),(nu,j+1,k+1),(nu,j,k+1))
            end
        end
    end
    return acc
end

function _plugin_exec_cutbox(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    tag,view=_plugin_view(numbers,caller)
    nu=_mesh_query_integer(Int(numbers["NumPointsU"]),caller,"NumPointsU")
    nv=_mesh_query_integer(Int(numbers["NumPointsV"]),caller,"NumPointsV")
    nw=_mesh_query_integer(Int(numbers["NumPointsW"]),caller,"NumPointsW")
    (nu<=0 || nv<=0 || nw<=0) && return nothing
    connect=!iszero(numbers["ConnectPoints"])
    boundary=!iszero(numbers["Boundary"])
    p0=(numbers["X0"],numbers["Y0"],numbers["Z0"])
    p1=(numbers["X1"],numbers["Y1"],numbers["Z1"])
    p2=(numbers["X2"],numbers["Y2"],numbers["Z2"])
    p3=(numbers["X3"],numbers["Y3"],numbers["Z3"])
    pnts=Array{NTuple{3,Float64}}(undef,nu,nv,nw)
    for i in 1:nu,j in 1:nv,k in 1:nw
        u=nu>1 ? (i-1)/(nu-1) : 0.0
        v=nv>1 ? (j-1)/(nv-1) : 0.0
        w=nw>1 ? (k-1)/(nw-1) : 0.0
        pnts[i,j,k]=ntuple(a->p0[a]+u*(p1[a]-p0[a])+v*(p2[a]-p0[a])+
                               w*(p3[a]-p0[a]),3)
    end
    elements,nsteps,times=_plugin_elements_all_steps(view,caller)
    sorted=_plugin_sorted_elements(elements)
    acc=_plugin_grid_blocks()
    for ncomp in _plugin_view_classes(elements)
        vals=Array{Vector{Float64}}(undef,nu,nv,nw)
        for idx in eachindex(vals)
            vals[idx]=_plugin_class_sample(view,sorted,ncomp,pnts[idx],
                                           nsteps,caller)
        end
        _plugin_cutbox_addinview!(acc,connect,boundary,nu,nv,nw,
                                  pnts,vals,ncomp,nsteps)
    end
    _plugin_publish_list_view!(view.name*"_CutBox",acc.order,acc.counts,
                               acc.blocks,times;caller)
    return nothing
end

# ---- CutGrid ---------------------------------------------------------------------
#
# Upstream `Plugin(CutGrid)` samples `View` on the Nu×Nv grid of the
# parametric quad through P0,P1,P2. `ConnectPoints` emits U-lines when Nv==1,
# V-lines when Nu==1, quadrangles for a full grid, and plain points otherwise.

function _plugin_exec_cutgrid(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    tag,view=_plugin_view(numbers,caller)
    nu=_mesh_query_integer(Int(numbers["NumPointsU"]),caller,"NumPointsU")
    nv=_mesh_query_integer(Int(numbers["NumPointsV"]),caller,"NumPointsV")
    (nu<=0 || nv<=0) && return nothing
    connect=!iszero(numbers["ConnectPoints"])
    p0=(numbers["X0"],numbers["Y0"],numbers["Z0"])
    p1=(numbers["X1"],numbers["Y1"],numbers["Z1"])
    p2=(numbers["X2"],numbers["Y2"],numbers["Z2"])
    pnts=Matrix{NTuple{3,Float64}}(undef,nu,nv)
    for i in 1:nu,j in 1:nv
        u=nu>1 ? (i-1)/(nu-1) : 0.0
        v=nv>1 ? (j-1)/(nv-1) : 0.0
        pnts[i,j]=ntuple(a->p0[a]+u*(p1[a]-p0[a])+v*(p2[a]-p0[a]),3)
    end
    elements,nsteps,times=_plugin_elements_all_steps(view,caller)
    sorted=_plugin_sorted_elements(elements)
    acc=_plugin_grid_blocks()
    for ncomp in _plugin_view_classes(elements)
        vals=Matrix{Vector{Float64}}(undef,nu,nv)
        for idx in eachindex(vals)
            vals[idx]=_plugin_class_sample(view,sorted,ncomp,pnts[idx],
                                           nsteps,caller)
        end
        if !connect || (nu==1 && nv==1)
            for i in 1:nu,j in 1:nv
                _plugin_grid_push!(acc,:point,ncomp,pnts,((i,j),),vals,nsteps)
            end
        elseif nu==1
            for j in 1:nv-1
                _plugin_grid_push!(acc,:line,ncomp,pnts,((1,j),(1,j+1)),
                                   vals,nsteps)
            end
        elseif nv==1
            for i in 1:nu-1
                _plugin_grid_push!(acc,:line,ncomp,pnts,((i,1),(i+1,1)),
                                   vals,nsteps)
            end
        else
            for i in 1:nu-1,j in 1:nv-1
                _plugin_grid_push!(acc,:quadrangle,ncomp,pnts,
                                   ((i,j),(i+1,j),(i+1,j+1),(i,j+1)),
                                   vals,nsteps)
            end
        end
    end
    _plugin_publish_list_view!(view.name*"_CutGrid",acc.order,acc.counts,
                               acc.blocks,times;caller)
    return nothing
end

# ---- CutParametric ------------------------------------------------------------------
#
# Upstream `Plugin(CutParametric)` evaluates the X/Y/Z expressions over the
# MinU..MaxU × MinV..MaxV grid (i*NbV+j row-major flattening), probes `View`
# at each point and emits per-class points, consecutive lines when
# `ConnectPoints` on a degenerate grid, or quadrangles on a full grid.
# The Tessella evaluator binds `x`,`y`,`z`, so `u`,`v` rewrite to `x`,`y`.

function _plugin_uv_rewrite(expr::AbstractString)
    out=IOBuffer();i=firstindex(expr);n=lastindex(expr)
    while i<=n
        c=expr[i]
        if _plugin_ident_char(c) && !(c in ('+','-'))
            j=i
            while j<=n && _plugin_ident_char(expr[j]);j=nextind(expr,j);end
            ident=expr[i:prevind(expr,j)]
            write(out,ident=="u" ? "x" : ident=="v" ? "y" : ident)
            i=j;continue
        end
        write(out,c);i=nextind(expr,i)
    end
    return String(take!(out))
end

function _plugin_uv_eval(prog,u,v)
    return try
        SizeField._eval_matheval_entity(prog,u,v,0.0,(),nothing)
    catch err
        err isa InterruptException && rethrow()
        return nothing
    end
end

function _plugin_exec_cutparametric(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers;strings=spec.strings
    tag,view=_plugin_view(numbers,caller)
    nu=_mesh_query_integer(Int(numbers["NumPointsU"]),caller,"NumPointsU")
    nv=_mesh_query_integer(Int(numbers["NumPointsV"]),caller,"NumPointsV")
    connect=!iszero(numbers["ConnectPoints"])
    (nu<2 && nv<2) && (connect=false)
    progs=[try
        SizeField.parse_matheval(_plugin_uv_rewrite(strings[name]))
    catch err
        err isa InterruptException && rethrow()
        nothing
    end for name in ("X","Y","Z")]
    nu>=0 && nv>=0 || return nothing
    minu=numbers["MinU"];maxu=numbers["MaxU"]
    minv=numbers["MinV"];maxv=numbers["MaxV"]
    pnts=Matrix{NTuple{3,Float64}}(undef,nu,nv)
    # Upstream `mathEvaluator::eval` is all-or-nothing per grid point and keeps
    # the previous `res` on failure (zero-initialized on first run).
    res=(0.0,0.0,0.0)
    for i in 1:nu
        u=nu>1 ? minu+(maxu-minu)*(i-1)/(nu-1) : minu
        for j in 1:nv
            v=nv>1 ? minv+(maxv-minv)*(j-1)/(nv-1) : minv
            x=progs[1]===nothing ? nothing : _plugin_uv_eval(progs[1],u,v)
            y=progs[2]===nothing ? nothing : _plugin_uv_eval(progs[2],u,v)
            z=progs[3]===nothing ? nothing : _plugin_uv_eval(progs[3],u,v)
            if x!==nothing && y!==nothing && z!==nothing
                res=(x,y,z)
            end
            pnts[i,j]=res
        end
    end
    elements,nsteps,times=_plugin_elements_all_steps(view,caller)
    sorted=_plugin_sorted_elements(elements)
    acc=_plugin_grid_blocks()
    for ncomp in _plugin_view_classes(elements)
        vals=Matrix{Vector{Float64}}(undef,nu,nv)
        for idx in eachindex(vals)
            vals[idx]=_plugin_class_sample(view,sorted,ncomp,pnts[idx],
                                           nsteps,caller)
        end
        if nu==1 || nv==1 || !connect
            for i in 1:nu,j in 1:nv
                idx=(i-1)*nv+j
                if connect
                    # Upstream emits no record for i==0; it only seeds res0.
                    idx>1 || continue
                    _plugin_grid_push!(acc,:line,ncomp,pnts,
                                       (CartesianIndices(pnts)[idx-1],
                                        CartesianIndices(pnts)[idx]),
                                       vals,nsteps)
                else
                    _plugin_grid_push!(acc,:point,ncomp,pnts,
                                       (CartesianIndices(pnts)[idx],),
                                       vals,nsteps)
                end
            end
        else
            for i in 1:nu-1,j in 1:nv-1
                _plugin_grid_push!(acc,:quadrangle,ncomp,pnts,
                                   ((i,j),(i,j+1),(i+1,j+1),(i+1,j)),
                                   vals,nsteps)
            end
        end
    end
    _plugin_publish_list_view!(view.name*"_CutParametric",acc.order,
                               acc.counts,acc.blocks,times;caller)
    return nothing
end

# ---- StreamLines / Particles ------------------------------------------------------
#
# Upstream `Plugin(StreamLines)` seeds `Nu×Nv` particles on the P0,P1,P2 grid
# and integrates `dX/dt = V` with the 4-stage scheme
#   X1 = X + a1·DT·V(X)    X2 = X + a2·DT·V(X1)    X3 = X + a3·DT·V(X2)
#   X4 = X + a4·DT·V(X3)
#   X += b1(X1-X) + b2(X2-X) + b3(X3-X) + b4(X4-X)
# with a1=a2=0.5, a3=a4=1, b1=1/3, b2=2/3, b3=1/3, b4=1/6. Without `OtherView`
# each seed emits one VP record whose `MaxIter` vector slots hold the running
# displacement; with a scalar `OtherView` each step emits an SL segment valued
# by the other view. `TimeStep < 0` tracks the current physical time through
# the view's step table.
#
# `Plugin(Particles)` solves `A2 X'' + A1 X' + A0 X = V` with the constant-
# acceleration Newmark scheme and emits the same displacement VP records.

function _plugin_seed_point(numbers,i,j,nu,nv)
    u=nu>1 ? (i-1)/(nu-1) : 0.0
    v=nv>1 ? (j-1)/(nv-1) : 0.0
    p0=(numbers["X0"],numbers["Y0"],numbers["Z0"])
    p1=(numbers["X1"],numbers["Y1"],numbers["Z1"])
    p2=(numbers["X2"],numbers["Y2"],numbers["Z2"])
    return ntuple(a->p0[a]+u*(p1[a]-p0[a])+v*(p2[a]-p0[a]),3)
end

# Vector-field sample of `view` at `p` restricted to 3-component elements at
# 0-based step index `stepnum` (the upstream `searchVector(..., step)` scan).
function _plugin_vector_sample(view,elements,p,stepnum,caller)
    if view.model_kind===nothing
        for element in elements
            element.ncomp==3 || continue
            N=_plugin_local_weights(element.geo,element.kind,p)
            N===nothing && continue
            va=size(element.vcoords,2)
            acc=(0.0,0.0,0.0)
            for node in 1:min(va,length(N)),comp in 1:3
                acc=Base.setindex(acc,
                    muladd(N[node],
                        element.vals[stepnum*va*3+(node-1)*3+comp],
                        acc[comp]),comp)
            end
            return acc
        end
        return (0.0,0.0,0.0)
    end
    isempty(elements) && return (0.0,0.0,0.0)
    first(elements).ncomp==3 || return (0.0,0.0,0.0)
    accs,step_nums=_plugin_probe_model(view,p,caller)
    i=findfirst(==(stepnum),step_nums)
    i===nothing && return (0.0,0.0,0.0)
    acc=accs[i]
    return (acc[1],acc[2],acc[3])
end

function _plugin_exec_streamlines(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    tag,view=_plugin_view(numbers,caller)
    nu=_mesh_query_integer(Int(numbers["NumPointsU"]),caller,"NumPointsU")
    nv=_mesh_query_integer(Int(numbers["NumPointsV"]),caller,"NumPointsV")
    (nu<=0 || nv<=0) && return nothing
    dt=numbers["DT"]
    maxiter=_mesh_query_integer(Int(numbers["MaxIter"]),caller,"MaxIter")
    maxiter>=0 || throw(ArgumentError("$caller: MaxIter must be non-negative"))
    timestep=_mesh_query_integer(Int(numbers["TimeStep"]),caller,"TimeStep")
    elements,nsteps,times=_plugin_elements_all_steps(view,caller)
    if timestep>=nsteps
        throw(ArgumentError(
            "$caller: Invalid time step ($timestep) in view $tag"))
    end
    sorted=_plugin_sorted_elements(elements)
    # Optional scalar `OtherView` probes its own per-step time table.
    oview=nothing;osteps=0
    if numbers["OtherView"]>=0
        _,oview=_plugin_view_tag(numbers,"OtherView",caller)
        _,osteps,_=_plugin_elements_all_steps(oview,caller)
    end
    b1,b2,b3,b4=1/3,2/3,1/3,1/6
    a1,a2,a3,a4=0.5,0.5,1.0,1.0
    flat=Float64[]
    if oview===nothing
        records=Vector{NTuple{3,Float64}}(undef,0)
        payload=Float64[]
        for i in 1:nu,j in 1:nv
            x=collect(_plugin_seed_point(numbers,i,j,nu,nv))
            xinit=copy(x)
            push!(records,Tuple(xinit))
            for iter in 1:maxiter
                cur=timestep<0 ? let
                    t0=isempty(times) ? 0.0 : first(times)
                    current=t0+dt*(iter-1)
                    s=0
                    while s<nsteps-1 && current>
                        0.5*(_plugin_step_time(view,s)+
                             _plugin_step_time(view,s+1))
                        s+=1
                    end
                    s
                end : timestep
                v=_plugin_vector_sample(view,sorted,(x[1],x[2],x[3]),cur,caller)
                x1=x.+dt.*a1.*collect(v)
                v=_plugin_vector_sample(view,sorted,(x1[1],x1[2],x1[3]),cur,caller)
                x2=x.+dt.*a2.*collect(v)
                v=_plugin_vector_sample(view,sorted,(x2[1],x2[2],x2[3]),cur,caller)
                x3=x.+dt.*a3.*collect(v)
                v=_plugin_vector_sample(view,sorted,(x3[1],x3[2],x3[3]),cur,caller)
                x4=x.+dt.*a4.*collect(v)
                x.+=b1.*(x1.-x).+b2.*(x2.-x).+b3.*(x3.-x).+b4.*(x4.-x)
                append!(payload,x.-xinit)
            end
        end
        for (k,p) in enumerate(records)
            push!(flat,p[1],p[2],p[3])
            base=(k-1)*3*maxiter
            append!(flat,@view(payload[base+1:base+3*maxiter]))
        end
        # Upstream pushes `currentT = T0 + DT*iter` into the view's Time table
        # once per iteration per seed only when `TimeStep < 0`.
        t0=isempty(times) ? 0.0 : _plugin_step_time(view,0)
        outtimes=timestep<0 ?
            Float64[t0+dt*(iter-1) for _ in 1:nu*nv for iter in 1:maxiter] :
            Float64[]
        _plugin_publish_list_view!(view.name*"_StreamLines",["VP"],
            Dict("VP"=>length(records)),Dict("VP"=>flat),outtimes;caller)
        return nothing
    end
    oelements=_plugin_sorted_elements(
        _plugin_elements_all_steps(oview,caller)[1])
    for i in 1:nu,j in 1:nv
        x=collect(_plugin_seed_point(numbers,i,j,nu,nv))
        val2=_plugin_class_sample(oview,oelements,1,
                                  (x[1],x[2],x[3]),osteps,caller)
        for iter in 1:maxiter
            xprev=copy(x)
            cur=timestep<0 ? let
                t0=isempty(times) ? 0.0 : first(times)
                current=t0+dt*(iter-1)
                s=0
                while s<nsteps-1 && current>
                    0.5*(_plugin_step_time(view,s)+_plugin_step_time(view,s+1))
                    s+=1
                end
                s
            end : timestep
            v=_plugin_vector_sample(view,sorted,(x[1],x[2],x[3]),cur,caller)
            x1=x.+dt.*a1.*collect(v)
            v=_plugin_vector_sample(view,sorted,(x1[1],x1[2],x1[3]),cur,caller)
            x2=x.+dt.*a2.*collect(v)
            v=_plugin_vector_sample(view,sorted,(x2[1],x2[2],x2[3]),cur,caller)
            x3=x.+dt.*a3.*collect(v)
            v=_plugin_vector_sample(view,sorted,(x3[1],x3[2],x3[3]),cur,caller)
            x4=x.+dt.*a4.*collect(v)
            x.+=b1.*(x1.-x).+b2.*(x2.-x).+b3.*(x3.-x).+b4.*(x4.-x)
            # SL record: x-block y-block z-block, then ALL steps of the
            # scalar at the previous point followed by ALL steps at the new
            # point (the upstream `for k` / `for k` push order).
            append!(flat,(xprev[1],x[1],xprev[2],x[2],xprev[3],x[3]))
            valnew=_plugin_class_sample(oview,oelements,1,
                                        (x[1],x[2],x[3]),osteps,caller)
            append!(flat,val2)
            append!(flat,valnew)
            val2=valnew
        end
    end
    nseg=nu*nv*maxiter
    t0=isempty(times) ? 0.0 : _plugin_step_time(view,0)
    outtimes=timestep<0 ?
        Float64[t0+dt*(iter-1) for _ in 1:nu*nv for iter in 1:maxiter] :
        Float64[]
    _plugin_publish_list_view!(view.name*"_StreamLines",["SL"],
        Dict("SL"=>nseg),Dict("SL"=>flat),outtimes;caller)
    return nothing
end

function _plugin_exec_particles(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    tag,view=_plugin_view(numbers,caller)
    nu=_mesh_query_integer(Int(numbers["NumPointsU"]),caller,"NumPointsU")
    nv=_mesh_query_integer(Int(numbers["NumPointsV"]),caller,"NumPointsV")
    (nu<=0 || nv<=0) && return nothing
    a2=numbers["A2"];a1=numbers["A1"];a0=numbers["A0"]
    dt=numbers["DT"]
    maxiter=_mesh_query_integer(Int(numbers["MaxIter"]),caller,"MaxIter")
    maxiter>=0 || throw(ArgumentError("$caller: MaxIter must be non-negative"))
    timestep=_mesh_query_integer(Int(numbers["TimeStep"]),caller,"TimeStep")
    elements,nsteps,times=_plugin_elements_all_steps(view,caller)
    if timestep>=nsteps
        _native_diagnostic_emit_locked(:error,
            "Invalid time step ($timestep) in view[$tag]: using 0")
        timestep=0
    end
    sorted=_plugin_sorted_elements(elements)
    # Constant-acceleration Newmark coefficients (γ=1/2, β=1/4).
    γ=0.5;β=0.25
    c1=a2+γ*dt*a1+β*dt*dt*a0
    c2=2a2-(1-2γ)*dt*a1-(0.5+γ-2β)*dt*dt*a0
    c3=-a2-(γ-1)*dt*a1-(0.5-γ+β)*dt*dt*a0
    c4=dt*dt*(β+(0.5+γ-2β)+(0.5-γ+β))
    flat=Float64[]
    nrec=0
    for i in 1:nu,j in 1:nv
        xinit=collect(_plugin_seed_point(numbers,i,j,nu,nv))
        x0=copy(xinit)
        x1=copy(xinit)
        # Per-record layout like upstream `data2->VP`: the seed point then
        # `MaxIter` displacement vectors (X - XINIT).
        push!(flat,xinit[1],xinit[2],xinit[3]);nrec+=1
        for iter in 1:maxiter
            f=_plugin_vector_sample(view,sorted,(x1[1],x1[2],x1[3]),
                                    timestep,caller)
            x=(c2.*x1.+c3.*x0.+c4.*collect(f))./c1
            append!(flat,(x[1]-xinit[1],x[2]-xinit[2],x[3]-xinit[3]))
            x0,x1=x1,x
        end
    end
    _plugin_publish_list_view!(view.name*"_Particles",["VP"],
        Dict("VP"=>nrec),Dict("VP"=>flat),Float64[];caller)
    return nothing
end

# ---- Triangulate / Tetrahedralize -------------------------------------------------
#
# Upstream `Plugin(Triangulate)` gathers the 1-node records, projects them on
# the SVD mean plane (`computeMeanPlaneSimple`: the least-singular right
# vector is the normal, the in-plane axes come from normal × e_i on the first
# vanishing component), jitters by `1e-12*lc*rand()`, 2-D Delaunay, and emits
# the triangles at the projected coordinates carrying all steps. Tessella
# uses the deterministic `delaunay2d` kernel instead of `rand()`.
#
# `Plugin(Tetrahedralize)` gathers 1-node records (≥4), runs 3-D Delaunay and
# emits real tets carrying all steps.

function _plugin_point_records(view,caller)
    elements,steps,_=_plugin_elements_all_steps(view,caller)
    points=_PluginElement[]
    for element in elements
        element.kind===:point || continue
        push!(points,element)
    end
    return points,steps
end

function _plugin_exec_triangulate(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    tag,view=_plugin_view(numbers,caller)
    _plugin_check_single_mesh(view,"Triangulate") || return nothing
    points,steps=_plugin_point_records(view,caller)
    length(points)>=3 || throw(ArgumentError(
        "$caller: Need at least 3 points to triangulate"))
    n=length(points)
    xm=ym=zm=0.0
    for e in points
        xm+=e.geo[1,1];ym+=e.geo[2,1];zm+=e.geo[3,1]
    end
    xm/=n;ym/=n;zm/=n
    U=Matrix{Float64}(undef,n,3)
    for (i,e) in enumerate(points)
        U[i,1]=e.geo[1,1]-xm;U[i,2]=e.geo[2,1]-ym;U[i,3]=e.geo[3,1]-zm
    end
    svd=LinearAlgebra.svd(U)
    Vt=svd.Vt
    # Upstream `computeMeanPlaneSimple` takes `min = 0` iff σ0 is strictly
    # smallest, `1` iff σ1 is strictly smallest, otherwise 2 — ties resolve to
    # the last index.
    s=svd.S
    minidx=(abs(s[1])<abs(s[2]) && abs(s[1])<abs(s[3])) ? 1 :
           (abs(s[2])<abs(s[1]) && abs(s[2])<abs(s[3])) ? 2 : 3
    normal=collect(Vt[minidx,:])
    nrm=_plugin_norm3(Tuple(normal))
    nrm>0 && (normal./=nrm)
    ex=normal[1]==0 ? (1.0,0.0,0.0) : normal[2]==0 ? (0.0,1.0,0.0) : (0.0,0.0,1.0)
    t1=_plugin_cross3(Tuple(normal),ex)
    t1n=_plugin_norm3(t1)
    t1=t1n>0 ? t1./t1n : (1.0,0.0,0.0)
    t2=_plugin_cross3(t1,Tuple(normal))
    t2n=_plugin_norm3(t2)
    t2=t2n>0 ? t2./t2n : (0.0,1.0,0.0)
    # `fillMeanPlane` anchors the projected frame at the dominant-axis
    # intercept of the plane `res·p = d` (d = res·centroid), not the centroid.
    d=normal[1]*xm+normal[2]*ym+normal[3]*zm
    anchor=if abs(normal[1])>=abs(normal[2]) && abs(normal[1])>=abs(normal[3])
        (d/normal[1],0.0,0.0)
    elseif abs(normal[2])>=abs(normal[3])
        (0.0,d/normal[2],0.0)
    else
        (0.0,0.0,d/normal[3])
    end
    us=Float64[LinearAlgebra.dot(
                   (e.geo[1,1]-anchor[1],e.geo[2,1]-anchor[2],
                    e.geo[3,1]-anchor[3]),t1) for e in points]
    vs=Float64[LinearAlgebra.dot(
                   (e.geo[1,1]-anchor[1],e.geo[2,1]-anchor[2],
                    e.geo[3,1]-anchor[3]),t2) for e in points]
    # `delaunay2d` dedups vertices; keep the unique→first-original remap so the
    # emitted triangles carry the original record payloads.
    _,_,remap=Mesh2D.dedup_points(us,vs)
    first_of=Dict{Int,Int}()
    for i in 1:n
        get!(first_of,Int(remap[i]),i)
    end
    T=Mesh2D.delaunay2d(us,vs)
    tris=Tuple{Int,Int,Int}[]
    for t in eachindex(T.alive)
        T.alive[t] || continue
        Mesh2D._is_ghost_tri(T,t) && continue
        push!(tris,ntuple(k->first_of[Int(Mesh2D._vert(T,t,k))],3))
    end
    times=[_plugin_step_time(view,s) for s in 0:steps-1]
    acc=_plugin_grid_blocks()
    # Upstream picks the record class per triangle from the vertices' stored
    # width (9→TT, 3→VT, else ST); mixed vertex widths fall back to scalar.
    for (j,tri) in enumerate(tris)
        widths=[points[tri[k]].ncomp for k in 1:3]
        ncomp=allequal(widths) && widths[1] in (3,9) ? widths[1] : 1
        pmat=Matrix{NTuple{3,Float64}}(undef,3,1)
        vmat=Matrix{Vector{Float64}}(undef,3,1)
        for k in 1:3
            e=points[tri[k]]
            # Upstream emits the original 3-D positions (the 2-D projection
            # only drives the Delaunay connectivity).
            pmat[k,1]=(e.geo[1,1],e.geo[2,1],e.geo[3,1])
            if e.ncomp==ncomp
                vmat[k,1]=e.vals
            else
                # Scalar fallback reads the flat `v[3 + step]` slots like
                # upstream's `v[3 + numComp * step + comp]` (numComp=1).
                vmat[k,1]=Float64[e.vals[min(s,length(e.vals))] for s in 1:steps]
            end
        end
        _plugin_grid_push!(acc,:triangle,ncomp,pmat,
                           ((1,1),(2,1),(3,1)),vmat,steps)
    end
    _plugin_publish_list_view!(view.name*"_Triangulate",acc.order,acc.counts,
                               acc.blocks,times;caller)
    return nothing
end

function _plugin_exec_tetrahedralize(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    tag,view=_plugin_view(numbers,caller)
    _plugin_check_single_mesh(view,"Tetrahedralize") || return nothing
    points,steps=_plugin_point_records(view,caller)
    length(points)>=4 || throw(ArgumentError(
        "$caller: Need at least 4 points to tetrahedralize"))
    xs=Float64[e.geo[1,1] for e in points]
    ys=Float64[e.geo[2,1] for e in points]
    zs=Float64[e.geo[3,1] for e in points]
    _,_,_,remap=Mesh3D._dedup3(xs,ys,zs)
    first_of=Dict{Int,Int}()
    for i in 1:length(points)
        get!(first_of,Int(remap[i]),i)
    end
    T=Mesh3D.delaunay3d(xs,ys,zs)
    tets=Tuple{Int,Int,Int,Int}[]
    for t in eachindex(T.alive)
        (T.alive[t] && !Mesh3D._is_ghost_tet(T,t)) || continue
        # `to_mesh3` convention: swap v3,v4 so exported tets have positive
        # signed volume.
        push!(tets,(first_of[Int(Mesh3D._vert(T,t,1))],
                    first_of[Int(Mesh3D._vert(T,t,2))],
                    first_of[Int(Mesh3D._vert(T,t,4))],
                    first_of[Int(Mesh3D._vert(T,t,3))]))
    end
    times=[_plugin_step_time(view,s) for s in 0:steps-1]
    acc=_plugin_grid_blocks()
    # Upstream decides the record class per tet: TS when all four vertices
    # carry 9*steps values, VS for 3*steps, SS for steps — any other mix is
    # skipped with `Msg::Warning("Skipping unknown type of data")`.
    warned=false
    for tet in tets
        widths=[points[tet[k]].ncomp for k in 1:4]
        ncomp=allequal(widths) && widths[1] in (1,3,9) ? widths[1] : 0
        if ncomp==0
            warned || (_native_diagnostic_emit_locked(:warning,
                "Skipping unknown type of data");warned=true)
            continue
        end
        pmat=Matrix{NTuple{3,Float64}}(undef,4,1)
        vmat=Matrix{Vector{Float64}}(undef,4,1)
        for k in 1:4
            e=points[tet[k]]
            pmat[k,1]=(e.geo[1,1],e.geo[2,1],e.geo[3,1])
            vmat[k,1]=e.vals
        end
        _plugin_grid_push!(acc,:tetrahedron,ncomp,pmat,
                           ((1,1),(2,1),(3,1),(4,1)),vmat,steps)
    end
    _plugin_publish_list_view!(view.name*"_Tetrahedralize",acc.order,
                               acc.counts,acc.blocks,times;caller)
    return nothing
end

# ---- NearestNeighbor ---------------------------------------------------------------
#
# Upstream `Plugin(NearestNeighbor)` rewrites each 1-node record's scalar
# value in place with the distance to its nearest *other* point (ANN k=2
# search, index 1). Tessella scans point records directly.

function _plugin_exec_nearestneighbor(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    tag,view=_plugin_view(numbers,caller)
    view.model_kind===nothing || throw(ArgumentError(
        "$caller: NearestNeighbor requires a list-based view"))
    steps=_plugin_list_steps(view,caller)
    npts=0
    for type in view.list_types
        layout=get(Post._POS_RECORD_TYPES,type,nothing)
        layout===nothing && continue
        layout[2]==1 && (npts+=view.list_counts[type])
    end
    npts>0 || throw(ArgumentError(
        "$caller: View[$tag] contains no points"))
    for type in view.list_types
        layout=get(Post._POS_RECORD_TYPES,type,nothing)
        layout===nothing && continue
        kind,arity,nc=layout
        arity==1 || continue
        data=view.list_data[type];count=view.list_counts[type]
        width=3+nc*steps
        for element in 0:count-1
            base=element*width
            x,y,z=data[base+1],data[base+2],data[base+3]
            best=Inf
            for type2 in view.list_types
                layout2=get(Post._POS_RECORD_TYPES,type2,nothing)
                layout2===nothing && continue
                layout2[2]==1 || continue
                data2=view.list_data[type2];count2=view.list_counts[type2]
                width2=3+layout2[3]*steps
                for element2 in 0:count2-1
                    base2=element2*width2
                    (type2==type && element2==element) && continue
                    dx=data2[base2+1]-x;dy=data2[base2+2]-y
                    dz=data2[base2+3]-z
                    d=dx*dx+dy*dy+dz*dz
                    d<best && (best=d)
                end
            end
            # Upstream `setValue(step=0, nod=0, comp=0)`: only the first
            # component of step 0 receives the distance.
            data[base+3+1]=sqrt(best)
        end
    end
    # Upstream renames the modified view in place.
    view.name=view.name*"_NearestNeighbor"
    return nothing
end

# ---- Scal2Tens ---------------------------------------------------------------------
#
# Upstream `Plugin(Scal2Tens)` packs component 0 of up to `NumberOfComponents`
# scalar views into a `NameNewView` tensor list view at the first set view's
# (element, node, step) slots — the `Scal2Vec` pattern with 9 slots.

function _plugin_exec_scal2tens(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers;strings=spec.strings
    name=strings["NameNewView"]
    ncomp=_mesh_query_integer(Int(numbers["NumberOfComponents"]),
                              caller,"NumberOfComponents")
    1<=ncomp<=9 || throw(ArgumentError(
        "$caller: NumberOfComponents must be in 1:9"))
    cviews=Tuple{Int,Any}[]
    ref=nothing
    # Upstream reads `View comp` options only for comp < numComp.
    for comp in 0:ncomp-1
        raw=numbers["View$comp"]
        if raw<0
            push!(cviews,(0,nothing));continue
        end
        ctag,cview=_plugin_view_tag(numbers,"View$comp",caller)
        push!(cviews,(ctag,cview))
        ref===nothing && (ref=(ctag,cview))
    end
    ref===nothing && throw(ArgumentError(
        "$caller: Scal2Tens plugin could not find any view."))
    _,rview=ref
    elements,steps,times=_plugin_elements_all_steps(rview,caller)
    merged=Vector{Vector{_PluginElement}}(undef,ncomp)
    for comp in 1:ncomp
        merged[comp]=cviews[comp][2]===nothing ? _PluginElement[] :
            _plugin_elements_all_steps(cviews[comp][2],caller)[1]
    end
    out=_PluginElement[]
    for (i,element) in enumerate(elements)
        va=size(element.vcoords,2);ga=size(element.geo,2)
        va==ga || continue
        _plugin_record_type(element.kind,ncomp)===nothing && continue
        vals=zeros(va*ncomp*steps)
        for comp in 1:ncomp
            cviews[comp][2]===nothing && continue
            oe=merged[comp]
            i<=length(oe) || continue
            oelement=oe[i];oa=size(oelement.vcoords,2);onc=oelement.ncomp
            onc>=1 || continue
            osteps=length(oelement.vals)÷(oa*onc)
            # Upstream `getValue` failures contribute the zero-initialized
            # `val` — out-of-range (step, node) reads stay zero.
            for step in 1:min(steps,osteps),node in 1:min(va,oa)
                vals[(step-1)*va*ncomp+(node-1)*ncomp+comp]=
                    oelement.vals[(step-1)*oa*onc+(node-1)*onc+1]
            end
        end
        push!(out,_PluginElement(element.kind,element.dim,element.vcoords,
                                 element.vcoords,vals,ncomp,element.etag))
    end
    _plugin_extract_emit(name,out,-1,steps,times,caller)
    return nothing
end

# ---- Eigenvectors --------------------------------------------------------------------
#
# Upstream `Plugin(Eigenvectors)` decomposes each 9-component record per node
# (real `mat.eig`, ascending eigenvalue order) into the `*_MinEigenvectors`,
# `*_MidEigenvectors` and `*_MaxEigenvectors` vector views, scaled by the
# eigenvalue when `ScaleByEigenvalues` is set. Complex spectra are reported.

function _plugin_exec_eigenvectors(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    tag,view=_plugin_view(numbers,caller)
    _plugin_check_single_mesh(view,"Eigenvectors") || return nothing
    scale=!iszero(numbers["ScaleByEigenvalues"])
    elements,nsteps,times=_plugin_elements_all_steps(view,caller)
    outmin=_PluginElement[];outmid=_PluginElement[];outmax=_PluginElement[]
    nbcomplex=0
    for element in elements
        element.ncomp==9 || continue
        va=size(element.vcoords,2);ga=size(element.geo,2)
        va==ga || continue
        _plugin_record_type(element.kind,3)===nothing && continue
        vmin=Vector{Float64}(undef,va*3*nsteps)
        vmid=Vector{Float64}(undef,va*3*nsteps)
        vmax=Vector{Float64}(undef,va*3*nsteps)
        ok=true
        for s in 1:nsteps,node in 1:va
            m=Matrix{Float64}(undef,3,3)
            for i in 1:3,j in 1:3
                m[i,j]=element.vals[(s-1)*va*9+(node-1)*9+3*(i-1)+j]
            end
            F=try
                LinearAlgebra.eigen(m)
            catch err
                err isa InterruptException && rethrow()
                ok=false;break
            end
            vals=F.values;vecs=F.vectors
            if any(v->imag(v)!=0,vals)
                nbcomplex+=1
            end
            order=sortperm(vals;by=real)
            for (slot,pos) in enumerate(order)
                ev=scale ? real(vals[pos]) : 1.0
                vec=real.(vecs[:,pos])
                n=_plugin_norm3(Tuple(vec))
                # Upstream keeps the raw eigenvector; normalize for a stable
                # direction like the upstream unit right-eigenvectors.
                n>0 && (vec=vec./n)
                target=slot==1 ? vmin : slot==2 ? vmid : vmax
                for a in 1:3
                    target[(s-1)*va*3+(node-1)*3+a]=ev*vec[a]
                end
            end
        end
        ok || continue
        etag=element.etag
        push!(outmin,_PluginElement(element.kind,element.dim,element.vcoords,
                                    element.vcoords,vmin,3,etag))
        push!(outmid,_PluginElement(element.kind,element.dim,element.vcoords,
                                    element.vcoords,vmid,3,etag))
        push!(outmax,_PluginElement(element.kind,element.dim,element.vcoords,
                                    element.vcoords,vmax,3,etag))
    end
    nbcomplex>0 && _native_diagnostic_emit_locked(:error,
        "$nbcomplex tensors have complex eigenvalues")
    name=view.name
    _plugin_extract_emit(name*"_MinEigenvectors",outmin,-1,nsteps,times,caller)
    _plugin_extract_emit(name*"_MidEigenvectors",outmid,-1,nsteps,times,caller)
    _plugin_extract_emit(name*"_MaxEigenvectors",outmax,-1,nsteps,times,caller)
    return nothing
end

# ---- Summation ------------------------------------------------------------------------
#
# Upstream `Plugin(Summation)` sums component-wise values across up to 8 views
# sharing the same grid (same entity/element structure); empty steps are
# skipped and each populated step of the reference view sums all others at the
# same step index.

function _plugin_exec_summation(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers;strings=spec.strings
    gathered=Tuple{Int,Any}[]
    for i in 0:7
        raw=numbers["View $i"]
        if raw<0
            i==0 && throw(ArgumentError(
                "$caller: Summation plugin could not find view 0"))
            continue
        end
        ctag,cview=_plugin_view_tag(numbers,"View $i",caller)
        _plugin_check_single_mesh(cview,"Summation") || return nothing
        push!(gathered,(ctag,cview))
    end
    _,rview=gathered[1]
    elements,nsteps,times=_plugin_elements_all_steps(rview,caller)
    others=[_plugin_elements_all_steps(v,caller)[1]
            for (_,v) in gathered[2:end]]
    # Upstream reports the grid mismatch once per offending view, then keeps
    # summing whatever `getValue` can reach (failures leave zeros).
    for oe in others
        length(oe)==length(elements) || _native_diagnostic_emit_locked(:error,
            "Summation plugin: views based on different grid.")
    end
    out=_PluginElement[]
    for (i,element) in enumerate(elements)
        va=size(element.vcoords,2);nc=element.ncomp
        vals=copy(element.vals)
        for oe in others
            i<=length(oe) || continue
            oelement=oe[i];oa=size(oelement.vcoords,2);onc=oelement.ncomp
            osteps=length(oelement.vals)÷(oa*onc)
            for s in 1:min(nsteps,osteps),node in 1:min(va,oa),
                comp in 1:min(nc,onc)
                vals[(s-1)*va*nc+(node-1)*nc+comp]+=
                    oelement.vals[(s-1)*oa*onc+(node-1)*onc+comp]
            end
        end
        push!(out,_PluginElement(element.kind,element.dim,element.geo,
                                 element.vcoords,vals,nc,element.etag))
    end
    # The upstream default name is `<view0 name>_Summation`.
    name=strings["Resuling View Name"]
    name=="default" && (name=rview.name*"_Summation")
    _plugin_extract_emit(name,out,-1,nsteps,times,caller)
    return nothing
end

# ---- MeshVolume ----------------------------------------------------------------------
#
# Upstream `Plugin(MeshVolume)` prints the summed element measure over the
# requested `Dimension`/`PhysicalGroup` and emits one `MeshVolume` SP view at
# the model bounding-box centre carrying the volume.

function _plugin_exec_meshvolume(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    dim=_mesh_query_integer(Int(numbers["Dimension"]),caller,"Dimension")
    physical=_mesh_query_integer(Int(numbers["PhysicalGroup"]),
                                 caller,"PhysicalGroup")
    _slot_capture_locked!()
    index=_current_slot_index_locked()
    index==0 && throw(ArgumentError("$caller: no current model"))
    slot=MODEL_SLOTS[index]
    mesh=slot.mesh
    mesh===nothing && throw(ArgumentError(
        "$caller: Plugin(MeshVolume) needs a mesh"))
    class=slot.class===nothing ? _cached_classification_locked(mesh) :
        slot.class
    public=_mesh_public_tags(mesh)
    _,_,total=_mesh_element_offsets(mesh)
    total>0 || throw(ArgumentError("$caller: Plugin(MeshVolume) needs a mesh"))
    vol=0.0;nelem=0
    for index in 1:total
        rec=mesh_element_record(mesh,index)
        # Upstream `getEntities(entities, dim)` with dim=-1 returns entities
        # of every dimension.
        (dim<0 || rec.dimension==dim) || continue
        # Upstream collects `groups[physical]` entities unless physical == -1.
        if physical != -1
            entity=_plugin_element_entity(mesh,class,index)
            entity===nothing && continue
            phys=get(class.physicals,entity,Int32[])
            isempty(phys) && continue
            Int32(physical) in phys || continue
        end
        spec_rec=msh_spec(Int(rec.element_type))
        kind=get(_PLUGIN_FAMILY_KIND,spec_rec.family,nothing)
        kind===nothing && continue
        ga=length(rec.node_tags)
        geo=Matrix{Float64}(undef,3,ga)
        for axis in 1:3,node in 1:ga
            geo[axis,node]=mesh.coords[axis,Int(rec.node_tags[node])]
        end
        vol+=_plugin_element_measure(geo,kind)
        nelem+=1
    end
    # Upstream warns but still publishes the (zero-volume) view.
    nelem==0 && _native_diagnostic_emit_locked(:warning,
        "The specified domain is empty")
    _native_diagnostic_emit_locked(:info,
        "Mesh volume (physical $physical | dimension $dim): $vol")
    # Upstream writes the point at `GModel::current()->bounds().center()`.
    lo=fill(Inf,3);hi=fill(-Inf,3)
    if slot.model!==nothing
        try
            bbox=model_bounding_box(slot.model,-1,-1)
            for a in 1:3
                lo[a]=bbox[a];hi[a]=bbox[3+a]
            end
        catch err
            err isa InterruptException && rethrow()
        end
    end
    if any(!isfinite,lo)
        lo=fill(Inf,3);hi=fill(-Inf,3)
        for i in axes(mesh.coords,2),a in 1:3
            v=mesh.coords[a,i]
            v<lo[a] && (lo[a]=v)
            v>hi[a] && (hi[a]=v)
        end
    end
    center=any(!isfinite,lo) ? (0.0,0.0,0.0) : Tuple((lo.+hi)./2)
    _plugin_publish_list_view!("MeshVolume",["SP"],Dict("SP"=>1),
        Dict("SP"=>[center...,vol]),Float64[];caller)
    return nothing
end

# ---- LongitudeLatitude ---------------------------------------------------------------
#
# Upstream `Plugin(LongitudeLatitude)` rewrites node coordinates in place:
# `x = atan2(y,x)`, `y = asin(z/r)`, `z = 0` in degrees — restricted to
# list-based views like the other in-place writers.

function _plugin_exec_longitude_latitude(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    tag,view=_plugin_view(numbers,caller)
    view.model_kind===nothing || throw(ArgumentError(
        "$caller: LongitudeLatitude requires a list-based view"))
    steps=_plugin_list_steps(view,caller)
    # Upstream iterates steps outermost and rewrites coordinates through the
    # shared list slots once per step — for multi-step views the transform is
    # applied repeatedly (upstream quirk preserved here).
    for step in 1:steps
        for type in view.list_types
            layout=get(Post._POS_RECORD_TYPES,type,nothing)
            layout===nothing && continue
            kind,arity,nc=layout
            data=view.list_data[type];count=view.list_counts[type]
            width=3arity+arity*nc*steps
            for element in 0:count-1
                base=element*width
                xmin=Inf;xmax=-Inf
                for node in 1:arity
                    x=data[base+node];y=data[base+arity+node]
                    z=data[base+2arity+node]
                    r=sqrt(x*x+y*y+z*z)
                    x2=atan(y,x)
                    y2=asin(z/r)
                    xmin=min(xmin,x2);xmax=max(xmax,x2)
                    if nc==3
                        vbase=base+3arity+(step-1)*arity*nc+(node-1)*nc
                        vin=(data[vbase+1],data[vbase+2],data[vbase+3])
                        sxy,cxy=sin(x2),cos(x2)
                        syy,cyy=sin(y2),cos(y2)
                        vout=(-sxy*vin[1]+cxy*vin[2],
                              -syy*(cxy*vin[1]+sxy*vin[2])+cyy*vin[3],
                              cyy*(cxy*vin[1]+sxy*vin[2])+syy*vin[3])
                        for a in 1:3
                            data[vbase+a]=vout[a]
                        end
                    end
                    data[base+node]=rad2deg(x2)
                    data[base+arity+node]=rad2deg(y2)
                    data[base+2arity+node]=0.0
                end
                # Antimeridian unwrap: when the element spans more than 180°
                # push the wrapped nodes into [xmax-180°, xmax+180°).
                if xmax-xmin>pi
                    for node in 1:arity
                        x=data[base+node]
                        xmax*180/pi-x>180 && (data[base+node]=x+360)
                    end
                end
            end
        end
    end
    return nothing
end

# ---- ShowNeighborElements ---------------------------------------------------------------
#
# Upstream `Plugin(ShowNeighborElements)` keeps `Element1..5` (public element
# tags) visible plus, for each dim-2/dim-3 entity, hidden elements sharing
# more than `3 - NumLayers` vertices with the seed set. Everything else is
# hidden. Tessella writes `ELEMENT_VISIBILITY` (public tag → flag) — the same
# per-element visibility slot `API.mesh.set_visibility` manipulates.

function _plugin_exec_showneighbor(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    nlayers=_mesh_query_integer(Int(numbers["NumLayers"]),caller,"NumLayers")
    seeds=Set{UInt64}()
    for i in 1:5
        raw=numbers["Element$i"]
        raw==0 && continue
        isinteger(raw) || throw(ArgumentError(
            "$caller: Element$i tag must be an integer"))
        push!(seeds,UInt64(Int(raw)))
    end
    _slot_capture_locked!()
    index=_current_slot_index_locked()
    index==0 && throw(ArgumentError("$caller: no current model"))
    slot=MODEL_SLOTS[index]
    mesh=slot.mesh
    mesh===nothing && throw(ArgumentError("$caller: mesh is empty"))
    class=slot.class===nothing ? _cached_classification_locked(mesh) :
        slot.class
    public=_mesh_public_tags(mesh)
    _,_,total=_mesh_element_offsets(mesh)
    # Per-entity vertex→element map for dim-2/3 elements plus the seed set,
    # mirroring `_init`/`_showLayers`. `_vertices` is an upstream member: it
    # accumulates seed vertices across every entity, so neighbor layers leak
    # over entity boundaries exactly like upstream.
    entitymap=Dict{Tuple{Int,Int},Vector{Int}}()
    for index in 1:total
        rec=mesh_element_record(mesh,index)
        rec.dimension<2 && continue
        entity=_plugin_element_entity(mesh,class,index)
        entity===nothing && continue
        push!(get!(entitymap,entity,Int[]),index)
    end
    entities=sort!(collect(keys(entitymap));by=ent->(ent[1],ent[2]))
    vertices=Set{Int}()
    for entity in entities
        elements=entitymap[entity]
        vert2elem=Dict{Int,Vector{Int}}()
        for index in elements
            rec=mesh_element_record(mesh,index)
            etag=public===nothing ? UInt64(index) :
                get(public.element_tags,index,UInt64(index))
            if etag in seeds
                ELEMENT_VISIBILITY[][etag]=Int32(1)
                for t in rec.node_tags;push!(vertices,Int(t));end
            else
                ELEMENT_VISIBILITY[][etag]=Int32(0)
                for t in rec.node_tags
                    push!(get!(vert2elem,Int(t),Int[]),index)
                end
            end
        end
        (isempty(vertices) || nlayers<1) && continue
        counts=Dict{Int,Int}()
        for v in vertices
            for index in get(vert2elem,v,Int[])
                counts[index]=get(counts,index,0)+1
            end
        end
        for (index,count) in counts
            count>3-nlayers || continue
            etag=public===nothing ? UInt64(index) :
                get(public.element_tags,index,UInt64(index))
            ELEMENT_VISIBILITY[][etag]=Int32(1)
        end
    end
    return nothing
end

# ---- SimplePartition --------------------------------------------------------------------
#
# Upstream `Plugin(SimplePartition)` clears existing partitions, evaluates
# `MappingX/Y/Z` (`t` over [0,1]) into per-axis slice abscissae over the model
# bounding box and assigns each element — barycenter centre — partition
# `kx*Ny*Nz + ky*Nz + kz + 1`. Slices are half-open `(pp[k], pp[k+1]]` except
# the lowest slice which also takes exactly `pp[0]`; a degenerate axis
# (`pp[0]==pp[N]`) accepts every coordinate. Tessella drives the same
# assignment through the explicit `_partition` path so
# `LAST_MESH_PARTITION`/ghost/topology state stays authoritative.

function _plugin_simple_partition_exec(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers;strings=spec.strings
    nsx=_mesh_query_integer(Int(numbers["NumSlicesX"]),caller,"NumSlicesX")
    nsy=_mesh_query_integer(Int(numbers["NumSlicesY"]),caller,"NumSlicesY")
    nsz=_mesh_query_integer(Int(numbers["NumSlicesZ"]),caller,"NumSlicesZ")
    (nsx<1 || nsy<1 || nsz<1) && throw(ArgumentError(
        "$caller: Number of slices should be strictly positive"))
    _slot_capture_locked!()
    index=_current_slot_index_locked()
    index==0 && throw(ArgumentError("$caller: no current model"))
    slot=MODEL_SLOTS[index]
    mesh=slot.mesh
    mesh===nothing && throw(ArgumentError(
        "$caller: Plugin(SimplePartition) requires a mesh"))
    _,_,total=_mesh_element_offsets(mesh)
    total>0 || throw(ArgumentError(
        "$caller: Plugin(SimplePartition) requires a mesh"))
    # Model bounds like upstream `m->bounds()`; pure meshes fall back to the
    # mesh coordinate extents.
    lo=fill(Inf,3);hi=fill(-Inf,3)
    bounds=try
        model_bounding_box(slot.model,-1,-1)
    catch err
        err isa InterruptException && rethrow()
        nothing
    end
    if bounds!==nothing
        # Flat (xmin,ymin,zmin,xmax,ymax,zmax) like upstream `m->bounds()`.
        for a in 1:3
            lo[a]=bounds[a];hi[a]=bounds[3+a]
        end
    else
        for i in axes(mesh.coords,2),a in 1:3
            v=mesh.coords[a,i]
            v<lo[a] && (lo[a]=v);v>hi[a] && (hi[a]=v)
        end
    end
    # Upstream `mathEvaluator` binds variable `t`; the Tessella evaluator
    # binds `x`, so `t` rewrites to `x` (`u`/`v` stay illegal like upstream).
    progs=NTuple{3,Any}(SizeField.parse_matheval(
        _plugin_tvar_rewrite(strings[("MappingX","MappingY","MappingZ")[i]]))
        for i in 1:3)
    cuts=(Vector{Float64}(undef,nsx+1),Vector{Float64}(undef,nsy+1),
          Vector{Float64}(undef,nsz+1))
    counts=(nsx,nsy,nsz)
    for axis in 1:3
        pp=cuts[axis];n=counts[axis]
        for p in 0:n
            t=p/n
            r=_plugin_uv_eval(progs[axis],t,0.0)
            r!==nothing && (t=r)
            pp[p+1]=lo[axis]+t*(hi[axis]-lo[axis])
        end
    end
    emptyaxes=ntuple(a->cuts[a][1]==cuts[a][end],3)
    tags=Int64[];parts=Int64[]
    public=_mesh_public_tags(mesh)
    for index in 1:total
        rec=mesh_element_record(mesh,index)
        center=(sum(mesh.coords[1,rec.node_tags])/length(rec.node_tags),
                sum(mesh.coords[2,rec.node_tags])/length(rec.node_tags),
                sum(mesh.coords[3,rec.node_tags])/length(rec.node_tags))
        part=0
        for kx in 0:nsx-1
            part!=0 && break
            for ky in 0:nsy-1
                part!=0 && break
                for kz in 0:nsz-1
                    part!=0 && break
                    inx=emptyaxes[1] || (kx==0 && cuts[1][1]==center[1]) ||
                        (cuts[1][kx+1]<center[1]<=cuts[1][kx+2])
                    iny=emptyaxes[2] || (ky==0 && cuts[2][1]==center[2]) ||
                        (cuts[2][ky+1]<center[2]<=cuts[2][ky+2])
                    inz=emptyaxes[3] || (kz==0 && cuts[3][1]==center[3]) ||
                        (cuts[3][kz+1]<center[3]<=cuts[3][kz+2])
                    if inx && iny && inz
                        part=kx*nsy*nsz+ky*nsz+kz+1
                    end
                end
            end
        end
        # Elements whose barycenter lands outside every slice stay
        # unpartitioned upstream — the explicit `_partition` path reproduces
        # that as the "All elements are not partitioned" failure, impossible
        # in practice for monotone mappings over the model bounds.
        part==0 && continue
        etag=public===nothing ? Int64(index) :
            Int64(get(public.element_tags,index,index))
        push!(tags,etag);push!(parts,part)
    end
    _partition(nsx*nsy*nsz,tags,parts)
    return nothing
end

# `MappingX/Y/Z` use `t` — the evaluator binds x; rewrite `t`→`x` only
# (u/v stay illegal as upstream).
function _plugin_tvar_rewrite(expr::AbstractString)
    out=IOBuffer();i=firstindex(expr);n=lastindex(expr)
    while i<=n
        c=expr[i]
        if _plugin_ident_char(c)
            j=i
            while j<=n && _plugin_ident_char(expr[j]);j=nextind(expr,j);end
            ident=expr[i:prevind(expr,j)]
            write(out,ident=="t" ? "x" : ident)
            i=j;continue
        end
        write(out,c);i=nextind(expr,i)
    end
    return String(take!(out))
end

# ---- GaussPoints ---------------------------------------------------------------------
#
# Upstream `Plugin(GaussPoints)` walks all `Dimension` entities (optionally
# restricted to `PhysicalGroup`) and emits one SP record per Gauss point of
# `getIntegrationPoints(Order)` — the element's real coordinates, value set
# to the public element number. Tessella uses the same Gauss tables through
# `mesh_integration_points(element_type, "GaussN")`.

function _plugin_exec_gausspoints(spec::_PluginSpec,caller::AbstractString)
    numbers=spec.numbers
    order=_mesh_query_integer(Int(numbers["Order"]),caller,"Order")
    dim=_mesh_query_integer(Int(numbers["Dimension"]),caller,"Dimension")
    physical=_mesh_query_integer(Int(numbers["PhysicalGroup"]),
                                 caller,"PhysicalGroup")
    _slot_capture_locked!()
    index=_current_slot_index_locked()
    index==0 && throw(ArgumentError("$caller: no current model"))
    slot=MODEL_SLOTS[index]
    mesh=slot.mesh
    mesh===nothing && throw(ArgumentError("$caller: No entities"))
    class=slot.class===nothing ? _cached_classification_locked(mesh) :
        slot.class
    public=_mesh_public_tags(mesh)
    _,_,total=_mesh_element_offsets(mesh)
    total>0 || throw(ArgumentError("$caller: No entities"))
    flat=Float64[]
    for index in 1:total
        rec=mesh_element_record(mesh,index)
        rec.dimension==dim || continue
        # Upstream filters to `groups[physical]` whenever `physical` is set.
        if physical!=0
            entity=_plugin_element_entity(mesh,class,index)
            entity===nothing && continue
            phys=get(class.physicals,entity,Int32[])
            isempty(phys) && continue
            Int32(physical) in phys || continue
        end
        etag=public===nothing ? UInt64(index) :
            get(public.element_tags,index,UInt64(index))
        spec_rec=msh_spec(Int(rec.element_type))
        spec_rec.family===:trih && continue
        gauss=try
            mesh_integration_points(Int(rec.element_type),"Gauss$order")
        catch err
            err isa InterruptException && rethrow()
            continue
        end
        # `pnt` points arrive flattened `[u1 v1 w1 u2 ...]`; evaluate the
        # element's own nodal basis there (upstream `e->pnt(u,v,w,p)`).
        pnt=gauss[1]
        npts=length(pnt)÷3
        nn=length(rec.node_tags)
        for p in 1:npts
            u=pnt[3p-2];v=pnt[3p-1];w=pnt[3p]
            _,N,_=mesh_basis_functions(Int(rec.element_type),[u,v,w],
                                       "Lagrange")
            x=y=z=0.0
            for node in 1:min(nn,length(N))
                x+=N[node]*mesh.coords[1,Int(rec.node_tags[node])]
                y+=N[node]*mesh.coords[2,Int(rec.node_tags[node])]
                z+=N[node]*mesh.coords[3,Int(rec.node_tags[node])]
            end
            push!(flat,x,y,z,Float64(etag))
        end
    end
    if isempty(flat)
        _native_diagnostic_emit_locked(:error,"No entities")
        return nothing
    end
    _plugin_publish_list_view!("GaussPoints",["SP"],Dict("SP"=>length(flat)÷4),
                               Dict("SP"=>flat),Float64[];caller)
    return nothing
end

# ---- AnalyseMeshQuality -----------------------------------------------------------------
#
# Upstream `Plugin(AnalyseMeshQuality)` always computes per-element min/max
# Jacobian determinants over `DimensionOfElements` (-1 → model dimension, 4 →
# dims 2 and 3, else that dimension), plus minIGE/minICN when flagged; prints
# the (min, avg, max) statistics, creates `ElementData` views per dimension
# when `CreateView`, and hides elements whose measure crosses
# `HidingThreshold` (ICN preferred over IGE over J). Tessella delegates the
# measures to `mesh_element_qualities` (`minDetJac`/`maxDetJac`/`minSIGE`/
# `minSICN`) and emits the same views as ElementData on the owning model.

function _plugin_publish_element_view!(name,tags,values,slot,mesh;caller)
    view=_SessionView(name)
    view.model_kind="ElementData"
    view.model_step_ranges=[_NativeViewStepRange(0,0,slot,1)]
    counts=fill(Int32(1),length(tags))
    flat=Float64.(values)
    extent=_native_model_view_extent(tags,
        _native_view_model_count(slot,"ElementData"))
    view.model_data[(0,"ElementData")]=(UInt64.(tags),flat,1)
    view.model_contexts[(0,"ElementData")]=(
        model_slot=slot,time=0.0,counts=counts,extent=extent,
        partitions=Set{Int}([0]),
        interpolation_msh=_plugin_mesh_msh_types(mesh))
    view.numbers["NbTimeStep"]=1.0
    tag=_fresh_view_tag()
    SESSION_VIEWS[tag]=view
    VIEW_TAG_MAX[]=max(VIEW_TAG_MAX[],tag)
    return tag
end

function _plugin_exec_analyzemeshquality(spec::_PluginSpec,
                                         caller::AbstractString)
    numbers=spec.numbers
    compute_jac=!iszero(numbers["JacobianDeterminant"])
    compute_ige=!iszero(numbers["IGEMeasure"])
    compute_icn=!iszero(numbers["ICNMeasure"])
    threshold=numbers["HidingThreshold"]
    threshold_greater=!iszero(numbers["ThresholdGreater"])
    create_view=!iszero(numbers["CreateView"])
    asked=_mesh_query_integer(Int(numbers["DimensionOfElements"]),
                              caller,"DimensionOfElements")
    _slot_capture_locked!()
    index=_current_slot_index_locked()
    index==0 && throw(ArgumentError("$caller: no current model"))
    slot=MODEL_SLOTS[index]
    mesh=slot.mesh
    mesh===nothing && throw(ArgumentError(
        "$caller: Plugin(AnalyseMeshQuality) needs a mesh"))
    _,_,total=_mesh_element_offsets(mesh)
    total>0 || throw(ArgumentError(
        "$caller: Plugin(AnalyseMeshQuality) needs a mesh"))
    model_dim=Int(rec_dimension_max(mesh))
    (asked<0 || asked>4) && (asked=model_dim)
    dims=asked==4 ? [2,3] : [asked]
    filter!(d->d<=model_dim,dims)
    public=_mesh_public_tags(mesh)
    # Per-dim tag vectors preserving mesh order.
    per_dim=Dict{Int,Vector{UInt64}}(d=>UInt64[] for d in dims)
    for index in 1:total
        rec=mesh_element_record(mesh,index)
        rec.dimension in dims || continue
        etag=public===nothing ? UInt64(index) :
            get(public.element_tags,index,UInt64(index))
        push!(per_dim[rec.dimension],etag)
    end
    allmin=Float64[];allmax=Float64[];allige=Float64[];allicn=Float64[]
    allelem=UInt64[]
    for d in dims
        tags=per_dim[d]
        isempty(tags) && continue
        # Upstream `jacobianBasedQuality::minMaxJacobianDeterminant` on a
        # straight segment gives the constant det J = length/2 (reference
        # length 2); a zero-length segment collapses to 0.
        if d==1
            lens=mesh_element_qualities(mesh,tags,"minEdge")
            minj=lens./2;maxj=copy(minj)
            ige=compute_ige ?
                mesh_element_qualities(mesh,tags,"minSIGE") : Float64[]
            icn=compute_icn ?
                mesh_element_qualities(mesh,tags,"minSICN") : Float64[]
        else
            minj=mesh_element_qualities(mesh,tags,"minDetJac")
            maxj=mesh_element_qualities(mesh,tags,"maxDetJac")
            ige=compute_ige ?
                mesh_element_qualities(mesh,tags,"minSIGE") : Float64[]
            icn=compute_icn ?
                mesh_element_qualities(mesh,tags,"minSICN") : Float64[]
        end
        if create_view
            # The minJ/maxJ view exists only when `JacobianDeterminant` set;
            # the J measures themselves are always computed upstream.
            if compute_jac
                ratio=[maxj[i]>0 ? minj[i]/maxj[i] :
                       maxj[i]<0 ? maxj[i]/minj[i] : 0.0
                       for i in eachindex(tags)]
                _plugin_publish_element_view!("minJ/maxJ $(d)D",tags,ratio,
                                              slot,mesh;caller)
            end
            compute_ige && _plugin_publish_element_view!("IGE $(d)D",tags,
                                                         ige,slot,mesh;caller)
            compute_icn && _plugin_publish_element_view!("ICN $(d)D",tags,
                                                         icn,slot,mesh;caller)
        end
        append!(allelem,tags);append!(allmin,minj);append!(allmax,maxj)
        compute_ige && append!(allige,ige)
        compute_icn && append!(allicn,icn)
    end
    # Upstream `_printStat*` formats: `%8.3g` triples.
    _fmt(x)=string(round(x,sigdigits=4))
    function _print3(name,vec,tail)
        isempty(vec) && return
        _native_diagnostic_emit_locked(:info,
            "$name = $(_fmt(minimum(vec))), $(_fmt(sum(vec)/length(vec))), " *
            "$(_fmt(maximum(vec))) $tail")
    end
    if !isempty(allelem)
        _print3("minJ     ",allmin,"(min, avg, max)")
        ratio=[allmax[i]>0 ? allmin[i]/allmax[i] :
               allmax[i]<0 ? allmax[i]/allmin[i] : 0.0
               for i in eachindex(allmin)]
        nconst=count(<(1-1e-5),ratio)
        if 0<nconst<length(ratio)
            avg=sum(ratio[i] for i in eachindex(ratio) if ratio[i]<1-1e-5)/
                nconst
            _native_diagnostic_emit_locked(:info,
                "minJ/maxJ =           $(_fmt(avg))           (avg on the " *
                "$nconst non-constant elements)")
        end
        _print3("minJ/maxJ",ratio,"(worst, avg, best)")
        compute_ige && _print3("IGE      ",allige,"(worst, avg, best)")
        compute_icn && _print3("ICN      ",allicn,"(worst, avg, best)")
    else
        _native_diagnostic_emit_locked(:info,"No stat to print")
    end
    # `whichMeasure = ICN ? 2 : IGE ? 1 : Jac ? 0 : -1` — hide (and re-show
    # passing) elements of the scanned dims.
    which=compute_icn ? 2 : compute_ige ? 1 : compute_jac ? 0 : -1
    if threshold<99 && which>=0
        measure=which==2 ? allicn : which==1 ? allige :
            Float64[allmax[i]>0 ? allmin[i]/allmax[i] :
                    allmax[i]<0 ? allmax[i]/allmin[i] : 0.0
                    for i in eachindex(allmin)]
        for (i,etag) in enumerate(allelem)
            i>length(measure) && break
            hide=threshold_greater ? measure[i]>threshold :
                measure[i]<threshold
            ELEMENT_VISIBILITY[][etag]=Int32(hide ? 0 : 1)
        end
    end
    return nothing
end

# `rec_dimension_max(mesh)` — the model's highest scanned dimension.
function rec_dimension_max(mesh)
    total=last(_mesh_element_offsets(mesh))
    dmax=0
    for index in 1:total
        dmax=max(dmax,Int(mesh_element_record(mesh,index).dimension))
    end
    return dmax
end

# ---- MeshSizeFieldView --------------------------------------------------------------------
#
# Upstream `Plugin(MeshSizeFieldView)` evaluates mesh-size field
# `MeshSizeField` at every node of `View` and rewrites every component in
# place (`Field::putOnView`), or publishes a fresh list view when `View < 0`
# and no model view exists. Tessella resolves the field through the session
# size-field machinery, restricted to evaluable contexts.

# Build a session field spec into an evaluable `AbstractField`, mirroring
# `_session_size_field_locked` but rooted at an arbitrary tag. Only the
# natively-built catalog kinds plus `PostView` resolve here; the remaining
# context-backed kinds (`BoundaryLayer`, `Extend`, `AttractorAnisoCurve`,
# `AutomaticMeshSizeField`) fail like an unavailable context.
function _plugin_session_field_build(slot,fieldnum::Int,
                                     caller::AbstractString)
    spec0=get(FIELD_SPECS[],fieldnum,nothing)
    spec0===nothing && throw(ArgumentError(
        "$caller: unknown field tag $fieldnum"))
    specs=Dict{Int,GeoFieldSpec}(
        t=>GeoFieldSpec(t,f.kind,Dict(f.options),copy(f.order))
        for (t,f) in FIELD_SPECS[])
    params=GeoParams(OPTIONS["Mesh.MeshSizeMin"],OPTIONS["Mesh.MeshSizeMax"],
        OPTIONS["Mesh.MeshSizeFactor"],0,Dict{Tuple{Int,Int},String}(),
        specs,BACKGROUND_FIELD[],copy(BOUNDARY_LAYER_FIELDS[]))
    bounds=try
        model_bounding_box(slot.model,-1,-1)
    catch err
        err isa InterruptException && rethrow()
        nothing
    end
    bbox=bounds===nothing ? nothing :
        (bounds[1],bounds[4],bounds[2],bounds[5],bounds[3],bounds[6])
    entities=(d,name)->_field_entity_mesh(slot.model,d,name)
    function context_fields(spec,config,entities_,_params)
        if lowercase(spec.kind)=="postview"
            function selector(name,default)
                haskey(spec.options,name) || return default
                numeric=tryparse(Float64,strip(spec.options[name]))
                numeric!==nothing && isfinite(numeric) && isinteger(numeric) ||
                    throw(ArgumentError(
                        "$caller: Field[$(spec.tag)].$name must be a " *
                        "finite integer"))
                return Int(numeric)
            end
            view_tag=selector("ViewTag",-1)
            selector("ViewIndex",0);selector("IView",0)
            if view_tag<0
                tags=sort!(collect(keys(SESSION_VIEWS)))
                0<=config.view_index<length(tags) || throw(ArgumentError(
                    "$caller: Field[$(spec.tag)] references unknown view " *
                    "index $(config.view_index)"))
                view_tag=tags[config.view_index+1]
            end
            ref=get(SESSION_VIEWS,view_tag,nothing)
            ref===nothing && throw(ArgumentError(
                "$caller: Field[$(spec.tag)] references unknown view " *
                "$view_tag"))
            return _post_view_field(ref,spec.tag;
                crop_negative=config.crop_negative_values,
                use_closest=config.use_closest)
        end
        throw(ArgumentError(
            "$caller: Field[$(spec.tag)] kind $(spec.kind) requires model " *
            "context unavailable to plugin(MeshSizeFieldView)"))
    end
    return SizeField._build_geo_field(params,entities,fieldnum;
        model_bbox=bbox,context_fields=context_fields)
end

function _plugin_exec_meshsizefieldview(spec::_PluginSpec,
                                        caller::AbstractString)
    numbers=spec.numbers
    fieldnum=_mesh_query_integer(Int(numbers["MeshSizeField"]),
                                 caller,"MeshSizeField")
    _mesh_query_integer(Int(numbers["Component"]),caller,"Component")
    tag,view=_plugin_view(numbers,caller)
    _slot_capture_locked!()
    index=_current_slot_index_locked()
    index==0 && throw(ArgumentError("$caller: no current model"))
    field=_plugin_session_field_build(MODEL_SLOTS[index],fieldnum,caller)
    # `Field::putOnView` iterates step-0 elements only and its inner
    # `for(comp...)` shadows the `Component` argument: every component of
    # step 0 is rewritten, and the view name becomes "Field <id>".
    if view.model_kind===nothing
        steps=_plugin_list_steps(view,caller)
        for type in view.list_types
            layout=get(Post._POS_RECORD_TYPES,type,nothing)
            layout===nothing && continue
            kind,arity,nc=layout
            data=view.list_data[type];count=view.list_counts[type]
            width=3arity+arity*nc*steps
            for element in 0:count-1
                base=element*width
                for node in 1:arity
                    x=data[base+node];y=data[base+arity+node]
                    z=data[base+2arity+node]
                    v=Float64(SizeField.field_value(field,x,y,z))
                    for c in 1:nc
                        data[base+3arity+(node-1)*nc+c]=v
                    end
                end
            end
        end
    else
        # `Field::putOnView` iterates the model's mesh elements and writes
        # `getNumNodes()` (primary-count) slots: NodeData rows get the field
        # at each primary node, ElementData rows the value at the last
        # written primary node, ElementNodeData rows one value per primary
        # node slot — all components rewritten.
        for key in _plugin_model_keys(view)
            step,kind=key
            step==0 || continue
            ctx=view.model_contexts[key]
            mesh=_plugin_model_mesh(ctx.model_slot,caller)
            public=_mesh_public_tags(mesh)
            entity_tags,values,nc=view.model_data[key]
            counts=get(ctx,:counts,Int32[])
            offsets=Vector{Int}(undef,length(entity_tags)+1);offsets[1]=1
            for row in eachindex(entity_tags)
                width=row<=length(counts) ? Int(counts[row]) :
                    length(values)÷length(entity_tags)
                offsets[row+1]=offsets[row]+width
            end
            row_of=Dict{UInt64,Int}(t=>r for (r,t) in enumerate(entity_tags))
            _,_,etotal=_mesh_element_offsets(mesh)
            for index in 1:etotal
                rec=mesh_element_record(mesh,index)
                spec_rec=msh_spec(Int(rec.element_type))
                np=min(length(rec.node_tags),
                       Elements._msh_primary_nodes(spec_rec.family))
                if kind=="NodeData"
                    for nod in 1:np
                        ntag=_cache_node_tag(mesh,Int(rec.node_tags[nod]))
                        row=get(row_of,ntag,0)
                        row==0 && continue
                        p=Tuple(mesh.coords[:,Int(rec.node_tags[nod])])
                        v=Float64(SizeField.field_value(field,p...))
                        for c in 1:nc
                            values[offsets[row]+c-1]=v
                        end
                    end
                    continue
                end
                etag=public===nothing ? UInt64(index) :
                    get(public.element_tags,index,UInt64(index))
                row=get(row_of,etag,0)
                row==0 && continue
                if kind=="ElementData"
                    # `for comp` writes `data[comp]` per primary node — the
                    # last write wins: field at the last primary node.
                    p=Tuple(mesh.coords[:,Int(rec.node_tags[np])])
                    v=Float64(SizeField.field_value(field,p...))
                    for c in 1:nc
                        values[offsets[row]+c-1]=v
                    end
                else
                    # ElementNodeData: `nc*nod+comp` per primary node.
                    width=offsets[row+1]-offsets[row]
                    for nod in 1:np
                        (nod-1)*nc+nc<=width || break
                        p=Tuple(mesh.coords[:,Int(rec.node_tags[nod])])
                        v=Float64(SizeField.field_value(field,p...))
                        for c in 1:nc
                            values[offsets[row]+(nod-1)*nc+c-1]=v
                        end
                    end
                end
            end
        end
    end
    view.name="Field $fieldnum"
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
    "Gradient"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"View"=>-1.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_gradient),
    "CutPlane"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"A"=>1.0,"B"=>0.0,"C"=>0.0,"D"=>0.0,
        "ExtractVolume"=>0.0,"View"=>-1.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_cutplane),
    "CutSphere"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"Xc"=>0.0,"Yc"=>0.0,"Zc"=>0.0,"R"=>0.25,
        "ExtractVolume"=>0.0,"RecurLevel"=>3.0,"TargetError"=>1e-4,
        "View"=>-1.0,"Visible"=>1.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_cutsphere),
    "Isosurface"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"Value"=>0.0,"ExtractVolume"=>0.0,"RecurLevel"=>3.0,
        "TargetError"=>1e-4,"View"=>-1.0,"Visible"=>1.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_isosurface),
    "Warp"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"Factor"=>1.0,"TimeStep"=>0.0,"SmoothingAngle"=>180.0,
        "View"=>-1.0,"OtherView"=>-1.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_warp),
    "MathEval"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"TimeStep"=>-1.0,"View"=>-1.0,"OtherTimeStep"=>-1.0,
        "OtherView"=>-1.0,"ForceInterpolation"=>0.0,"PhysicalGroup"=>-1.0,
        "Dimension"=>-1.0),
        string_options=Dict{String,String}(
            "Expression0"=>"Sqrt(v0^2+v1^2+v2^2)",
            "Expression1"=>"","Expression2"=>"","Expression3"=>"",
            "Expression4"=>"","Expression5"=>"","Expression6"=>"",
            "Expression7"=>"","Expression8"=>""),
        run_fn=_plugin_exec_matheval),
    "Scal2Vec"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"ViewX"=>-1.0,"ViewY"=>-1.0,"ViewZ"=>-1.0),
        string_options=Dict{String,String}("NameNewView"=>"NewView"),
        run_fn=_plugin_exec_scal2vec),
    "Divergence"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"View"=>-1.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_divergence),
    "Curl"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"View"=>-1.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_curl),
    "Lambda2"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"Eigenvalue"=>2.0,"View"=>-1.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_lambda2),
    "SphericalRaise"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"Xc"=>0.0,"Yc"=>0.0,"Zc"=>0.0,"Raise"=>1.0,"Offset"=>0.0,
        "TimeStep"=>0.0,"View"=>-1.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_spherical_raise),
    "Annotate"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"X"=>50.0,"Y"=>30.0,"Z"=>0.0,"ThreeD"=>0.0,
        "FontSize"=>14.0,"View"=>-1.0),
        string_options=Dict{String,String}("Text"=>"My Text",
            "Font"=>"Helvetica","Align"=>"Left"),
        run_fn=_plugin_exec_annotate),
    "Eigenvalues"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"View"=>-1.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_eigenvalues),
    "ModulusPhase"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"RealPart"=>0.0,"ImaginaryPart"=>1.0,"View"=>-1.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_modulus_phase),
    "HarmonicToTime"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"RealPart"=>0.0,"ImaginaryPart"=>1.0,"NumSteps"=>20.0,
        "TimeSign"=>-1.0,"Frequency"=>1.0,"NumPeriods"=>1.0,"View"=>-1.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_harmonic_to_time),
    "Density"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"View"=>-1.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_density),
    "MakeSimplex"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"View"=>-1.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_makesimplex),
    # Upstream registers `ModifyComponents` (plural); there is no singular
    # `ModifyComponent` in Gmsh 4.15.2.
    "ModifyComponents"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"TimeStep"=>-1.0,"View"=>-1.0,"OtherTimeStep"=>-1.0,
        "OtherView"=>-1.0,"ForceInterpolation"=>0.0),
        string_options=Dict{String,String}(
            "Expression0"=>"v0 * Sin(x)",
            "Expression1"=>"","Expression2"=>"","Expression3"=>"",
            "Expression4"=>"","Expression5"=>"","Expression6"=>"",
            "Expression7"=>"","Expression8"=>""),
        run_fn=_plugin_exec_modifycomponents),
    "FieldFromAmplitudePhase"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"Wavenumber"=>5.0,"AmplitudeView"=>0.0,"PhaseView"=>1.0),
        string_options=Dict{String,String}("MeshFile"=>"fine.msh"),
        run_fn=_plugin_exec_field_from_amplitude_phase),
    "ExtractEdges"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"Angle"=>40.0,"IncludeBoundary"=>1.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_extractedges),
    # Tessella extension: Gmsh 4.15.2 registers no `Vorticity` plugin.
    "Vorticity"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"View"=>-1.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_vorticity),
    "CutBox"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"X0"=>0.0,"Y0"=>0.0,"Z0"=>0.0,
        "X1"=>1.0,"Y1"=>0.0,"Z1"=>0.0,
        "X2"=>0.0,"Y2"=>1.0,"Z2"=>0.0,
        "X3"=>0.0,"Y3"=>0.0,"Z3"=>1.0,
        "NumPointsU"=>20.0,"NumPointsV"=>20.0,"NumPointsW"=>20.0,
        "ConnectPoints"=>1.0,"Boundary"=>1.0,"View"=>-1.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_cutbox),
    "CutGrid"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"X0"=>0.0,"Y0"=>0.0,"Z0"=>0.0,
        "X1"=>1.0,"Y1"=>0.0,"Z1"=>0.0,
        "X2"=>0.0,"Y2"=>1.0,"Z2"=>0.0,
        "NumPointsU"=>20.0,"NumPointsV"=>20.0,"ConnectPoints"=>1.0,
        "View"=>-1.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_cutgrid),
    "CutParametric"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"MinU"=>0.0,"MaxU"=>2*3.1416,"NumPointsU"=>180.0,
        "MinV"=>0.0,"MaxV"=>2*3.1416,"NumPointsV"=>180.0,
        "ConnectPoints"=>0.0,"View"=>-1.0),
        string_options=Dict{String,String}(
            "X"=>"2 * Cos(u) * Sin(v)","Y"=>"4 * Sin(u) * Sin(v)",
            "Z"=>"0.1 + 0.5 * Cos(v)"),
        run_fn=_plugin_exec_cutparametric),
    "StreamLines"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"X0"=>0.0,"Y0"=>0.0,"Z0"=>0.0,
        "X1"=>1.0,"Y1"=>0.0,"Z1"=>0.0,
        "X2"=>0.0,"Y2"=>1.0,"Z2"=>0.0,
        "NumPointsU"=>10.0,"NumPointsV"=>1.0,"DT"=>0.1,"MaxIter"=>100.0,
        "TimeStep"=>0.0,"View"=>-1.0,"OtherView"=>-1.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_streamlines),
    "Particles"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"X0"=>0.0,"Y0"=>0.0,"Z0"=>0.0,
        "X1"=>1.0,"Y1"=>0.0,"Z1"=>0.0,
        "X2"=>0.0,"Y2"=>1.0,"Z2"=>0.0,
        "NumPointsU"=>10.0,"NumPointsV"=>1.0,"A2"=>1.0,"A1"=>0.0,"A0"=>0.0,
        "DT"=>0.1,"MaxIter"=>100.0,"TimeStep"=>0.0,"View"=>-1.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_particles),
    "Triangulate"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"View"=>-1.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_triangulate),
    "Tetrahedralize"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"View"=>-1.0),
        string_options=Dict{String,String}(),
        run_fn=_plugin_exec_tetrahedralize),
    "NearestNeighbor"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"View"=>-1.0),
        string_options=Dict{String,String}(),
        run_fn=_plugin_exec_nearestneighbor),
    "Scal2Tens"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"NumberOfComponents"=>9.0,
        "View0"=>-1.0,"View1"=>-1.0,"View2"=>-1.0,"View3"=>-1.0,
        "View4"=>-1.0,"View5"=>-1.0,"View6"=>-1.0,"View7"=>-1.0,
        "View8"=>-1.0),
        string_options=Dict{String,String}("NameNewView"=>"NewView"),
        run_fn=_plugin_exec_scal2tens),
    "Eigenvectors"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"ScaleByEigenvalues"=>1.0,"View"=>-1.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_eigenvectors),
    "Summation"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"View 0"=>-1.0,"View 1"=>-1.0,"View 2"=>-1.0,
        "View 3"=>-1.0,"View 4"=>-1.0,"View 5"=>-1.0,"View 6"=>-1.0,
        "View 7"=>-1.0),
        string_options=Dict{String,String}("Resuling View Name"=>"default"),
        run_fn=_plugin_exec_summation),
    "MeshVolume"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"PhysicalGroup"=>-1.0,"Dimension"=>3.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_meshvolume),
    "LongitudeLatitude"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"View"=>-1.0),
        string_options=Dict{String,String}(),
        run_fn=_plugin_exec_longitude_latitude),
    "ShowNeighborElements"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"NumLayers"=>1.0,"Element1"=>0.0,"Element2"=>0.0,
        "Element3"=>0.0,"Element4"=>0.0,"Element5"=>0.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_showneighbor),
    "GaussPoints"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"Order"=>0.0,"Dimension"=>2.0,"PhysicalGroup"=>0.0),
        string_options=Dict{String,String}(),run_fn=_plugin_exec_gausspoints),
    "AnalyseMeshQuality"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"JacobianDeterminant"=>0.0,"IGEMeasure"=>0.0,
        "ICNMeasure"=>0.0,"HidingThreshold"=>99.0,"ThresholdGreater"=>1.0,
        "CreateView"=>0.0,"Recompute"=>0.0,"DimensionOfElements"=>-1.0),
        string_options=Dict{String,String}(),
        run_fn=_plugin_exec_analyzemeshquality),
    "MeshSizeFieldView"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"MeshSizeField"=>0.0,"View"=>-1.0,"Component"=>0.0),
        string_options=Dict{String,String}(),
        run_fn=_plugin_exec_meshsizefieldview),
    "SimplePartition"=>(number_options=Dict{String,Float64}(
        "Run"=>0.0,"NumSlicesX"=>4.0,"NumSlicesY"=>1.0,"NumSlicesZ"=>1.0),
        string_options=Dict{String,String}(
            "MappingX"=>"t","MappingY"=>"t","MappingZ"=>"t"),
        run_fn=_plugin_simple_partition_exec),
)
