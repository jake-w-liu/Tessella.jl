# Generation owns a temporary coarse Curve view. Original lower products and
# their identities remain detached from the temporary Line constructors.
struct _NativeFullQuadMiddle
    tag::UInt64
    owner::Tuple{Int,Int32}
    parameter::Union{Nothing,Vector{Float64}}
    world::NTuple{3,Float64}
    chord_world::NTuple{3,Float64}
    # `born` marks a middle minted for an uncoarsened Curve cell edge: it is
    # a fresh identity absent from the lower producer, not a retained
    # halved vertex restoring onto the coarse boundary.
    born::Bool
end
_NativeFullQuadMiddle(tag,owner,parameter,world)=
    _NativeFullQuadMiddle(tag,owner,parameter,world,world,false)
_NativeFullQuadMiddle(tag,owner,parameter,world,chord)=
    _NativeFullQuadMiddle(tag,owner,parameter,world,chord,false)

# A supplied lower Point can split a Curve cell. Transport that actual UID
# into the detached Curve product before pairing cells: otherwise the Point
# splits a temporary super-edge after its stored middle has been detached,
# and that middle can no longer restore. No vertex identity is inferred from
# an output coordinate; the source contains only the consumed lower products.
function _gmsh_fullquad_transport_lower!(m,tag,authority,high_water,options,caller)
    changed=Set{Int32}()
    (_surface_type(m,tag)===:plane && !isempty(get(m.embeds,(2,tag),NTuple{2,Int}[]))) ||
        return changed,high_water
    frame=_gmsh_native_surface_frame(m,tag,options,caller)
    source=_gmsh_native_surface_source(m,tag,authority,frame,caller;ctx_lc=options.ctx_lc)
    _gmsh_native_surface_split_constraints!(source)
    data=authority.mesh.entity_data
    rows=Dict(uid=>Int32(row) for (row,uid) in enumerate(data.external_node_tags))
    for (curve,chain) in _gmsh_native_surface_curve_chains(source)
        owner=(1,curve);original=authority.products[owner]
        labels=Dict{Tuple{UInt64,UInt64},UInt64}()
        members=Set{UInt64}()
        for (block,tags) in zip(original.blocks,original.element_tags),column in axes(block.nodes,2)
            a,b=original.rows[block.nodes[:,column]]
            u,v=data.external_node_tags[a],data.external_node_tags[b]
            labels[(u,v)]=tags[column];push!(members,u);push!(members,v)
        end
        all(in(members),chain) && continue
        global_rows=Int32[rows[uid] for uid in unique(chain)]
        local_rows=Dict(data.external_node_tags[row]=>Int32(i) for (i,row) in enumerate(global_rows))
        connectivity=Matrix{Int32}(undef,2,length(chain)-1)
        tags=UInt64[]
        for i in 1:length(chain)-1
            a,b=chain[i],chain[i+1]
            connectivity[:,i].=(local_rows[a],local_rows[b])
            label=get(labels,(a,b),UInt64(0))
            if label==0
                high_water<typemax(UInt64) || throw(ArgumentError(
                    "$caller: full-quad transported Curve cells exhaust UInt64 labels"))
                label=(high_water+=UInt64(1))
            end
            push!(tags,label)
        end
        authority.products[owner]=_ModelMeshProduct(owner,global_rows,
            ElementBlock[ElementBlock(1,connectivity)],Vector{UInt64}[tags])
        push!(changed,curve)
    end
    return changed,high_water
end

function _gmsh_fullquad_input(m,tag,input,element_high_water,options,caller)
    middles=Dict{Tuple{UInt64,UInt64},_NativeFullQuadMiddle}()
    compound=any(pair->first(pair)==2 && tag in last(pair),m.meshing.compounds)
    # Mesh-periodic slaves copy an already generated master. This Plane
    # producer has no CAD periodic parametrization or seam to reject.
    if compound
        return (authority=input,middles,element_high_water,admitted=false,
            ctx_lc=options.ctx_lc,chains=Dict{Int32,Vector{UInt64}}(),
            parameters=Dict{UInt64,Float64}(),node_high_water=UInt64(0))
    end
    authority=deepcopy(input)
    source=authority.mesh
    data=source.entity_data
    _,embedded=_model_surface_embedding_tags(m,tag,caller)
    curves=vcat(embedded,_gmsh_native_surface_curves(m,tag))
    high_water=UInt64(element_high_water)
    transported,high_water=_gmsh_fullquad_transport_lower!(m,tag,authority,high_water,options,caller)
    # The committed dim-1 write-back restores each consumed Curve's own
    # chain — the temporary paired Line constructors never leave this
    # producer. `chains` records that pre-coarse ordered node-UID chain and
    # `parameters` its stored Curve parameters so a genuinely grown
    # recovered chain can splice its new vertices back at the right span.
    chains=Dict{Int32,Vector{UInt64}}()
    parameters=Dict{UInt64,Float64}()
    records=NamedTuple[]
    proposed=Dict{Tuple{UInt64,UInt64},Int}()
    for curve in curves
        haskey(m.meshing.transfinite_curves,curve) && _model_mesh_warn(options,
            "Full-quad recombination only compatible with transfinite meshes if those are performed first")
        curve in m.meshing.degenerated && continue
        owner=(1,Int32(curve))
        original=get(authority.products,owner,nothing)
        original!==nothing || throw(ArgumentError(
            "$caller: full-quad input lacks actual Curve[$curve] identities"))
        cells=NTuple{2,Int32}[]
        for block in original.blocks
            block.msh==1 || throw(ArgumentError(
                "$caller: full-quad Curve[$curve] needs linear Line constructors"))
            append!(cells,Tuple(block.nodes[:,column]) for column in axes(block.nodes,2))
        end
        adjacency=Dict{Int32,Int32}();incoming=Set{Int32}();consistent=true
        for (a,b) in cells
            haskey(adjacency,a) && (consistent=false)
            adjacency[a]=b;push!(incoming,b)
        end
        heads=consistent ? Int32[a for (a,_) in cells if !(a in incoming)] : Int32[]
        (consistent && !isempty(cells) && length(heads)<=1) || (consistent=false)
        if consistent
            vertex=isempty(heads) ? first(cells)[1] : only(heads)
            start=vertex;chain=Int32[vertex]
            for _ in eachindex(cells)
                vertex=get(adjacency,vertex,Int32(0))
                vertex==0 && break
                push!(chain,vertex);vertex==start && break
            end
            length(chain)==length(cells)+1 &&
                (chains[Int32(curve)]=UInt64[data.external_node_tags[original.rows[row]]
                                             for row in chain])
        end
        for row in original.rows
            data.node_entities[row]==owner || continue
            parameter=data.node_parametric[row]
            (parameter===nothing || isempty(parameter)) && continue
            parameters[data.external_node_tags[row]]=parameter[1]
        end
        if isodd(length(cells)) && !(Int32(curve) in transported)
            _model_mesh_error(options,"1D mesh cannot be divided by 2")
        end
        # Keep supplied Point rows as coarse endpoints. A transported Point
        # can leave an unmatched literal cell; it receives a fresh Curve
        # middle at subdivision, just like an uncoarsened collision edge.
        groups=NTuple{3,Int32}[]
        index=1
        while index<=length(cells)
            first,last=cells[index]
            if index<length(cells) &&
                    data.node_entities[original.rows[last]][1]!=0
                next,finish=cells[index+1]
                last==next || throw(ArgumentError(
                    "$caller: Curve[$curve] consecutive Line constructors do not share their middle identity"))
                push!(groups,(first,last,finish));index+=2
            else
                push!(groups,(first,Int32(0),last));index+=1
            end
        end
        keys=Tuple{UInt64,UInt64}[]
        for (first,middle,last) in groups
            middle==0 && continue
            a,c=original.rows[first],original.rows[last]
            key=minmax(data.external_node_tags[a],data.external_node_tags[c])
            push!(keys,key)
            proposed[key]=get(proposed,key,0)+1
        end
        # Every literal cell edge is a coarse-mesh constrained edge when its
        # Curve is left uncoarsened; `middles` (and the triangulation edge
        # set itself) cannot hold two entries for one node pair, so a
        # proposed super-edge colliding with another proposed super-edge —
        # a four-cell closed Curve, or two Curves spanning the same two
        # vertices — or with another Curve's literal cell edge leaves its
        # Curve uncoarsened instead.
        edges=Set{Tuple{UInt64,UInt64}}(
            minmax(data.external_node_tags[original.rows[cell[1]]],
                   data.external_node_tags[original.rows[cell[2]]]) for cell in cells)
        push!(records,(curve=curve,cells=cells,groups=groups,keys=keys,edges=edges,original=original))
    end
    node_high_water=maximum(data.external_node_tags;init=UInt64(0))
    for record in records
        (;curve,cells,keys,edges,original)=record
        owner=(1,Int32(curve))
        keyset=Set(keys)
        collide=any(key->proposed[key]>1,keys) ||
            any(other->!isdisjoint(keyset,other.edges),
                other for other in records if other!==record)
        if collide && haskey(chains,Int32(curve))
            # The Curve keeps its literal cells in the coarse view. Each
            # cell edge still needs its subdivision middle — born as a fresh
            # Curve-owned vertex on the curve itself — and the committed
            # chain records the grown sequence so the dim-1 write-back keeps
            # pace with the bisected constraint edges.
            chain=chains[Int32(curve)]
            born=Dict{Tuple{UInt64,UInt64},UInt64}()
            for cell in cells
                a,b=original.rows[cell[1]],original.rows[cell[2]]
                key=minmax(data.external_node_tags[a],data.external_node_tags[b])
                wa=Tuple(source.coords[:,a]);wb=Tuple(source.coords[:,b])
                chord=ntuple(axis->0.5*(wa[axis]+wb[axis]),3)
                parameter=only(model_parametrization(m,1,Int(curve),collect(chord)))
                isfinite(parameter) || throw(ArgumentError(
                    "$caller: full-quad Curve[$curve] keeps a non-finite edge parameter"))
                world=_gmsh_native_curve_sample(m,Int(curve),parameter,caller)
                node_high_water<typemax(UInt64) || throw(ArgumentError(
                    "$caller: full-quad uncoarsened Curve middles exhaust UInt64 labels"))
                uid=(node_high_water+=UInt64(1))
                middles[key]=_NativeFullQuadMiddle(uid,owner,Float64[parameter],
                    world,chord,true)
                parameters[uid]=parameter
                born[key]=uid
            end
            grown=UInt64[]
            for index in 1:length(chain)-1
                push!(grown,chain[index])
                push!(grown,born[minmax(chain[index],chain[index+1])])
            end
            push!(grown,chain[end])
            chains[Int32(curve)]=grown
            continue
        end
        coarse=NTuple{2,Int32}[]
        labels=UInt64[]
        born=Dict{Tuple{UInt64,UInt64},UInt64}()
        for (first,middle,last) in record.groups
            if middle==0
                a,c=original.rows[first],original.rows[last]
                key=minmax(data.external_node_tags[a],data.external_node_tags[c])
                chord=ntuple(axis->0.5*(source.coords[axis,a]+source.coords[axis,c]),3)
                parameter=only(model_parametrization(m,1,Int(curve),collect(chord)))
                isfinite(parameter) || throw(ArgumentError(
                    "$caller: full-quad Curve[$curve] keeps a non-finite edge parameter"))
                node_high_water<typemax(UInt64) || throw(ArgumentError(
                    "$caller: full-quad literal Curve middles exhaust UInt64 labels"))
                uid=(node_high_water+=UInt64(1))
                middles[key]=_NativeFullQuadMiddle(uid,owner,Float64[parameter],
                    _gmsh_native_curve_sample(m,Int(curve),parameter,caller),chord,true)
                parameters[uid]=parameter;born[key]=uid
                push!(coarse,(first,last))
                high_water<typemax(UInt64) || throw(ArgumentError(
                    "$caller: full-quad temporary Line constructors exhaust UInt64 labels"))
                push!(labels,high_water+=UInt64(1))
                continue
            end
            a,b,c=original.rows[first],original.rows[middle],original.rows[last]
            point=ntuple(axis->0.5*(source.coords[axis,a]+source.coords[axis,c]),3)
            parameter=data.node_parametric[b]
            restored=if !m.meshing.second_order_linear && data.node_entities[b][1]==1 &&
                        parameter!==nothing && !isempty(parameter)
                _gmsh_native_curve_sample(m,Int(data.node_entities[b][2]),parameter[1],caller)
            else
                point
            end
            payload=_NativeFullQuadMiddle(data.external_node_tags[b],data.node_entities[b],
                parameter===nothing ? nothing : copy(parameter),restored,point)
            middles[minmax(data.external_node_tags[a],data.external_node_tags[c])]=payload
            source.coords[:,b].=point
            push!(coarse,(first,last))
            high_water<typemax(UInt64) || throw(ArgumentError(
                "$caller: full-quad temporary Line constructors exhaust UInt64 labels"))
            high_water+=UInt64(1)
            push!(labels,high_water)
        end
        if !isempty(born) && haskey(chains,Int32(curve))
            chain=chains[Int32(curve)];grown=UInt64[]
            for index in 1:length(chain)-1
                push!(grown,chain[index])
                middle=get(born,minmax(chain[index],chain[index+1]),UInt64(0))
                middle==0 || push!(grown,middle)
            end
            push!(grown,chain[end]);chains[Int32(curve)]=grown
        end
        connectivity=Matrix{Int32}(undef,2,length(coarse))
        for (column,cell) in enumerate(coarse)
            connectivity[:,column].=cell
        end
        authority.products[owner]=_ModelMeshProduct(owner,original.rows,
            ElementBlock[ElementBlock(1,connectivity)],Vector{UInt64}[labels])
    end
    return (authority,middles,element_high_water=high_water,admitted=!isempty(middles),
        ctx_lc=options.ctx_lc,chains,parameters,node_high_water)
end

# The halved Curve middles restore through `middles` — they are not coarse
# boundary vertices, so their stored parameters must not force a fresh
# source vertex inside a paired segment. Keep any other forced parameter
# (periodic synchronizations, size-at-parameter points) for the ordinary
# subdivision path to consume.
function _gmsh_fullquad_forced(forced,quad_input)
    consumed=Dict{Int,Vector{Float64}}()
    for payload in values(quad_input.middles)
        payload.owner[1]==1 || continue
        parameter=payload.parameter
        if parameter!==nothing && !isempty(parameter)
            push!(get!(()->Float64[],consumed,Int(payload.owner[2])),parameter[1])
        end
    end
    isempty(consumed) && return forced
    output=Dict{Int,Vector{Float64}}()
    for (curve,list) in forced
        drop=get(consumed,curve,nothing)
        if drop===nothing
            output[curve]=list;continue
        end
        kept=Float64[]
        for value in list
            any(d->isequal(value,d) ||
                abs(value-d)<=128eps(Float64)*max(1.0,abs(value),abs(d)),drop) &&
                continue
            push!(kept,value)
        end
        isempty(kept) || (output[curve]=kept)
    end
    return output
end

# The committed dim-1 product is the restored (pre-coarse) chain — the
# temporary paired Line constructors never leave the producer. A recovered
# chain that genuinely grew past its consumed cells (a forced periodic or
# size-at-parameter vertex minted inside a coarse segment) splices its new
# node UIDs back into the consumed chain at their Curve parameters. A curve
# whose pre-coarse chain could not be reconstructed stays untouched — the
# coarse view must never overwrite the committed product.
function _gmsh_fullquad_restored_chains(m,source,grown,quad_input,caller)
    restored=Dict{Int32,Vector{UInt64}}()
    for (curve,chain) in grown
        original=get(quad_input.chains,curve,nothing)
        original===nothing && continue
        members=Set{UInt64}(original)
        inserts=UInt64[uid for uid in chain if !(uid in members)]
        isempty(inserts) && (restored[curve]=original;continue)
        curve_id=Int(curve)
        t0,t1=_model_curve_param_bounds(m,curve_id,caller)
        a,b=m.curves[curve_id]
        closed=a==b
        forward=true
        if closed
            samples=Float64[]
            for index in 2:length(original)-1
                value=get(quad_input.parameters,original[index],nothing)
                value===nothing || push!(samples,value)
            end
            forward=length(samples)<2 || issorted(samples)
        end
        function parameter_at(uid,position)
            closed && position==1 && return forward ? t0 : t1
            closed && position==length(original) && return forward ? t1 : t0
            stored=get(quad_input.parameters,uid,nothing)
            stored===nothing || return stored
            row=get(source.positions,uid,0)
            if row!=0
                parameter=source.parameters[row]
                if parameter!==nothing && !isempty(parameter)
                    return parameter[1]
                end
                if source.owners[row][1]==0
                    row==get(source.point_positions,Int(a),0) && return t0
                    row==get(source.point_positions,Int(b),0) && return t1
                end
                return model_parametrization(m,1,curve_id,
                    collect(source.world[row]))[1]
            end
            throw(ArgumentError("$caller: full-quad restored Curve[$curve] " *
                "vertex $uid lacks its Curve parameter"))
        end
        merged=UInt64[]
        pending=copy(inserts)
        for index in 1:length(original)-1
            first,last=original[index],original[index+1]
            push!(merged,first)
            pa=parameter_at(first,index);pb=parameter_at(last,index+1)
            lo,hi=minmax(pa,pb)
            inside=UInt64[]
            for uid in pending
                value=parameter_at(uid,0)
                lo<value<hi && push!(inside,uid)
            end
            isempty(inside) && continue
            sort!(inside;by=uid->parameter_at(uid,0),rev=pa>pb)
            append!(merged,inside)
            filter!(uid->!(uid in inside),pending)
        end
        push!(merged,original[end])
        isempty(pending) || throw(ArgumentError(
            "$caller: full-quad restored Curve[$curve] loses a grown vertex"))
        restored[curve]=merged
    end
    return restored
end

function _gmsh_fullquad_size_options(options::_ModelMesh1DOptions)
    names=fieldnames(_ModelMesh1DOptions)
    return _ModelMesh1DOptions((name===:lc_factor ? 2.0*getfield(options,name) :
        getfield(options,name) for name in names)...)
end

function _gmsh_fullquad_append_node!(state,tag,owner,parameter,point,uv,
                                    parameters,retained,caller;keep=owner[1]<2)
    length(state.node_tags)<typemax(Int32) || throw(ArgumentError(
        "$caller: full-quad subdivision exhausts dense node indices"))
    row=length(state.node_tags)+1
    if row>size(state.world,2)
        state.world=hcat(state.world,collect(point))
        state.uv=hcat(state.uv,collect(uv))
    else
        state.world[:,row].=point
        state.uv[:,row].=uv
    end
    push!(state.node_tags,tag)
    push!(state.owners,owner)
    push!(state.alive,true)
    push!(parameters,parameter)
    push!(retained,keep)
    return Int32(length(state.node_tags))
end

# Owned-boundary middle UV resolution for the planar/OCC fronts — the
# reparamOnFace transport `_gmsh_native_surface_add_source!` applies to the
# same Curve/Point payloads: a Curve-owned middle maps through its stored
# Curve parameter (`_occ_edge_reparam_on_face` on an OCC frame, the curve
# sample projected onto the mean Plane otherwise), a Point-owned middle
# projects its restored world position. Parametric fronts supply their own
# `middle_uv_of` instead (e.g. the ruled patch's `_ruled_curve_parameters`).
function _gmsh_fullquad_boundary_middle_uv(m,tag,frame,payload,point,ctx_lc,caller)
    if payload.owner[1]==1
        value=payload.parameter===nothing || isempty(payload.parameter) ?
            0.0 : payload.parameter[1]
        return frame.occ ? _occ_edge_reparam_on_face(m,m.surface_geometry[tag],
            Int(payload.owner[2]),value,1,ctx_lc,caller) :
            _gmsh_native_surface_parameters(frame,
                _gmsh_native_curve_sample(m,Int(payload.owner[2]),value,caller))
    end
    return _gmsh_native_surface_parameters(frame,point)
end

function _gmsh_fullquad_subdivide!(state,m,tag,middles,parameters,retained,
                                  node_high_water,frame,caller;ctx_lc=0.0,
                                  middle_uv_of=nothing,point_lift=nothing)
    high_water=Ref(UInt64(node_high_water))
    owner=(2,Int32(tag))
    edge_rows=Dict{Tuple{UInt64,UInt64},Int32}()
    positions=Dict(uid=>Int32(row) for (row,uid) in enumerate(state.node_tags))
    quads=NTuple{4,Int32}[]
    labels=UInt64[]
    # Reserve exactly the potentially new shared edge and cell-center rows.
    # Appending each coordinate column separately would copy O(N^2) storage.
    edges=Set{Tuple{UInt64,UInt64}}()
    for cells in (state.triangles,state.quads),cell in cells,index in eachindex(cell)
        push!(edges,minmax(state.node_tags[cell[index]],
            state.node_tags[cell[mod1(index+1,length(cell))]]))
    end
    known=Set(state.node_tags);new_rows=length(state.triangles)+length(state.quads)
    for key in edges
        payload=get(middles,key,nothing)
        if payload===nothing
            new_rows+=1
        elseif !(payload.tag in known)
            push!(known,payload.tag);new_rows+=1
        end
    end
    count=length(state.node_tags)+new_rows
    count<=typemax(Int32) || throw(ArgumentError(
        "$caller: full-quad subdivision exhausts dense node indices"))
    previous=length(state.node_tags)
    world=Matrix{Float64}(undef,3,count);world[:,1:previous].=state.world
    uv=Matrix{Float64}(undef,2,count);uv[:,1:previous].=state.uv
    state.world=world;state.uv=uv
    function birth()
        high_water[]<typemax(UInt64) || throw(ArgumentError(
            "$caller: full-quad subdivision exhausts UInt64 node labels"))
        high_water[]+=UInt64(1)
        return high_water[]
    end
    function edge_middle(a,b)
        key=minmax(state.node_tags[a],state.node_tags[b])
        row=get(edge_rows,key,Int32(0))
        row!=0 && return row
        payload=get(middles,key,nothing)
        middle_uv=ntuple(axis->0.5*(state.uv[axis,a]+state.uv[axis,b]),2)
        if payload!==nothing && payload.owner[1] in (0,1)
            middle_uv=middle_uv_of===nothing ?
                _gmsh_fullquad_boundary_middle_uv(m,tag,frame,payload,
                    payload.world,ctx_lc,caller) :
                middle_uv_of(payload,payload.world)
        end
        # Interior edge middles and parametric middles store `middle_uv` on
        # the face; the world row is the chord midpoint on a planar frame and
        # the on-surface `point_lift` on a curved patch (upstream's
        # MFaceVertex at the averaged parameters).
        point=payload===nothing ?
            (point_lift===nothing ?
                ntuple(axis->0.5*(state.world[axis,a]+state.world[axis,b]),3) :
                point_lift(middle_uv...)) :
            payload.world
        uid=payload===nothing ? birth() : payload.tag
        row=get(positions,uid,Int32(0))
        if row==0
            row=_gmsh_fullquad_append_node!(state,uid,
                payload===nothing ? owner : payload.owner,
                payload===nothing ? collect(middle_uv) : deepcopy(payload.parameter),point,middle_uv,
                parameters,retained,caller;
                keep=payload===nothing ? false : !payload.born)
            positions[uid]=row
        else
            state.world[:,row].=point
            state.uv[:,row].=middle_uv
        end
        edge_rows[key]=row
        return row
    end
    for cells in (state.triangles,state.quads),cell in cells
        width=length(cell)
        edge=ntuple(index->edge_middle(cell[index],cell[mod1(index+1,width)]),width)
        scale=width==3 ? 1.0/3.0 : 0.25
        center_uv=ntuple(axis->begin
            value=state.uv[axis,cell[1]]+state.uv[axis,cell[2]]+state.uv[axis,cell[3]]
            width==4 && (value+=state.uv[axis,cell[4]])
            value*scale
        end,2)
        # The center's world row is the edge-middle mean on a planar frame —
        # identical to `point_lift(center_uv)` under the affine lift — and the
        # on-surface lift of the averaged parameters on a curved patch.
        point=point_lift===nothing ? ntuple(axis->begin
            value=state.world[axis,edge[1]]+state.world[axis,edge[2]]+state.world[axis,edge[3]]
            width==4 && (value+=state.world[axis,edge[4]])
            value*scale
        end,3) : point_lift(center_uv...)
        center=_gmsh_fullquad_append_node!(state,birth(),owner,collect(center_uv),point,center_uv,
            parameters,retained,caller)
        for index in 1:width
            push!(quads,(cell[index],edge[index],center,edge[mod1(index-1,width)]))
            push!(labels,_native_recombine_birth!(state))
        end
    end
    empty!(state.triangles)
    empty!(state.triangle_tags)
    state.quads=quads
    state.quad_tags=labels
    length(state.node_tags)==count || throw(ArgumentError(
        "$caller: full-quad reserved rows disagree with actual shared identities"))
    return high_water[]
end

function _gmsh_native_surface_fullquad(m,tag,product,world,protected,frame,input,caller;
                                      reverse_initial::Bool=false,lift=nothing,
                                      middle_uv_of=nothing,point_lift=nothing)
    # The planar/OCC fronts relocate through the frame's lift; parametric
    # fronts (the ruled patch's own chart) hand in `lift` plus the boundary
    # maps `middle_uv_of`/`point_lift` since no plane frame exists.
    state=lift===nothing ? _native_recombine_state(product,world,frame,caller) :
        _native_recombine_state_lift(product,world,lift,caller)
    # The initial-only producer stores recovered UV triangles CCW. Native
    # subdivision consumes the already constructed GFace winding, so perform
    # that conversion before pairing and subdivision can allocate new cells.
    reverse_initial && (state.triangles=[(a,c,b) for (a,b,c) in state.triangles])
    _native_recombine_process!(state,m,tag,protected;minimum_quality=0.1)
    parameters=deepcopy(product.node_parameters)
    retained=copy(product.retained_nodes)
    # `input.node_high_water` sits above any identities reserved for
    # uncoarsened-Curve middles, so fresh births cannot reuse them.
    node_high_water=_gmsh_fullquad_subdivide!(state,m,tag,input.middles,parameters,retained,
        max(product.node_high_water,input.node_high_water),frame,caller;
        ctx_lc=input.ctx_lc,middle_uv_of,point_lift)
    # Restored Curve line segments remain protected in the second stage.
    restored=Set{Tuple{Int32,Int32}}()
    positions=Dict(uid=>Int32(row) for (row,uid) in enumerate(state.node_tags))
    for (ends,middle) in input.middles
        a=get(positions,ends[1],Int32(0));b=get(positions,ends[2],Int32(0))
        center=get(positions,middle.tag,Int32(0))
        a>0 && center>0 && push!(restored,minmax(a,center))
        b>0 && center>0 && push!(restored,minmax(center,b))
    end
    _native_recombine_process!(state,m,tag,restored;minimum_quality=1e-3)
    output=_native_recombine_export(product,state,tag,caller;retained,parameters,node_high_water)
    final_positions=Dict(uid=>row for (row,uid) in enumerate(state.node_tags))
    # Born middles are fresh identities absent from the lower producer — the
    # source-update map only carries retained halved middles.
    updates=Dict(middle.tag=>begin
        row=get(final_positions,middle.tag,0)
        row==0 ? middle.chord_world : ntuple(axis->state.world[axis,row],3)
    end for middle in values(input.middles) if !middle.born)
    output=_ModelSurfaceProduction(output.mesh,output.node_tags,output.node_owners,
        output.node_parameters,output.retained_nodes,output.node_high_water,
        output.element_tags,output.element_high_water,output.birth_rows,output.constructor_mesh,updates)
    return output,state.world[:,findall(state.alive)]
end
