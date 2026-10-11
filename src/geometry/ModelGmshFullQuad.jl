# Generation owns a temporary coarse Curve view. Original lower products and
# their identities remain detached from the temporary Line constructors.
struct _NativeFullQuadMiddle
    tag::UInt64
    owner::Tuple{Int,Int32}
    parameter::Union{Nothing,Vector{Float64}}
    world::NTuple{3,Float64}
    chord_world::NTuple{3,Float64}
end
_NativeFullQuadMiddle(tag,owner,parameter,world)=
    _NativeFullQuadMiddle(tag,owner,parameter,world,world)

function _gmsh_fullquad_input(m,tag,input,element_high_water,options,caller)
    middles=Dict{Tuple{UInt64,UInt64},_NativeFullQuadMiddle}()
    compound=any(pair->first(pair)==2 && tag in last(pair),m.meshing.compounds)
    # Mesh-periodic slaves copy an already generated master. This Plane
    # producer has no CAD periodic parametrization or seam to reject.
    if compound
        return (authority=input,middles,element_high_water,admitted=false,ctx_lc=options.ctx_lc)
    end
    authority=deepcopy(input)
    source=authority.mesh
    data=source.entity_data
    _,embedded=_model_surface_embedding_tags(m,tag,caller)
    curves=vcat(embedded,_gmsh_native_surface_curves(m,tag))
    high_water=UInt64(element_high_water)
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
        if isodd(length(cells))
            _model_mesh_error(options,"1D mesh cannot be divided by 2")
        end
        coarse=NTuple{2,Int32}[]
        labels=UInt64[]
        for index in 1:2:length(cells)-1
            first,middle=cells[index]
            next,last=cells[index+1]
            middle==next || throw(ArgumentError(
                "$caller: Curve[$curve] consecutive Line constructors do not share their middle identity"))
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
        connectivity=Matrix{Int32}(undef,2,length(coarse))
        for (column,cell) in enumerate(coarse)
            connectivity[:,column].=cell
        end
        authority.products[owner]=_ModelMeshProduct(owner,original.rows,
            ElementBlock[ElementBlock(1,connectivity)],Vector{UInt64}[labels])
    end
    return (authority,middles,element_high_water=high_water,admitted=!isempty(middles),ctx_lc=options.ctx_lc)
end

function _gmsh_fullquad_size_options(options::_ModelMesh1DOptions)
    names=fieldnames(_ModelMesh1DOptions)
    return _ModelMesh1DOptions((name===:lc_factor ? 2.0*getfield(options,name) :
        getfield(options,name) for name in names)...)
end

function _gmsh_fullquad_append_node!(state,tag,owner,parameter,point,uv,
                                    parameters,retained,caller)
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
    push!(retained,owner[1]<2)
    return Int32(length(state.node_tags))
end

function _gmsh_fullquad_subdivide!(state,m,tag,middles,parameters,retained,
                                  node_high_water,frame,caller;ctx_lc=0.0)
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
        point=payload===nothing ? ntuple(axis->0.5*(state.world[axis,a]+state.world[axis,b]),3) :
            payload.world
        middle_uv=ntuple(axis->0.5*(state.uv[axis,a]+state.uv[axis,b]),2)
        if payload!==nothing && payload.owner[1]==1
            value=payload.parameter===nothing || isempty(payload.parameter) ? 0.0 : payload.parameter[1]
            middle_uv=frame.occ ? _occ_edge_reparam_on_face(m,m.surface_geometry[tag],
                Int(payload.owner[2]),value,1,ctx_lc,caller) :
                _gmsh_native_surface_parameters(frame,
                    _gmsh_native_curve_sample(m,Int(payload.owner[2]),value,caller))
        elseif payload!==nothing && payload.owner[1]==0
            middle_uv=_gmsh_native_surface_parameters(frame,point)
        end
        uid=payload===nothing ? birth() : payload.tag
        row=get(positions,uid,Int32(0))
        if row==0
            row=_gmsh_fullquad_append_node!(state,uid,
                payload===nothing ? owner : payload.owner,
                payload===nothing ? collect(middle_uv) : deepcopy(payload.parameter),point,middle_uv,
                parameters,retained,caller)
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
        point=ntuple(axis->begin
            value=state.world[axis,edge[1]]+state.world[axis,edge[2]]+state.world[axis,edge[3]]
            width==4 && (value+=state.world[axis,edge[4]])
            value*scale
        end,3)
        center_uv=ntuple(axis->begin
            value=state.uv[axis,cell[1]]+state.uv[axis,cell[2]]+state.uv[axis,cell[3]]
            width==4 && (value+=state.uv[axis,cell[4]])
            value*scale
        end,2)
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
                                      reverse_initial::Bool=false)
    state=_native_recombine_state(product,world,frame,caller)
    # The initial-only producer stores recovered UV triangles CCW. Native
    # subdivision consumes the already constructed GFace winding, so perform
    # that conversion before pairing and subdivision can allocate new cells.
    reverse_initial && (state.triangles=[(a,c,b) for (a,b,c) in state.triangles])
    _native_recombine_process!(state,m,tag,protected;minimum_quality=0.1)
    parameters=deepcopy(product.node_parameters)
    retained=copy(product.retained_nodes)
    node_high_water=_gmsh_fullquad_subdivide!(state,m,tag,input.middles,parameters,retained,
        product.node_high_water,frame,caller;ctx_lc=input.ctx_lc)
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
    updates=Dict(middle.tag=>begin
        row=get(final_positions,middle.tag,0)
        row==0 ? middle.chord_world : ntuple(axis->state.world[axis,row],3)
    end for middle in values(input.middles))
    output=_ModelSurfaceProduction(output.mesh,output.node_tags,output.node_owners,
        output.node_parameters,output.retained_nodes,output.node_high_water,
        output.element_tags,output.element_high_water,output.birth_rows,output.constructor_mesh,updates)
    return output,state.world[:,findall(state.alive)]
end
