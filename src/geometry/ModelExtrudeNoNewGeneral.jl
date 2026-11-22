# The general QuadTriNoNewVerts planner: Gmsh's QuadToTriEdgeGenerator and
# QuadToTriCreateElements (src/structured/QuadToTri3D.jl) on a lone region
# whose source the bounded catalogs do not cover — mixed triangle and
# quadrangle sources, sources with cells on a rotation axis, and every other
# source those catalogs refuse.
#
# Gmsh orders vertices by MVertex address, which differs between runs (12
# runs of one input give two element-type mixes). Here the order is Gmsh's
# creation order: points, curves, surfaces, the volume, each by tag, and
# within an entity the order its mesher creates them.

function _extrude_nonew_general_candidate(m::GeoModel,t::Int,caller)
    params=_extrude_gate(m,3,t)
    _extrude_is_nonew(params) || return false
    spec=get(m.meshing.extrude_specs,(3,t),nothing)
    link=get(m.meshing.extrude_sources,(3,t),nothing)
    spec!==nothing && link!==nothing && link[1]==2 || return false
    spec.type===:rotate && abs(spec.angle)>=2pi && return false
    source=abs(link[2])
    haskey(m.surfaces,source) && !haskey(m.meshing.extrude_sources,(2,source)) || return false
    isempty(get(m.embeds,(3,t),NTuple{2,Int}[])) || return false
    for surface in abs.(_extrude_nonew_boundary_surfaces(m,t,caller))
        isempty(get(m.embeds,(2,surface),NTuple{2,Int}[])) || return false
        _extrude_quadtri_regions(m,surface,caller)==[t] || return false
    end
    _extrude_nonew_unmerged_b4_incidence(m,t,caller) && return false
    return true
end

function _extrude_nonew_general_plan(m::GeoModel,t::Int,caller::AbstractString;
        min_angle_deg::Real=25.0,max_periodic_passes=8,
        size_field::Union{Nothing,AbstractSizeField}=nothing,_mesh_input=nothing)
    params=_extrude_gate(m,3,t)
    spec=m.meshing.extrude_specs[(3,t)]
    source=abs(m.meshing.extrude_sources[(3,t)][2])
    boundary=abs.(_extrude_nonew_boundary_surfaces(m,t,caller))
    top=0;laterals=Int[]
    for surface in boundary
        surface==source && continue
        link=get(m.meshing.extrude_sources,(2,surface),nothing)
        link===nothing && throw(ArgumentError(
            "$caller: QuadTriNoNewVerts boundary surface lacks extrusion provenance"))
        if link[1]==2 && abs(link[2])==source
            top==0 || throw(ArgumentError("$caller: QuadTriNoNewVerts has multiple copied caps"))
            top=surface
        elseif link[1]==1
            push!(laterals,surface)
        end
    end
    top!=0 || throw(ArgumentError("$caller: QuadTriNoNewVerts requires one copied cap"))
    sort!(laterals)
    produced_source=_mesh_input isa _ModelMeshAuthority &&
        haskey(_mesh_input.producer_rows,(2,Int32(source))) ?
        _model_mesh_product(_mesh_input,2,source;caller,producer_order=true) : nothing
    source_mesh=produced_source===nothing ?
        mesh_model_surface(m,source;min_angle_deg,max_periodic_passes,size_field,
            _mesh_input) : produced_source
    tris=NTuple{3,Int}[];quads=NTuple{4,Int}[]
    if source_mesh isa Mesh
        for column in axes(source_mesh.tris,2)
            push!(tris,ntuple(k->Int(source_mesh.tris[k,column]),3))
        end
    else
        for block in source_mesh.blocks
            block isa ElementBlock && msh_dimension(block.msh)==2 || continue
            block.msh in (2,3) || throw(ArgumentError(
                "$caller: QuadTriNoNewVerts needs a linear source mesh"))
            for column in axes(block.nodes,2)
                if block.msh==2
                    push!(tris,ntuple(k->Int(block.nodes[k,column]),3))
                else
                    push!(quads,ntuple(k->Int(block.nodes[k,column]),4))
                end
            end
        end
    end
    nsource=nnodes(source_mesh)
    levels,refs=_extrude_nonew_levels(params,caller;source_nodes=nsource,extra_nodes=0,
        cells_per_interval=Base.checked_mul(16,max(1,length(tris)+length(quads))))
    nlev=length(levels)
    cols=_extrude_volume_columns(m,t,source,source_mesh,params,spec,levels,caller)
    all(p->all(isfinite,p),cols) || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts columns are not finite"))
    position_node=Dict{NTuple{3,Float64},Int}(cols[1,i]=>i for i in axes(cols,2))
    # Bounding curves: their chains (source rows) and the lateral of each.
    curves=_model_projection_surface_curves(m,source)
    point_row=Dict{Int,Int}();curve_rank=Dict{Int,Tuple{Int,Int}}()
    boundary_curves=Int[];boundary_chains=Vector{Int32}[];boundary_widths=Int[]
    lateral_of=Dict{Int,Int}()
    for s in laterals
        link=m.meshing.extrude_sources[(2,s)]
        lateral_of[abs(link[2])]=s
    end
    lateral_lines=Vector{Vector{Tuple{Int,Int}}}()
    lateral_tags=Int[]
    for curve in curves
        chain=_extrude_curve_nodes(m,curve,caller)
        rows=Int[get(position_node,p,0) for p in chain]
        all(>(0),rows) || throw(ArgumentError(
            "$caller: QuadTriNoNewVerts source boundary leaves the source mesh"))
        push!(boundary_curves,curve);push!(boundary_chains,Int32.(rows));push!(boundary_widths,length(rows))
        a,b=m.curves[curve]
        point_row[a]=rows[1];point_row[b]=rows[end]
        for (index,row) in enumerate(rows)
            (index==1 || index==length(rows)) && continue
            curve_rank[row]=(curve,index)
        end
        s=get(lateral_of,curve,0)
        if s!=0
            closed=a==b
            push!(lateral_lines,[(rows[i],rows[closed ? mod1(i+1,length(rows)) : i+1])
                for i in 1:(closed ? length(rows) : length(rows)-1)])
            push!(lateral_tags,s)
        end
    end
    row_point=Dict(row=>point for (point,row) in point_row)
    # Entities the swept vertices land on.
    region_curves=Set{Int}()
    for surface in boundary,loop in m.surfaces[surface],signed in m.loops[loop]
        push!(region_curves,abs(signed))
    end
    connector=Dict{Int,Int}();copied_curve=Dict{Int,Int}()
    for c in region_curves
        link=get(m.meshing.extrude_sources,(1,c),nothing)
        link===nothing && continue
        link[1]==0 && (connector[abs(link[2])]=c)
        link[1]==1 && (copied_curve[abs(link[2])]=c)
    end
    # Creation key of source row i at level l (0-based).
    function creation_key(i,l)
        if haskey(row_point,i)
            point=row_point[i]
            l==0 && return (0,point,0,0)
            c=get(connector,point,0)
            c==0 && return (0,point,0,0)   # on the axis: the point itself
            l==nlev-1 && return (0,m.curves[c][2],0,0)
            return (1,c,l,0)
        elseif haskey(curve_rank,i)
            curve,index=curve_rank[i]
            l==0 && return (1,curve,index,0)
            l==nlev-1 && return (1,get(copied_curve,curve,curve),index,0)
            return (2,get(lateral_of,curve,0),index,l)
        end
        l==0 && return (2,source,i,0)
        l==nlev-1 && return (2,top,i,0)
        return (3,t,i,l)
    end
    ids=Dict{NTuple{3,Float64},Int}();keys_of=Dict{NTuple{3,Float64},NTuple{4,Int}}()
    for i in axes(cols,2),l in 0:nlev-1
        p=cols[l+1,i]
        key=creation_key(i,l)
        old=get(keys_of,p,nothing)
        (old===nothing || key<old) && (keys_of[p]=key)
    end
    ordered=sort!(collect(keys(keys_of));by=p->keys_of[p])
    point_of=Vector{NTuple{3,Float64}}(undef,length(ordered))
    for (id,p) in enumerate(ordered)
        ids[p]=id;point_of[id]=p
    end
    vid=(v,l)->ids[cols[l+1,v]]
    # Source rows as ids: the port reads elements by vertex id.
    source_id=[ids[cols[1,i]] for i in axes(cols,2)]
    row_of=Dict(id=>i for (i,id) in enumerate(source_id))
    region=QuadToTri3D.QTRegion(
        [map(r->source_id[r],c) for c in tris],[map(r->source_id[r],c) for c in quads],
        collect(Int,params.layers),
        (v,l)->vid(row_of[v],l),id->(id<=length(point_of) ? point_of[id] :
            error("unknown QuadToTri vertex $id")),
        Set(source_id[r] for r in values(point_row)) ∪ Set(source_id[r] for r in keys(curve_rank)),
        spec.type===:rotate || spec.type===:twist,
        [[(source_id[a],source_id[b]) for (a,b) in lines] for lines in lateral_lines],
        fill(params.recomb_laterals ? 0 : 2,length(lateral_lines)),
        fill(params.recomb_laterals,length(lateral_lines)),
        false,true,length(point_of)+1)
    result=QuadToTri3D.qt_mesh_region(region)
    out=result.mesh
    coordinate(id)=id<=length(point_of) ? point_of[id] :
        out.new_coords[findfirst(==(id),out.new_ids)]
    P(c)=map(coordinate,c)
    sweep=_ExtrudeVolumeSweep(t,true,cols,
        NTuple{6,NTuple{3,Float64}}[],NTuple{4,NTuple{3,Float64}}[P(c) for c in out.tets],
        NTuple{8,NTuple{3,Float64}}[P(c) for c in out.hexes],
        NTuple{6,NTuple{3,Float64}}[P(c) for c in out.prisms],
        NTuple{5,NTuple{3,Float64}}[P(c) for c in out.pyramids],
        copy(out.new_coords),Bool[false])
    edges=Set{NTuple{2,NTuple{3,Float64}}}(
        _extrude_ekey(coordinate(a),coordinate(b)) for (a,b) in result.state.quadToTri_edges
        if a<=length(point_of) && b<=length(point_of))
    top_cells=NTuple{3,NTuple{3,Float64}}[P(c) for c in result.top]
    cells=Tuple{Vararg{Int32}}[Int32.(c) for c in tris]
    append!(cells,Tuple{Vararg{Int32}}[Int32.(c) for c in quads])
    catalog=_ExtrudeNoNewMultiCellCatalog(cells,boundary_curves,boundary_chains,
        boundary_widths,levels,refs,top_cells)
    return (sweep=sweep,edges=edges,source_mesh=source_mesh,source_tag=source,
        top_tag=top,lateral_tags=sort!(lateral_tags),catalog=catalog)
end

# The general planner's cap: its own triangles (MeshQuadToTriTopSurface).
function _extrude_nonew_general_top(cells)
    lookup=Dict{NTuple{3,Float64},Int32}();coords=NTuple{3,Float64}[]
    nodes=Matrix{Int32}(undef,3,length(cells))
    for (column,cell) in enumerate(cells),k in 1:3
        nodes[k,column]=_extrude_node!(lookup,coords,cell[k])
    end
    out=Matrix{Float64}(undef,3,length(coords))
    for (i,p) in enumerate(coords),k in 1:3
        out[k,i]=p[k]
    end
    return Mesh(out;tris=nodes)
end
