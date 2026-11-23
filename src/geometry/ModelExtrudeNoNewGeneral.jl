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
        owners=_extrude_quadtri_regions(m,surface,caller)
        owners==[t] && continue
        # A Source swept by several NoNewVerts Regions: Gmsh runs QuadToTri
        # on each Region alone (the Source's cells are never edited), as two
        # sweeps of one Source mesh show (both Regions = the lone result).
        surface==source && all(owners) do owner
            entry=_extrude_entity_params(m,3,owner,caller)
            entry!==nothing && _extrude_is_nonew(entry[1]) &&
                entry[3][1]==2 && abs(entry[3][2])==source
        end && continue
        # A lateral of two NoNewVerts Regions with the same layer record:
        # the first Region meshes it (its quadrangles free), the second reads
        # its cells (QuadToTriGetRegionDiags).
        _extrude_nonew_general_shared_lateral(m,t,surface,owners,caller) && continue
        # A lateral of an unrecombined layered Region: Gmsh's global
        # SubdivideExtrudedMesh splits it first and the QuadToTri Region is
        # meshed again on those diagonals (meshQuadToTriRegionAfterGlobalSubdivide).
        _extrude_nonew_general_ordinary_partner(m,t,surface,owners,caller) || return false
    end
    _extrude_nonew_unmerged_b4_incidence(m,t,caller) && return false
    return true
end

function _extrude_nonew_general_shared_lateral(m::GeoModel,t::Int,surface::Int,owners,caller)
    length(owners)==2 && t in owners || return false
    s_link=get(m.meshing.extrude_sources,(2,surface),nothing)
    s_link!==nothing && s_link[1]==1 || return false
    t_entry=_extrude_entity_params(m,3,t,caller)
    for owner in owners
        entry=_extrude_entity_params(m,3,owner,caller)
        entry!==nothing && _extrude_is_nonew(entry[1]) && entry[3][1]==2 || return false
        isequal(entry[1].layers,t_entry[1].layers) && isequal(entry[1].heights,t_entry[1].heights) &&
            isequal(entry[2],t_entry[2]) || return false
        abs(entry[3][2])==surface && return false
        abs(s_link[2]) in _model_projection_surface_curves(m,abs(entry[3][2])) || return false
    end
    return true
end

function _extrude_nonew_general_ordinary_partner(m::GeoModel,t::Int,surface::Int,owners,caller)
    length(owners)==2 && t in owners || return false
    s_link=get(m.meshing.extrude_sources,(2,surface),nothing)
    s_link!==nothing && s_link[1]==1 || return false
    t_entry=_extrude_entity_params(m,3,t,caller)
    partner=only(o for o in owners if o!=t)
    entry=_extrude_entity_params(m,3,partner,caller)
    entry===nothing && return false
    params,spec,link=entry
    (_extrude_is_nonew(params) || params.recombine || params.quad_to_tri!==:none) && return false
    link[1]==2 && abs(link[2])!=surface || return false
    # A recombined Source cannot be swept without Recombine; Gmsh stops there.
    haskey(m.meshing.recombine,(2,abs(link[2]))) && return false
    isequal(params.layers,t_entry[1].layers) && isequal(params.heights,t_entry[1].heights) &&
        isequal(spec,t_entry[2]) || return false
    return abs(s_link[2]) in _model_projection_surface_curves(m,abs(link[2]))
end

# The diagonal set of Gmsh's global SubdivideExtrudedMesh over every
# unrecombined layered Region of the model, run before any seed exists.
function _extrude_nonew_general_subdivision(m::GeoModel,caller;_mesh_input=nothing)
    sweeps=_ExtrudeVolumeSweep[]
    for tag in sort!(collect(keys(m.volumes)))
        _extrude_gate(m,3,tag)===nothing && continue
        entry=_extrude_entity_params(m,3,tag,caller;_mesh_input)
        entry===nothing && continue
        params,spec,link=entry
        (_extrude_is_nonew(params) || params.recombine) && continue
        push!(sweeps,_extrude_volume_sweep(m,tag,params,spec,link[2],caller;_mesh_input))
    end
    return _extrude_subdivide!(sweeps,caller)
end

# The wall an unrecombined partner's subdivision leaves on lateral `s`:
# two triangles per layer quadrangle, split on the subdivision diagonal.
function _extrude_nonew_general_subdivided_wall(cols,lines,edges,caller)
    coords=NTuple{3,Float64}[];index=Dict{NTuple{3,Float64},Int32}()
    node(p)=get!(index,p) do
        push!(coords,p);Int32(length(coords))
    end
    tris=NTuple{3,Int32}[]
    has(p,q)=(p,q) in edges || (q,p) in edges
    for (a,b) in lines,l in 1:size(cols,1)-1
        A0=cols[l,a];B0=cols[l,b];A1=cols[l+1,a];B1=cols[l+1,b]
        if has(A0,B1)
            push!(tris,(node(A0),node(B0),node(B1)),(node(A0),node(B1),node(A1)))
        elseif has(B0,A1)
            push!(tris,(node(A0),node(B0),node(A1)),(node(B0),node(B1),node(A1)))
        else
            throw(ArgumentError("$caller: QuadTriNoNewVerts shared wall quadrangle " *
                "has no subdivision diagonal"))
        end
    end
    xyz=Matrix{Float64}(undef,3,length(coords))
    for (k,p) in enumerate(coords);xyz[:,k].=p;end
    nodes=Matrix{Int32}(undef,3,length(tris))
    for (k,c) in enumerate(tris);nodes[:,k].=c;end
    return MixedMesh(xyz,[ElementBlock(2,nodes)])
end

function _extrude_nonew_general_plan(m::GeoModel,t::Int,caller::AbstractString;
        min_angle_deg::Real=25.0,max_periodic_passes=8,
        size_field::Union{Nothing,AbstractSizeField}=nothing,_mesh_input=nothing,
        _given_walls=nothing)
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
        Int[_extrude_quadtri_regions(m,s,caller)==[t] ? (params.recomb_laterals ? 0 : 2) : 1
            for s in lateral_tags],
        Bool[_extrude_quadtri_regions(m,s,caller)==[t] ? params.recomb_laterals : true
            for s in lateral_tags],
        false,true,length(point_of)+1)
    # Laterals shared with the other NoNewVerts Region: given cells when that
    # Region was meshed first, else free quadrangles this Region splits.
    given=Dict{Int,Any}();free=Set{Int}();initial=Dict{Int,Any}()
    subdivision=nothing
    for (index,s) in enumerate(lateral_tags)
        owners=_extrude_quadtri_regions(m,s,caller)
        owners==[t] && continue
        if _extrude_nonew_general_ordinary_partner(m,t,s,owners,caller)
            subdivision===nothing &&
                (subdivision=_extrude_nonew_general_subdivision(deepcopy(m),caller;_mesh_input))
            wall=_extrude_nonew_general_subdivided_wall(cols,lateral_lines[index],subdivision,caller)
            given[index]=_extrude_nonew_general_wall_cells(wall,ids,lateral_lines[index],
                source_id,region,caller)
            continue
        end
        part=_given_walls===nothing ? nothing : get(_given_walls,s,nothing)
        if part===nothing
            push!(free,index)
            # A transfinite wall is meshed before either Region; its
            # triangles' diagonals are adjustable, not fixed.
            if haskey(m.meshing.transfinite_surfaces,s)
                wall=_transfinite_surface_mesh(m,s,Dict{Tuple{Int,Float64},Float64}(),caller;size_field)
                initial[index]=_extrude_nonew_general_wall_cells(wall,ids,lateral_lines[index],
                    source_id,region,caller)
            end
        else
            given[index]=_extrude_nonew_general_wall_cells(part,ids,lateral_lines[index],
                source_id,region,caller)
        end
    end
    # A transfinite cap is meshed before the Region (meshGFace's order); the
    # Region reads its triangles and their diagonals.
    top_override=nothing
    if haskey(m.meshing.transfinite_surfaces,top)
        actual=_mesh_input isa _ModelMeshAuthority ?
            _model_mesh_product(_mesh_input,2,top;caller) : nothing
        cap=actual!==nothing && any(block->!isempty(block.nodes),actual.blocks) ? actual :
            mesh_model_surface(m,top;min_angle_deg,max_periodic_passes,size_field,_mesh_input)
        top_override=_extrude_nonew_general_cap(cap,ids,region,caller)
    elseif _mesh_input isa _ModelMeshAuthority
        # A cap the model already holds as triangles (a retained product being
        # projected) is read as it stands, as Gmsh reads a meshed face.
        actual=_model_mesh_product(_mesh_input,2,top;caller)
        if (actual isa Mesh && size(actual.tris,2)>0) || (actual isa MixedMesh &&
           any(b->b isa ElementBlock && b.msh==2 && !isempty(b.nodes),actual.blocks) &&
           all(b->!(b isa ElementBlock) || b.msh in (1,2) || isempty(b.nodes),actual.blocks))
            top_override=_extrude_nonew_general_cap(actual,ids,region,caller)
        end
    end
    result=QuadToTri3D.qt_mesh_region(region;top=top_override,given,free,initial)
    out=result.mesh
    # Gmsh logs these and meshes on; its API call then fails with the text.
    isempty(out.errors) || throw(ArgumentError("$caller: "*first(out.errors)))
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

# The cells of a shared lateral the other Region already emitted, in this
# planner's vertex ids, grouped per source line as createQuaTri made them.
function _extrude_nonew_general_wall_cells(part,ids,lines,source_id,region,caller)
    nlev=sum(region.layers)
    blocks=part isa Mesh ? [part.tris] : [b.nodes for b in part.blocks if b isa ElementBlock && b.msh in (2,3)]
    cells=Vector{Tuple{Vararg{Int}}}()
    for nodes in blocks,column in axes(nodes,2)
        cell=Tuple(begin
            p=(part.coords[1,v],part.coords[2,v],part.coords[3,v])
            id=get(ids,p,0)
            id==0 && throw(ArgumentError(
                "$caller: QuadTriNoNewVerts shared wall node is off this Region's columns"))
            id
        end for v in @view nodes[:,column])
        push!(cells,cell)
    end
    per_line=Vector{Vector{Tuple{Vararg{Int}}}}()
    used=falses(length(cells))
    for (a,b) in lines
        column=Set{Int}(region.vid(v,l) for v in (source_id[a],source_id[b]) for l in 0:nlev)
        mine=Tuple{Vararg{Int}}[]
        for (k,cell) in enumerate(cells)
            used[k] && continue
            all(v->v in column,cell) || continue
            used[k]=true;push!(mine,cell)
        end
        push!(per_line,mine)
    end
    all(used) || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts shared wall cell spans two source lines"))
    return per_line
end

# A cap meshed ahead of the Region, in the planner's vertex ids: its
# triangles, and as diagonals the triangle edges that cross a swept source
# quadrangle. A cap with nodes off the swept columns, quadrangles, or a
# triangle across two source cells cannot be followed; Gmsh stops there too.
function _extrude_nonew_general_cap(cap,ids,region,caller)
    top=sum(region.layers)
    blocks=cap isa Mesh ? [(2,cap.tris)] : [(Int(b.msh),b.nodes) for b in cap.blocks if b isa ElementBlock]
    cells=NTuple{3,Int}[]
    for (msh,nodes) in blocks
        msh==1 && continue
        msh==2 || throw(ArgumentError(
            "$caller: QuadTriNoNewVerts cap mesh cannot be followed by the Region"))
        for column in axes(nodes,2)
            cell=ntuple(3) do k
                p=(cap.coords[1,nodes[k,column]],cap.coords[2,nodes[k,column]],cap.coords[3,nodes[k,column]])
                id=get(ids,p,0)
                id==0 && throw(ArgumentError(
                    "$caller: QuadTriNoNewVerts cap mesh has a node off the swept columns"))
                id
            end
            push!(cells,cell)
        end
    end
    quads=[map(v->region.vid(v,top),q) for q in region.quads]
    edge_quads=Dict{QuadToTri3D.Edge,Int}()
    diagonal_of=Dict{QuadToTri3D.Edge,Int}()
    for (k,q) in enumerate(quads)
        for i in 1:4
            edge_quads[QuadToTri3D._ekey(q[i],q[mod1(i+1,4)])]=k
        end
        diagonal_of[QuadToTri3D._ekey(q[1],q[3])]=k
        diagonal_of[QuadToTri3D._ekey(q[2],q[4])]=k
    end
    diagonals=QuadToTri3D.Edge[]
    for cell in cells,i in 1:3
        e=QuadToTri3D._ekey(cell[i],cell[mod1(i+1,3)])
        haskey(diagonal_of,e) && !(e in diagonals) && push!(diagonals,e)
    end
    length(diagonals)==length(quads) || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts finalized surfaces do not exactly cover the typed volume boundary"))
    return cells,diagonals
end

# The general planner's cap: its own triangles (MeshQuadToTriTopSurface),
# its nodes in the Source's node order (copyMesh copies the Source's
# vertices in order before the cap is triangulated).
function _extrude_nonew_general_top(cells,cols=nothing)
    lookup=Dict{NTuple{3,Float64},Int32}();coords=NTuple{3,Float64}[]
    if cols!==nothing
        for i in axes(cols,2)
            _extrude_node!(lookup,coords,cols[end,i])
        end
    end
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
