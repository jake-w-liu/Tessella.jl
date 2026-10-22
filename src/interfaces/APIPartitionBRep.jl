# Actual lower BRep transport across every complete cell family at any
# catalogued order (Point15; Line surfaces; Tri/Quad surfaces;
# Tet/Hex/Prism/Pyr volumes). Interface seams are derived
# from literal node UID incidence, with one authoritative original cell per
# UID. Facet and edge corner tables come from the canonical
# `MeshEntityTopology` outward cycles, which share the MSH pointsGenerators
# order, so an order-p cell carries n_corners+(edge_ordinal-1)*(p-1)+i edge
# nodes along the table edge direction followed by per-face interior nodes in
# `_FACES` order, and the face-interior offsets accumulate triangle interiors
# before quadrangle interiors exactly as pointsGenerators appends them.
const _PARTITION_SURFACE_TYPES=Tuple(sort!(Int[msh for (msh,spec) in
    MSH_CATALOG if spec.family in (:tri,:qua) &&
    !spec.serendipity && spec.order>=1]))
const _PARTITION_VOLUME_TYPES=Tuple(sort!(Int[msh for (msh,spec) in
    MSH_CATALOG if spec.family in (:tet,:hex,:pri,:pyr) &&
    !spec.serendipity && spec.order>=1]))
const _PARTITION_CURVE_TYPES=Tuple(sort!(Int[msh for (msh,spec) in
    MSH_CATALOG if spec.family==:lin && !spec.serendipity &&
    spec.order>=1]))
# Non-volume admitted types — the bare-Tet4 producer routes to the BRep path
# whenever any declared lower cell exists at any complete order.
const _PARTITION_LOWER_TYPES=(15,_PARTITION_CURVE_TYPES...,
    _PARTITION_SURFACE_TYPES...)
_partition_support(uids)=Tuple(sort!(collect(uids)))
# Linear-family base type for any admitted cell type, resolved through the
# MSH catalog family rather than a numeric offset so third-order types map
# the same way second-order types do.
function _partition_linear_base(msh)
    family=Elements.msh_family(msh)
    family==:lin && return 1;family==:tri && return 2;family==:qua && return 3
    family==:tet && return 4;family==:hex && return 5;family==:pri && return 6
    family==:pyr && return 7
    return Int(msh)
end
# Corner-edge index pairs for any admitted cell type: the canonical tables are
# keyed on linear types, so higher-order cells resolve through their family.
_partition_edge_patterns(msh)=_simplex_edge_patterns(_partition_linear_base(msh))
# `(corners, edge node count per edge)` for an admitted surface facet type.
function _partition_surface_shape(msh)
    spec=Elements.msh_spec(msh)
    corners=spec.family==:tri ? 3 : spec.family==:qua ? 4 : 0
    return corners,spec.order-1
end
# Reverse a linear corner cycle, preserving the first corner.
_partition_reversed_corners(b)=ntuple(i->i==1 ? b[1] : b[length(b)+2-i],length(b))
# Source-index permutation `src` for a corner map σ (new slot σ(k) holds old
# vertex k) on a Tri/Quad cell of order p, so `transformed[i]=u[src[i]]`.
# Every admitted surface node sits on the integer monomial lattice that the
# pointsGenerators translation enumerates — the first nc monomials are the
# corner lattice points — and σ induces the unique affine map T sending
# old vertex k's lattice point to new vertex σ(k)'s. Each old slot's image
# coordinate resolves to exactly one new slot, so interior nodes of any
# order permute correctly without per-order tables.
function _partition_facet_perm(msh,σ)
    spec=msh_spec(msh);p=spec.order
    L=spec.family==:tri ? _tri_monomials(p) :
      spec.family==:qua ? _qua_monomials(p) :
      throw(ArgumentError("partition: facet permutations need Tri/Quad cells"))
    nc=length(σ);verts=Tuple(L[i] for i in 1:nc)
    c=verts[σ[1]]
    e1=Tuple(verts[σ[2]][d]-c[d] for d in 1:2)
    e2=Tuple(verts[σ[nc]][d]-c[d] for d in 1:2)
    lut=Dict(m=>i for (i,m) in enumerate(L))
    src=Vector{Int}(undef,length(L))
    for (i,m) in enumerate(L)
        t=(c[1]+(m[1]*e1[1]+m[2]*e2[1])÷p,
           c[2]+(m[1]*e1[2]+m[2]*e2[2])÷p)
        src[lut[t]]=i
    end
    return Tuple(src)
end
# Perms depend only on the cell type; every caller runs serially, but the
# cache is locked anyway so a future threaded caller stays safe.
const _PARTITION_FACET_PERM_LOCK=ReentrantLock()
const _PARTITION_FACET_PERM_CACHE=Dict{Tuple{Int,Tuple{Vararg{Int}}},
                                       Tuple{Vararg{Int}}}()
function _partition_facet_perm_cached(msh,σ)
    return lock(_PARTITION_FACET_PERM_LOCK) do
        get!(_PARTITION_FACET_PERM_CACHE,(msh,σ)) do
            _partition_facet_perm(msh,σ)
        end
    end
end
# Reverse an oriented cell ordering while preserving the family's node layout:
# the reversed corner cycle maps old vertex k to slot 1 or nc+2-k, and the
# affine lattice map carries edge and interior nodes to their slots at any
# order.
function _partition_reversed(msh,u)
    nc=Elements.msh_family(msh)==:tri ? 3 : 4
    src=_partition_facet_perm_cached(msh,ntuple(k->k==1 ? 1 : nc+2-k,nc))
    return ntuple(i->u[src[i]],length(u))
end
# Cycle-equality on oriented layouts: a equals some cyclic rotation of b,
# with edge and interior nodes rotating consistently through the same
# affine lattice map.
function _partition_cycle_equal(amsh,a,b)
    length(a)==length(b) || return false
    nc=Elements.msh_family(amsh)==:tri ? 3 : 4
    for shift in 0:nc-1
        src=_partition_facet_perm_cached(amsh,ntuple(k->mod1(k-shift,nc),nc))
        all(i->a[i]==b[src[i]],eachindex(a)) && return true
    end
    return false
end
_partition_cycle_reverse(amsh,a,bmsh,b)=
    _partition_cycle_equal(amsh,a,_partition_reversed(bmsh,b))
# `(facet surface msh type, cell node tuple)` for each outward facet of a
# volume family at its canonical order. An order-p cell appends p-1 edge
# nodes along every consecutive corner pair — taken in the facet's own
# direction, so a facet edge running against the canonical edge table picks
# the volume's edge nodes in reverse — followed by the facet's interior
# slots, whose offsets accumulate in `_FACES` order.
function _partition_volume_facets(msh)
    spec=Elements.msh_spec(msh);p=spec.order
    linear=_partition_linear_base(msh)
    corners=Elements.msh_spec(linear).nnodes
    edges=_simplex_edge_patterns(linear)
    edge_index=Dict(minmax(a,b)=>i for (i,(a,b)) in enumerate(edges))
    facets=Tuple{Int,Tuple{Vararg{Int}}}[]
    offset=corners+(p-1)*length(edges)
    for (family,patterns) in
            ((:tri,_simplex_face_patterns(linear,3)),
             (:qua,_simplex_face_patterns(linear,4)))
        interior=family==:tri ? (p-1)*(p-2)÷2 : (p-1)^2
        isempty(patterns) && continue
        surface_type=Elements._msh_type_exact(family,p,false,"partition")
        for pattern in patterns
            full=collect(Int,pattern)
            for k in eachindex(pattern)
                pair=(pattern[k],pattern[mod1(k+1,length(pattern))])
                ordinal=edge_index[minmax(pair...)]
                base=corners+(ordinal-1)*(p-1)
                forward=edges[ordinal]==pair
                for i in 1:p-1
                    push!(full,forward ? base+i : base+p-i)
                end
            end
            for i in 1:interior;push!(full,offset+i);end
            offset+=interior
            push!(facets,(surface_type,Tuple(full)))
        end
    end
    return facets
end

function _partition_brep_surface_borders(faces,caller;workers=1)
    face_edges=Vector{Vector{Tuple{Tuple{UInt64,UInt64},Tuple{UInt64,UInt64},
                                  Union{Nothing,Tuple{Vararg{UInt64}}}}}}(undef,length(faces))
    if workers>1
        Threads.@threads for index in eachindex(faces)
            face=faces[index];corners,e=_partition_surface_shape(face.msh)
            entries=Tuple{Tuple{UInt64,UInt64},Tuple{UInt64,UInt64},
                          Union{Nothing,Tuple{Vararg{UInt64}}}}[]
            for (order,(a,b)) in enumerate(_partition_edge_patterns(face.msh))
                edge=(face.uids[a],face.uids[b]);key=minmax(edge...)
                mid=e==0 ? nothing : ntuple(i->face.uids[corners+(order-1)*e+i],e)
                push!(entries,(key,edge,mid))
            end
            face_edges[index]=entries
        end
    else
        for index in eachindex(faces)
            face=faces[index];corners,e=_partition_surface_shape(face.msh)
            entries=Tuple{Tuple{UInt64,UInt64},Tuple{UInt64,UInt64},
                          Union{Nothing,Tuple{Vararg{UInt64}}}}[]
            for (order,(a,b)) in enumerate(_partition_edge_patterns(face.msh))
                edge=(face.uids[a],face.uids[b]);key=minmax(edge...)
                mid=e==0 ? nothing : ntuple(i->face.uids[corners+(order-1)*e+i],e)
                push!(entries,(key,edge,mid))
            end
            face_edges[index]=entries
        end
    end
    incidences=Dict{Int32,Dict{Tuple{UInt64,UInt64},Vector{Tuple{UInt64,UInt64}}}}()
    mid_incidences=Dict{Int32,Dict{Tuple{UInt64,UInt64},Vector{Union{Nothing,Tuple{Vararg{UInt64}}}}}}()
    for (face,entries) in zip(faces,face_edges),(key,edge,mid) in entries
        edges=get!(()->Dict{Tuple{UInt64,UInt64},Vector{Tuple{UInt64,UInt64}}}(),incidences,face.parent)
        push!(get!(()->Tuple{UInt64,UInt64}[],edges,key),edge)
        mids=get!(()->Dict{Tuple{UInt64,UInt64},Vector{Union{Nothing,Tuple{Vararg{UInt64}}}}}(),mid_incidences,face.parent)
        push!(get!(()->Union{Nothing,Tuple{Vararg{UInt64}}}[],mids,key),mid)
    end
    borders=Dict{Int32,Dict{Tuple{UInt64,UInt64},Tuple{UInt64,UInt64}}}()
    mid_borders=Dict{Int32,Dict{Tuple{UInt64,UInt64},Union{Nothing,Tuple{Vararg{UInt64}}}}}()
    for (parent,edges) in incidences
        border=Dict{Tuple{UInt64,UInt64},Tuple{UInt64,UInt64}}()
        mid_border=Dict{Tuple{UInt64,UInt64},Union{Nothing,Tuple{Vararg{UInt64}}}}()
        for (key,hits) in edges
            length(hits)<=2 || throw(ArgumentError("$caller: original Surface has a nonmanifold edge"))
            if length(hits)==2
                hits[1]==reverse(hits[2]) || throw(ArgumentError("$caller: original Surface edge winding disagrees"))
                # Co-bordant cells traverse the shared edge oppositely, so
                # their interior edge-node sequences must be mutual reversals.
                first,second=mid_incidences[parent][key]
                (first===nothing ? second===nothing : Tuple(first)==reverse(Tuple(second))) ||
                    throw(ArgumentError("$caller: original Surface edge node disagrees"))
            else
                border[key]=only(hits);mid_border[key]=only(mid_incidences[parent][key])
            end
        end
        borders[parent]=border;mid_borders[parent]=mid_border
    end
    return (;pairs=borders,mids=mid_borders)
end

# The native Region sign under which one declared Surface's stored cells
# already match this volume's own facet carrier; a Surface with no carrier
# keeps its stored sign (the caller's coverage check then reports the gap).
function _partition_region_surface_sign(parent,surface,signed,faces,facets)
    for face in faces
        face.parent==surface || continue
        hits=get(facets,_partition_support(face.uids),NamedTuple[])
        index=findfirst(hit->hit.volume.parent==parent,hits)
        index===nothing && continue
        return _partition_cycle_equal(face.msh,face.uids,hits[index].uids) ?
            Int(surface) : -Int(surface)
    end
    return Int(signed)
end

function _partition_brep_snapshot(model,cached,class,assignment,caller)
    class.public_tags!==nothing || throw(ArgumentError("$caller: lower BRep needs actual UID authority"))
    _high_order_overlay(cached)===nothing || throw(ArgumentError("$caller: lower BRep partitioning is unsupported for overlay-backed high-order caches"))
    catalog=[row for row in _cache_catalog(cached,class) if !isempty(row[4])]
    all(row->row[1] in (15,) || row[1] in _PARTITION_SURFACE_TYPES ||
        row[1] in _PARTITION_CURVE_TYPES || row[1] in _PARTITION_VOLUME_TYPES,catalog) &&
    any(row->row[1] in _PARTITION_VOLUME_TYPES,catalog) || throw(ArgumentError(
        "$caller: lower BRep partitioning supports complete cell families only"))
    cached isa MixedMesh && cached.partition_data!==nothing &&
        throw(ArgumentError("$caller: existing partition topology is unsupported"))
    public=class.public_tags
    length(public.element_tags)<=typemax(Int32) || throw(ArgumentError("$caller: lower BRep cell capacity exceeded"))
    for (_,record) in Model._discrete_mesh_records_model(model)
        all(uid->haskey(public.node_indices,uid),record.node_tags) &&
        all(uid->haskey(public.element_indices,uid),record.element_tags) || throw(ArgumentError(
            "$caller: lower BRep partitioning cannot omit caller-owned records outside the cache"))
    end
    assembly=_Dim01Assembly(_Dim01Node[],_Dim01Cell[],Dict{UInt64,Int32}(),
        NODE_TAG_MAX[],ELEMENT_TAG_MAX[],typemax(Int32),typemax(Int32),caller)
    _dim01_ingest!(assembly,model,cached,class;preserve_source_order=true)
    length(assembly.nodes)==nnodes(cached) && length(assembly.cells)==length(public.element_tags) ||
        throw(ArgumentError("$caller: cached and record BRep rows do not agree"))
    entities,_=_dim01_entities(model,model_entities(model),assembly,eachindex(assembly.nodes))
    workers=_partition_workers()
    blocks=ElementBlock[];owners=Vector{Int32}[];uids=Vector{UInt64}[]
    faces=NamedTuple[];volumes=NamedTuple[];lines=NamedTuple[];points=NamedTuple[]
    for (msh,dim,offset,cells,classified) in catalog
        physical=Int32[_dim01_physical(model,(dim,Int(owner))) for owner in classified]
        push!(blocks,ElementBlock(msh,cells,physical));push!(owners,copy(classified))
        push!(uids,copy(public.element_tags[offset+1:offset+size(cells,2)]))
        for column in axes(cells,2)
            rows=Tuple(cells[:,column]);labels=Tuple(public.node_tags[collect(rows)]);index=offset+column
            record=(;rows,uids=labels,msh=Int(msh),index,parent=classified[column],part=assignment[index])
            push!(msh in _PARTITION_SURFACE_TYPES ? faces :
                  msh in _PARTITION_VOLUME_TYPES ? volumes :
                  msh in _PARTITION_CURVE_TYPES ? lines : points,record)
        end
        if msh in _PARTITION_VOLUME_TYPES
            linear=_partition_linear_base(msh)
            positive=Vector{Bool}(undef,size(cells,2))
            if workers>1
                Threads.@threads for column in axes(cells,2)
                    positive[column]=Model._volume_cell_signed_volume(
                        cached.coords,cells,column,linear)>0
                end
            else
                for column in axes(cells,2)
                    positive[column]=Model._volume_cell_signed_volume(
                        cached.coords,cells,column,linear)>0
                end
            end
            all(positive) || throw(ArgumentError(
                "$caller: lower BRep needs positive oriented volume cells"))
        end
    end
    parents=sort!(unique(volume.parent for volume in volumes))
    surfaces=sort!(unique(face.parent for face in faces))
    curve_parents=sort!(unique(line.parent for line in lines))
    point_parents=sort!(unique(point.parent for point in points))
    all(parent->parent>0 && _model_entity_known(model,3,Int(parent)),parents) &&
    all(parent->parent>0 && _model_entity_known(model,2,Int(parent)),surfaces) &&
    all(parent->parent>0 && _model_entity_known(model,1,Int(parent)),curve_parents) &&
    all(parent->parent>0 && _model_entity_known(model,0,Int(parent)),point_parents) || throw(ArgumentError(
        "$caller: lower BRep requires Region, Surface, Curve and Point parents"))
    line_support=Dict{Tuple{UInt64,UInt64},NamedTuple}();curve_edges=Dict{Int32,Vector{NamedTuple}}()
    curve_nodes=Dict{Int32,Set{Int32}}();point_cells=Dict{Int32,NamedTuple}();point_rows=Set{Int32}()
    for line in lines
        key=minmax(line.uids[1],line.uids[2])
        haskey(line_support,key) && throw(ArgumentError(
            "$caller: original BRep repeats a Curve cell support"))
        line_support[key]=line
        ends=(line.uids[1],line.uids[2])
        order=Elements.msh_order(line.msh)
        push!(get!(()->NamedTuple[],curve_edges,line.parent),
            (uids=ends,rows=line.rows,mid=order==1 ? nothing :
                ntuple(i->line.uids[2+i],order-1)))
        nodes=get!(()->Set{Int32}(),curve_nodes,line.parent)
        for row in line.rows;push!(nodes,row);end
    end
    for point in points
        haskey(point_cells,point.parent) && throw(ArgumentError(
            "$caller: each declared lower Point needs exactly one actual Point15 cell"))
        point_cells[point.parent]=point
        point.rows[1] in point_rows && throw(ArgumentError("$caller: original BRep repeats a Point15 vertex"))
        push!(point_rows,point.rows[1])
    end
    curve_ends=Dict{Int32,NTuple{2,UInt64}}();used_points=Set{Int32}()
    curve_declared=Dict{Int32,Vector{NTuple{2,Int}}}()
    for parent in curve_parents
        edges=curve_edges[parent];outgoing=Dict{UInt64,Int}();incoming=Dict{UInt64,Int}()
        for (index,edge) in enumerate(edges)
            (haskey(outgoing,edge.uids[1]) || haskey(incoming,edge.uids[2])) && throw(ArgumentError(
                "$caller: declared Curve $parent is not one nonbranching directed chain"))
            outgoing[edge.uids[1]]=index;incoming[edge.uids[2]]=index
        end
        seed=something(findfirst(index->!haskey(incoming,edges[index].uids[1]),eachindex(edges)),1)
        seen=falses(length(edges));index=seed;last=seed
        while index!=0 && !seen[index]
            seen[index]=true;last=index;index=get(outgoing,edges[index].uids[2],0)
        end
        all(seen) || throw(ArgumentError(
            "$caller: declared Curve $parent is not one nonbranching directed chain"))
        boundary=_partition_declared_boundary(model,1,Int(parent),caller)
        if index==seed
            isempty(boundary) || throw(ArgumentError(
                "$caller: a closed declared Curve $parent cannot have Point boundaries"))
            curve_ends[parent]=(UInt64(0),UInt64(0))
            curve_declared[parent]=boundary
            continue
        end
        curve_ends[parent]=(edges[seed].uids[1],edges[last].uids[2])
        if !haskey(model.discrete,(1,Int(parent))) && length(boundary)==2
            # A native Curve's signed ends describe the actual stored chain:
            # `Reverse` or an opposite OCC edge use can leave the segments
            # running against the declared endpoint order.
            firstp,lastp=Int32(abs(boundary[1][2])),Int32(abs(boundary[2][2]))
            if haskey(point_cells,firstp) && haskey(point_cells,lastp)
                firstp,lastp=point_cells[firstp].uids[1]==curve_ends[parent][1] ?
                    (firstp,lastp) : (lastp,firstp)
                boundary=Tuple{Int,Int}[(0,firstp),(0,-lastp)]
            end
        end
        curve_declared[parent]=boundary
        length(boundary)==2 || throw(ArgumentError(
            "$caller: declared Curve $parent needs signed start and end Point boundaries"))
        covered_signs=Set{Int}()
        for (bdim,bsigned) in boundary
            bdim==0 || throw(ArgumentError(
                "$caller: declared Curve $parent bounds require Point entities"))
            point=Int32(abs(bsigned))
            haskey(point_cells,point) || throw(ArgumentError(
                "$caller: declared boundary Point $point has no actual Point15 cell"))
            side=bsigned>0 ? 1 : 2
            side in covered_signs && throw(ArgumentError(
                "$caller: declared Curve $parent repeats a boundary orientation"))
            push!(covered_signs,side);push!(used_points,point)
            node=side==1 ? curve_ends[parent][1] : curve_ends[parent][2]
            point_cells[point].uids[1]==node || throw(ArgumentError(
                "$caller: declared boundary Point $point does not lie on the actual Curve end"))
        end
    end
    surface_borders=_partition_brep_surface_borders(faces,caller;workers)
    used_curves=Set{Int32}()
    for parent in surfaces
        boundary=_partition_declared_boundary(model,2,Int(parent),caller)
        isempty(boundary) && continue
        all(entry->entry[1]==1,boundary) && allunique(Int32(abs(entry[2])) for entry in boundary) ||
            throw(ArgumentError(
                "$caller: each declared Surface boundary must list distinct Curve entities"))
        border=surface_borders.pairs[parent];covered=Set{Tuple{UInt64,UInt64}}()
        native=!haskey(model.discrete,(2,Int(parent)))
        for (_,signed) in boundary
            curve=Int32(abs(signed));haskey(curve_edges,curve) || throw(ArgumentError(
                "$caller: declared boundary Curve $curve has no actual Curve cells"))
            if native
                # A native Surface's Curve signs describe the actual border
                # traversal: the OCC face's use flag, its parametric frame
                # and Reverse/orientation flips can each leave the stored
                # cells running against the stored loop text.
                oriented=get(border,minmax(curve_edges[curve][1].uids...),nothing)
                oriented===nothing ||
                    (signed=oriented==curve_edges[curve][1].uids ? curve : -curve)
            end
            push!(used_curves,curve)
            for edge in curve_edges[curve]
                key=minmax(edge.uids...);oriented=get(border,key,nothing)
                oriented===nothing && throw(ArgumentError(
                    "$caller: declared Curve $curve cells must lie on the actual Surface border"))
                expected=signed>0 ? edge.uids : (edge.uids[2],edge.uids[1])
                oriented==expected || throw(ArgumentError(
                    "$caller: declared Curve $curve orientation disagrees with its Surface border"))
                # A higher-order Curve cell must share the Surface border's
                # actual edge nodes; a negatively declared cell carries them
                # along the border in reversed sequence. Linear cells carry
                # none.
                if edge.mid!==nothing
                    aligned=signed>0 ? edge.mid : reverse(edge.mid)
                    surface_borders.mids[parent][key]==aligned || throw(ArgumentError(
                        "$caller: declared Curve $curve edge node disagrees with its Surface border"))
                end
                key in covered && throw(ArgumentError(
                    "$caller: declared Surface boundary repeats a border Curve cell"))
                push!(covered,key)
            end
        end
        covered==Set(keys(border)) || throw(ArgumentError(
            "$caller: declared lower Surface Curves do not cover its complete border"))
    end
    Set(curve_parents)==used_curves || throw(ArgumentError(
        "$caller: unbound lower Curves are unsupported"))
    Set(point_parents)==used_points || throw(ArgumentError(
        "$caller: unbound lower Points are unsupported"))
    all(row->begin
            dim,tag=class.node_entities[row]
            dim==0 ? (get(point_cells,tag,nothing)!==nothing && point_cells[tag].rows[1]==row) :
            dim==1 ? (haskey(curve_nodes,tag) && row in curve_nodes[tag]) :
            dim==2 ? tag in surfaces : dim==3 && tag in parents
        end,
        unique(vcat([collect(cell.rows) for cell in vcat(faces,volumes,lines,points)]...))) ||
        throw(ArgumentError(
            "$caller: referenced lower BRep vertices have inconsistent entity classification"))
    facet_tables=Dict{Int,Vector{Tuple{Int,Tuple{Vararg{Int}}}}}()
    for msh in unique(volume.msh for volume in volumes)
        facet_tables[msh]=_partition_volume_facets(msh)
    end
    volume_facets=Vector{Vector{Tuple{Tuple{Vararg{UInt64}},NamedTuple}}}(undef,length(volumes))
    if workers>1
        Threads.@threads for index in eachindex(volumes)
            volume=volumes[index];entries=Tuple{Tuple{Vararg{UInt64}},NamedTuple}[]
            for (face_msh,pattern) in facet_tables[volume.msh]
                rows=ntuple(i->volume.rows[pattern[i]],length(pattern))
                labels=ntuple(i->volume.uids[pattern[i]],length(pattern))
                push!(entries,(_partition_support(labels),
                    (;rows,uids=labels,msh=face_msh,volume)))
            end
            volume_facets[index]=entries
        end
    else
        for index in eachindex(volumes)
            volume=volumes[index];entries=Tuple{Tuple{Vararg{UInt64}},NamedTuple}[]
            for (face_msh,pattern) in facet_tables[volume.msh]
                rows=ntuple(i->volume.rows[pattern[i]],length(pattern))
                labels=ntuple(i->volume.uids[pattern[i]],length(pattern))
                push!(entries,(_partition_support(labels),
                    (;rows,uids=labels,msh=face_msh,volume)))
            end
            volume_facets[index]=entries
        end
    end
    facets=Dict{Tuple{Vararg{UInt64}},Vector{NamedTuple}}()
    facet_order=Tuple{Vararg{UInt64}}[]
    for index in eachindex(volumes),(key,hit) in volume_facets[index]
        haskey(facets,key) || push!(facet_order,key)
        push!(get!(()->NamedTuple[],facets,key),hit)
    end
    for hits in values(facets)
        length(hits)<=2 || throw(ArgumentError("$caller: lower BRep facet has too many volume carriers"))
        length(hits)==2 && !_partition_cycle_reverse(hits[1].msh,hits[1].uids,
            hits[2].msh,hits[2].uids) &&
            throw(ArgumentError("$caller: shared volume facet has inconsistent winding"))
    end
    face_edges=Vector{Vector{Tuple{Tuple{UInt64,UInt64},NamedTuple}}}(undef,length(faces))
    if workers>1
        Threads.@threads for index in eachindex(faces)
            face=faces[index];corners,e=_partition_surface_shape(face.msh)
            entries=Tuple{Tuple{UInt64,UInt64},NamedTuple}[]
            for (order,(a,b)) in enumerate(_partition_edge_patterns(face.msh))
                mid=e==0 ? nothing : ntuple(i->(face.uids[corners+(order-1)*e+i],
                    face.rows[corners+(order-1)*e+i]),e)
                push!(entries,(minmax(face.uids[a],face.uids[b]),
                    (uids=(face.uids[a],face.uids[b]),rows=(face.rows[a],face.rows[b]),face,mid)))
            end
            face_edges[index]=entries
        end
    else
        for index in eachindex(faces)
            face=faces[index];corners,e=_partition_surface_shape(face.msh)
            entries=Tuple{Tuple{UInt64,UInt64},NamedTuple}[]
            for (order,(a,b)) in enumerate(_partition_edge_patterns(face.msh))
                mid=e==0 ? nothing : ntuple(i->(face.uids[corners+(order-1)*e+i],
                    face.rows[corners+(order-1)*e+i]),e)
                push!(entries,(minmax(face.uids[a],face.uids[b]),
                    (uids=(face.uids[a],face.uids[b]),rows=(face.rows[a],face.rows[b]),face,mid)))
            end
            face_edges[index]=entries
        end
    end
    face_support=Dict{Tuple{Vararg{UInt64}},NamedTuple}()
    lower_edges=Dict{Tuple{UInt64,UInt64},Vector{NamedTuple}}()
    for (face,entries) in zip(faces,face_edges)
        key=_partition_support(face.uids)
        haskey(face_support,key) && throw(ArgumentError(
            "$caller: original BRep repeats a Surface cell support"))
        face_support[key]=face
        hits=get(facets,key,NamedTuple[])
        !isempty(hits) && any(hit->hit.volume.part==face.part,hits) || throw(ArgumentError(
            "$caller: each lower Surface cell partition must match an actual volume facet carrier"))
        for (edge_key,entry) in entries
            push!(get!(()->NamedTuple[],lower_edges,edge_key),entry)
        end
    end
    declared=Dict{Int32,Vector{Int32}}();used=Set{Int32}()
    for parent in parents
        boundary=_partition_declared_boundary(model,3,Int(parent),caller)
        !isempty(boundary) && all(entry->entry[1]==2 && Int32(abs(entry[2])) in surfaces,boundary) &&
            allunique(abs(entry[2]) for entry in boundary) || throw(ArgumentError(
                "$caller: each Region needs a complete declared Surface boundary"))
        haskey(model.discrete,(3,Int(parent))) ||
            (boundary=NTuple{2,Int}[(d,_partition_region_surface_sign(
                parent,Int32(abs(signed)),signed,faces,facets))
                for (d,signed) in boundary])
        expected=Set(key for (key,hits) in facets if
            count(hit->hit.volume.parent==parent,hits)==1)
        covered=Set{Tuple{Vararg{UInt64}}}();declared[parent]=Int32[value for (_,value) in boundary]
        for (_,signed) in boundary
            surface=Int32(abs(signed));push!(used,surface)
            for face in faces
                face.parent==surface || continue
                key=_partition_support(face.uids);key in expected && !(key in covered) ||
                    throw(ArgumentError("$caller: declared Surface does not uniquely cover its Region facets"))
                hit=only(filter(hit->hit.volume.parent==parent,facets[key]))
                oriented=signed>0 ? face.uids : _partition_reversed(face.msh,face.uids)
                _partition_cycle_equal(face.msh,oriented,hit.uids) || throw(ArgumentError(
                    "$caller: declared Surface sign disagrees with its actual volume facet winding"))
                push!(covered,key)
            end
        end
        covered==expected || throw(ArgumentError("$caller: declared lower BRep omits Region facets"))
    end
    used==Set(surfaces) || throw(ArgumentError("$caller: unbound lower Surfaces are unsupported"))
    return (;entities,blocks,owners,uids,public,faces,volumes,lines,points,
        parents,surfaces,curve_parents,point_parents,curve_edges,point_cells,
        facets,facet_order,face_support,lower_edges,declared,surface_borders,
        curve_declared)
end

# Directed-traversal balance on one shared edge: an even family partitions
# into mutually reversed sheet pairs exactly when every direction count
# equals its reverse's — the local condition for a union of manifold sheets.
# A partition interface may therefore carry a junction edge several sheets
# cross through (valence 4, 6, …) while odd or lopsided incidence stays a
# genuine nonmanifold defect.
function _partition_edge_balanced(directions)
    counts=Dict{Tuple{UInt64,UInt64},Int}()
    for d in directions
        counts[d]=get(counts,d,0)+1
    end
    return all(d->counts[d]==get(counts,reverse(d),0),keys(counts))
end

function _partition_brep_boundary_certificate(source,children,part_entities,boundaries,
                                               interfaces,curves,caller;topology)
    surface_faces=Dict{Int32,Vector{NamedTuple}}()
    for face in source.faces
        child=children[(2,face.parent,face.part)]
        push!(get!(()->NamedTuple[],surface_faces,child),face)
    end
    merge!(surface_faces,interfaces)
    for ((dim,tag),_) in part_entities
        dim==3 || continue
        expected=Dict{Tuple{Vararg{UInt64}},NamedTuple}()
        for (key,hits) in source.facets
            selected=[hit for hit in hits if children[(3,hit.volume.parent,hit.volume.part)]==tag]
            length(selected)==1 && (expected[key]=only(selected))
        end
        covered=Set{Tuple{Vararg{UInt64}}}()
        for signed in boundaries[(dim,tag)],face in surface_faces[abs(signed)]
            key=_partition_support(face.uids)
            haskey(expected,key) && !(key in covered) || throw(ArgumentError(
                "$caller: child Region boundary repeats or includes a foreign volume facet"))
            oriented=signed>0 ? face.uids : _partition_reversed(face.msh,face.uids)
            _partition_cycle_equal(face.msh,oriented,expected[key].uids) || throw(ArgumentError(
                "$caller: child Region boundary has inconsistent facet winding"))
            push!(covered,key)
        end
        topology && covered!=Set(keys(expected)) && throw(ArgumentError(
            "$caller: child Region topology does not cover its complete volume boundary"))
    end
    topology || return nothing
    for (tag,faces) in surface_faces
        edges=Dict{Tuple{UInt64,UInt64},Vector{Tuple{UInt64,UInt64}}}()
        for face in faces,(a,b) in _partition_edge_patterns(face.msh)
            pair=(face.uids[a],face.uids[b]);push!(get!(()->Tuple{UInt64,UInt64}[],edges,minmax(pair...)),pair)
        end
        expected=Dict{Tuple{UInt64,UInt64},Tuple{UInt64,UInt64}}()
        for (key,hits) in edges
            if length(hits)==1
                expected[key]=only(hits);continue
            end
            _partition_edge_balanced(hits) || throw(ArgumentError(length(hits)==2 ?
                "$caller: child Surface internal edge winding disagrees" :
                "$caller: child Surface has a nonmanifold edge"))
        end
        covered=Set{Tuple{UInt64,UInt64}}()
        for signed in boundaries[(2,Int(tag))],edge in curves[abs(signed)]
            key=minmax(edge.uids...)
            # A declared Curve can span sibling children; only its edges on
            # this child's actual border bind here and must not repeat.
            haskey(expected,key) || continue
            oriented=signed>0 ? edge.uids : reverse(edge.uids)
            expected[key]==oriented && !(key in covered) || throw(ArgumentError(
                "$caller: Surface Curve incidence repeats or disagrees with its actual border"))
            push!(covered,key)
        end
        parent=part_entities[(2,Int(tag))].parent
        inherited=parent[1]==2 ? source.surface_borders.pairs[parent[2]] :
            Dict{Tuple{UInt64,UInt64},Tuple{UInt64,UInt64}}()
        all(key->key in covered || get(inherited,key,nothing)==expected[key],keys(expected)) ||
            throw(ArgumentError("$caller: child Surface omits a new seam or changes its original border"))
    end
    return nothing
end

# Seam-edge decomposition into Curve chains. Vertices of degree two pass a
# chain through; every other vertex — an open endpoint (degree 1) or a
# junction where three or more seam edges meet — terminates one. A junction
# vertex therefore splits the incident network into the maximal degree-2
# arcs its Point15 entity then binds, instead of rejecting the fan, while a
# pure cycle of degree-2 vertices still emits one closed chain.
function _partition_brep_chains(edges,caller)
    adjacent=Dict{UInt64,Vector{Int}}()
    for (index,edge) in enumerate(edges),uid in edge.uids
        push!(get!(()->Int[],adjacent,uid),index)
    end
    manifold(uid)=length(get(adjacent,uid,Int[]))==2
    # Edge directions are arbitrary — a junction arc can mix carrier and
    # interface-derived orientations — so each chain also records whether an
    # edge is traversed against its literal direction.
    starts=[index for (index,edge) in enumerate(edges)
        if !manifold(edge.uids[1]) || !manifold(edge.uids[2])]
    seen=falses(length(edges));chains=NamedTuple[]
    for seed in vcat(starts,collect(eachindex(edges)))
        seen[seed] && continue
        edge=edges[seed]
        entry=manifold(edge.uids[1]) && !manifold(edge.uids[2]) ?
            edge.uids[2] : edge.uids[1]
        indices=Int[];flips=Bool[];index=seed;vertex=entry
        while true
            seen[index]=true;current=edges[index];forward=current.uids[1]==vertex
            push!(indices,index);push!(flips,!forward)
            vertex=forward ? current.uids[2] : current.uids[1]
            manifold(vertex) || break
            nexts=filter(!=(index),adjacent[vertex]);next=only(nexts)
            if seen[next]
                next==seed || throw(ArgumentError(
                    "$caller: partition seam joins another chain"))
                break
            end
            index=next
        end
        closed=vertex==entry
        closed && length(indices)<3 && throw(ArgumentError("$caller: partition seam is not a closed polygon"))
        push!(chains,(indices=indices,flips=flips,closed=closed))
    end
    return chains
end

# Face-adjacency links for one partition-interface pair. Manifold edges
# (exactly two incident faces) link their pair outright; a junction edge
# with three or more incident faces is split into sheets: two faces link
# only when they are the two facets of the *same* incident volume carrier —
# the unique pairing whose co-oriented normals traverse the shared edge
# oppositely after the interface orientation is unified. Each face carries
# at most one partner per edge (canonical side first, then the other), so
# every emitted component stays edge-manifold for the junction-fan
# certificate; faces left unmatched simply bound their sheet there.
function _partition_interface_links(faces,pair,children,facets)
    supports=Dict{Tuple{Vararg{UInt64}},Int}()
    incidence=Dict{Tuple{UInt64,UInt64},Vector{Int}}()
    for (index,face) in enumerate(faces)
        supports[_partition_support(face.uids)]=index
        for (a,b) in _partition_edge_patterns(face.msh)
            push!(get!(()->Int[],incidence,minmax(face.uids[a],face.uids[b])),index)
        end
    end
    links=Tuple{Int,Int}[]
    for (key,hits) in incidence
        if length(hits)==2
            push!(links,(hits[1],hits[2]));continue
        end
        length(hits)>2 || continue
        # Candidate partner per face: the face the same incident carrier
        # pairs it with on this edge — searched on each side in `pair`
        # order. Two faces only ever share one carrier's corner here, so
        # the first found candidate already defines a mutual link; a face
        # keeps at most one partner.
        partners=Dict{Int,Int}()
        for i in hits
            haskey(partners,i) && continue
            face=faces[i];fsup=_partition_support(face.uids)
            fhits=get(facets,fsup,NamedTuple[])
            for side in pair
                carrier=findfirst(hit->children[(3,hit.volume.parent,hit.volume.part)]==side,
                                  fhits)
                carrier===nothing && continue
                volume=fhits[carrier].volume;partner=0
                for (_,pattern) in _partition_volume_facets(volume.msh)
                    fu=ntuple(slot->volume.uids[pattern[slot]],length(pattern))
                    (key[1] in fu && key[2] in fu) || continue
                    candidate=get(supports,_partition_support(fu),nothing)
                    candidate===nothing && continue
                    (candidate==i || haskey(partners,candidate)) && continue
                    partner=candidate;break
                end
                partner==0 && continue
                partners[i]=partner;partners[partner]=i
                break
            end
        end
        for (i,j) in partners
            i<j && push!(links,(i,j))
        end
    end
    return links
end

function _partition_brep_production_plan(model,cached,class,count,assignment,caller;
                                         topology::Bool,ghosts::Bool)
    source=_partition_brep_snapshot(model,cached,class,assignment,caller)
    volume_parts=Int32[volume.part for volume in source.volumes]
    length(unique(volume_parts))==count || throw(ArgumentError(
        "$caller: every partition must own a volume cell"))
    # A declared Surface cell on a cross-part facet already carries the
    # interface role: its partition child keeps the cells, binds both Region
    # children through the declared boundary with opposite signs and only
    # extends its partition membership to the incident set. No duplicate
    # interface entity or cell is emitted for that support.
    coincident=Set{Tuple{Vararg{UInt64}}}()
    for key in source.facet_order
        hits=source.facets[key]
        topology && length(hits)==2 && hits[1].volume.part!=hits[2].volume.part &&
            haskey(source.face_support,key) && push!(coincident,key)
    end
    staged=Model._model_mesh_transaction_copy(model)
    part_entities=Dict{Tuple{Int,Int},Elements.MixedEntity}()
    children=Dict{Tuple{Int,Int32,Int32},Int32}()
    interfaces=Dict{Int32,Vector{NamedTuple}}();curves=Dict{Int32,Vector{NamedTuple}}()
    for dim in (2,3)
        cells=dim==2 ? source.faces : source.volumes
        for (parent,part) in sort!(unique((cell.parent,cell.part) for cell in cells))
            tag=Int32(_partition_entity_birth!(staged,dim,caller));children[(dim,parent,part)]=tag
            groups=Model._model_projection_physical_tags(staged,dim,Int(parent))
            staged.entity_physicals[(dim,Int(tag))]=copy(groups)
            for group in groups
                members=staged.physical[(dim,Int(group))];Int(tag) in members || push!(members,Int(tag))
            end
            selected=[cell for cell in cells if cell.parent==parent && cell.part==part]
            rows=unique(vcat([collect(cell.rows) for cell in selected]...))
            part_entities[(dim,Int(tag))]=Elements.MixedEntity(dim,tag,_partition_bounds(cached.coords,rows);
                physical_tags=groups,partitions=Int32[part],parent=(dim,parent))
        end
    end
    curve_children=Dict{Int32,Int32}();point_children=Dict{Int32,Int32}()
    curve_partitions=Dict{Int32,Vector{Int32}}()
    for parent in source.curve_parents
        tag=Int32(_partition_entity_birth!(staged,1,caller));curve_children[parent]=tag
        groups=Model._model_projection_physical_tags(staged,1,Int(parent))
        staged.entity_physicals[(1,Int(tag))]=copy(groups)
        for group in groups
            members=staged.physical[(1,Int(group))];Int(tag) in members || push!(members,Int(tag))
        end
        edges=source.curve_edges[parent]
        parts=sort!(unique(Int32[entry.face.part for edge in edges
            for entry in source.lower_edges[minmax(edge.uids...)]]))
        curve_partitions[parent]=parts
        rows=unique(reduce(vcat,[collect(edge.rows) for edge in edges]))
        part_entities[(1,Int(tag))]=Elements.MixedEntity(1,tag,_partition_bounds(cached.coords,rows);
            physical_tags=groups,partitions=parts,parent=(1,Int(parent)))
        curves[tag]=edges
    end
    point_declaring=Dict{Int32,Vector{Int32}}()
    for parent in source.curve_parents,(_,signed) in source.curve_declared[parent]
        push!(get!(()->Int32[],point_declaring,Int32(abs(signed))),parent)
    end
    for parent in source.point_parents
        tag=Int32(_partition_entity_birth!(staged,0,caller));point_children[parent]=tag
        groups=Model._model_projection_physical_tags(staged,0,Int(parent))
        staged.entity_physicals[(0,Int(tag))]=copy(groups)
        for group in groups
            members=staged.physical[(0,Int(group))];Int(tag) in members || push!(members,Int(tag))
        end
        cell=source.point_cells[parent]
        parts=sort!(unique(reduce(vcat,[curve_partitions[curve] for curve in point_declaring[parent]])))
        part_entities[(0,Int(tag))]=Elements.MixedEntity(0,tag,_partition_bounds(cached.coords,[cell.rows[1]]);
            physical_tags=groups,partitions=parts,parent=(0,Int(parent)))
    end
    flat_owners=Int32[owner for owners in source.owners for owner in owners]
    for cell in source.faces;flat_owners[cell.index]=children[(2,cell.parent,cell.part)];end
    for cell in source.volumes;flat_owners[cell.index]=children[(3,cell.parent,cell.part)];end
    for cell in source.lines;flat_owners[cell.index]=curve_children[cell.parent];end
    for cell in source.points;flat_owners[cell.index]=point_children[cell.parent];end
    blocks=copy(source.blocks);owners=Vector{Int32}[];offset=0
    for block in blocks
        n=size(block.nodes,2);push!(owners,flat_owners[offset+1:offset+n]);offset+=n
    end
    labels=copy(source.uids);cell_count=length(source.public.element_tags)
    maximum_uid=max(ELEMENT_TAG_MAX[],maximum(source.public.element_tags;init=UInt64(0)))
    boundaries=Dict(key=>Int32[] for key in keys(part_entities))
    for volume in source.volumes
        child=children[(3,volume.parent,volume.part)]
        for signed in source.declared[volume.parent],face in source.faces
            face.parent==abs(signed) || continue
            hits=source.facets[_partition_support(face.uids)]
            any(hit->hit.volume.index==volume.index,hits) || continue
            face_child=children[(2,face.parent,face.part)]*sign(signed)
            face_child in boundaries[(3,Int(child))] ||
                push!(boundaries[(3,Int(child))],face_child)
        end
    end
    for parent in source.curve_parents
        boundaries[(1,Int(curve_children[parent]))]=Int32[
            sign(signed)*point_children[Int32(abs(signed))]
            for (_,signed) in source.curve_declared[parent]]
    end
    seams=Dict{Tuple{Tuple{Int,Int32},Tuple{Vararg{Int32}}},Vector{NamedTuple}}()
    seam_order=Tuple{Tuple{Int,Int32},Tuple{Vararg{Int32}}}[]
    if topology
        edge_curve=Dict{Tuple{UInt64,UInt64},NamedTuple}()
        for parent in source.curve_parents,edge in source.curve_edges[parent]
            edge_curve[minmax(edge.uids...)]=(uids=edge.uids,tag=curve_children[parent])
        end
        seam_lookup=Dict{Tuple{UInt64,UInt64},NamedTuple}()
        for (key,entry) in edge_curve
            seam_lookup[key]=entry
        end
        grouped=Dict{Tuple{Int32,Int32},Vector{NamedTuple}}()
        for key in source.facet_order
            hits=source.facets[key];length(hits)==2 || continue
            left,right=hits;left.volume.part==right.volume.part && continue
            key in coincident && continue
            child=children[(3,left.volume.parent,left.volume.part)]
            other=children[(3,right.volume.parent,right.volume.part)]
            push!(get!(()->NamedTuple[],grouped,minmax(child,other)),
                (rows=left.rows,uids=left.uids,msh=left.msh,child=child,other=other))
        end
        for key in coincident
            face=source.face_support[key];hits=source.facets[key]
            tag=children[(2,face.parent,face.part)];entity=part_entities[(2,Int(tag))]
            merged=sort!(unique!(vcat(entity.partitions,Int32[hit.volume.part for hit in hits])))
            merged!=entity.partitions &&
                (part_entities[(2,Int(tag))]=Elements.MixedEntity(2,tag,entity.bbox;
                    physical_tags=entity.physical_tags,partitions=merged,parent=entity.parent))
        end
        # Every interface cell edge is indexed globally, across components:
        # an edge internal to one component stays silent, while any other
        # edge must close a junction fan — every volume incident on the edge
        # contributes exactly two of its own facets (interfaces, same-part
        # interior facets, or lower Surface carriers).
        global_edges=Dict{Tuple{UInt64,UInt64},Vector{NamedTuple}}()
        global_order=Tuple{UInt64,UInt64}[]
        edge_volumes=Dict{Tuple{UInt64,UInt64},Vector{NamedTuple}}()
        for key in source.facet_order,hits in (source.facets[key],)
            for (a,b) in _partition_edge_patterns(hits[1].msh),hit in hits
                push!(get!(()->NamedTuple[],edge_volumes,minmax(hit.uids[a],hit.uids[b])),hit.volume)
            end
        end
        for pair in sort!(collect(keys(grouped)))
            faces=grouped[pair]
            links=_partition_interface_links(faces,pair,children,source.facets)
            for component in _partition_interface_components(faces,links)
                primary=faces[first(component)].child;tag=Int32(_partition_entity_birth!(staged,2,caller))
                cell_count<=typemax(Int32)-length(component) || throw(ArgumentError("$caller: interface cell capacity exceeded"))
                cell_count+=length(component)
                members=sort!(Int32[part_entities[(3,Int(child))].partitions[1] for child in pair])
                interface_faces=NamedTuple[]
                oriented_rows=Vector{Tuple{Vararg{Int32}}}(undef,length(component))
                ids=UInt64[]
                for (slot,index) in enumerate(component)
                    face=faces[index]
                    rows=face.child==primary ? face.rows : _partition_reversed(face.msh,face.rows)
                    uids=Tuple(source.public.node_tags[collect(rows)])
                    oriented_rows[slot]=rows
                    maximum_uid<typemax(UInt64) || throw(ArgumentError("$caller: interface UID capacity exceeded"))
                    maximum_uid+=UInt64(1);push!(ids,maximum_uid)
                    push!(interface_faces,(;rows,uids,msh=face.msh))
                end
                rows_flat=Int32[row for rows in oriented_rows for row in rows]
                part_entities[(2,Int(tag))]=Elements.MixedEntity(2,tag,_partition_bounds(cached.coords,unique(rows_flat));
                    partitions=members,parent=part_entities[(3,Int(primary))].parent)
                boundaries[(2,Int(tag))]=Int32[];interfaces[tag]=interface_faces
                for face_msh in sort!(unique(face.msh for face in interface_faces))
                    slots=findall(face->face.msh==face_msh,interface_faces)
                    nodes=Matrix{Int32}(undef,length(oriented_rows[first(slots)]),length(slots))
                    for (column,slot) in enumerate(slots)
                        nodes[:,column]=collect(oriented_rows[slot])
                    end
                    push!(blocks,ElementBlock(Int(face_msh),nodes))
                    push!(owners,fill(tag,length(slots)));push!(labels,ids[slots])
                end
                for child in pair;push!(boundaries[(3,Int(child))],child==primary ? tag : -tag);end
                edge_hits=Dict{Tuple{UInt64,UInt64},Vector{NamedTuple}}()
                for face in interface_faces,(order,(a,b)) in enumerate(_partition_edge_patterns(face.msh))
                    key=minmax(face.uids[a],face.uids[b])
                    corners,e=_partition_surface_shape(face.msh)
                    mid=e==0 ? nothing : ntuple(i->(face.uids[corners+(order-1)*e+i],
                        face.rows[corners+(order-1)*e+i]),e)
                    push!(get!(()->NamedTuple[],edge_hits,key),
                        (uids=(face.uids[a],face.uids[b]),rows=(face.rows[a],face.rows[b]),mid))
                    haskey(global_edges,key) || push!(global_order,key)
                    push!(get!(()->NamedTuple[],global_edges,key),
                        (uids=(face.uids[a],face.uids[b]),rows=(face.rows[a],face.rows[b]),mid,
                         interface=tag))
                end
                for (key,hits) in edge_hits
                    length(hits)==1 && continue
                    _partition_edge_balanced(Tuple{UInt64,UInt64}[hit.uids for hit in hits]) ||
                        throw(ArgumentError(length(hits)==2 ?
                            "$caller: interface edge winding disagrees" :
                            "$caller: nonmanifold interface edge"))
                    # Sheet pairs traverse the shared edge in mutually
                    # reversed directions, so normalized to corner order
                    # every sheet must carry the identical node sequence.
                    norm(hit)=hit.mid===nothing ? nothing :
                        (hit.uids[1]<hit.uids[2] ? Tuple(hit.mid) : reverse(Tuple(hit.mid)))
                    allequal(norm(hit) for hit in hits) ||
                        throw(ArgumentError("$caller: interface edge node disagrees"))
                end
            end
        end
        for key in global_order
            hits=global_edges[key]
            length(hits)==2 && allequal(hit.interface for hit in hits) && continue
            lower=get(source.lower_edges,key,NamedTuple[])
            volume_cells=Dict{NamedTuple,Int}()
            for volume in get(edge_volumes,key,NamedTuple[])
                volume_cells[volume]=get(volume_cells,volume,0)+1
            end
            all(==(2),values(volume_cells)) || throw(ArgumentError(
                "$caller: interface seam must meet actual lower Surface carriers closing its junction fan"))
            # Edge-node sequences follow each incident cell's own traversal
            # direction; align them to the canonical corner order so cells
            # covering the edge in opposite directions compare consistently.
            norm_mid=(uids,mid)->mid===nothing ? nothing :
                (uids[1]<uids[2] ? Tuple(mid) : reverse(Tuple(mid)))
            mids=Union{Nothing,Tuple{Vararg{Tuple{UInt64,Int32}}}}[
                norm_mid(hit.uids,hit.mid) for hit in hits]
            append!(mids,Union{Nothing,Tuple{Vararg{Tuple{UInt64,Int32}}}}[
                norm_mid(edge.uids,edge.mid) for edge in lower])
            allequal(mids) || throw(ArgumentError(
                "$caller: interface seam edge node disagrees across its incident cells"))
            parts=sort!(unique!(Int32[volume.part for volume in keys(volume_cells)]))
            curve_edge=get(edge_curve,key,nothing)
            if curve_edge!==nothing
                entity=part_entities[(1,Int(curve_edge.tag))]
                merged=sort!(unique!(vcat(entity.partitions,parts)))
                merged!=entity.partitions &&
                    (part_entities[(1,Int(curve_edge.tag))]=Elements.MixedEntity(1,curve_edge.tag,entity.bbox;
                        physical_tags=entity.physical_tags,partitions=merged,parent=entity.parent))
                seam_lookup[key]=(uids=curve_edge.uids,tag=curve_edge.tag)
                continue
            end
            if isempty(lower)
                parent=(3,Int32(minimum(Int(volume.parent) for volume in keys(volume_cells))))
                edge=(;uids=hits[1].uids,rows=hits[1].rows,mid=hits[1].mid)
            else
                reference=argmin([children[(2,edge.face.parent,edge.face.part)] for edge in lower])
                edge=lower[reference];parent=(2,edge.face.parent)
            end
            group=(parent,Tuple(parts))
            haskey(seams,group) || push!(seam_order,group)
            current=get!(()->NamedTuple[],seams,group)
            any(existing->_partition_support(existing.uids)==_partition_support(edge.uids),current) ||
                push!(current,(uids=edge.uids,rows=edge.rows,mid=edge.mid))
        end
        # Junction vertices — nodes where three or more seam edges meet —
        # materialize as dim-0 entities carrying Point15 cells; their parent is
        # the lowest-dimension parent among incident seam groups. A vertex
        # already bound by a declared Point child reuses that child instead.
        vertex_degree=Dict{UInt64,Int}()
        vertex_groups=Dict{UInt64,Set{Tuple{Tuple{Int,Int32},Tuple{Vararg{Int32}}}}}()
        for group in seam_order,edge in seams[group],uid in edge.uids
            vertex_degree[uid]=get(vertex_degree,uid,0)+1
            push!(get!(()->Set{Tuple{Tuple{Int,Int32},Tuple{Vararg{Int32}}}}(),vertex_groups,uid),group)
        end
        for (_,edge) in edge_curve,uid in edge.uids
            vertex_degree[uid]=get(vertex_degree,uid,0)+1
        end
        point_nodes=Dict{UInt64,Int32}()
        for (parent,child) in point_children
            cell=source.point_cells[parent];uid=source.public.node_tags[cell.rows[1]]
            point_nodes[uid]=min(get(point_nodes,uid,Int32(typemax(Int32))),child)
        end
        junctions=Dict{UInt64,Int32}()
        for uid in sort!(collect(uid for (uid,degree) in vertex_degree if degree>=3))
            groups=get(vertex_groups,uid,nothing)
            groups===nothing && continue
            parts=sort!(unique!(Int32[part for group in collect(groups) for part in group[2]]))
            existing=get(point_nodes,uid,nothing)
            if existing!==nothing
                entity=part_entities[(0,Int(existing))]
                merged=sort!(unique!(vcat(entity.partitions,parts)))
                part_entities[(0,Int(existing))]=Elements.MixedEntity(0,existing,entity.bbox;
                    physical_tags=entity.physical_tags,partitions=merged,parent=entity.parent)
                junctions[uid]=existing;continue
            end
            parent=minimum(group[1] for group in groups)
            tag=Int32(_partition_entity_birth!(staged,0,caller))
            cell_count<=typemax(Int32) || throw(ArgumentError("$caller: junction cell capacity exceeded"))
            cell_count+=1
            maximum_uid<typemax(UInt64) || throw(ArgumentError("$caller: junction UID capacity exceeded"))
            maximum_uid+=UInt64(1)
            row=source.public.node_indices[uid]
            part_entities[(0,Int(tag))]=Elements.MixedEntity(0,tag,_partition_bounds(cached.coords,Int32[row]);
                partitions=parts,parent=parent)
            boundaries[(0,Int(tag))]=Int32[]
            junctions[uid]=tag
            push!(blocks,ElementBlock(15,reshape(Int32[row],1,1)));push!(owners,[tag]);push!(labels,[maximum_uid])
        end
        for group in seam_order
            edges=seams[group]
            for chain in _partition_brep_chains(edges,caller)
                oriented=map(zip(chain.indices,chain.flips)) do (index,flip)
                    edge=edges[index]
                    flip ? (uids=(edge.uids[2],edge.uids[1]),
                            rows=(edge.rows[2],edge.rows[1]),
                            mid=edge.mid===nothing ? nothing : reverse(edge.mid)) : edge
                end
                cell_count<=typemax(Int32)-length(oriented) || throw(ArgumentError("$caller: seam cell capacity exceeded"))
                cell_count+=length(oriented)
                arity=all(edge->edge.mid===nothing,oriented) ? 2 :
                    2+length(oriented[findfirst(edge->edge.mid!==nothing,oriented)].mid)
                all(edge->(edge.mid===nothing)==(arity==2) &&
                    (edge.mid===nothing || length(edge.mid)==arity-2),oriented) ||
                    throw(ArgumentError("$caller: partition seam mixes cell orders"))
                tag=Int32(_partition_entity_birth!(staged,1,caller));nodes=Matrix{Int32}(undef,arity,length(oriented));ids=UInt64[]
                for (column,edge) in enumerate(oriented)
                    nodes[:,column]=arity==2 ? collect(edge.rows) :
                        Int32[edge.rows[1],edge.rows[2],(edge.mid[i][2] for i in 1:arity-2)...]
                    maximum_uid<typemax(UInt64) || throw(ArgumentError("$caller: seam UID capacity exceeded"))
                    maximum_uid+=UInt64(1);push!(ids,maximum_uid)
                    seam_lookup[minmax(edge.uids...)]=(uids=edge.uids,tag=tag)
                end
                part_entities[(1,Int(tag))]=Elements.MixedEntity(1,tag,_partition_bounds(cached.coords,unique(vec(nodes)));
                    partitions=collect(group[2]),parent=group[1])
                bounds=Int32[]
                if !chain.closed
                    start=get(junctions,oriented[1].uids[1],Int32(0))
                    stop=get(junctions,oriented[end].uids[2],Int32(0))
                    start!=0 && push!(bounds,start)
                    stop!=0 && push!(bounds,-stop)
                end
                boundaries[(1,Int(tag))]=bounds
                curves[tag]=oriented
                line_type=Elements._msh_type_exact(:lin,arity-1,false,caller)
                push!(blocks,ElementBlock(line_type,nodes));push!(owners,fill(tag,length(ids)));push!(labels,ids)
            end
        end
        for ((dim,tag),entity) in part_entities
            dim==2 || continue
            faces=haskey(interfaces,Int32(tag)) ? interfaces[Int32(tag)] :
                [face for face in source.faces if children[(2,face.parent,face.part)]==tag]
            # A junction edge several sheets of one Surface cross through is
            # traversed in both directions; those canceling incidences are
            # internal to the child and vote no boundary sign, while every
            # remaining border incidence must still agree on one sign.
            tally=Dict{Tuple{Int32,Tuple{UInt64,UInt64}},Int}()
            for face in faces,(a,b) in _partition_edge_patterns(face.msh)
                key=minmax(face.uids[a],face.uids[b]);seam=get(seam_lookup,key,nothing);seam===nothing && continue
                vote=(face.uids[a],face.uids[b])==seam.uids ? 1 : -1
                slot=(seam.tag,key);tally[slot]=get(tally,slot,0)+vote
            end
            signs=Dict{Int32,Int}()
            for ((tag,key),net) in tally
                net==0 && continue
                sign=net>0 ? 1 : -1
                haskey(signs,tag) && signs[tag]!=sign && throw(ArgumentError("$caller: Surface seam has inconsistent orientation"))
                signs[tag]=sign
            end
            boundaries[(dim,tag)]=Int32[Int32(signs[curve])*curve for curve in sort!(collect(keys(signs)))]
        end
    end
    _partition_brep_boundary_certificate(source,children,part_entities,boundaries,interfaces,curves,caller;topology)
    node_owners=copy(class.node_entities)
    sentinel=(4,Int32(typemax(Int32)))
    # Per-block best (lowest-dimension incident entity, parametric-eligible
    # incident entity) merge by `isless`, so the result is order-independent.
    partials=Vector{Dict{Int32,Tuple{Tuple{Int,Int32},Tuple{Int,Int32}}}}(undef,length(blocks))
    incident_block=function (b)
        block=blocks[b];classified=owners[b];dim=Elements.msh_dimension(block.msh)
        d=Dict{Int32,Tuple{Tuple{Int,Int32},Tuple{Int,Int32}}}()
        for column in axes(block.nodes,2),row in block.nodes[:,column]
            candidate=(dim,classified[column])
            old=get(d,row,(sentinel,sentinel))
            best=isless(candidate,old[1]) ? candidate : old[1]
            parbest=old[2]
            parameters=class.node_parametric[row]
            if parameters!==nothing && candidate[1]==length(parameters) &&
                    isless(candidate,old[2])
                entry=get(part_entities,(candidate[1],Int(candidate[2])),nothing)
                owner=(class.node_entities[row][1],Int(class.node_entities[row][2]))
                entry!==nothing && entry.parent==owner && (parbest=candidate)
            end
            d[row]=(best,parbest)
        end
        return d
    end
    workers=_partition_workers()
    if workers>1
        Threads.@threads for b in eachindex(blocks);partials[b]=incident_block(b);end
    else
        for b in eachindex(blocks);partials[b]=incident_block(b);end
    end
    incident=fill(sentinel,nnodes(cached))
    incident_par=fill(sentinel,nnodes(cached))
    for d in partials,(row,(c,pc)) in pairs(d)
        isless(c,incident[row]) && (incident[row]=c)
        isless(pc,incident_par[row]) && (incident_par[row]=pc)
    end
    for row in eachindex(incident)
        pick=class.node_parametric[row]===nothing ? incident[row] : incident_par[row]
        pick[1]<4 && (node_owners[row]=pick)
    end
    for ((dim,tag),entity) in collect(part_entities)
        boundary=boundaries[(dim,tag)]
        staged.discrete[(dim,tag)].boundary=Tuple{Int,Int}[(dim-1,Int(value)) for value in boundary]
        part_entities[(dim,tag)]=Elements.MixedEntity(dim,tag,entity.bbox;
            physical_tags=entity.physical_tags,boundaries=boundary,partitions=entity.partitions,parent=entity.parent)
    end
    global_entities=copy(source.entities);ghost_entities=Tuple{Int32,Int32}[];ghost_elements=Elements.MixedGhostElement[]
    if ghosts
        node_parts=[Set{Int32}() for _ in 1:nnodes(cached)]
        for volume in source.volumes,row in volume.rows;push!(node_parts[row],volume.part);end
        destinations=Vector{Set{Int32}}(undef,length(source.volumes))
        if workers>1
            Threads.@threads for index in eachindex(source.volumes)
                parts=Set{Int32}()
                for row in source.volumes[index].rows;union!(parts,node_parts[row]);end
                destinations[index]=parts
            end
        else
            for index in eachindex(source.volumes)
                parts=Set{Int32}()
                for row in source.volumes[index].rows;union!(parts,node_parts[row]);end
                destinations[index]=parts
            end
        end
        for (index,volume) in enumerate(source.volumes)
            delete!(destinations[index],volume.part);isempty(destinations[index]) && continue
            push!(ghost_elements,Elements.MixedGhostElement(volume.index,volume.part,
                sort!(collect(destinations[index]))))
        end
        for part in 1:count
            tag=Int32(_partition_entity_birth!(staged,3,caller));push!(ghost_entities,(tag,Int32(part)))
            selected=[volume for (index,volume) in enumerate(source.volumes) if Int32(part) in destinations[index]]
            rows=isempty(selected) ? Int32[] : unique(vcat([collect(volume.rows) for volume in selected]...))
            global_entities[(3,Int(tag))]=Elements.MixedEntity(3,tag,_partition_bounds(cached.coords,rows))
        end
    end
    sum(size(block.nodes,2) for block in blocks)<=typemax(Int32) || throw(ArgumentError("$caller: partition cell capacity exceeded"))
    data=Elements._native_mixed_entity_data(global_entities;node_entities=node_owners,
        node_parametric=deepcopy(class.node_parametric),external_node_tags=copy(source.public.node_tags),
        block_entities=owners,external_element_tags=labels)
    partition=Elements.MixedPartitionData(count;entities=part_entities,ghost_entities,ghost_elements)
    target=MixedMesh(cached.coords,blocks;entity_data=data,elementary_entities=owners,
        physical_names=staged.physical_names,partition_data=partition,
        periodic_links=cached isa MixedMesh ? cached.periodic_links : Elements.MixedPeriodicLink[],
        ancillary_sections=cached isa MixedMesh ? cached.ancillary_sections : [],
        data_sections=cached isa MixedMesh ? cached.data_sections : Elements.MshDataSection[])
    validate(target).ok || throw(ArgumentError("$caller: lower BRep partition produced invalid or duplicate cells"))
    _geo_import_entity_receipts!(staged,target)
    identities=Dict(uid=>uid for uid in source.public.node_tags)
    declarations=union(Set(keys(global_entities)),Set(keys(part_entities)))
    _dim01_reconcile_records!(staged,declarations,target,identities;runtime_parameters=data.node_parametric)
    previous=model.meshing.mesh_authority
    ranks=previous!==nothing && previous.mesh===cached ? previous.private_rank : collect(UInt64,1:nnodes(target))
    staged.meshing.mesh_authority=Model._model_mesh_authority(target;caller,private_rank=ranks)
    cell_map=Dict(uid=>uid for uid in source.public.element_tags)
    return (;model=staged,mesh=target,class=_higher_classification_from_data(target,caller),
        authority=Set{Tuple{Int,Int32}}((dim,Int32(tag)) for (dim,tag) in declarations),
        max_node_tag=NODE_TAG_MAX[],max_element_tag=maximum_uid,
        node_tag_map=identities,element_tag_map=cell_map,visibility_tag_map=cell_map)
end
