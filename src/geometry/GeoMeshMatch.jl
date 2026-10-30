# ─── `Geometry.MatchGeomAndMesh` — `GeomMeshMatcher::match` parity ───────────
#
# Upstream (OpenFile.cpp + GeomMeshMatcher.cpp): when `matchGeomAndMesh` is set
# and the current model is non-empty, an opened `.msh` is read into a temporary
# model, scaled by `MatchMeshScaleFactor`, then matched onto the current
# model's entities. Vertex entities pair by position within
# `MatchMeshTolerance` (max-abs coordinate distance); curves pair through the
# matched endpoints they share (closed curves admit any closed candidate);
# faces pair through their matched boundary-curve sets; regions require all of
# the matched boundary faces. Multi-candidate ties break on an oriented-
# bounding-box distance. Point/curve counts compare against the *mesh* totals,
# surface counts are advisory only, and region counts compare against the
# *model* total — all mirrored here. On success the mesh is grafted onto the
# model entities: element data keeps the geometry's tags, physical names
# transfer wholesale, point/curve entities keep the geometry's physical
# memberships and face/region entities adopt the mesh file's. On failure
# nothing is copied.
#
# Tessella keeps the same observable contract natively: the imported actual
# mesh is re-keyed to the matched model entities, the previous mesh on those
# entities is replaced (upstream `deleteMesh` per matched entity), physical
# records and names follow upstream's direction rules, and every mutation
# happens only after a complete match — a failed match leaves the model,
# mesh state and receipts untouched.
#
# Bounded scope: `.msh` inputs, simplex/quad/prism/pyramid/hex cell families
# as imported data (matching needs only connectivity and positions), no
# partitioned inputs, and no geometric node reprojection — matched nodes keep
# their file coordinates (upstream projects edge vertices onto the CAD curve;
# genuine matches differ by less than the tolerance).

function _geo_match_active(m::GeoModel,context)
    on=something(_geo_option_number(context,"Geometry",0,"MatchGeomAndMesh"),0.0)
    iszero(on) && return false
    return !isempty(Model.model_entities(m))
end

# Canonical corner count — MSH records list corner nodes first.
function _geo_match_cell_corners(msh::Int)
    msh_spec(msh).dim==0 && return 1
    family=msh_family(msh)
    family===:lin && return 2
    family===:tri && return 3
    family===:quad && return 4
    family===:tet && return 4
    family===:hex && return 8
    family===:prism && return 6
    family===:pyr && return 5
    return 0
end

_geo_match_edge_key(a,b)=a<b ? (Int32(a),Int32(b)) : (Int32(b),Int32(a))

function _geo_match_pad4(rows)
    t=Tuple(sort!(Int32[r for r in rows]))
    return ntuple(i->i<=length(t) ? t[i] : Int32(0),4)
end

# Extreme nodes of an entity's 1-D cells: corner nodes whose seg degree
# differs from 2. Closed curves have none — the empty result is what the
# caller treats as "closed".
function _geo_match_curve_extremes(incoming,data,owner)
    degree=Dict{Int32,Int}();seen=false
    for (bi,block) in enumerate(incoming.blocks)
        block isa SpecialElementBlock && continue
        msh_dimension(block.msh)==1 || continue
        n=_geo_match_cell_corners(block.msh);n<2 && continue
        for c in axes(block.nodes,2)
            Int(data.block_entities[bi][c])==owner || continue
            seen=true
            for k in 1:2
                r=block.nodes[k,c];degree[r]=get(degree,r,0)+1
            end
        end
    end
    seen || return Int32[]
    return Int32[r for (r,d) in degree if d!=2]
end

# Incoming-entity boundary tag sets in mesh-file numbering. Declared
# `$Entities` boundary records win; records-less entities fall back to
# connectivity incidence — curves by seg-extreme nodes owned by point
# entities, faces by any shared cell edge, volumes by any shared cell face.
function _geo_match_inboundaries(incoming,data)
    entities=data.entities
    cells=Dict{Tuple{Int,Int},Vector{NTuple{2,Int}}}()
    for (bi,block) in enumerate(incoming.blocks)
        block isa SpecialElementBlock && continue
        dim=msh_dimension(block.msh)
        for c in axes(block.nodes,2)
            key=(dim,Int(data.block_entities[bi][c]))
            push!(get!(()->NTuple{2,Int}[],cells,key),(bi,c))
        end
    end
    bounds=Dict{Tuple{Int,Int},Set{Int32}}()
    for (key,entity) in entities
        dim,tag=key
        dim==0 && continue
        if !isempty(entity.boundaries)
            bounds[key]=Set{Int32}(abs.(entity.boundaries));continue
        end
        out=Set{Int32}()
        if dim==1
            for r in _geo_match_curve_extremes(incoming,data,tag)
                owner=data.node_entities[r]
                owner[1]==0 && push!(out,Int32(owner[2]))
            end
        elseif dim==2
            edges=Set{NTuple{2,Int32}}()
            for (bi,c) in get(cells,key,NTuple{2,Int}[])
                block=incoming.blocks[bi];n=_geo_match_cell_corners(block.msh)
                n<3 && continue
                for k in 1:n
                    push!(edges,_geo_match_edge_key(
                        block.nodes[k,c],block.nodes[mod1(k+1,n),c]))
                end
            end
            isempty(edges) && (bounds[key]=out;continue)
            for ((cd,ct),clist) in cells
                cd==1 || continue
                hit=false
                for (bi,c) in clist
                    block=incoming.blocks[bi]
                    _geo_match_edge_key(block.nodes[1,c],block.nodes[2,c]) in edges ||
                        continue
                    hit=true;break
                end
                hit && push!(out,Int32(ct))
            end
        elseif dim==3
            faces=Set{NTuple{4,Int32}}()
            for (bi,c) in get(cells,key,NTuple{2,Int}[])
                block=incoming.blocks[bi];family=msh_family(block.msh)
                rows=view(block.nodes,:,c)
                if family===:tet
                    for f in ((1,2,3),(1,2,4),(1,3,4),(2,3,4))
                        push!(faces,_geo_match_pad4(rows[k] for k in f))
                    end
                elseif family===:hex
                    for f in ((1,2,3,4),(5,6,7,8),(1,2,6,5),(2,3,7,6),
                              (3,4,8,7),(4,1,5,8))
                        push!(faces,_geo_match_pad4(rows[k] for k in f))
                    end
                elseif family===:prism
                    for f in ((1,2,3),(4,5,6),(1,2,5,4),(2,3,6,5),(3,1,4,6))
                        push!(faces,_geo_match_pad4(rows[k] for k in f))
                    end
                elseif family===:pyr
                    for f in ((1,2,3,4),(1,2,5),(2,3,5),(3,4,5),(4,1,5))
                        push!(faces,_geo_match_pad4(rows[k] for k in f))
                    end
                end
            end
            isempty(faces) && (bounds[key]=out;continue)
            for ((fd,ft),flist) in cells
                fd==2 || continue
                hit=false
                for (bi,c) in flist
                    block=incoming.blocks[bi];n=_geo_match_cell_corners(block.msh)
                    n<3 && continue
                    _geo_match_pad4(block.nodes[k,c] for k in 1:n) in faces || continue
                    hit=true;break
                end
                hit && push!(out,Int32(ft))
            end
        end
        bounds[key]=out
    end
    return bounds
end

# Axis-aligned box distance used for upstream's OBB tie-break: per-axis
# lower/upper endpoint deviations, summed in quadrature. Deterministic;
# equal scores fall back to the lower mesh tag.
function _geo_match_box_score(a::NTuple{6,Float64},b::NTuple{6,Float64})
    s=0.0
    for i in 1:3
        s+=(a[i]-b[i])^2+(a[i+3]-b[i+3])^2
    end
    return sqrt(s)
end

function _geo_match_model_box(m::GeoModel,dim::Int,tag::Int)
    return try
        Model._model_entity_bounding_box(m,dim,tag,"GeomMeshMatcher")
    catch err
        err isa InterruptException && rethrow();nothing
    end
end

function _geo_match_model_boundary(m::GeoModel,dim::Int,tag::Int)
    return try
        Model.model_boundary(m,[(dim,tag)])
    catch err
        err isa InterruptException && rethrow();nothing
    end
end

# Model point position — native points store their coordinate; discrete
# points store their first node.
function _geo_match_point_pos(m::GeoModel,tag::Int)
    haskey(m.points,tag) &&
        return (m.points[tag][1],m.points[tag][2],m.points[tag][3])
    rec=get(m.discrete,(0,tag),nothing)
    rec===nothing || size(rec.node_coords,2)>=1 || return nothing
    rec===nothing && return nothing
    return (rec.node_coords[1,1],rec.node_coords[2,1],rec.node_coords[3,1])
end

# Lowest-scoring candidate, or the lowest tag on a score tie — deterministic
# even when boxes are unavailable.
function _geo_match_pick(candidates::Set{Int32},m::GeoModel,dim::Int,gtag::Int,
                         entities)
    length(candidates)==1 && return first(candidates)
    box=_geo_match_model_box(m,dim,gtag)
    best=Inf;pick=Int32(0)
    for mtag in sort!(collect(candidates))
        score=box===nothing ? 0.0 :
            _geo_match_box_score(box,entities[(dim,mtag)].bbox)
        (score<best || (score==best && (pick==0 || mtag<pick))) &&
            (best=score;pick=mtag)
    end
    return pick==0 ? Int32(0) : pick
end

# The `(dim, mesh_tag) => model_tag` correspondence, or a loud bounded
# failure carrying upstream's matched/total accounting. Sorted traversal at
# every level keeps matches deterministic.
function _geo_match_entities(m::GeoModel,incoming::MixedMesh,tol::Float64,
                             context,caller::AbstractString)
    data=incoming.entity_data
    data===nothing && throw(ArgumentError(
        "$caller: -match requires an MSH with entity classification"))
    entities=data.entities
    bydim=[sort!(Int32[tag for ((edim,tag),_) in entities if edim==d]) for d in 0:3]
    model_baskets=[sort!(Int[tag for (edim,tag) in Model.model_entities(m)
        if edim==d]) for d in 0:3]
    pairs=Dict{Tuple{Int,Int32},Int}()
    # ── vertices: position within tol, best max-abs score ──
    v2v=Dict{Int,Int32}()
    for gtag in model_baskets[1]
        pos=_geo_match_point_pos(m,gtag)
        pos===nothing && continue
        best=Inf;choice=Int32(0)
        for mtag in bydim[1]
            e=entities[(0,mtag)]
            score=max(abs(pos[1]-e.bbox[1]),abs(pos[2]-e.bbox[2]),
                      abs(pos[3]-e.bbox[3]))
            (score<tol && (score<best || (score==best && (choice==0 || mtag<choice)))) &&
                (best=score;choice=mtag)
        end
        choice==0 && continue
        v2v[gtag]=choice;pairs[(0,choice)]=gtag
    end
    length(pairs)!=length(bydim[1]) && throw(ArgumentError(
        "$caller: Could only match $(length(pairs)) points out of " *
        "$(length(bydim[1])): check mesh/CAD or increase " *
        "Geometry.MatchMeshTolerance (currently $tol)"))
    inbound=_geo_match_inboundaries(incoming,data)
    # ── curves: shared matched endpoints; closed model curves admit any
    # closed mesh curve; multi-candidates resolve on the bbox score ──
    incident=Dict{Int32,Set{Int32}}()       # mesh point tag -> incident curves
    closed_curves=Set{Int32}()
    for mtag in bydim[2]
        b=get(inbound,(1,mtag),Set{Int32}())
        if isempty(b)
            push!(closed_curves,mtag)
        else
            for p in b
                push!(get!(()->Set{Int32}(),incident,p),mtag)
            end
        end
    end
    e2e=Dict{Int,Int32}();matched_curves=0
    for gtag in model_baskets[2]
        boundary=_geo_match_model_boundary(m,1,gtag)
        boundary===nothing && continue
        ep=Int[abs(p[2]) for p in boundary if p[1]==0]
        candidates=Set{Int32}()
        if isempty(ep) || allequal(ep)
            candidates=closed_curves
        elseif length(ep)>=2 && haskey(v2v,ep[1]) && haskey(v2v,ep[2])
            candidates=intersect(get(incident,v2v[ep[1]],Set{Int32}()),
                                 get(incident,v2v[ep[2]],Set{Int32}()))
        end
        isempty(candidates) && continue
        choice=_geo_match_pick(candidates,m,1,gtag,entities)
        choice==0 && continue
        e2e[gtag]=choice;pairs[(1,choice)]=gtag;matched_curves+=1
    end
    matched_curves!=length(bydim[2]) && throw(ArgumentError(
        "$caller: Could only match $matched_curves curves out of " *
        "$(length(bydim[2])): stopping match"))
    # ── faces: candidates must contain every matched boundary curve;
    # unmatched surfaces are advisory upstream (elements are simply not
    # transferred) — mirror with a warning ──
    f2f=Dict{Int,Int32}();matched_faces=0
    for gtag in model_baskets[3]
        boundary=_geo_match_model_boundary(m,2,gtag)
        boundary===nothing && continue
        required=Set{Int32}();complete=true
        for (bd,bt) in boundary
            bd==1 || continue
            mc=get(e2e,bt,nothing)
            mc===nothing && (complete=false;break)
            push!(required,mc)
        end
        complete || continue
        candidates=Set{Int32}()
        for mtag in bydim[3]
            issubset(required,get(inbound,(2,mtag),Set{Int32}())) &&
                push!(candidates,mtag)
        end
        isempty(candidates) && continue
        choice=_geo_match_pick(candidates,m,2,gtag,entities)
        choice==0 && continue
        f2f[gtag]=choice;pairs[(2,choice)]=gtag;matched_faces+=1
    end
    matched_faces!=length(bydim[3]) && _geo_msg_warn!(context,
        "GeomMeshMatcher: could only match $matched_faces surfaces out of " *
        "$(length(bydim[3])) — unmatched face elements are not transferred")
    # ── regions: every model boundary face must land on the candidate ──
    matched_regions=0
    for gtag in model_baskets[4]
        boundary=_geo_match_model_boundary(m,3,gtag)
        boundary===nothing && continue
        required=Set{Int32}();complete=true
        for (bd,bt) in boundary
            bd==2 || continue
            mf=get(f2f,bt,nothing)
            mf===nothing && (complete=false;break)
            push!(required,mf)
        end
        complete || continue
        candidates=Set{Int32}()
        for mtag in bydim[4]
            issubset(required,get(inbound,(3,mtag),Set{Int32}())) &&
                push!(candidates,mtag)
        end
        isempty(candidates) && continue
        choice=_geo_match_pick(candidates,m,3,gtag,entities)
        choice==0 && continue
        pairs[(3,choice)]=gtag;matched_regions+=1
    end
    matched_regions!=length(model_baskets[4]) && throw(ArgumentError(
        "$caller: Volumes matched: $matched_regions / " *
        "$(length(model_baskets[4])) — failed to match mesh to CAD"))
    return pairs
end

# Rebuild the imported mesh bound to the matched model entities: keys move
# to model tags, unmatched entities drop entirely, cells survive only when
# every node's owning entity matched (upstream copies elements through the
# per-entity vertex map), and entity physicals follow upstream's direction —
# faces/regions take the mesh's, points/curves keep the model's.
function _geo_match_rebind(m::GeoModel,incoming::MixedMesh,
                           pairs::Dict{Tuple{Int,Int32},Int})
    data=incoming.entity_data
    matched=Set{Tuple{Int,Int32}}(keys(pairs))
    keep=[(data.node_entities[r][1],data.node_entities[r][2]) in matched
          for r in 1:nnodes(incoming)]
    rows=Int32.(findall(keep))
    remap=zeros(Int32,nnodes(incoming));remap[rows]=Int32.(1:length(rows))
    # New entity records on model tags.
    new_entities=Dict{Tuple{Int,Int},Elements.MixedEntity}()
    for ((dim,mtag),entity) in data.entities
        gtag=get(pairs,(dim,mtag),nothing)
        gtag===nothing && continue
        phys=dim>=2 ? copy(entity.physical_tags) :
            Int32.(abs.(get(m.entity_physicals,(dim,gtag),Int[])))
        bounds=Int32[]
        for b in entity.boundaries
            partner=get(pairs,(dim-1,Int32(abs(b))),nothing)
            partner===nothing || push!(bounds,Int32(b<0 ? -partner : partner))
        end
        new_entities[(dim,gtag)]=Elements._native_mixed_entity(dim,gtag,entity.bbox;
            physical_tags=phys,boundaries=bounds,
            embedded_curves=copy(entity.embedded_curves),
            partitions=copy(entity.partitions),parent=entity.parent)
    end
    # Cells survive when their owner matched and every node row survives.
    blocks=ElementBlock[];block_entities=Vector{Int32}[]
    external_element_tags=Vector{UInt64}[]
    for (bi,block) in enumerate(incoming.blocks)
        block isa SpecialElementBlock && continue
        dim=msh_dimension(block.msh)
        owners=data.block_entities[bi]
        keepcols=[c for c in axes(block.nodes,2) if
            (dim,owners[c]) in matched &&
            all(r->remap[r]>0,view(block.nodes,:,c))]
        isempty(keepcols) && continue
        tags=Int32[]
        sizehint!(tags,length(keepcols))
        for c in keepcols
            e=get(new_entities,(dim,Int(pairs[(dim,owners[c])])),nothing)
            push!(tags,(e===nothing || isempty(e.physical_tags)) ?
                Int32(0) : first(e.physical_tags))
        end
        push!(blocks,ElementBlock(block.msh,
            remap[block.nodes[:,keepcols]],tags))
        push!(block_entities,Int32[pairs[(dim,owners[c])] for c in keepcols])
        push!(external_element_tags,data.external_element_tags[bi][keepcols])
    end
    links=Elements.MixedPeriodicLink[]
    for link in incoming.periodic_links
        sp=get(pairs,(link.dim,Int32(link.slave_entity)),nothing)
        mp=get(pairs,(link.dim,Int32(link.master_entity)),nothing)
        (sp===nothing || mp===nothing) && continue
        all(r->remap[r]>0,link.slave_nodes) &&
            all(r->remap[r]>0,link.master_nodes) || continue
        push!(links,Elements.MixedPeriodicLink(link.dim,sp,mp,
            remap[link.slave_nodes],remap[link.master_nodes];affine=link.affine))
    end
    node_entities=Tuple{Int,Int32}[
        (data.node_entities[r][1],
         Int32(pairs[(data.node_entities[r][1],data.node_entities[r][2])]))
        for r in rows]
    new_data=Elements._native_mixed_entity_data(new_entities;
        node_entities=node_entities,
        node_parametric=[data.node_parametric[r] for r in rows],
        external_node_tags=data.external_node_tags[rows],
        block_entities=block_entities,
        external_element_tags=external_element_tags)
    return MixedMesh(incoming.coords[:,rows],blocks;
        physical_names=copy(incoming.physical_names),entity_data=new_data,
        elementary_entities=new_data.block_entities,periodic_links=links,
        ancillary_sections=incoming.ancillary_sections)
end

# Model-side physical records after a successful match: the file's
# `$PhysicalNames` replace the model's name table (upstream
# `setPhysicalNames` assigns the whole map), and dim>=2 entities adopt the
# mesh's physical memberships while dim<=1 keep their own.
function _geo_match_model_physicals!(m::GeoModel,incoming::MixedMesh,
                                     pairs::Dict{Tuple{Int,Int32},Int})
    empty!(m.physical_names)
    for (key,name) in incoming.physical_names
        m.physical_names[(Int(key[1]),Int(key[2]))]=name
    end
    data=incoming.entity_data
    for ((dim,mtag),gtag) in pairs
        dim>=2 || continue
        entity=data.entities[(dim,mtag)]
        phys=Int.(entity.physical_tags)
        for (pd,pt) in collect(keys(m.physical))
            pd==dim && filter!(t->t!=gtag,m.physical[(pd,pt)])
        end
        for pt in phys
            members=get!(()->Int[],m.physical,(dim,pt))
            gtag in members || push!(members,gtag)
            m.physical_tag_max=max(m.physical_tag_max,pt)
        end
        m.entity_physicals[(dim,gtag)]=phys
    end
    return nothing
end

# Remove previous actual cells/nodes owned by matched model entities — the
# upstream `deleteMesh` on each matched entity — and drop cells that
# reference removed nodes. Returns `(mesh, remap, dropped_uids)` where
# `remap[old_row]` is the stripped row (0 when dropped) and `dropped_uids`
# the removed elements' external tags.
function _geo_match_strip_actual(prev::MixedMesh,
                                 matched::Set{Tuple{Int,Int32}},
                                 caller::AbstractString)
    data=prev.entity_data
    data===nothing && return prev,collect(Int32,1:nnodes(prev)),Set{UInt64}()
    drop=[(data.node_entities[r][1],data.node_entities[r][2]) in matched
          for r in 1:nnodes(prev)]
    keep_cols=Vector{Int}[]
    for (bi,block) in enumerate(prev.blocks)
        block isa SpecialElementBlock &&
            (push!(keep_cols,collect(eachindex(block.tags)));continue)
        dim=msh_dimension(block.msh);owners=data.block_entities[bi]
        push!(keep_cols,[c for c in axes(block.nodes,2) if
            !((dim,owners[c]) in matched) &&
            all(r->!drop[r],view(block.nodes,:,c))])
    end
    referenced=falses(nnodes(prev))
    for (bi,block) in enumerate(prev.blocks),c in keep_cols[bi],
            r in view(block.nodes,:,c)
        referenced[r]=true
    end
    rows=Int32[r for r in 1:nnodes(prev) if !drop[r] && referenced[r]]
    remap=zeros(Int32,nnodes(prev));remap[rows]=Int32.(1:length(rows))
    dropped_uids=Set{UInt64}()
    for (bi,block) in enumerate(prev.blocks),c in axes(block.nodes,2)
        c in keep_cols[bi] ||
            push!(dropped_uids,data.external_element_tags[bi][c])
    end
    blocks=ElementBlock[]
    for (bi,block) in enumerate(prev.blocks)
        isempty(keep_cols[bi]) && continue
        push!(blocks,ElementBlock(block.msh,
            remap[block.nodes[:,keep_cols[bi]]],block.tags[keep_cols[bi]]))
    end
    links=[Elements.MixedPeriodicLink(link.dim,link.slave_entity,
        link.master_entity,remap[link.slave_nodes],remap[link.master_nodes];
        affine=link.affine) for link in prev.periodic_links if
        all(r->remap[r]>0,link.slave_nodes) &&
        all(r->remap[r]>0,link.master_nodes)]
    new_data=Elements._native_mixed_entity_data(copy(data.entities);
        node_entities=data.node_entities[rows],
        node_parametric=data.node_parametric[rows],
        external_node_tags=data.external_node_tags[rows],
        block_entities=[data.block_entities[bi][keep_cols[bi]] for
            bi in eachindex(prev.blocks) if !isempty(keep_cols[bi])],
        external_element_tags=[data.external_element_tags[bi][keep_cols[bi]]
            for bi in eachindex(prev.blocks) if !isempty(keep_cols[bi])])
    return MixedMesh(prev.coords[:,rows],blocks;
        physical_names=copy(prev.physical_names),entity_data=new_data,
        elementary_entities=new_data.block_entities,periodic_links=links,
        ancillary_sections=prev.ancillary_sections,
        partition_data=prev.partition_data),remap,dropped_uids
end

# Entity union for the rebound join: on matched keys the rebound record
# already encodes the final physicals, so it wins; other keys merge like the
# ordinary import union.
function _geo_match_entity_union(a,b)
    output=copy(a)
    for (key,incoming) in b
        previous=get(output,key,nothing)
        previous===nothing && (output[key]=incoming;continue)
        box=previous.dim==0 ? previous.bbox : ntuple(i->i<=3 ?
            min(previous.bbox[i],incoming.bbox[i]) :
            max(previous.bbox[i],incoming.bbox[i]),6)
        output[key]=Elements._native_mixed_entity(previous.dim,previous.tag,box;
            physical_tags=copy(incoming.physical_tags),
            boundaries=unique(vcat(previous.boundaries,incoming.boundaries)),
            embedded_curves=unique(vcat(previous.embedded_curves,
                incoming.embedded_curves)),
            partitions=previous.partitions,parent=previous.parent)
    end
    return output
end

function _geo_match_join_actual(a,b,caller)
    a===nothing && return b
    ad=a.entity_data;bd=b.entity_data
    total=Base.checked_add(nnodes(a),nnodes(b))
    total<=typemax(Int32) || throw(ArgumentError(
        "$caller: merged node rows exceed Int32 indexing"))
    node_offset=Int32(nnodes(a))
    entities=_geo_match_entity_union(ad.entities,bd.entities)
    blocks=ElementBlock[]
    for (mesh,data,offset) in ((a,ad,Int32(0)),(b,bd,node_offset)),
            (block,owners) in zip(mesh.blocks,data.block_entities)
        projected=Vector{Int32}(undef,length(owners))
        for i in eachindex(owners)
            entity=get(entities,(msh_dimension(block.msh),Int(owners[i])),nothing)
            projected[i]=(entity===nothing || isempty(entity.physical_tags)) ?
                Int32(0) : first(entity.physical_tags)
        end
        push!(blocks,ElementBlock(block.msh,block.nodes.+offset,projected))
    end
    links=vcat(a.periodic_links,[Elements.MixedPeriodicLink(link.dim,
        link.slave_entity,link.master_entity,link.slave_nodes.+node_offset,
        link.master_nodes.+node_offset;affine=link.affine)
        for link in b.periodic_links])
    data=Elements._native_mixed_entity_data(entities;
        node_entities=vcat(ad.node_entities,bd.node_entities),
        node_parametric=vcat(ad.node_parametric,bd.node_parametric),
        external_node_tags=vcat(ad.external_node_tags,bd.external_node_tags),
        block_entities=vcat(ad.block_entities,bd.block_entities),
        external_element_tags=vcat(ad.external_element_tags,
            bd.external_element_tags))
    names=merge(copy(b.physical_names),a.physical_names)
    return MixedMesh(hcat(a.coords,b.coords),blocks;physical_names=names,
        entity_data=data,elementary_entities=data.block_entities,
        periodic_links=links,
        ancillary_sections=vcat(a.ancillary_sections,b.ancillary_sections),
        partition_data=a.partition_data)
end

# Strip previous carrier cells bound to dropped actual elements. A simplex
# `Mesh` carrier uses the authority's `carrier_cells` UID binding; a
# `MixedMesh` carrier filters owners directly. Returns
# `(carrier, kept_old_rows)` — `kept_old_rows` indexes the old carrier.
function _geo_match_strip_carrier(authority,carrier,dropped_uids,matched)
    carrier===nothing && return nothing,Int32[]
    if carrier isa Mesh
        fams=_model_mesh_carrier_blocks_for_match(carrier)
        bindings=authority===nothing ? Union{Nothing,UInt64}[] :
            authority.carrier_cells
        keep=Vector{Int}[]
        for (fi,cells) in enumerate(fams)
            bound=fi<=length(bindings) ? bindings[fi] : Union{Nothing,UInt64}[]
            push!(keep,[c for c in axes(cells,2) if
                !(c<=length(bound) && bound[c]!==nothing &&
                  bound[c] in dropped_uids)])
        end
        used=Set{Int32}()
        for (fi,cells) in enumerate(fams),c in keep[fi],r in view(cells,:,c)
            push!(used,r)
        end
        rows=sort!(collect(used))
        remap=zeros(Int32,nnodes(carrier));remap[rows]=Int32.(1:length(rows))
        return Mesh(carrier.coords[:,rows];
            segs=remap[fams[1][:,keep[1]]],tris=remap[fams[2][:,keep[2]]],
            tets=remap[fams[3][:,keep[3]]],
            seg_tag=carrier.seg_tag[keep[1]],tri_tag=carrier.tri_tag[keep[2]],
            tet_tag=carrier.tet_tag[keep[3]]),rows
    elseif carrier isa MixedMesh
        data=carrier.entity_data
        data===nothing && return carrier,collect(Int32,1:nnodes(carrier))
        keep_cols=Vector{Int}[]
        for (bi,block) in enumerate(carrier.blocks)
            block isa SpecialElementBlock &&
                (push!(keep_cols,collect(eachindex(block.tags)));continue)
            dim=msh_dimension(block.msh);owners=data.block_entities[bi]
            push!(keep_cols,[c for c in axes(block.nodes,2) if
                !((dim,owners[c]) in matched)])
        end
        used=Set{Int32}()
        for (bi,block) in enumerate(carrier.blocks),c in keep_cols[bi],
                r in view(block.nodes,:,c)
            push!(used,r)
        end
        rows=sort!(collect(used))
        remap=zeros(Int32,nnodes(carrier));remap[rows]=Int32.(1:length(rows))
        blocks=ElementBlock[ElementBlock(block.msh,
            remap[block.nodes[:,keep_cols[bi]]],block.tags[keep_cols[bi]])
            for (bi,block) in enumerate(carrier.blocks) if !isempty(keep_cols[bi])]
        new_data=Elements._native_mixed_entity_data(copy(data.entities);
            node_entities=data.node_entities[rows],
            node_parametric=data.node_parametric[rows],
            external_node_tags=data.external_node_tags[rows],
            block_entities=[data.block_entities[bi][keep_cols[bi]] for
                bi in eachindex(carrier.blocks) if !isempty(keep_cols[bi])],
            external_element_tags=[data.external_element_tags[bi][keep_cols[bi]]
                for bi in eachindex(carrier.blocks) if !isempty(keep_cols[bi])])
        return MixedMesh(carrier.coords[:,rows],blocks;
            physical_names=copy(carrier.physical_names),entity_data=new_data,
            elementary_entities=new_data.block_entities,
            periodic_links=carrier.periodic_links,
            ancillary_sections=carrier.ancillary_sections,
            partition_data=carrier.partition_data),rows
    end
    return carrier,collect(Int32,1:nnodes(carrier))
end

_model_mesh_carrier_blocks_for_match(mesh::Mesh)=(mesh.segs,mesh.tris,mesh.tets)

# The match tail mirroring `_geo_import_mesh_file!`: join the stripped
# previous actual with the rebound import, rebuild authority/products/
# context, then publish — everything runs only after the match succeeded.
function _geo_match_install!(m::GeoModel,context,incoming::MixedMesh,
                             rebound::MixedMesh,pairs,identities,
                             ghost_dimension,caller::AbstractString)
    previous=m.meshing.mesh_authority
    matched_model=Set{Tuple{Int,Int32}}((d,Int32(t)) for ((d,_),t) in pairs)
    if context.mesh!==nothing && (previous===nothing || previous.mesh===nothing ||
            previous.carrier!==context.mesh)
        working=_geo_mesh_working_model(m,context)
        plan=Model._model_mesh_generate_higher(working,context.mesh,
            context.mesh_parts,context.mesh_part_rows,context.mesh_node_owner,
            0,nothing,caller;operation=:bind_completed,
            initial_max_node_tag=context.mesh_node_max,
            initial_max_element_tag=context.mesh_element_max)
        _geo_publish_actual_plan!(m,context,plan,caller;carrier=context.mesh)
        previous=m.meshing.mesh_authority
    end
    prev_actual=previous===nothing ? nothing : previous.mesh
    stripped_prev=prev_actual;actual_remap=nothing;dropped_uids=Set{UInt64}()
    if prev_actual!==nothing && !isempty(matched_model)
        stripped_prev,actual_remap,dropped_uids=
            _geo_match_strip_actual(prev_actual,matched_model,caller)
    end
    actual=_geo_match_join_actual(stripped_prev,rebound,caller)
    actual,node_mapping=_geo_import_order_nodes(actual)
    imported_carrier=_geo_mixed_simplex_mesh(rebound,caller)
    old_carrier_rows=previous===nothing ? Int32[] : previous.carrier_rows
    stripped_carrier=context.mesh;kept_carrier_rows=Int32[]
    if context.mesh!==nothing && previous!==nothing
        stripped_carrier,kept_carrier_rows=_geo_match_strip_carrier(
            previous,context.mesh,dropped_uids,matched_model)
    elseif context.mesh!==nothing
        kept_carrier_rows=collect(Int32,1:nnodes(context.mesh))
    end
    carrier=stripped_carrier===nothing ? imported_carrier :
        _geo_concat_meshes(stripped_carrier,imported_carrier)
    incoming_cells=_geo_import_simplex_cell_bindings(rebound)
    cells=if previous===nothing
        incoming_cells
    else
        retained=_geo_import_retained_cell_bindings(previous)
        kept=[Union{Nothing,UInt64}[uid for uid in fam if
            uid===nothing || !(uid in dropped_uids)] for fam in retained]
        if previous.carrier isa Mesh
            [vcat(old,new) for (old,new) in zip(kept,incoming_cells)]
        else
            vcat(kept,[tags for (tags,(_,nodes)) in zip(incoming_cells,
                Model._model_mesh_carrier_blocks(imported_carrier))
                if size(nodes,2)>0])
        end
    end
    done=previous===nothing ? Set{Tuple{Int,Int32}}() : copy(previous.done)
    union!(done,rebound.entity_data.node_entities)
    for (block,owners) in zip(rebound.blocks,rebound.entity_data.block_entities)
        union!(done,((msh_dimension(block.msh),owner) for owner in owners))
    end
    ranks=previous===nothing ? UInt64[] : copy(previous.private_rank)
    if actual_remap===nothing
        next_rank=maximum(ranks;init=UInt64(0))
        UInt64(nnodes(rebound))<=typemax(UInt64)-next_rank || throw(ArgumentError(
            "$caller: imported private vertex rank overflows UInt64"))
        append!(ranks,[next_rank+UInt64(row) for row in 1:nnodes(rebound)])
    else
        kept=[r for r in eachindex(ranks) if actual_remap[r]>0]
        base=maximum(ranks;init=UInt64(0))
        ranks=vcat(ranks[kept],[base+UInt64(i) for i in 1:nnodes(rebound)])
    end
    ordered_ranks=similar(ranks);ordered_ranks[node_mapping]=ranks
    # Carrier rows: kept old carrier rows map through the actual-strip and
    # join into combined space; rebound rows follow the stripped actual.
    stripped_count=stripped_prev===nothing ? 0 : nnodes(stripped_prev)
    if previous===nothing
        rows=node_mapping[Int32.(1:nnodes(rebound))]
    else
        old_rows=Int32[]
        sizehint!(old_rows,length(kept_carrier_rows))
        for k in kept_carrier_rows
            arow=old_carrier_rows[k]
            push!(old_rows,actual_remap===nothing ? arow : actual_remap[arow])
        end
        combined=vcat(old_rows,
            Int32.(stripped_count .+ (1:nnodes(rebound))))
        rows=node_mapping[combined]
    end
    authority=Model._model_mesh_authority(actual,done;caller,carrier,
        carrier_rows=rows,carrier_cells=cells,private_rank=ordered_ranks)
    if previous!==nothing
        native_inputs=copy(previous.native_inputs)
        producer_rows=Dict(owner=>begin
                mapped=actual_remap===nothing ? positions :
                    actual_remap[positions]
                node_mapping[Int32[p for p in mapped if p>0]]
            end for (owner,positions) in previous.producer_rows)
        for owner in matched_model
            delete!(native_inputs,owner);delete!(producer_rows,owner)
        end
        authority=Model._ModelMeshAuthority(authority.mesh,authority.products,
            authority.done,authority.private_rank,native_inputs,carrier,rows,
            producer_rows,authority.carrier_cells)
        constructors=Dict(owner=>constructor for (owner,constructor) in
            previous.surface_constructor_inputs if !(owner in matched_model))
        authority=Model._model_mesh_with_surface_constructor_inputs(authority,
            constructors;caller)
    end
    parts,part_rows=_geo_actual_product_parts(authority,caller)
    carrier_owners=_geo_actual_carrier_owners(authority,carrier)
    geometry_receipt=_geo_geometry_rebind_receipt(authority,carrier,
        context.mesh_geometry_parts,nothing,caller)
    if previous!==nothing && previous.emitted_geometry!==nothing &&
            geometry_receipt!==nothing
        authority=Model._model_mesh_with_emitted_receipt(authority,
            geometry_receipt.parts,geometry_receipt.rows;caller)
    end
    _geo_merge_discrete!(m,rebound,caller)
    _geo_import_entity_receipts!(m,rebound;ghost_dimension)
    m.meshing.mesh_authority=authority
    context.mesh=carrier
    context.mesh_parts=parts;context.mesh_part_rows=part_rows
    context.mesh_geometry_parts=geometry_receipt
    context.mesh_node_owner=carrier_owners
    context.mesh_node_parametric=
        Model._model_mesh_runtime_parameters(m,actual)[rows]
    context.mesh_node_max=max(context.mesh_node_max,
        maximum(identities.node_tags;init=UInt64(0)))
    context.mesh_element_max=max(context.mesh_element_max,
        maximum((maximum(tags;init=UInt64(0)) for tags in identities.element_tags);
            init=UInt64(0)))
    _geo_match_element_lookup!(context,rebound)
    return nothing
end

# Rebound-aware element-lookup cache — `_geo_import_element_lookup!` zips
# raw identity blocks, which no longer align with the filtered rebound
# blocks; rebuild the same cache shape from the rebound authority metadata.
function _geo_match_element_lookup!(context,rebound)
    receipt=context.native_element_lookup_receipt
    if receipt===nothing
        receipt=(vector=Dict{UInt64,Tuple{Int32,Int,Int}}(),
            vector_max=nothing,map=Dict{UInt64,Tuple{Int32,Int,Int}}(),
            node_vector=Dict{UInt64,Bool}(),node_vector_max=nothing,
            node_map=Dict{UInt64,Bool}())
    end
    data=rebound.entity_data
    vector=copy(receipt.vector);map=copy(receipt.map)
    for (block,tags,owners) in zip(rebound.blocks,
            data.external_element_tags,data.block_entities),
            column in eachindex(tags)
        uid=tags[column]
        get!(map,uid) do
            (Int32(block.msh),msh_dimension(block.msh),Int(owners[column]))
        end
    end
    node_map=copy(receipt.node_map);node_vector=copy(receipt.node_vector)
    for uid in data.external_node_tags
        node_map[uid]=true
    end
    context.native_element_lookup_receipt=(vector=vector,
        vector_max=receipt.vector_max,map=map,node_vector=node_vector,
        node_vector_max=receipt.node_vector_max,node_map=node_map)
    return nothing
end

# `GeomMeshMatcher::match` — read the MSH, scale it, pair every importable
# entity with a model entity, bind physicals, then install the rebound mesh
# onto the model. Any failed dimension count aborts before the model or
# context is touched.
function _geo_match_mesh_file!(m::GeoModel,context,file,caller)
    tol=Float64(something(
        _geo_option_number(context,"Geometry",0,"MatchMeshTolerance"),1e-6))
    scale=Float64(something(
        _geo_option_number(context,"Geometry",0,"MatchMeshScaleFactor"),1.0))
    switch_tags=_geo_mesh_switch_tags(context)
    mixed,identities=Elements._read_mixed_msh_with_identities(file;switch_tags)
    file_dimension=mixed.entity_data===nothing ? -1 :
        maximum(first,keys(mixed.entity_data.entities);init=-1)
    ghost_dimension=max(file_dimension,
        maximum(first,Model.model_entities(m);init=-1))
    incoming=_geo_import_actual_mesh(mixed,identities,caller;ghost_dimension)
    incoming.partition_data!==nothing && throw(ArgumentError(
        "$caller: -match on a partitioned MSH is not supported"))
    if !(scale==1.0)
        entities=Dict{Tuple{Int,Int},Elements.MixedEntity}()
        for (key,e) in incoming.entity_data.entities
            lo=Tuple(sort!(Float64[e.bbox[i]*scale for i in 1:3]))[[1,2,3]]
            hi=Tuple(sort!(Float64[e.bbox[i+3]*scale for i in 1:3]))[[1,2,3]]
            entities[key]=Elements._native_mixed_entity(e.dim,e.tag,
                (lo[1],lo[2],lo[3],hi[1],hi[2],hi[3]);
                physical_tags=e.physical_tags,boundaries=e.boundaries,
                embedded_curves=e.embedded_curves,partitions=e.partitions,
                parent=e.parent)
        end
        data=incoming.entity_data
        incoming=MixedMesh(incoming.coords.*scale,incoming.blocks;
            physical_names=incoming.physical_names,
            entity_data=Elements._native_mixed_entity_data(entities;
                node_entities=data.node_entities,
                node_parametric=data.node_parametric,
                external_node_tags=data.external_node_tags,
                block_entities=data.block_entities,
                external_element_tags=data.external_element_tags),
            elementary_entities=incoming.elementary_entities,
            periodic_links=incoming.periodic_links,
            ancillary_sections=incoming.ancillary_sections,
            data_sections=incoming.data_sections)
    end
    pairs=_geo_match_entities(m,incoming,tol,context,caller)
    rebound=_geo_match_rebind(m,incoming,pairs)
    _geo_match_model_physicals!(m,incoming,pairs)
    _geo_match_install!(m,context,incoming,rebound,pairs,identities,
        ghost_dimension,caller)
    _geo_msg_warn!(context,"Successfully matched mesh to CAD")
    return nothing
end
