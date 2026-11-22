# Gmsh 4.15.2's surface classification (geo/GModelParametrize.cpp:
# classifyFaces(gm, angle, includeBoundary, forParametrization, curveAngle)
# and classifyFaces(gm, curveAngle), mesh/meshGFaceOptimize.cpp's edge
# angles, geo/MEdge.cpp's SortEdgeConsecutive and
# GModel::pruneMeshVertexAssociations), ported step by step so the patches,
# curves and points come out with Gmsh's tags, element numbers, node
# ownership and element order.
#
# Element numbers follow GModel's element counter: every MLine Gmsh builds
# draws the next number, including the temporary ones the classification
# builds and drops (the sharp, boundary and non-manifold edge lines, and the
# stack MLine it builds for every triangle edge it looks up), so the new
# points and curves start well above the mesh's last element (t13: 8997).

const _CF_Edge = NTuple{2,UInt64}

_cf_key(a::UInt64, b::UInt64) = a < b ? (a, b) : (b, a)

# Triangle and quadrangle edges in Gmsh's local order.
_cf_local_edges(n::Int) = n == 3 ? ((1, 2), (2, 3), (3, 1)) : ((1, 2), (2, 3), (3, 4), (4, 1))

# edge_angle: the angle between the two elements' normals, from their first
# three vertices (meshGFaceOptimize.cpp).
function _cf_edge_angle(coords, e1::Vector{UInt64}, e2::Vector{UInt64})
    function normal(e)
        p1 = coords[e[1]]; p2 = coords[e[2]]; p3 = coords[e[3]]
        a = (p1[1] - p2[1], p1[2] - p2[2], p1[3] - p2[3])
        b = (p1[1] - p3[1], p1[2] - p3[2], p1[3] - p3[3])
        c = [a[2] * b[3] - a[3] * b[2], -a[1] * b[3] + a[3] * b[1], a[1] * b[2] - a[2] * b[1]]
        mod = sqrt(c[1] * c[1] + c[2] * c[2] + c[3] * c[3])
        if mod != 0.0
            inv = 1.0 / mod
            c[1] *= inv; c[2] *= inv; c[3] *= inv
        end
        return c
    end
    c1 = normal(e1); c2 = normal(e2)
    c3 = (c1[2] * c2[3] - c1[3] * c2[2], -c1[1] * c2[3] + c1[3] * c2[1], c1[1] * c2[2] - c1[2] * c2[1])
    cosa = c1[1] * c2[1] + c1[2] * c2[2] + c1[3] * c2[3]
    sina = sqrt(c3[1] * c3[1] + c3[2] * c3[2] + c3[3] * c3[3])
    return _gm_atan2(sina, cosa)
end

# breakForLargeAngle (GModelParametrize.cpp).
function _cf_break_for_large_angle(coords, vprev::UInt64, vmid::UInt64, vpos::UInt64, threshold::Float64)
    threshold >= pi - 1e-12 && return false
    threshold <= 0 && return true
    p = coords[vprev]; q = coords[vmid]; r = coords[vpos]
    v1 = (q[1] - p[1], q[2] - p[2], q[3] - p[3]); v2 = (r[1] - q[1], r[2] - q[2], r[3] - q[3])
    c = (v1[2] * v2[3] - v1[3] * v2[2], v1[3] * v2[1] - v1[1] * v2[3], v1[1] * v2[2] - v1[2] * v2[1])
    a = _gm_atan2(sqrt(c[1] * c[1] + c[2] * c[2] + c[3] * c[3]), v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3])
    return (a > threshold && a < pi - threshold) || (a > pi + threshold && a < 2pi - threshold)
end

# SortEdgeConsecutive (MEdge.cpp): chains of vertices from a list of edges,
# walking from the smallest-numbered free end (or the smallest vertex of a
# loop). Where Gmsh returns false the chains found so far are returned, as
# classifyFaces uses them regardless.
function _cf_sort_edge_consecutive(edges::Vector{_CF_Edge})
    vs = Vector{Vector{UInt64}}()
    isempty(edges) && return vs
    c = Dict{UInt64,Tuple{UInt64,UInt64}}()      # 0 stands for nullptr
    for (v0, v1) in edges
        if haskey(c, v0)
            c[v0][2] == 0 || return vs
            c[v0] = (c[v0][1], v1)
        else
            c[v0] = (v1, UInt64(0))
        end
        if haskey(c, v1)
            c[v1][2] == 0 || return vs
            c[v1] = (c[v1][1], v0)
        else
            c[v1] = (v0, UInt64(0))
        end
    end
    while !isempty(c)
        v = UInt64[]
        keys_sorted = sort!(collect(keys(c)))
        start = keys_sorted[1]
        for k in keys_sorted
            if c[k][2] == 0
                start = k; break
            end
        end
        first_, second_ = c[start]
        prev = second_ == start ? first_ : second_
        current = start
        while true
            isempty(c) && return vs
            push!(v, current)
            haskey(c, current) || return vs
            v1, v2 = c[current]
            delete!(c, current)
            temp = current
            if v1 == prev
                current = v2
            elseif v2 == prev
                current = v1
            else
                break
            end
            prev = temp
            current == start && push!(v, current)
            (current != start && current != 0) || break
        end
        length(v) > 2 && v[end - 1] == v[end] && pop!(v)
        push!(vs, v)
    end
    return vs
end

# GEdgeLoop(edges) (GEdgeLoop.cpp), the curve loops GFace::setBoundEdges
# makes of a face's curves: signed tags in loop order. `ends[c]` is curve
# c's (begin point, end point).
function _cf_edge_loop(curves::Vector{Int}, ends::Dict{Int,Tuple{Int,Int}})
    wire = copy(curves)
    loop = Int[]
    begin_of(s) = s > 0 ? ends[s][1] : ends[-s][2]
    end_of(s) = s > 0 ? ends[s][2] : ends[-s][1]
    function next_one(previous)
        if previous == 0
            ge = popfirst!(wire)
            return ge
        end
        gv = end_of(previous)
        choices = [ge for ge in wire if ends[ge][1] == gv || ends[ge][2] == gv]
        for ge in choices
            count(==(ge), choices) == 2 || continue
            filter!(!=(ge), wire); push!(wire, ge)
            ends[ge][1] == gv && return ge
            ends[ge][2] == gv && return -ge
        end
        for ge in choices
            ge == abs(previous) && continue
            filter!(!=(ge), wire)
            ends[ge][1] == gv && return ge
            ends[ge][2] == gv && return -ge
        end
        return 0
    end
    while !isempty(wire)
        previous = 0
        while !isempty(wire)
            s = next_one(previous)
            s == 0 && break
            push!(loop, s); previous = s
        end
    end
    return loop
end

function _cf_max_entity_tag(m::GeoModel, dim::Int)
    result = 0
    for tag in keys(_model_entity_dictionary(m, dim)); result = max(result, tag); end
    for (d, tag) in keys(m.discrete); d == dim && (result = max(result, tag)); end
    return result
end

# Replace a record's nodes by `tags` (coordinates from `coords`), keeping
# the mesh-vertex parameters a node carried.
function _cf_set_record_nodes!(record::DiscreteEntity, tags::Vector{UInt64}, coords, native)
    record.node_tags = tags
    record.node_coords = Matrix{Float64}(undef, 3, length(tags))
    for (j, t) in enumerate(tags)
        p = coords[t]
        record.node_coords[1, j] = p[1]; record.node_coords[2, j] = p[2]; record.node_coords[3, j] = p[3]
    end
    record.node_params = Matrix{Float64}(undef, 0, 0)
    record.node_native_params = Dict{UInt64,Vector{Float64}}(t => get(native, t, Float64[]) for t in tags)
    return record
end

function _classify_faces_gmsh!(m::GeoModel, angle::Float64, include_boundary::Bool, curve_angle::Float64)
    records = _discrete_mesh_records_model(m)
    coords = Dict{UInt64,NTuple{3,Float64}}()
    native = Dict{UInt64,Vector{Float64}}()
    for (_, record) in records
        for (j, t) in enumerate(record.node_tags)
            coords[t] = (record.node_coords[1, j], record.node_coords[2, j], record.node_coords[3, j])
            p = get(record.node_native_params, t, nothing)
            p === nothing || isempty(p) || (native[t] = p)
        end
    end
    face_keys = sort!([key for (key, _) in records if key[1] == 2]; by=last)
    record_of = Dict(key => record for (key, record) in records)
    element_max = Int(_model_fresh_element_tag(m)) - 1
    next_element() = (element_max += 1; UInt64(element_max))

    # Outer classifyFaces: the sharp and boundary edges become lines of one
    # temporary curve.
    max1_outer = _cf_max_entity_tag(m, 1)
    elements = Vector{Vector{UInt64}}()
    for key in face_keys       # each face's triangles, then its quadrangles
        record = record_of[key]
        for (i, type) in enumerate(record.element_types)
            type == 2 && push!(elements, record.element_nodes[i][1:3])
        end
        for (i, type) in enumerate(record.element_types)
            type == 3 && push!(elements, record.element_nodes[i][1:4])
        end
    end
    adj = Dict{_CF_Edge,Tuple{_CF_Edge,Int,Int}}()   # key => (first MEdge, t1, t2)
    for (i, e) in enumerate(elements), (a, b) in _cf_local_edges(length(e))
        k = _cf_key(e[a], e[b])
        entry = get(adj, k, nothing)
        adj[k] = entry === nothing ? ((e[a], e[b]), i, 0) : (entry[1], entry[2], i)
    end
    detected = Tuple{_CF_Edge,Float64}[]; lonely = _CF_Edge[]
    for k in sort!(collect(keys(adj)))
        edge, t1, t2 = adj[k]
        t2 == 0 ? push!(lonely, edge) : push!(detected, (edge, _cf_edge_angle(coords, elements[t1], elements[t2])))
    end
    sort!(detected; by=last, rev=true)
    cut = _CF_Edge[]
    for (edge, a) in detected
        a <= angle && break
        push!(cut, edge); next_element()
    end
    if include_boundary
        for edge in lonely; push!(cut, edge); next_element(); end
    end
    # computeNonManifoldEdges(gm, cut, true): triangles only, boundary always.
    count = Dict{_CF_Edge,Tuple{_CF_Edge,Int}}()
    for key in face_keys, (i, type) in enumerate(record_of[key].element_types)
        type == 2 || continue
        e = record_of[key].element_nodes[i]
        for (a, b) in _cf_local_edges(3)
            k = _cf_key(e[a], e[b]); entry = get(count, k, nothing)
            count[k] = entry === nothing ? ((e[a], e[b]), 1) : (entry[1], entry[2] + 1)
        end
    end
    for k in sort!(collect(keys(count)))
        edge, n = count[k]
        n > 2 && (push!(cut, edge); next_element())
        n == 1 && (push!(cut, edge); next_element())
    end
    lines = Dict{_CF_Edge,_CF_Edge}()          # makeMLinesUnique, then the set
    for edge in sort(cut; by=e -> _cf_key(e...))
        k = _cf_key(edge...); haskey(lines, k) || (lines[k] = edge)
    end

    # Inner classifyFaces(gm, curveAngle).
    max0 = _cf_max_entity_tag(m, 0)
    max1 = max(_cf_max_entity_tag(m, 1), max1_outer + 1)
    max2 = _cf_max_entity_tag(m, 2)
    for (key, record) in records
        key[1] == 1 || continue
        for (i, type) in enumerate(record.element_types)
            type == 1 || continue
            e = record.element_nodes[i]; k = _cf_key(e[1], e[2])
            haskey(lines, k) || (lines[k] = (e[1], e[2]))
        end
    end
    triangles = Tuple{UInt64,Vector{UInt64},Tuple{Int,Int}}[]   # (tag, nodes, old face)
    tris = Dict{_CF_Edge,Vector{Int}}()
    for key in face_keys, (i, type) in enumerate(record_of[key].element_types)
        type == 2 || continue
        push!(triangles, (record_of[key].element_tags[i], record_of[key].element_nodes[i], key))
        t = length(triangles); e = triangles[t][2]
        for (a, b) in _cf_local_edges(3)
            push!(get!(Vector{Int}, tris, _cf_key(e[a], e[b])), t)
        end
    end
    isempty(triangles) && return false
    touched = sort!(collect(eachindex(triangles)); by=t -> triangles[t][1])
    touched_set = Set(touched)
    order = Dict(t => i for (i, t) in enumerate(touched))
    face_of = zeros(Int, length(triangles))
    new_faces = Tuple{Int,Vector{Int}}[]
    pointer = 1
    while !isempty(touched_set)
        while !(touched[pointer] in touched_set); pointer += 1; end
        seed = touched[pointer]; delete!(touched_set, seed)
        max2 += 1; face = max2; members = Int[]
        stack = [seed]
        while !isempty(stack)
            t = pop!(stack)
            push!(members, t); face_of[t] = face
            e = triangles[t][2]
            for (a, b) in _cf_local_edges(3)
                k = _cf_key(e[a], e[b])
                owners = tris[k]
                next_element()                      # the stack MLine ll
                haskey(lines, k) && continue
                tt = owners[1] == t ? owners[2] : owners[1]
                if tt in touched_set
                    push!(stack, tt); delete!(touched_set, tt)
                end
            end
        end
        push!(new_faces, (face, members))
    end
    replaced_by = Dict{Tuple{Int,Int},Set{Int}}()
    for (face, members) in new_faces, t in members
        push!(get!(Set{Int}, replaced_by, triangles[t][3]), face)
    end

    # Curves: one per set of adjacent faces, in edge order.
    new_edges = Tuple{Int,Vector{Int},Vector{_CF_Edge}}[]
    temp_tag = max1
    for k in sort!(collect(keys(tris)))
        next_element()                              # the stack MLine ml
        haskey(lines, k) || continue
        faces = [face_of[t] for t in tris[k]]
        length(faces) == 2 && faces[1] == faces[2] && continue
        found = 0
        for (j, (_, fs, _)) in enumerate(new_edges)
            length(fs) == length(faces) || continue
            if all(f -> f in faces, fs) && all(f -> f in fs, faces)
                found = j; break
            end
        end
        if found == 0
            temp_tag += 1
            push!(new_edges, (temp_tag, faces, _CF_Edge[]))
            found = length(new_edges)
        end
        push!(new_edges[found][3], lines[k])
    end
    max1 = temp_tag
    force_split = Set{UInt64}()
    embedded = Dict{Int,Int}()                     # temporary curve => face
    for (j, (_, fs, ls)) in enumerate(new_edges)
        for a in eachindex(fs), b in a+1:length(fs)
            fs[a] == fs[b] || continue
            embedded[j] = fs[a]
            vs_ = _cf_sort_edge_consecutive(ls)
            (vs_ === nothing || isempty(vs_)) && continue
            push!(force_split, vs_[1][1]); push!(force_split, vs_[1][end])
        end
    end
    point_of = Dict{UInt64,Int}()
    point_elements = Dict{Int,UInt64}()
    curves = Tuple{Int,Int,Int,Vector{Tuple{UInt64,_CF_Edge}}}[]   # tag, begin, end, lines
    face_topology = Dict{Int,Vector{Int}}()
    embedded_new = Dict{Int,Int}()
    for (j, (_, fs, ls)) in enumerate(new_edges)
        emb = get(embedded, j, 0)
        vs_ = _cf_sort_edge_consecutive(ls)
        vs_ === nothing && (vs_ = Vector{UInt64}[])
        pieces = Vector{Vector{UInt64}}()
        for i in eachindex(vs_)
            chain = vs_[i]
            periodic = chain[end] == chain[1]
            if periodic
                for q in 1:length(chain)-1
                    v0 = chain[q == 1 ? length(chain) - 1 : q - 1]
                    v1 = chain[q]; v2 = chain[q + 1]
                    if v1 in force_split || _cf_break_for_large_angle(coords, v0, v1, v2, curve_angle)
                        chain = [chain[mod1(r, length(chain))] for r in q:length(chain)+q-1]
                        break
                    end
                end
            end
            cuts = [1]
            for q in 2:length(chain)-1
                (chain[q] in force_split ||
                 _cf_break_for_large_angle(coords, chain[q-1], chain[q], chain[q+1], curve_angle)) && push!(cuts, q)
            end
            push!(cuts, length(chain))
            first_vertex = chain[cuts[1]]
            for c in 2:length(cuts)
                vv = UInt64[]
                for q in cuts[c-1]:cuts[c]
                    (q == cuts[c-1] || chain[q] != chain[q-1]) && push!(vv, chain[q])
                end
                periodic && c == length(cuts) && length(cuts) > 2 && push!(vv, first_vertex)
                push!(pieces, vv)
            end
        end
        for vv in pieces
            ends = Int[]
            for v in (vv[1], vv[end])
                if !haskey(point_of, v)
                    max0 += 1; point_of[v] = max0
                    point_elements[max0] = next_element()
                end
                push!(ends, point_of[v])
            end
            max1 += 1
            segment = [(next_element(), (vv[q-1], vv[q])) for q in 2:length(vv)]
            push!(curves, (max1, ends[1], ends[2], segment))
            emb == 0 || (embedded_new[max1] = emb)
            for f in fs; push!(get!(Vector{Int}, face_topology, f), max1); end
        end
    end

    # Replace the model's curves, points and faces.
    for key in [key for (key, _) in records if key[1] <= 1]
        haskey(m.discrete, key) && delete!(m.discrete, key)
        haskey(m.meshing.attached, key) && delete!(m.meshing.attached, key)
    end
    for key in face_keys
        record = record_of[key]
        keep = findall(t -> t != 2, record.element_types)
        record.element_types = record.element_types[keep]
        record.element_tags = record.element_tags[keep]
        record.element_nodes = record.element_nodes[keep]
    end
    for (v, tag) in point_of
        record = DiscreteEntity()
        _record_append_element!(record, 15, point_elements[tag], UInt64[v])
        m.discrete[(0, tag)] = record
    end
    curve_ends = Dict{Int,Tuple{Int,Int}}(tag => (b, e) for (tag, b, e, _) in curves)
    for (tag, b, e, segment) in curves
        record = DiscreteEntity()
        record.boundary = b == e ? [(0, b)] : [(0, b), (0, e)]
        for (element, (v1, v2)) in segment
            _record_append_element!(record, 1, element, UInt64[v1, v2])
        end
        m.discrete[(1, tag)] = record
    end
    for (face, members) in new_faces
        record = DiscreteEntity()
        for t in members
            _record_append_element!(record, 2, triangles[t][1], triangles[t][2])
        end
        bound = Int[]; emb_curves = Int[]
        for c in get(face_topology, face, Int[])
            get(embedded_new, c, 0) == face ? push!(emb_curves, c) : push!(bound, c)
        end
        record.boundary = [(1, s) for s in _cf_edge_loop(bound, curve_ends)]
        m.discrete[(2, face)] = record
        for c in emb_curves
            embed!(m, 1, [c], 2, face)
        end
    end
    for (key, record) in records       # a region's faces become their pieces
        key[1] == 3 || continue
        faces = Set{Int}()
        for (d, t) in record.boundary
            d == 2 && haskey(replaced_by, (2, abs(t))) && union!(faces, replaced_by[(2, abs(t))])
        end
        record.boundary = [(2, f) for f in sort!(collect(faces))]
    end
    for key in face_keys
        record = record_of[key]
        isempty(record.element_types) || continue
        haskey(m.discrete, key) && delete!(m.discrete, key)
        haskey(m.meshing.attached, key) && delete!(m.meshing.attached, key)
    end

    # pruneMeshVertexAssociations: each node belongs to the lowest-dimension
    # entity whose elements use it; entities list their nodes by number.
    owner = Dict{UInt64,Tuple{Int,Int}}()
    all_records = _discrete_mesh_records_model(m)
    for dim in 3:-1:0
        for (key, record) in sort!([(k, r) for (k, r) in all_records if k[1] == dim]; by=x -> x[1][2])
            for nodes in record.element_nodes, v in nodes
                owner[v] = key
            end
        end
    end
    lists = Dict{Tuple{Int,Int},Vector{UInt64}}()
    for v in sort!(collect(keys(owner)))
        push!(get!(Vector{UInt64}, lists, owner[v]), v)
    end
    for (key, record) in all_records
        nodes = get(lists, key, UInt64[])
        # Nodes no element uses stay where they were (Gmsh drops them).
        for v in record.node_tags
            haskey(owner, v) || push!(nodes, v)
        end
        _cf_set_record_nodes!(record, nodes, coords, native)
    end
    return true
end

# The points' and curves' new MPoint/MLine elements join the retained mesh
# (the mesh snapshot the API reads), one block per entity, in entity order;
# elements of the removed points and curves leave it.
function _cf_publish_new_elements!(m::GeoModel)
    authority = m.meshing.mesh_authority
    (authority === nothing || authority.mesh === nothing) && return nothing
    mesh = authority.mesh
    mesh isa MixedMesh || return nothing
    data = mesh.entity_data
    data === nothing && return nothing
    row = Dict{UInt64,Int32}(t => Int32(i) for (i, t) in enumerate(data.external_node_tags))
    live = Set{UInt64}()
    for (_, record) in _discrete_mesh_records_model(m); union!(live, record.element_tags); end
    # Drop point and line blocks whose elements no record holds any more.
    keep = Int[]
    for (bi, block) in enumerate(mesh.blocks)
        stale = block isa ElementBlock && msh_dimension(block.msh) <= 1 &&
            bi <= length(data.external_element_tags) &&
            !any(t -> t in live, data.external_element_tags[bi])
        stale || push!(keep, bi)
    end
    if length(keep) != length(mesh.blocks)
        blocks = mesh.blocks[keep]; empty!(mesh.blocks); append!(mesh.blocks, blocks)
        entities = data.block_entities[keep]; empty!(data.block_entities); append!(data.block_entities, entities)
        tags = data.external_element_tags[keep]; empty!(data.external_element_tags); append!(data.external_element_tags, tags)
        if mesh.elementary_entities !== nothing
            el = mesh.elementary_entities[keep]; empty!(mesh.elementary_entities); append!(mesh.elementary_entities, el)
        end
    end
    # A surface lists its triangles in the order the classification collected
    # them: reorder each triangle block's columns to the records' order.
    position = Dict{UInt64,Int}()
    for (key, record) in sort!([(k, r) for (k, r) in _discrete_mesh_records_model(m) if k[1] == 2]; by=x -> x[1][2])
        for t in record.element_tags; position[t] = length(position) + 1; end
    end
    for (bi, block) in enumerate(mesh.blocks)
        (block isa ElementBlock && block.msh == 2 && bi <= length(data.external_element_tags)) || continue
        tags = data.external_element_tags[bi]
        all(t -> haskey(position, t), tags) || continue
        perm = sortperm(tags; by=t -> position[t])
        perm == eachindex(tags) && continue
        mesh.blocks[bi] = ElementBlock(block.msh, block.nodes[:, perm], block.tags[perm])
        data.external_element_tags[bi] = tags[perm]
        data.block_entities[bi] = data.block_entities[bi][perm]
        mesh.elementary_entities === nothing ||
            (mesh.elementary_entities[bi] = mesh.elementary_entities[bi][perm])
    end
    known = Set{UInt64}()
    for tags in data.external_element_tags; union!(known, tags); end
    record_of = Dict(_discrete_mesh_records_model(m))
    for dim in 0:1
        for key in sort!([k for k in keys(record_of) if k[1] == dim]; by=last)
            record = record_of[key]
            columns = [i for i in eachindex(record.element_tags) if !(record.element_tags[i] in known)]
            isempty(columns) && continue
            msh = dim == 0 ? 15 : 1
            all(i -> record.element_types[i] == msh, columns) || continue
            nodes = Matrix{Int32}(undef, dim + 1, length(columns))
            ok = true
            for (j, i) in enumerate(columns), (k, v) in enumerate(record.element_nodes[i])
                r = get(row, v, Int32(0)); r == 0 && (ok = false); nodes[k, j] = r
            end
            ok || continue
            push!(mesh.blocks, ElementBlock(msh, nodes, zeros(Int32, length(columns))))
            push!(data.block_entities, fill(Int32(key[2]), length(columns)))
            push!(data.external_element_tags, record.element_tags[columns])
            mesh.elementary_entities === nothing || push!(mesh.elementary_entities, fill(Int32(key[2]), length(columns)))
        end
    end
    return nothing
end
