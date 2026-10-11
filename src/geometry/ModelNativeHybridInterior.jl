# The embedded native recovery policy has a zero orient4d objective. The
# segment, face, physical-pool and inverse phases remain operation-owned;
# the generic convex Arena admission policy is independent of this route.
module ModelNativeHybridInterior
const A=parentmodule(@__MODULE__).ModelNativeTetgenArena
const R=parentmodule(@__MODULE__).ModelNativeHybridRecovery
const U=parentmodule(@__MODULE__).ModelNativeHybridUnflip
const L=parentmodule(@__MODULE__).ModelNativeHybridLawson
const C=parentmodule(@__MODULE__).ModelNativeHybridCarve

function currentedge(a,record)
    first,last=record.vertices[1:2]
    h=record.handle
    if a.cells[h.cell].alive
        (A.org(a,h),A.dest(a,h))==(first,last) && return h
        (A.org(a,h),A.dest(a,h))==(last,first) && return A.esym(h)
    end
    h,kind=A.finddirection!(a,A.point2tetorg(a,first),last)
    kind===:across_vertex && A.dest(a,h)==last ? h : A.NULL
end

function unchanged(a,saved)
    length(a.cells)==length(saved.cells) && a.dead==saved.dead &&
        a.hints==saved.hints || return false
    return all(i->a.cells[i].vertices==saved.cells[i].vertices &&
        a.cells[i].alive==saved.cells[i].alive && a.cells[i].links==saved.cells[i].links,
        eachindex(a.cells))
end

function quality!(a,hulls)
    queue=Tuple{Int32,NTuple{4,Int32},Vector{Float64}}[]
    work=A.QualityWorkspace()
    for (slot,cell) in enumerate(a.cells)
        cell.alive && !(Int32(slot) in hulls) || continue
        values=collect(A.quality_dihedrals(ntuple(i->a.points[cell.vertices[i]],4),work))
        minimum(values)<A.QUALITY_MAX_COS && push!(queue,(Int32(slot),cell.vertices,values))
    end
    flips=0
    while !isempty(queue)
        nextqueue=Tuple{Int32,NTuple{4,Int32},Vector{Float64}}[]
        for (slot,vertices,values) in queue
            a.cells[slot].alive && a.cells[slot].vertices==vertices || continue
            for edge in 1:6
                values[edge]<A.QUALITY_MAX_COS || continue
                handle=A.Handle(slot,A.QUALITY_EDGE_VERSIONS[edge])
                ring=try
                    A.spin(a,handle)
                catch error
                    error isa A.Unsupported || rethrow()
                    continue
                end
                length(ring)<=10 && all(h->!(h.cell in hulls),ring) || continue
                saved=deepcopy(a)
                constraint=R.QualityConstraint(values[edge],Ref(0.0))
                removed=length(ring)<=4 ? U.removeedge!(a,handle,constraint).removed :
                    R.removeedge!(a,handle,constraint)
                if removed
                    flips+=1
                    if constraint.output[]<A.QUALITY_MAX_COS
                        for (index,cell) in enumerate(a.cells)
                            cell.alive && !(Int32(index) in hulls) || continue
                            old=index<=length(saved.cells) ? saved.cells[index] : nothing
                            old!==nothing && old.alive && old.vertices==cell.vertices && continue
                            newvalues=collect(A.quality_dihedrals(
                                ntuple(i->a.points[cell.vertices[i]],4),work))
                            minimum(newvalues)<A.QUALITY_MAX_COS &&
                                push!(nextqueue,(Int32(index),cell.vertices,newvalues))
                        end
                    end
                    break
                elseif length(ring)>4 && !unchanged(a,saved)
                    throw(A.Unsupported("hybrid quality refusal requires a recorded larger-star inverse"))
                end
            end
        end
        queue=nextqueue
    end
    return flips
end

function recover(kernel,facets;seed=1,component=nothing,embedded_segments=NTuple{2,Int32}[])
    all(segment->length(segment)==2 && allunique(segment) &&
        all(row->0<row<=kernel.count,segment),embedded_segments) || throw(A.Unsupported(
        "hybrid embedded segments contain missing or duplicate current rows"))
    faces=NTuple{3,Int32}[Tuple(Int32.(face)) for face in facets]
    isempty(faces) && throw(A.Unsupported("hybrid recovery requires explicit constructor facets"))
    wanted=Set{NTuple{3,Int32}}();edges=Dict{Tuple{Int32,Int32},Int}()
    for face in faces
        allunique(face) && all(row->0<row<=kernel.count,face) || throw(A.Unsupported(
            "hybrid constructor facets contain missing or duplicate current rows"))
        key=A.facekey(face)
        key in wanted && throw(A.Unsupported("hybrid constructor facets repeat an actual panel"))
        push!(wanted,key)
        for side in 1:3
            edge=A.edgekey(face[side],face[mod1(side+1,3)])
            edges[edge]=get(edges,edge,0)+1
        end
    end
    multiplicities=if component===nothing
        all(==(2),values(edges)) || throw(A.Unsupported(
            "hybrid constructor facets are not a closed two-face manifold"))
        Dict(face=>1 for face in wanted)
    else
        parentmodule(@__MODULE__)._model_native_hybrid_recovery_multiplicities(component,kernel,faces)
    end
    a=A.Arena(kernel;seed)
    initial_hulls=C.initial_hulls!(a)
    protocol=R.recoverboundary!(a,faces;embedded_segments)
    A.consistency(a) || throw(A.Unsupported("hybrid recovered face bonds are not reciprocal"))
    carving=C.carve!(a,initial_hulls)
    hulls=Set(carving.hulls)
    lawson=L.lawson!(a,hulls)
    for record in lawson.unflippable
        handle=currentedge(a,record)
        handle==A.NULL && continue
        # Native's embedded orient4d is zero, so every root objective removal
        # is rejected and unflipped. An accepted removal violates that policy.
        result=U.removeedge!(a,handle,0.0)
        result.removed && throw(A.Unsupported(
            "hybrid embedded zero objective accepted a root Delaunay removal"))
    end
    quality=quality!(a,hulls)
    A.consistency(a) || throw(A.Unsupported("hybrid final face bonds are not reciprocal"))
    cells=NTuple{4,Int32}[];slots=Int32[]
    counts=Dict{NTuple{3,Int32},Int}()
    for (slot,cell) in enumerate(a.cells)
        cell.alive && !(Int32(slot) in hulls) || continue
        all(row->0<row<=kernel.count,cell.vertices) && allunique(cell.vertices) ||
            throw(A.Unsupported("hybrid final cell retains an enclosing-box or duplicate row"))
        A.orient(a,cell.vertices...)<0 || throw(A.Unsupported(
            "hybrid temporary cell loses native positive orientation"))
        push!(cells,cell.vertices);push!(slots,Int32(slot))
        for omit in 1:4
            face=A.facekey(Tuple(cell.vertices[i] for i in 1:4 if i!=omit))
            counts[face]=get(counts,face,0)+1
        end
    end
    all(count->count<=2,values(counts)) || throw(A.Unsupported(
        "hybrid final constructor cells are not a manifold"))
    observed=Set(face for (face,count) in counts if count==1)
    exterior=Set(face for (face,count) in multiplicities if count==1)
    observed==exterior && all(get(counts,face,0)==count for (face,count) in multiplicities) || throw(A.Unsupported(
        "hybrid final cell exterior differs from the owned constructor facets"))
    for segment in embedded_segments
        A.edgekey(segment...) in a.protected_segments &&
            any(cell->all(row->row in cell,segment),cells) || throw(A.Unsupported(
            "hybrid recovery lost a protected embedded Curve segment"))
    end
    return (cells=cells,pool_slots=slots,arena=a,protocol=protocol,carving=carving,
        lawson=lawson,quality_flips=quality)
end
end
