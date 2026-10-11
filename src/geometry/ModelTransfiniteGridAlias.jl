# Only a compact TF constructor may expose its newly evaluated diagonal rows
# on the exterior. The literal lower Surface and actual Region retain separate
# UIDs and coordinates. This receipt maps their recorded logical grid slots;
# it does not merge either record or select identity by coordinate proximity.
struct _ModelTransfiniteGridAlias
    row::Int32
    source_uid::UInt64
    point::NTuple{3,Float64}
end

struct _ModelTransfiniteGridAliasReceipt
    counts::NTuple{3,Int}
    owner::Int32
    events::Vector{_ModelTransfiniteGridAlias}
end

function _model_transfinite_compact_alias_events(mesh,ids,counts,owner,caller)
    n,nv,nw=counts
    n==nv && 0<n<=size(mesh.coords,2) && 0<nw<=size(mesh.coords,2) &&
        owner>0 && Int128(length(ids))==(Int128(n)+1)*(Int128(nw)+1) ||
        throw(ArgumentError("$caller: compact diagonal alias layout is inconsistent"))
    lattice=(Int128(n)+1)*(Int128(n)+2)÷2;block=1+Int128(n)*(n+1)
    events=_ModelTransfiniteGridAlias[]
    for i in 1:n-1,k in 1:nw-1
        event=lattice+(k-1)*block+2+Int128(i-1)*(n+1)+i
        event<=min(size(mesh.coords,2),typemax(Int32)) || throw(ArgumentError(
            "$caller: compact diagonal event has no actual constructor row"))
        row=Int32(event)
        uid=UInt64(ids[i+1+k*(n+1)])
        uid>0 || throw(ArgumentError("$caller: compact diagonal event lacks its lower grid UID"))
        push!(events,_ModelTransfiniteGridAlias(row,uid,ntuple(axis->mesh.coords[axis,row],3)))
    end
    return _ModelTransfiniteGridAliasReceipt(counts,Int32(owner),events)
end

function _model_transfinite_compact_alias_binding(receipt,faces,mesh,source_tags,caller)
    receipt===nothing && return Dict{Int32,UInt64}()
    receipt isa _ModelTransfiniteGridAliasReceipt && mesh!==nothing &&
        size(mesh.coords,2)==length(source_tags) && length(faces)==5 || throw(ArgumentError(
        "$caller: compact diagonal aliases lack their actual five-face constructor"))
    found=findall(face->face.owner==receipt.owner,faces)
    length(found)==1 || throw(ArgumentError("$caller: compact diagonal alias owner is absent or repeated"))
    face=faces[only(found)];n,nv,nw=receipt.counts
    face.extents==(n,nw) || throw(ArgumentError(
        "$caller: compact diagonal alias face changes its grid extents"))
    expected=_model_transfinite_compact_alias_events(mesh,face.ids,receipt.counts,receipt.owner,caller)
    length(receipt.events)==length(expected.events) || throw(ArgumentError(
        "$caller: compact diagonal alias event list is incomplete"))
    aliases=Dict{Int32,UInt64}()
    for (actual,current) in zip(receipt.events,expected.events)
        actual.row==current.row && actual.source_uid==current.source_uid &&
            isequal(actual.point,current.point) && iszero(source_tags[actual.row]) || throw(ArgumentError(
            "$caller: compact diagonal alias changed its recorded row, UID or evaluated point"))
        haskey(aliases,actual.row) && throw(ArgumentError(
            "$caller: compact diagonal alias repeats an evaluated event"))
        aliases[actual.row]=actual.source_uid
    end
    return aliases
end
