"""
    Periodic

Periodic identification of equal-count master/slave node sets related by a
finite translation. Identified slave nodes are rewritten onto the translated
master coordinates; connectivity is remapped. Embedded constraints are the
existing CDT segment/triangle cells already stored on a `Mesh`.
"""
module Periodic

using ..MeshTypes: Mesh, nnodes, validate

export periodic_identify

function _finite3(raw, caller)
    length(raw)==3 || throw(ArgumentError("$caller: translation must have 3 components"))
    p=try (Float64(raw[1]),Float64(raw[2]),Float64(raw[3])) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: translation must be Float64-representable"))
    end
    all(isfinite,p) || throw(ArgumentError("$caller: translation must be finite"))
    return p
end

function periodic_identify(mesh::Mesh, translation, master::AbstractVector{<:Integer},
                           slave::AbstractVector{<:Integer}; atol::Real=1e-12)
    caller="periodic_identify"
    t=_finite3(translation,caller)
    hypot(t...)>0 || throw(ArgumentError("$caller: translation must be nonzero"))
    tol=try Float64(atol) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: atol must be Float64-representable"))
    end
    (isfinite(tol) && tol>=0) || throw(ArgumentError("$caller: atol must be non-negative"))
    length(master)==length(slave) || throw(ArgumentError(
        "$caller: master and slave node counts differ"))
    isempty(master) && throw(ArgumentError("$caller: need at least one identified pair"))
    nn=nnodes(mesh)
    coords=copy(mesh.coords)
    used=falses(nn)
    for (a,b) in zip(master,slave)
        ia,ib=Int(a),Int(b)
        (1<=ia<=nn && 1<=ib<=nn) || throw(ArgumentError("$caller: node index out of range"))
        ia==ib && throw(ArgumentError("$caller: master and slave must be distinct"))
        used[ib] && throw(ArgumentError("$caller: slave node $ib is paired twice"))
        used[ib]=true
        expected=(coords[1,ia]+t[1], coords[2,ia]+t[2], coords[3,ia]+t[3])
        got=(coords[1,ib],coords[2,ib],coords[3,ib])
        hypot(expected[1]-got[1],expected[2]-got[2],expected[3]-got[3])<=tol ||
            throw(ArgumentError("$caller: slave $ib is not master $ia + translation"))
        coords[1,ib]=expected[1]; coords[2,ib]=expected[2]; coords[3,ib]=expected[3]
    end
    out=Mesh(coords; segs=copy(mesh.segs), tris=copy(mesh.tris), tets=copy(mesh.tets),
             seg_tag=copy(mesh.seg_tag), tri_tag=copy(mesh.tri_tag), tet_tag=copy(mesh.tet_tag))
    diag=validate(out)
    diag.ok || throw(ErrorException("$caller: invalid mesh — "*join(diag.messages,"; ")))
    return out
end

end # module
