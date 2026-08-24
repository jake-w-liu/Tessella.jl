"""
    Periodic

Certification and exact coordinate snapping for equal-count master/slave node
sets related by a finite translation. Node numbering, connectivity, and tags are
preserved. The caller retains the node-pair map; this compact `Mesh` operation
does not claim Gmsh model-level periodic entity metadata.
"""
module Periodic

using ..MeshTypes: Mesh, nnodes, validate

export periodic_identify

function _periodic_translation(raw,caller)
    (raw isa Tuple || raw isa AbstractVector) || throw(ArgumentError(
        "$caller: translation must be a tuple or vector with 3 components"))
    length(raw)==3 || throw(ArgumentError(
        "$caller: translation must have 3 components"))
    p=ntuple(3) do i
        value=raw[i]
        value isa Bool && throw(ArgumentError(
            "$caller: translation component $i must not be Bool"))
        value isa Real || throw(ArgumentError(
            "$caller: translation component $i must be real"))
        result=try
            Float64(value)
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError(
                "$caller: translation component $i must be Float64-representable"))
        end
        isfinite(result) || throw(ArgumentError(
            "$caller: translation component $i must be finite"))
        result
    end
    return p
end

function _periodic_index(value,nn::Int,caller,name,index)
    value isa Bool && throw(ArgumentError(
        "$caller: $name index $index must not be Bool"))
    converted=try
        Int(value)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "$caller: $name index $index exceeds the platform Int range"))
    end
    1<=converted<=nn || throw(ArgumentError(
        "$caller: $name node $converted is outside 1:$nn"))
    return converted
end

"""
    periodic_identify(mesh, translation, master, slave; atol=1e-12) -> Mesh

Validate a one-to-one correspondence between disjoint `master` and `slave` node
sets satisfying `slave ≈ master + translation`, then snap each slave coordinate
to that translated master coordinate. The input mesh is validated and left
unchanged; node numbering, simplex connectivity, and cell tags are copied
without remapping. `atol` is an absolute Euclidean-distance tolerance.

The returned `Mesh` stores the corrected geometry, not a persistent periodic
entity relation. Callers that need solver constraints must retain the supplied
node pairs separately.
"""
function periodic_identify(mesh::Mesh, translation, master::AbstractVector{<:Integer},
                           slave::AbstractVector{<:Integer}; atol::Real=1e-12)
    caller="periodic_identify"
    input_diagnostic=validate(mesh)
    input_diagnostic.ok || throw(ArgumentError(
        "$caller: input mesh is invalid — "*join(input_diagnostic.messages,"; ")))
    t=_periodic_translation(translation,caller)
    magnitude=hypot(t...)
    (isfinite(magnitude) && magnitude>0) || throw(ArgumentError(
        "$caller: translation must have finite positive length"))
    atol isa Bool && throw(ArgumentError("$caller: atol must not be Bool"))
    tol=try Float64(atol) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: atol must be Float64-representable"))
    end
    (isfinite(tol) && tol>=0) || throw(ArgumentError("$caller: atol must be non-negative"))
    length(master)==length(slave) || throw(ArgumentError(
        "$caller: master and slave node counts differ"))
    isempty(master) && throw(ArgumentError("$caller: need at least one identified pair"))
    nn=nnodes(mesh)
    masters=Vector{Int}(undef,length(master));slaves=similar(masters)
    used_master=falses(nn);used_slave=falses(nn)
    for (i,(master_node,slave_node)) in enumerate(zip(master,slave))
        ia=_periodic_index(master_node,nn,caller,"master",i)
        ib=_periodic_index(slave_node,nn,caller,"slave",i)
        used_master[ia] && throw(ArgumentError(
            "$caller: master node $ia is paired more than once"))
        used_slave[ib] && throw(ArgumentError(
            "$caller: slave node $ib is paired more than once"))
        used_master[ia]=true;used_slave[ib]=true
        masters[i]=ia;slaves[i]=ib
    end
    overlap=findfirst(used_master .& used_slave)
    overlap===nothing || throw(ArgumentError(
        "$caller: node $overlap belongs to both master and slave sets"))

    expected_coordinates=Vector{NTuple{3,Float64}}(undef,length(masters))
    for i in eachindex(masters)
        ia=masters[i];ib=slaves[i]
        expected=(mesh.coords[1,ia]+t[1],mesh.coords[2,ia]+t[2],
                  mesh.coords[3,ia]+t[3])
        all(isfinite,expected) || throw(ArgumentError(
            "$caller: translated master node $ia is not Float64-representable"))
        got=(mesh.coords[1,ib],mesh.coords[2,ib],mesh.coords[3,ib])
        hypot(expected[1]-got[1],expected[2]-got[2],expected[3]-got[3])<=tol ||
            throw(ArgumentError("$caller: slave $ib is not master $ia + translation"))
        expected_coordinates[i]=expected
    end
    coords=copy(mesh.coords)
    for i in eachindex(slaves)
        ib=slaves[i];expected=expected_coordinates[i]
        coords[1,ib]=expected[1];coords[2,ib]=expected[2];coords[3,ib]=expected[3]
    end
    out=Mesh(coords; segs=copy(mesh.segs), tris=copy(mesh.tris), tets=copy(mesh.tets),
             seg_tag=copy(mesh.seg_tag), tri_tag=copy(mesh.tri_tag), tet_tag=copy(mesh.tet_tag))
    diag=validate(out)
    diag.ok || throw(ErrorException("$caller: invalid mesh — "*join(diag.messages,"; ")))
    return out
end

end # module
