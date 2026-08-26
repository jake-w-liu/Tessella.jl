"""
    Periodic

Certification and exact coordinate snapping for equal-count master/slave node
sets related by a finite translation or nonsingular affine transformation. Node
numbering, connectivity, and tags are preserved. The caller retains the
node-pair map and transformation. These compact `Mesh` operations do not add
entity metadata; [`Tessella.Model.set_periodic!`](@ref) owns supported native
straight-curve relations separately.
"""
module Periodic

using ..MeshTypes: Mesh, nnodes, validate
using ..Transform: _affine_coordinate, _transform_homogeneous

export periodic_identify, periodic_identify_affine

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

function _periodic_tolerance(atol,caller::AbstractString)
    atol isa Bool && throw(ArgumentError("$caller: atol must not be Bool"))
    atol isa Real || throw(ArgumentError("$caller: atol must be real"))
    tolerance=try
        Float64(atol)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: atol must be Float64-representable"))
    end
    (isfinite(tolerance) && tolerance>=0) || throw(ArgumentError(
        "$caller: atol must be finite and non-negative"))
    return tolerance
end

function _periodic_pairs(mesh::Mesh,master::AbstractVector{<:Integer},
                         slave::AbstractVector{<:Integer},caller::AbstractString)
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
    return masters,slaves
end

function _periodic_output(mesh::Mesh,slaves::Vector{Int},
                          expected_coordinates::Vector{NTuple{3,Float64}},
                          caller::AbstractString)
    coords=copy(mesh.coords)
    for i in eachindex(slaves)
        ib=slaves[i];expected=expected_coordinates[i]
        coords[1,ib]=expected[1];coords[2,ib]=expected[2];coords[3,ib]=expected[3]
    end
    out=Mesh(coords; segs=copy(mesh.segs), tris=copy(mesh.tris), tets=copy(mesh.tets),
             seg_tag=copy(mesh.seg_tag), tri_tag=copy(mesh.tri_tag),
             tet_tag=copy(mesh.tet_tag))
    diagnostic=validate(out)
    diagnostic.ok || throw(ErrorException(
        "$caller: invalid mesh — "*join(diagnostic.messages,"; ")))
    return out
end

function _periodic_affine(raw,caller::AbstractString)
    coefficients,translation,_=_transform_homogeneous(raw,caller)
    return coefficients,translation
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
    tol=_periodic_tolerance(atol,caller)
    masters,slaves=_periodic_pairs(mesh,master,slave,caller)

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
    return _periodic_output(mesh,slaves,expected_coordinates,caller)
end

"""
    periodic_identify_affine(mesh, affine, master, slave; atol=1e-12) -> Mesh

Validate a one-to-one correspondence between disjoint `master` and `slave` node
sets satisfying `slave ≈ affine(master)`, then snap every slave coordinate to
the transformed master coordinate. `affine` is either a finite 4×4 matrix or a
16-entry tuple/vector in Gmsh row-major order. Its homogeneous row must be
`(0, 0, 0, 1)` and its 3×3 linear part must be nonsingular.

The input mesh is validated and left unchanged; node numbering, simplex
connectivity, and cell tags are copied without remapping. `atol` is an absolute
Euclidean-distance tolerance. The returned `Mesh` stores corrected geometry,
not persistent periodic entity metadata; callers that need solver constraints
must retain the supplied pairs and affine transform.
"""
function periodic_identify_affine(mesh::Mesh,affine,
                                  master::AbstractVector{<:Integer},
                                  slave::AbstractVector{<:Integer};
                                  atol::Real=1e-12)
    caller="periodic_identify_affine"
    input_diagnostic=validate(mesh)
    input_diagnostic.ok || throw(ArgumentError(
        "$caller: input mesh is invalid — "*join(input_diagnostic.messages,"; ")))
    coefficients,translation=_periodic_affine(affine,caller)
    tolerance=_periodic_tolerance(atol,caller)
    masters,slaves=_periodic_pairs(mesh,master,slave,caller)
    a11,a21,a31,a12,a22,a32,a13,a23,a33=coefficients
    expected_coordinates=Vector{NTuple{3,Float64}}(undef,length(masters))
    for i in eachindex(masters)
        ia=masters[i];ib=slaves[i]
        x=mesh.coords[1,ia];y=mesh.coords[2,ia];z=mesh.coords[3,ia]
        expected=(
            _affine_coordinate(0.0,translation[1],a11,a12,a13,x,y,z,
                               0.0,0.0,0.0,caller,ia),
            _affine_coordinate(0.0,translation[2],a21,a22,a23,x,y,z,
                               0.0,0.0,0.0,caller,ia),
            _affine_coordinate(0.0,translation[3],a31,a32,a33,x,y,z,
                               0.0,0.0,0.0,caller,ia),
        )
        got=(mesh.coords[1,ib],mesh.coords[2,ib],mesh.coords[3,ib])
        hypot(expected[1]-got[1],expected[2]-got[2],expected[3]-got[3])<=tolerance ||
            throw(ArgumentError(
                "$caller: slave $ib is not affine(master $ia)"))
        expected_coordinates[i]=expected
    end
    return _periodic_output(mesh,slaves,expected_coordinates,caller)
end

end # module
