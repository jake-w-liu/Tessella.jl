# OCC shapes and GModel entities have separate identities and lifetimes.
# A public Point mutation does not change an OCC control pole. Shape IDs are
# operation-owned references, including separate IDs for coincident vertices.
struct _ModelOCCShape
    dimension::Int
    geometry::NamedTuple
    children::Vector{UInt64}
end

Base.:(==)(a::_ModelOCCShape,b::_ModelOCCShape)=
    a.dimension==b.dimension && _model_occ_geometry_equal(a.geometry,b.geometry) &&
    a.children==b.children

mutable struct _ModelNativeOCCKernel
    bindings::Dict{Tuple{Int,Int},UInt64}
    shape_tags::Dict{UInt64,Tuple{Int,Int}}
    shapes::Dict{UInt64,_ModelOCCShape}
    entity_maps::NTuple{6,_ModelOCCEntityMap}
    maxima::Vector{Int}
    next_shape::UInt64
    removed::Set{Tuple{Int,Int}}
    changed::Bool
    # Synchronize is idempotent: a bound shell maps to one GModel surface
    # loop for its whole lifetime, so repeated syncs reuse the slot.
    shell_loops::Dict{UInt64,Int}
    # tag => the kernel point value last imported into `m.points` by
    # synchronize. Lets the dim-0 arm distinguish a kernel-side rewrite
    # (OCC transform — push to public) from a public `set_coordinates`
    # move (kernel record untouched — the public value wins), mirroring
    # upstream where both names land on the same vertex object.
    mirrored_points::Dict{Int,NTuple{3,Float64}}
end

_ModelNativeOCCKernel()=_ModelNativeOCCKernel(
    Dict{Tuple{Int,Int},UInt64}(),Dict{UInt64,Tuple{Int,Int}}(),
    Dict{UInt64,_ModelOCCShape}(),ntuple(_->_ModelOCCEntityMap(),6),zeros(Int,6),UInt64(1),
    Set{Tuple{Int,Int}}(),false,Dict{UInt64,Int}(),
    Dict{Int,NTuple{3,Float64}}())

Base.:(==)(a::_ModelNativeOCCKernel,b::_ModelNativeOCCKernel)=
    all(name->getfield(a,name)==getfield(b,name),fieldnames(_ModelNativeOCCKernel))
