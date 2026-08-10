"""
    SizeField

Mesh-size (edge-length target) fields (PLAN.md §3 "SizeField"). A size field maps
a point to a desired local element size `h`. Stage 2 needs only the basics used to
drive 1-D and surface meshing; Stage 4 adds curvature/distance/background fields.

`size_at(sf, x, y, z)` (and the 3-tuple form) returns `h > 0`.
"""
module SizeField

export AbstractSizeField, ConstantSize, FunctionSize, MinSize, size_at

abstract type AbstractSizeField end

"""Uniform target size `h` everywhere."""
struct ConstantSize <: AbstractSizeField
    h::Float64
    function ConstantSize(h::Real)
        h > 0 || throw(ArgumentError("ConstantSize: h must be positive"))
        new(Float64(h))
    end
end

"""Size from an arbitrary callable `f(x,y,z) -> h`."""
struct FunctionSize{F} <: AbstractSizeField
    f::F
end

"""Pointwise minimum of several fields (the meshing-conservative combination)."""
struct MinSize <: AbstractSizeField
    fields::Vector{AbstractSizeField}
    function MinSize(fields::AbstractVector{<:AbstractSizeField})
        isempty(fields) && throw(ArgumentError("MinSize: need at least one field"))
        new(collect(AbstractSizeField, fields))
    end
end

@inline size_at(sf::ConstantSize, x, y, z) = sf.h
@inline function size_at(sf::FunctionSize, x, y, z)
    h = sf.f(x, y, z)
    (h > 0 && isfinite(h)) || throw(ArgumentError("FunctionSize: f returned non-positive/NaN size $h at ($x,$y,$z)"))
    return Float64(h)
end
function size_at(sf::MinSize, x, y, z)
    m = Inf
    @inbounds for f in sf.fields
        h = size_at(f, x, y, z); h < m && (m = h)
    end
    return m
end

@inline size_at(sf::AbstractSizeField, p) = size_at(sf, p[1], p[2], length(p) >= 3 ? p[3] : 0.0)

end # module SizeField
