"""
    Post

List-based post-processing views (Gmsh `.pos` analogue). A view is a named
scalar field sampled at mesh nodes. Plugins are registered functions of a view.
"""
module Post

using ..MeshTypes: Mesh, nnodes

export View, view_value, add_plugin!, apply_plugin

struct View
    name::String
    coords::Matrix{Float64}
    values::Vector{Float64}
end

function View(name::AbstractString, mesh::Mesh, values::AbstractVector{<:Real})
    nnodes(mesh)==length(values) || throw(ArgumentError("View: value count must match nodes"))
    V=Float64[Float64(v) for v in values]
    all(isfinite,V) || throw(ArgumentError("View: values must be finite"))
    return View(String(name), copy(mesh.coords), V)
end

function view_value(v::View, i::Integer)
    1<=i<=length(v.values) || throw(ArgumentError("view_value: index out of range"))
    return v.values[i]
end

const PLUGINS = Dict{String,Function}()

function add_plugin!(name::AbstractString, f::Function)
    isempty(name) && throw(ArgumentError("add_plugin!: empty plugin name"))
    PLUGINS[String(name)]=f
    return name
end

function apply_plugin(name::AbstractString, v::View)
    f=get(PLUGINS,String(name),nothing)
    f===nothing && throw(ArgumentError("apply_plugin: unknown plugin $name"))
    return f(v)
end

# Built-in: Scale plugin multiplies all samples by a constant stored as the name
# "Scale". The factor is passed by creating a closure via add_plugin!.
add_plugin!("Abs", v -> View(v.name*"_abs", v.coords, abs.(v.values)))
add_plugin!("IsosurfaceZeroCount", v -> count(iszero, v.values))

end # module
