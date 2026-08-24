"""
    Post

List-based post-processing views (Gmsh `.pos` analogue). A view is an owned,
named scalar field sampled at mesh nodes. Plugins are registered functions of a
view; this compact interface does not yet claim the full Gmsh view data model.
"""
module Post

using ..MeshTypes: Mesh, nnodes

export View, view_value, add_plugin!, apply_plugin

"""
    View(name, mesh, values)
    View(name, coords, values)

Owned scalar samples at three-dimensional coordinates. `coords` must be a
finite `3 × n` real matrix and `values` a finite length-`n` real vector. Inputs
are converted to `Float64` and copied, so later mutation of caller arrays does
not change the view.
"""
struct View
    name::String
    coords::Matrix{Float64}
    values::Vector{Float64}
    function View(name::AbstractString,coords::AbstractMatrix{<:Real},
                  values::AbstractVector{<:Real})
        caller="View"
        size(coords,1)==3 || throw(ArgumentError(
            "$caller: coords must be a 3 × n matrix"))
        size(coords,2)==length(values) || throw(ArgumentError(
            "$caller: value count must match coordinate count"))
        coordinates=try
            Matrix{Float64}(coords)
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError(
                "$caller: coordinates must be Float64-representable: "*
                sprint(showerror,err)))
        end
        samples=try
            Vector{Float64}(values)
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError(
                "$caller: values must be Float64-representable: "*
                sprint(showerror,err)))
        end
        @inbounds for node in axes(coordinates,2),component in 1:3
            isfinite(coordinates[component,node]) || throw(ArgumentError(
                "$caller: coordinate $component of node $node must be finite"))
        end
        @inbounds for sample in eachindex(samples)
            isfinite(samples[sample]) || throw(ArgumentError(
                "$caller: value $sample must be finite"))
        end
        new(String(name),coordinates,samples)
    end
end

function View(name::AbstractString, mesh::Mesh, values::AbstractVector{<:Real})
    nnodes(mesh)==length(values) || throw(ArgumentError(
        "View: value count must match mesh nodes"))
    return View(name,mesh.coords,values)
end

"""Return scalar sample `i` from `v`, rejecting Boolean and out-of-range indices."""
function view_value(v::View, i::Integer)
    i isa Bool && throw(ArgumentError("view_value: index must not be Bool"))
    index=try
        Int(i)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("view_value: index exceeds the platform Int range"))
    end
    1<=index<=length(v.values) || throw(ArgumentError(
        "view_value: index $index is outside 1:$(length(v.values))"))
    return @inbounds v.values[index]
end

const PLUGINS = Dict{String,Function}()
const PLUGIN_LOCK = ReentrantLock()

"""
    add_plugin!(name, callback) -> String

Register `callback(view)` under a nonempty `name`, replacing any callback with
the same name. Registration is synchronized; callback execution is not held
under the registry lock.
"""
function add_plugin!(name::AbstractString, f::Function)
    isempty(name) && throw(ArgumentError("add_plugin!: empty plugin name"))
    plugin_name=String(name)
    lock(PLUGIN_LOCK) do
        PLUGINS[plugin_name]=f
    end
    return plugin_name
end

"""Apply the plugin registered as `name` to `v`, propagating its return value or error."""
function apply_plugin(name::AbstractString, v::View)
    plugin_name=String(name)
    f=lock(PLUGIN_LOCK) do
        get(PLUGINS,plugin_name,nothing)
    end
    f===nothing && throw(ArgumentError("apply_plugin: unknown plugin $name"))
    return f(v)
end

# Built-in: Scale plugin multiplies all samples by a constant stored as the name
# "Scale". The factor is passed by creating a closure via add_plugin!.
add_plugin!("Abs", v -> View(v.name*"_abs", v.coords, abs.(v.values)))
add_plugin!("IsosurfaceZeroCount", v -> count(iszero, v.values))

end # module
