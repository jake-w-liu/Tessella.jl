"""
    GUI

Headless GUI/state machine for Tessella. Commands mutate a `GuiState`
(model, mesh, camera, selection, log). The interactive window launcher is
optional; its logic is unit-tested here. Pixel FLTK chrome is a non-goal.
"""
module GUI

using ..API
using ..MeshTypes: Mesh, nnodes, ntets

export GuiState, gui_command!, gui_new!, camera_reset!

mutable struct GuiState
    initialized::Bool
    camera::NTuple{3,Float64}
    selected::Vector{Int}
    log::Vector{String}
    mesh_nodes::Int
    mesh_tets::Int
end

GuiState() = GuiState(false,(0.0,0.0,1.0),Int[],String[],0,0)

function _log!(s::GuiState, msg::AbstractString)
    push!(s.log,String(msg))
    return nothing
end

function gui_new!(s::GuiState)
    API.initialize()
    s.initialized=true
    s.selected=Int[]
    s.mesh_nodes=0; s.mesh_tets=0
    camera_reset!(s)
    _log!(s,"new model")
    return s
end

function camera_reset!(s::GuiState)
    s.camera=(0.0,0.0,1.0)
    return s.camera
end

function gui_command!(s::GuiState, cmd::AbstractString)
    parts=split(strip(cmd))
    isempty(parts) && throw(ArgumentError("gui_command!: empty command"))
    head=lowercase(parts[1])
    head=="new" && return gui_new!(s)
    s.initialized || throw(ArgumentError("gui_command!: call new first"))
    if head=="box" && length(parts)==8
        nums=parse.(Float64, parts[2:7]); tag=parse(Int,parts[8])
        API.model.add_box(nums[1],nums[2],nums[3],nums[4],nums[5],nums[6]; tag=tag)
        _log!(s,"box $tag")
    elseif head=="mesh" && length(parts)==2
        dim=parse(Int,parts[2])
        mesh=API.mesh.generate(dim)
        s.mesh_nodes=nnodes(mesh); s.mesh_tets=ntets(mesh)
        _log!(s,"mesh dim=$dim nodes=$(s.mesh_nodes) tets=$(s.mesh_tets)")
    elseif head=="select" && length(parts)>=2
        s.selected=parse.(Int, parts[2:end])
        _log!(s,"select $(s.selected)")
    elseif head=="camera" && length(parts)==4
        s.camera=(parse(Float64,parts[2]),parse(Float64,parts[3]),parse(Float64,parts[4]))
        _log!(s,"camera $(s.camera)")
    elseif head=="reset"
        camera_reset!(s); _log!(s,"camera reset")
    else
        throw(ArgumentError("gui_command!: unknown command $cmd"))
    end
    return s
end

end # module
