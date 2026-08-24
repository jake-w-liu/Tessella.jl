"""
    GUI

Headless GUI/state machine for Tessella. Commands mutate a `GuiState`
(model, mesh, camera, selection, log). The interactive window launcher is
optional; its logic is unit-tested here. Pixel FLTK chrome is a non-goal.
"""
module GUI

using ..API
using ..MeshTypes: nnodes, ntets

export GuiState, gui_command!, gui_new!, camera_reset!

const _MAX_COMMAND_BYTES=1_000_000
const _MAX_COMMAND_TOKENS=10_000
const _MAX_SELECTION=10_000
const _MAX_LOG_ENTRIES=100_000

"""
    GuiState()
    GuiState(initialized, camera, selected, log, mesh_nodes, mesh_tets)

Owned state for the headless GUI command machine. Camera values must be finite,
selection tags must be unique positive integers, and mesh counts must be
nonnegative. Caller-provided selection and log vectors are copied.
"""
mutable struct GuiState
    initialized::Bool
    camera::NTuple{3,Float64}
    selected::Vector{Int}
    log::Vector{String}
    mesh_nodes::Int
    mesh_tets::Int

    function GuiState(initialized::Bool,camera::Tuple{Vararg{Real,3}},
                      selected::AbstractVector{<:Integer},
                      log::AbstractVector{<:AbstractString},
                      mesh_nodes::Integer,mesh_tets::Integer)
        any(x->x isa Bool,camera) && throw(ArgumentError(
            "GuiState: camera values must not be Bool"))
        cam=try
            (Float64(camera[1]),Float64(camera[2]),Float64(camera[3]))
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("GuiState: camera values must be Float64-representable"))
        end
        all(isfinite,cam) || throw(ArgumentError("GuiState: camera values must be finite"))
        length(selected)<=_MAX_SELECTION || throw(ArgumentError(
            "GuiState: selection exceeds $_MAX_SELECTION tags"))
        ids=Int[];sizehint!(ids,length(selected))
        for raw in selected
            raw isa Bool && throw(ArgumentError("GuiState: selection tags must not be Bool"))
            id=try
                Int(raw)
            catch err
                err isa InterruptException && rethrow()
                throw(ArgumentError("GuiState: selection tag exceeds the platform Int range"))
            end
            id>0 || throw(ArgumentError("GuiState: selection tags must be positive"))
            push!(ids,id)
        end
        length(unique(ids))==length(ids) || throw(ArgumentError(
            "GuiState: selection tags must be unique"))
        length(log)<=_MAX_LOG_ENTRIES || throw(ArgumentError(
            "GuiState: log exceeds $_MAX_LOG_ENTRIES entries"))
        entries=String.(log)
        nodes=_nonnegative_count(mesh_nodes,"mesh_nodes")
        tets=_nonnegative_count(mesh_tets,"mesh_tets")
        new(initialized,cam,ids,entries,nodes,tets)
    end
end

GuiState() = GuiState(false,(0.0,0.0,1.0),Int[],String[],0,0)

function _nonnegative_count(raw::Integer,what::AbstractString)
    raw isa Bool && throw(ArgumentError("GuiState: $what must not be Bool"))
    value=try
        Int(raw)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("GuiState: $what exceeds the platform Int range"))
    end
    value>=0 || throw(ArgumentError("GuiState: $what must be nonnegative"))
    return value
end

function _validate_state(s::GuiState)
    all(isfinite,s.camera) || throw(ArgumentError("GuiState: camera values must be finite"))
    (s.mesh_nodes>=0 && s.mesh_tets>=0) || throw(ArgumentError(
        "GuiState: mesh counts must be nonnegative"))
    length(s.selected)<=_MAX_SELECTION || throw(ArgumentError(
        "GuiState: selection exceeds $_MAX_SELECTION tags"))
    all(>(0),s.selected) || throw(ArgumentError("GuiState: selection tags must be positive"))
    length(unique(s.selected))==length(s.selected) || throw(ArgumentError(
        "GuiState: selection tags must be unique"))
    length(s.log)<=_MAX_LOG_ENTRIES || throw(ArgumentError(
        "GuiState: log exceeds $_MAX_LOG_ENTRIES entries"))
    return nothing
end

function _reserve_log(s::GuiState)
    length(s.log)<_MAX_LOG_ENTRIES || throw(ArgumentError(
        "GuiState: log limit of $_MAX_LOG_ENTRIES entries reached"))
    return nothing
end

function _log!(s::GuiState, msg::AbstractString)
    push!(s.log,String(msg))
    return nothing
end

"""Reset `s` and the process-global [`API`](@ref) session to an empty model."""
function gui_new!(s::GuiState)
    _reserve_log(s)
    API.initialize()
    s.initialized=true
    s.selected=Int[]
    s.mesh_nodes=0; s.mesh_tets=0
    camera_reset!(s)
    _log!(s,"new model")
    return s
end

"""Reset the GUI camera to `(0.0, 0.0, 1.0)` and return that tuple."""
function camera_reset!(s::GuiState)
    s.camera=(0.0,0.0,1.0)
    return s.camera
end

function _parse_float(token::AbstractString,what::AbstractString)
    value=try
        parse(Float64,token)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("gui_command!: invalid $what $(repr(String(token)))"))
    end
    isfinite(value) || throw(ArgumentError("gui_command!: $what must be finite"))
    return value
end

function _parse_int(token::AbstractString,what::AbstractString)
    return try
        parse(Int,token)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("gui_command!: invalid $what $(repr(String(token)))"))
    end
end

function _arity(parts,expected::Int,head::AbstractString)
    length(parts)==expected || throw(ArgumentError(
        "gui_command!: $head expects $(expected-1) arguments, got $(length(parts)-1)"))
    return nothing
end

"""
    gui_command!(state, command) -> GuiState

Apply one bounded command: `new`, `box x y z dx dy dz tag`, `mesh 2|3`,
`select tag...`, `camera x y z`, or `reset`. A command validates all of its
arguments before changing GUI state; geometry and meshing delegate to [`API`](@ref).
"""
function gui_command!(s::GuiState, cmd::AbstractString)
    ncodeunits(cmd)<=_MAX_COMMAND_BYTES || throw(ArgumentError(
        "gui_command!: command exceeds $_MAX_COMMAND_BYTES bytes"))
    parts=split(strip(cmd);limit=_MAX_COMMAND_TOKENS+1)
    isempty(parts) && throw(ArgumentError("gui_command!: empty command"))
    length(parts)<=_MAX_COMMAND_TOKENS || throw(ArgumentError(
        "gui_command!: command exceeds $_MAX_COMMAND_TOKENS tokens"))
    head=lowercase(parts[1])
    if head=="new"
        _arity(parts,1,head)
        return gui_new!(s)
    end
    s.initialized || throw(ArgumentError("gui_command!: call new first"))
    _validate_state(s)
    if head=="box"
        _arity(parts,8,head);_reserve_log(s)
        nums=[_parse_float(parts[i],"box value") for i in 2:7]
        tag=_parse_int(parts[8],"box tag")
        actual=API.model.add_box(nums[1],nums[2],nums[3],nums[4],nums[5],nums[6];tag=tag)
        _log!(s,"box $actual")
    elseif head=="mesh"
        _arity(parts,2,head);_reserve_log(s)
        dim=_parse_int(parts[2],"mesh dimension")
        mesh=API.mesh.generate(dim)
        s.mesh_nodes=nnodes(mesh);s.mesh_tets=ntets(mesh)
        _log!(s,"mesh dim=$dim nodes=$(s.mesh_nodes) tets=$(s.mesh_tets)")
    elseif head=="select"
        length(parts)>=2 || throw(ArgumentError("gui_command!: select expects at least one tag"))
        length(parts)-1<=_MAX_SELECTION || throw(ArgumentError(
            "gui_command!: selection exceeds $_MAX_SELECTION tags"))
        _reserve_log(s)
        selected=[_parse_int(parts[i],"selection tag") for i in 2:length(parts)]
        all(>(0),selected) || throw(ArgumentError(
            "gui_command!: selection tags must be positive"))
        length(unique(selected))==length(selected) || throw(ArgumentError(
            "gui_command!: selection tags must be unique"))
        s.selected=selected
        _log!(s,"select $(s.selected)")
    elseif head=="camera"
        _arity(parts,4,head);_reserve_log(s)
        s.camera=(_parse_float(parts[2],"camera coordinate"),
                  _parse_float(parts[3],"camera coordinate"),
                  _parse_float(parts[4],"camera coordinate"))
        _log!(s,"camera $(s.camera)")
    elseif head=="reset"
        _arity(parts,1,head);_reserve_log(s)
        camera_reset!(s);_log!(s,"camera reset")
    else
        throw(ArgumentError("gui_command!: unknown command $cmd"))
    end
    return s
end

end # module
