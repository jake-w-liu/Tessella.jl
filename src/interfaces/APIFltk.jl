# `fltk` module: the FLTK user-interface surface, matching the gmsh 4.15.2
# `gmshFltk*` API names. Tessella ships no FLTK backend, so the headless
# contract is deterministic:
#   * `is_available()` returns `0` — the upstream build-time capability flag —
#     which is honest both for the missing backend and for headless runs.
#   * All event-loop, window and status entry points are no-ops; a headless
#     `run()` returns immediately instead of failing.
#   * Selection entry points return empty selections — nothing is displayed,
#     so nothing can be selected.
# Session rules mirror the rest of the API: most calls require
# `API.initialize()` first (`_model_locked` inside `STATE_LOCK`).

function _fltk_check(action,caller::AbstractString)
    lock(STATE_LOCK) do
        _model_locked()
        return nothing
    end
    return nothing
end

function _fltk_wait(time)
    caller="API.fltk.wait"
    t=try Float64(time) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: time must be a number"))
    end
    isfinite(t) || throw(ArgumentError("$caller: time must be finite"))
    return _fltk_check(nothing,caller)
end

function _fltk_split_current_window(how,ratio)
    caller="API.fltk.split_current_window"
    how isa AbstractString || throw(ArgumentError("$caller: how must be a string"))
    String(how) in ("v","h") || throw(ArgumentError(
        "$caller: how must be \"v\" or \"h\" (got '$(String(how))')"))
    r=try Float64(ratio) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: ratio must be a number"))
    end
    0<r<1 || throw(ArgumentError("$caller: ratio must be inside (0,1)"))
    return _fltk_check(nothing,caller)
end

function _fltk_set_current_window(index)
    caller="API.fltk.set_current_window"
    (index isa Integer && !(index isa Bool)) || throw(ArgumentError(
        "$caller: window index must be an integer"))
    index>=0 || throw(ArgumentError("$caller: window index must be non-negative"))
    return _fltk_check(nothing,caller)
end

function _fltk_set_status_message(message,graphics)
    caller="API.fltk.set_status_message"
    message isa AbstractString || throw(ArgumentError(
        "$caller: message must be a string"))
    graphics isa Bool || throw(ArgumentError(
        "$caller: graphics must be a Bool"))
    return _fltk_check(nothing,caller)
end

function _fltk_show_context_window(dim,tag)
    caller="API.fltk.show_context_window"
    (dim isa Integer && !(dim isa Bool) && 0<=dim<=3) || throw(ArgumentError(
        "$caller: dimension must be an integer in 0:3"))
    (tag isa Integer && !(tag isa Bool) && tag>0) || throw(ArgumentError(
        "$caller: tag must be a positive integer"))
    return _fltk_check(nothing,caller)
end

function _fltk_tree_item(name,caller)
    name isa AbstractString || throw(ArgumentError(
        "$caller: name must be a string"))
    occursin('\0',name) && throw(ArgumentError(
        "$caller: name contains a NUL byte"))
    return _fltk_check(nothing,caller)
end

module fltk

using ..API: _fltk_check,_fltk_wait,_fltk_split_current_window,
             _fltk_set_current_window,_fltk_set_status_message,
             _fltk_show_context_window,_fltk_tree_item

"Initialize the FLTK user interface — no-op: no FLTK backend exists."
initialize()=_fltk_check(nothing,"API.fltk.initialize")
"Finalize the FLTK user interface — no-op: no FLTK backend exists."
finalize()=_fltk_check(nothing,"API.fltk.finalize")
"Wait `time` seconds for events — returns immediately in the headless build."
wait(time=-1.0)=_fltk_wait(time)
"Update the user interface — no-op: no FLTK backend exists."
update()=_fltk_check(nothing,"API.fltk.update")
"Signal an event `action` — no-op: no FLTK backend exists."
awake(action="")=(action isa AbstractString || throw(ArgumentError(
    "API.fltk.awake: action must be a string"));
    _fltk_check(nothing,"API.fltk.awake"))
"Lock the interface — no-op: no FLTK backend exists."
lock()=_fltk_check(nothing,"API.fltk.lock")
"Unlock the interface — no-op: no FLTK backend exists."
unlock()=_fltk_check(nothing,"API.fltk.unlock")
"Run the event loop — returns immediately: no FLTK backend exists."
run(optionFileName="")=(optionFileName isa AbstractString || throw(ArgumentError(
    "API.fltk.run: optionFileName must be a string"));
    _fltk_check(nothing,"API.fltk.run"))
"""
Return `1` when an FLTK backend is available. Always `0`: Tessella ships no
FLTK backend, which is honest both for the missing implementation and for
headless runs without a display.
"""
is_available()=Int32(0)
"Select entities interactively — returns an empty selection: nothing is displayed."
select_entities(dim=-1)=( _fltk_check(nothing,"API.fltk.select_entities");
                          Tuple{Int,Int32}[] )
"Select elements interactively — returns an empty selection: nothing is displayed."
select_elements()=( _fltk_check(nothing,"API.fltk.select_elements");
                    Int[] )
"Select views interactively — returns an empty selection: nothing is displayed."
select_views()=( _fltk_check(nothing,"API.fltk.select_views");
                 Int[] )
"Split the current window — no-op: no FLTK backend exists."
split_current_window(how="v",ratio=0.5)=
    _fltk_split_current_window(how,ratio)
"Set the current window — no-op: no FLTK backend exists."
set_current_window(windowIndex=0)=_fltk_set_current_window(windowIndex)
"Set a status message — no-op: no FLTK backend exists."
set_status_message(message,graphics=false)=
    _fltk_set_status_message(message,graphics)
"Show the context window for entity `(dim, tag)` — no-op: no FLTK backend exists."
show_context_window(dim,tag)=_fltk_show_context_window(dim,tag)
"Open the tree item `name` — no-op: no FLTK backend exists."
open_tree_item(name)=_fltk_tree_item(name,"API.fltk.open_tree_item")
"Close the tree item `name` — no-op: no FLTK backend exists."
close_tree_item(name)=_fltk_tree_item(name,"API.fltk.close_tree_item")

end # module fltk
