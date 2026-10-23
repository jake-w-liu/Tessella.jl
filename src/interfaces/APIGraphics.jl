# `graphics` module: the graphics redraw surface, matching the gmsh 4.15.2
# `gmshGraphics*` API. With no rendering backend, `draw()` is a deterministic
# no-op — upstream `gmsh.graphics.draw()` also returns without error when no
# graphics context exists.

module graphics

using ..API: _fltk_check

"Redraw all views — no-op: Tessella ships no rendering backend."
draw()=_fltk_check(nothing,"API.graphics.draw")

end # module graphics
