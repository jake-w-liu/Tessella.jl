# Boundary-recovery conformity check (uses the SHIPPED recover_boundary).
#
# Independent geometric conformity oracle: the recovered tet mesh's BOUNDARY area
# equals the input surface area (the boundary IS the input surface, triangulation-
# independently), the domain volume is exact, the mesh is valid and closed-manifold.
# Covers convex, non-convex genus-1 (through-tunnel), and hollow-shell PLCs.
#
# Run:  julia --project=. validation/boundary_recovery/conformity_check.jl
using Pkg; Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella, Tessella.MeshTypes, Tessella.Geometry
using Printf

meshvol(m)=sum(tet_volume(node(m,m.tets[1,t]),node(m,m.tets[2,t]),node(m,m.tets[3,t]),node(m,m.tets[4,t])) for t in 1:ntets(m); init=0.0)
surfarea(s)=sum(triangle_area(node(s,s.tris[1,t]),node(s,s.tris[2,t]),node(s,s.tris[3,t])) for t in 1:ntris(s); init=0.0)
bndarea(m)=(bf=first(boundary_faces(m.tets)); sum(triangle_area(node(m,f[1]),node(m,f[2]),node(m,f[3])) for f in bf; init=0.0))

function check(name, s, exactvol)
    sa=surfarea(s); t0=time()
    m=recover_boundary(s); dt=time()-t0
    ba=bndarea(m); v=meshvol(m)
    @printf("%-14s ntets=%-5d vol=%.4f(exact %.4f:%s) bnd_area=%.4f==surf %.4f:%s valid=%s closed=%s %.1fs\n",
        name, ntets(m), v, exactvol, isapprox(v,exactvol;rtol=1e-6),
        ba, sa, isapprox(ba,sa;rtol=1e-6), validate(m).ok, is_closed_manifold(m), dt)
end

check("box",         box_surface(0,4,0,4,0,4), 64.0)
check("box_tunnel",  box_tunnel_surface(0,6, 0,6, 0,6, 2,4, 2,4), 216.0-24.0)   # genus-1
check("box_shell",   box_shell_surface(0,6,0,6,0,6, 1,5,1,5,1,5), 216.0-64.0)   # hollow
# Observed: every case boundary-area == surface-area, exact volume, valid, closed — conforming.
