# Stage-4 RESOLUTION (shipped): mesh_box — Kuhn/Freudenthal structured subdivision.
#
# The three failures in this folder pointed at the fix: emit EXPLICIT connectivity
# instead of Delaunay-ing degenerate points. mesh_box subdivides a box grid (spacing
# hmax/√3) into 6 path-tetrahedra per cube, all sharing the main diagonal. Result:
#   • maxedge ≤ hmax by construction (cube main diagonal = (hmax/√3)·√3 = hmax)
#   • all-positive tets (validate.ok), watertight (boundary Euler χ = 2)
#   • exact volume, sliver-free (min dihedral 45°/≥42°, radius-edge < 0.9)
#   • terminating (finite grid), deterministic
# Covers axis-aligned box regions (the enclosure case/air cavities).
#
# Run:  julia --project=. validation/stage4_size_refinement/works_mesh_box_kuhn.jl
using Pkg; Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella, Tessella.MeshTypes, Tessella.Geometry
using Printf

maxedge(m)=maximum(begin vs=(m.tets[1,t],m.tets[2,t],m.tets[3,t],m.tets[4,t]); e=0.0
    for i in 1:4,j in i+1:4; p=node(m,vs[i]);q=node(m,vs[j]); e=max(e,hypot(p[1]-q[1],p[2]-q[2],p[3]-q[3])); end; e end for t in 1:ntets(m))
meshvol(m)=sum(tet_volume(node(m,m.tets[1,t]),node(m,m.tets[2,t]),node(m,m.tets[3,t]),node(m,m.tets[4,t])) for t in 1:ntets(m); init=0.0)

for (bx,hmax) in (((0.0,4,0,4,0,4),3.0), ((0.0,4,0,4,0,4),1.0), ((0.0,4,0,4,0,4),0.5), ((-1.0,2,0,5,1,3),0.8))
    m=mesh_box(bx...; hmax=hmax)
    dmin=minimum(tet_dihedral_extrema(node(m,m.tets[1,t]),node(m,m.tets[2,t]),node(m,m.tets[3,t]),node(m,m.tets[4,t]))[1] for t in 1:ntets(m))
    boxvol=(bx[2]-bx[1])*(bx[4]-bx[3])*(bx[6]-bx[5])
    @printf("box=%s hmax=%.1f ntets=%d maxedge=%.4f (≤%.1f) vol=%.6f (exact %.1f) valid=%s χ_bnd=%d mindih=%.2f°\n",
            bx, hmax, ntets(m), maxedge(m), hmax, meshvol(m), boxvol, validate(m).ok, boundary_euler(m), rad2deg(dmin))
end
# Observed: every case maxedge ≤ hmax, exact volume, valid=true, χ_bnd=2 (sphere),
#           min dihedral ≥ 42° — a correct, sliver-free, size-controlled box mesh.
