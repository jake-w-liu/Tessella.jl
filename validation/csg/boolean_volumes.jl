# Native mesh-CSG check (uses the SHIPPED mesh_boolean) — no OpenCASCADE.
#
# Independent oracle: divergence-theorem volume of each Boolean RESULT surface
# equals the analytic value, the result is watertight, and it fills into a valid
# tet mesh via recover_boundary. Covers axis-aligned boxes (exact plane-arrangement
# path, incl. coplanar shared faces) and box×cylinder (exact tri-tri + Mesh2D-CDT).
#
# Run:  julia --project=. validation/csg/boolean_volumes.jl
using Pkg; Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella, Tessella.MeshTypes, Tessella.Geometry
using Printf

divvol(s)=abs(sum(begin a=node(s,s.tris[1,t]);b=node(s,s.tris[2,t]);c=node(s,s.tris[3,t])
    a[1]*(b[2]*c[3]-b[3]*c[2])-a[2]*(b[1]*c[3]-b[3]*c[1])+a[3]*(b[1]*c[2]-b[2]*c[1]) end
    for t in 1:ntris(s); init=0.0)/6)
wt(s)=isempty(first(boundary_edges(s.tris)))

function trio(name, A, B, vu, vi, vd)
    for (op,ex) in ((:union,vu),(:intersection,vi),(:difference,vd))
        R=mesh_boolean(A,B,op); m=recover_boundary(R)
        @printf("%-22s %-13s vol=%.4f(exact %.4f:%s) watertight=%s fill_valid=%s\n",
            name, op, divvol(R), ex, isapprox(divvol(R),ex;rtol=1e-9), wt(R), validate(m).ok)
    end
end

trio("box∩box (x-offset)",  box_surface(0,4,0,4,0,4),  box_surface(2,6,0,4,0,4),   96.0, 32.0, 32.0)
trio("box∩box (xy-offset)", box_surface(0,4,0,4,0,4),  box_surface(2,6,2,6,0,4),  112.0, 16.0, 48.0)
trio("box∩box (gen-pos)",   box_surface(0,10,0,10,0,10), box_surface(5,15,5,15,5,15), 1875.0, 125.0, 875.0)

A=box_surface(0,4,0,4,0,4); cyl=cylinder_surface((2.,2.,-1.),(0.,0.,1.),1.0,6.0; nθ=16)
removed=0.5*16*1.0^2*sin(2pi/16)*4.0; R=mesh_boolean(A,cyl,:difference)
@printf("box − cylinder         difference    vol=%.4f(exact %.4f:%s) watertight=%s fill_valid=%s\n",
    divvol(R), 64.0-removed, isapprox(divvol(R),64.0-removed;rtol=1e-9), wt(R), validate(recover_boundary(R)).ok)
# Observed: every Boolean volume exact, watertight, and recover_boundary-fillable.
