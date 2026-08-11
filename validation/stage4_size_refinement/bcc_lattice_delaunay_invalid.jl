# Stage-4 measured finding (3/3): BCC-lattice + Delaunay is robust only for the
# trivial even-spacing case; general spacing yields INVALID (inverted) tets.
#
# A body-centered-cubic (BCC) point set is the classic high-quality structured tet
# mesh (congruent tets, 60°/90° dihedrals). Placing BCC points in a box at spacing
# a = hmax/√2 (so face diagonals ≤ hmax) and taking their Delaunay works PERFECTLY
# when the box is an even multiple of the cell (hmax=3.0 → 2×2×2, a=2: valid, 60°
# min dihedral, exact volume, max edge bounded) — but for general spacing the mesh
# is INVALID: inverted tets (min dihedral 0°) and even wrong total volume (64→64.44).
#
# Mechanism: a regular lattice is *maximally* cospherical-degenerate. Delaunay of
# cospherical points is ambiguous; the exact-predicate + SoS kernel resolves ties
# deterministically but not always into a *valid* (non-inverted, non-overlapping)
# tetrahedralization. A correct BCC mesher must emit the KNOWN BCC connectivity
# explicitly (a structured-meshing subsystem), not Delaunay-of-the-lattice.
#
# Run:  julia --project=. validation/stage4_size_refinement/bcc_lattice_delaunay_invalid.jl
using Pkg; Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella.MeshTypes, Tessella.Geometry
import Tessella.Mesh3D as M3
using Printf

function bcc_box(x0,x1,y0,y1,z0,z1, hmax)
    a0 = hmax/sqrt(2)
    nx = max(1, ceil(Int,(x1-x0)/a0)); ny = max(1, ceil(Int,(y1-y0)/a0)); nz = max(1, ceil(Int,(z1-z0)/a0))
    ax = (x1-x0)/nx; ay=(y1-y0)/ny; az=(z1-z0)/nz
    xs=Float64[]; ys=Float64[]; zs=Float64[]
    for i in 0:nx, j in 0:ny, k in 0:nz
        push!(xs, x0+i*ax); push!(ys, y0+j*ay); push!(zs, z0+k*az)
    end
    for i in 0:nx-1, j in 0:ny-1, k in 0:nz-1
        push!(xs, x0+(i+0.5)*ax); push!(ys, y0+(j+0.5)*ay); push!(zs, z0+(k+0.5)*az)
    end
    return xs,ys,zs
end
maxedge(m)=begin mx=0.0
    for t in 1:size(m.tets,2); vs=(m.tets[1,t],m.tets[2,t],m.tets[3,t],m.tets[4,t])
        for i in 1:4,j in i+1:4; p=node(m,vs[i]);q=node(m,vs[j]); d=sqrt((p[1]-q[1])^2+(p[2]-q[2])^2+(p[3]-q[3])^2); d>mx&&(mx=d); end
    end; mx; end
meshvol(m)=sum(tet_volume(node(m,m.tets[1,t]),node(m,m.tets[2,t]),node(m,m.tets[3,t]),node(m,m.tets[4,t])) for t in 1:size(m.tets,2); init=0.0)
function mindihedral(m)
    mn=180.0
    for t in 1:size(m.tets,2)
        p=(node(m,m.tets[1,t]),node(m,m.tets[2,t]),node(m,m.tets[3,t]),node(m,m.tets[4,t]))
        for (a,b,c,d) in ((1,2,3,4),(1,3,2,4),(1,4,2,3),(2,3,1,4),(2,4,1,3),(3,4,1,2))
            u=(p[b][1]-p[a][1],p[b][2]-p[a][2],p[b][3]-p[a][3])
            v=(p[c][1]-p[a][1],p[c][2]-p[a][2],p[c][3]-p[a][3])
            w=(p[d][1]-p[a][1],p[d][2]-p[a][2],p[d][3]-p[a][3])
            n1=(u[2]*v[3]-u[3]*v[2],u[3]*v[1]-u[1]*v[3],u[1]*v[2]-u[2]*v[1])
            n2=(u[2]*w[3]-u[3]*w[2],u[3]*w[1]-u[1]*w[3],u[1]*w[2]-u[2]*w[1])
            l1=sqrt(n1[1]^2+n1[2]^2+n1[3]^2); l2=sqrt(n2[1]^2+n2[2]^2+n2[3]^2)
            (l1==0||l2==0) && continue
            cang=clamp((n1[1]*n2[1]+n1[2]*n2[2]+n1[3]*n2[3])/(l1*l2),-1,1)
            mn=min(mn,180.0-acosd(cang))
        end
    end
    mn
end

for hmax in (3.0, 2.0, 1.0)
    xs,ys,zs = bcc_box(0,4,0,4,0,4, hmax)
    T=M3.delaunay3d(xs,ys,zs; perturb=false)   # exact coords: flat faces stay on the boundary
    m=M3.to_mesh3(T); d=validate(m)
    @printf("hmax=%.1f npts=%d ntets=%d maxedge=%.4f (bound %.1f) vol=%.4f valid=%s mindih=%.2f°\n",
            hmax, length(xs), size(m.tets,2), maxedge(m), hmax, meshvol(m), d.ok, mindihedral(m))
end
# Observed: hmax=3.0 (2×2×2, a=2) → valid, 60° dihedral, exact vol, bounded edge.
#           hmax=2.0, 1.0 → valid=false / wrong volume: Delaunay of the maximally
#           cospherical lattice yields inverted tets. Needs explicit BCC connectivity.
