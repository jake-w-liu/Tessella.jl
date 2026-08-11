# Stage-4 RESOLUTION (shipped): mesh_sized_conforming — interior size control for
# curved domains, via an inset interior Steiner lattice + the exact conformity gate.
#
# For a THICK curved domain (sphere) it returns a conforming mesh with interior edges
# ≤ hmax (genuine size reduction vs the no-interior recover_boundary baseline). For a
# THIN cylinder the inset removes the interior lattice, so it degrades SAFELY to
# conforming-only (never a silent invalid mesh; the gate blocks otherwise).
#
# Run:  julia --project=. validation/stage4_size_refinement/works_sized_conforming.jl
using Pkg; Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella, Tessella.MeshTypes, Tessella.Geometry
using Printf

function icosphere(R,k)
    t=(1+sqrt(5))/2
    V=[(-1.,t,0.),(1.,t,0.),(-1.,-t,0.),(1.,-t,0.),(0.,-1.,t),(0.,1.,t),(0.,-1.,-t),(0.,1.,-t),(t,0.,-1.),(t,0.,1.),(-t,0.,-1.),(-t,0.,1.)]
    F=[(1,12,6),(1,6,2),(1,2,8),(1,8,11),(1,11,12),(2,6,10),(6,12,5),(12,11,3),(11,8,7),(8,2,9),(4,10,5),(4,5,3),(4,3,7),(4,7,9),(4,9,10),(5,10,6),(3,5,12),(7,3,11),(9,7,8),(10,9,2)]
    coords=[v for v in V]; idx=Dict{NTuple{3,Float64},Int}(); for (i,v) in enumerate(coords); idx[v]=i; end; faces=F
    for _ in 1:k
        vid(p)=get!(idx,p) do; push!(coords,p); length(coords); end; nf=NTuple{3,Int}[]
        for (a,b,c) in faces; pa=coords[a];pb=coords[b];pc=coords[c]; md(x,y)=((x[1]+y[1])/2,(x[2]+y[2])/2,(x[3]+y[3])/2)
            ab=vid(md(pa,pb));bc=vid(md(pb,pc));ca=vid(md(pc,pa)); push!(nf,(a,ab,ca));push!(nf,(ab,b,bc));push!(nf,(ca,bc,c));push!(nf,(ab,bc,ca)); end
        faces=nf; end
    pr=[(v[1]/hypot(v...)*R,v[2]/hypot(v...)*R,v[3]/hypot(v...)*R) for v in coords]
    C=Matrix{Float64}(undef,3,length(pr)); for (i,p) in enumerate(pr); C[:,i]=[p...]; end
    Tm=Matrix{Int32}(undef,3,length(faces)); for (t,f) in enumerate(faces); Tm[:,t]=Int32[f...]; end
    Mesh(C; tris=Tm)
end
sfa(s)=sum(triangle_area(node(s,s.tris[1,t]),node(s,s.tris[2,t]),node(s,s.tris[3,t])) for t in 1:ntris(s); init=0.0)
bfa(m)=(bf=first(boundary_faces(m.tets)); sum(triangle_area(node(m,f[1]),node(m,f[2]),node(m,f[3])) for f in bf; init=0.0))
intmax(m)=begin isb=falses(nnodes(m)); for f in first(boundary_faces(m.tets)); isb[f[1]]=isb[f[2]]=isb[f[3]]=true; end; mx=0.0
    for t in 1:ntets(m); vs=(m.tets[1,t],m.tets[2,t],m.tets[3,t],m.tets[4,t])
        for i in 1:4,j in i+1:4; (isb[vs[i]]||isb[vs[j]])&&continue; p=node(m,vs[i]);q=node(m,vs[j]); mx=max(mx,hypot(p[1]-q[1],p[2]-q[2],p[3]-q[3])); end; end; mx; end

s=icosphere(3.0,2); base=recover_boundary(s); m=mesh_sized_conforming(s; hmax=1.5)
@printf("sphere R=3 hmax=1.5: baseline_nt=%d → sized_nt=%d int_maxedge=%.3f(≤1.5) conf=%s valid=%s closed=%s\n",
    ntets(base), ntets(m), intmax(m), isapprox(bfa(m),sfa(s);rtol=1e-9), validate(m).ok, is_closed_manifold(m))
cyl=cylinder_surface((0.,0,0),(0.,0,1),2.0,4.0; nθ=16, nz=4); mc=mesh_sized_conforming(cyl; hmax=1.5)
@printf("thin cylinder hmax=1.5: nt=%d (degrades safely) conf=%s valid=%s closed=%s\n",
    ntets(mc), isapprox(bfa(mc),sfa(cyl);rtol=1e-9), validate(mc).ok, is_closed_manifold(mc))
# Observed: sphere genuinely refined (interior ≤ hmax, conforming); thin cylinder conforming+valid.
