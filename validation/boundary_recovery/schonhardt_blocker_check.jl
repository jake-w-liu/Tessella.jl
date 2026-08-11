# Schönhardt blocker check (uses the SHIPPED recover_boundary).
#
# The Schönhardt polyhedron cannot be tetrahedralized using only its 6 vertices
# (needs a Steiner point). recover_boundary must therefore raise its EXPLICIT
# blocker rather than return a silently non-conforming mesh. The convex
# triangulation of the SAME 6 points is tetrahedralizable and meshes cleanly —
# so the blocker is exactly discriminating, not blanket.
#
# Run:  julia --project=. validation/boundary_recovery/schonhardt_blocker_check.jl
using Pkg; Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella, Tessella.MeshTypes, Tessella.Geometry

# twisted triangular prism; diag=:aibn = the REFLEX diagonal ⇒ Schönhardt (blocked);
# diag=:anbi = the convex diagonal ⇒ tetrahedralizable (meshes).
function schon(θ, diag)
    C=Matrix{Float64}(undef,3,6)
    for k in 0:2; C[:,k+1]=[cos(2pi*k/3),sin(2pi*k/3),0.0]; end
    for k in 0:2; C[:,k+4]=[cos(2pi*k/3+θ),sin(2pi*k/3+θ),1.0]; end
    tris=NTuple{3,Int}[(1,3,2),(4,5,6)]
    for i in 1:3
        ai=i; an=i%3+1; bi=i+3; bn=(i%3)+1+3
        if diag==:anbi; push!(tris,(ai,an,bi)); push!(tris,(an,bn,bi))
        else;           push!(tris,(ai,an,bn)); push!(tris,(ai,bn,bi)); end
    end
    T=Matrix{Int32}(undef,3,length(tris)); for (t,f) in enumerate(tris); T[:,t]=Int32[f...]; end
    Mesh(C; tris=T)
end

for (θ,diag) in ((pi/6,:aibn),(pi/6,:anbi),(pi/4,:aibn),(pi/4,:anbi))
    s=schon(θ,diag)
    print("θ=$(round(Int,rad2deg(θ)))° diag=$diag → ")
    try
        m=recover_boundary(s; max_seeds=16)
        println("MESHED ntets=$(ntets(m)) (tetrahedralizable)")
    catch e
        println("BLOCKED ✓ (explicit): ", sprint(showerror,e)[1:min(end,60)])
    end
end
# Observed: reflex diagonal (:aibn) ⇒ BLOCKED; convex diagonal (:anbi) ⇒ MESHED.
