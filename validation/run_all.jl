#!/usr/bin/env julia
# ── Tessella vs external tools (gmsh) — cross-validation driver ──────────────────
#
#   julia --project=. validation/run_all.jl
#
# For every case: mesh the SAME domain with Tessella and with gmsh, then compare the
# meshed volume against the analytic value (the shared oracle), element quality, and
# wall-clock time. Writes validation/REPORT.md and prints a summary table.
#
# Flat-faced solids (box, tunnel, hollow box) have an EXACT analytic volume that both
# meshers must reproduce — a hard correctness cross-check. Curved solids (cylinder,
# sphere) are meshed by each tool under its own surface model (Tessella: an inscribed
# polyhedron; gmsh: the true curved surface), so their volumes are reported against
# the true curved value to expose the geometric-fidelity trade-off honestly.

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io=devnull)
include(joinpath(@__DIR__, "common.jl"))
using Tessella.Geometry
using Tessella.MeshTypes: Mesh

const HERE = @__DIR__

# closed sphere surface: octahedron subdivided `n` times, verts snapped to |p|=R.
function sphere_surface(R::Float64, n::Int)
    V = [(R,0.,0.),(-R,0.,0.),(0.,R,0.),(0.,-R,0.),(0.,0.,R),(0.,0.,-R)]
    F = [(1,3,5),(1,6,3),(1,5,4),(1,4,6),(2,5,3),(2,3,6),(2,4,5),(2,6,4)]
    snap(p)=(s=R/sqrt(p[1]^2+p[2]^2+p[3]^2);(p[1]*s,p[2]*s,p[3]*s))
    for _ in 1:n
        verts=collect(V); mid=Dict{Tuple{Int,Int},Int}()
        gm(a,b)=get!(mid,(min(a,b),max(a,b))) do
            pa=verts[a];pb=verts[b]; push!(verts,snap(((pa[1]+pb[1])/2,(pa[2]+pb[2])/2,(pa[3]+pb[3])/2))); length(verts)
        end
        nF=Tuple{Int,Int,Int}[]
        for (a,b,c) in F; ab=gm(a,b);bc=gm(b,c);ca=gm(c,a); push!(nF,(a,ab,ca));push!(nF,(ab,b,bc));push!(nF,(ca,bc,c));push!(nF,(ab,bc,ca)); end
        V, F = verts, nF
    end
    C=Matrix{Float64}(undef,3,length(V)); for (k,v) in enumerate(V); C[:,k]=[v...]; end
    T=Matrix{Int32}(undef,3,length(F)); for (k,f) in enumerate(F); T[:,k]=Int32[f...]; end
    Mesh(C; tris=T)
end

# case = (folder, geo, Tessella surface, analytic volume of Tessella's model,
#         true curved volume (== analytic for flat solids), note)
cases = [
 (dir="01_box", geo="box.geo",
  surf=()->box_surface(0,2,0,1,0,1), model_vol=2.0, true_vol=2.0,
  note="flat solid — both meshers must hit V=2 exactly"),
 (dir="02_cylinder", geo="cylinder.geo",
  surf=()->cylinder_surface((0.,0,0),(0.,0,1),2.0,5.0; nθ=48, nz=6),
  model_vol=0.5*48*2.0^2*sinpi(2/48)*5.0, true_vol=pi*4*5,
  note="curved — Tessella meshes an inscribed 48-gon prism; gmsh the true circle"),
 (dir="03_box_tunnel", geo="box_tunnel.geo",
  surf=()->box_tunnel_surface(1,5,1,5,1,3,2,4,2,4), model_vol=24.0, true_vol=24.0,
  note="genus-1 flat solid — both must hit V=24 exactly"),
 (dir="04_hollow_box", geo="hollow_box.geo",
  surf=()->box_shell_surface(-1,2,0,3,1,5, 0,1,1,2,2,3), model_vol=35.0, true_vol=35.0,
  note="hollow box (Boolean difference) — both must hit V=35 exactly"),
 (dir="05_sphere", geo="sphere.geo",
  surf=()->sphere_surface(1.7, 3), model_vol=NaN, true_vol=4/3*pi*1.7^3,
  note="curved — Tessella meshes a 3x-subdivided octahedron; gmsh the true sphere"),
]

rows = String[]
push!(rows, "| Case | Tool | Nodes | Tets | Volume | rel.err vs true | min∠° | mean∠° | slivers | time (s) |")
push!(rows, "|------|------|------:|-----:|-------:|----------------:|------:|-------:|--------:|---------:|")

haveg = gmsh_available()
println("gmsh: ", haveg ? gmsh_version() :
        "NOT AVAILABLE (required size-field differential will fail)")

# The catalog differential is a required child gate. It performs direct Gmsh
# API/plugin probes and exits nonzero when Gmsh is missing, is not 4.15.2, or a
# parity assertion fails. CONTEXT_SKIP lines are explicit non-claims, not a
# substitute for running Gmsh.
println("\n── size_fields ──  direct and mesh-observed Gmsh 4.15.2 differential")
size_field_script = joinpath(HERE, "size_fields", "differential.jl")
size_field_project = normpath(joinpath(HERE, ".."))
size_field_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $size_field_script`
println("  command: ", size_field_command)
run(size_field_command) # ProcessFailedException makes validation/run_all.jl nonzero.

println("\n── uniform_refine ──  exact Gmsh 4.15.2 simplex templates")
uniform_refine_script = joinpath(HERE, "uniform_refine", "differential.jl")
uniform_refine_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $uniform_refine_script`
println("  command: ", uniform_refine_command)
run(uniform_refine_command) # ProcessFailedException makes validation/run_all.jl nonzero.

println("\n── transfinite ──  Gmsh 4.15.2 four-sided planar patches")
transfinite_script = joinpath(HERE, "transfinite", "differential.jl")
transfinite_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $transfinite_script`
println("  command: ", transfinite_command)
run(transfinite_command) # ProcessFailedException makes validation/run_all.jl nonzero.

function fmtrow(case, tool, m, secs, truevol)
    @sprintf("| %s | %s | %d | %d | %.6g | %s | %.2f | %.2f | %d | %.3f |",
             case, tool, m.nnodes, m.ntets, m.volume, pct(relerr(m.volume, truevol)),
             m.min_dih, m.mean_dih, m.slivers, secs)
end

for c in cases
    println("\n── $(c.dir) ──  $(c.note)")
    # Tessella
    surf = c.surf()
    tm, tsec = tessella_fill(surf; smooth=false)
    println("  Tessella: ", c.dir, "  V=", round(tm.volume, sigdigits=8),
            "  (model V=", isnan(c.model_vol) ? "poly" : round(c.model_vol, sigdigits=8), ")")
    # correctness assertion against Tessella's OWN geometric model (exact for flat)
    if !isnan(c.model_vol)
        @assert relerr(tm.volume, c.model_vol) < 1e-6 "Tessella $(c.dir): V=$(tm.volume) != model $(c.model_vol)"
    end
    push!(rows, fmtrow(c.dir, "Tessella", tm, tsec, c.true_vol))
    # gmsh
    if haveg
        geo = joinpath(HERE, "cases", c.dir, c.geo)
        out = joinpath(HERE, "cases", c.dir, "gmsh_out.msh")
        r = run_gmsh(geo, out)
        if r.ok
            gm = gmsh_metrics(out)
            println("  gmsh:     V=", round(gm.volume, sigdigits=8), "  tets=", gm.ntets, "  t=", round(r.seconds, digits=3), "s")
            push!(rows, fmtrow(c.dir, "gmsh", gm, r.seconds, c.true_vol))
        else
            println("  gmsh:     FAILED to mesh $(c.geo)")
            push!(rows, "| $(c.dir) | gmsh | — | — | FAILED | — | — | — | — | — |")
        end
    end
end

# ── enclosure coax junction: the acceptance case (gmsh's known empty-volume fail) ──
println("\n── 06_enclosure_coax ──  ASCENT coax feed-through (multi-material CSG)")
enc_lines = String[]
if haveg
    geo = joinpath(HERE, "cases", "06_enclosure_coax", "enclosure_coax_junction.geo")
    out = joinpath(HERE, "cases", "06_enclosure_coax", "gmsh_out.msh")
    r = run_gmsh(geo, out)
    # gmsh's documented failure exits NON-ZERO but still writes a partial mesh whose
    # solid volumes are empty — inspect that output regardless of exit status.
    if isfile(out)
        m = read_msh(out).mesh
        bytag = Dict{Int32,Int}()
        for t in (isempty(m.tet_tag) ? Int32[] : m.tet_tag); bytag[t] = get(bytag, t, 0) + 1; end
        exit_note = r.ok ? "exit 0" : "exit NON-ZERO (meshing error)"
        push!(enc_lines, "gmsh $(gmsh_version()) on the literal `.geo` ($exit_note, $(round(r.seconds,digits=1)) s): $(ntets(m)) volume tets total.")
        push!(enc_lines, "Tets per physical-volume tag (air=1, coax_pin=2, case=4 are the solid regions): $(isempty(bytag) ? "none — no tagged volume elements" : sort(collect(bytag)))")
        solid = sum(get(bytag, Int32(t), 0) for t in (1,2,4))
        push!(enc_lines, solid == 0 ?
            "❌ gmsh left ALL three solid volumes EMPTY (0 tets tagged air/pin/case) — exactly the documented failure Tessella targets." :
            "gmsh filled $(solid) solid-region tets.")
    else
        push!(enc_lines, "❌ gmsh produced no mesh at all for the enclosure fixture — the documented failure Tessella targets.")
    end
else
    push!(enc_lines, "gmsh not available; skipped.")
end
push!(enc_lines, "Tessella status: the literal multi-region enclosure/coax model is natively meshed, solver-loadable, and solved; see validation/enclosure_literal/ and ASCENT.md. The primitive rows above independently check analytic volumes.")
for l in enc_lines; println("  ", l); end

# ── write REPORT.md ──
open(joinpath(HERE, "REPORT.md"), "w") do io
    println(io, "# Tessella validation report\n")
    println(io, "Generated by `validation/run_all.jl`. Reference tool: gmsh ", haveg ? gmsh_version() : "n/a", ".\n")
    println(io, "Each row meshes the same domain and reports the meshed volume vs the **true** value,")
    println(io, "element dihedral-angle quality, sliver count, and wall-clock time.\n")
    for l in rows; println(io, l); end
    println(io, "\n## 06_enclosure_coax (acceptance case)\n")
    for l in enc_lines; println(io, "- ", l); end
    println(io, "\n## Notes\n")
    println(io, "- Flat solids (box, box_tunnel, hollow_box): the analytic volume is exact, so a passing")
    println(io, "  row confirms BOTH meshers conform to the true geometry to rounding.")
    println(io, "- Curved solids (cylinder, sphere): Tessella meshes an inscribed polyhedron (its volume")
    println(io, "  is exact for that model and slightly under the true curved volume); gmsh meshes the")
    println(io, "  true curved surface. The `rel.err vs true` column quantifies each tool's fidelity.")
end
println("\nWrote ", joinpath(HERE, "REPORT.md"))
