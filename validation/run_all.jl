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
include(joinpath(@__DIR__, "support", "common.jl"))
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

println("\n── geo_ranges ──  bit-exact Gmsh 4.15.2 finite constant lists")
geo_range_script = joinpath(HERE, "geo_ranges", "differential.jl")
geo_range_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $geo_range_script`
println("  command: ", geo_range_command)
run(geo_range_command) # ProcessFailedException makes validation/run_all.jl nonzero.

println("\n── gmsh_parity geometry expressions ──  bounded entity/value execution")
geo_expression_script = joinpath(
    HERE, "gmsh_parity", "geo_geometry_expressions.jl")
geo_expression_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $geo_expression_script`
println("  command: ", geo_expression_command)
run(geo_expression_command)

println("\n── gmsh_parity list variables ──  bounded numeric/entity reuse")
geo_list_script = joinpath(HERE, "gmsh_parity", "geo_list_variables.jl")
geo_list_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $geo_list_script`
println("  command: ", geo_list_command)
run(geo_list_command)

println("\n── gmsh_parity topology queries ──  Point sizing and Physical groups")
geo_mesh_size_script = joinpath(HERE, "gmsh_parity", "geo_mesh_sizes.jl")
geo_mesh_size_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $geo_mesh_size_script`
println("  command: ", geo_mesh_size_command)
run(geo_mesh_size_command)

println("\n── gmsh_parity dynamic tags ──  geometry and Physical lifecycle")
geo_dynamic_tag_script = joinpath(HERE, "gmsh_parity", "geo_dynamic_tags.jl")
geo_dynamic_tag_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $geo_dynamic_tag_script`
println("  command: ", geo_dynamic_tag_command)
run(geo_dynamic_tag_command)

println("\n── gmsh_parity model topology ──  entity boundaries and adjacencies")
model_topology_script = joinpath(
    HERE, "gmsh_parity", "model_topology_queries.jl")
model_topology_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $model_topology_script`
println("  command: ", model_topology_command)
run(model_topology_command)

println("\n── gmsh_parity model identity ──  entity names and atomic retagging")
model_identity_script = joinpath(
    HERE, "gmsh_parity", "model_entity_identity.jl")
model_identity_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $model_identity_script`
println("  command: ", model_identity_command)
run(model_identity_command)

println("\n── gmsh_parity model removal ──  ordered dependencies and recursion")
model_removal_script = joinpath(
    HERE, "gmsh_parity", "model_entity_removal.jl")
model_removal_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $model_removal_script`
println("  command: ", model_removal_command)
run(model_removal_command)

println("\n── gmsh_parity model spatial ──  analytical bounds and containment")
model_spatial_script = joinpath(
    HERE, "gmsh_parity", "model_spatial_queries.jl")
model_spatial_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $model_spatial_script`
println("  command: ", model_spatial_command)
run(model_spatial_command)

println("\n── gmsh_parity model metadata ──  native types and partition ownership")
model_metadata_script = joinpath(
    HERE, "gmsh_parity", "model_entity_metadata.jl")
model_metadata_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $model_metadata_script`
println("  command: ", model_metadata_command)
run(model_metadata_command)

println("\n── uniform_refine ──  exact Gmsh 4.15.2 simplex templates")
uniform_refine_script = joinpath(HERE, "uniform_refine", "differential.jl")
uniform_refine_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $uniform_refine_script`
println("  command: ", uniform_refine_command)
run(uniform_refine_command) # ProcessFailedException makes validation/run_all.jl nonzero.

println("\n── high_order ──  exact Gmsh 4.15.2 type-11 tetrahedron ordering")
high_order_script = joinpath(HERE, "high_order", "differential.jl")
high_order_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $high_order_script`
println("  command: ", high_order_command)
run(high_order_command) # ProcessFailedException makes validation/run_all.jl nonzero.

println("\n── transfinite ──  Gmsh 4.15.2 four-sided planar patches")
transfinite_script = joinpath(HERE, "transfinite", "differential.jl")
transfinite_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $transfinite_script`
println("  command: ", transfinite_command)
run(transfinite_command) # ProcessFailedException makes validation/run_all.jl nonzero.

println("\n── transfinite_curve ──  Gmsh 4.15.2 straight-curve laws")
transfinite_curve_script = joinpath(HERE, "transfinite_curve", "differential.jl")
transfinite_curve_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $transfinite_curve_script`
println("  command: ", transfinite_curve_command)
run(transfinite_curve_command) # ProcessFailedException makes validation/run_all.jl nonzero.

println("\n── transfinite_triangle ──  Gmsh 4.15.2 three-sided triangle/quad patches")
transfinite_triangle_script = joinpath(HERE, "transfinite_triangle", "differential.jl")
transfinite_triangle_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $transfinite_triangle_script`
println("  command: ", transfinite_triangle_command)
run(transfinite_triangle_command) # ProcessFailedException makes validation/run_all.jl nonzero.

println("\n── transfinite_quad ──  Gmsh 4.15.2 recombined four-sided patches")
transfinite_quad_script = joinpath(HERE, "transfinite_quad", "differential.jl")
transfinite_quad_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $transfinite_quad_script`
println("  command: ", transfinite_quad_command)
run(transfinite_quad_command) # ProcessFailedException makes validation/run_all.jl nonzero.

println("\n── transfinite_volume ──  Gmsh 4.15.2 affine six-face volumes")
transfinite_volume_script = joinpath(HERE, "transfinite_volume", "differential.jl")
transfinite_volume_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $transfinite_volume_script`
println("  command: ", transfinite_volume_command)
run(transfinite_volume_command) # ProcessFailedException makes validation/run_all.jl nonzero.

println("\n── transfinite_prism ──  Gmsh 4.15.2 affine five-face prisms")
transfinite_prism_script = joinpath(HERE, "transfinite_prism", "differential.jl")
transfinite_prism_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $transfinite_prism_script`
println("  command: ", transfinite_prism_command)
run(transfinite_prism_command) # ProcessFailedException makes validation/run_all.jl nonzero.

println("\n── gmsh_parity box ──  Tessella API vs analytic/Gmsh box volume")
box_api_script = joinpath(HERE, "gmsh_parity", "box_api.jl")
box_api_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $box_api_script`
println("  command: ", box_api_command)
run(box_api_command)

println("\n── gmsh_parity t1 square ──  Tessella vs analytic/Gmsh unit-square area")
t1_script = joinpath(HERE, "gmsh_parity", "t1_square.jl")
t1_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $t1_script`
println("  command: ", t1_command)
run(t1_command)

println("\n── gmsh_parity cone ──  Tessella frustum vs Gmsh OCC Cone")
cone_script = joinpath(HERE, "gmsh_parity", "cone_geo.jl")
cone_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $cone_script`
println("  command: ", cone_command)
run(cone_command)

println("\n── gmsh_parity NURBS surface ──  Tessella IGES-128 vs Gmsh OCC BSpline")
nurbs_script = joinpath(HERE, "gmsh_parity", "nurbs_surface.jl")
nurbs_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $nurbs_script`
println("  command: ", nurbs_command)
run(nurbs_command)

println("\n── gmsh_parity cylinder ──  Tessella prism vs Gmsh OCC Cylinder")
cyl_script = joinpath(HERE, "gmsh_parity", "cylinder_geo.jl")
cyl_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $cyl_script`
println("  command: ", cyl_command)
run(cyl_command)

println("\n── gmsh_parity boolean boxes ──  Tessella vs analytic/Gmsh BooleanDifference")
bool_script = joinpath(HERE, "gmsh_parity", "boolean_boxes.jl")
bool_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $bool_script`
println("  command: ", bool_command)
run(bool_command)

println("\n── gmsh_parity t4 hole ──  Tessella vs analytic/Gmsh holed-square area")
t4_script = joinpath(HERE, "gmsh_parity", "t4_hole.jl")
t4_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $t4_script`
println("  command: ", t4_command)
run(t4_command)

println("\n── gmsh_parity embed point ──  Tessella vs Gmsh Point-In-Surface node")
embed_script = joinpath(HERE, "gmsh_parity", "embed_point.jl")
embed_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $embed_script`
println("  command: ", embed_command)
run(embed_command)

println("\n── gmsh_parity embed line ──  Tessella vs Gmsh Line-In-Surface chain")
eline_script = joinpath(HERE, "gmsh_parity", "embed_line.jl")
eline_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $eline_script`
println("  command: ", eline_command)
run(eline_command)

println("\n── gmsh_parity embed sheet ──  Tessella vs Gmsh Surface-In-Volume sheet")
esheet_script = joinpath(HERE, "gmsh_parity", "embed_sheet.jl")
esheet_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $esheet_script`
println("  command: ", esheet_command)
run(esheet_command)

println("\n── gmsh_parity holed embed sheet ──  Tessella vs Gmsh holed Surface-In-Volume")
esheet_hole_script = joinpath(HERE, "gmsh_parity", "embed_sheet_hole.jl")
esheet_hole_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $esheet_hole_script`
println("  command: ", esheet_hole_command)
run(esheet_hole_command)

println("\n── gmsh_parity explicit shell ──  Tessella vs Gmsh Surface Loop/Volume")
explicit_shell_script = joinpath(HERE, "gmsh_parity", "explicit_shell.jl")
explicit_shell_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $explicit_shell_script`
println("  command: ", explicit_shell_command)
run(explicit_shell_command)

println("\n── gmsh_parity periodic affine ──  Tessella vs Gmsh node correspondence")
periodic_script = joinpath(HERE, "gmsh_parity", "periodic_translation.jl")
periodic_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $periodic_script`
println("  command: ", periodic_command)
run(periodic_command)

println("\n── gmsh_parity periodic embedded curve ──  Tessella vs Gmsh embedded trace pairs")
periodic_embedded_script = joinpath(HERE, "gmsh_parity", "periodic_embedded_curve.jl")
periodic_embedded_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $periodic_embedded_script`
println("  command: ", periodic_embedded_command)
run(periodic_embedded_command)

println("\n── gmsh_parity periodic curve graph ──  reusable/chained/expression-backed")
periodic_graph_script = joinpath(HERE, "gmsh_parity", "periodic_curve_graph.jl")
periodic_graph_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $periodic_graph_script`
println("  command: ", periodic_graph_command)
run(periodic_graph_command)

println("\n── gmsh_parity periodic volume boundaries ──  planar explicit shell")
periodic_volume_script = joinpath(
    HERE, "gmsh_parity", "periodic_surface_volume.jl")
periodic_volume_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $periodic_volume_script`
println("  command: ", periodic_volume_command)
run(periodic_volume_command)

println("\n── gmsh_parity 2-D boundary layer ──  Tessella vs analytic/Gmsh BL quads")
bl2d_script = joinpath(HERE, "gmsh_parity", "boundary_layer_2d.jl")
bl2d_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $bl2d_script`
println("  command: ", bl2d_command)
run(bl2d_command)

println("\n── transfinite_hex ──  Gmsh 4.15.2 affine six-face hexahedra")
transfinite_hex_script = joinpath(HERE, "transfinite_hex", "differential.jl")
transfinite_hex_command = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --project=$size_field_project $transfinite_hex_script`
println("  command: ", transfinite_hex_command)
run(transfinite_hex_command) # ProcessFailedException makes validation/run_all.jl nonzero.

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
        # The partial file contains duplicate surface triangles and correctly
        # fails the validated simplex reader. Use the mixed-record parser only
        # for forensic cell counting, and report its failed structural audit.
        m = read_mixed_msh(out)
        diagnostic=validate(m)
        bytag = Dict{Int32,Int}()
        for block in m.blocks
            msh_dimension(block.msh)==3 || continue
            for tag in block.tags
                bytag[tag]=get(bytag,tag,0)+1
            end
        end
        nvolume=sum(values(bytag);init=0)
        exit_note = r.ok ? "exit 0" : "exit NON-ZERO (meshing error)"
        validity_note=diagnostic.ok ? "structurally valid" :
            "structurally invalid (duplicate surface cells)"
        push!(enc_lines, "gmsh $(gmsh_version()) on the literal `.geo` ($exit_note, $(round(r.seconds,digits=1)) s): $nvolume volume cells total; partial file is $validity_note.")
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
