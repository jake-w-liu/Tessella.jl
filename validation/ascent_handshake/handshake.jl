# ASCENT side of the mesh handshake: parse the Tessella-written .msh exactly as
# ASCENT.load_mesh does (GmshDiscreteModel + top-dimensional physical-group check).
# Run in the ASCENT env (GridapGmsh is an ASCENT dependency, not a Tessella one):
#   julia --project=<2026_066/ASCENT> validation/ascent_handshake/handshake.jl
using GridapGmsh: GmshDiscreteModel
using Gridap.Geometry: get_face_labeling, num_tags, get_tag_name, num_cells

path = joinpath(@__DIR__, "ascent_coax.msh")
isfile(path) || error("run generate.jl (Tessella env) first to produce $path")

# mirrors ASCENT/src/core/mesh.jl:_gmsh_discrete_model (v4.x dense-tag → renumber=false)
model = try
    GmshDiscreteModel(path; renumber=false)
catch
    GmshDiscreteModel(path)
end

# mirrors ASCENT/src/core/mesh.jl:load_mesh — keep only tags on top-dimensional entities
labels = get_face_labeling(model)
nT = num_tags(labels)
top = Set(Int.(labels.d_to_dface_to_entity[end]))
voltags = String[]
for i in 1:nT
    any(Int(e) in top for e in labels.tag_to_entities[i]) || continue
    push!(voltags, String(get_tag_name(labels, i)))
end

println("RESULT num_cells=", num_cells(model), " num_tags=", nT, " volume_groups=", voltags)
if isempty(voltags)
    error("HANDSHAKE_FAIL: ASCENT.load_mesh would reject — no top-dimensional physical groups")
end
println("HANDSHAKE_OK — ASCENT.load_mesh would return a valid MeshData with tags $voltags")
