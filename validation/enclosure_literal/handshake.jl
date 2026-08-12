# Confirm the natively-reconstructed LITERAL ENC-COAX mesh loads into ASCENT's parser
# with its three physical volumes. Run in the ASCENT env (GridapGmsh is an ASCENT dep):
#   julia --project=<2026_066/ASCENT> validation/enclosure_literal/handshake.jl
using GridapGmsh: GmshDiscreteModel
using Gridap.Geometry: get_face_labeling, num_tags, get_tag_name, num_cells

path = joinpath(@__DIR__, "enclosure_literal.msh")
isfile(path) || error("run reconstruct.jl (Tessella env) first to produce $path")
model = try; GmshDiscreteModel(path; renumber=false); catch; GmshDiscreteModel(path); end
labels = get_face_labeling(model)
names = [String(get_tag_name(labels, i)) for i in 1:num_tags(labels)]
vols = intersect(names, ["coax_pin", "air", "case"])
println("LITERAL ENC-COAX → ASCENT: num_cells=", num_cells(model), " volumes=", vols)
length(vols) == 3 || error("INCOMPLETE: volumes=$vols")
println("LITERAL_HANDSHAKE_OK — the literal enclosure meshes natively and loads in ASCENT")
