# ASCENT side of the BC handshake: confirm a Tessella mesh carrying material volumes
# AND boundary-condition surfaces loads fully into ASCENT's parser (GridapGmsh), with
# every physical group — 3-D materials and 2-D BCs — visible to the solver.
#   julia --project=<2026_066/ASCENT> validation/ascent_handshake/bc_handshake.jl
using GridapGmsh: GmshDiscreteModel
using Gridap.Geometry: get_face_labeling, num_tags, get_tag_name, num_cells

path = joinpath(@__DIR__, "enc_coax_bc.msh")
isfile(path) || error("run generate_bc.jl (Tessella env) first to produce $path")
model = try; GmshDiscreteModel(path; renumber=false); catch; GmshDiscreteModel(path); end

labels = get_face_labeling(model); nT = num_tags(labels)
names = [String(get_tag_name(labels, i)) for i in 1:nT]
vols  = intersect(names, ["coax_pin", "air", "case"])
surfs = intersect(names, ["radiation", "coax_pin_pec"])
println("num_cells=", num_cells(model), " num_tags=", nT)
println("materials (volumes): ", vols)
println("boundary conditions (surfaces): ", surfs)
if length(vols) == 3 && length(surfs) == 2
    println("BC_HANDSHAKE_OK — ASCENT sees all 3 material volumes + 2 BC surfaces")
else
    error("INCOMPLETE: volumes=$vols surfaces=$surfs")
end
