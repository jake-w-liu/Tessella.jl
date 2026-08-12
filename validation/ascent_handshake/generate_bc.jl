# Tessella side of the ASCENT *boundary-condition* handshake: build the ENC-COAX
# conforming mesh (pin/air/case) AND tag its boundary + material-interface faces as
# surface physical groups, then write gmsh MSH v4.1 carrying BOTH 3-D material volumes
# and 2-D BC surfaces — the full solver-ready structure ASCENT needs.
#   julia --project=<Tessella.jl> validation/ascent_handshake/generate_bc.jl
using Pkg; Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella, Tessella.Mesh3D, Tessella.MeshTypes, Tessella.Geometry, Tessella.IO

outer  = box_surface(0.,6., 0.,6., 0.,6.)                       # case outer
cavity = box_surface(1.,5., 1.,5., 1.,5.)                       # air+pin cavity
pin    = cylinder_surface((3.,3.,1.5), (0.,0.,1.), 0.7, 3.0; nθ=8, nz=2)
m = tetrahedralize_conforming([pin, cavity, outer])             # pin(1)/air(2)/case(3)
validate(m).ok || error("generate_bc.jl: conforming mesh invalid")

# classify each interior face by the region tags incident to it
faceregs = Dict{NTuple{3,Int32},Vector{Int32}}()
for t in 1:ntets(m), k in 1:4
    vs = (m.tets[1,t],m.tets[2,t],m.tets[3,t],m.tets[4,t])
    f = Tuple(sort(Int32[vs[j] for j in 1:4 if j != k]))
    push!(get!(faceregs, f, Int32[]), m.tet_tag[t])
end
rad = NTuple{3,Int32}[]      # domain boundary (incident to one tet) -> radiation BC
pec = NTuple{3,Int32}[]      # pin(1)/air(2) interface -> coax_pin PEC skin
for (f, rs) in faceregs
    length(rs) == 1 && push!(rad, f)
    (length(rs) == 2 && Set(rs) == Set(Int32[1,2])) && push!(pec, f)
end

tris = vcat(rad, pec)
tritag = vcat(fill(Int32(10), length(rad)), fill(Int32(11), length(pec)))
Tm = Matrix{Int32}(undef, 4, ntets(m))
for t in 1:ntets(m); Tm[:,t] = Int32[m.tets[1,t],m.tets[2,t],m.tets[3,t],m.tets[4,t]]; end
trm = Matrix{Int32}(undef, 3, length(tris)); for (i,f) in enumerate(tris); trm[:,i] = Int32[f...]; end
mm = Mesh(m.coords; tets=Tm, tet_tag=m.tet_tag, tris=trm, tri_tag=tritag)

out = joinpath(@__DIR__, "enc_coax_bc.msh")
write_msh(out, mm; version=4.1, physical_names=Dict(
    (3,Int32(1))=>"coax_pin", (3,Int32(2))=>"air", (3,Int32(3))=>"case",
    (2,Int32(10))=>"radiation", (2,Int32(11))=>"coax_pin_pec"))
println("wrote $out  (tets=$(ntets(m)), 3 material volumes + 2 BC surfaces: ",
        "$(length(rad)) radiation, $(length(pec)) coax_pin_pec)")
