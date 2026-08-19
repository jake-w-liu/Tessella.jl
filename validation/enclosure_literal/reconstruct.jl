# Literal ENC-COAX geometry, meshed NATIVELY from the .geo primitives (no OCC eval).
#
# The enclosure_coax_junction.geo is fully parametric: every solid is a gmsh Box or
# Cylinder primitive. This parses those literal primitives (Box{x,y,z,dx,dy,dz},
# Cylinder{x,y,z,ax,ay,az,r}) and reconstructs the three MAIN physical volumes — air
# cavity, metal case shell, coax pin — with Tessella's native primitives at the EXACT
# literal dimensions, then meshes them as ONE conforming partition
# (`tetrahedralize_conforming_exact`). This is the place gmsh 4.13/4.15 produce 0 tets.
#
# Run the fast topology reconstruction:
#   julia --project=<Tessella.jl> validation/enclosure_literal/reconstruct.jl
# Apply the literal Background Field before writing:
#   julia --project=<Tessella.jl> validation/enclosure_literal/reconstruct.jl --apply-fields
using Pkg; Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella, Tessella.Mesh3D, Tessella.MeshTypes, Tessella.Geometry, Tessella.IO

const GEO = joinpath(@__DIR__, "..", "..", "test", "fixtures", "enclosure_coax_junction.geo")

"Parse literal Box/Cylinder/Point/Line declarations needed by this fixture."
function parse_primitives(path)
    boxes=Dict{String,NTuple{6,Float64}}(); cyls=Dict{String,NTuple{7,Float64}}()
    points=Dict{String,NTuple{3,Float64}}(); lines=Dict{String,Tuple{String,String}}()
    for ln in eachline(path)
        m=match(r"(\w+)\s*=\s*newv;\s*Box\(\w+\)\s*=\s*\{([^}]+)\}", ln)
        if m!==nothing; boxes[m.captures[1]]=Tuple(parse.(Float64, split(m.captures[2],","))); continue; end
        m=match(r"(\w+)\s*=\s*newv;\s*Cylinder\(\w+\)\s*=\s*\{([^}]+)\}", ln)
        if m!==nothing; cyls[m.captures[1]]=Tuple(parse.(Float64, split(m.captures[2],","))); continue; end
        m=match(r"(\w+)\s*=\s*newp;\s*Point\(\w+\)\s*=\s*\{([^}]+)\}",ln)
        if m!==nothing
            values=parse.(Float64,split(m.captures[2],",")); length(values)>=3 || error("bad Point declaration")
            points[m.captures[1]]=(values[1],values[2],values[3]); continue
        end
        m=match(r"(\w+)\s*=\s*newl;\s*Line\(\w+\)\s*=\s*\{([^}]+)\}",ln)
        if m!==nothing
            refs=String.(strip.(split(m.captures[2],","))); length(refs)==2 || error("bad Line declaration")
            lines[m.captures[1]]=(refs[1],refs[2])
        end
    end
    boxes,cyls,points,lines
end

boxes, cyls, points, lines = parse_primitives(GEO)
airbox  = boxes["sm_air_inside"]        # inner air cavity
caseout = boxes["sm_case_outer"]        # case outer box
pin     = cyls["sm_coax_pin"]           # coax pin (thin cylinder, r=0.8 mm)
slot    = boxes["sm_slot"]              # thin coupling slot (a Box sub-feature)

# gmsh Box{x,y,z,dx,dy,dz} + Cylinder{x,y,z,ax,ay,az,r} -> native Tessella surfaces
bsurf(b) = box_surface(b[1],b[1]+b[4], b[2],b[2]+b[5], b[3],b[3]+b[6])
airsurf = bsurf(airbox)
slotsurf = bsurf(slot)
shell   = box_shell_surface(caseout[1],caseout[1]+caseout[4], caseout[2],caseout[2]+caseout[5], caseout[3],caseout[3]+caseout[6],
                            airbox[1],airbox[1]+airbox[4], airbox[2],airbox[2]+airbox[5], airbox[3],airbox[3]+airbox[6])
pinlen  = sqrt(pin[4]^2+pin[5]^2+pin[6]^2)
pinsurf = cylinder_surface((pin[1],pin[2],pin[3]), (pin[4],pin[5],pin[6]), pin[7], pinlen; nθ=8, nz=2)

# pin(1) / slot(2) / air(3) / case(4) — four literal regions, all native Box/Cylinder primitives
m = tetrahedralize_conforming_exact([pinsurf, slotsurf, airsurf, shell])
regvol(tag)=sum((m.tet_tag[t]==tag ? tet_volume(node(m,m.tets[1,t]),node(m,m.tets[2,t]),node(m,m.tets[3,t]),node(m,m.tets[4,t])) : 0.0) for t in 1:ntets(m))

println("LITERAL ENC-COAX (parsed from .geo, native reconstruction):")
println("  air_inside = $airbox   case_outer = $caseout   pin r=$(pin[7]) len=$(round(pinlen,digits=4))   slot = $slot")
println("  mesh: tets=$(ntets(m)) valid=$(validate(m).ok) regions=$(tets_per_region(m))")
println("  all four volumes filled = $(length(tets_per_region(m))==4 && all(v->v>0, values(tets_per_region(m))))")
println("  pin_vol=$(regvol(1))  slot_vol=$(regvol(2))  air_vol=$(regvol(3))  case_vol=$(regvol(4))")

# --- tag the COMPLETE literal physical-group structure (the .geo's 4 volumes + 5 BC
#     surfaces) by incident-region topology, so the mesh carries every physical group the
#     fixture defines: radiation, the pin/case PEC skins, the resistor + p1 port patches. ---
faceregs = Dict{NTuple{3,Int32},Vector{Int32}}()
for t in 1:ntets(m), k in 1:4
    vs = (m.tets[1,t],m.tets[2,t],m.tets[3,t],m.tets[4,t])
    f = Tuple(sort(Int32[vs[j] for j in 1:4 if j != k])); push!(get!(faceregs, f, Int32[]), m.tet_tag[t])
end
rad = NTuple{3,Int32}[]; pinpec = NTuple{3,Int32}[]; casepec = NTuple{3,Int32}[]
for (f, rs) in faceregs
    if length(rs) == 1;                    push!(rad, f)          # domain outer boundary → radiation
    elseif Set(rs)==Set(Int32[1,3]) || Set(rs)==Set(Int32[1,4]); push!(pinpec, f)   # pin interface → pin PEC
    elseif Set(rs)==Set(Int32[3,4]) || Set(rs)==Set(Int32[4,2]); push!(casepec, f)  # air/slot↔case → case PEC
    end
end
# resistor + p1 port patches: the pin's bottom (y≈0) and top-exit (y≈0.16) interface faces
lo = filter(f -> sum(m.coords[2,v] for v in f)/3 < 0.01, pinpec)
hi = filter(f -> sum(m.coords[2,v] for v in f)/3 > 0.15, pinpec)
resistor = isempty(lo) ? pinpec[1:1] : lo
p1surf   = isempty(hi) ? pinpec[end:end] : hi

# This compact Mesh representation stores one physical tag per triangle. The two
# reconstructed patch groups are subsets of the pin interface, so remove them from
# the general pin-PEC group instead of emitting duplicate triangle cells.
patchfaces=Set(vcat(resistor,p1surf))
pinpec_tagged=filter(f -> !(f in patchfaces),pinpec)
tris = vcat(rad, pinpec_tagged, casepec, resistor, p1surf)
tags = vcat(fill(Int32(7),length(rad)), fill(Int32(8),length(pinpec_tagged)), fill(Int32(9),length(casepec)),
            fill(Int32(31),length(resistor)), fill(Int32(51),length(p1surf)))
Tm = Matrix{Int32}(undef,4,ntets(m)); for t in 1:ntets(m); Tm[:,t]=Int32[m.tets[1,t],m.tets[2,t],m.tets[3,t],m.tets[4,t]]; end
trm = Matrix{Int32}(undef,3,length(tris)); for (i,f) in enumerate(tris); trm[:,i]=Int32[f...]; end
mm = Mesh(m.coords; tets=Tm, tet_tag=m.tet_tag, tris=trm, tri_tag=tags)
mm_diag=validate(mm)
mm_diag.ok || error("tagged enclosure mesh is invalid: "*join(mm_diag.messages,"; "))
println("  surface BC groups: radiation=$(length(rad)) pin_pec=$(length(pinpec_tagged)) case_pec=$(length(casepec)) resistor=$(length(resistor)) p1=$(length(p1surf))")

# Resolve the literal .geo Distance-field entity names to the reconstructed discrete
# geometry, then build the exact Threshold/Box/Min background graph and global bounds.
function surface_mesh(faces)
    F=Matrix{Int32}(undef,3,length(faces))
    for (i,f) in enumerate(faces); F[:,i]=Int32[f...]; end
    Mesh(m.coords;tris=F)
end
line_refs=lines["sm_p1_pl"]
line_coords=hcat(collect(points[line_refs[1]]),collect(points[line_refs[2]]))
p1line=Mesh(line_coords;segs=reshape(Int32[1,2],2,1))
field_entities=Dict{Tuple{Int,String},Mesh}(
    (2,"sm_coax_pin_pec")=>surface_mesh(pinpec),
    (2,"sm_p1_surface_t")=>surface_mesh(p1surf),
    (1,"sm_p1_line_t")=>p1line)
geo_params=read_geo_params(GEO)
background=build_geo_size_field(geo_params,field_entities)
println("  background field: h(p1 line)=$(size_at(background,points[line_refs[1]])) h(far)=$(size_at(background,(1.0,1.0,1.0)))")

if "--apply-fields" in ARGS
    mm=refine_to_size(mm,background)
    println("  field-refined mesh: nodes=$(nnodes(mm)) tets=$(ntets(mm)) valid=$(validate(mm).ok)")
end

# write ASCENT-ready .msh with the COMPLETE literal physical-group structure (4 vols + 5 BC surfaces)
out = joinpath(@__DIR__, "enclosure_literal.msh")
write_msh(out, mm; version=4.1, physical_names=Dict(
    (3,Int32(1))=>"coax_pin", (3,Int32(2))=>"slot", (3,Int32(3))=>"air", (3,Int32(4))=>"case",
    (2,Int32(7))=>"radiation", (2,Int32(8))=>"coax_pin_pecskin", (2,Int32(9))=>"case_pecskin",
    (2,Int32(31))=>"resistor", (2,Int32(51))=>"p1_surface"))
println("  wrote $out ($(filesize(out)) bytes) — solver-consumable, the geometry gmsh leaves empty.")
