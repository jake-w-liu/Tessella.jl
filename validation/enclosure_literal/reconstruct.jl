# Literal ENC-COAX geometry, meshed NATIVELY from the .geo primitives (no OCC eval).
#
# The enclosure_coax_junction.geo is fully parametric: every solid is a gmsh Box or
# Cylinder primitive. This parses those literal primitives (Box{x,y,z,dx,dy,dz},
# Cylinder{x,y,z,ax,ay,az,r}) and reconstructs the three MAIN physical volumes — air
# cavity, metal case shell, coax pin — with Tessella's native primitives at the EXACT
# literal dimensions, then meshes them as ONE conforming partition
# (`tetrahedralize_conforming_exact`). This is the place gmsh 4.13/4.15 produce 0 tets.
#
# Run:  julia --project=<Tessella.jl> validation/enclosure_literal/reconstruct.jl
using Pkg; Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella, Tessella.Mesh3D, Tessella.MeshTypes, Tessella.Geometry, Tessella.IO

const GEO = joinpath(@__DIR__, "..", "..", "test", "fixtures", "enclosure_coax_junction.geo")

"parse gmsh Box/Cylinder primitives from a .geo file into (boxes, cylinders) dicts."
function parse_primitives(path)
    boxes=Dict{String,NTuple{6,Float64}}(); cyls=Dict{String,NTuple{7,Float64}}()
    for ln in eachline(path)
        m=match(r"(\w+)\s*=\s*newv;\s*Box\(\w+\)\s*=\s*\{([^}]+)\}", ln)
        if m!==nothing; boxes[m.captures[1]]=Tuple(parse.(Float64, split(m.captures[2],","))); continue; end
        m=match(r"(\w+)\s*=\s*newv;\s*Cylinder\(\w+\)\s*=\s*\{([^}]+)\}", ln)
        m!==nothing && (cyls[m.captures[1]]=Tuple(parse.(Float64, split(m.captures[2],","))))
    end
    boxes,cyls
end

boxes, cyls = parse_primitives(GEO)
airbox  = boxes["sm_air_inside"]        # inner air cavity
caseout = boxes["sm_case_outer"]        # case outer box
pin     = cyls["sm_coax_pin"]           # coax pin (thin cylinder, r=0.8 mm)

# gmsh Box{x,y,z,dx,dy,dz} + Cylinder{x,y,z,ax,ay,az,r} -> native Tessella surfaces
airsurf = box_surface(airbox[1],airbox[1]+airbox[4], airbox[2],airbox[2]+airbox[5], airbox[3],airbox[3]+airbox[6])
shell   = box_shell_surface(caseout[1],caseout[1]+caseout[4], caseout[2],caseout[2]+caseout[5], caseout[3],caseout[3]+caseout[6],
                            airbox[1],airbox[1]+airbox[4], airbox[2],airbox[2]+airbox[5], airbox[3],airbox[3]+airbox[6])
pinlen  = sqrt(pin[4]^2+pin[5]^2+pin[6]^2)
pinsurf = cylinder_surface((pin[1],pin[2],pin[3]), (pin[4],pin[5],pin[6]), pin[7], pinlen; nθ=8, nz=2)

m = tetrahedralize_conforming_exact([pinsurf, airsurf, shell])   # pin(1)/air(2)/case(3)
regvol(tag)=sum((m.tet_tag[t]==tag ? tet_volume(node(m,m.tets[1,t]),node(m,m.tets[2,t]),node(m,m.tets[3,t]),node(m,m.tets[4,t])) : 0.0) for t in 1:ntets(m))

println("LITERAL ENC-COAX (parsed from .geo, native reconstruction):")
println("  air_inside = $airbox   case_outer = $caseout   pin r=$(pin[7]) len=$(round(pinlen,digits=4))")
println("  mesh: tets=$(ntets(m)) valid=$(validate(m).ok) regions=$(tets_per_region(m))")
println("  all three volumes filled = $(all(v->v>0, values(tets_per_region(m))))")
println("  air_vol=$(regvol(2))  case_vol=$(regvol(3))  pin_vol=$(regvol(1))")

# write ASCENT-ready .msh with the three literal physical volume names
out = joinpath(@__DIR__, "enclosure_literal.msh")
write_msh(out, m; version=4.1, physical_names=Dict(
    (3,Int32(1))=>"coax_pin", (3,Int32(2))=>"air", (3,Int32(3))=>"case"))
println("  wrote $out ($(filesize(out)) bytes) — solver-consumable, the geometry gmsh leaves empty.")
