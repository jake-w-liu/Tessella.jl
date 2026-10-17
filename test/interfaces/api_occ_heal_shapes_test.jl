# `model.occ.heal_shapes` — tolerance healing: coincident-but-distinct
# vertices weld within `tolerance`, line edges sharing both endpoints
# collapse to one curve entity, restoring the two-surfaces-per-curve
# shell rule for cracked solids (distinct-but-coincident geometry).
using Test
using Tessella
const A=Tessella.API

function _errtext(f)
    try
        f()
    catch err
        return sprint(showerror,err)
    end
    return ""
end

@testset "occ.heal_shapes" begin
    A.initialize(); occ=A.model.occ
    try
        # Unit cube with EVERY face authored over its own private
        # vertices — 24 coincident-but-distinct corner entities, 24
        # private edges; no curve is shared so the shell cannot
        # validate until healing welds them.
        faces=Int[]
        quads=[((0,0,0),(0,1,0),(1,1,0),(1,0,0)),   # -z
               ((0,0,1),(1,0,1),(1,1,1),(0,1,1)),   # +z
               ((0,0,0),(0,0,1),(0,1,1),(0,1,0)),   # -x
               ((1,0,0),(1,1,0),(1,1,1),(1,0,1)),   # +x
               ((0,0,0),(1,0,0),(1,0,1),(0,0,1)),   # -y
               ((0,1,0),(0,1,1),(1,1,1),(1,1,0))]   # +y
        for q in quads
            pts=[occ.add_point(Float64(p[1]),Float64(p[2]),Float64(p[3]))
                 for p in q]
            ls=[occ.add_line(pts[j],pts[mod1(j+1,4)]) for j in 1:4]
            w=occ.add_wire(ls)
            push!(faces,occ.add_plane_surface([w]))
        end
        occ.heal_shapes([(2,f) for f in faces])
        loop=occ.add_surface_loop(faces)
        v=occ.add_volume([loop])
        @test occ.get_mass(3,v)≈1.0 atol=1e-6
        c=_errtext(()->occ.get_center_of_mass(3,v))
        @test c==""
        # Welded cube: 24 private edges collapse to the 12 unique
        # geometric edges; 24 private vertices weld to 8. Count bound
        # OCC entities per dimension (welded duplicates are unbound).
        k=A._with_model(m->hasproperty(m,:occ_kernel) ? m.occ_kernel : m)
        cnt(d)=count(kv->first(kv.first)==d,k.bindings)
        @test cnt(1)==12
        @test cnt(0)==8

        # Second cracked cube with sewing disabled — edges stay
        # private entities (vertex welding alone cannot merge them).
        faces2=Int[]
        for q in quads
            pts=[occ.add_point(Float64(p[1]+2),Float64(p[2]),
                               Float64(p[3])) for p in q]
            ls=[occ.add_line(pts[j],pts[mod1(j+1,4)]) for j in 1:4]
            w=occ.add_wire(ls)
            push!(faces2,occ.add_plane_surface([w]))
        end
        occ.heal_shapes([(2,f) for f in faces2];
                        fixSmallEdges=false,sewFaces=false)
        @test cnt(1)==12+24

        # unknown entities reject; no-op on nothing is fine
        @test occursin("unknown",
                       _errtext(()->occ.heal_shapes([(2,99999)])))
    finally
        A.finalize()
    end
end
