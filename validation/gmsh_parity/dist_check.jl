# model_distance differential vs gmsh 4.15.2 occ.getDistance oracle values.
# Each case mirrors dist_oracle.py; compare d within tolerance and verify the
# returned points attain it.
using Tessella
using Tessella.Model: model_distance

fails = String[]

function check(key, got, expected; atol=1e-9, pt_atol=1e-7)
    d, pa, pb = got
    de = expected
    rel = abs(d - de) <= atol * max(1.0, abs(de))
    attain = hypot(pa[1] - pb[1], pa[2] - pb[2], pa[3] - pb[3])
    rel && abs(attain - d) <= pt_atol * max(1.0, d) ||
        push!(fails, "$key: d=$d want=$de attain=$attain pa=$pa pb=$pb")
    nothing
end

# pt-pt
m = GeoModel()
add_point!(m, 0, 0, 0; tag=1); add_point!(m, 3, 4, 0; tag=2)
check("pt-pt", model_distance(m, 0, 1, 0, 2), 5.0)

# pt-line
m = GeoModel()
add_point!(m, 0, 0, 0; tag=1); add_point!(m, 2, 0, 0; tag=2)
add_line!(m, 1, 2; tag=1); add_point!(m, 1, 4, 0; tag=3)
check("pt-line", model_distance(m, 0, 3, 1, 1), 4.0)

# pt-arc: quarter circle about origin, pt at (0,0,3) -> sqrt(10)
m = GeoModel()
add_point!(m, 1, 0, 0; tag=1); add_point!(m, 0, 0, 0; tag=2)
add_point!(m, 0, 1, 0; tag=3); add_point!(m, 0, 0, 3; tag=4)
add_circle_arc!(m, 1, 2, 3; tag=1)
check("pt-arc", model_distance(m, 0, 4, 1, 1), 3.1622776601683795)

# pt-plane
function unit_square!(m, z, tbase)
    for (t, x, y) in ((tbase + 1, 0.0, 0.0), (tbase + 2, 1.0, 0.0),
                      (tbase + 3, 1.0, 1.0), (tbase + 4, 0.0, 1.0))
        add_point!(m, x, y, z; tag=t)
    end
    for (t, a, b) in ((tbase + 1, tbase + 1, tbase + 2),
                      (tbase + 2, tbase + 2, tbase + 3),
                      (tbase + 3, tbase + 3, tbase + 4),
                      (tbase + 4, tbase + 4, tbase + 1))
        add_line!(m, a, b; tag=t)
    end
    add_curve_loop!(m, [tbase + 1, tbase + 2, tbase + 3, tbase + 4]; tag=tbase + 1)
    add_plane_surface!(m, [tbase + 1]; tag=tbase + 1)
    return tbase + 1
end
m = GeoModel()
sf = unit_square!(m, 0.0, 0)
add_point!(m, 0.5, 0.5, 2; tag=9)
check("pt-plane", model_distance(m, 0, 9, 2, sf), 2.0)

# pt-box out/in
m = GeoModel()
add_box!(m, 0, 0, 0, 1, 1, 1; tag=1)
add_point!(m, 3, 0, 0; tag=50)
check("pt-box-out", model_distance(m, 0, 50, 3, 1), 2.0)
m = GeoModel()
add_box!(m, 0, 0, 0, 1, 1, 1; tag=1)
add_point!(m, 0.2, 0.5, 0.5; tag=50)
check("pt-box-in", model_distance(m, 0, 50, 3, 1), 0.0)

# line-line skew / crossing
m = GeoModel()
add_point!(m, 0, 0, 0; tag=1); add_point!(m, 1, 0, 0; tag=2)
add_point!(m, 0, 1, 3; tag=3); add_point!(m, 0, 1, 4; tag=4)
add_line!(m, 1, 2; tag=1); add_line!(m, 3, 4; tag=2)
check("line-line-skew", model_distance(m, 1, 1, 1, 2), 3.1622776601683795)
m = GeoModel()
add_point!(m, -1, 0, 0; tag=1); add_point!(m, 1, 0, 0; tag=2)
add_point!(m, 0, -1, 0; tag=3); add_point!(m, 0, 1, 0; tag=4)
add_line!(m, 1, 2; tag=1); add_line!(m, 3, 4; tag=2)
check("line-line-x", model_distance(m, 1, 1, 1, 2), 0.0)

# line-plane
m = GeoModel()
sf = unit_square!(m, 0.0, 0)
add_point!(m, 0.5, 0.5, 2; tag=9); add_point!(m, 0.5, 0.6, 2; tag=10)
add_line!(m, 9, 10; tag=9)
check("line-plane", model_distance(m, 1, 9, 2, sf), 2.0)

# face-face
m = GeoModel()
s1 = unit_square!(m, 0.0, 0)
s2 = unit_square!(m, 2.0, 10)
check("face-face", model_distance(m, 2, s1, 2, s2), 2.0)

# box-box
m = GeoModel()
add_box!(m, 0, 0, 0, 1, 1, 1; tag=1)
add_box!(m, 4, 0, 0, 1, 1, 1; tag=2)
check("box-box", model_distance(m, 3, 1, 3, 2), 3.0)

# box-sphere (gap 2)
m = GeoModel()
add_box!(m, 0, 0, 0, 1, 1, 1; tag=1)
add_sphere!(m, 4, 0.5, 0.5, 1.0; tag=2)
check("box-sphere", model_distance(m, 3, 1, 3, 2), 2.0; atol=1e-8)

# sphere-sphere (gap 1)
m = GeoModel()
add_sphere!(m, 0, 0, 0, 1.0; tag=1)
add_sphere!(m, 3, 0, 0, 1.0; tag=2)
check("sphere-sphere", model_distance(m, 3, 1, 3, 2), 1.0; atol=1e-8)

# overlapping box-sphere -> 0
m = GeoModel()
add_box!(m, 0, 0, 0, 2, 2, 2; tag=1)
add_sphere!(m, 1, 1, 3, 1.0; tag=2)
check("box-sphere-ovl", model_distance(m, 3, 1, 3, 2), 0.0; atol=1e-8)

# cyl-box
m = GeoModel()
add_cylinder!(m, 0, 0, 0, 0, 0, 3, 1.0; tag=1)
add_box!(m, 4, 0, 0, 1, 1, 1; tag=2)
check("cyl-box", model_distance(m, 3, 1, 3, 2), 3.0; atol=1e-8)

# pt-torus: R=3 r=1 at origin, center pt -> inner equator distance 2
m = GeoModel()
add_torus!(m, 0, 0, 0, 3.0, 1.0; tag=1)
add_point!(m, 0, 0, 0; tag=50)
check("pt-torus", model_distance(m, 0, 50, 3, 1), 2.0; atol=1e-7)

println("=== differential results ===")
if isempty(fails)
    println("ALL-DIFFERENTIAL-OK")
else
    for f in fails
        println("FAIL ", f)
    end
end
