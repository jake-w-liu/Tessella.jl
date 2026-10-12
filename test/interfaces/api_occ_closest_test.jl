# `model.occ.get_closest_entities` — `OCC_Internals::getClosestEntities`
# semantics on the pending-shape kernel (no `synchronize`): per input
# dim-tag the closest point is the query vertex's `BRepExtrema` answer
# (solid inputs short-circuit to zero distance at the query point when
# inside/on), results sort through the upstream `std::set` key
# `(distance, (dim, tag, x, y, z))` — duplicate inputs collapse — and
# `n <= 0` still emits the first entry like upstream's `nn >= n` loop.
# The kernel materializes each queried shape (with its wire/edge/vertex
# closure) into a private scratch `GeoModel`, so the verified `_dist_*`
# projectors apply verbatim and public model state stays untouched.
# Layout mirrors `OCC_Internals::getClosestEntities` on Gmsh 4.15.2.

using Test
using Tessella
const A=Tessella.API
const M=Tessella.Model

function _occ_setup()
    A.initialize()
    return A.model.occ
end

_errtext(f)=(try;f();catch err;return sprint(showerror,err);end;"")

@testset "occ.get_closest_entities analytic distances" begin
    occ=_occ_setup()
    b=occ.add_box(0.0,0.0,0.0, 2.0,3.0,4.0)
    s=occ.add_sphere(10.0,0.0,0.0, 1.5)
    ents=[(3,b),(3,s)]
    # Box face x=2 wins; closest point is the face projection.
    tags,dist,coord=occ.get_closest_entities(4.0,1.0,2.0,ents)
    @test tags==[(3,b)] && dist≈[2.0] && coord≈[2.0,1.0,2.0]
    # Both entities sorted by distance; the sphere's closest point is the
    # radial projection c + r·(q−c)/|q−c| = (10−9/√41, 1.5/√41, 3/√41).
    tags,dist,coord=occ.get_closest_entities(4.0,1.0,2.0,ents,2)
    @test tags==[(3,b),(3,s)]
    @test dist≈[2.0,sqrt(41)-1.5]
    @test coord≈[2.0,1.0,2.0,10-9/sqrt(41),1.5/sqrt(41),3/sqrt(41)]
    # Inside the box → zero distance at the query point (upstream
    # classifier IN/ON semantics), not a boundary projection.
    tags,dist,coord=occ.get_closest_entities(1.0,1.0,1.0,ents,2)
    @test tags[1]==(3,b) && dist[1]==0.0 && coord[1:3]≈[1.0,1.0,1.0]
    @test dist[2]≈sqrt(83)-1.5
    # On the boundary also reports zero.
    tags,dist,_=occ.get_closest_entities(2.0,1.0,2.0,[(3,b)])
    @test dist[1]==0.0
end

@testset "occ.get_closest_entities ordering and n" begin
    occ=_occ_setup()
    b=occ.add_box(0.0,0.0,0.0, 2.0,3.0,4.0)
    s=occ.add_sphere(10.0,0.0,0.0, 1.5)
    ents=[(3,b),(3,s)]
    # Duplicate inputs collapse through the upstream std::set key.
    tags,_,_=occ.get_closest_entities(4.0,1.0,2.0,[(3,b),(3,b),(3,s)],2)
    @test tags==[(3,b),(3,s)]
    # n <= 0 still returns the single closest (upstream nn>=n quirk).
    for n in (0,-3)
        tags,_,_=occ.get_closest_entities(4.0,1.0,2.0,ents,n)
        @test tags==[(3,b)]
    end
    # n beyond the candidate count returns all of them.
    tags,_,_=occ.get_closest_entities(4.0,1.0,2.0,ents,9)
    @test tags==[(3,b),(3,s)]
    # Empty input yields empty outputs.
    tags,dist,coord=occ.get_closest_entities(4.0,1.0,2.0,Tuple{Int,Int}[])
    @test isempty(tags) && isempty(dist) && isempty(coord)
    # Ties order by distance first, then (dim, tag): a vertex and two
    # box edges all at distance √2 from (3,0,1).
    p=occ.add_point(2.0,0.0,0.0)
    tags,dist,_=occ.get_closest_entities(3.0,0.0,1.0,
        [(1,10),(0,p),(1,2),(1,1)],4)
    @test tags[1]==(1,10) && dist[1]≈1.0
    @test (0,p) in tags && all(≈(sqrt(2.0)),dist[2:4])
end

@testset "occ.get_closest_entities wire, shell and face dims" begin
    occ=_occ_setup()
    p1=occ.add_point(5.0,0.0,0.0);p2=occ.add_point(5.0,3.0,0.0)
    p3=occ.add_point(5.0,3.0,3.0);p4=occ.add_point(5.0,0.0,3.0)
    c1=occ.add_line(p1,p2);c2=occ.add_line(p2,p3)
    c3=occ.add_line(p3,p4);c4=occ.add_line(p4,p1)
    w=occ.add_wire([c1,c2,c3,c4])
    # A bound wire resolves through its member edges; the output keeps
    # the wire's own (-1, tag).
    tags,dist,coord=occ.get_closest_entities(4.0,0.5,0.5,[(-1,w)])
    @test tags==[(-1,w)] && dist≈[sqrt(1.25)]
    @test coord[2]≈0.5 && coord[1]≈5.0
    # Standalone plane face: the interior projection wins over boundary.
    f=occ.add_plane_surface([w])
    tags,dist,coord=occ.get_closest_entities(4.0,0.5,0.5,[(2,f)])
    @test tags==[(2,f)] && dist≈[1.0] && coord≈[5.0,0.5,0.5]
    # A bound shell resolves through its member faces.
    b=occ.add_box(0.0,0.0,0.0, 2.0,3.0,4.0)
    loops,_=occ.get_surface_loops(b)
    tags,dist,coord=occ.get_closest_entities(4.0,1.0,2.0,[(-2,loops[1])])
    @test tags==[(-2,loops[1])] && dist≈[2.0] && coord≈[2.0,1.0,2.0]
end

@testset "occ.get_closest_entities trimmed faces and curved records" begin
    occ=_occ_setup()
    # Plane face with a square hole — the interior projection of a point
    # over the hole is rejected and the hole's boundary edge wins.
    pts=Int[]
    for (x,y) in [(0,0),(4,0),(4,4),(0,4),
                  (1.5,1.5),(2.5,1.5),(2.5,2.5),(1.5,2.5)]
        push!(pts,occ.add_point(Float64(x),Float64(y),0.0))
    end
    outer=occ.add_wire([occ.add_line(pts[i],pts[mod(i,4)+1]) for i in 1:4])
    inner=occ.add_wire(
        [occ.add_line(pts[i],pts[mod(i-4,4)+5]) for i in 5:8])
    f=occ.add_plane_surface([outer,inner])
    tags,dist,coord=occ.get_closest_entities(2.0,2.0,1.0,[(2,f)])
    @test tags==[(2,f)] && dist≈[sqrt(1.25)]
    @test coord[3]==0.0 && (coord[1]≈2.0 && coord[2]≈1.5 ||
                            coord[1]≈2.0 && coord[2]≈2.5 ||
                            coord[1]≈1.5 && coord[2]≈2.0 ||
                            coord[1]≈2.5 && coord[2]≈2.0)
    # Torus face: closest point is the analytic tube projection
    # (3,0,0)+ (1,0,3)/√10 → distance √10−1.
    occ2=_occ_setup()
    t=occ2.add_torus(0.0,0.0,0.0, 3.0,1.0)
    fe=A.model.occ.get_entities(2)
    tags,dist,coord=occ2.get_closest_entities(4.0,0.0,3.0,fe)
    @test dist[1]≈sqrt(10)-1
    @test coord[1:3]≈[3.0+1/sqrt(10),0.0,3/sqrt(10)]
    # Copied-curve families (Bezier interpolant) project correctly.
    occ3=_occ_setup()
    b1=occ3.add_point(0.0,0.0,40.0);b2=occ3.add_point(1.0,2.0,40.0)
    b3=occ3.add_point(2.0,0.0,40.0)
    bz=occ3.add_bezier([b1,b2,b3])
    tags,dist,coord=occ3.get_closest_entities(1.0,0.0,40.0,[(1,bz)])
    @test dist[1]≈sqrt(0.75)
end

@testset "occ.get_closest_entities transform and state isolation" begin
    occ=_occ_setup()
    b=occ.add_box(0.0,0.0,0.0, 2.0,3.0,4.0)
    # Kernel-side transforms are visible without synchronize.
    occ.translate([(3,b)],10.0,0.0,0.0)
    tags,dist,coord=occ.get_closest_entities(4.0,1.0,2.0,[(3,b)])
    @test dist==[6.0] && coord≈[10.0,1.0,2.0]
    # The query mutates nothing: bindings, model entities and maxima are
    # identical before and after.
    m=A._with_model(identity)
    before=(copy(m.occ_kernel.bindings),copy(m.occ_kernel.shapes),
            copy(m.points),copy(m.curves),copy(m.surfaces),
            copy(m.volumes),copy(m.occ_kernel.maxima))
    occ.get_closest_entities(1.0,1.0,1.0,
        vcat(A.model.occ.get_entities(0),A.model.occ.get_entities(1),
             A.model.occ.get_entities(2),A.model.occ.get_entities(3)),99)
    @test m.occ_kernel.bindings==before[1]
    @test m.occ_kernel.shapes==before[2]
    @test m.points==before[3] && m.curves==before[4]
    @test m.surfaces==before[5] && m.volumes==before[6]
    @test m.occ_kernel.maxima==before[7]
    # Errors before any resolution leave state untouched as well.
    _errtext(()->occ.get_closest_entities(0.0,0.0,0.0,[(3,99)]))
    @test m.occ_kernel.bindings==before[1]
    @test m.occ_kernel.shapes==before[2]
end

@testset "occ.get_distance entity pairs" begin
    occ=_occ_setup()
    b=occ.add_box(0.0,0.0,0.0, 2.0,3.0,4.0)
    s=occ.add_sphere(10.0,0.0,0.0, 1.5)
    p=occ.add_point(6.0,0.0,0.0)
    # Upstream out-param layout: (distance, x1,y1,z1, x2,y2,z2).
    @test collect(occ.get_distance(0,p,3,b))≈[4.0,6.0,0.0,0.0,2.0,0.0,0.0]
    @test collect(occ.get_distance(0,p,3,s))≈[2.5,6.0,0.0,0.0,8.5,0.0,0.0]
    @test collect(occ.get_distance(3,b,3,s))≈[6.5,2.0,0.0,0.0,8.5,0.0,0.0]
    # Argument order maps to shape1/shape2 points.
    d=occ.get_distance(3,b,0,p)
    @test d[1]≈4.0 && d[2:4]==(2.0,0.0,0.0) && d[5:7]==(6.0,0.0,0.0)
    # Inside a volume reports zero.
    pin=occ.add_point(1.0,1.0,1.0)
    d=occ.get_distance(0,pin,3,b)
    @test d[1]==0.0 && d[2:4]==(1.0,1.0,1.0)
    # A point on a face's support plane still answers the boundary
    # minimum (upstream's extrema never classifies).
    p2=occ.add_point(3.0,0.0,1.0)
    @test collect(occ.get_distance(0,p2,3,b))≈[1.0,3.0,0.0,1.0,2.0,0.0,1.0]
    # Shared subshapes between the inputs: a wire against its own
    # member edge, a shell against its own face.
    loops,curves=occ.get_curve_loops(1)
    @test occ.get_distance(-1,loops[1],1,curves[1][1])[1]==0.0
    sl,faces=occ.get_surface_loops(b)
    @test occ.get_distance(-2,sl[1],2,faces[1][1])[1]==0.0
    # Segment-segment closed form and curved-record pairs.
    l1p1=occ.add_point(0.0,0.0,10.0);l1p2=occ.add_point(2.0,0.0,10.0)
    l2p1=occ.add_point(1.0,-1.0,10.0);l2p2=occ.add_point(1.0,1.0,12.0)
    l1=occ.add_line(l1p1,l1p2);l2=occ.add_line(l2p1,l2p2)
    d=occ.get_distance(1,l1,1,l2)
    @test d[1]≈sqrt(0.5) && d[2:4]==(1.0,0.0,10.0) &&
          d[5:7]==(1.0,-0.5,10.5)
    c=occ.add_circle(0.0,0.0,15.0,1.0)
    lp1=occ.add_point(0.0,0.0,20.0);lp2=occ.add_point(1.0,0.0,20.0)
    d=occ.get_distance(1,c,1,occ.add_line(lp1,lp2))
    @test d[1]≈5.0
    # Self distance is zero; member units cover a bound wire.
    @test occ.get_distance(3,b,3,b)[1]==0.0
    wpts=[occ.add_point(x,y,20.0) for (x,y) in
          [(0,0),(2,0),(2,2),(0,2)]]
    wes=[occ.add_line(wpts[1],wpts[2]),occ.add_line(wpts[2],wpts[3]),
         occ.add_line(wpts[3],wpts[4]),occ.add_line(wpts[4],wpts[1])]
    w=occ.add_wire(wes)
    d=occ.get_distance(-1,w,0,p)
    # Wire rim at z=20 to (6,0,0): closest member point (2,0,20).
    @test d[1]≈sqrt(416.0) && d[2:4]==(2.0,0.0,20.0) &&
          d[5:7]==(6.0,0.0,0.0)
end

@testset "occ.get_distance errors and isolation" begin
    occ=_occ_setup()
    occ.add_box(0.0,0.0,0.0, 1.0,1.0,1.0)
    @test occursin("Unknown OpenCASCADE entity of dimension 3 with tag 99",
        _errtext(()->occ.get_distance(3,99,3,1)))
    @test occursin("Unknown OpenCASCADE entity of dimension 2 with tag 99",
        _errtext(()->occ.get_distance(3,1,2,99)))
    @test occursin("OCC dimension must be in -2:3",
        _errtext(()->occ.get_distance(4,1,3,1)))
    @test occursin("tag must be an integer",
        _errtext(()->occ.get_distance(3,1.5,3,1)))
    # Queries leave registry and model state untouched, errors included.
    m=A._with_model(identity)
    before=(copy(m.occ_kernel.bindings),copy(m.occ_kernel.shapes),
            copy(m.points),copy(m.occ_kernel.maxima))
    occ.get_distance(3,1,3,1)
    _errtext(()->occ.get_distance(3,1,0,99))
    @test m.occ_kernel.bindings==before[1]
    @test m.occ_kernel.shapes==before[2]
    @test m.points==before[3] && m.occ_kernel.maxima==before[4]
end

@testset "occ.get_closest_entities errors" begin
    occ=_occ_setup()
    occ.add_box(0.0,0.0,0.0, 1.0,1.0,1.0)
    @test occursin("Unknown OpenCASCADE entity of dimension 3 with tag 99",
        _errtext(()->occ.get_closest_entities(0.0,0.0,0.0,[(3,99)])))
    @test occursin("Unknown OpenCASCADE entity of dimension 1 with tag 99",
        _errtext(()->occ.get_closest_entities(0.0,0.0,0.0,[(1,99)])))
    @test occursin("OCC dimension must be in -2:3",
        _errtext(()->occ.get_closest_entities(0.0,0.0,0.0,[(4,1)])))
    @test occursin("coordinates must be finite",
        _errtext(()->occ.get_closest_entities(NaN,0.0,0.0,[(3,1)])))
    @test occursin("coordinates must be finite",
        _errtext(()->occ.get_closest_entities(0.0,0.0,Inf,[(3,1)])))
    @test occursin("n must be an integer",
        _errtext(()->occ.get_closest_entities(0.0,0.0,0.0,[(3,1)],1.5)))
    @test occursin("tag must be an integer",
        _errtext(()->occ.get_closest_entities(0.0,0.0,0.0,[(3,1.5)])))
end
