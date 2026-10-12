# `model.occ` kernel-side queries: `get_bounding_box`,
# `get_entities_in_bounding_box`, `get_curve_loops`, `get_surface_loops`
# and `occ.mesh.set_size`. Upstream queries the bound OCC registry
# without synchronizing — these read `kernel.shapes`/`kernel.bindings`
# directly, so rebound transforms and fresh shapes are visible before
# `synchronize`. Boxes are exact (line endpoints, circle/ellipse/torus
# extrema, the sphere's ball box, NURBS hulls padded by 1e-7 like
# `Bnd_Box`, boundary unions elsewhere); loop traversal order comes from
# the wire/shell `uses`, matching `BRepTools_WireExplorer` ordering.
# Layout mirrors `OCC_Internals::getBoundingBox`/`getEntitiesInBoundingBox`/
# `getCurveLoops`/`getSurfaceLoops`/`setMeshSize` on Gmsh 4.15.2.

using Test
using Tessella
const A=Tessella.API
const M=Tessella.Model

function _occ_setup()
    A.initialize()
    return A.model.occ
end

_current_model()=A._with_model(identity)

function _errtext(f)
    try
        f()
    catch err
        return sprint(showerror,err)
    end
    return ""
end

@testset "occ.get_bounding_box per dimension" begin
    occ=_occ_setup()
    v=occ.add_box(0.0,0.0,0.0, 2.0,1.0,3.0)
    occ.synchronize()
    # Solid, face, wire and vertex boxes recurse through the kernel graph.
    @test collect(occ.get_bounding_box(3,v))≈[0,0,0,2,1,3]
    @test collect(occ.get_bounding_box(2,1))≈[0,0,0,2,1,0]
    @test collect(occ.get_bounding_box(-1,1))≈[0,0,0,2,1,0]
    @test collect(occ.get_bounding_box(0,1))≈[0,0,0,0,0,0]
    # The upstream dispatch accepts every bound dimension.
    @test collect(occ.get_bounding_box(-2,1))≈[0,0,0,2,1,3]
    # Exact curved-record boxes without synchronization.
    occ2=_occ_setup()
    s=occ2.add_sphere(1.0,2.0,3.0, 2.0)
    @test collect(occ2.get_bounding_box(3,s))≈[-1,0,1,3,4,5]
    occ2.synchronize()
    # The sphere's single face carries the analytic ball box — its
    # meridian wire alone would not cover the poles.
    @test collect(occ2.get_bounding_box(2,1))≈[-1,0,1,3,4,5]
    t=occ2.add_torus(0.0,0.0,0.0, 4.0,1.0)
    @test collect(occ2.get_bounding_box(3,t))≈[-5,-5,-1,5,5,1] atol=1e-9
    c=occ2.add_circle(0.0,0.0,0.0,2.0)
    @test collect(occ2.get_bounding_box(1,c))≈[-2,-2,0,2,2,0]
    e=occ2.add_ellipse(0.0,0.0,0.0,3.0,1.0)
    @test collect(occ2.get_bounding_box(1,e))≈[-3,-1,0,3,1,0]
    # Unknown entity and bad dims reject with the upstream message.
    @test occursin("Unknown OpenCASCADE entity of dimension 3 with tag 99",
                   _errtext(()->occ.get_bounding_box(3,99)))
    @test occursin("dimension",_errtext(()->occ.get_bounding_box(4,1)))
end

@testset "occ.get_bounding_box sees transforms before sync" begin
    occ=_occ_setup()
    v=occ.add_cylinder(0.0,0.0,0.0, 0.0,0.0,2.0, 1.0)
    occ.synchronize()
    occ.translate([(3,v)],5.0,0.0,0.0)
    # Kernel-side query — no synchronize needed, matching upstream's
    # direct `_find` on the bound registry.
    @test collect(occ.get_bounding_box(3,v))≈[4,-1,0,6,1,2]
    occ.synchronize()
    @test collect(occ.get_bounding_box(3,v))≈[4,-1,0,6,1,2]
end

@testset "occ.get_entities_in_bounding_box" begin
    occ=_occ_setup()
    v=occ.add_box(0.0,0.0,0.0, 1.0,1.0,1.0)
    p=occ.add_point(10.0,0.0,0.0)
    q=occ.add_point(11.0,0.0,0.0)
    l=occ.add_line(p,q)
    occ.synchronize()
    # Containment semantics: an entity qualifies only when its box is
    # inside the query box — the standalone line sits outside.
    out=occ.get_entities_in_bounding_box(-0.5,-0.5,-0.5, 1.5,1.5,1.5)
    @test length(out)==27
    @test (1,l) ∉ out && (0,p) ∉ out
    @test (3,v) in out
    # Dimension filters keep the same ordering structure.
    @test occ.get_entities_in_bounding_box(-0.5,-0.5,-0.5, 1.5,1.5,1.5;dim=3)==
        [(3,v)]
    @test occ.get_entities_in_bounding_box(-0.5,-0.5,-0.5, 1.5,1.5,1.5;dim=0)==
        [(0,t) for t in 1:8]
    # A boundary-inclusive query picks up only the exact-touching entities.
    @test occ.get_entities_in_bounding_box(10.0,0.0,0.0, 11.0,0.0,0.0;dim=1)==
        [(1,l)]
    # Reversed boxes select nothing.
    @test isempty(occ.get_entities_in_bounding_box(5.0,0.0,0.0, 4.0,0.0,0.0))
    @test occursin("dimension",_errtext(()->
        occ.get_entities_in_bounding_box(0.,0.,0.,1.,1.,1.;dim=9)))
end

@testset "occ.get_curve_loops" begin
    occ=_occ_setup()
    # Plane face with a hole: two wires, outer first, traversal order.
    outer_p=[occ.add_point(x,y,0.0) for (x,y) in [(0,0),(4,0),(4,4),(0,4)]]
    outer_e=[occ.add_line(outer_p[a],outer_p[b])
             for (a,b) in [(1,2),(2,3),(3,4),(4,1)]]
    w1=occ.add_wire(outer_e)
    hole_p=[occ.add_point(x,y,0.0) for (x,y) in [(1,1),(3,1),(3,3),(1,3)]]
    hole_e=[occ.add_line(hole_p[a],hole_p[b])
            for (a,b) in [(1,2),(2,3),(3,4),(4,1)]]
    w2=occ.add_wire(hole_e)
    f=occ.add_plane_surface([w1,w2])
    loops,curves=occ.get_curve_loops(f)
    @test loops==[w1,w2]
    @test curves[1]==outer_e
    @test curves[2]==hole_e
    # The query is read-only — repeated calls return the same bound tags.
    loops2,curves2=occ.get_curve_loops(f)
    @test loops2==loops && curves2==curves
    # Unknown surface rejects with the upstream message.
    @test occursin("Unknown OpenCASCADE surface with tag 7",
                   _errtext(()->occ.get_curve_loops(7)))
end

@testset "occ.get_surface_loops" begin
    occ=_occ_setup()
    outer=occ.add_box(0.0,0.0,0.0, 4.0,4.0,4.0)
    loops,surfaces=occ.get_surface_loops(outer)
    @test loops==[1]
    @test sort(surfaces[1])==collect(1:6)
    # A standalone loop/volume pair reports its bound shell and faces.
    w=occ.add_wire([occ.add_line(a,b) for (a,b) in
        [(occ.add_point(0.,0.,-2.),occ.add_point(1.,0.,-2.)),
         (occ.add_point(1.,0.,-2.),occ.add_point(1.,1.,-2.)),
         (occ.add_point(1.,1.,-2.),occ.add_point(0.,1.,-2.)),
         (occ.add_point(0.,1.,-2.),occ.add_point(0.,0.,-2.))]])
    f=occ.add_plane_surface([w])
    s2=occ.add_surface_loop([f];sewing=false)
    v2=occ.add_volume([s2])
    loops2,surfaces2=occ.get_surface_loops(v2)
    @test loops2==[s2]
    @test surfaces2==[[f]]
    # Unknown volume rejects with the upstream message.
    @test occursin("Unknown OpenCASCADE volume with tag 42",
                   _errtext(()->occ.get_surface_loops(42)))
end

@testset "occ.mesh.set_size" begin
    occ=_occ_setup()
    p=occ.add_point(0.0,0.0,0.0)
    q=occ.add_point(1.0,0.0,0.0)
    occ.mesh.set_size([(0,p)],0.25)
    occ.synchronize()
    m=_current_model()
    @test m.point_size[p]≈0.25
    # Non-point dims and unbound points are silently ignored, matching
    # `OCC_Internals::setMeshSize`'s early return and `IsBound` guard.
    occ.mesh.set_size([(3,999),(0,4242)],7.0)
    occ.synchronize()
    @test m.point_size[p]≈0.25
    # A later call overwrites the constraint.
    occ.mesh.set_size([(0,p)],0.5)
    occ.synchronize()
    @test m.point_size[p]≈0.5
    @test m.point_size[q]==0.0
    @test occursin("finite",_errtext(()->occ.mesh.set_size([(0,p)],NaN)))
end
