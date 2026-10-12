# `model.occ` ellipse endpoints: `add_ellipse` (center/radii/axis form, full
# and trimmed), `add_ellipse_arc` (start/center/major/end four-point form),
# and `add_disk` (full-ellipse wire + plane face). Ellipse edges carry the
# native `(occ=:ellipse)` record — center, frame axes, radii, trim — and
# plane faces store `ellips2d` pcurves, matching
# `OCC_Internals::addEllipse`/`addEllipseArc`/`makeDisk` semantics observed
# on Gmsh 4.15.2.

using Test
using LinearAlgebra
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

@testset "OCC add_ellipse full" begin
    occ=_occ_setup()
    e=occ.add_ellipse(1.0,2.0,3.0,4.0,2.0)
    @test e==1
    @test sort!(collect(occ.get_entities(1)))==[(1,1)]
    occ.synchronize()
    m=_current_model()
    # A full ellipse collapses to a single seam vertex.
    a,b=m.curves[1]
    @test a==b
    @test m.curve_types[1]==:ellipse
    g=m.curve_geometry[1]
    @test g.occ==:ellipse && g.rx==4.0 && g.ry==2.0
    @test g.center==(1.0,2.0,3.0) && g.n==(0.0,0.0,1.0)
    @test g.t0==0.0 && g.t1==2π
    @test collect(m.points[a])≈[5.0,2.0,3.0] atol=1e-12
    # Evaluation, derivative, curvature at the cardinal parameters.
    @test A.model.get_value(1,1,[0.0])≈[5.0,2.0,3.0] atol=1e-12
    @test A.model.get_value(1,1,[π/2])≈[1.0,4.0,3.0] atol=1e-12
    @test A.model.get_value(1,1,[π])≈[-3.0,2.0,3.0] atol=1e-12
    @test A.model.get_derivative(1,1,[0.0])≈[0.0,2.0,0.0] atol=1e-12
    @test A.model.get_derivative(1,1,[π/2])≈[-4.0,0.0,0.0] atol=1e-12
    # κ = a/b² at the major ends, b/a² at the minor ends.
    @test A.model.get_curvature(1,1,[0.0])[1]≈1.0 atol=1e-9
    @test A.model.get_curvature(1,1,[π/2])[1]≈0.125 atol=1e-9
    lo,hi=A.model.get_parametrization_bounds(1,1)
    @test lo==[0.0] && hi[1]==2π
    @test collect(A.model.get_bounding_box(1,1))≈[-3.0,0.0,3.0,5.0,4.0,3.0] atol=1e-12
    # An explicit tag is honored; a duplicate is rejected.
    @test occ.add_ellipse(0.0,0.0,0.0,3.0,1.5;tag=9)==9
    @test occursin("already exists",_errtext(()->occ.add_ellipse(0,0,0,3,1;tag=9)))
end

@testset "OCC add_ellipse trimmed" begin
    occ=_occ_setup()
    e=occ.add_ellipse(0.0,0.0,0.0,4.0,2.0;angle1=π/4,angle2=3π/4)
    occ.synchronize()
    m=_current_model()
    a,b=m.curves[1]
    @test a!=b
    g=m.curve_geometry[1]
    @test g.t0==π/4 && g.t1==3π/4
    s2=sqrt(0.5)
    @test collect(m.points[a])≈[4s2,2s2,0.0] atol=1e-12
    @test collect(m.points[b])≈[-4s2,2s2,0.0] atol=1e-12
    @test A.model.get_value(1,1,[g.t0])≈collect(m.points[a]) atol=1e-12
    @test A.model.get_value(1,1,[g.t1])≈collect(m.points[b]) atol=1e-12
    mid=A.model.get_value(1,1,[(g.t0+g.t1)/2])
    @test mid≈[0.0,2.0,0.0] atol=1e-12
    # Degenerate ranges reject like Geom_TrimmedCurve's U1==U2.
    @test occursin("U1 == U2",_errtext(()->occ.add_ellipse(0,0,0,4,2;angle1=0.5,angle2=0.5)))
end

@testset "OCC add_ellipse axis forms" begin
    occ=_occ_setup()
    # Normal-only: the major direction is whatever gp_Ax2 picks; check the
    # plane and radii, not the frame.
    e1=occ.add_ellipse(0.0,0.0,0.0,4.0,2.0;zAxis=[1.0,0.0,0.0])
    # Normal + x direction: major along y in the y-z plane.
    e2=occ.add_ellipse(0.0,0.0,0.0,4.0,2.0;zAxis=[1.0,0.0,0.0],xAxis=[0.0,1.0,0.0])
    occ.synchronize()
    m=_current_model()
    g1=m.curve_geometry[e1]; g2=m.curve_geometry[e2]
    @test abs(g1.n[1])≈1.0
    @test g2.n==(1.0,0.0,0.0) && g2.X==(0.0,1.0,0.0)
    @test A.model.get_value(1,e2,[0.0])≈[0.0,4.0,0.0] atol=1e-12
    @test A.model.get_value(1,e2,[π/2])≈[0.0,0.0,2.0] atol=1e-12
    for u in range(g1.t0,g1.t1,length=9)
        p=A.model.get_value(1,e1,[u])
        @test p[1]≈0.0 atol=1e-12
        rel=p.-collect(g1.center)
        @test (dot(rel,collect(g1.X))/g1.rx)^2+
              (dot(rel,collect(g1.Y))/g1.ry)^2≈1.0 atol=1e-12
    end
    @test occursin("parallel",_errtext(()->occ.add_ellipse(0,0,0,4,2;zAxis=[0,0,1],xAxis=[0,0,2])))
    @test occursin("null modulus",_errtext(()->occ.add_ellipse(0,0,0,4,2;zAxis=[0,0,0])))
    @test occursin("3-element",_errtext(()->occ.add_ellipse(0,0,0,4,2;zAxis=[0,0])))
end

@testset "OCC add_ellipse_arc" begin
    occ=_occ_setup()
    s=occ.add_point(4.0,0.0,0.0)
    c=occ.add_point(0.0,0.0,0.0)
    mj=occ.add_point(4.0,0.0,0.0)
    e=occ.add_point(0.0,2.0,0.0)
    a=occ.add_ellipse_arc(s,c,mj,e)
    @test a==1
    @test sort!(collect(occ.get_entities(1)))==[(1,1)]
    occ.synchronize()
    m=_current_model()
    # The edge runs start→end; center and major are not edge vertices.
    @test m.curves[a]==(s,e)
    @test m.curve_types[a]==:ellipse
    g=m.curve_geometry[a]
    @test g.occ==:ellipse && g.rx==4.0 && g.ry==2.0
    @test g.t1-g.t0≈π/2 atol=1e-12
    @test collect(M._occ_ellipse_point(g,g.t0))≈[4.0,0.0,0.0] atol=1e-12
    @test collect(M._occ_ellipse_point(g,g.t1))≈[0.0,2.0,0.0] atol=1e-12
    # The reversed traversal (start/end swapped) takes the same quarter —
    # upstream keeps the <= pi branch, re-framed so the range mirrors.
    b=occ.add_ellipse_arc(e,c,mj,s)
    occ.synchronize()
    g2=m.curve_geometry[b]
    @test g2.t1-g2.t0≈π/2 atol=1e-12
    @test collect(M._occ_ellipse_point(g2,g2.t0))≈[0.0,2.0,0.0] atol=1e-12
    @test collect(M._occ_ellipse_point(g2,g2.t1))≈[4.0,0.0,0.0] atol=1e-12
    @test collect(M._occ_ellipse_point(g2,(g2.t0+g2.t1)/2))≈
        collect(M._occ_ellipse_point(g,(g.t0+g.t1)/2)) atol=1e-12
    # Axis swap: the supplied major direction solves to the smaller radius;
    # OCCT warns and re-derives the frame.
    s2=occ.add_point(0.433,0.5,0.0)
    mj2=occ.add_point(1.0,0.0,0.0)
    e2=occ.add_point(-0.25,0.9682458366,0.0)
    sw=occ.add_ellipse_arc(s2,c,mj2,e2)
    occ.synchronize()
    g3=m.curve_geometry[sw]
    @test g3.rx>g3.ry
    @test collect(M._occ_ellipse_point(g3,g3.t0))≈[0.433,0.5,0.0] atol=1e-9
    @test collect(M._occ_ellipse_point(g3,g3.t1))≈[-0.25,0.9682458366,0.0] atol=1e-9
    # Every sampled point satisfies the conic in the stored frame.
    for t in range(g3.t0,g3.t1,length=9)
        p=collect(M._occ_ellipse_point(g3,t))
        ux=dot(p.-collect(g3.center),collect(g3.X))/g3.rx
        uy=dot(p.-collect(g3.center),collect(g3.Y))/g3.ry
        @test ux*ux+uy*uy≈1.0 atol=1e-12
    end
end

@testset "OCC add_ellipse_arc guards" begin
    occ=_occ_setup()
    s=occ.add_point(4.0,0.0,0.0)
    c=occ.add_point(0.0,0.0,0.0)
    mj=occ.add_point(4.0,0.0,0.0)
    e=occ.add_point(0.0,2.0,0.0)
    a=occ.add_ellipse_arc(s,c,mj,e)
    @test occursin("already exists",_errtext(()->occ.add_ellipse_arc(s,c,mj,e;tag=a)))
    @test occursin("unknown OCC Point",_errtext(()->occ.add_ellipse_arc(999,c,mj,e)))
    @test occursin("unknown OCC Point",_errtext(()->occ.add_ellipse_arc(s,c,mj,999)))
    # Both endpoints on the major axis.
    w=occ.add_point(-4.0,0.0,0.0)
    @test occursin("on the major axis",_errtext(()->occ.add_ellipse_arc(w,c,mj,s)))
    # Symmetric with respect to the major axis (same x, mirrored y).
    q1=occ.add_point(3.0,1.0,0.0); q2=occ.add_point(3.0,-1.0,0.0)
    @test occursin("symmetric",_errtext(()->occ.add_ellipse_arc(q1,c,mj,q2)))
    # Coincident center and major points -> gp_Dir on the null vector.
    @test occursin("null modulus",_errtext(()->occ.add_ellipse_arc(s,c,c,e)))
    # Failure atomicity: every rejection above left the registry untouched.
    @test sort!(collect(occ.get_entities(1)))==[(1,1)]
    @test sort!(collect(occ.get_entities(0)))==[(0,i) for i in 1:7]
end

@testset "OCC add_disk" begin
    occ=_occ_setup()
    d=occ.add_disk(1.0,2.0,0.0,4.0,2.0)
    @test d==1
    occ.synchronize()
    m=_current_model()
    @test m.surface_types[d]==:plane
    sg=m.surface_geometry[d]
    @test sg.axis==(0.0,0.0,1.0) && sg.center==(1.0,2.0,0.0)
    @test sg.uvb==([-4.0,-2.0],[4.0,2.0])
    # Topology: one loop of one ellipse edge, closed on its seam vertex.
    @test length(m.surfaces[d])==1
    edge=abs(only(m.loops[only(m.surfaces[d])]))
    @test m.curve_types[edge]==:ellipse
    va,vb=m.curves[edge]
    @test va==vb
    # The stored pcurve round-trips: pcurve(t) -> surface eval -> 3-D edge point.
    eg=m.curve_geometry[edge]
    pc=sg.pcurves[edge].fwd
    @test hasproperty(pc,:ellips2d)
    for t in range(eg.t0,eg.t1,length=13)
        p3=M._occ_ellipse_point(eg,t)
        uv=M._occ_pcurve_d0(pc,t)
        back=M._occ_surface_point(sg,uv[1],uv[2])
        @test collect(back)≈collect(p3) atol=1e-12
    end
    # Face classification: interior, near-boundary inside, outside.
    @test M._occ_surface_contains(m,d,sg,(1.0,2.0,0.0))
    @test M._occ_surface_contains(m,d,sg,(4.9,2.0,0.0))
    @test !M._occ_surface_contains(m,d,sg,(1.0,4.5,0.0))
    @test !M._occ_surface_contains(m,d,sg,(5.01,2.0,0.0))
    # Axis forms.
    d2=occ.add_disk(0.0,0.0,0.0,3.0,1.5;zAxis=[1.0,0.0,0.0],xAxis=[0.0,1.0,0.0])
    occ.synchronize()
    sg2=m.surface_geometry[d2]
    @test sg2.axis==(1.0,0.0,0.0)
    eg2=m.curve_geometry[abs(only(m.loops[only(m.surfaces[d2])]))]
    @test collect(M._occ_ellipse_point(eg2,0.0))≈[0.0,3.0,0.0] atol=1e-12
    @test collect(M._occ_ellipse_point(eg2,π/2))≈[0.0,0.0,1.5] atol=1e-12
    # Guards.
    @test occursin("positive",_errtext(()->occ.add_disk(0,0,0,0.0,2.0)))
    @test occursin("smaller",_errtext(()->occ.add_disk(0,0,0,4.0,5.0)))
    @test occursin("already exists",_errtext(()->occ.add_disk(0,0,0,4.0,2.0;tag=d)))
end

@testset "OCC ellipse transform and hole" begin
    occ=_occ_setup()
    d=occ.add_disk(0.0,0.0,0.0,4.0,2.0)
    occ.synchronize()
    m=_current_model()
    # Anisotropic dilation stretches each semi-axis independently.
    M.transform_entities!(m,M._affine_dilation((0.0,0.0,0.0),(2.0,3.0,1.0),
        "test"),[(2,d)];caller="test")
    eg=m.curve_geometry[abs(only(m.loops[only(m.surfaces[d])]))]
    @test eg.rx≈8.0 && eg.ry≈6.0
    @test collect(M._occ_ellipse_point(eg,0.0))≈[8.0,0.0,0.0] atol=1e-12
    @test collect(M._occ_ellipse_point(eg,π/2))≈[0.0,6.0,0.0] atol=1e-12
    # A shear would warp the conic: non-representable transforms reject.
    shear=M._AffineTransform((1.0,0.5,0.0,0.0,1.0,0.0,0.0,0.0,1.0),
        (0.0,0.0,0.0),NTuple{16,Float64}[],(linear=(1.0,0.5,0.0,0.0,1.0,0.0,0.0,0.0,1.0),loc=(0.0,0.0,0.0)))
    @test occursin("not representable",
        _errtext(()->M.transform_entities!(m,shear,[(2,d)];caller="test")))
    @test m.curve_geometry[abs(only(m.loops[only(m.surfaces[d])]))].rx≈8.0
    # Ellipse as a hole loop: the inner wire reverses on sync.
    occ=_occ_setup()
    p=[occ.add_point(x,y,0.0) for (x,y) in ((-10.0,-10.0),(10.0,-10.0),(10.0,10.0),(-10.0,10.0))]
    ls=[occ.add_line(p[1],p[2]),occ.add_line(p[2],p[3]),
        occ.add_line(p[3],p[4]),occ.add_line(p[4],p[1])]
    wo=occ.add_wire(ls)
    ei=occ.add_ellipse(0.0,0.0,0.0,4.0,2.0)
    wi=occ.add_wire([ei])
    f=occ.add_plane_surface([wo,wi])
    occ.synchronize()
    m=_current_model()
    sg=m.surface_geometry[f]
    # Inside the ellipse hole -> outside the face; between -> inside.
    @test M._occ_surface_contains(m,f,sg,(9.0,0.0,0.0))
    @test !M._occ_surface_contains(m,f,sg,(2.0,0.0,0.0))
    @test !M._occ_surface_contains(m,f,sg,(11.0,0.0,0.0))
end

@testset "OCC ellipse projection and bbox" begin
    occ=_occ_setup()
    e=occ.add_ellipse(0.0,0.0,0.0,4.0,2.0)
    t=occ.add_ellipse(0.0,0.0,0.0,4.0,2.0;angle1=0.0,angle2=π/2)
    occ.synchronize()
    m=_current_model()
    # Closest-point projection lands on the curve and is orthogonal.
    par=A.model.get_parametrization(1,e,[5.0,0.0,0.0])
    @test A.model.get_value(1,e,par)≈[4.0,0.0,0.0] atol=1e-9
    par=A.model.get_parametrization(1,e,[0.0,3.0,0.0])
    @test A.model.get_value(1,e,par)≈[0.0,2.0,0.0] atol=1e-9
    # Interior query: the foot stays on the ellipse.
    par=A.model.get_parametrization(1,e,[3.9,1.9,0.0])
    p=A.model.get_value(1,e,par)
    @test (p[1]/4)^2+(p[2]/2)^2≈1.0 atol=1e-9
    # Trimmed arc: an outside-range query clamps to the nearer endpoint —
    # (-4,0.5) is closest to the (0,2) terminus of the [0,pi/2] arc.
    par=A.model.get_parametrization(1,t,[-4.0,0.5,0.0])
    p=A.model.get_value(1,t,par)
    @test p≈[0.0,2.0,0.0] atol=1e-9
    # Trimmed bounding box tightens to the arc.
    @test collect(A.model.get_bounding_box(1,t))≈[0.0,0.0,0.0,4.0,2.0,0.0] atol=1e-9
end

@testset "OCC ellipse coherence merge" begin
    occ=_occ_setup()
    # Identical arcs on the same conic merge under explicit coherence.
    occ.add_ellipse(0.0,0.0,0.0,2.0,1.0;angle1=0.0,angle2=π/2)
    occ.add_ellipse(0.0,0.0,0.0,2.0,1.0;angle1=0.0,angle2=π/2)
    occ.synchronize()
    m=_current_model()
    M.coherence!(m)
    @test length(A.model.get_entities(1))==1
    # Complementary arcs share both endpoints but differ in span.
    occ=_occ_setup()
    occ.add_ellipse(0.0,0.0,0.0,2.0,1.0;angle1=0.0,angle2=π/2)
    occ.add_ellipse(0.0,0.0,0.0,2.0,1.0;angle1=π/2,angle2=2π)
    occ.synchronize()
    m=_current_model()
    M.coherence!(m)
    @test length(A.model.get_entities(1))==2
    # Same center and radii with a rotated major axis is a different edge.
    occ=_occ_setup()
    occ.add_ellipse(0.0,0.0,0.0,2.0,1.0;angle1=0.0,angle2=π/2)
    occ.add_ellipse(0.0,0.0,0.0,2.0,1.0;angle1=0.0,angle2=π/2,
                    zAxis=[0.0,0.0,1.0],xAxis=[0.8660254037844386,0.5,0.0])
    occ.synchronize()
    m=_current_model()
    M.coherence!(m)
    @test length(A.model.get_entities(1))==2
end

@testset "OCC ellipse meshing" begin
    occ=_occ_setup()
    d=occ.add_disk(0.0,0.0,0.0,4.0,2.0)
    occ.synchronize()
    A.mesh.generate(2)
    nt,xyz=A.mesh.get_nodes(2,d,true)
    @test length(nt)>4
    for i in 1:3:3length(nt)
        x,y=xyz[i],xyz[i+1]
        @test x*x/16+y*y/4<=1.0+1e-6
    end
end
