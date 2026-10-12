# `model.occ` circle endpoints: `add_circle_arc` (center and through-point
# forms) and `add_circle` (center/radius/axis form, full and trimmed).
# Shapes bind into the OCC kernel with copied geometry, `synchronize`
# imports them as :circle GModel curves, and the OCC geometry records
# drive evaluation and 1-D meshing — matching `OCC_Internals::addCircleArc`
# / `addCircle` semantics observed on Gmsh 4.15.2.

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

@testset "OCC add_circle_arc center form" begin
    occ=_occ_setup()
    s=occ.add_point(1.0,0.0,0.0)
    c=occ.add_point(0.0,0.0,0.0)
    e=occ.add_point(0.0,1.0,0.0)
    @test (s,c,e)==(1,2,3)
    @test occ.add_circle_arc(s,c,e)==1
    @test sort!(collect(occ.get_entities(1)))==[(1,1)]
    occ.synchronize()
    m=_current_model()
    # The edge runs start→end; the center point is not an edge vertex.
    @test m.curves[1]==(s,e)
    @test m.curve_types[1]==:circle
    lower,upper=A.model.get_parametrization_bounds(1,1)
    @test lower==[0.0] || lower[1]>=0.0
    # Whatever frame OCCT picks, the trim evaluates onto the unit circle
    # and lands on both vertex positions at the bounds.
    @test A.model.get_value(1,1,[lower[1]])≈collect(m.points[s]) atol=1e-12
    @test A.model.get_value(1,1,[upper[1]])≈collect(m.points[e]) atol=1e-12
    mid=A.model.get_value(1,1,[(lower[1]+upper[1])/2])
    @test sqrt(mid[1]^2+mid[2]^2+mid[3]^2)≈1.0 atol=1e-12
    # Explicit tag honored; duplicate rejected.
    @test occ.add_circle_arc(s,c,e;tag=7)==7
    @test occursin("already exists",_errtext(()->occ.add_circle_arc(s,c,e;tag=7)))
end

@testset "OCC add_circle_arc through-point form" begin
    occ=_occ_setup()
    p1=occ.add_point(1.0,0.0,0.0)
    p2=occ.add_point(0.0,-1.0,0.0)   # middle point ON the arc, long way
    p3=occ.add_point(0.0,1.0,0.0)
    @test occ.add_circle_arc(p1,p2,p3;center=false)==1
    occ.synchronize()
    m=_current_model()
    @test m.curves[1]==(p1,p3)
    @test m.curve_types[1]==:circle
    # Trim must span 3π/2 — the arc containing the middle point.
    lower,upper=A.model.get_parametrization_bounds(1,1)
    @test upper[1]-lower[1]≈3π/2 atol=1e-12
    # Every evaluated point stays on the unit circle and the middle
    # control point lies on the trim (its parameter is interior).
    for u in range(lower[1],upper[1];length=7)
        p=A.model.get_value(1,1,[u])
        @test sqrt(p[1]^2+p[2]^2+p[3]^2)≈1.0 atol=1e-12
    end
    middle_param=A.model.get_parametrization(1,1,[0.0,-1.0,0.0])
    @test lower[1]<middle_param[1]<upper[1]
    @test A.model.get_value(1,1,middle_param)≈[0.0,-1.0,0.0] atol=1e-12
end

@testset "OCC add_circle_arc guards" begin
    occ=_occ_setup()
    s=occ.add_point(1.0,0.0,0.0)
    c=occ.add_point(0.0,0.0,0.0)
    e=occ.add_point(0.0,1.0,0.0)
    @test occursin("unknown OCC Point[99]",_errtext(()->occ.add_circle_arc(99,c,e)))
    @test occursin("center must be a boolean",_errtext(()->occ.add_circle_arc(s,c,e;center=1)))
    far=occ.add_point(0.0,2.0,0.0)
    @test occursin("Could not create circle arc",_errtext(()->occ.add_circle_arc(s,c,far)))
    a=occ.add_point(2.0,0.0,0.0)
    b=occ.add_point(3.0,0.0,0.0)
    @test occursin("Could not build circle through three points",
                   _errtext(()->occ.add_circle_arc(s,a,b;center=false)))
    # Failures bind nothing: the curve list is still empty.
    @test occ.get_entities(1)==[]
    @test occursin("U1 == U2",_errtext(()->occ.add_circle_arc(s,c,s)))
end

@testset "OCC add_circle full and trimmed" begin
    occ=_occ_setup()
    @test occ.add_circle(0.0,0.0,0.0,2.0)==1
    @test occ.add_circle(0.0,0.0,0.0,1.0;tag=9,angle1=0.0,angle2=π/2)==9
    occ.synchronize()
    m=_current_model()
    # Full circle: single bound vertex at center+r·X, edge closed on it.
    v=m.curves[1][1]
    @test m.curves[1]==(v,v)
    @test collect(m.points[v])≈[2.0,0.0,0.0] atol=1e-12
    @test sort!(collect(occ.get_entities(1)))==[(1,1),(1,9)]
    lower,upper=A.model.get_parametrization_bounds(1,1)
    @test lower==[0.0] && upper[1]≈2π
    for u in (0.0,π/3,π,4.7)
        p=A.model.get_value(1,1,[u])
        @test sqrt(p[1]^2+p[2]^2+p[3]^2)≈2.0 atol=1e-12
    end
    # Trimmed arc: two distinct bound vertices at the angle endpoints.
    v1,v2=m.curves[9]
    @test v1!=v2
    @test collect(m.points[v1])≈[1.0,0.0,0.0] atol=1e-12
    @test collect(m.points[v2])≈[0.0,1.0,0.0] atol=1e-12
    @test A.model.get_parametrization_bounds(1,9)==([0.0],[π/2])
end

@testset "OCC add_circle axis forms" begin
    occ=_occ_setup()
    # Reversed z-axis: trim runs clockwise seen from +z.
    @test occ.add_circle(0,0,0,1;zAxis=[0,0,-1])==1
    # Custom plane + reference direction.
    @test occ.add_circle(0,0,0,1;zAxis=[0,1,0],xAxis=[1,0,0])==2
    # xAxis without zAxis is ignored upstream — default z/x frame.
    @test occ.add_circle(0,0,0,1;xAxis=[0,1,0])==3
    occ.synchronize()
    m=_current_model()
    g1,g2,g3=m.curve_geometry[1],m.curve_geometry[2],m.curve_geometry[3]
    @test g1.n==(0.0,0.0,-1.0)
    @test g2.n==(0.0,1.0,0.0)
    @test collect(g2.X)≈[1.0,0.0,0.0] atol=1e-12
    @test g3.n==(0.0,0.0,1.0)
    @test collect(g3.X)≈[1.0,0.0,0.0] atol=1e-12
    # Circle on the y-normal plane has y=0 everywhere.
    @test A.model.get_value(1,2,[π/2])[2]≈0.0 atol=1e-12
    # Reversed axis arc still spans π/2.
    a=occ.add_circle(0,0,0,1;angle1=0.0,angle2=π/2,zAxis=[0,0,-1])
    occ.synchronize()
    @test A.model.get_parametrization_bounds(1,a)==([0.0],[π/2])
end

@testset "OCC add_circle guards" begin
    occ=_occ_setup()
    @test occursin("radius should be positive",_errtext(()->occ.add_circle(0,0,0,0.0)))
    @test occursin("U1 == U2",_errtext(()->occ.add_circle(0,0,0,1;angle1=1.0,angle2=1.0)))
    @test occursin("gp_Dir",_errtext(()->occ.add_circle(0,0,0,1;zAxis=[0,0,0])))
    @test occursin("parallel",_errtext(()->occ.add_circle(0,0,0,1;zAxis=[0,0,1],xAxis=[0,0,1])))
    @test occursin("3-element vector",_errtext(()->occ.add_circle(0,0,0,1;zAxis=[0,1])))
    @test occursin("3-element vector",_errtext(()->occ.add_circle(0,0,0,1;zAxis="z")))
    @test occursin("must be finite",_errtext(()->occ.add_circle(0,0,0,1;angle1=NaN)))
    @test occ.get_entities(1)==[]
    # Angle swap: [max,min] input trims CCW min→max like upstream's
    # Geom_TrimmedCurve parameter swap.
    @test occ.add_circle(0,0,0,1;angle1=3π/2,angle2=π/2)==1
    occ.synchronize()
    @test A.model.get_parametrization_bounds(1,1)==([π/2],[3π/2])
end

@testset "OCC circle end-to-end mesh" begin
    occ=_occ_setup()
    s=occ.add_point(1.0,0.0,0.0)
    c=occ.add_point(0.0,0.0,0.0)
    e=occ.add_point(0.0,1.0,0.0)
    occ.add_circle_arc(s,c,e)
    occ.add_circle(5.0,0.0,0.0,2.0)
    occ.synchronize()
    A.mesh.generate(1)
    _,coords,_=A.mesh.get_nodes()
    pts=[[coords[i],coords[i+1],coords[i+2]] for i in 1:3:length(coords)]
    on_arc=count(p->abs(sqrt(p[1]^2+p[2]^2+p[3]^2)-1.0)<1e-9,pts)
    on_full=count(p->abs(sqrt((p[1]-5.0)^2+p[2]^2+p[3]^2)-2.0)<1e-9,pts)
    # The arc center stays a standalone vertex entity and meshes to a node.
    on_center=count(p->sqrt(p[1]^2+p[2]^2+p[3]^2)<1e-12,pts)
    @test on_arc>=4
    @test on_full>=4
    @test on_center==1
    @test on_arc+on_full+on_center==length(pts)
end
