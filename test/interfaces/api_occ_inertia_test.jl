# `gmsh.model.occ.get_matrix_of_inertia` — upstream
# `OCC_Internals::getMatrixOfInertia` reports the centroidal inertia tensor
# (Gmsh 4.15.2 public control: a unit box is diag(1/6), not the origin-moment
# diag(2/3) with -1/4 products) as nine row-major values. The tensor is
# I = tr(Q_c)I₃ − Q_c with Q_c = ∫(r-c)⊗(r-c)dμ about the entity centroid.
# Curves and faces integrate the centered tensor per record — line/circle
# closed forms, the plane fan tensor, and Gauss-Legendre panel quadrature for
# the conic/spline/NURBS/quadric families; solids reuse the compact
# `primitive` records and the divergence theorem
# ∭(r-c)_i(r-c)_j dV = n_k∬(r-c)_i(r-c)_j(r-c)_k dA/5. Conic/spline plane
# wires use boundary 1-forms; full/partial torus sectors use Pappus moments;
# trimmed quadrics use the centered CDT divergence. Unparameterized shells
# reject rather than approximating. Dims outside 1:3 report the empty-GProp
# zero matrix like upstream's silent switch fall-through.
using Test
using LinearAlgebra: Diagonal

using Tessella
const A=Tessella.API

function _occ_setup()
    A.initialize()
    return A.model.occ
end

function _errtext(f)
    try
        f()
    catch err
        return sprint(showerror,err)
    end
    return ""
end

_symmetric(v)=isapprox(v[2],v[4];atol=1e-12) &&
    isapprox(v[3],v[7];atol=1e-12) && isapprox(v[6],v[8];atol=1e-12)

@testset "occ.get_matrix_of_inertia — curves" begin
    occ=_occ_setup()
    # Segment on the x-axis t∈[1,3]: centroid (2,0,0); ∫(x-2)²dx=2/3.
    a=occ.add_point(1.0,0.0,0.0); b=occ.add_point(3.0,0.0,0.0)
    l=occ.add_line(a,b)
    @test occ.get_matrix_of_inertia(1,l)≈
        [0.0,0.0,0.0, 0.0,2/3,0.0, 0.0,0.0,2/3]
    # Quarter arc r=1: centroid (2/π,2/π,0). Centered Qxx=Qyy=π/4−2/π,
    # Qxy=1/2−2/π, Qzz=0.
    c0=occ.add_point(0.0,0.0,0.0)
    p1=occ.add_point(1.0,0.0,0.0); p2=occ.add_point(0.0,1.0,0.0)
    arc=occ.add_circle_arc(p1,c0,p2)
    @test occ.get_matrix_of_inertia(1,arc)≈
        [π/4-2/π,2/π-1/2,0.0, 2/π-1/2,π/4-2/π,0.0,
         0.0,0.0,π/2-4/π] atol=1e-9
    # Full circle r=2 centered at the centroid: unchanged from origin moments.
    ci=occ.add_circle(0.0,0.0,0.0,2.0)
    @test occ.get_matrix_of_inertia(1,ci)≈
        [8π,0.0,0.0, 0.0,8π,0.0, 0.0,0.0,16π] atol=1e-9
    # Full ellipse a=2,b=1 (centered): exact arc-length moments.
    e=occ.add_ellipse(0.0,0.0,0.0,2.0,1.0)
    Ie=occ.get_matrix_of_inertia(1,e)
    @test Ie[5]≈16.27956516493 atol=1e-6   # ∮x²ds
    @test Ie[1]≈5.61855692932 atol=1e-6    # ∮y²ds
    @test Ie[9]≈21.89812209425 atol=1e-6   # ∮(x²+y²)ds
    @test _symmetric(Ie)
    # A bspline's quadrature stays symmetric about the plane of the
    # control polygon (z=0 → Ixz=Iyz=0) and its centroid lies on x=3/2.
    pts=[occ.add_point(x,y,0.0) for (x,y) in ((0,0),(1,2),(2,2),(3,0))]
    bs=occ.add_bspline(pts)
    Ib=occ.get_matrix_of_inertia(1,bs)
    @test _symmetric(Ib)
    @test Ib[3]≈0.0 atol=1e-12
    @test Ib[6]≈0.0 atol=1e-12
    @test Ib[1]>0.0 && Ib[9]>Ib[1]
end

@testset "occ.get_matrix_of_inertia — faces" begin
    occ=_occ_setup()
    # Unit square z=0: centroid (1/2,1/2); centered Q=diag(1/12,1/12,0).
    p=[occ.add_point(x,y,0.0) for (x,y) in ((0,0),(1,0),(1,1),(0,1))]
    ls=[occ.add_line(p[i],p[mod1(i+1,4)]) for i in 1:4]
    w=occ.add_wire(ls)
    f=occ.add_plane_surface([w])
    @test occ.get_matrix_of_inertia(2,f)≈
        [1/12,0.0,0.0, 0.0,1/12,0.0, 0.0,0.0,1/6] atol=1e-12
    # Rectangle [1,3]×[2,6] at z=1: centroid (2,4,1);
    # centered Q=diag(8/3,32/3,0) → I=diag(32/3,8/3,40/3).
    p2=[occ.add_point(x,y,1.0) for (x,y) in ((1,2),(3,2),(3,6),(1,6))]
    ls2=[occ.add_line(p2[i],p2[mod1(i+1,4)]) for i in 1:4]
    f2=occ.add_plane_surface([occ.add_wire(ls2)])
    @test occ.get_matrix_of_inertia(2,f2)≈
        [32/3,0.0,0.0, 0.0,8/3,0.0, 0.0,0.0,40/3] atol=1e-12
    # Cylinder lateral face r=1 h=2: centroid (0,0,1);
    # centered Q=diag(2π,2π,4π/3) → I=diag(10π/3,10π/3,4π).
    occ.add_cylinder(0.0,0.0,0.0,0.0,0.0,2.0,1.0)
    lateral=only([tag for (dim,tag) in occ.get_entities(2)
        if isapprox(occ.get_mass(2,tag),4π;rtol=1e-10)])
    @test occ.get_matrix_of_inertia(2,lateral)≈
        [10π/3,0.0,0.0, 0.0,10π/3,0.0, 0.0,0.0,4π] atol=1e-8
    # Unit sphere shell: centered at the origin; Q=4π/3·I.
    occ.add_sphere(0.0,0.0,0.0,1.0)
    shell=maximum(tag for (dim,tag) in occ.get_entities(2))
    @test occ.get_matrix_of_inertia(2,shell)≈
        [8π/3,0.0,0.0, 0.0,8π/3,0.0, 0.0,0.0,8π/3] atol=1e-7
    # Disk face r=1 at z=0: conic boundary moments; centered
    # I=diag(π/4,π/4,π/2) (Gmsh 4.15.2 public control).
    d=occ.add_disk(0.0,0.0,0.0,1.0,1.0)
    @test collect(occ.get_center_of_mass(2,d))≈[0.,0.,0.] atol=1e-12
    @test occ.get_matrix_of_inertia(2,d)≈
        [π/4,0.0,0.0, 0.0,π/4,0.0, 0.0,0.0,π/2] atol=1e-9
end

@testset "occ.get_matrix_of_inertia — solids" begin
    occ=_occ_setup()
    # Unit box at origin: centroid (1/2,1/2,1/2); centered Q=(1/12)I.
    b=occ.add_box(0.0,0.0,0.0,1.0,1.0,1.0)
    @test occ.get_matrix_of_inertia(3,b)≈
        [1/6,0.0,0.0, 0.0,1/6,0.0, 0.0,0.0,1/6] atol=1e-12
    # Box [0,2]×[0,1]×[0,3]: centroid (1,1/2,3/2);
    # centered Q=diag(2,1/2,9/2) → I=diag(5,13/2,5/2).
    b2=occ.add_box(0.0,0.0,0.0,2.0,1.0,3.0)
    @test occ.get_matrix_of_inertia(3,b2)≈
        [5.0,0.0,0.0, 0.0,6.5,0.0, 0.0,0.0,2.5] atol=1e-12
    # Gmsh/OpenCASCADE wedge dx=4,dy=3,dz=2: triangular xy profile
    # extruded along z; centroid (4/3,1,1). Public Gmsh 4.15.2 control:
    # I=[10,4,0; 4,44/3,0; 0,0,50/3].
    wg=occ.add_wedge(0.0,0.0,0.0,4.0,3.0,2.0;ltx=0.0)
    @test occ.get_matrix_of_inertia(3,wg)≈
        [10.0,4.0,0.0, 4.0,44/3,0.0, 0.0,0.0,50/3] atol=1e-12
    # Sphere r=1 at origin: I=(2/5)Vr²·I₃=8π/15·I₃.
    s=occ.add_sphere(0.0,0.0,0.0,1.0)
    @test occ.get_matrix_of_inertia(3,s)≈
        [8π/15,0.0,0.0, 0.0,8π/15,0.0, 0.0,0.0,8π/15] atol=1e-12
    # Cylinder r=2 h=3 on the base origin: centroid (0,0,3/2);
    # centered Q=diag(12π,12π,9π) → I=diag(21π,21π,24π).
    cy=occ.add_cylinder(0.0,0.0,0.0,0.0,0.0,3.0,2.0)
    @test occ.get_matrix_of_inertia(3,cy)≈
        [21π,0.0,0.0, 0.0,21π,0.0, 0.0,0.0,24π] atol=1e-12
    # Apex cone r=2 h=3 about the base centre: centered
    # Ixx=Iyy=15π/4 and Izz=24π/5 (Gmsh 4.15.2 control).
    co=occ.add_cone(0.0,0.0,0.0,0.0,0.0,3.0,2.0,0.0)
    @test occ.get_matrix_of_inertia(3,co)≈
        [15π/4,0.0,0.0, 0.0,15π/4,0.0, 0.0,0.0,24π/5] atol=1e-12
    # Axisymmetric/spherical sector primitives retain their analytic
    # centroidal moments (Gmsh 4.15.2 controls for π/2 sectors).
    cs=occ.add_cylinder(0.0,0.0,0.0,0.0,0.0,3.0,2.0;angle=π/2)
    @test collect(occ.get_center_of_mass(3,cs))≈
        collect((0.848826363156775,0.848826363156775,1.5)) atol=1e-10
    @test occ.get_matrix_of_inertia(3,cs)≈
        [9.702750526092,0.790610905254,0.0, 0.790610905254,9.702750526092,0.0,
         0.0,0.0,5.26833411103] rtol=1e-10
    cos=occ.add_cone(0.0,0.0,0.0,0.0,0.0,3.0,2.0,0.5;angle=π/2)
    @test occ.get_matrix_of_inertia(3,cos)≈
        [2.973887818387,0.112714975951,0.461383928571,
         0.112714975951,2.973887818387,0.461383928571,
         0.461383928571,0.461383928571,1.599334555321] rtol=1e-10
    sps=occ.add_sphere(0.0,0.0,0.0,2.0;angle1=-π/4,angle2=π/4,angle3=π/2)
    @test occ.get_matrix_of_inertia(3,sps)≈
        [5.958878968496,0.501865106343,0.0, 0.501865106343,5.958878968496,0.0,
         0.0,0.0,3.62437645243] rtol=1e-10
    # General trimmed quadric shell: sphere cut by z>=0 (lower hemisphere).
    hs=occ.add_sphere(0.0,0.0,0.0,1.0;tag=50)
    hb=occ.add_box(-2.0,-2.0,0.0,4.0,4.0,2.0;tag=51)
    hout,_=occ.cut([(3,50)],[(3,51)];tag=52)
    @test hout==[(3,52)]
    occ.synchronize()
    @test occ.get_mass(3,52)≈2π/3 rtol=1e-6
    @test collect(occ.get_center_of_mass(3,52))≈[0.,0.,-0.375] atol=1e-6
    @test occ.get_matrix_of_inertia(3,52)≈
        [0.543233729683,0.0,0.0, 0.0,0.543233729683,0.0,
         0.0,0.0,0.837758040957] rtol=1e-6 atol=1e-7
    t=occ.add_torus(0.0,0.0,0.0,3.0,1.0)
    @test occ.get_matrix_of_inertia(3,t)≈
        [303.490335333498,0.0,0.0, 0.0,303.490335333498,0.0,
         0.0,0.0,577.371857463727] rtol=1e-10
    tp=occ.add_torus(0.0,0.0,0.0,3.0,1.0;angle=π/2)
    @test collect(occ.get_center_of_mass(3,tp))≈
        [1.9629109648000427,1.962910964800042,0.0] atol=1e-10
    @test occ.get_matrix_of_inertia(3,tp)≈
        [18.830917166708,11.095874107916,0.0,
         11.095874107916,18.830917166708,0.0,
         0.0,0.0,30.259631032599] rtol=1e-8
end

@testset "occ.get_matrix_of_inertia — misc dims + lifecycle" begin
    occ=_occ_setup()
    z=zeros(9)
    # Points, wires, shells and an empty volume report the zero matrix
    # upstream's unfired switch leaves behind.
    p=occ.add_point(1.0,2.0,3.0)
    @test occ.get_matrix_of_inertia(0,p)≈z
    p2=occ.add_point(2.0,0.0,0.0)
    l=occ.add_line(p,p2); w=occ.add_wire([l])
    @test occ.get_matrix_of_inertia(-1,w)≈z
    b=occ.add_box(0.0,0.0,0.0,1.0,1.0,1.0)
    @test occ.get_matrix_of_inertia(-2,1)≈z
    v=occ.add_volume(Int[])
    @test occ.get_matrix_of_inertia(3,v)≈z
    # Degenerate pole edges are zero-length → zero matrix.
    occ.add_sphere(0.0,0.0,0.0,2.0)
    m=A._with_model(identity)
    deg=sort([t for ((d,t),id) in m.occ_kernel.bindings
              if d==1 && m.occ_kernel.shapes[id].geometry.geometry.occ===
                  :degenerate])
    @test !isempty(deg)
    @test occ.get_matrix_of_inertia(1,deg[1])≈z
    # Unknown entity and out-of-range dims reject.
    @test occursin("Unknown OpenCASCADE entity of dimension 3 with tag 99",
                   _errtext(()->occ.get_matrix_of_inertia(3,99)))
    @test occursin("dimension",_errtext(()->occ.get_matrix_of_inertia(4,1)))
    # Centroidal moments are translation invariant: a unit box shifted by
    # (10,20,30) keeps diag(1/6).
    shifted=occ.add_box(0.0,0.0,0.0,1.0,1.0,1.0)
    occ.translate([(3,shifted)],10.0,20.0,30.0)
    Ie=vec(Diagonal(fill(1/6,3)))
    @test occ.get_matrix_of_inertia(3,shifted)≈Ie atol=1e-9
    occ.synchronize()
    @test occ.get_matrix_of_inertia(3,shifted)≈Ie atol=1e-9
end
