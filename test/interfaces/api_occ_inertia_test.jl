# `gmsh.model.occ.get_matrix_of_inertia` — upstream
# `OCC_Internals::getMatrixOfInertia` reports `GProps::MatrixOfInertia`
# under the per-dim Linear/Surface/VolumeProperties switch as nine
# row-major values: the inertia tensor I = tr(Q)I₃ − Q about the world
# origin with Q = ∫r⊗r dμ. Curves and faces integrate r⊗r per record —
# line/circle closed forms, the plane fan tensor, and Gauss-Legendre
# panel quadrature for the conic/spline/NURBS/quadric families; solids
# reuse the compact `primitive` records and the divergence theorem
# ∭x_i x_j dV = n_k∬x_i x_j x_k dA/5 with an exact cubic Stroud fan over
# all-plane shells. Conic wires, quadric-bounded shells and other
# unsupported boundaries reject loudly; dims outside 1:3 report the
# empty-GProp zero matrix like upstream's silent switch fall-through.
using Test

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
    # Segment on the x-axis, t∈[1,3]: Qxx=∫x²dx=26/3 → Iyy=Izz=26/3.
    a=occ.add_point(1.0,0.0,0.0); b=occ.add_point(3.0,0.0,0.0)
    l=occ.add_line(a,b)
    @test occ.get_matrix_of_inertia(1,l)≈
        [0.0,0.0,0.0, 0.0,26/3,0.0, 0.0,0.0,26/3]
    # Quarter arc r=1 on [0,π/2]: Qxx=Qyy=∫cos²=π/4, Qxy=∫sin·cos=1/2.
    c0=occ.add_point(0.0,0.0,0.0)
    p1=occ.add_point(1.0,0.0,0.0); p2=occ.add_point(0.0,1.0,0.0)
    arc=occ.add_circle_arc(p1,c0,p2)
    @test occ.get_matrix_of_inertia(1,arc)≈
        [π/4,-0.5,0.0, -0.5,π/4,0.0, 0.0,0.0,π/2] atol=1e-9
    # Full circle r=2: Izz=2πr³=16π, Ixx=Iyy=8π.
    ci=occ.add_circle(0.0,0.0,0.0,2.0)
    @test occ.get_matrix_of_inertia(1,ci)≈
        [8π,0.0,0.0, 0.0,8π,0.0, 0.0,0.0,16π] atol=1e-9
    # Full ellipse a=2,b=1: GL8 panels match a dense independent
    # quadrature to ~1e-6.
    e=occ.add_ellipse(0.0,0.0,0.0,2.0,1.0)
    Ie=occ.get_matrix_of_inertia(1,e)
    @test Ie[5]≈16.2795654 atol=1e-6   # Iyy = ∮x²ds
    @test Ie[1]≈5.6185567 atol=1e-6    # Ixx = ∮y²ds
    @test Ie[9]≈21.8981221 atol=1e-6   # Izz = ∮(x²+y²)ds
    @test _symmetric(Ie)
    # A bspline's quadrature stays symmetric about the plane of the
    # control polygon (z=0 → Ixz=Iyz=0).
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
    # Unit square z=0: ∬x²=∬y²=1/3, ∬xy=1/4.
    p=[occ.add_point(x,y,0.0) for (x,y) in ((0,0),(1,0),(1,1),(0,1))]
    ls=[occ.add_line(p[i],p[mod1(i+1,4)]) for i in 1:4]
    w=occ.add_wire(ls)
    f=occ.add_plane_surface([w])
    @test occ.get_matrix_of_inertia(2,f)≈
        [1/3,-0.25,0.0, -0.25,1/3,0.0, 0.0,0.0,2/3] atol=1e-12
    # A rectangle off the axes mixes products: [1,3]×[2,6] at z=1.
    p2=[occ.add_point(x,y,1.0) for (x,y) in ((1,2),(3,2),(3,6),(1,6))]
    ls2=[occ.add_line(p2[i],p2[mod1(i+1,4)]) for i in 1:4]
    f2=occ.add_plane_surface([occ.add_wire(ls2)])
    # A=8; ∬x²=(9-1)/3·(6-2)=32/3·4/4? → ∫x²=26/3? compute: ∫₁³x²dx=26/3
    # times Δy=4 → 104/3; ∬y²=(216-8)/3·2=416/3; ∬z²=8; ∬xy=4·16/2·8/8:
    # ∫x=4,∫y=16 → Qxy=4·16=64; ∬xz=∫x·1·Δy=16; ∬yz=∫y·1·Δx=32.
    Q=[104/3,64.0,16.0, 64.0,416/3,32.0, 16.0,32.0,8.0]
    tr=(104/3+416/3+8.0)
    @test occ.get_matrix_of_inertia(2,f2)≈
        [tr-Q[1],-Q[2],-Q[3], -Q[4],tr-Q[5],-Q[6], -Q[7],-Q[8],tr-Q[9]]
    # Cylinder lateral face r=1 h=2: Izz=r²·A=4π; ∬y²=r²πh=2π;
    # ∬z²=(8/3)·2π=16π/3 → Ixx=2π+16π/3=22π/3.
    occ.add_cylinder(0.0,0.0,0.0,0.0,0.0,2.0,1.0)
    @test occ.get_matrix_of_inertia(2,3)≈
        [22π/3,0.0,0.0, 0.0,22π/3,0.0, 0.0,0.0,4π] atol=1e-8
    # Unit sphere shell: ∬xᵢ²dA=4π/3 → I=8π/3·I₃.
    occ.add_sphere(0.0,0.0,0.0,1.0)
    @test occ.get_matrix_of_inertia(2,6)≈
        [8π/3,0.0,0.0, 0.0,8π/3,0.0, 0.0,0.0,8π/3] atol=1e-7
    # Disk faces carry conic wires → the plane tensor path rejects.
    d=occ.add_disk(0.0,0.0,0.0,1.0,1.0)
    @test occursin("conic",_errtext(()->occ.get_matrix_of_inertia(2,d)))
end

@testset "occ.get_matrix_of_inertia — solids" begin
    occ=_occ_setup()
    # Unit box at origin: ∭x²=∭y²=∭z²=1/3, ∭xy=∭xz=∭yz=1/4.
    b=occ.add_box(0.0,0.0,0.0,1.0,1.0,1.0)
    @test occ.get_matrix_of_inertia(3,b)≈
        [2/3,-0.25,-0.25, -0.25,2/3,-0.25, -0.25,-0.25,2/3] atol=1e-12
    # Box [0,2]×[0,1]×[0,3]: V=6; ∭x²=8,∭y²=2,∭z²=18; ∭xy=3,∭xz=9,
    # ∭yz=4.5 — the divergence fan lands on the closed form.
    b2=occ.add_box(0.0,0.0,0.0,2.0,1.0,3.0)
    @test occ.get_matrix_of_inertia(3,b2)≈
        [20.0,-3.0,-9.0, -3.0,26.0,-4.5, -9.0,-4.5,10.0] atol=1e-12
    # Triangular wedge dx=4,dy=3,dz=2 (profile legs 4,2): ∭x²=32,
    # ∭y²=36, ∭z²=8, ∭xy=24, ∭xz=8, ∭yz=12 — non-rectangular fan.
    wg=occ.add_wedge(0.0,0.0,0.0,4.0,3.0,2.0;ltx=0.0)
    @test occ.get_matrix_of_inertia(3,wg)≈
        [44.0,-24.0,-8.0, -24.0,40.0,-12.0, -8.0,-12.0,68.0] atol=1e-12
    # Sphere r=1 at origin: I=(2/5)Vr²·I₃=8π/15·I₃.
    s=occ.add_sphere(0.0,0.0,0.0,1.0)
    @test occ.get_matrix_of_inertia(3,s)≈
        [8π/15,0.0,0.0, 0.0,8π/15,0.0, 0.0,0.0,8π/15] atol=1e-12
    # Cylinder r=2 h=3 on the base origin: Izz=Vr²/2=24π;
    # Ixx=V(r²/4+h²/3)=48π.
    cy=occ.add_cylinder(0.0,0.0,0.0,0.0,0.0,3.0,2.0)
    @test occ.get_matrix_of_inertia(3,cy)≈
        [48π,0.0,0.0, 0.0,48π,0.0, 0.0,0.0,24π] atol=1e-12
    # Apex cone r=2 h=3 about the base centre: Izz=3Vr²/10=24π/5,
    # Ixx=∭y²+∭z²=12π/5+18π/5=6π.
    co=occ.add_cone(0.0,0.0,0.0,0.0,0.0,3.0,2.0,0.0)
    @test occ.get_matrix_of_inertia(3,co)≈
        [6π,0.0,0.0, 0.0,6π,0.0, 0.0,0.0,24π/5] atol=1e-12
    # Torus-bounded shells are outside the all-plane/primitive support
    # and reject rather than approximate.
    t=occ.add_torus(0.0,0.0,0.0,3.0,1.0)
    @test occursin("unsupported solid boundary",
                   _errtext(()->occ.get_matrix_of_inertia(3,t)))
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
    # Kernel-side visibility across transforms, pre-synchronization:
    # unit box shifted by (10,20,30) → Q=J_c+V·c⊗c, J_c=diag(1/12).
    occ.translate([(3,b)],10.0,20.0,30.0)
    c=(10.5,20.5,30.5)
    Q=[(i==j ? 1/12 : 0.0)+c[i]*c[j] for i in 1:3,j in 1:3]
    trq=sum(Q[i,i] for i in 1:3)
    Ie=[i==j ? trq-Q[i,j] : -Q[i,j] for i in 1:3,j in 1:3]
    @test occ.get_matrix_of_inertia(3,b)≈vec(Ie') atol=1e-9
    occ.synchronize()
    @test occ.get_matrix_of_inertia(3,b)≈vec(Ie') atol=1e-9
end
