module OCCBooleanExactConicChecks
using Test,Tessella
const M=Tessella.Model

# A plane cuts a sphere in a circle, a cylinder obliquely in an ellipse, and
# two spheres meet in a circle — exact conics even when their trace on the
# curved face's chart is neither a parallel nor a meridian. These sections
# used to be traced numerically (1e-7): the edges were anonymous splines, the
# plane faces' areas were off by up to 1e-3, volumes by 1e-5, and the plane
# faces could not be meshed ("Point[N] is not coplanar with Surface[M]").
# The 3-D conic is kept exact now and only its chart trace is fitted.
#
# The section builders are called directly: a Boolean that walks a fitted
# trace through trimming and region building costs minutes interpreted (the
# three Booleans this file used to run took 2300 s). With the exact sections
# Box(-1..1)³ minus Sphere(.8,.55,.35,.6) is 7.3636776202, the closed form
# (Gmsh 4.15.2: 7.3636776470).
pcs(pc)=pc isa AbstractVector ? pc : [pc]
function lifts_onto(sg,pc,curve)
    worst=0.0
    for θ in range(0.1,6.2;length=23)
        uv=M._brep_pc_eval(pc,θ)
        worst=max(worst,M._b3dist(M._brep_surface_eval(sg,uv[1],uv[2]),
                                  M._brep_curve_eval(curve,θ)))
    end
    return worst
end

@testset "OCC Boolean sections that are exact conics" begin
    ball=(occ=:sphere,center=(0.0,0.0,0.0),radius=1.0,axis=(0.0,0.0,1.0),
          X=(1.0,0.0,0.0),Y=(0.0,1.0,0.0))
    # the cutting plane x=0.4 is off the centre: a circle of radius √0.84
    plane=(occ=:plane,center=(0.4,0.0,0.0),axis=(1.0,0.0,0.0),
           X=(0.0,1.0,0.0),Y=(0.0,0.0,1.0))
    section=only(M._brep_sec_plane_sphere(plane,ball,true,"test"))
    @test section.closed && section.curve.kind===:circle
    @test section.curve.r≈sqrt(0.84) rtol=1e-14
    @test all(pc->haskey(pc,:circ2d),pcs(section.pcA))
    trace=only(pcs(section.pcB))
    @test haskey(trace,:per2d)
    @test lifts_onto(ball,trace,section.curve)<2e-9

    # two unit spheres 1.2 apart meet in the circle x=0.6, radius 0.8
    other=merge(ball,(center=(1.2,0.0,0.0),))
    section=only(M._brep_sec_sphere_sphere(ball,other,"test"))
    @test section.curve.kind===:circle
    @test section.curve.r≈0.8 rtol=1e-14
    @test collect(section.curve.c)≈[0.6,0.0,0.0] atol=1e-14
    @test lifts_onto(other,only(pcs(section.pcB)),section.curve)<2e-9

    # a plane tilted 0.5 rad about y cuts a cylinder of radius 0.5 in an
    # ellipse with semi-axes 0.5/cos(0.5) and 0.5
    rod=(occ=:cylinder,center=(0.0,0.0,-1.0),axis=(0.0,0.0,1.0),
         X=(1.0,0.0,0.0),Y=(0.0,1.0,0.0),radius=0.5,height=2.0)
    tilted=(occ=:plane,center=(0.0,0.0,0.0),axis=(sin(0.5),0.0,cos(0.5)),
            X=(cos(0.5),0.0,-sin(0.5)),Y=(0.0,1.0,0.0))
    section=only(M._brep_sec_plane_cylinder(tilted,rod,true,"test"))
    @test section.curve.kind===:ellipse
    @test sort([section.curve.rx,section.curve.ry])≈[0.5,0.5/cos(0.5)] rtol=1e-14
    @test lifts_onto(rod,only(pcs(section.pcB)),section.curve)<2e-9
end
end
