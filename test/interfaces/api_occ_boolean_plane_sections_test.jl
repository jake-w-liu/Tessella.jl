module OCCBooleanPlaneSectionsChecks
using Test,Tessella
const A=Tessella.API
const M=Tessella.Model

# A plane meets a cone in a conic and a torus, when it is parallel to the
# axis, in a spiric section. Both were followed by the level-set tracer on
# the curved face's chart to 1e-7, with the plane's pcurve fitted after it:
# the edge left the plane, so the plane faces could not be meshed ("Point[21]
# is not coplanar with Surface[13]"), `fuse` of a box and a cone threw
# "degenerate zero-area wire", and volumes were 5e-6 off. The sections now
# come from their closed-form parametrization: both traces lift onto the
# exact curve within 1e-10 and the curve is the plane's own lift.
function worst_lift(section,plane,curved,npoints)
    lo,hi=M._brep_pc_domain(section.pcA)
    off_plane=0.0;off_curved=0.0;apart=0.0
    for i in 0:npoints
        t=lo+(hi-lo)*i/npoints
        P=M._brep_curve_eval(section.curve,t)
        off_plane=max(off_plane,abs(M._b3dot(M._b3(P,plane.center),plane.axis)))
        uv=M._brep_pc_eval(section.pcB,t)
        Q=M._brep_surface_eval(curved,uv[1],uv[2])
        apart=max(apart,M._b3dist(P,Q))
        st=M._brep_pc_eval(section.pcA,t)
        off_curved=max(off_curved,M._b3dist(P,M._brep_surface_eval(plane,st[1],st[2])))
    end
    return off_plane,apart,off_curved
end

@testset "OCC plane sections of cones and tori" begin
    # x = 1 through a cone with its axis along z, 0.2 from the axis where
    # the cone's top radius is 0.24: the two arms of a hyperbola, each
    # running from one rim to the other
    cone=(occ=:cone,center=(0.8,0.55,-0.25),axis=(0.0,0.0,1.0),
          X=(1.0,0.0,0.0),Y=(0.0,1.0,0.0),r1=0.6,r2=0.24,height=1.2)
    wall=(occ=:plane,center=(1.0,0.0,0.0),axis=(1.0,0.0,0.0),
          X=(0.0,1.0,0.0),Y=(0.0,0.0,1.0))
    arms=M._brep_sec_cone_conic(wall,cone,true,"test")
    @test length(arms)==2 && !any(s->s.closed,arms)
    slant=hypot(1.2,0.36)
    for arm in arms
        off_plane,apart,_=worst_lift(arm,wall,cone,57)
        @test off_plane<1e-14
        @test apart<1e-9
        lo,hi=M._brep_pc_domain(arm.pcB)
        @test sort([M._brep_pc_eval(arm.pcB,t)[2] for t in (lo,hi)])≈
            [0.0,slant] atol=1e-9
    end

    # a tilted plane that cuts every ruling: an ellipse, one closed branch
    tilted=(occ=:plane,center=(0.8,0.55,0.3),axis=(sin(0.2),0.0,cos(0.2)),
            X=(cos(0.2),0.0,-sin(0.2)),Y=(0.0,1.0,0.0))
    ellipse=M._brep_sec_cone_conic(tilted,cone,true,"test")
    @test length(ellipse)==1 && ellipse[1].closed
    off_plane,apart,_=worst_lift(ellipse[1],tilted,cone,57)
    @test off_plane<1e-14 && apart<1e-9

    # spiric sections of a torus: one oval when the plane passes between
    # the inner and the outer equator, two loops around the tube inside
    ring=(occ=:torus,center=(0.8,0.55,0.35),axis=(0.0,0.0,1.0),
          X=(1.0,0.0,0.0),Y=(0.0,1.0,0.0),r1=0.39,r2=0.21)
    oval=M._brep_sec_spiric(wall,ring,true,"test")
    @test length(oval)==1 && oval[1].closed
    off_plane,apart,_=worst_lift(oval[1],wall,ring,57)
    @test off_plane<1e-14 && apart<1e-9
    near=merge(wall,(center=(0.9,0.0,0.0),))
    loops=M._brep_sec_spiric(near,ring,true,"test")
    @test length(loops)==2
    @test all(s->worst_lift(s,near,ring,57)[2]<1e-9,loops)
end

# The Booleans themselves walk fitted traces through trimming and region
# building, which costs minutes interpreted; they run compiled only. Values
# are closed forms (the cone's are also Gmsh 4.15.2's to every digit; the
# torus's agree with Gmsh to its 1e-6, and Gmsh cannot mesh that cut:
# "The 1D mesh seems not to be forming a closed loop").
if Base.JLOptions().compile_enabled==1
@testset "OCC Booleans with conic and spiric sections" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        box=occ.add_box(-1,-1,-1,2,2,2)
        cone=occ.add_cone(0.8,0.55,-0.25,0,0,1.2,0.6,0.24)
        out,_=occ.fragment([(3,box)],[(3,cone)])
        occ.synchronize()
        volumes=sort([occ.get_mass(3,t) for (_,t) in out])
        @test volumes≈[0.174129917,0.531597457,7.468402543] rtol=2e-9
        @test sum(volumes[1:2])≈π*1.2*(0.36+0.144+0.0576)/3 rtol=1e-11
        @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(16,27,15,3)
    finally
        A.finalize()
    end

    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        box=occ.add_box(-1,-1,-1,2,2,2)
        ring=occ.add_torus(0.8,0.55,0.35,0.39,0.21)
        out,_=occ.cut([(3,box)],[(3,ring)])
        occ.synchronize()
        @test occ.get_mass(3,out[1][2])≈7.789031646176062 rtol=1e-9
        @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(14,21,7,1)
        A.option("Mesh.MeshSizeMax",0.3)
        A.mesh.generate(2)
        @test sum(length(A.mesh.get_elements(2,s)[2][1])
                  for (_,s) in A.model.get_entities(2))>500
    finally
        A.finalize()
    end
end
end
end
