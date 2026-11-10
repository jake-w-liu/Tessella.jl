module OCCConicRationalChartChecks
using Test,Tessella
const A=Tessella.API
const M=Tessella.Model

# Independent inverse tangent-half-angle chart; opaque gmsh 4.15.2 SDK
# controls agree at identical public parameters for circles and quadrics.
function chart(parameter,lower,upper)
    count=max(1,ceil(Int,(upper-lower)/(2π/3)))
    width=(upper-lower)/count
    span=clamp(floor(Int,(parameter-lower)/width),0,count-1)
    start=lower+span*width
    q=tan(width/4);z=(2*(parameter-start)/width-1)*q
    first=4q/(width*(1+z*z))
    return start+width/2+2atan(z),first,-first*first*z
end

@testset "Converted conics expose rational public parameters" begin
    A.initialize()
    try
        for (lower,upper,rx,ry) in ((.1,1.1,1.,1.),(0.,2π,1.,1.),
                (.4,.4+2π,1.,1.),(.1,3.7,3.,1.),
                (-2π,0.,1.,1.),(-1.7,.9,3.,1.))
            tag=rx==ry ? A.model.occ.add_circle(0.,0.,0.,rx;angle1=lower,angle2=upper) :
                A.model.occ.add_ellipse(0.,0.,0.,rx,ry;angle1=lower,angle2=upper)
            A.model.occ.synchronize()
            probe=lower+.23*(upper-lower)
            angular=A.model.get_value(1,tag,[probe])
            A.model.occ.convert_to_nurbs([(1,tag)]);A.model.occ.synchronize()
            @test A.model.get_type(1,tag)=="BSpline"
            @test maximum(abs,A.model.get_value(1,tag,[probe]).-angular)>1e-4
            lo,hi=A.model.get_parametrization_bounds(1,tag)
            @test lo≈[lower] && hi≈[upper]
            geometry=A.CURRENT[].curve_geometry[tag]
            # Avoid exact interior rational-knot parameters for the analytic
            # internal D2 comparison: at a double knot the quadratic rational
            # spline is C1 but its pointwise second derivative jumps, while the
            # public finite-difference query is separately checked below.
            for parameter in (lower,probe,lower+.4*(upper-lower),upper)
                theta,d,dd=chart(parameter,lower,upper)
                point=[rx*cos(theta),ry*sin(theta),0.]
                tangent=[-rx*sin(theta),ry*cos(theta),0.]
                @test A.model.get_value(1,tag,[parameter])≈point atol=2e-14
                @test A.model.get_derivative(1,tag,[parameter])≈d.*tangent atol=3e-14
                # Internal rational-curve D2 is analytic and must match the
                # exact conic independently of the public query contract.
                @test collect(M._model_occ_curve_d2(geometry.curve,parameter)[3])≈
                    dd.*tangent.-d*d.*point atol=8e-14
                # OCC's public GEdge::secondDer is a 1e-3 difference of the
                # analytic first derivative (opaque gmsh 4.15.2 controls match
                # it exactly, including one-sided differences at the bounds).
                eps=1e-3
                public=A.model.get_second_derivative(1,tag,[parameter])
                first=A.model.get_derivative(1,tag,[parameter])
                expected=if parameter-eps<=lo[1]
                    second=A.model.get_derivative(1,tag,[parameter+eps])
                    ntuple(i->1000*(second[i]-first[i]),3)
                elseif parameter+eps>=hi[1]
                    second=A.model.get_derivative(1,tag,[parameter-eps])
                    ntuple(i->1000*(first[i]-second[i]),3)
                else
                    second=A.model.get_derivative(1,tag,[parameter+eps])
                    first=A.model.get_derivative(1,tag,[parameter-eps])
                    ntuple(i->500*(second[i]-first[i]),3)
                end
                @test collect(public)≈collect(expected) atol=1e-12
            end
            # Interior rational knots are valid public query points.
            count=max(1,ceil(Int,(upper-lower)/(2π/3)))
            if count>1
                parameter=lower+(upper-lower)/count
                eps=1e-3
                public=A.model.get_second_derivative(1,tag,[parameter])
                first=A.model.get_derivative(1,tag,[parameter+eps])
                second=A.model.get_derivative(1,tag,[parameter-eps])
                expected=ntuple(i->500*(first[i]-second[i]),3)
                @test collect(public)≈collect(expected) atol=1e-12
            end
            expected_length=rx==ry ? rx*(upper-lower) : nothing
            expected_length===nothing || @test A.model.occ.get_mass(1,tag)≈expected_length rtol=2e-12
            @test M._brep_edge_curve(A.CURRENT[],tag)[1].kind===:bspline
            A.model.occ.convert_to_nurbs([(1,tag)]);A.model.occ.synchronize()
            @test A.CURRENT[].curve_geometry[tag]==geometry
        end
    finally
        A.finalize()
    end
end

@testset "Quadric rational charts transport retained and converted boundaries" begin
    for (kind,create) in ((:cylinder,()->A.model.occ.add_cylinder(0.,0.,0.,0.,0.,2.,1.)),
            (:cone,()->A.model.occ.add_cone(0.,0.,0.,0.,0.,2.,1.,.5)),
            (:sphere,()->A.model.occ.add_sphere(0.,0.,0.,1.)),
            (:torus,()->A.model.occ.add_torus(0.,0.,0.,3.,1.)))
        A.initialize()
        try
            volume=create();A.model.occ.synchronize()
            model=A.CURRENT[]
            face=only(tag for (tag,g) in model.surface_geometry if g.occ===kind)
            source=deepcopy(model.surface_geometry[face])
            lower,upper=A.model.get_parametrization_bounds(2,face)
            uv=lower .+ .17 .* (upper .- lower)
            angular=A.model.get_value(2,face,uv)
            mass=A.model.occ.get_mass(2,face)
            for selection in ([(2,face)],[(3,volume)])
                A.model.occ.convert_to_nurbs(selection);A.model.occ.synchronize()
                model=A.CURRENT[];g=model.surface_geometry[face]
                @test A.model.get_type(2,face)=="BSpline surface"
                @test maximum(abs,A.model.get_value(2,face,uv).-angular)>1e-4
                for fraction in (0.,.17,.37,.63,1.)
                    u,v=lower .+ fraction .* (upper .- lower)
                    a,da,dda=chart(u,lower[1],upper[1])
                    b,db,ddb=kind in (:sphere,:torus) ? chart(v,lower[2],upper[2]) : (v,1.,0.)
                    p,du,dv,duu,dvv,duv=M._occ_surface_d2(source,a,b)
                    expected=(p,da.*du,db.*dv,da^2 .* duu .+ dda .* du,
                        db^2 .* dvv .+ ddb .* dv,(da*db).*duv)
                    actual=M._occ_surface_d2(g,u,v)
                    for order in 1:6
                        @test collect(actual[order])≈collect(expected[order]) atol=8e-13
                    end
                end
                @test A.model.occ.get_mass(2,face)≈mass rtol=2e-12
                for loop in model.surfaces[face],signed in model.loops[loop]
                    entry=M._occ_face_pcurve(model,g,signed)
                    entry===nothing && continue
                    pc,t0,t1=entry
                    for parameter in (t0,t0+.23*(t1-t0),t1)
                        point=M._periodic_curve_point(model,abs(signed),parameter,"test")
                        trace=M._occ_pcurve_d0(pc,parameter)
                        @test collect(M._occ_surface_point(g,trace...))≈collect(point) atol=2e-12
                        @test all(isfinite,M._occ_pcurve_d1(pc,parameter)[2])
                    end
                end
            end
            # A pole transform preserves the existing rational knot coordinates.
            original=M._occ_surface_point(A.CURRENT[].surface_geometry[face],uv...)
            A.model.occ.dilate([(3,volume)],0.,0.,0.,2.,2.,2.)
            A.model.occ.synchronize()
            model=A.CURRENT[]
            moved=only(tag for (tag,g) in model.surface_geometry if g.occ===kind)
            g=model.surface_geometry[moved]
            @test collect(M._occ_surface_point(g,uv...))≈2 .*collect(original) atol=2e-12
            @test A.model.occ.get_mass(2,moved)≈4mass rtol=2e-12
            for loop in model.surfaces[moved],signed in model.loops[loop]
                entry=M._occ_face_pcurve(model,g,signed)
                entry===nothing && continue
                pc,t0,t1=entry;parameter=t0+.37*(t1-t0)
                point=M._periodic_curve_point(model,abs(signed),parameter,"test")
                @test collect(M._occ_surface_point(g,M._occ_pcurve_d0(pc,parameter)...))≈collect(point) atol=3e-12
            end
        finally
            A.finalize()
        end
    end
end

@testset "Reflected converted conics keep the rational parameter domain" begin
    A.initialize()
    try
        # Trimmed/full circles and an ellipse keep p'(t)=T·p(t).
        for (create,lower,upper) in (
                (()->A.model.occ.add_circle(0.,0.,0.,1.;angle1=.2,angle2=2.1),.2,2.1),
                (()->A.model.occ.add_circle(0.,0.,0.,1.;angle1=-2π,angle2=0.),-2π,0.),
                (()->A.model.occ.add_ellipse(0.,0.,0.,3.,1.;angle1=.1,angle2=3.7),.1,3.7))
            tag=create();A.model.occ.synchronize()
            A.model.occ.convert_to_nurbs([(1,tag)]);A.model.occ.synchronize()
            probe=lower+.31*(upper-lower)
            before=A.model.get_value(1,tag,[probe])
            A.model.occ.mirror([(1,tag)],1.,0.,0.,0.);A.model.occ.synchronize()
            lo,hi=A.model.get_parametrization_bounds(1,tag)
            @test lo≈[lower] && hi≈[upper]
            after=A.model.get_value(1,tag,[probe])
            @test collect(after)≈[-before[1],before[2],before[3]] atol=2e-13
            A.model.occ.convert_to_nurbs([(1,tag)]);A.model.occ.synchronize()
        end
        # Native-backed center-form arcs reflect through their native
        # geometry too, not only through the rational poles.
        p1=A.model.occ.add_point(1.,0.,0.)
        p2=A.model.occ.add_point(0.,0.,0.)
        p3=A.model.occ.add_point(0.,1.,0.)
        arc=A.model.occ.add_circle_arc(p1,p2,p3)
        A.model.occ.synchronize()
        A.model.occ.convert_to_nurbs([(1,arc)]);A.model.occ.synchronize()
        lo,hi=A.model.get_parametrization_bounds(1,arc)
        probe=lo[1]+.37*(hi[1]-lo[1])
        before=A.model.get_value(1,arc,[probe])
        A.model.occ.mirror([(1,arc)],1.,0.,0.,0.);A.model.occ.synchronize()
        @test collect(A.model.get_value(1,arc,[probe]))≈
            [-before[1],before[2],before[3]] atol=2e-13
        # A determinant-negative nonuniform affine warps a converted circle
        # into an exact ellipse on the same rational parameter domain.
        tag=A.model.occ.add_circle(0.,0.,0.,1.;angle1=0.,angle2=2π)
        A.model.occ.synchronize()
        A.model.occ.convert_to_nurbs([(1,tag)]);A.model.occ.synchronize()
        A.model.occ.affine_transform([(1,tag)],
            [-2.,0.,0.,0., 0.,1.,0.,0., 0.,0.,1.,0.])
        A.model.occ.synchronize()
        for parameter in (0.,.7,2.4,5.9)
            theta,_,_=chart(parameter,0.,2π)
            @test collect(A.model.get_value(1,tag,[parameter]))≈
                [-2cos(theta),sin(theta),0.] atol=3e-13
        end
        # Extreme anisotropic scales keep the ellipse branch and stable axes.
        small=A.model.occ.add_circle(0.,0.,0.,1.;angle1=0.,angle2=2π)
        A.model.occ.synchronize()
        A.model.occ.convert_to_nurbs([(1,small)]);A.model.occ.synchronize()
        A.model.occ.affine_transform([(1,small)],
            [-1e-15,0.,0.,0., 0.,2e-15,0.,0., 0.,0.,1e-15,0.])
        A.model.occ.synchronize()
        for parameter in (0.,.7,2.4,5.9)
            theta,_,_=chart(parameter,0.,2π)
            @test collect(A.model.get_value(1,small,[parameter]))≈
                [-1e-15*cos(theta),2e-15*sin(theta),0.] atol=2e-28
        end
        # A converted cylinder face reflects with exact value and mass.
        volume=A.model.occ.add_cylinder(0.,0.,0.,0.,0.,2.,1.)
        A.model.occ.synchronize()
        face=only(tag for (tag,g) in A.CURRENT[].surface_geometry if g.occ===:cylinder)
        A.model.occ.convert_to_nurbs([(3,volume)]);A.model.occ.synchronize()
        lo,hi=A.model.get_parametrization_bounds(2,face)
        uv=lo .+ .37 .* (hi .- lo)
        before=A.model.get_value(2,face,uv)
        mass=A.model.occ.get_mass(2,face)
        A.model.occ.mirror([(3,volume)],1.,0.,0.,0.);A.model.occ.synchronize()
        moved=only(tag for (tag,g) in A.CURRENT[].surface_geometry if g.occ===:cylinder)
        @test collect(A.model.get_value(2,moved,uv))≈
            [-before[1],before[2],before[3]] atol=2e-12
        @test A.model.occ.get_mass(2,moved)≈mass rtol=2e-12
    finally
        A.finalize()
    end
end
end
