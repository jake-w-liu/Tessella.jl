using Test, Tessella, LinearAlgebra

const _LP2B=Tessella.API

function _lp2b_mesh(dimension,beta;scale=1.0,reflected=false)
    reference=dimension==2 ?
        ((0.,0.,0.),(1.,0.,0.),(0.,1.,0.),(.5,0.,0.),(.5,.5,0.),(0.,.5,0.)) :
        ((0.,0.,0.),(1.,0.,0.),(0.,1.,0.),(0.,0.,1.),
         (.5,0.,0.),(.5,.5,0.),(0.,.5,0.),(0.,0.,.5),(0.,.5,.5),(.5,0.,.5))
    points=reduce(hcat,(begin
        x=u+beta/2*((u-.23)^2-(v-.31)^2)
        y=v-beta*(u-.23)*(v-.31)
        scale .* [reflected ? -x : x,y,w]
    end for (u,v,w) in reference))
    # The represented maps below have positive exact Jacobians. Structural
    # construction permits the independent reflected/extreme-scale variants.
    return dimension==2 ?
        Tessella.HighOrder.P2TriMesh(points,reshape(Int32.(1:6),6,1);_check_geometry=false) :
        Tessella.HighOrder.P2Mesh(points,reshape(Int32.(1:10),10,1);
            require_positive_tets=false,_check_geometry=false)
end

function _lp2b_quality(mesh,name)
    vertices,scale=_LP2B._legacy_p2_quality_vertices(mesh,1,"bounds regression")
    normal=mesh isa Tessella.HighOrder.P2TriMesh ?
        first(_LP2B._legacy_p2_primary_normal(mesh,1,scale)) : nothing
    return _LP2B._legacy_p2_bounded_quality(mesh,vertices,normal,scale,name,"bounds regression")
end

@testset "Actual simplex P2 Bernstein quality bounds" begin
    # Independent polynomial:
    # x=u+b/2*((u-a)^2-(v-c)^2), y=v-b*(u-a)*(v-c), z=w.
    # Its determinant is 1-b^2*((u-a)^2+(v-c)^2). The strict interior
    # maximum is 1 at (a,c)=(.23,.31); all original corners are lower.
    # The fixed values are from Gmsh 4.15.2 on these exact represented nodes.
    # Upstream reports converged Bernstein bounds, including its 1e-3
    # relative gap, and the documented ICN lower/corner-bound blend.
    pinned=(
        (2,.8,(.55904,1.0006399999999998,.30299162315884115)),
        (2,1.1,(.16630999999999896,1.0008318750000003,.09648934802496811)),
        (3,.8,(.5590399999999989,1.0006399999999995,.4334655234738617)),
        (3,1.1,(.16630999999999763,1.000642812500001,.1808937919447706)))
    for (dimension,beta,expected) in pinned
        mesh=_lp2b_mesh(dimension,beta)
        observed=ntuple(i->_lp2b_quality(mesh,
            ("minDetJac","maxDetJac","minIsotropy")[i]),3)
        @test collect(observed)≈collect(expected) atol=2e-11 rtol=2e-11
        exact_min=1-beta^2*max(.23^2+.31^2,.77^2+.31^2,.23^2+.69^2)
        @test observed[1]≈exact_min atol=2e-13 rtol=2e-13
        @test 1<=observed[2]<=1.0011
        @test 0<observed[3]<1
        # The initial corner samples cannot establish the interior maximum.
        vertices,scale=_LP2B._legacy_p2_quality_vertices(mesh,1,"bounds regression")
        normal=dimension==2 ? first(_LP2B._legacy_p2_primary_normal(mesh,1,scale)) : nothing
        initial=_LP2B._legacy_p2_jac_domain(vertices,normal)
        @test initial.max_l*scale^dimension<1
        @test !_LP2B._legacy_p2_jac_bounds_ok(initial,initial.min_l,initial.max_l)
        # Uniform rescaling changes determinants and preserves isotropy.
        scaled=_lp2b_mesh(dimension,beta;scale=8.)
        @test _lp2b_quality(scaled,"minDetJac")≈8.0^dimension*observed[1] atol=2e-10 rtol=2e-12
        @test _lp2b_quality(scaled,"maxDetJac")≈8.0^dimension*observed[2] atol=2e-10 rtol=2e-12
        @test _lp2b_quality(scaled,"minIsotropy")≈observed[3] atol=2e-12 rtol=2e-12
        reflected=_lp2b_mesh(dimension,beta;reflected=true)
        if dimension==2
            # The primary surface normal follows its reflected winding.
            @test _lp2b_quality(reflected,"minDetJac")≈observed[1] atol=2e-12 rtol=2e-12
            @test _lp2b_quality(reflected,"minIsotropy")≈observed[3] atol=2e-12 rtol=2e-12
        else
            @test _lp2b_quality(reflected,"minDetJac")≈-observed[2] atol=2e-12 rtol=2e-12
            @test _lp2b_quality(reflected,"maxDetJac")≈-observed[1] atol=2e-12 rtol=2e-12
            @test _lp2b_quality(reflected,"minIsotropy")==0.
        end
    end
    for dimension in (2,3)
        affine=_lp2b_mesh(dimension,0.)
        @test _lp2b_quality(affine,"minDetJac")≈1. atol=2e-14 rtol=2e-14
        @test _lp2b_quality(affine,"maxDetJac")≈1. atol=2e-14 rtol=2e-14
        expected_isotropy=dimension==2 ? sqrt(3.)/2 : 2.0^(4/3)/3
        @test _lp2b_quality(affine,"minIsotropy")≈expected_isotropy atol=2e-14 rtol=2e-14
        near=_lp2b_mesh(dimension,2.0^-20)
        @test abs(_lp2b_quality(near,"minDetJac")-1)<2e-12
        @test abs(_lp2b_quality(near,"maxDetJac")-1)<2e-12
        @test abs(_lp2b_quality(near,"minIsotropy")-expected_isotropy)<1e-5
        for scale in (2.0^-400,2.0^400)
            extreme=_lp2b_mesh(dimension,.8;scale)
            @test _lp2b_quality(extreme,"minIsotropy")≈
                _lp2b_quality(_lp2b_mesh(dimension,.8),"minIsotropy") atol=2e-12 rtol=2e-12
            @test _lp2b_quality(extreme,"minDetJac")>0
            determinant=_lp2b_quality(extreme,"maxDetJac")
            # Existing quality conversion preserves sign on underflow and
            # reports infinity for a quantity exceeding Float64's range.
            @test !isnan(determinant) && determinant>0
        end
        frame=((0.,0.,0.),(0.,0.,0.),(0.,0.,0.))
        zero_vertices=ntuple(_->frame,dimension+1)
        normal=dimension==2 ? (0.,0.,1.) : nothing
        for name in ("minDetJac","maxDetJac","minIsotropy")
            @test _LP2B._legacy_p2_bounded_quality(affine,zero_vertices,normal,1.,name,"zero regression")==0.
        end
    end
end
