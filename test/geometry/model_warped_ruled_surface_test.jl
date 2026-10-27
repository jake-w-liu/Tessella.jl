module ModelWarpedRuledSurfaceChecks
using Test,Tessella,LinearAlgebra
const MT=Tessella.MeshTypes
const M2=Tessella.Mesh2D

# A genuinely non-coplanar four-wire patch: the fourth corner lifts out of
# the z=0 plane, so no planar projection exists and the surface must mesh in
# its own ruled parameter domain with the first-fundamental metric.
const WARPED_SMALL="""
Point(1) = {0, 0, 0, 0.5};
Point(2) = {1, 0, 0, 0.5};
Point(3) = {1, 1, 0, 0.5};
Point(4) = {0, 1, 0.5, 0.5};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Surface(1) = {1};
"""
const WARPED_LARGE="""
Point(1) = {0, 0, 0, 0.5};
Point(2) = {10, 0, 0, 0.5};
Point(3) = {10, 10, 0, 0.5};
Point(4) = {0, 10, 5.0, 0.5};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Surface(1) = {1};
"""

function execute_source(source;extra="")
    return mktemp() do path,io
        write(io,source*extra)
        close(io)
        execute_geo(path;mesh_dim=2)
    end
end

# The patch is the bilinear map of the unit square: u=x/scale, v=y/scale,
# z = z4·(1-u)·v. Every lifted node must sit on it to machine precision.
function surface_residual(mesh,scale,z4)
    worst=0.0
    for column in axes(mesh.coords,2)
        x,y,z=mesh.coords[1,column],mesh.coords[2,column],mesh.coords[3,column]
        u=clamp(x/scale,0.0,1.0);v=clamp(y/scale,0.0,1.0)
        worst=max(worst,abs(z-z4*(1.0-u)*v))
    end
    return worst
end

function min_angle(mesh)
    worst=Inf
    for cell in eachcol(mesh.tris)
        p=(mesh.coords[:,cell[1]],mesh.coords[:,cell[2]],mesh.coords[:,cell[3]])
        for k in 1:3
            a=p[k];b=p[mod1(k+1,3)];c=p[mod1(k+2,3)]
            v1=b-a;v2=c-a
            cosine=dot(v1,v2)/(norm(v1)*norm(v2))
            worst=min(worst,acosd(clamp(cosine,-1.0,1.0)))
        end
    end
    return worst
end

function max_edge(mesh)
    worst=0.0
    for cell in eachcol(mesh.tris)
        p=(mesh.coords[:,cell[1]],mesh.coords[:,cell[2]],mesh.coords[:,cell[3]])
        for k in 1:3
            worst=max(worst,norm(p[mod1(k+1,3)]-p[k]))
        end
    end
    return worst
end

@testset "warped ruled surface fills its own parameter domain" begin
    # Frontal-Delaunay (the default Mesh.Algorithm=6) on the unit warped
    # patch. Gmsh 4.15.2 reports 14 nodes for this fixture; the element
    # count folds the ten boundary segments into the total. Pinned exactly:
    # the construction is deterministic.
    first_run=execute_source(WARPED_SMALL)
    mesh=first_run.mesh
    @test mesh isa Mesh
    @test MT.validate(mesh).ok
    @test MT.nnodes(mesh)==14
    @test MT.ntris(mesh)==16
    @test MT.nsegs(mesh)==10
    @test surface_residual(mesh,1.0,0.5)<1e-12
    @test min_angle(mesh)>20.0
    @test max_edge(mesh)<0.75
    # The four declared corner coordinates survive at their exact positions.
    corners=Set((mesh.coords[1,i],mesh.coords[2,i],mesh.coords[3,i])
                for i in axes(mesh.coords,2))
    for corner in ((0.0,0.0,0.0),(1.0,0.0,0.0),(1.0,1.0,0.0),(0.0,1.0,0.5))
        @test corner in corners
    end
    # Deterministic: a second execution reproduces the identical product.
    second_run=execute_source(WARPED_SMALL)
    @test second_run.mesh.coords==mesh.coords
    @test second_run.mesh.tris==mesh.tris

    # The Delaunay arm (Mesh.Algorithm=5) and AUTO resolution (=2 → the
    # non-planar arm) both fill the same warped parameter domain on-surface.
    delaunay=execute_source(WARPED_SMALL;extra="Mesh.Algorithm = 5;\n")
    @test delaunay.mesh isa Mesh && MT.validate(delaunay.mesh).ok
    @test MT.ntris(delaunay.mesh)>0
    @test surface_residual(delaunay.mesh,1.0,0.5)<1e-12
    auto=execute_source(WARPED_SMALL;extra="Mesh.Algorithm = 2;\n")
    @test auto.mesh isa Mesh && MT.validate(auto.mesh).ok
    @test surface_residual(auto.mesh,1.0,0.5)<1e-12
end

@testset "metric ground floor: repeated eigenvalue and eigenvector sign" begin
    # The warped patch's first fundamental form and the bgm=15 metric at
    # (0.9,0.5): metric = G/225, so (metric, G/h²) is a proportional pair.
    # det(M−λN) cancels to ~−3e−11 in Float64; the floor must clamp the
    # discriminant and return the pulled-back bound G/h², not I/h².
    ground=(106.25,-1.25,100.25)
    metric=ground./225.0
    floored=M2._gmsh_metric_ground_floor(metric,ground,0.5)
    @test floored[1]≈4*ground[1] rtol=1e-10
    @test floored[2]≈4*ground[2] rtol=1e-10
    @test floored[3]≈4*ground[3] rtol=1e-10
    @test floored[1]>100.0   # pulled back, not the 4·I world-aligned floor

    # A straddling generalized spectrum (λmin<1<λmax) with n12≠0 exercises
    # the second null-space candidate — the eigenvector whose y component is
    # b−λmin·n12. The returned floor's generalized spectrum must be exactly
    # (1, λmax): the weak direction is lifted to the bound and the strong
    # direction is untouched.
    n=(1.0,0.4,1.0)
    a,b,d=0.5,0.0,3.0
    detn=n[1]*n[3]-n[2]*n[2]
    cB=a*n[3]+d*n[1]-2.0*b*n[2]
    detm=a*d-b*b
    disc=sqrt(cB*cB-4.0*detn*detm)
    lmax=(cB+disc)/(2.0*detn)
    result=M2._gmsh_metric_ground_floor((a,b,d),n,1.0)
    ra,rb,rd=result
    rcB=ra*n[3]+rd*n[1]-2.0*rb*n[2]
    rdetm=ra*rd-rb*rb
    rdisc=sqrt(max(rcB*rcB-4.0*detn*rdetm,0.0))
    rlmin=(rcB-rdisc)/(2.0*detn)
    rlmax=(rcB+rdisc)/(2.0*detn)
    @test rlmin≈1.0 atol=1e-10
    @test rlmax≈lmax rtol=1e-10
end

@testset "warped ruled surface at reference scale" begin
    # The 10×10 patch at lc=0.5. Gmsh 4.15.2's Frontal-Delaunay produces
    # 573 nodes and 1148 elements (≈1058 surface triangles); the native
    # construction lands at the same density. Ranges rather than exact
    # counts: the Newton surface-circle iterate can drift a node across
    # platforms without changing the coverage contract.
    run=execute_source(WARPED_LARGE)
    mesh=run.mesh
    @test mesh isa Mesh
    @test MT.validate(mesh).ok
    @test 500<=MT.nnodes(mesh)<=650
    @test 900<=MT.ntris(mesh)<=1250
    @test surface_residual(mesh,10.0,5.0)<1e-11
    @test min_angle(mesh)>20.0
    @test max_edge(mesh)<1.0
end
end
