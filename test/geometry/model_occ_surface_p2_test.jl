# Order-2 placement on OCC analytic revolution faces (HO01 slice).
#
# High-order edge/face supports on a curved face must evaluate the actual
# OCC surface record — Gmsh's high-order pass places the support at the
# projected (u,v) instead of the straight first-order interpolation.  OCC
# faces are consumed inside the volume PLC during `generate(3)`, so the
# observable route is a mesh attached to the face entity through
# `add_nodes`/`add_elements_by_type` (matching `add_discrete_*` semantics on
# native entities) followed by `set_order(2)`.
#
# Fixture chords are short enough for a unique closest point: on the unit
# cylinder a chord across a 25° arc sits at rho=cos(25°)≈0.906 while the
# projected support returns to rho=1.  Corners share no common circle (the
# third corner is off the z-level of the first edge), so the curved Tri6
# keeps a strictly positive Jacobian — an all-arc-inscribed quadratic
# triangle is the classic borderline-invalid element and is rightly
# refused by the exact certificate.

using Test
using Tessella
using LinearAlgebra: norm
const API=Tessella.API

# Cylinder r=1 about +z; lateral Surface[1] is the analytic :cylinder face.
# Corners (u,z): (0°,0.3), (50°,0.3), (20°,0.8).
function _occ_p2_cylinder_corners()
    return Float64[cosd(0.),sind(0.),0.3,
                   cosd(50.),sind(50.),0.3,
                   cosd(20.),sind(20.),0.8]
end

function _occ_p2_cylinder_fixture(api)
    api.model.add_cylinder(0.,0.,0.,0.,0.,1.,1.;tag=1)
    api.mesh.add_nodes(2,1,[101,102,103],_occ_p2_cylinder_corners())
    api.mesh.add_elements_by_type(1,2,[201],[101,102,103])
    return nothing
end

# Elevated supports are the nodes carrying a stored (u,v) parameter.
_occ_p2_supports(api,tags)=[t for t in tags if length(api.mesh.get_node(t)[2])==2]

@testset "P2 supports on an OCC cylinder face project to the surface" begin
    api=API
    api.initialize()
    try
        _occ_p2_cylinder_fixture(api)
        api.mesh.set_order(2)
        tags,xyz,uv=api.mesh.get_nodes(2,1,false,true)
        @test length(tags)==6
        @test api.mesh.get_element_types(2,1)==Int32[9]
        supports=_occ_p2_supports(api,tags)
        @test length(supports)==3
        for t in tags
            p,par,dim,ent=api.mesh.get_node(t)
            @test (dim,ent)==(2,1)
            # Corners were attached exactly on rho=1; projected supports
            # land back on the face instead of the rho<1 chord midpoint.
            @test hypot(p[1],p[2])≈1.0 atol=1e-12
            @test 0.29<=p[3]<=0.81
            if t in supports
                @test length(par)==2
                # Stored (u,v) evaluates back to the published coordinate.
                @test collect(api.model.get_value(2,1,par))≈p atol=1e-12
            end
        end
        # Order-1 round-trip restores the three corner nodes.
        api.mesh.set_order(1)
        tags1,xyz1,_=api.mesh.get_nodes(2,1,false,true)
        @test length(tags1)==3
        @test reshape(xyz1,3,3)≈reshape(_occ_p2_cylinder_corners(),3,3)
    finally
        api.finalize()
    end
end

@testset "P2 supports on an OCC sphere face project to the sphere" begin
    api=API
    api.initialize()
    try
        api.model.add_sphere(0.,0.,0.,1.;tag=1)
        # Corners on three different great-circle arcs: no common circle.
        pts=[(1.,0.,0.),(0.,1.,0.),(cosd(40.)*sind(50.),sind(40.)*sind(50.),cosd(50.))]
        api.mesh.add_nodes(2,1,[101,102,103],collect(Iterators.flatten(pts)))
        api.mesh.add_elements_by_type(1,2,[201],[101,102,103])
        api.mesh.set_order(2)
        tags,xyz,uv=api.mesh.get_nodes(2,1,false,true)
        @test length(tags)==6
        supports=_occ_p2_supports(api,tags)
        @test length(supports)==3
        for t in supports
            p,par,dim,ent=api.mesh.get_node(t)
            # Chord midpoints of these arcs sit strictly inside |p|<1.
            @test norm(p)≈1.0 atol=1e-12
            @test collect(api.model.get_value(2,1,par))≈p atol=1e-12
        end
    finally
        api.finalize()
    end
end

@testset "SecondOrderLinear keeps straight placement on OCC faces" begin
    api=API
    api.initialize()
    try
        _occ_p2_cylinder_fixture(api)
        api.option("Mesh.SecondOrderLinear",1)
        api.mesh.set_order(2)
        tags,xyz,_=api.mesh.get_nodes(2,1,false,true)
        @test length(tags)==6
        # Interpolated edge supports stay on the chords (rho<1) and carry
        # no stored face parameters.
        @test all(t->hypot(api.mesh.get_node(t)[1][1],
                           api.mesh.get_node(t)[1][2])<1.0-1e-3 &&
                      isempty(api.mesh.get_node(t)[2]),tags[4:6])
    finally
        api.option("Mesh.SecondOrderLinear",0)
        api.finalize()
    end
end

@testset "Mixed-cache quadratic elevation projects OCC face supports" begin
    # The cache elevation path used by `_native_mixed_set_order_labels_plan`
    # applies the same projection for cells already classified on the face.
    api=API
    api.initialize()
    try
        api.model.add_cylinder(0.,0.,0.,0.,0.,1.,1.;tag=1)
        model=api.CURRENT[]
        linear=MixedMesh(reshape(_occ_p2_cylinder_corners(),3,3),
            [ElementBlock(2,Int32[1;2;3;;])])
        owner=(2,Int32(1))
        class=api._mixed_classification(linear,owner,[owner],fill(owner,3),
            Dict(owner=>Int32[]),Dict(2=>Int32[1]))
        quadratic,_=api._mixed_quadratic_cache(linear,class,"test";model)
        @test quadratic.blocks[1].msh==9
        @test validate(quadratic).ok
        for k in 4:6
            @test hypot(quadratic.coords[1,k],quadratic.coords[2,k])≈1.0 atol=1e-12
        end
    finally
        api.finalize()
    end
end
