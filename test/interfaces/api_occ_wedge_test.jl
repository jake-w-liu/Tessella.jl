# `model.occ.add_wedge` — `BRepPrimAPI_MakeWedge`: the xz trapezoid
# (0,0)→(dx,0)→(ltx,dz)→(0,dz) in the `gp_Ax2(p,zAxis)` frame extruded
# `dy` along +y. `ltx==0` collapses the top edge into the ridge vertex
# (triangular profile, 5 faces); `ltx==dx` reproduces the box layout.
# Signed extents wind outward consistently (verified through the
# divergence-path `get_mass`), zero extents reject "Degenerate wedge",
# and `zAxis` follows the strict 3-vector contract like `add_torus`.

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

@testset "occ.add_wedge" begin
    occ=_occ_setup()
    # Standard ramp: triangular profile (0,0),(2,0),(0,4) extruded 3.
    v=occ.add_wedge(0.0,0.0,0.0,2.0,3.0,4.0)
    @test v==1
    occ.synchronize()
    m=_current_model()
    # 6 vertices, 9 edges, 5 wires, 5 plane faces, 1 shell, 1 solid.
    @test length(m.points)==6
    @test length(m.curves)==9
    @test length(m.surfaces)==5
    @test length(m.surface_loops)==1
    @test m.volumes==Dict(1=>[1])
    @test sort(m.surface_loops[1])==collect(1:5)
    # Two triangular caps, three quad sides.
    wire_sizes=sort!([length(m.loops[only(m.surfaces[t])])
                      for t in keys(m.surfaces)])
    @test wire_sizes==[3,3,4,4,4]
    # Volume = dz·dy·(dx+ltx)/2 = 4·3·2/2 = 12; prism centroid at the
    # profile triangle centroid plus dy/2.
    @test occ.get_mass(3,v)≈12.0
    @test collect(occ.get_center_of_mass(3,v))≈collect((2/3,1.5,4/3))
    xmin,ymin,zmin,xmax,ymax,zmax=A.model.get_bounding_box(3,v)
    @test collect((xmin,ymin,zmin))≈collect((0.0,0.0,0.0)) atol=1e-12
    @test collect((xmax,ymax,zmax))≈collect((2.0,3.0,4.0))
    # Outward normals: caps ±y, bottom −z, back −x, slope (dz,0,dx)/√.
    expected=((0.0,-1.0,0.0),(0.0,1.0,0.0),(0.0,0.0,-1.0),
              (-1.0,0.0,0.0),(4.0,0.0,2.0)./sqrt(20.0))
    axes=[m.surface_geometry[t].axis for t in 1:5]
    @test length(axes)==5
    for e in expected
        @test any(a->isapprox(collect(a),collect(e);atol=1e-9),axes)
    end
    # Manifold shell: every edge appears in exactly two face wires.
    edge_faces=Dict{Int,Vector{Int}}()
    for (ftag,loops) in m.surfaces
        for ctag in m.loops[only(loops)]
            push!(get!(edge_faces,abs(ctag),Int[]),ftag)
        end
    end
    @test sort!(collect(keys(edge_faces)))==collect(1:9)
    @test all(length(unique(fs))==2 for fs in values(edge_faces))
end

@testset "occ.add_wedge variants" begin
    occ=_occ_setup()
    # Truncated wedge: trapezoid (0,0),(3,0),(1,4),(0,4), area 8, V=16;
    # polygon centroid (13/12, 5/3) in xz, dy/2 in y.
    v=occ.add_wedge(0.0,0.0,0.0,3.0,2.0,4.0;ltx=1.0)
    occ.synchronize()
    m=_current_model()
    @test length(m.points)==8
    @test length(m.surfaces)==6
    @test sort!([length(m.loops[only(m.surfaces[t])])
                 for t in keys(m.surfaces)])==[4,4,4,4,4,4]
    @test occ.get_mass(3,v)≈16.0
    @test collect(occ.get_center_of_mass(3,v))≈
        collect((13/12,1.0,5/3))
    # The flat top adds a +z normal to the set.
    @test (0.0,0.0,1.0) in [m.surface_geometry[t].axis for t in keys(m.surfaces)]
    A.initialize(); occ=_occ_setup()
    # ltx==dx reproduces the box: six rectangles, V=dx·dy·dz.
    vb=occ.add_wedge(1.0,2.0,3.0,2.0,2.0,2.0;ltx=2.0)
    @test occ.get_mass(3,vb)≈8.0
    @test collect(occ.get_center_of_mass(3,vb))≈collect((2.0,3.0,4.0))
    # ltx>dx is a valid overhanging trapezoid: area = dz·(dx+ltx)/2 = 7.
    vo=occ.add_wedge(0.0,0.0,0.0,2.0,1.0,2.0;ltx=5.0)
    @test occ.get_mass(3,vo)≈7.0
    xmin,_,_,xmax,_,_=occ.get_bounding_box(3,vo)
    @test xmin≈0.0
    @test xmax≈5.0
    # All-negative extents stay outward-consistent (positive mass).
    vn=occ.add_wedge(0.0,0.0,0.0,-2.0,-3.0,-4.0)
    @test occ.get_mass(3,vn)≈12.0
    xmin,ymin,zmin,xmax,ymax,zmax=occ.get_bounding_box(3,vn)
    @test collect((xmin,ymin,zmin))≈collect((-2.0,-3.0,-4.0)) atol=1e-12
    @test collect((xmax,ymax,zmax))≈collect((0.0,0.0,0.0)) atol=1e-12
    # zAxis frame: local z→x, so dz maps onto world x and volume keeps.
    vz=occ.add_wedge(0.0,0.0,0.0,2.0,3.0,4.0;zAxis=[1.0,0.0,0.0])
    @test occ.get_mass(3,vz)≈12.0
    xmin,ymin,zmin,xmax,ymax,zmax=occ.get_bounding_box(3,vz)
    @test xmax≈4.0
    @test ymax≈0.0 atol=1e-12
    @test ymin≈-3.0
    @test zmax≈2.0
    @test zmin≈0.0 atol=1e-12
end

@testset "occ.add_wedge validation + lifecycle" begin
    occ=_occ_setup()
    v=occ.add_wedge(0.0,0.0,0.0,1.0,1.0,1.0)
    before=length(occ.get_entities())
    # Zero extents reject before any sub-entity allocation.
    for ext in ((0.0,1.0,1.0),(1.0,0.0,1.0),(1.0,1.0,0.0))
        @test occursin("Degenerate wedge",
            _errtext(()->occ.add_wedge(0.0,0.0,0.0,ext...)))
    end
    @test length(occ.get_entities())==before
    # Non-finite inputs reject on each slot.
    @test occursin("must be finite",
        _errtext(()->occ.add_wedge(NaN,0.0,0.0,1.0,1.0,1.0)))
    @test occursin("must be finite",
        _errtext(()->occ.add_wedge(0.0,0.0,0.0,1.0,1.0,Inf)))
    @test occursin("ltx must be finite",
        _errtext(()->occ.add_wedge(0.0,0.0,0.0,1.0,1.0,1.0;ltx=NaN)))
    # zAxis contract matches add_torus.
    @test occursin("zAxis must be a 3-vector",
        _errtext(()->occ.add_wedge(0.0,0.0,0.0,1.0,1.0,1.0;zAxis=[1.0,0.0])))
    @test occursin("zAxis must have positive length",
        _errtext(()->occ.add_wedge(0.0,0.0,0.0,1.0,1.0,1.0;
                                   zAxis=[0.0,0.0,0.0])))
    # Tag honoring + collision; failed calls consume nothing.
    @test length(occ.get_entities())==before
    @test occ.add_wedge(9.0,9.0,9.0,1.0,1.0,1.0;tag=7)==7
    @test occursin("OCC Volume[7] already exists",
        _errtext(()->occ.add_wedge(0.0,0.0,0.0,1.0,1.0,1.0;tag=7)))
    # Auto tags are getMaxTag(3)+1 — 8 after 1 and 7 exist.
    @test occ.add_wedge(0.0,0.0,0.0,1.0,1.0,1.0)==8
    # Recursive removal carries the shell's GModel loop away.
    occ.synchronize()
    m=_current_model()
    @test haskey(m.volumes,v)
    @test length(m.surface_loops)==3
    occ.remove([(3,v)],true)
    occ.synchronize()
    m=_current_model()
    @test !haskey(m.volumes,v)
    @test length(m.surface_loops)==2
end

@testset "occ.add_wedge transforms + meshing" begin
    occ=_occ_setup()
    v=occ.add_wedge(0.0,0.0,0.0,2.0,3.0,4.0)
    # Kernel-side queries see the transform before synchronize.
    occ.translate([(3,v)],10.0,0.0,0.0)
    xmin,_,_,xmax,_,_=occ.get_bounding_box(3,v)
    @test xmin≈10.0
    @test xmax≈12.0
    @test occ.get_mass(3,v)≈12.0
    occ.translate([(3,v)],-10.0,0.0,0.0)
    occ.synchronize()
    A.mesh.generate(3)
    node_tags,coords=A.mesh.get_nodes(3,v,true)
    @test length(node_tags)>0
    for i in 1:3:3length(node_tags)
        x,y,z=coords[i],coords[i+1],coords[i+2]
        @test -1e-9<=x<=2.0+1e-9
        @test -1e-9<=y<=3.0+1e-9
        @test -1e-9<=z<=4.0+1e-9
    end
    element_types,element_tags,_=A.mesh.get_elements(3,v)
    @test element_types==Int32[4]
    @test length(element_tags[1])>0
end
