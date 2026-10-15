# `model.occ` primitive solid endpoints: `add_cylinder`, `add_sphere`,
# `add_cone`, `add_torus` — the `BRepPrimAPI_Make*` solids bound into the OCC
# shape registry. Each publishes its full `BRepPrim_OneAxis` entity graph
# atomically (rim/pole vertices, circle/degenerate/seam edges, wires, the
# analytic face plus Plane caps, shell, solid) and `synchronize` imports it
# through the same arms as `add_box`, additionally retaining the compact
# `m.cylinders`/`m.spheres`/`m.cones` encoding so volume meshing flows through
# the native analytic tessellation like the geo-pipeline `add_*!` producers.
# Partial sweeps (`angle`/`angle1`/`angle2`/`angle3` away from the
# full-rotation defaults) reject loudly — the side-face layout has no
# verified reference yet. Layout mirrors `OCC_Internals::makeCylinder`/
# `makeSphere`/`makeCone`/`makeTorus` on Gmsh 4.15.2.

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

@testset "occ.add_cylinder" begin
    occ=_occ_setup()
    v=occ.add_cylinder(0.0,0.0,0.0, 0.0,0.0,2.0, 1.5)
    @test v==1
    occ.synchronize()
    m=_current_model()
    @test length(m.points)==2
    @test length(m.curves)==3
    @test length(m.surfaces)==3
    @test length(m.surface_loops)==1
    @test m.volumes==Dict(1=>[1])
    @test m.surface_types[1]==:cylinder
    @test m.surface_types[2]==:plane
    @test m.surface_types[3]==:plane
    # Two rim circles + one seam line; circles closed on their rim vertices.
    kinds=[m.curve_types[t] for t in 1:3]
    @test kinds==[:circle,:line,:circle]
    @test m.curves[1]==(1,1)
    @test m.curves[2]==(2,1)
    @test m.curves[3]==(2,2)
    # Wire [-top,-seam,+bottom,+seam]; shell [lateral,+top,-bottom].
    @test m.loops[m.surfaces[1][1]]==[-1,-2,3,2]
    @test m.surface_loops[1]==[1,2,-3]
    # Analytic record.
    g=m.surface_geometry[1]
    @test g.occ===:cylinder
    @test collect(g.center)≈[0.0,0.0,0.0]
    @test collect(g.axis)≈[0.0,0.0,1.0]
    @test g.radius==1.5
    @test g.height==2.0
    # Surface evaluation: u azimuth, v along axis.
    @test M.model_value(m,2,1,[0.0,0.0])≈[1.5,0.0,0.0]
    @test M.model_value(m,2,1,[π/2,1.0])≈[0.0,1.5,1.0] atol=1e-12
    @test M.model_value(m,2,1,[π,2.0])≈[-1.5,0.0,2.0] atol=1e-12
    # Parametrization bounds.
    lo,hi=M.model_parametrization_bounds(m,2,1)
    @test lo≈[0.0,0.0] && hi≈[2π,2.0]
    lo,hi=M.model_parametrization_bounds(m,2,2)
    @test lo≈[-1.5,-1.5] && hi≈[1.5,1.5]
    # Outward normals.
    @test A.model.get_normal(1,[0.0,1.0])≈[1.0,0.0,0.0] atol=1e-12
    @test A.model.get_normal(2,[0.0,0.0])≈[0.0,0.0,1.0]
    @test A.model.get_normal(3,[0.0,0.0])≈[0.0,0.0,-1.0] atol=1e-12
    # Pcurve -> 3-D round-trip on the lateral face.
    worst=0.0
    for (ct,e) in g.pcurves
        cg=m.curve_geometry[ct]
        for pc in (e.fwd,e.rev)
            pc===nothing && continue
            for t in range(cg.t0,cg.t1,length=5)
                uvv=M._occ_pcurve_d0(pc,t)
                srf=M.model_value(m,2,1,[uvv[1],uvv[2]])
                crv=M.model_value(m,1,ct,[t])
                worst=max(worst,maximum(abs.(srf.-crv)))
            end
        end
    end
    @test worst<1e-12
    # Compact encoding retained for the native analytic volume path.
    @test haskey(m.cylinders,v)
    @test m.cylinders[v].radius==1.5
    @test M._implicit_volume_surface(m,v)
    @test M._materialized_curved_consistent(m,v)
    # Bounding box through the encoding.
    xmin,ymin,zmin,xmax,ymax,zmax=A.model.get_bounding_box(3,v)
    @test xmin≈-1.5 && xmax≈1.5 && zmin≈0.0 && zmax≈2.0
    # Arbitrary axis: consistency and frame hold for a non-z direction.
    occ2=_occ_setup()
    v2=occ2.add_cylinder(0.0,0.0,0.0, 1.0,1.0,1.0, 0.5)
    occ2.synchronize()
    m2=_current_model()
    @test M._materialized_curved_consistent(m2,v2)
    xmin,ymin,zmin,xmax,ymax,zmax=A.model.get_bounding_box(3,v2)
    d=0.5*sqrt(2.0/3.0)
    @test xmin≈-d atol=1e-12
    @test xmax≈1.0+d atol=1e-12
    # Guards.
    occ3=_occ_setup()
    @test occursin("radius must be positive",
        _errtext(()->occ3.add_cylinder(0.0,0.0,0.0,0.0,0.0,1.0,-1.0)))
    @test occursin("radius must be positive",
        _errtext(()->occ3.add_cylinder(0.0,0.0,0.0,0.0,0.0,1.0,0.0)))
    @test occursin("Cannot build cylinder of zero height",
        _errtext(()->occ3.add_cylinder(0.0,0.0,0.0,0.0,0.0,0.0,1.0)))
    @test occursin("Cannot build cylinder with angle <= 0 or angle > 2*Pi",
        _errtext(()->occ3.add_cylinder(0.0,0.0,0.0,0.0,0.0,1.0,1.0;angle=0.0)))
    @test occursin("Cannot build cylinder with angle <= 0 or angle > 2*Pi",
        _errtext(()->occ3.add_cylinder(0.0,0.0,0.0,0.0,0.0,1.0,1.0;angle=3π)))
    # A half-sector cylinder materializes: lateral wall, two meridian
    # planes, two cap planes; mass is the half-cylinder volume.
    pc=occ3.add_cylinder(0.0,0.0,0.0,0.0,0.0,1.0,1.0;angle=π)
    @test pc==1
    occ3.synchronize()
    m3=_current_model()
    @test occ3.get_mass(3,pc)≈π/2
    @test count(==(:plane),values(m3.surface_types))==4
    @test count(==(:cylinder),values(m3.surface_types))==1
    @test occursin("coordinates must be finite",
        _errtext(()->occ3.add_cylinder(NaN,0.0,0.0,0.0,0.0,1.0,1.0)))
    @test occursin("coordinates must be finite",
        _errtext(()->occ3.add_cylinder(0.0,0.0,0.0,0.0,0.0,Inf,1.0)))
    # Failed calls consume no volume tag: the next valid call gets 2.
    @test occ3.add_cylinder(0.0,0.0,0.0,0.0,0.0,1.0,1.0)==2
end

@testset "occ.add_sphere" begin
    occ=_occ_setup()
    v=occ.add_sphere(1.0,2.0,3.0,2.5)
    @test v==1
    occ.synchronize()
    m=_current_model()
    @test length(m.points)==2
    @test length(m.curves)==3
    @test length(m.surfaces)==1
    @test length(m.surface_loops)==1
    @test m.volumes==Dict(1=>[1])
    @test m.surface_types[1]==:sphere
    kinds=[m.curve_types[t] for t in 1:3]
    @test kinds==[:degenerate,:circle,:degenerate]
    # Wire [-degN,-meridian,+degS,+meridian]; shell [face].
    @test m.loops[m.surfaces[1][1]]==[-1,-2,3,2]
    @test m.surface_loops[1]==[1]
    g=m.surface_geometry[1]
    @test g.occ===:sphere
    @test collect(g.center)≈[1.0,2.0,3.0]
    @test g.radius==2.5
    # Equator and pole evaluation.
    @test M.model_value(m,2,1,[0.0,0.0])≈[3.5,2.0,3.0]
    @test M.model_value(m,2,1,[0.0,π/2])≈[1.0,2.0,5.5] atol=1e-12
    @test M.model_value(m,2,1,[0.0,-π/2])≈[1.0,2.0,0.5] atol=1e-12
    lo,hi=M.model_parametrization_bounds(m,2,1)
    @test lo≈[0.0,-π/2] && hi≈[2π,π/2]
    @test A.model.get_normal(1,[0.0,0.0])≈[1.0,0.0,0.0] atol=1e-12
    # Degenerate edges evaluate to their pole via the bound pcurve.
    for e in (1,3)
        p=M.model_value(m,1,e,[0.0])
        @test hypot(p[1]-1.0,p[2]-2.0,abs(p[3]-3.0)-2.5)<1e-12
    end
    @test haskey(m.spheres,v)
    @test m.spheres[v].radius==2.5
    @test M._implicit_volume_surface(m,v)
    @test M._materialized_curved_consistent(m,v)
    xmin,ymin,zmin,xmax,ymax,zmax=A.model.get_bounding_box(3,v)
    @test collect((xmin,ymin,zmin,xmax,ymax,zmax))≈
          [-1.5,-0.5,0.5,3.5,4.5,5.5]
    # Guards.
    occ2=_occ_setup()
    @test occursin("Sphere radius should be positive",
        _errtext(()->occ2.add_sphere(0.0,0.0,0.0,-1.0)))
    @test occursin("Sphere radius should be positive",
        _errtext(()->occ2.add_sphere(0.0,0.0,0.0,0.0)))
    @test occursin("Cannot build sphere with angle <= 0 or angle > 2*Pi",
        _errtext(()->occ2.add_sphere(0.0,0.0,0.0,1.0;angle3=0.0)))
    @test occursin("Cannot build sphere with angle <= 0 or angle > 2*Pi",
        _errtext(()->occ2.add_sphere(0.0,0.0,0.0,1.0;angle3=3π)))
    # Partial sweeps materialize: a polar band (hemisphere cap) and an
    # azimuthal wedge each measure half the unit sphere's volume.
    ps1=occ2.add_sphere(0.0,0.0,0.0,1.0;angle1=0.0)
    @test ps1==1
    ps2=occ2.add_sphere(0.0,0.0,0.0,1.0;angle3=π)
    @test ps2==2
    occ2.synchronize()
    @test occ2.get_mass(3,ps1)≈2π/3
    @test occ2.get_mass(3,ps2)≈2π/3
    @test occursin("coordinates must be finite",
        _errtext(()->occ2.add_sphere(NaN,0.0,0.0,1.0)))
    @test occ2.add_sphere(0.0,0.0,0.0,1.0)==3
end

@testset "occ.add_cone" begin
    occ=_occ_setup()
    # Frustum: lateral + both caps.
    v=occ.add_cone(0.0,0.0,10.0, 0.0,0.0,3.0, 2.0,1.0)
    occ.synchronize()
    m=_current_model()
    @test length(m.surfaces)==3
    @test m.surface_types[1]==:cone
    @test m.surface_types[2]==:plane
    @test m.surface_types[3]==:plane
    @test m.surface_loops[1]==[1,2,-3]
    g=m.surface_geometry[1]
    @test g.occ===:cone && g.r1==2.0 && g.r2==1.0 && g.height==3.0
    len=sqrt(9.0+1.0)
    @test M.model_value(m,2,1,[0.0,0.0])≈[2.0,0.0,10.0]
    @test M.model_value(m,2,1,[0.0,len])≈[1.0,0.0,13.0] atol=1e-12
    @test M._materialized_curved_consistent(m,v)
    # Apex at top (r2=0): degenerate top edge, no top cap.
    va=occ.add_cone(0.0,0.0,0.0, 0.0,0.0,3.0, 2.0,0.0)
    occ.synchronize()
    m=_current_model()
    @test m.surface_loops[2]==[4,-6]
    degs=[t for (t,k) in m.curve_types if k==:degenerate]
    @test length(degs)==1
    apex=M.model_value(m,1,only(degs),[0.0])
    @test apex≈[0.0,0.0,3.0] atol=1e-12
    @test M._materialized_curved_consistent(m,va)
    # Apex at bottom (r1=0): degenerate bottom edge, no bottom cap.
    vb=occ.add_cone(0.0,0.0,0.0, 0.0,0.0,3.0, 0.0,2.0)
    occ.synchronize()
    m=_current_model()
    loop=m.surface_loops[3]
    @test length(loop)==2 && loop[2]>0
    @test M._materialized_curved_consistent(m,vb)
    # Compact encodings for all three.
    @test length(m.cones)==3
    @test m.cones[va].r2==0.0
    # Guards.
    occ2=_occ_setup()
    @test occursin("Cannot build cone of zero height",
        _errtext(()->occ2.add_cone(0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.5)))
    @test occursin("Cone angle should be positive",
        _errtext(()->occ2.add_cone(0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.5;angle=0.0)))
    @test occursin("radii must be non-negative with at least one positive",
        _errtext(()->occ2.add_cone(0.0,0.0,0.0,0.0,0.0,1.0,-1.0,0.5)))
    @test occursin("radii must be non-negative with at least one positive",
        _errtext(()->occ2.add_cone(0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0)))
    # A half-sector cone materializes and measures the half-frustum volume.
    pcn=occ2.add_cone(0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.5;angle=π)
    @test pcn==1
    @test occ2.get_mass(3,pcn)≈(π/3)*(1.0+0.5+0.25)/2
    @test occursin("coordinates must be finite",
        _errtext(()->occ2.add_cone(0.0,0.0,0.0,0.0,NaN,1.0,1.0,0.5)))
    @test occ2.add_cone(0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.5)==2
end

@testset "occ.add_torus" begin
    occ=_occ_setup()
    v=occ.add_torus(0.0,0.0,0.0,5.0,1.0)
    occ.synchronize()
    m=_current_model()
    @test length(m.points)==1
    @test length(m.curves)==2
    @test length(m.surfaces)==1
    @test m.surface_types[1]==:torus
    # Rim vertex, equator + meridian circles both closed on it.
    @test m.curves[1]==(1,1)
    @test m.curves[2]==(1,1)
    # Wire [-equator,+meridian,+equator,-meridian]; shell [face].
    @test m.loops[m.surfaces[1][1]]==[-1,2,1,-2]
    @test m.surface_loops[1]==[1]
    g=m.surface_geometry[1]
    @test g.occ===:torus && g.r1==5.0 && g.r2==1.0
    @test M.model_value(m,2,1,[0.0,0.0])≈[6.0,0.0,0.0]
    @test M.model_value(m,2,1,[π,0.0])≈[-6.0,0.0,0.0] atol=1e-12
    @test A.model.get_normal(1,[0.0,0.0])≈[1.0,0.0,0.0] atol=1e-12
    # Partial torus: end/start rim vertices, trimmed equator arc, two
    # meridian circles, torus face + two Plane caps.
    vp=occ.add_torus(0.0,0.0,0.0,5.0,1.0;angle=π/2)
    occ.synchronize()
    m=_current_model()
    @test m.surface_types[2]==:torus
    @test m.surface_types[3]==:plane
    @test m.surface_types[4]==:plane
    @test m.surface_loops[2]==[2,3,-4]
    gp=m.surface_geometry[2]
    @test gp.occ===:torus && gp.angle==π/2
    @test M.model_value(m,2,2,[π/2,0.0])≈[0.0,6.0,0.0] atol=1e-12
    # The equator arc spans [0,angle]; meridians stay full circles.
    @test m.curve_geometry[3].t1==π/2
    # zAxis: the revolve axis rotates the whole construction.
    vz=occ.add_torus(0.0,0.0,0.0,5.0,1.0;zAxis=[1.0,0.0,0.0])
    occ.synchronize()
    m=_current_model()
    gz=m.surface_geometry[5]
    @test collect(gz.axis)≈[1.0,0.0,0.0]
    @test M.model_value(m,2,5,[0.0,0.0])≈[0.0,0.0,6.0] atol=1e-12
    # Guards.
    occ2=_occ_setup()
    @test occursin("Torus radii should be positive",
        _errtext(()->occ2.add_torus(0.0,0.0,0.0,-5.0,1.0)))
    @test occursin("Torus radii should be positive",
        _errtext(()->occ2.add_torus(0.0,0.0,0.0,5.0,0.0)))
    @test occursin("angle must lie in (0, 2*Pi]",
        _errtext(()->occ2.add_torus(0.0,0.0,0.0,5.0,1.0;angle=0.0)))
    @test occursin("angle must lie in (0, 2*Pi]",
        _errtext(()->occ2.add_torus(0.0,0.0,0.0,5.0,1.0;angle=3π)))
    @test occursin("zAxis must be a 3-vector",
        _errtext(()->occ2.add_torus(0.0,0.0,0.0,5.0,1.0;zAxis=[1.0,0.0])))
    @test occursin("zAxis must have positive length",
        _errtext(()->occ2.add_torus(0.0,0.0,0.0,5.0,1.0;zAxis=[0.0,0.0,0.0])))
    @test occursin("r2 must be finite",
        _errtext(()->occ2.add_torus(0.0,0.0,0.0,5.0,NaN)))
    @test occ2.add_torus(0.0,0.0,0.0,5.0,1.0)==1
end

@testset "occ primitives tags and lifecycle" begin
    occ=_occ_setup()
    # Explicit tags honored; collisions reject against the OCC namespace.
    v=occ.add_cylinder(0.0,0.0,0.0,0.0,0.0,1.0,0.5;tag=7)
    @test v==7
    @test occursin("OCC Volume[7] already exists",
        _errtext(()->occ.add_sphere(0.0,0.0,0.0,1.0;tag=7)))
    occ.synchronize()
    @test A.model.occ.get_max_tag(3)==7
    @test (3,7) in A.model.occ.get_entities(3)
    # Auto tags continue from the explicit maximum.
    @test occ.add_sphere(0.0,0.0,0.0,1.0)==8
    # Repeated synchronize is idempotent.
    occ.synchronize(); occ.synchronize()
    m=_current_model()
    @test length(m.volumes)==2
    @test length(m.surface_loops)==2
    @test sort!(collect(keys(m.cylinders)))==[7]
    @test sort!(collect(keys(m.spheres)))==[8]
    # Recursive removal cleans the volume and every bound sub-entity plus
    # the compact encoding.
    occ.remove([(3,7)],true)
    occ.synchronize()
    m=_current_model()
    @test isempty(m.cylinders)
    @test !haskey(m.volumes,7)
    @test (3,7) ∉ A.model.occ.get_entities(3)
    # The sphere's boundary stays — only the removed volume's graph left.
    @test haskey(m.volumes,8)
    @test haskey(m.spheres,8)
end

@testset "occ primitive volume meshing" begin
    occ=_occ_setup()
    v=occ.add_cylinder(0.0,0.0,0.0, 0.0,0.0,2.0, 1.0)
    vs=occ.add_sphere(10.0,0.0,0.0,2.0)
    occ.synchronize()
    A.mesh.generate(3)
    node_tags,coords=A.mesh.get_nodes(3,v,true)
    @test length(node_tags)>0
    for i in 1:3:3length(node_tags)
        x,y,z=coords[i],coords[i+1],coords[i+2]
        @test hypot(x,y)<=1.0+1e-9
        @test -1e-9<=z<=2.0+1e-9
    end
    element_types,element_tags,_=A.mesh.get_elements(3,v)
    @test element_types==Int32[4]
    @test length(element_tags[1])>0
    node_tags,coords=A.mesh.get_nodes(3,vs,true)
    @test length(node_tags)>0
    for i in 1:3:3length(node_tags)
        x,y,z=coords[i],coords[i+1],coords[i+2]
        @test hypot(x-10.0,y,z)<=2.0+1e-9
    end
end
