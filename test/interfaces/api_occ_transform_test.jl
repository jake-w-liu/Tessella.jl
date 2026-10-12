# `model.occ` transform endpoints: `translate`, `rotate`, `dilate`, `mirror`
# (`symmetrize` alias), `affine_transform` and `copy` — the shared
# `OCC_Internals::_transform`/`_multiBind` machinery rebound through the
# kernel: the transformed top-level shape keeps its tag while every
# re-created sub-entity takes a fresh `getMaxTag+1` tag, ownership-protected
# recursive unbind detaches shared sub-shapes, and `synchronize`
# removes-then-reimports through `kernel.removed`. Record families use the
# verified `_transform_occ_*` evaluators — circles/ellipses keep handedness
# under reflections, spheres/tori gate on isotropy, and unrepresentable
# transforms reject atomically before any bind.
# Layout mirrors `OCC_Internals::_transform`/`_bind`/`_unbind`/`_multiBind` on
# Gmsh 4.15.2.

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

@testset "occ.translate geometry and tag churn" begin
    occ=_occ_setup()
    v=occ.add_box(0.0,0.0,0.0, 2.0,1.0,3.0)
    p1=occ.add_point(10.0,0.0,0.0); p2=occ.add_point(20.0,0.0,0.0)
    l=occ.add_line(p1,p2)
    occ.synchronize()
    m=_current_model()
    @test sort(collect(keys(m.volumes)))==[1]
    @test sort(collect(keys(m.curves)))==collect(1:13)
    occ.translate([(3,v)],5.0,0.0,0.0)
    occ.synchronize()
    bb=A.model.get_bounding_box(3,v)
    @test collect(bb)≈[5.0,0.0,0.0,7.0,1.0,3.0] atol=1e-12
    # The standalone line keeps tag 13; the re-created box curves take fresh
    # tags after it — upstream `_bind(out, tag, true)` getMaxTag+1 dynamics.
    @test sort(collect(keys(m.curves)))==[13;collect(14:25)]
    @test m.curves[l]==(9,10)
    # Faces recycled their freed tags 1–6 but carry translated records.
    @test sort(collect(keys(m.surfaces)))==collect(1:6)
    # Repeated synchronize is idempotent on the rebound graph.
    occ.synchronize()
    @test sort(collect(keys(m.curves)))==[13;collect(14:25)]
    @test collect(A.model.get_bounding_box(3,v))≈[5.0,0.0,0.0,7.0,1.0,3.0] atol=1e-12
end

@testset "occ.translate shared sub-shape detach" begin
    occ=_occ_setup()
    p1=occ.add_point(0.0,0.0,0.0)
    p2=occ.add_point(1.0,0.0,0.0)
    p3=occ.add_point(1.0,1.0,0.0)
    c1=occ.add_line(p1,p2); c2=occ.add_line(p2,p3)
    occ.synchronize()
    m=_current_model()
    occ.translate([(1,c1)],0.0,0.0,5.0)
    occ.synchronize()
    # The shared vertex stays bound for the untransformed curve — the
    # ownership-protected unbind refuses it, so Curve[1] detaches onto two
    # fresh points exactly like upstream's `_unbind` containment check.
    @test collect(m.points[2])≈[1.0,0.0,0.0]
    @test m.curves[c2]==(2,3)
    @test m.curves[c1]==(4,5)
    @test collect(m.points[4])≈[0.0,0.0,5.0]
    @test collect(m.points[5])≈[1.0,0.0,5.0]
end

@testset "occ.rotate rigid motion" begin
    occ=_occ_setup()
    v=occ.add_box(0.0,0.0,0.0, 4.0,1.0,2.0)
    occ.synchronize()
    occ.rotate([(3,v)],0.0,0.0,0.0, 0.0,0.0,1.0, pi/2)
    occ.synchronize()
    @test collect(A.model.get_bounding_box(3,v))≈
        [-1.0,0.0,0.0,0.0,4.0,2.0] atol=1e-12
    # Cylinder axis rotates with the record.
    occ2=_occ_setup()
    c=occ2.add_cylinder(0.0,0.0,0.0, 0.0,0.0,2.0, 1.0)
    occ2.synchronize()
    occ2.rotate([(3,c)],0.0,0.0,0.0, 1.0,0.0,0.0, pi/2)
    occ2.synchronize()
    m=_current_model()
    @test collect(m.cylinders[c].axis)≈[0.0,-2.0,0.0] atol=1e-12
    @test collect(A.model.get_bounding_box(3,c))≈
        [-1.0,-2.0,-1.0,1.0,0.0,1.0] atol=1e-12
end

@testset "occ.dilate and representability" begin
    occ=_occ_setup()
    s=occ.add_sphere(1.0,0.0,0.0, 0.5)
    occ.synchronize()
    m=_current_model()
    occ.dilate([(3,s)],0.0,0.0,0.0, 2.0,2.0,2.0)
    occ.synchronize()
    @test collect(m.spheres[s].center)≈[2.0,0.0,0.0]
    @test m.spheres[s].radius≈1.0
    # Anisotropic scale warps the sphere — rejects before any bind.
    err=_errtext(()->occ.dilate([(3,s)],0.0,0.0,0.0, 2.0,1.0,1.0))
    @test occursin("not representable",err)
    occ.synchronize()
    @test collect(m.spheres[s].center)≈[2.0,0.0,0.0]
    @test m.spheres[s].radius≈1.0
end

@testset "occ.mirror reflection" begin
    occ=_occ_setup()
    v=occ.add_box(0.0,0.0,0.0, 1.0,1.0,1.0)
    occ.synchronize()
    occ.mirror([(3,v)],1.0,0.0,0.0,0.0)
    occ.synchronize()
    @test collect(A.model.get_bounding_box(3,v))≈
        [-1.0,0.0,0.0,0.0,1.0,1.0] atol=1e-12
    # Degenerate pole edges track the reflected sphere through the embedded
    # face record.
    occ2=_occ_setup()
    s=occ2.add_sphere(1.0,2.0,3.0,2.0)
    occ2.synchronize()
    occ2.symmetrize([(3,s)],0.0,0.0,1.0,0.0)
    occ2.synchronize()
    m=_current_model()
    @test collect(m.spheres[s].center)≈[1.0,2.0,-3.0]
    degs=sort([t for (t,k) in m.curve_types if k===:degenerate])
    @test length(degs)==2
    poles=sort([M.model_value(m,1,t,[0.0]) for t in degs],by=x->x[3])
    @test poles[1]≈[1.0,2.0,-5.0] atol=1e-12
    @test poles[2]≈[1.0,2.0,-1.0] atol=1e-12
end

@testset "occ.affine_transform" begin
    occ=_occ_setup()
    v=occ.add_cylinder(0.0,0.0,0.0, 0.0,0.0,3.0, 1.0)
    occ.synchronize()
    occ.affine_transform([(3,v)],Float64[1,0,0,3, 0,1,0,4, 0,0,1,5])
    occ.synchronize()
    @test collect(A.model.get_bounding_box(3,v))≈
        [2.0,3.0,5.0,4.0,5.0,8.0] atol=1e-12
    # A shear is still a parallelepiped — the plane records absorb it.
    occ2=_occ_setup()
    b=occ2.add_box(0.0,0.0,0.0, 1.0,1.0,1.0)
    occ2.synchronize()
    occ2.affine_transform([(3,b)],Float64[1,0.5,0,0, 0,1,0,0, 0,0,1,0])
    occ2.synchronize()
    bb=A.model.get_bounding_box(3,b)
    @test collect(bb)≈[0.0,0.0,0.0,1.5,1.0,1.0] atol=1e-12
    n=A.model.get_normal(6,[0.5,0.5])
    @test n[1]>0.8 && n[2]<-0.4   # tilted +x face normal
end

@testset "occ.copy" begin
    occ=_occ_setup()
    v=occ.add_cylinder(0.0,0.0,0.0, 0.0,0.0,3.0, 1.0)
    occ.synchronize()
    out=occ.copy([(3,v)])
    @test out==[(3,2)]
    occ.synchronize()
    m=_current_model()
    @test sort(collect(keys(m.volumes)))==[1,2]
    # The compact encoding duplicates with the copied solid.
    @test m.cylinders[1]==m.cylinders[2]
    @test collect(A.model.get_bounding_box(3,2))≈
        collect(A.model.get_bounding_box(3,1))
    # Copying a wire binds its edges at the highest present dimension —
    # `_multiBind` never returns the wire itself.
    p1=occ.add_point(0.0,0.0,10.0); p2=occ.add_point(1.0,0.0,10.0)
    p3=occ.add_point(1.0,1.0,10.0)
    e1=occ.add_line(p1,p2); e2=occ.add_line(p2,p3)
    w=occ.add_wire([e1,e2])
    out2=occ.copy([(-1,w)])
    @test all(d==1 for (d,_) in out2)
    @test length(out2)==2
    # A point copy lands under a fresh tag.
    out3=occ.copy([(0,p1)])
    @test length(out3)==1 && out3[1][1]==0 && out3[1][2]!=p1
    occ.synchronize()
    @test collect(m.points[out3[1][2]])≈collect(m.points[p1])
    # The same entity listed twice produces one copy (compound dedup).
    out4=occ.copy([(3,v),(3,v)])
    @test length(out4)==1
end

@testset "occ transform guards and atomicity" begin
    occ=_occ_setup()
    v=occ.add_box(0.0,0.0,0.0, 1.0,1.0,1.0)
    occ.synchronize()
    m=_current_model()
    before=deepcopy(m.points)
    # Unknown entity.
    @test occursin("Unknown OpenCASCADE entity",
                   _errtext(()->occ.translate([(3,99)],1.0,0.0,0.0)))
    # Duplicate input.
    @test occursin("listed twice",
                   _errtext(()->occ.translate([(3,v),(3,v)],1.0,0.0,0.0)))
    # Non-finite delta.
    @test occursin("finite",
                   _errtext(()->occ.translate([(3,v)],NaN,0.0,0.0)))
    # Zero rotation axis.
    @test occursin("nonzero",
                   _errtext(()->occ.rotate([(3,v)],0.,0.,0., 0.,0.,0., 1.0)))
    occ.synchronize()
    @test m.points==before
    @test sort(collect(keys(m.volumes)))==[1]
end

@testset "occ transform lifecycle" begin
    occ=_occ_setup()
    v=occ.add_cylinder(0.0,0.0,0.0, 0.0,0.0,2.0, 1.0)
    occ.synchronize()
    m=_current_model()
    occ.translate([(3,v)],1.0,0.0,0.0)
    occ.synchronize()
    # One loop, rebound content — the dead shell's surface_loop is pruned.
    @test length(m.surface_loops)==1
    # The transformed volume still removes cleanly.
    occ.remove([(3,v)],true)
    occ.synchronize()
    @test isempty(m.volumes)
    @test isempty(m.cylinders)
    @test isempty(m.surface_loops)
end

@testset "occ transform volume meshing" begin
    occ=_occ_setup()
    v=occ.add_cylinder(0.0,0.0,0.0, 0.0,0.0,3.0, 1.0)
    occ.synchronize()
    occ.translate([(3,v)],3.0,4.0,5.0)
    occ.synchronize()
    A.mesh.generate(3)
    node_tags,coords=A.mesh.get_nodes(3,v,true)
    @test length(node_tags)>0
    for i in 1:3:3length(node_tags)
        x,y,z=coords[i],coords[i+1],coords[i+2]
        @test (x-3.0)^2+(y-4.0)^2<=1.0+1e-9
        @test 5.0-1e-9<=z<=8.0+1e-9
    end
end
