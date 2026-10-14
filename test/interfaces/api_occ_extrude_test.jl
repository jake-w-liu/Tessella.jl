# `model.occ.extrude` — `OCC_Internals::extrude` translation sweeps
# (`BRepPrimAPI_MakePrism`) on the pending-shape kernel. Bounded
# coverage: points sweep to line edges, line edges to ruled
# parallelogram faces, line-bounded plane faces (with holes) to prism
# solids, and a full-circle disk along its normal to the compact
# cylinder. Output ordering follows upstream's `getReturnedShapes`
# reorder — per input `[top (dim-1), body (dim), laterals …]` —
# grouped by ascending dimension; dim-0 inputs emit just the swept
# edge. Swept vertices/lateral edges/top edges deduplicate across
# inputs like upstream's single compound. `numElements`/`heights`/
# `recombine` transport onto the Julia-native Layers machinery (the
# upstream `ExtrudeMesh` attachment); curved profiles, solid or wire
# inputs and in-plane directions reject loudly in phase 1 so
# rejections leave the registry untouched. Face records bake their
# intended outward normal into the stored `axis` (the solid measure
# path reads `σ·axis`), so `get_mass` verifies orientation end-to-end.

using Test
using Tessella
const A=Tessella.API
const M=Tessella.Model

function _occ_setup()
    A.initialize()
    return A.model.occ
end

_current_model()=A._with_model(identity)
_kernel()=A._with_model(m->m.occ_kernel)

_errtext(f)=(try;f();catch err;return sprint(showerror,err);end;"")

function _rect_face(occ,x0,y0,dx,dy)
    p=[occ.add_point(x,y,0.0)
       for (x,y) in [(x0,y0),(x0+dx,y0),(x0+dx,y0+dy),(x0,y0+dy)]]
    es=[occ.add_line(p[1],p[2]),occ.add_line(p[2],p[3]),
        occ.add_line(p[3],p[4]),occ.add_line(p[4],p[1])]
    return p,es,occ.add_plane_surface([occ.add_wire(es)])
end

@testset "occ.extrude point and line inputs" begin
    occ=_occ_setup()
    p=occ.add_point(0.0,0.0,0.0)
    out=occ.extrude([(0,p)],1.0,2.0,3.0)
    @test length(out)==1
    @test out[1][1]==1
    edge_tag=out[1][2]
    kernel=_kernel()
    eid=kernel.bindings[1,edge_tag]
    eg=kernel.shapes[eid].geometry.geometry
    @test eg.occ===:line
    @test collect(eg.origin)≈collect((0.0,0.0,0.0))
    @test collect(eg.dir)≈collect((1.0,2.0,3.0)./sqrt(14.0))
    @test eg.t1-eg.t0≈sqrt(14.0)
    # The swept top vertex is a fresh dim-0 entity translated by d.
    a_id,b_id=kernel.shapes[eid].children
    @test collect(kernel.shapes[b_id].geometry.point)≈collect((1.0,2.0,3.0))

    # Line extrusion reports (top edge, face, two laterals) in
    # profile order; the face is the ruled parallelogram.
    q1=occ.add_point(0.0,0.0,0.0)
    q2=occ.add_point(2.0,0.0,0.0)
    e=occ.add_line(q1,q2)
    out=occ.extrude([(1,e)],0.0,0.0,4.0)
    @test [d for (d,_) in out]==[1,2,1,1]
    top_tag,face_tag,la_tag,lb_tag=(out[1][2],out[2][2],
                                    out[3][2],out[4][2])
    kernel=_kernel()
    top=kernel.shapes[kernel.bindings[1,top_tag]]
    @test collect(top.geometry.geometry.origin)≈collect((0.0,0.0,4.0))
    @test collect(top.geometry.geometry.dir)≈collect((1.0,0.0,0.0))
    # Laterals connect source to top vertices.
    la=kernel.shapes[kernel.bindings[1,la_tag]]
    @test collect(la.geometry.geometry.origin)≈collect((0.0,0.0,0.0))
    @test collect(la.geometry.geometry.dir)≈collect((0.0,0.0,1.0))
    lb=kernel.shapes[kernel.bindings[1,lb_tag]]
    @test collect(lb.geometry.geometry.origin)≈collect((2.0,0.0,0.0))
    # The swept face is a plane through the four corners whose axis
    # points to the exterior of the profile direction (here +y or -y
    # depending on traversal: the quad a→b→b'→a'→a yields -y).
    fs=kernel.shapes[kernel.bindings[2,face_tag]]
    @test fs.geometry.kind===:plane
    @test collect(fs.geometry.geometry.axis)≈collect((0.0,-1.0,0.0))
    @test occ.get_mass(2,face_tag)≈8.0
    @test collect(occ.get_center_of_mass(2,face_tag))≈collect((1.0,0.0,2.0))
end

@testset "occ.extrude planar face to prism" begin
    occ=_occ_setup()
    _,_,f=_rect_face(occ,0.0,0.0,3.0,2.0)
    out=occ.extrude([(2,f)],0.0,0.0,5.0)
    @test [d for (d,_) in out]==[2,3,2,2,2,2]
    cap_tag,solid_tag=out[1][2],out[2][2]
    lateral_tags=[t for (d,t) in out[3:end]]
    @test solid_tag==1
    @test length(lateral_tags)==4
    @test occ.get_mass(3,solid_tag)≈30.0
    @test collect(occ.get_center_of_mass(3,solid_tag))≈collect((1.5,1.0,2.5))
    @test collect(occ.get_bounding_box(3,solid_tag))≈
        collect((0.0,0.0,0.0,3.0,2.0,5.0))
    # The translated cap face sits at z=5 with axis +z.
    kernel=_kernel()
    cap=kernel.shapes[kernel.bindings[2,cap_tag]]
    @test cap.geometry.kind===:plane
    @test collect(cap.geometry.geometry.axis)≈collect((0.0,0.0,1.0))
    @test collect(cap.geometry.geometry.center)≈collect((1.5,1.0,5.0))
    # Lateral axes point outward horizontally.
    for t in lateral_tags
        lg=kernel.shapes[kernel.bindings[2,t]].geometry.geometry
        @test abs(lg.axis[3])<1e-12
    end
    # Solid topology: one shell of six faces.
    sid=kernel.bindings[3,solid_tag]
    shell=kernel.shapes[only(kernel.shapes[sid].children)]
    @test length(shell.geometry.uses)==6
end

@testset "occ.extrude face with hole" begin
    occ=_occ_setup()
    po,eo,wo=(let p=[occ.add_point(x,y,0.0)
                    for (x,y) in [(0,0),(4,0),(4,4),(0,4)]]
              e=[occ.add_line(p[1],p[2]),occ.add_line(p[2],p[3]),
                 occ.add_line(p[3],p[4]),occ.add_line(p[4],p[1])]
              (p,e,occ.add_wire(e)) end)
    ph=[occ.add_point(x,y,0.0)
        for (x,y) in [(1.5,1.5),(2.5,1.5),(2.5,2.5),(1.5,2.5)]]
    eh=[occ.add_line(ph[1],ph[2]),occ.add_line(ph[2],ph[3]),
        occ.add_line(ph[3],ph[4]),occ.add_line(ph[4],ph[1])]
    wh=occ.add_wire(eh)
    f=occ.add_plane_surface([wo,wh])
    out=occ.extrude([(2,f)],0.0,0.0,3.0)
    solid_tag=out[2][2]
    # (4·4 − 1·1)·3 = 45 — the hole laterals subtract correctly.
    @test occ.get_mass(3,solid_tag)≈45.0
    @test collect(occ.get_center_of_mass(3,solid_tag))≈collect((2.0,2.0,1.5))
    @test collect(occ.get_bounding_box(3,solid_tag))≈
        collect((0.0,0.0,0.0,4.0,4.0,3.0))
    @test length(out)==2+8
    # Cap keeps the hole: two wires.
    kernel=_kernel()
    cap=kernel.shapes[kernel.bindings[2,out[1][2]]]
    @test length(cap.children)==2

    # Anti-normal sweep mirrors the cap and hole-wall traversals; the
    # solid still reports a positive outward-oriented volume.
    occ2=_occ_setup()
    _,_,f2=_rect_face(occ2,0.0,0.0,4.0,4.0)
    out2=occ2.extrude([(2,f2)],0.0,0.0,-3.0)
    @test occ2.get_mass(3,out2[2][2])≈48.0
    @test collect(occ2.get_bounding_box(3,out2[2][2]))≈
        collect((0.0,0.0,-3.0,4.0,4.0,0.0))
end

@testset "occ.extrude disk to cylinder" begin
    occ=_occ_setup()
    d=occ.add_disk(0.0,0.0,0.0,2.0,2.0)
    out=occ.extrude([(2,d)],0.0,0.0,5.0)
    @test [x for (x,_) in out]==[2,3,2]
    top_tag,solid_tag,lateral_tag=(out[1][2],out[2][2],out[3][2])
    @test occ.get_mass(3,solid_tag)≈pi*4.0*5.0 rtol=1e-10
    @test collect(occ.get_bounding_box(3,solid_tag))≈
        collect((-2.0,-2.0,0.0,2.0,2.0,5.0)) atol=1e-9
    kernel=_kernel()
    @test kernel.shapes[kernel.bindings[2,lateral_tag]].geometry.kind===:cylinder
    top=kernel.shapes[kernel.bindings[2,top_tag]]
    @test top.geometry.kind===:plane
    @test collect(top.geometry.geometry.center)≈collect((0.0,0.0,5.0))
end

@testset "occ.extrude shared topology dedup" begin
    occ=_occ_setup()
    p=occ.add_point(0.0,0.0,0.0)
    q=occ.add_point(1.0,0.0,0.0)
    e=occ.add_line(p,q)
    out=occ.extrude([(0,p),(1,e)],0.0,2.0,0.0)
    # Point output groups first, then the line's (top, face, laterals);
    # the shared generatrix is one edge — the line's lateral at `p`
    # equals the point's swept edge.
    @test out[1][1]==1
    shared=out[1][2]
    @test (1,shared) in out[2:end]
    kernel=_kernel()
    # Only one edge record carries vertex p's generatrix: exactly one
    # duplicate entry in the reported output.
    @test length(unique(out))==length(out)-1
end

@testset "occ.extrude rejection atomicity" begin
    occ=_occ_setup()
    p,es,f=_rect_face(occ,0.0,0.0,2.0,1.0)
    v=occ.add_box(10.0,0.0,0.0,1.0,1.0,1.0)
    c=occ.add_circle(0.0,10.0,0.0,1.0)
    disk=occ.add_disk(0.0,20.0,0.0,1.0,1.0)
    ell=occ.add_disk(0.0,30.0,0.0,2.0,1.0)
    function snapshot()
        kernel=_kernel()
        return (length(kernel.bindings),kernel.next_shape,
                copy(kernel.maxima))
    end
    cases=[
        (()->occ.extrude([(2,999)],0.0,0.0,1.0),
         "Unknown OpenCASCADE entity of dimension 2 with tag 999"),
        (()->occ.extrude([(3,v)],0.0,0.0,1.0),
         "extrusion of OCC volumes is not yet supported"),
        (()->occ.extrude([(-1,1)],0.0,0.0,1.0),
         "Wrong input dimension in extrusion"),
        (()->occ.extrude([(4,f)],0.0,0.0,1.0),
         "OCC dimension must be in -2:3"),
        (()->occ.extrude([(0,p[1])],0.0,0.0,0.0),
         "extrusion direction must have positive length"),
        (()->occ.extrude([(1,es[1])],1.0,0.0,0.0),
         "extrusion direction is parallel to the curve"),
        (()->occ.extrude([(2,f)],1.0,1.0,0.0),
         "coplanar with the surface"),
        (()->occ.extrude([(1,c)],0.0,0.0,1.0),
         "extrusion of OCC curve kind circle is not yet supported"),
        (()->occ.extrude([(2,disk)],1.0,0.0,0.0),
         "extrusion of a disk must follow its surface normal"),
        (()->occ.extrude([(2,ell)],0.0,0.0,1.0),
         "extrusion of elliptical profiles is not yet supported"),
        (()->occ.extrude([(2,f)],0.0,0.0,1.0;heights=[0.5,1.0]),
         "heights requires numElements"),
        (()->occ.extrude([(2,f)],0.0,0.0,1.0;numElements=[0]),
         "numElements entries must be positive integers"),
        (()->occ.extrude([(2,f)],0.0,0.0,1.0;
            numElements=[2,2],heights=[0.5]),
         "one cumulative entry per numElements layer"),
        (()->occ.extrude([(2,f)],0.0,0.0,1.0;
            numElements=[2,2],heights=[0.7,0.5]),
         "strictly increasing"),
        (()->occ.extrude([(2,f)],0.0,0.0,1.0;
            numElements=[2],heights=[0.5]),
         "normalized to 1"),
        # A trailing invalid input rejects the whole call — nothing binds.
        (()->occ.extrude([(2,f),(2,999)],0.0,0.0,1.0),
         "Unknown OpenCASCADE entity of dimension 2 with tag 999"),
        (()->occ.extrude([(2,f),(2,disk)],1.0,0.0,1.0),
         "extrusion of a disk must follow its surface normal"),
    ]
    for (fnc,frag) in cases
        before=snapshot()
        text=_errtext(fnc)
        @test occursin(frag,text)
        @test snapshot()==before
    end
end

@testset "occ.extrude prism mesh" begin
    occ=_occ_setup()
    _,_,f=_rect_face(occ,0.0,0.0,3.0,2.0)
    out=occ.extrude([(2,f)],0.0,0.0,4.0)
    solid_tag=out[2][2]
    occ.synchronize()
    m=_current_model()
    @test haskey(m.volumes,solid_tag)
    @test length(m.surfaces)==6
    A.mesh.generate(3)
    node_tags,coords=A.mesh.get_nodes(3,solid_tag,true)
    @test length(node_tags)>0
    ok=true
    for i in 1:3:3length(node_tags)
        x,y,z=coords[i],coords[i+1],coords[i+2]
        ok &= -1e-9<=x<=3.0+1e-9 && -1e-9<=y<=2.0+1e-9 &&
              -1e-9<=z<=4.0+1e-9
    end
    @test ok
    element_types,element_tags,_=A.mesh.get_elements(3,solid_tag)
    @test element_types==Int32[4]
    @test length(element_tags[1])>0
end

@testset "occ.extrude Layers transport" begin
    # numElements/heights/recombine flow into the native Layers records
    # on every swept entity. A rect face swept 4 units with layers
    # [2,3] at cumulative boundaries 0.4/1.0 grades the volume at
    # z = 0.2k*4 within each layer: {0,.8,1.6} + {2.4,3.2,4.0}.
    occ=_occ_setup()
    _,_,f=_rect_face(occ,0.0,0.0,3.0,2.0)
    out=occ.extrude([(2,f)],0.0,0.0,4.0;numElements=[2,3],
                    heights=[0.4,1.0])
    solid_tag=out[2][2]
    occ.synchronize()
    A.mesh.generate(3)
    _,coords=A.mesh.get_nodes(3,solid_tag,true)
    zs=sort!(unique!(round.([coords[i] for i in 3:3:length(coords)],
                            digits=9)))
    @test zs≈[0.0,0.8,1.6,2.4,3.2,4.0] atol=1e-9
    element_types,_,_=A.mesh.get_elements(3,solid_tag)
    @test element_types==Int32[4]

    # Recombined: the lateral strips recombine to quads, so each
    # triangulated-source column becomes a prism (type 6).
    occ=_occ_setup()
    _,_,f=_rect_face(occ,0.0,0.0,3.0,2.0)
    out=occ.extrude([(2,f)],0.0,0.0,4.0;numElements=[4],
                    recombine=true)
    solid_tag=out[2][2]
    occ.synchronize()
    A.mesh.generate(3)
    element_types,_,_=A.mesh.get_elements(3,solid_tag)
    @test 6 in element_types
    @test !(4 in element_types)

    # A line sweep with recombine emits quads (type 3) on the lateral.
    occ=_occ_setup()
    q1=occ.add_point(0.0,0.0,0.0)
    q2=occ.add_point(2.0,0.0,0.0)
    e=occ.add_line(q1,q2)
    out=occ.extrude([(1,e)],0.0,0.0,4.0;numElements=[3],
                    recombine=true)
    face_tag=out[2][2]
    occ.synchronize()
    A.mesh.generate(2)
    element_types,_,_=A.mesh.get_elements(2,face_tag)
    @test element_types==Int32[3]

    # The compact disk -> cylinder path transports the same records:
    # layers [2,2] at boundaries 0.4/1.0 grade the sweep at
    # z = {0,1,2} + {3.5,5}.
    occ=_occ_setup()
    d=occ.add_disk(0.0,0.0,0.0,2.0,2.0)
    out=occ.extrude([(2,d)],0.0,0.0,5.0;numElements=[2,2],
                    heights=[0.4,1.0])
    solid_tag=out[2][2]
    occ.synchronize()
    A.mesh.generate(3)
    _,coords=A.mesh.get_nodes(3,solid_tag,true)
    zs=sort!(unique!(round.([coords[i] for i in 3:3:length(coords)],
                            digits=9)))
    @test zs≈[0.0,1.0,2.0,3.5,5.0] atol=1e-9
    element_types,_,_=A.mesh.get_elements(3,solid_tag)
    @test element_types==Int32[4]

    # Recombined disk -> prisms (a triangulated source face sweeps
    # tri-column cells; type 6).
    occ=_occ_setup()
    d=occ.add_disk(0.0,0.0,0.0,2.0,2.0)
    out=occ.extrude([(2,d)],0.0,0.0,5.0;numElements=[3],
                    recombine=true)
    solid_tag=out[2][2]
    occ.synchronize()
    A.mesh.generate(3)
    element_types,_,_=A.mesh.get_elements(3,solid_tag)
    @test 6 in element_types
end
