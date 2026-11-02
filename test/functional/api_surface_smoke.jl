module ApiSurfaceSmokeChecks
using Test, Tessella
import Sockets

# Public entry points the rest of the suite does not reach:
# model.boolean_union/boolean_intersection, model.add_bezier_surface,
# occ.add_bezier_surface, occ.set_max_tag, mesh.import_stl,
# mesh.set_transfinite_automatic, mesh.set_size_at_parametric_points,
# mesh.set_compound, mesh.add_elements, field.set_numbers/get_numbers,
# field.set_as_boundary_layer, onelab.client_warning/client_define_string/
# client_preprocess, and the top-level rec_dimension_max.
let api=Tessella.API
    api.initialize()
    try
        # Legacy factory boolean operators return a fresh volume tag and
        # remove both operands. (Result queries route through the operand
        # snapshot — see the cross-cutting audit note on explicit volumes.)
        occ=api.model.occ
        a=occ.add_box(0.0,0.0,0.0,1.0,1.0,1.0)
        b=occ.add_box(0.5,0.5,0.5,1.0,1.0,1.0)
        occ.synchronize()
        u=api.model.boolean_union(a,b)
        @test u>0
        @test (3,u) in api.model.get_entities(3)
        @test (3,a) ∉ api.model.get_entities(3) &&
              (3,b) ∉ api.model.get_entities(3)
        @test_throws ArgumentError api.model.boolean_union(0,a)

        api.initialize()
        occ=api.model.occ
        a=occ.add_box(0.0,0.0,0.0,1.0,1.0,1.0)
        b=occ.add_box(0.5,0.5,0.5,1.0,1.0,1.0)
        occ.synchronize()
        i=api.model.boolean_intersection(a,b)
        @test i>0
        @test (3,i) in api.model.get_entities(3)

        # Factory and OCC tensor-product Bezier surfaces over a 2×2 public
        # Point grid; occ.set_max_tag raises the allocation floor.
        api.initialize()
        geo=api.model.geo
        g1=geo.add_point(0.0,0.0,0.0);g2=geo.add_point(1.0,0.0,0.0)
        g3=geo.add_point(0.0,1.0,0.0);g4=geo.add_point(1.0,1.0,0.0)
        geo.synchronize()
        fs=api.model.add_bezier_surface([g1,g2,g3,g4],2)
        @test (2,fs) in api.model.get_entities(2)

        api.initialize()
        occ=api.model.occ
        occ.set_max_tag(0,100)
        p1=occ.add_point(0.0,0.0,0.0)
        @test p1>=101
        p2=occ.add_point(1.0,0.0,0.0)
        p3=occ.add_point(0.0,1.0,0.0)
        p4=occ.add_point(1.0,1.0,0.0)
        os=occ.add_bezier_surface([p1,p2,p3,p4],2)
        @test os>0
        occ.synchronize()
        @test (2,os) in api.model.get_entities(2)

        # mesh.import_stl adds one discrete Surface per retained STL file to
        # the session model at the lowest free tag.
        mktempdir() do dir
            stl=joinpath(dir,"tri.stl")
            open(stl,"w") do io
                print(io,"solid t\nfacet normal 0 0 1\nouter loop\n",
                    "vertex 0 0 0\nvertex 1 0 0\nvertex 0 1 0\nendloop\nendfacet\nendsolid t\n")
            end
            api.initialize()
            api.open(stl)
            first_surfaces=sort!([t for (d,t) in api.model.get_entities(2)])
            @test !isempty(first_surfaces)
            api.mesh.import_stl()
            after=sort!([t for (d,t) in api.model.get_entities(2)])
            @test length(after)==length(first_surfaces)+1
        end

        # Transfinite-automatic on a unit quad with explicit corner sizes:
        # recombining surfaces yields a conforming 5x5 grid; parametric size
        # points insert a pinned node; a curve compound meshes two collinear
        # curves as one. (Corner sizes are required because unsized points —
        # point_size==0.0, which is what geo/occ synchronization writes for
        # every default vertex — currently trip an InexactError inside
        # _transfinite_automatic_surface!; see the audit report.)
        api.initialize()
        geo=api.model.geo
        t1=geo.add_point(0.0,0.0,0.0);t2=geo.add_point(1.0,0.0,0.0)
        t3=geo.add_point(1.0,1.0,0.0);t4=geo.add_point(0.0,1.0,0.0)
        e1=geo.add_line(t1,t2);e2=geo.add_line(t2,t3)
        e3=geo.add_line(t3,t4);e4=geo.add_line(t4,t1)
        s=geo.add_plane_surface([geo.add_curve_loop([e1,e2,e3,e4])])
        geo.synchronize()
        api.mesh.set_size([(0,t1),(0,t2),(0,t3),(0,t4)],0.25)
        api.mesh.set_transfinite_automatic()
        api.mesh.generate(2)
        mesh=api.mesh.get()
        @test api.rec_dimension_max(mesh)==2
        ttags,_,_=api.mesh.get_nodes()
        @test length(ttags)==25
        types,elements,nodes=api.mesh.get_elements()
        @test any(==(3),types)

        api.initialize()
        occ=api.model.occ
        c1=occ.add_point(0.0,0.0,0.0);c2=occ.add_point(1.0,0.0,0.0)
        l1=occ.add_line(c1,c2);occ.synchronize()
        api.mesh.set_size_at_parametric_points(1,l1,[0.5],[0.1])
        api.mesh.generate(1)
        _,xyz,_=api.mesh.get_nodes(1,l1,true,true)
        @test any(x->isapprox(x,0.5;atol=1e-9),xyz[1:3:end])

        api.initialize()
        occ=api.model.occ
        q1=occ.add_point(0.0,0.0,0.0);q2=occ.add_point(1.0,0.0,0.0)
        q3=occ.add_point(2.0,0.0,0.0)
        m1=occ.add_line(q1,q2);m2=occ.add_line(q2,q3)
        occ.synchronize()
        api.mesh.set_compound(1,[m1,m2])
        api.mesh.generate(1)
        tags,_,_=api.mesh.get_nodes()
        @test !isempty(tags)
        @test_throws ArgumentError api.mesh.set_compound(0,[q1,q2])

        # Discrete entity + manual node/element publication.
        api.initialize()
        api.model.add_discrete_entity(2,7)
        api.mesh.add_nodes(2,7,[101,102,103],[0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0])
        api.mesh.add_elements(2,7,[2],[[50]],[[101,102,103]])
        ft,fe,fn=api.mesh.get_elements(2,7)
        @test ft==[2] && fe==[[50]] && fn==[[101,102,103]]

        # Field option round-trip and boundary-layer flag lifecycle.
        api.initialize()
        ft=api.mesh.field.add("MathEval")
        api.mesh.field.set_numbers(ft,"MyList",[1.5,2,3])
        @test api.mesh.field.get_numbers(ft,"MyList")==[1.5,2.0,3.0]
        @test api.mesh.field.get_numbers(ft,"Absent")==Float64[]
        api.mesh.field.set_as_boundary_layer(ft)
        api.mesh.field.remove(ft)
        @test_throws ArgumentError api.mesh.field.set_as_boundary_layer(ft)
    finally
        api.finalize()
    end
end

# ONELAB client entry points not exercised by the interface suite:
# client_warning (fire-and-forget WARNING frame), client_define_string and
# client_preprocess (OLPARSE round-trip) — against a small in-process
# fixture socket speaking the documented framing.
let api=Tessella.API
    function read_exact(sock,n)
        buf=Vector{UInt8}(undef,n);read!(sock,buf);return buf
    end
    function recv_frame(sock)
        hdr=read_exact(sock,8)
        type=reinterpret(Int32,hdr[1:4])[1]
        len=reinterpret(Int32,hdr[5:8])[1]
        return type,String(read_exact(sock,Int(len)))
    end
    function send_frame(sock,type,payload="")
        write(sock,Int32(type));write(sock,Int32(ncodeunits(payload)));write(sock,payload)
        flush(sock)
    end
    function fixture(sock)
        commands=String[]
        while true
            type,msg=recv_frame(sock)
            if type==1
                continue                        # START
            elseif type==2
                break                           # STOP
            elseif type==24                     # PARAMETER_QUERY
                send_frame(sock,23,msg)         # answer with the server copy
            elseif type==28                     # OLPARSE
                send_frame(sock,28,"true")
            elseif type==23 || type==32         # PARAMETER / PARAMETER_UPDATE
                continue                        # define/update publication
            elseif type in (10,11,12,13)
                push!(commands,"msg$(type):"*msg)
            end
        end
        return commands
    end
    server=Sockets.listen(Sockets.localhost,0)
    _,port=Sockets.getsockname(server)
    done=Channel{Vector{String}}(1)
    listener=Threads.@spawn begin
        sock=Sockets.accept(server)
        put!(done,fixture(sock))
    end
    try
        c=api.onelab.client_connect("127.0.0.1",port;name="coverage client")
        api.onelab.client_warning(c,"watch out")
        values=api.onelab.client_define_string(c,"s/mode";values="check")
        @test values==["check"]
        @test api.onelab.client_preprocess(c,"pre","fixture.ol")==true
        api.onelab.client_disconnect(c)
        commands=take!(done)
        @test "msg11:watch out" in commands
    finally
        close(server)
        istaskfailed(listener) || wait(listener)
    end
end
end
