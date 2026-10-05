module PublicRecombineLabelTests
# Public API controls verified against pinned Gmsh 4.15.2 first-order pairing.
using Test,Tessella
const API=Tessella.API
function recombine_square(;seed=true)
    API.initialize();API.option("Mesh.Renumber",0)
    API.option("Mesh.RecombinationAlgorithm",0)
    for (tag,point) in enumerate(((0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.)))
        API.model.add_point(point...;tag)
    end
    for (tag,ends) in enumerate(((1,2),(2,3),(3,4),(4,1)))
        API.model.add_line(ends...;tag)
    end
    API.model.add_curve_loop([1,2,3,4];tag=1);API.model.add_plane_surface([1];tag=1)
    for (_,curve) in API.model.get_entities(1);API.mesh.set_transfinite_curve(curve,2);end
    if seed
        API.mesh.add_nodes(0,1,[111],[0.,0.,0.])
        API.mesh.add_elements_by_type(1,15,[777],[111])
    end
    API.mesh.set_transfinite_surface(1);mesh=API.mesh.generate(2)
    @test mesh isa Tessella.Mesh && API.LAST_MESH[] isa Tessella.Mesh
    @test API.mesh.get_element_types(2,1)==Int32[2]
    return nothing
end
function recombine_snapshot()
    model=API.CURRENT[]
    return (;model,repr=repr(model),cache=API.LAST_MESH[],class=API.LAST_MESH_CLASS[],
        overlay=API.LAST_MESH_HIGH_ORDER[],edges=API.LAST_MESH_EDGES[],faces=API.LAST_MESH_FACES[],
        payload=(API.mesh.get_nodes(),API.mesh.get_elements()),
        counters=(API.NODE_TAG_MAX[],API.ELEMENT_TAG_MAX[]),visibility=copy(API.ELEMENT_VISIBILITY[]))
end
@testset "Native recombine preserves public identities and allocator history" begin
    try
        recombine_square()
        old=API.mesh.get_nodes()[1];mapped=UInt64[10003+7i for i in 0:length(old)-1]
        API.mesh.renumber_nodes(old,mapped)
        elements=API.mesh.get_elements(2,1)[2][1]
        API.mesh.renumber_elements(elements,UInt64[30001+11i for i in 0:length(elements)-1])
        points=Dict(node=>API.mesh.get_node(node)[1] for node in mapped)
        for (entity,node,element,point) in ((41000,20001,60001,[10.,10.,10.]),(41001,20003,90001,[12.,12.,12.]))
            API.model.add_discrete_entity(0,entity);API.mesh.add_nodes(0,entity,[node],point)
            API.mesh.add_elements_by_type(entity,15,[element],[node])
        end
        API.mesh.remove_elements(0,41001,[90001])
        API.model.add_discrete_entity(2,52000)
        API.mesh.add_nodes(2,52000,[20005],[10.,10.,10.])
        API.mesh.add_elements_by_type(52000,2,[65001],[mapped[1],mapped[2],20005])
        referenced=API.mesh.get_elements(2,52000)
        foreign=(API.mesh.get_nodes(0,41000),API.mesh.get_elements(0,41000))
        before=API.mesh.get_nodes();@test API.mesh.get_max_element_tag()==90001
        API.mesh.set_recombine(2,1);API.mesh.recombine()
        types,ids,nodes=API.mesh.get_elements(2,1)
        @test types==Int32[3] && ids==[UInt64[90003]] && length(nodes[1])==4
        @test Set(nodes[1])==Set(mapped)
        @test API.mesh.get_nodes()==before
        @test (API.mesh.get_nodes(0,41000),API.mesh.get_elements(0,41000))==foreign
        @test API.mesh.get_elements(2,52000)==referenced
        @test API.mesh.get_element(65001)==(2,UInt64[mapped[1],mapped[2],20005],2,52000)
        xyz=[points[node] for node in nodes[1]]
        # Literal shoelace area and vertex coordinates are independent of the
        # recombination kernel and accept only the unit Square boundary cycle.
        area=sum(xyz[i][1]*xyz[mod1(i+1,4)][2]-xyz[mod1(i+1,4)][1]*xyz[i][2] for i in 1:4)/2
        @test area==1 && Set(Tuple.(xyz))==Set(((0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.)))
        @test API.mesh.get_barycenters(3,1,false,false)==[.5,.5,0.]
        globalids=vcat(API.mesh.get_elements()[2]...)
        @test allunique(globalids) && API.mesh.get_max_element_tag()==90003
        @test API.mesh.get_element(90003)==(3,nodes[1],2,1)
    finally
        API.finalize()
    end
    # Pinned primary leaves MSH9 triangles in place. Preserve the actual edited
    # P2 map and labels; never recombine its P1 skeleton as if it were MSH2.
    try
        recombine_square();API.mesh.set_order(2)
        old=API.mesh.get_nodes()[1];API.mesh.renumber_nodes(old,UInt64[10003+7i for i in 0:length(old)-1])
        interior=API.mesh.get_nodes(2,1,false,false)[1]
        @test length(interior)==1
        API.mesh.set_node(only(interior),[.45,.55,.01],Float64[])
        @test API.mesh.get_node(only(interior))[1]==[.45,.55,.01]
        overlay=API.LAST_MESH_HIGH_ORDER[]
        before=(API.mesh.get_nodes(),API.mesh.get_elements())
        API.mesh.set_recombine(2,1);API.mesh.recombine()
        @test (API.mesh.get_nodes(),API.mesh.get_elements())==before
        @test API.LAST_MESH_HIGH_ORDER[]===overlay
    finally
        API.finalize()
    end
    # A real exhausted element-label history rejects before changing model,
    # consumed triangles, support caches or counters.
    try
        recombine_square()
        API.model.add_discrete_entity(0,41000);API.mesh.add_nodes(0,41000,[50001],[10.,10.,10.])
        API.mesh.add_elements_by_type(41000,15,[typemax(Int32)],[50001])
        API.mesh.set_recombine(2,1);API.mesh.create_edges();API.mesh.create_faces()
        before=recombine_snapshot()
        @test_throws r"recombined element tags exceed Int32" API.mesh.recombine()
        @test API.CURRENT[]===before.model && repr(API.CURRENT[])==before.repr
        @test API.LAST_MESH[]===before.cache && API.LAST_MESH_CLASS[]===before.class
        @test API.LAST_MESH_HIGH_ORDER[]===before.overlay
        @test API.LAST_MESH_EDGES[]===before.edges && API.LAST_MESH_FACES[]===before.faces
        @test (API.mesh.get_nodes(),API.mesh.get_elements())==before.payload
        @test (API.NODE_TAG_MAX[],API.ELEMENT_TAG_MAX[])==before.counters
        @test API.ELEMENT_VISIBILITY[]==before.visibility
    finally
        API.finalize()
    end
end

# Real native CAD embedding controls retain explicitly Point-owned vertices.
# The discrete declared-boundary diagonal is separately a structural control.
# The pinned filter depends on actual line connectivity; an empty curve entity
# excludes no pair.
@testset "Pinned first-pass greedy admission and actual Curve-edge exclusions" begin
    for (exclusion,angle,expected_type,maximum) in
            ((:none,.5,2,90002),(:none,1.0,2,90002),(:none,1.01,3,90003),
             (:embedded,45.,2,90001),(:boundary,45.,2,90001),
             (:empty_embedded,45.,3,90003))
        try
            API.initialize();API.option("Mesh.Renumber",0)
            API.option("Mesh.RecombinationAlgorithm",0)
            tags=UInt64[10003,10010,10017,10024]
            coordinates=[0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.]
            if exclusion in (:embedded,:empty_embedded)
                for (tag,point) in enumerate(((0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.)))
                    API.model.add_point(point...;tag)
                    API.mesh.add_nodes(0,tag,[tags[tag]],collect(point))
                end
                for (tag,ends) in enumerate(((1,2),(2,3),(3,4),(4,1)))
                    API.model.add_line(ends...;tag)
                end
                API.model.add_line(1,3;tag=11)
                API.model.add_curve_loop([1,2,3,4];tag=1)
                API.model.add_plane_surface([1];tag=1)
            else
                exclusion===:none || API.model.add_discrete_entity(1,11)
                API.model.add_discrete_entity(2,1,exclusion===:boundary ? [11] : Int[])
                API.mesh.add_nodes(2,1,tags,coordinates)
            end
            triangles=UInt64[10003,10010,10017,10003,10017,10024]
            API.mesh.add_elements_by_type(1,2,[30001,30012],triangles)
            if exclusion in (:embedded,:boundary)
                API.mesh.add_elements_by_type(11,1,[40001],[10003,10017])
            end
            exclusion in (:embedded,:empty_embedded) && API.mesh.embed(1,[11],2,1)
            API.model.add_discrete_entity(0,41001)
            API.mesh.add_nodes(0,41001,[20003],[12.,12.,12.])
            API.mesh.add_elements_by_type(41001,15,[90001],[20003])
            API.mesh.remove_elements(0,41001,[90001])
            API.mesh.set_recombine(2,1,angle)
            API.mesh.recombine()
            types,elements,nodes=API.mesh.get_elements(2,1)
            @test types==Int32[expected_type]
            @test API.mesh.get_max_element_tag()==maximum
            @test all(API.mesh.get_node(tags[i])[1]==coordinates[3i-2:3i] for i in 1:4)
            if expected_type==2
                @test elements==[UInt64[30001,30012]] && nodes==[triangles]
            else
                @test elements==[UInt64[90003]] && Set(only(nodes))==Set(tags)
            end
        finally
            API.finalize()
        end
    end
    # The strict admission is a greedy branch contract in pinned Gmsh.
    # Successful Blossom still forms the square at angle0.5.
    try
        recombine_square();API.option("Mesh.RecombinationAlgorithm",1)
        API.mesh.set_recombine(2,1,.5);API.mesh.recombine()
        @test API.mesh.get_element_types(2,1)==Int32[3]
    finally
        API.finalize()
    end
end

@testset "Virgin native node-first sparse identity uses real public plan" begin
    for order in (1,2)
        try
            recombine_square(;seed=false)
            order==2 && API.mesh.set_order(2)
            old,coordinates,_=API.mesh.get_nodes()
            _,elements,_=API.mesh.get_elements()
            element_ids=vcat(elements...)
            mapped=UInt64[10003+7i for i in 0:length(old)-1]
            API.mesh.renumber_nodes(old,mapped)
            @test API.LAST_MESH[] isa Tessella.Mesh
            @test API.mesh.get_nodes()[1:2]==(mapped,coordinates)
            @test Set(vcat(API.mesh.get_elements()[2]...))==Set(element_ids)
            overlay=API.LAST_MESH_HIGH_ORDER[]
            actual=(API.mesh.get_nodes(),API.mesh.get_elements())
            API.mesh.set_recombine(2,1);API.mesh.recombine()
            if order==1
                @test API.mesh.get_element_types(2,1)==Int32[3]
                @test Set(only(API.mesh.get_elements(2,1)[3]))==Set(mapped)
            else
                @test (API.mesh.get_nodes(),API.mesh.get_elements())==actual
                @test API.LAST_MESH_HIGH_ORDER[]===overlay
            end
        finally
            API.finalize()
        end
    end
end

end
