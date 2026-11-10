module NativeLookupLifecycleProbeV2
using Tessella,Test,SHA,TOML,Serialization
const A=Tessella.API
const M=Tessella.Model
include(joinpath(@__DIR__,"native_lookup_lifecycle_fixtures_v1.jl"))
bits(xs)=[string(reinterpret(UInt64,Float64(x));base=16,pad=16) for x in xs]
digest(value)=let stream=IOBuffer();serialize(stream,value);bytes2hex(sha256(take!(stream)));end

function bulk()
    nt,xyz,uv=A.mesh.get_nodes()
    types,tags,nodes=A.mesh.get_elements()
    owners=Any[]
    for (dim,tag) in A.model.get_entities()
        typ,uid,conn=A.mesh.get_elements(dim,tag)
        push!(owners,Dict("dimension"=>Int(dim),"tag"=>Int(tag),"types"=>Int.(typ),
            "element_tags"=>[Int.(r) for r in uid],"connectivity"=>[Int.(r) for r in conn]))
    end
    return Dict("node_tags"=>Int.(nt),"xyz_bits"=>bits(xyz),"uv_bits"=>bits(uv),
        "element_types"=>Int.(types),"element_tags"=>[Int.(r) for r in tags],
        "connectivity"=>[Int.(r) for r in nodes],"owners"=>owners,
        "max_node_tag"=>Int(A.mesh.get_max_node_tag()),"max_element_tag"=>Int(A.mesh.get_max_element_tag()))
end

function query(uid)
    try
        typ,nodes,dim,owner=A.mesh.get_element(uid)
        return Dict("uid"=>Int(uid),"ok"=>true,"type"=>Int(typ),"nodes"=>Int.(nodes),
            "dimension"=>Int(dim),"owner"=>Int(owner))
    catch err
        err isa InterruptException && rethrow()
        return Dict("uid"=>Int(uid),"ok"=>false,"error"=>sprint(showerror,err))
    end
end

function check_query(actual,expected)
    @test actual["uid"]==expected["uid"]
    @test actual["ok"]==expected["ok"]
    if expected["ok"]
        @test actual==expected
    else
        @test occursin("unknown element",lowercase(actual["error"]))
        @test occursin(string(expected["uid"]),actual["error"])
    end
end

function inspect(children,expected)
    @test Int.(A.mesh.get_visibility(children))==expected["visibility"]
    @test length(children)==length(expected["queries"])
    for (uid,q) in zip(children,expected["queries"]);check_query(query(uid),q);end
    @test bulk()==expected["bulk"]
end

function seed()
    A.initialize();A.model.add("A")
    for tag in (49,50);A.model.add_discrete_entity(2,tag);end
    A.mesh.add_nodes(2,49,[1001,1002,1003,1004],[0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.])
    A.mesh.add_elements_by_type(49,3,[7001],[1001,1002,1003,1004])
    A.mesh.add_nodes(2,50,[2001,2002,2003],[3.,0.,0.,4.,0.,0.,3.,1.,0.])
    A.mesh.add_elements_by_type(50,2,[9001],[2001,2002,2003])
    A.mesh.set_visibility([7001],1);A.mesh.split_quadrangles(1.,49)
    return Int.(only(A.mesh.get_elements(2,49)[2]))
end

function cad_seed()
    A.initialize();A.model.add("native CAD lifetime");A.option("Mesh.Renumber",0)
    for (tag,(x,y)) in zip(101:104,((0.,0.),(1.,0.),(1.,1.),(0.,1.)))
        A.model.add_point(x,y,0.;meshSize=.2,tag)
    end
    for (tag,a,b) in ((11,101,102),(12,102,103),(13,103,104),(14,104,101))
        A.model.add_line(a,b;tag)
    end
    A.model.add_curve_loop([11,12,13,14];tag=31)
    A.model.add_plane_surface([31];tag=49);A.model.synchronize()
    for tag in (11,12,13,14);A.mesh.set_transfinite_curve(tag,3);end
    A.mesh.set_transfinite_surface(49,"Left",[101,102,103,104])
    A.mesh.set_recombine(2,49);A.mesh.generate(2)
    parents=Int.(only(A.mesh.get_elements(2,49)[2]))
    A.mesh.set_visibility(parents,1);A.mesh.split_quadrangles(1.,49)
    return Int.(only(A.mesh.get_elements(2,49)[2]))
end

function mutate(event,children)
    if event=="rebuild_only_if_necessary";A.mesh.rebuild_element_cache(true)
    elseif event=="rebuild_force";A.mesh.rebuild_element_cache(false)
    elseif event=="add_nodes";A.mesh.add_nodes(2,50,[4001],[5.,0.,0.])
    elseif event=="add_nodes_empty";A.mesh.add_nodes(2,50,[],[])
    elseif event=="add_nodes_invalid";A.mesh.add_nodes(2,50,[4001],[5.,0.])
    elseif event=="add_elements_by_type";A.mesh.add_elements_by_type(50,2,[10001],[2001,2002,2003])
    elseif event=="add_elements_by_type_empty";A.mesh.add_elements_by_type(50,2,[],[])
    elseif event=="add_elements_empty";A.mesh.add_elements(2,50,[],[],[])
    elseif event=="reverse_unrelated";A.mesh.reverse([(2,50)])
    elseif event=="reverse_empty";A.mesh.reverse([])
    elseif event=="reverse_elements_unrelated";A.mesh.reverse_elements([9001])
    elseif event=="reverse_elements_empty";A.mesh.reverse_elements([])
    elseif event=="set_node";A.mesh.set_node(2001,[3.125,0.,0.],[])
    elseif event=="set_visibility_missing";A.mesh.set_visibility(children,0)
    elseif event=="affine_unrelated";A.mesh.affine_transform([1.,0.,0.,.125,0.,1.,0.,0.,0.,0.,1.,0.],[(2,50)])
    elseif event=="affine_identity";A.mesh.affine_transform([1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,0.],[(2,49)])
    elseif event=="reorder_children";A.mesh.reorder_elements(2,49,[1,0])
    elseif event=="retag_child_owner";A.model.set_tag(2,49,149)
    elseif event=="retag_unrelated_owner";A.model.set_tag(2,50,150)
    elseif event=="set_order_one";A.mesh.set_order(1)
    elseif event=="set_order_two";A.mesh.set_order(2)
    elseif event=="renumber_elements_identity";A.mesh.renumber_elements([children;9001],[children;9001])
    elseif event=="renumber_elements";A.mesh.renumber_elements([children;9001],[children;9001].+20000);children=children.+20000
    elseif event=="renumber_nodes_identity";old=Int.(A.mesh.get_nodes()[1]);A.mesh.renumber_nodes(old,old)
    elseif event=="renumber_nodes";old=Int.(A.mesh.get_nodes()[1]);A.mesh.renumber_nodes(old,old.+10000)
    elseif event=="remove_unrelated";A.mesh.remove_elements(2,50,[9001])
    elseif event=="remove_empty";A.mesh.remove_elements(2,50,[])
    elseif event=="reclassify_nodes";A.mesh.reclassify_nodes()
    elseif event=="pure_geo_synchronize";A.model.geo.synchronize()
    elseif event=="repeat_geo_synchronize";A.model.geo.synchronize();A.model.geo.synchronize()
    elseif event=="switch_model";A.model.add("B");A.model.set_current("A")
    elseif event=="remove_other_model";A.model.add("B");A.model.remove();A.model.set_current("A")
    else;error(event)
    end
    return children
end

function node_query(uid)
    xyz,uv,dim,owner=A.mesh.get_node(uid)
    return Dict("uid"=>uid,"xyz_bits"=>bits(xyz),"uv_bits"=>bits(uv),"dimension"=>Int(dim),"owner"=>Int(owner))
end

function normalize(value)
    if value isa AbstractFloat
        return Dict("float64_bits"=>only(bits((value,))))
    elseif value isa AbstractVector || value isa Tuple
        return Any[normalize(x) for x in value]
    elseif value isa AbstractDict
        return Dict(k=>normalize(v) for (k,v) in value)
    end
    return value
end

function primer_seed()
    A.initialize();A.model.add("A")
    for tag in (49,50);A.model.add_discrete_entity(2,tag);end
    A.mesh.add_nodes(2,49,[1001,1002,1003,1004],[0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.])
    A.mesh.add_elements_by_type(49,3,[7001],[1001,1002,1003,1004])
    A.mesh.add_nodes(2,50,[2001,2002,2003],[3.,0.,0.,4.,0.,0.,3.,1.,0.])
    A.mesh.add_elements_by_type(50,2,[9001],[2001,2002,2003])
    return nothing
end

function primer(event)
    if event=="get_element";return A.mesh.get_element(9001)
    elseif event=="get_jacobian";return A.mesh.get_jacobian(9001,[.2,.2,0.])
    elseif event=="get_local_coordinates";return A.mesh.get_local_coordinates_in_element(9001,3.2,.2,0.)
    elseif event=="get_qualities";return A.mesh.get_element_qualities([9001],"minSICN")
    elseif event=="get_qualities_empty";return A.mesh.get_element_qualities([],"minSICN")
    elseif event=="get_keys_for_element";return A.mesh.get_keys_for_element(9001,"Lagrange",true)
    elseif event=="get_node";return A.mesh.get_node(2001)
    elseif event=="set_node";return A.mesh.set_node(2001,[3.,0.,0.],[])
    elseif event=="get_nodes";return A.mesh.get_nodes()
    elseif event=="get_nodes_by_element_type";return A.mesh.get_nodes_by_element_type(2,50,false)
    elseif event=="get_elements";return A.mesh.get_elements()
    elseif event=="get_elements_by_type";return A.mesh.get_elements_by_type(2,50)
    elseif event=="get_jacobians";return A.mesh.get_jacobians(2,[.2,.2,0.],50)
    elseif event=="get_keys";return A.mesh.get_keys(2,"Lagrange",50,true)
    elseif event=="get_element_by_coordinates";return A.mesh.get_element_by_coordinates(3.2,.2,0.,2,false)
    elseif event=="get_elements_by_coordinates";return A.mesh.get_elements_by_coordinates(3.2,.2,0.,2,false)
    elseif event=="get_visibility_known";return A.mesh.get_visibility([9001])
    elseif event=="get_visibility_missing";return A.mesh.get_visibility([12345678])
    elseif event=="set_visibility_missing";return A.mesh.set_visibility([12345678],0)
    elseif event=="rebuild_element_necessary";return A.mesh.rebuild_element_cache(true)
    elseif event=="rebuild_node_necessary";return A.mesh.rebuild_node_cache(true)
    elseif event=="unknown_element_abort0"
        A.option("General.AbortOnError",0);return A.mesh.get_element(12345678)
    elseif event=="unknown_node_abort0"
        A.option("General.AbortOnError",0);return A.mesh.get_node(12345678)
    end
    error(event)
end

function authoritative_snapshot()
    m=A.CURRENT[]
    refs=(m,A.LAST_MESH[],A.LAST_MESH_CLASS[],A.LAST_MESH_HIGH_ORDER[],
        A.LAST_MESH_HIGH_ORDER_MESH[],m.meshing.mesh_authority,m.meshing,
        A.LAST_MESH_LOCATOR[],A.LAST_MESH_EDGES[],A.LAST_MESH_FACES[])
    fields=(m,A.LAST_MESH[],A.LAST_MESH_CLASS[],A.LAST_MESH_HIGH_ORDER[],
        A.OPTIONS,A.SESSION_VIEWS,A.NODE_TAG_MAX[],A.ELEMENT_TAG_MAX[])
    return (;refs,hash=digest(fields))
end

function check_snapshot(before)
    after=authoritative_snapshot()
    @test all(a===b for (a,b) in zip(before.refs,after.refs))
    @test after.hash==before.hash
end

function run(output;fresh_case=nothing)
    records=Dict{String,Any}[]
    @testset "Native lookup catalog lifetime" begin
        if fresh_case===nothing
            # Keep the immutable primary's same-process event order, including
            # its accumulated query errors. No deleted native pointer is read.
            for row in NATIVE_LIFETIME64
                @testset "event=$(row["event"]) repeat=$(row["repeat"])" begin
                    children=seed()
                    try
                        @test children==[q["uid"] for q in row["initial"]["queries"]]
                        inspect(children,row["initial"])
                        failed=false
                        try;children=mutate(row["event"],children)
                        catch err;err isa InterruptException && rethrow();failed=true;end
                        @test failed==(row["operation_error"]!==nothing)
                        @test children==row["children"]
                        inspect(Int.(row["retained_children"]),row["after"])
                        if row["unrelated_query"]!==nothing;check_query(query(9001),row["unrelated_query"]);end
                        A.mesh.rebuild_element_cache(false)
                        inspect(Int.(row["retained_children"]),row["rebuilt"])
                        push!(records,Dict("scope"=>"lifetime64","event"=>row["event"],"repeat"=>row["repeat"]))
                    finally;A.finalize();end
                end
            end
            for row in NATIVE_LIVE20
                @testset "live=$(row["kind"])/$(row["value"])/$(row["renumber"]) repeat=$(row["repeat"])" begin
                    children=row["kind"]=="cad" ? cad_seed() : seed()
                    try
                        @test children==row["original_child_tags"];inspect(children,row["initial"])
                        if row["kind"]=="order"
                            A.option("Mesh.Renumber",row["renumber"]);A.mesh.set_order(row["value"])
                        elseif row["kind"]=="node_retag"
                            uid=row["value"]==49 ? 1001 : 2001
                            @test node_query(uid)==row["node_before"]
                            A.model.set_tag(2,row["value"],row["value"]+100)
                            @test node_query(uid)==row["node_after"]
                            row["retained_cached_element"]===nothing || check_query(query(9001),row["retained_cached_element"])
                        elseif row["value"]=="pure_sync";A.model.geo.synchronize()
                        elseif row["value"]=="repeat_sync";A.model.geo.synchronize();A.model.geo.synchronize()
                        elseif row["value"]=="add_point_sync"
                            A.model.geo.add_point(3.,0.,0.;meshSize=.2,tag=201);A.model.geo.synchronize()
                        elseif row["value"]=="add_line_sync"
                            A.model.geo.add_point(3.,0.,0.;meshSize=.2,tag=201)
                            A.model.geo.add_line(101,201;tag=21);A.model.geo.synchronize()
                        end
                        owner=row["kind"]=="node_retag" && row["value"]==49 ? 149 : 49
                        current=Int.(only(A.mesh.get_elements(2,owner)[2]))
                        @test current==row["current_child_tags"];inspect(current,row["after"])
                        A.mesh.rebuild_element_cache(false);inspect(current,row["rebuilt"])
                        push!(records,Dict("scope"=>"live20","kind"=>row["kind"],"repeat"=>row["repeat"]))
                    finally;A.finalize();end
                end
            end
            for row in NATIVE_ERROR16
                order_case(row,records;fresh=false)
            end
            for row in NATIVE_PRIMERS46
                @testset "primer=$(row["event"]) repeat=$(row["repeat"])" begin
                    primer_seed()
                    try
                        @test bulk()==row["initial"]
                        actual=normalize(primer(row["event"]))
                        @test actual==row["result"]
                        @test bulk()==row["after_primer"]
                        A.mesh.split_quadrangles(1.,49)
                        children=Int.(only(A.mesh.get_elements(2,49)[2]))
                        @test children==row["children"]
                        inspect(children,row["after_split"])
                        A.mesh.rebuild_element_cache(false);inspect(children,row["rebuilt"])
                        push!(records,Dict("scope"=>"primers46","event"=>row["event"],"repeat"=>row["repeat"]))
                    finally;A.finalize();end
                end
            end
            @testset "Raw and empty rebuilds preserve authoritative state" begin
                A.initialize()
                try
                    before=authoritative_snapshot()
                    A.mesh.rebuild_node_cache(false);A.mesh.rebuild_element_cache(false)
                    A.mesh.rebuild_node_cache(true);A.mesh.rebuild_element_cache(true)
                    check_snapshot(before)
                    @test A.LAST_MESH[]===nothing
                finally;A.finalize();end
                children=seed()
                try
                    before=authoritative_snapshot()
                    A.mesh.rebuild_node_cache(false);A.mesh.rebuild_element_cache(false)
                    @test A.mesh.get_node(1001)[3:4]==(2,49)
                    @test A.mesh.get_element(9001)[3:4]==(2,50)
                    check_snapshot(before)
                    @test A.LAST_MESH[]===nothing
                finally;A.finalize();end
            end
        else
            # This entry runs in an actual fresh Julia process, exactly like
            # the independent native worker protocol. No private counter reset.
            @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==0
            order_case(NATIVE_FRESH32[fresh_case],records;fresh=true)
        end
    end
    open(output,"w") do io;TOML.print(io,Dict("julia_version"=>string(VERSION),"rows"=>records,
        "native_bits"=>true,"no_deleted_pointer_queries"=>true));end
    return records
end

function order_case(row,records;fresh)
    case=fresh ? row["case"] : row
    @testset "order=$(case["order"]) unknown=$(case["initial_unknown_query"]) abort=$(case["abort_option"]) fresh=$fresh" begin
        children=seed()
        try
            A.option("Mesh.Renumber",fresh ? case["renumber"] : 0)
            case["abort_option"]==0 && A.option("General.AbortOnError",0)
            @test A.option("General.AbortOnError")==row["effective_abort_option"]
            @test Int.(A.mesh.get_visibility(children))==row["initial_visibility"]
            expected=fresh ? row["initial_queries"] : row["initial_error_queries"]
            if case["initial_unknown_query"]
                for (uid,q) in zip(children,expected);check_query(query(uid),q);end
            end
            A.mesh.set_order(case["order"])
            current=Int.(only(A.mesh.get_elements(2,49)[2]))
            @test current==row["current_children"]
            inspect(current,row["after"])
            A.mesh.rebuild_element_cache(false);inspect(current,row["rebuilt"])
            push!(records,Dict("scope"=>fresh ? "fresh32" : "error16","order"=>case["order"],
                "initial_unknown_query"=>case["initial_unknown_query"],"abort_option"=>string(case["abort_option"])))
        finally;A.finalize();end
    end
end
end
