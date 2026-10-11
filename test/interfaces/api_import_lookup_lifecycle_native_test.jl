module APIImportLookupLifecycleNativeTests
using Tessella,Test,TOML,SHA
const A=Tessella.API
const ROOT=joinpath(@__DIR__,"..","..")
const ARTIFACT=joinpath(@__DIR__,"..","artifacts","import_lookup_lifecycle_native_v1.toml")
const PACKET=TOML.parsefile(ARTIFACT)
@testset "native MSH4 read-time element cache and generation lifetime" begin
    @test bytes2hex(sha256(replace(read(ARTIFACT,String),"\r\n"=>"\n")))==
        "366212ce21ac6c4512a5678a80477fe5b7f064dab2bccfaf92ed1e514a2c2a6e"
    @test PACKET["gmsh_version"]=="4.15.2" && length(PACKET["cases"])==7
    source=joinpath(ROOT,PACKET["source"])
    @test bytes2hex(sha256(read(source)))==PACKET["source_sha256"]
    capture=joinpath(ROOT,PACKET["capture_script"])
    @test bytes2hex(sha256(replace(read(capture,String),"\r\n"=>"\n")))==PACKET["capture_script_sha256"]
    for (name,expected) in sort!(collect(PACKET["cases"]);by=first)
        @testset "$name" begin
            A.initialize()
            try
                A.option("Mesh.Renumber",0)
                mktempdir() do directory
                    body=read(source,String)
                    ghost=findfirst("\$GhostElements",body)
                    if name in ("no_ghost","zero_ghost")
                        body=body[1:prevind(body,first(ghost))]
                        name=="zero_ghost" && (body*="\$GhostElements\n0\n\$EndGhostElements\n")
                    end
                    msh=joinpath(directory,"partition.msh");write(msh,body)
                    @test bytes2hex(sha256(read(msh)))==expected["source_msh_sha256"]
                    geo=joinpath(directory,"import.geo")
                    command=name=="geo_mesh" ? "Mesh.Renumber=0; Mesh 0;\n" :
                        name=="geo_renumber" ? "RenumberMeshNodes; RenumberMeshElements;\n" : ""
                    write(geo,"Merge \""*replace(msh,'\\'=>'/')*"\";\n"*command)
                    A.open_geo!(geo)
                end
                # No scalar or bulk element query warms a catalogue before
                # Unpartition: the native MSH4 reader itself primes it.
                if name=="api_mesh"
                    A.mesh.generate(0)
                elseif name=="api_failed_mesh"
                    A.mesh.add_homology_request("Homology",[99],[],[0])
                    failure=try A.mesh.generate(0);nothing catch error;error end
                    @test failure isa ErrorException
                    @test sprint(showerror,failure)==expected["error"]
                    @test A.logger.get_last_error()==expected["last_error"]
                end
                @test A.mesh.unpartition()===nothing
                @test A.mesh.get_elements()==(Int32.(expected["bulk_types"]),
                    [UInt64.(tags) for tags in expected["bulk_tags"]],
                    [UInt64.(nodes) for nodes in expected["bulk_nodes"]])
                scalar=(Int32(expected["scalar_type"]),UInt64.(expected["scalar_nodes"]),
                    expected["scalar_dimension"],expected["scalar_owner"])
                @test A.mesh.get_element(expected["scalar_uid"])==scalar
                @test A.mesh.get_max_node_tag()==expected["maximum_node"]
                @test A.mesh.get_max_element_tag()==expected["maximum_element"]
                @test A.mesh.unpartition()===nothing
                @test A.mesh.get_element(expected["scalar_uid"])==scalar
            finally A.finalize() end
        end
    end
end

bits(values)=string.(reinterpret(UInt64,Float64.(values));base=16,pad=16)
function check_retained_queries(expected)
    @test A.mesh.get_elements()==(Int32.(expected["bulk_types"]),
        [UInt64.(tags) for tags in expected["bulk_tags"]],
        [UInt64.(nodes) for nodes in expected["bulk_nodes"]])
    for (uid,query) in sort!(collect(expected["elements"]);by=row->parse(UInt64,first(row)))
        tag=parse(UInt64,uid)
        if query["status"]=="error"
            failure=try A.mesh.get_element(tag);nothing catch error;error end
            @test failure isa ArgumentError
            failure isa ArgumentError && (@test failure.msg==query["message"])
        else
            @test A.mesh.get_element(tag)==(Int32(query["type"]),UInt64.(query["nodes"]),
                query["dimension"],query["owner"])
        end
    end
    for (uid,query) in sort!(collect(expected["node_queries"]);by=row->parse(UInt64,first(row)))
        tag=parse(UInt64,uid)
        if query["status"]=="error"
            failure=try A.mesh.get_node(tag);nothing catch error;error end
            @test failure isa ArgumentError
            failure isa ArgumentError && (@test failure.msg==query["message"])
        else
            xyz,uv,dimension,owner=A.mesh.get_node(tag)
            @test bits(xyz)==query["xyz_bits"]
            @test bits(uv)==query["uv_bits"]
            @test (dimension,owner)==(query["dimension"],query["owner"])
        end
    end
end
@testset "native retained import cache misses and explicit rebuild" begin
    artifact=joinpath(@__DIR__,"..","artifacts","retained_lookup_native_v1.toml")
    packet=TOML.parsefile(artifact)
    @test bytes2hex(sha256(replace(read(artifact,String),"\r\n"=>"\n")))==
        "e776d0cf883e444610b969bf2fe4dc8a00374cfbb09eb643ba330f50b577fe83"
    @test packet["gmsh_version"]=="4.15.2" && length(packet["cases"])==4
    capture=joinpath(ROOT,packet["capture_script"])
    @test bytes2hex(sha256(replace(read(capture,String),"\r\n"=>"\n")))==packet["capture_script_sha256"]
    for (name,profiles) in sort!(collect(packet["cases"]);by=first)
        @testset "$name" begin
            A.initialize()
            try
                expected=profiles["cold"]
                paths=[joinpath(ROOT,"test","data","native_msh_import_uint64_v1",file)
                    for file in expected["files"]]
                @test bytes2hex.(sha256.(read.(paths)))==expected["source_sha256"]
                mktempdir() do directory
                    source=joinpath(directory,"import.geo")
                    write(source,join("Merge \"$(replace(path,'\\'=>'/'))\";\n" for path in paths))
                    A.open_geo!(source)
                end
                check_retained_queries(expected)
                # A forced native rebuild is the manifest's declared protocol;
                # the conditional form keeps its existing miss entries.
                @test A.mesh.rebuild_node_cache(true)===nothing
                @test A.mesh.rebuild_element_cache(true)===nothing
                check_retained_queries(expected)
                @test A.mesh.rebuild_node_cache(false)===nothing
                @test A.mesh.rebuild_element_cache(false)===nothing
                check_retained_queries(profiles["rebuilt"])
            finally A.finalize() end
        end
    end
end
function check_element_query(uid,expected)
    if expected["status"]=="error"
        failure=try A.mesh.get_element(uid);nothing catch error;error end
        @test failure isa ArgumentError
        failure isa ArgumentError && (@test failure.msg==expected["message"])
    else
        @test A.mesh.get_element(uid)==(Int32(expected["type"]),UInt64.(expected["nodes"]),
            expected["dimension"],expected["owner"])
    end
end
@testset "native import lookup lifetime under GEO conversions and Ghost rebuilds" begin
    artifact=joinpath(@__DIR__,"..","artifacts","import_mutation_lookup_native_v1.toml")
    packet=TOML.parsefile(artifact)
    @test bytes2hex(sha256(replace(read(artifact,String),"\r\n"=>"\n")))==
        "af70a5be98efe8cc4ea149dfa4579ad675bc4cffe254a1e86de1779f52dd0479"
    @test packet["gmsh_version"]=="4.15.2" && packet["isolated_runs"]==6
    @test bytes2hex(sha256(replace(read(joinpath(ROOT,packet["capture_script"]),String),
        "\r\n"=>"\n")))==packet["capture_script_sha256"]
    for (name,case) in sort!(collect(packet["cases"]);by=first)
        @testset "$name" begin
            source=joinpath(ROOT,case["source"])
            @test bytes2hex(sha256(read(source)))==case["source_sha256"]
            A.initialize()
            try
                mktempdir() do directory
                    path=joinpath(directory,"import.geo")
                    write(path,"Mesh.Renumber=0;\nMerge \""*replace(source,'\\'=>'/')*
                        "\";\n"*case["command"]*"\n")
                    A.open_geo!(path)
                end
                expected=case["payload"]
                if case["mutation"]
                    @test A.mesh.unpartition()===nothing
                    @test A.mesh.get_elements()==(Int32.(expected["bulk_types"]),
                        [UInt64.(tags) for tags in expected["bulk_tags"]],
                        [UInt64.(nodes) for nodes in expected["bulk_nodes"]])
                    check_element_query(101,expected["scalar"])
                    xyz,uv,dimension,owner=A.mesh.get_node(10)
                    @test bits(xyz)==expected["node"]["xyz_bits"]
                    @test bits(uv)==expected["node"]["uv_bits"]
                    @test (dimension,owner)==(expected["node"]["dimension"],expected["node"]["owner"])
                    @test A.model.get_number_of_partitions()==expected["partitions"]
                    @test A.mesh.unpartition()===nothing
                    check_element_query(101,expected["scalar"])
                else
                    bulk=A.mesh.get_elements()
                    stored=deepcopy(A.LAST_MESH[].entity_data.external_element_tags)
                    for stage in ("cold","conditional","rebuilt")
                        stage=="conditional" && (@test A.mesh.rebuild_element_cache(true)===nothing)
                        stage=="rebuilt" && (@test A.mesh.rebuild_element_cache(false)===nothing)
                        for uid in (101,102)
                            check_element_query(uid,expected[stage][string(uid)])
                        end
                        @test A.mesh.get_elements()==bulk
                        @test A.LAST_MESH[].entity_data.external_element_tags==stored
                    end
                end
                @test (A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag())==
                    Tuple(UInt64.(expected["maxima"]))
            finally A.finalize() end
        end
    end
end
end
