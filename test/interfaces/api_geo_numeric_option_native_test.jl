module APIGEONumericOptionNativeTests
using Tessella,Test,TOML,SHA
const A=Tessella.API
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
bits(values)=string.(reinterpret(UInt64,Float64.(values));base=16,pad=16)
function check_state(expected)
    @test A.model.get_entities()==Tuple{Int,Int}[Tuple(row) for row in expected["entities"]]
    ids,xyz,uv=A.mesh.get_nodes(-1,-1,false,true)
    @test ids==UInt64.(expected["global_nodes"]["tags"])
    @test bits(xyz)==expected["global_nodes"]["xyz_bits"]
    @test bits(uv)==expected["global_nodes"]["uv_bits"]
    @test A.mesh.get_elements()==(Int32.(expected["global_cells"]["types"]),
        [UInt64.(row) for row in expected["global_cells"]["tags"]],
        [UInt64.(row) for row in expected["global_cells"]["connectivity"]])
    for row in expected["scalar_nodes"]
        xyz,uv,dimension,owner=A.mesh.get_node(row["tag"])
        @test bits(xyz)==row["xyz_bits"]
        @test bits(uv)==row["uv_bits"]
        @test [dimension,owner]==row["owner"]
    end
    for row in expected["scalar_cells"]
        @test A.mesh.get_element(row["tag"])==(Int32(row["type"]),UInt64.(row["nodes"]),
            row["owner"][1],row["owner"][2])
    end
    @test [A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag()]==UInt64.(expected["maxima"])
end
@testset "GEO numeric options persist on success and roll back on import failure" begin
    artifact=joinpath(ROOT,"test","artifacts","numeric_geo_option_native_v1.toml")
    packet=TOML.parsefile(artifact)
    @test bytes2hex(sha256(replace(read(artifact,String),"\r\n"=>"\n")))==
        "6c117cfd8f5a560a4608090bb10c9e0507e79022c0af71957610f5aac1d7ad7b"
    @test packet["gmsh_version"]=="4.15.2" && packet["isolated_runs"]==5
    @test bytes2hex(sha256(replace(read(joinpath(ROOT,packet["capture_script"]),String),
        "\r\n"=>"\n")))==packet["capture_script_sha256"]
    for (name,expected) in sort!(collect(packet["cases"]);by=first)
        @testset "$name" begin
            A.initialize()
            try
                for (key,value) in expected["initial"];A.option(key,value);end
                if name=="aborting"
                    A.model.add_discrete_entity(0,77)
                    A.mesh.add_nodes(0,77,[777],[7.0,7.0,7.0])
                    A.mesh.add_elements_by_type(77,15,[888],[777])
                end
                before_options=Dict(key=>A.option(key) for key in keys(expected["options"]))
                before=(A.model.get_entities(),A.mesh.get_nodes(),A.mesh.get_elements(),
                    A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag())
                mktempdir() do directory
                    source=joinpath(directory,"options.geo");write(source,expected["scene"])
                    @test bytes2hex(sha256(read(source)))==expected["scene_sha256"]
                    failure=try A.open_geo!(source);nothing catch error;error end
                    if isempty(expected["failure"])
                        @test failure===nothing
                    else
                        @test failure isa ErrorException
                        failure isa ErrorException && (@test sprint(showerror,failure)==expected["failure"])
                    end
                end
                # Native failure snapshots retain file writes; API import
                # now stages those writes until successful publication.
                wanted_options=isempty(expected["failure"]) ? expected["options"] : before_options
                for (key,value) in wanted_options
                    @test bits([A.option(key)])==bits([value])
                end
                if name=="aborting"
                    @test (A.model.get_entities(),A.mesh.get_nodes(),A.mesh.get_elements(),
                        A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag())==before
                    @test A.mesh.get_element(888)==(Int32(15),UInt64[777],0,77)
                else
                    check_state(expected["after"])
                end
            finally A.finalize() end
        end
    end
end
@testset "native GEO option aliases share parser and API slots" begin
    artifact=joinpath(ROOT,"test","artifacts","numeric_geo_alias_native_v1.toml")
    packet=TOML.parsefile(artifact)
    @test bytes2hex(sha256(replace(read(artifact,String),"\r\n"=>"\n")))==
        "8f9d536375a55045300a941956ff7972b83663af46d2057351f02ad6ace76c6f"
    @test packet["gmsh_version"]=="4.15.2"
    @test bytes2hex(sha256(replace(read(joinpath(ROOT,packet["capture_script"]),String),
        "\r\n"=>"\n")))==packet["capture_script_sha256"]
    A.initialize()
    try
        mktempdir() do directory
            source=joinpath(directory,"aliases.geo");write(source,packet["scene"])
            @test bytes2hex(sha256(read(source)))==packet["scene_sha256"]
            # Pure GEO and API import consume the same native option aliases.
            @test bits(execute_geo(source).lists["copied"])==bits(packet["values"])
            @test bits(A.open_geo!(source).lists["copied"])==bits(packet["values"])
        end
        for (key,value) in packet["options"]
            @test bits([A.option(key)])==bits([value])
        end
    finally A.finalize() end
end
end
