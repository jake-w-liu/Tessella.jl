module APIPartitionGhostCollisionNativeTests
using Tessella,Test,TOML,SHA
const A=Tessella.API
include(joinpath(@__DIR__,"native_public_state.jl"))
using .NativePublicState:assert_state
const NATIVE=TOML.parsefile(joinpath(@__DIR__,"..","artifacts","partition_ghost_collision_native_v1.toml"))
@testset "native partition/Ghost declaration collision keeps bulk and scalar views" begin
    @test NATIVE["gmsh_version"]=="4.15.2" && length(NATIVE["cases"])==2
    for (name,expected) in sort!(collect(NATIVE["cases"]);by=first)
        @testset "$name" begin
            source=joinpath(@__DIR__,"..","..",expected["source"])
            @test bytes2hex(sha256(read(source)))==expected["source_sha256"]
            A.initialize()
            try
                previous_errors=A.NATIVE_LOOKUP_QUERY_ERRORS[]
                previous_last_error=A.logger.get_last_error()
                mktempdir() do directory
                    wrapper=joinpath(directory,"import.geo")
                    write(wrapper,"Merge \""*replace(source,'\\'=>'/')*"\";\n")
                    A.open_geo!(wrapper)
                end
                assert_state(expected["state"])
                # Native Ghost cells remain directly addressable even when
                # they are absent from both global and entity bulk views.
                @test isempty(A.mesh.get_elements()[1])
                @test isempty(A.mesh.get_elements(1,21)[1])
                for uid in (101,102)
                    native=expected["element$uid"]
                    kind,connected,dimension,owner=A.mesh.get_element(uid)
                    @test kind==native[1]
                    @test connected==UInt64.(native[2])
                    @test (dimension,owner)==(native[3],native[4])
                end
                tags,partitions=A.mesh.get_ghost_elements(1,21)
                @test tags==UInt64.(expected["ghost"][1])
                @test partitions==Int32.(expected["ghost"][2])
                @test A.NATIVE_LOOKUP_QUERY_ERRORS[]-previous_errors==length(expected["error_events"])
                @test A.logger.get_last_error()==previous_last_error
            finally A.finalize() end
        end
    end
end
@testset "zero-dimensional partition declaration keeps ordinary Point cells" begin
    expected=TOML.parsefile(joinpath(@__DIR__,"..","artifacts","partition_point_native_v1.toml"))
    source=joinpath(@__DIR__,"..","..",expected["source"])
    @test expected["gmsh_version"]=="4.15.2"
    @test bytes2hex(sha256(read(source)))==expected["source_sha256"]
    A.initialize()
    try
        mktempdir() do directory
            wrapper=joinpath(directory,"import.geo")
            write(wrapper,"Merge \""*replace(source,'\\'=>'/')*"\";\n")
            A.open_geo!(wrapper)
        end
        assert_state(expected["state"])
        @test A.model.get_number_of_partitions()==expected["num_partitions"]
        native=expected["element101"]
        @test A.mesh.get_element(101)==(native[1],UInt64.(native[2]),native[3],native[4])
        @test A.mesh.get_element_types(0,21)==Int32[15]
        @test A.mesh.get_elements_by_type(15,21)==(UInt64[101],UInt64[10])
        @test A.model.get_type(0,21)=="Partition point"
        parts=A.model.get_partitions(0,21)
        empty!(parts)
        @test A.model.get_partitions(0,21)==[1,2]
    finally A.finalize() end
end
@testset "native empty and assigned Ghost declarations" begin
    packet=TOML.parsefile(joinpath(@__DIR__,"..","artifacts","partition_followups_native_v1.toml"))
    @test packet["gmsh_version"]=="4.15.2" && length(packet["cases"])==5
    for (name,expected) in sort!(collect(packet["cases"]);by=first)
        @testset "$name" begin
            source=joinpath(@__DIR__,"..","..",expected["source"])
            @test bytes2hex(sha256(read(source)))==expected["source_sha256"]
            A.initialize()
            try
                mktempdir() do directory
                    wrapper=joinpath(directory,"import.geo")
                    write(wrapper,"Merge \""*replace(source,'\\'=>'/')*"\";\n")
                    A.open_geo!(wrapper)
                end
                assert_state(expected["state"])
                @test A.model.get_number_of_partitions()==expected["num_partitions"]
                @test A.mesh.get_ghost_elements(1,21)==
                    (UInt64.(expected["ghost"][1]),Int32.(expected["ghost"][2]))
            finally A.finalize() end
        end
    end
end
@testset "native Ghost alias keeps shared element occurrences" begin
    expected=TOML.parsefile(joinpath(@__DIR__,"..","artifacts","partition_alias_native_v1.toml"))
    source=joinpath(@__DIR__,"..","..",expected["source"])
    @test expected["gmsh_version"]=="4.15.2"
    @test bytes2hex(sha256(read(source)))==expected["source_sha256"]
    A.initialize()
    try
        mktempdir() do directory
            wrapper=joinpath(directory,"import.geo")
            write(wrapper,"Merge \""*replace(source,'\\'=>'/')*"\";\n")
            A.open_geo!(wrapper)
        end
        assert_state(expected["state"])
        @test A.model.get_number_of_partitions()==expected["num_partitions"]
        for tag in (-1,21,31)
            query=expected["queries"][string(tag)]
            @test A.mesh.get_element_types(1,tag)==Int32.(query["types"])
            @test A.mesh.get_elements_by_type(1,tag)==Tuple(UInt64.(v) for v in query["by_type"])
            nt,xyz,uv=A.mesh.get_nodes_by_element_type(1,tag,true)
            @test nt==UInt64.(query["nodes"])
            @test NativePublicState.bits(xyz)==query["xyz_bits"]
            @test NativePublicState.bits(uv)==query["uv_bits"]
            @test NativePublicState.bits(A.mesh.get_barycenters(1,tag,false,false))==query["barycenter_bits"]
            for task in (0,1,2)
                # Tessella's documented task result is compact. The native
                # packet retains full zero-padded buffers; select its actual
                # task positions by the native begin/end formula, preserving
                # zero UIDs if they occur inside the populated task range.
                count=length(query["by_type"][1])
                chosen=min(fld(task*count,2),count)+1:min(fld((task+1)*count,2),count)
                raw=query["task_by_type"][task+1]
                native_tags=UInt64.(raw[1]);native_nodes=reshape(UInt64.(raw[2]),2,count)
                @test length(native_tags)==count
                @test all(i->i in chosen || iszero(native_tags[i]),eachindex(native_tags))
                @test A.mesh.get_elements_by_type(1,tag,task,2)==
                    (native_tags[chosen],vec(native_nodes[:,chosen]))
            end
        end
        @test A.mesh.get_ghost_elements(1,31)==(UInt64[102],Int32[1])
        # Both bulk occurrences refer to one stored cell and one scalar UID.
        @test sum(length, A.LAST_MESH[].entity_data.external_element_tags)==2
        @test A.mesh.get_element(102)==(1,UInt64[30,20],1,21)
    finally A.finalize() end
end
@testset "native Ghost kinds follow model dimension at each import event" begin
    packet=TOML.parsefile(joinpath(@__DIR__,"..","artifacts","partition_model_order_native_v1.toml"))
    @test packet["gmsh_version"]=="4.15.2"
    files=Dict{String,String}()
    for (name,expected) in packet["files"]
        source=joinpath(@__DIR__,"..","..",expected["source"])
        @test bytes2hex(sha256(read(source)))==expected["sha256"]
        files[name]=source
    end
    for (name,expected) in sort!(collect(packet["cases"]);by=first)
        @testset "$name" begin
            A.initialize()
            try
                mktempdir() do directory
                    wrapper=joinpath(directory,"imports.geo")
                    write(wrapper,join(("Merge \""*replace(files[key],'\\'=>'/')*"\";\n"
                        for key in expected["order"])))
                    A.open_geo!(wrapper)
                end
                assert_state(expected["state"])
                @test A.model.get_number_of_partitions()==expected["num_partitions"]
            finally A.finalize() end
        end
    end
end
@testset "native imported unpartition transfers and removes owned objects" begin
    packet=TOML.parsefile(joinpath(@__DIR__,"..","artifacts","partition_unpartition_native_v1.toml"))
    @test packet["gmsh_version"]=="4.15.2" && length(packet["cases"])==6
    for (name,expected) in sort!(collect(packet["cases"]);by=first)
        @testset "$name" begin
            source=joinpath(@__DIR__,"..","..",expected["source"])
            @test bytes2hex(sha256(read(source)))==expected["source_sha256"]
            A.initialize()
            try
                mktempdir() do directory
                    wrapper=joinpath(directory,"import.geo")
                    write(wrapper,"Merge \""*replace(source,'\\'=>'/')*"\";\n")
                    A.open_geo!(wrapper)
                end
                expected["warm_before"] && assert_state(expected["before"])
                @test A.mesh.unpartition()===nothing
                assert_state(expected["after"])
                @test A.model.get_number_of_partitions()==expected["num_partitions"]
                @test A.LAST_MESH[].partition_data===nothing
                @test A.mesh.unpartition()===nothing
                assert_state(expected["after"])
            finally A.finalize() end
        end
    end
end
@testset "native Ghost membership survives shared partitions and label changes" begin
    root=normpath(joinpath(@__DIR__,"..",".."))
    artifact=joinpath(root,"test","artifacts","partition_ghost_membership_native_v2.toml")
    @test bytes2hex(sha256(replace(read(artifact,String),"\r\n"=>"\n")))==
        "d0dc2fd15a85cbd14ce3ef9f54191549223d57194cdacd0859c8f8b923f8f114"
    packet=TOML.parsefile(artifact)
    @test packet["gmsh_version"]=="4.15.2" && packet["isolated_runs"]==30
    @test packet["all_three_samples_exact"] && packet["cold_ghost_queries_before_rebuild"]
    @test bytes2hex(sha256(replace(read(joinpath(root,packet["capture_script"]),String),"\r\n"=>"\n")))==
        packet["capture_script_sha256"]
    for (name,expected) in sort!(collect(packet["cases"]);by=first)
        @testset "$name" begin
            @test length(expected["raw_sha256"])==3
            @test Set(expected["order"])==Set(keys(expected["files"]))
            for (file,digest) in expected["files"]
                @test bytes2hex(sha256(read(joinpath(root,file))))==digest
            end
            A.initialize()
            try
                operation=expected["operation"]
                mktempdir() do directory
                    wrapper=joinpath(directory,"import.geo")
                    tail=operation=="default_geo" ? "RenumberMeshElements;\n" :
                        operation=="mesh0_geo" ? "Mesh 0;\n" : ""
                    write(wrapper,join(("Merge \""*replace(joinpath(root,file),'\\'=>'/')*"\";\n"
                        for file in expected["order"]))*tail)
                    A.open_geo!(wrapper)
                end
                if operation=="explicit"
                    @test A.mesh.renumber_elements([102],[502])===nothing
                elseif operation=="default_api"
                    @test A.mesh.renumber_elements()===nothing
                elseif operation=="mesh0_api"
                    A.mesh.generate(0)
                end
                # Compare Ghost membership before any scalar cache rebuild,
                # exactly as the isolated native captures do.
                for query in expected["cold_ghost_queries"]
                    @test A.mesh.get_ghost_elements(query["owner"]...)==
                        (parse.(UInt64,query["uids"]),Int32.(query["partitions"]))
                end
                A.mesh.rebuild_node_cache(false)
                A.mesh.rebuild_element_cache(false)
                assert_state(expected["state"])
                @test A.model.get_number_of_partitions()==expected["num_partitions"]
            finally A.finalize() end
        end
    end
end
@testset "Ghost ordering matches a complete native fork throughout mutations" begin
    root=normpath(joinpath(@__DIR__,"..",".."))
    artifact=joinpath(root,"test","artifacts","partition_ghost_order_native_v1.toml")
    @test bytes2hex(sha256(replace(read(artifact,String),"\r\n"=>"\n")))==
        "c4158870b90fdfc3a0c9c10a18def2249f1c67e15a1e4047fcea825b50e86648"
    packet=TOML.parsefile(artifact)
    @test packet["gmsh_version"]=="4.15.2" && packet["isolated_runs"]==12
    @test packet["complete_coupled_lifetimes"] && packet["native_pointer_order_can_vary"]
    @test packet["deterministic_source_birth_order_extension"] && !packet["raw_order_exact_parity_claim"]
    for prefix in ("capture_script","snapshot_helper")
        @test bytes2hex(sha256(replace(read(joinpath(root,packet[prefix]),String),"\r\n"=>"\n")))==
            packet[prefix*"_sha256"]
    end
    for (name,expected) in sort!(collect(packet["cases"]);by=first)
        @testset "$name" begin
            @test length(expected["raw_sha256"])==6
            @test bytes2hex(sha256(read(joinpath(root,expected["source"]))))==expected["source_sha256"]
            A.initialize()
            try
                mktempdir() do directory
                    geo=joinpath(directory,"import.geo")
                    write(geo,"Merge \""*replace(joinpath(root,expected["source"]),'\\'=>'/')*"\";\n")
                    A.open_geo!(geo)
                end
                observed=Any[]
                for action in ("import","reorder","renumber")
                    action=="reorder" && A.mesh.reorder_elements(1,21,[1,0])
                    action=="renumber" && A.mesh.renumber_elements()
                    uids,parts=A.mesh.get_ghost_elements(1,31)
                    # Native's same-family UID and bulk orders agree within
                    # every captured process, including after a parent reorder.
                    @test uids==only(A.mesh.get_elements(1,31)[2])
                    A.mesh.rebuild_node_cache(false);A.mesh.rebuild_element_cache(false)
                    push!(observed,Dict("action"=>action,
                        "ghost_queries"=>[Dict("owner"=>[1,31],"uids"=>string.(uids),
                            "partitions"=>Int.(parts))],"state"=>NativePublicState.public_state()))
                end
                # Match the entire operation lifetime to one measured fork;
                # choosing unrelated variants at separate steps cannot pass.
                @test any(variant->variant==observed,expected["variants"])
            finally A.finalize() end
        end
    end
end
end
