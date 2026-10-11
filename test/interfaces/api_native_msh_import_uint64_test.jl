module APINativeMshImportUInt64Tests
using Tessella,Test,TOML,InteractiveUtils
const E=Tessella.Elements
const G=Tessella.GeoExec
const M=Tessella.Model
const A=Tessella.API
const DATA=joinpath(@__DIR__,"..","data","native_msh_import_uint64_v1")
const MANIFEST=TOML.parsefile(joinpath(DATA,"native_controls_v2.toml"))
@test MANIFEST["isolated_native_complete"]

bits(values)=join((string(reinterpret(UInt64,Float64(v));base=16,pad=16) for v in values),",")
uids(values)=join(string.(values),",")
node_signature(uid)=begin
    xyz,uv,dim,tag=A.mesh.get_node(uid)
    "$(uid)|$(dim)|$(tag)|$(bits(xyz))|$(bits(uv))"
end
cell_signature(uid)=begin
    msh,nodes,dim,tag=A.mesh.get_element(uid)
    "$(uid)|$(dim)|$(tag)|$(msh)|$(uids(nodes))"
end
function bulk_signature(dim=-1,tag=-1)
    types,tags,nodes=A.mesh.get_elements(dim,tag)
    return ["$(msh)|$(uids(labels))|$(uids(connected))" for
        (msh,labels,connected) in zip(types,tags,nodes)]
end
function storage_cells(mesh)
    return sort!([(block.msh,Tuple(block.nodes[:,column]),block.tags[column])
        for block in mesh.blocks for column in axes(block.nodes,2)])
end
function assert_public(expected)
    tags,xyz,uv=A.mesh.get_nodes()
    @test eltype(tags)===UInt64
    @test uids(tags)==expected["aggregate_nodes"]
    @test bits(xyz)==expected["aggregate_xyz_bits"]
    @test bits(uv)==expected["aggregate_uv_bits"]
    @test [node_signature(tag) for tag in tags]==expected["node_signatures"]
    @test bulk_signature()==expected["aggregate_cells"]
    @test sort!([cell_signature(parse(UInt64,first(split(signature,'|')))) for
        signature in expected["cell_signatures"]])==sort!(copy(expected["cell_signatures"]))
    @test A.mesh.get_max_node_tag()==parse(UInt64,expected["max_node_uid"])
    @test A.mesh.get_max_element_tag()==parse(UInt64,expected["max_element_uid"])
    entities=get(expected,"entities",Any[])
    @test sort!(Tuple{Int,Int}[(Int(dim),Int(tag)) for
        (dim,tag) in A.model.get_entities()])==
        sort!(Tuple{Int,Int}[(entity["dim"],entity["tag"]) for entity in entities])
    for entity in entities
        dim=entity["dim"];tag=entity["tag"]
        nt,co,pa=A.mesh.get_nodes(dim,tag,false,true)
        @test uids(nt)==entity["nodes"]
        @test bits(co)==entity["xyz_bits"]
        @test bits(pa)==entity["uv_bits"]
        @test bulk_signature(dim,tag)==entity["cells"]
    end
    for uid in get(expected,"missing_nodes",String[])
        @test_throws ErrorException A.mesh.get_node(parse(UInt64,uid))
    end
    # Owned records remain available to model producers; public cache coverage
    # prevents them from becoming a second aggregate node/cell collection.
    @test all(entry->isempty(entry[3].node_tags) && isempty(entry[3].element_tags),
        A._discrete_mesh_records(A.CURRENT[]))
end
function import_source(path;command="Merge",tail="")
    mktempdir() do directory
        geo=joinpath(directory,"import.geo")
        write(geo,"$command \"$(replace(abspath(path),'\\'=>'/'))\";\n"*tail)
        return A.open_geo!(geo)
    end
end
function with_session(f)
    A.initialize()
    try f() finally A.finalize() end
end

@testset "native full MSH import payloads and every supported format" begin
    for control in MANIFEST["cases"]
        @testset "$(control["file"])" begin
        path=joinpath(DATA,control["file"])
        with_session() do
            result=import_source(path)
            @test result.model.meshing.mesh_authority!==nothing
            @test A.LAST_MESH[] isa E.MixedMesh
            assert_public(control)
        end
        # Bulk storage/writer contracts are checked separately from the public
        # import receipt. MSH2 writing deliberately emits dense file labels.
        mesh=E.read_mixed_msh(path)
        @test eltype(mesh.coords)===Float64
        @test all(block->eltype(block.nodes)===Int32,mesh.blocks)
        mktempdir() do directory
            for version in (2.2,4.1),binary in (false,true)
                output=joinpath(directory,"roundtrip-$version-$binary.msh")
                if version==2.2 && mesh.partition_data!==nothing
                    @test_throws ArgumentError E.write_mixed_msh(output,mesh;version,binary)
                    @test !isfile(output)
                    continue
                end
                E.write_mixed_msh(output,mesh;version,binary)
                back=E.read_mixed_msh(output)
                @test reinterpret(UInt64,vec(back.coords))==reinterpret(UInt64,vec(mesh.coords))
                # MSH4 groups cells by entity. A legacy MSH2 input can have
                # those entities interleaved, so serialization changes column
                # order while preserving every oriented connectivity record.
                @test storage_cells(back)==storage_cells(mesh)
                if version==4.1 && mesh.entity_data!==nothing
                    @test back.entity_data.external_node_tags==mesh.entity_data.external_node_tags
                    @test back.entity_data.external_element_tags==mesh.entity_data.external_element_tags
                end
            end
        end
        end
    end
end

@testset "Include/Merge/MergeWithBoundingBox and Save/Print public routes" begin
    for control in filter(c->c["route_family"]=="ordinary_width",MANIFEST["cases"]),
            command in ("Include","Merge","MergeWithBoundingBox")
        path=joinpath(DATA,control["file"])
        with_session() do
            import_source(path;command)
            assert_public(control)
        end
        mktempdir() do directory
            for output_command in ("Save","Print")
                output=joinpath(directory,"export-$output_command.msh")
                with_session() do
                    import_source(path;command,tail="$output_command \"$(replace(output,'\\'=>'/'))\";\n")
                    assert_public(control)
                end
                @test isfile(output)
                back=E.read_mixed_msh(output)
                @test length(back.blocks)==1
                @test back.blocks[1].nodes==reshape(Int32[1,2],2,1)
                @test size(back.coords,2)==2
            end
        end
    end
end

@testset "full UInt64 integral contracts beyond native element import width" begin
    # The pinned DLL narrows external element labels through signed Int32;
    # these are lossless integral/serialization contracts, not DLL parity.
    for fixture in MANIFEST["uint64_contracts"]
        path=joinpath(DATA,fixture["file"])
        nodes=parse.(UInt64,fixture["node_uids"]);element=parse(UInt64,fixture["element_uid"])
        mixed,receipt=E._read_mixed_msh_with_identities(path)
        @test receipt.node_tags==nodes
        @test receipt.element_tags==[UInt64[element]]
        @test mixed.blocks[1].nodes==reshape(Int32[1,2],2,1)
        if fixture["ascii"]
            simplex=Tessella.IO.read_msh(path).mesh
            @test simplex.segs==reshape(Int32[1,2],2,1)
            @test reinterpret(UInt64,vec(simplex.coords))==reinterpret(UInt64,vec(mixed.coords))
            mktempdir() do directory
                for version in (2.2,4.1)
                    output=joinpath(directory,"simplex-$version.msh")
                    Tessella.IO.write_msh(output,simplex;version)
                    @test Tessella.IO.read_msh(output).mesh.segs==simplex.segs
                end
            end
        end
        for binary in (false,true)
            mktempdir() do directory
                output=joinpath(directory,"wide-roundtrip.msh")
                # MSH2 reader receipts enter v4 owned metadata through the same
                # public import normalization, retaining every literal UID.
                actual=G._geo_import_actual_mesh(mixed,receipt,"UInt64 integral contract")
                E.write_mixed_msh(output,actual;version=4.1,binary)
                back=E.read_mixed_msh(output)
                @test back.entity_data.external_node_tags==nodes
                @test back.entity_data.external_element_tags==[UInt64[element]]
                @test back.blocks[1].nodes==mixed.blocks[1].nodes
            end
        end
        with_session() do
            import_source(path)
            @test A.mesh.get_nodes()[1]==nodes
            @test A.mesh.get_elements()[2]==[UInt64[element]]
            @test A.mesh.get_elements()[3]==[nodes]
            @test A.mesh.get_max_node_tag()==maximum(nodes)
            @test A.mesh.get_max_element_tag()==element
            @test A.mesh.get_element(element)==(Int32(1),nodes,Int32(1),Int32(7))
        end
    end
end

@testset "all 125 fixed storage/import contracts retain their literal UID records" begin
    for fixture in MANIFEST["fixed_contracts"]
        path=joinpath(DATA,fixture["file"])
        mixed,receipt=E._read_mixed_msh_with_identities(path)
        nodes=parse.(UInt64,fixture["node_uids"])
        @test receipt.node_tags==nodes
        @test receipt.element_tags==[UInt64[701]]
        @test mixed.blocks[1].msh==fixture["msh"]
        @test mixed.blocks[1].nodes==reshape(collect(Int32,1:fixture["arity"]),:,1)
        with_session() do
            import_source(path)
            @test A.mesh.get_nodes()[1]==nodes
            @test A.mesh.get_elements()[1]==Int32[fixture["msh"]]
            @test A.mesh.get_elements()[2]==[UInt64[701]]
            @test A.mesh.get_elements()[3]==[nodes]
            @test A.mesh.get_element(UInt64(701))[1:2]==(Int32(fixture["msh"]),nodes)
            @test A.mesh.get_max_node_tag()==maximum(nodes)
            @test A.mesh.get_max_element_tag()==UInt64(701)
        end
        mktempdir() do directory
            for version in (2.2,4.1),binary in (false,true)
                output=joinpath(directory,"all-fixed-$version-$binary.msh")
                unsafe=fixture["msh"] in (version==2.2 ?
                    E.GMSH_4_15_2_MSH_READER_GAPS_V2 : E.GMSH_4_15_2_MSH_READER_GAPS_V4)
                if unsafe
                    @test_throws ArgumentError E.write_mixed_msh(output,mixed;version,binary)
                    @test !isfile(output)
                end
                E.write_mixed_msh(output,mixed;version,binary,gmsh_compatible=!unsafe)
                back=E.read_mixed_msh(output)
                @test back.blocks[1].msh==fixture["msh"]
                @test back.blocks[1].nodes==mixed.blocks[1].nodes
                @test reinterpret(UInt64,vec(back.coords))==reinterpret(UInt64,vec(mixed.coords))
            end
        end
    end
end

@testset "fixed-cell corner projection retains columns and constant cells" begin
    for (msh,width,coords,nodes,property) in (
        (8,2,Float64[0 1 2 3 .5 2.5;zeros(2,6)],Int32[1 3;2 4;5 6],:segs),
        (9,3,Float64[0 1 0 2 3 2 .5 .5 0 2.5 2.5 2;
                      0 0 1 0 0 1 0 .5 .5 0 .5 .5;zeros(1,12)],
            Int32[1 4;2 5;3 6;7 10;8 11;9 12],:tris),
        (11,4,Float64[0 1 0 0 2 3 2 2 .5 .5 0 0 .5 0 2.5 2.5 2 2 2.5 2;
                       0 0 1 0 0 0 1 0 0 .5 .5 0 0 .5 0 .5 .5 0 0 .5;
                       0 0 0 1 0 0 0 1 0 0 0 .5 .5 .5 0 0 0 .5 .5 .5],
            Int32[1 5;2 6;3 7;4 8;9 15;10 16;11 17;12 18;13 19;14 20],:tets))
        mixed=E.MixedMesh(coords,[E.ElementBlock(msh,nodes,Int32[0,0])])
        folded=G._geo_mixed_simplex_mesh(mixed,"two-column corner projection")
        @test getproperty(folded,property)==nodes[1:width,:]
    end
    for msh in (84,85,87)
        mixed=E.MixedMesh(zeros(3,1),[E.ElementBlock(msh,ones(Int32,1,1),Int32[0])])
        folded=G._geo_mixed_simplex_mesh(mixed,"constant-cell import")
        @test isempty(folded.segs) && isempty(folded.tris) && isempty(folded.tets)
    end
end

@testset "empty, node-only, retained, mixed, and partition metadata" begin
    for fixture in MANIFEST["retained_controls"]
        with_session() do
            mktempdir() do directory
                paths=[replace(abspath(joinpath(DATA,path)),'\\'=>'/') for path in fixture["files"]]
                geo=joinpath(directory,"retained.geo")
                write(geo,join(["Merge \"$path\";\n" for path in paths]))
                A.open_geo!(geo)
                # These sealed controls explicitly compare rebuilt lookup
                # caches. Import-time cold lookup results are retained in
                # their separate native lifetime controls.
                @test startswith(fixture["native_protocol"],"Explicit rebuild comparison;")
                A.mesh.rebuild_node_cache(false)
                A.mesh.rebuild_element_cache(false)
                assert_public(fixture)
                @test allunique(A.mesh.get_nodes()[1])
                mesh=A.LAST_MESH[]
                @test length(mesh.entity_data.external_node_tags)==length(fixture["node_signatures"])
            end
        end
    end
    for fixture in MANIFEST["partition_controls"]
        path=joinpath(DATA,fixture["file"])
        with_session() do
            result=import_source(path)
            actual=A.LAST_MESH[]
            @test actual.partition_data!==nothing
            @test actual.partition_data.num_partitions==2
            @test actual.partition_data.entities[(1,21)].parent==(1,11)
            @test actual.partition_data.entities[(1,21)].partitions==Int32[1,2]
            @test actual.partition_data.entities[(1,21)].boundaries==Int32[21,-22]
            @test actual.entity_data.external_element_tags==[UInt64[101,102]]
            # Native Ghost curve stores both scalar objects while omitting
            # them from ordinary bulk family collections.
            @test A.mesh.get_elements(1,21)==(Int32[],Vector{UInt64}[],Vector{UInt64}[])
            @test A.mesh.get_elements()==(Int32[],Vector{UInt64}[],Vector{UInt64}[])
            @test A.mesh.get_element(101)==(Int32(1),UInt64[10,30],Int32(1),Int32(21))
            @test A.mesh.get_element(102)==(Int32(1),UInt64[30,20],Int32(1),Int32(21))
            @test A.model.get_type(1,21)=="Ghost curve"
            @test A.model.get_type(0,21)=="Partition point"
            @test A.model.get_parent(1,21)==(-1,-1)
            @test A.model.get_parent(0,21)==(0,1)
            @test A.model.get_partitions(1,21)==[2]
            @test A.model.get_partitions(0,21)==[1]
            @test A.model.get_number_of_partitions()==2
            @test A.mesh.get_node(UInt64(10))[3:4]==(Int32(0),Int32(21))
            @test A.mesh.get_node(UInt64(20))[3:4]==(Int32(0),Int32(22))
            @test length(result.model.meshing.mesh_authority.products[(1,Int32(21))].element_tags[1])==2
            mktempdir() do directory
                for binary in (false,true)
                    output=joinpath(directory,"partition-roundtrip-$binary.msh")
                    E.write_mixed_msh(output,actual;version=4.1,binary)
                    back=E.read_mixed_msh(output)
                    @test back.partition_data.entities[(1,21)].parent==(1,11)
                    @test back.entity_data.external_element_tags==[UInt64[101,102]]
                end
            end
        end
    end
end

@testset "parser negatives, strict UID constructors and resource bounds" begin
    for fixture in MANIFEST["malformed_controls"]
        path=joinpath(DATA,fixture["file"])
        @test_throws ArgumentError E.read_mixed_msh(path)
        fixture["ascii"] && @test_throws ArgumentError Tessella.IO.read_msh(path)
    end
    ordinary=joinpath(DATA,MANIFEST["ordinary_file"])
    @test_throws ArgumentError E.read_mixed_msh(ordinary;max_nodes=1)
    @test_throws ArgumentError E.read_mixed_msh(ordinary;max_elements=0)
    @test_throws ArgumentError E.read_mixed_msh(ordinary;max_connectivity=1)
    @test_throws ArgumentError E.read_mixed_msh(ordinary;max_blocks=0)
    @test_throws ArgumentError E.read_mixed_msh(ordinary;max_entities=0)
    @test_throws ArgumentError E.read_mixed_msh(ordinary;max_file_bytes=filesize(ordinary)-1)
    @test_throws ArgumentError E.read_mixed_msh(ordinary;max_nodes=true)
    @test_throws ArgumentError E.read_mixed_msh(ordinary;max_elements=1.0)
    @test length(E.read_mixed_msh(ordinary;max_nodes=2,max_elements=1,max_connectivity=2,
        max_blocks=2,max_entities=1,max_file_bytes=filesize(ordinary)).coords)==6
    entity=E.MixedEntity(1,7,(0.,0.,0.,1.,0.,0.))
    for value in (true,-1,BigInt(typemax(UInt64))+1,1.0)
        @test_throws ArgumentError E.MixedEntityData(Dict((1,7)=>entity);
            node_entities=[(1,7)],node_parametric=[nothing],external_node_tags=[value])
    end
end

@testset "explicit independently unsupported guards preserve scope" begin
    for fixture in MANIFEST["special_controls"]
        path=joinpath(DATA,fixture["file"])
        if fixture["file"]=="special_parent70_ascii2.msh"
            # Preserve the recovered three-tag input as a negative control:
            # a third MSH2 tag is ambiguous with partition metadata. The
            # storage control supplies the explicit partition prefix.
            @test_throws ArgumentError E._read_mixed_msh_with_identities(path)
            path=joinpath(DATA,"special_parent70_partition_ascii2.msh")
        end
        mesh,receipt=E._read_mixed_msh_with_identities(path)
        @test any(block->block isa E.SpecialElementBlock,mesh.blocks)
        error=try G._geo_import_actual_mesh(mesh,receipt,"special authority requirement");nothing catch e;e end
        @test error isa ArgumentError
        @test occursin("actual special-cell authority",sprint(showerror,error))
        # Storage/serialization of the special cell remains independently
        # supported; this guard does not count as import parity completion.
        mktempdir() do directory
            output=joinpath(directory,"special-roundtrip.msh")
            E.write_mixed_msh(output,mesh;version=fixture["version"],binary=fixture["binary"])
            @test E.mixed_crc(E.read_mixed_msh(output)).sha==E.mixed_crc(mesh).sha
        end
    end
    with_session() do
        path=joinpath(DATA,MANIFEST["ordinary_file"])
        mktempdir() do directory
            geo=joinpath(directory,"duplicate.geo")
            fixed=replace(abspath(path),'\\'=>'/')
            write(geo,"Merge \"$fixed\";\nMerge \"$fixed\";\n")
            error=try A.open_geo!(geo);nothing catch e;e end
            @test error isa ArgumentError
            @test occursin("occurrence identities",sprint(showerror,error))
        end
    end
end

@testset "actual method bindings, UInt64 record methods and linear allocation growth" begin
    @test basename(String(which(G._geo_import_mesh_file!,
        Tuple{M.GeoModel,Tessella.IO._GeoNumericContext,String,String}).file))=="GeoMeshImport.jl"
    @test basename(String(which(M._record_append_element!,
        Tuple{M.DiscreteEntity,Int,UInt64,Vector{UInt64}}).file))=="ModelMeshingAttributes.jl"
    @test !applicable(M._record_append_element!,M.DiscreteEntity(),1,Int32(7),Int32[1,2])
    function boxed_control()
        value=0
        f=()->(value+=1)
        f();return f()
    end
    function has_box(method)
        info=Base.uncompressed_ast(method)
        any(statement->statement isa Expr && occursin("Core.Box",sprint(show,statement)),info.code)
    end
    @test has_box(which(boxed_control,Tuple{}))
    for function_name in (:_geo_import_actual_mesh,:_geo_join_import_actual,
            :_geo_import_order_nodes,:_geo_import_mesh_file!,:_geo_append_discrete_import_nodes!)
        for method in methods(getfield(G,function_name))
            @test !has_box(method)
        end
    end
    for function_name in (:_model_mesh_authority,:_model_mesh_product,
            :_model_mesh_product_entities)
        for method in methods(getfield(M,function_name))
            @test !has_box(method)
        end
    end
    allocated=Int[]
    for count in (256,512,1024)
        nodes=reshape(collect(Int32,1:count),2,:)
        coords=zeros(3,count);coords[1,:].=1:count
        mixed=E.MixedMesh(coords,[E.ElementBlock(1,nodes,zeros(Int32,count÷2))])
        function ingest()
            record=M.DiscreteEntity()
            G._geo_append_discrete_import_nodes!(record,mixed,collect(Int32,1:count),nothing)
            return record
        end
        ingest();GC.gc()
        push!(allocated,@allocated ingest())
    end
    @test allocated[2]<=2.8allocated[1]+65_536
    @test allocated[3]<=2.8allocated[2]+65_536
    println("IO_WIDTH_BATCH_NODE_ALLOCATIONS ",allocated)
    public_allocations=Dict{String,Vector{Int}}()
    for shape in ("one_curve","many_points")
        values=Int[]
        mktempdir() do directory
            for count in (128,256,512)
                path=joinpath(directory,"growth-$shape-$count.msh")
                open(path,"w") do stream
                    println(stream,"\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n\$Nodes\n",count)
                    for row in 1:count
                        println(stream,UInt64(2)^63+UInt64(row)," ",row," 0 0")
                    end
                    println(stream,"\$EndNodes\n\$Elements\n",shape=="many_points" ? count : count÷2)
                    for cell in 1:(shape=="many_points" ? count : count÷2)
                        uid=UInt64(2)^63+UInt64(cell)
                        if shape=="many_points"
                            println(stream,uid," 15 2 0 ",cell," ",UInt64(2)^63+UInt64(cell))
                        else
                            println(stream,uid," 1 2 0 7 ",UInt64(2)^63+UInt64(2cell-1),
                                " ",UInt64(2)^63+UInt64(2cell))
                        end
                    end
                    println(stream,"\$EndElements")
                end
                function public_import()
                    with_session() do
                        import_source(path)
                        @assert length(A.mesh.get_nodes()[1])==count
                        @assert A.mesh.get_max_node_tag()==UInt64(2)^63+UInt64(count)
                    end
                end
                public_import();GC.gc()
                push!(values,@allocated public_import())
            end
        end
        @test values[2]<=3.1values[1]+262_144
        @test values[3]<=3.1values[2]+262_144
        public_allocations[shape]=values
    end
    println("IO_WIDTH_FULL_PUBLIC_IMPORT_ALLOCATIONS ",public_allocations)
end
end
