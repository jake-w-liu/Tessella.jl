module NativeRecordPatternWidthTestsV1
using Tessella,Test,TOML
const A=Tessella.API
const DATA=TOML.parsefile(joinpath(@__DIR__,"..","artifacts","native_record_pattern_wide_bulk_v1.toml"))
word_bits(values)=string.(reinterpret(UInt64,Float64.(values));base=16,pad=16)
function nodes(dim=-1,tag=-1)
    ids,xyz,uv=A.mesh.get_nodes(dim,tag,false,true)
    return Dict("node_tags"=>string.(ids),"xyz_bits"=>word_bits(xyz),"uv_bits"=>word_bits(uv))
end
function cells(dim=-1,tag=-1)
    types,ids,connectivity=A.mesh.get_elements(dim,tag)
    return Dict("element_types"=>Int.(types),"element_tags"=>[string.(group) for group in ids],
        "connectivity"=>[string.(group) for group in connectivity])
end
function snapshot()
    return Dict("nodes"=>nodes(),"cells"=>cells(),
        "entities"=>[Dict("dimension"=>dim,"tag"=>tag,"nodes"=>nodes(dim,tag),
            "cells"=>cells(dim,tag)) for (dim,tag) in A.model.get_entities()],
        "physical_groups"=>[Dict("dimension"=>dim,"tag"=>tag,
            "entities"=>Int.(A.model.get_entities_for_physical_group(dim,tag)))
            for (dim,tag) in A.model.get_physical_groups()],
        "max_node_tag"=>string(A.mesh.get_max_node_tag()),
        "max_element_tag"=>string(A.mesh.get_max_element_tag()))
end
function query(row)
    if row["kind"]=="edge"
        return A.mesh.get_element_edge_nodes(row["element_type"],row["tag"],
            row["primary"],row["task"],row["num_tasks"])
    end
    return A.mesh.get_element_face_nodes(row["element_type"],row["face_type"],
        row["tag"],row["primary"],row["task"],row["num_tasks"])
end
function has_box(value)
    value===Core.Box && return true
    value isa GlobalRef && return value.mod===Core && value.name===:Box
    value isa QuoteNode && return has_box(value.value)
    value isa Core.CodeInfo && return any(has_box,value.code)
    value isa Expr && return any(has_box,value.args)
    return false
end
function positive_capture(flag)
    value=0
    callback=()->value
    value=flag ? 1 : 2
    return callback()
end
function allocation_calls(cells,patterns,count)
    return @allocated for _ in 1:count
        A._record_pattern_nodes(cells,patterns)
    end
end
@testset "Native bulk edge and face patterns retain UInt64 record labels" begin
    @test DATA["native_version"]=="4.15.2"
    @test length(DATA["queries"])==182
    A.initialize()
    try
        for entity in DATA["input_entities"]
            dim,tag=entity["dimension"],entity["tag"]
            @test A.model.add_discrete_entity(dim,tag)==tag
            isempty(entity["node_tags"]) && continue
            A.mesh.add_nodes(dim,tag,parse.(UInt64,entity["node_tags"]),
                entity["coordinates"],entity["parameters"])
            A.mesh.add_elements_by_type(tag,entity["element_type"],
                parse.(UInt64,entity["element_tags"]),parse.(UInt64,entity["connectivity"]))
        end
        A.model.add_physical_group(3,[31];tag=501)
        A.model.add_physical_group(2,[21];tag=502)
        @test snapshot()==DATA["before"]
        strict=compact=invalid=0
        for row in DATA["queries"]
            if row["role"]=="invalid_face_pattern" || row["role"]=="invalid_num_tasks"
                # The existing library validates these arguments. Native face
                # patterns return empty; native task zero reports Msg Error.
                @test_throws ArgumentError query(row)
                invalid+=1
            elseif row["tag"]==999
                # Native unknown-owner Error/Abort2 and the current empty record
                # selector differ; exact native rows stay in the fixture.
                @test row["native_error"]
            elseif row["role"]=="partition"
                # The library returns a compact slice. Native returns its full
                # preallocated output with zeros outside this exact slice.
                @test query(row)==parse.(UInt64,row["selected_result"])
                compact+=1
            else
                @test !row["native_error"]
                @test query(row)==parse.(UInt64,row["native_result"])
                strict+=1
            end
            @test snapshot()==DATA["before"]
        end
        @test (strict,compact,invalid)==(36,120,14)
        @test snapshot()==DATA["before"]
        patterns=((1,2),(2,3),(3,1),(4,1),(4,3),(4,2))
        input=first(filter(row->row["element_type"]==4 && !isempty(row["node_tags"]),DATA["input_entities"]))
        labels=parse.(UInt64,input["connectivity"])
        record_cells=[labels[index:index+3] for index in 1:4:length(labels)]
        reference=only(filter(row->row["kind"]=="edge" && row["element_type"]==4 &&
            row["tag"]==31 && !row["primary"] && row["role"]=="selector",DATA["queries"]))
        @test A._record_pattern_nodes(record_cells,patterns)==parse.(UInt64,reference["native_result"])
        signature=Tuple{typeof(record_cells),typeof(patterns)}
        method=which(A._record_pattern_nodes,signature)
        @test realpath(String(method.file))==realpath(joinpath(Base.pkgdir(Tessella),"src","interfaces","API.jl"))
        @test occursin("UInt64",string(method.sig))
        lowered=code_lowered(A._record_pattern_nodes,signature)
        @test !isempty(lowered) && !any(has_box,lowered)
        typed=code_typed(A._record_pattern_nodes,signature;optimize=false)
        @test !isempty(typed) && all(row->!has_box(first(row)),typed)
        @test any(has_box,code_lowered(positive_capture,Tuple{Bool}))
        @test positive_capture(true)==1
        @test positive_capture(false)==2
        A._record_pattern_nodes(record_cells,patterns)
        rows=[(calls=count,bytes=allocation_calls(record_cells,patterns,count)) for count in (1000,2000,4000)]
        @test all(row->row.bytes<=4096*row.calls+1048576,rows)
        @test rows[2].bytes<=2.5*rows[1].bytes+65536
        @test rows[3].bytes<=2.5*rows[2].bytes+65536
        for row in rows
            println("NATIVE_RECORD_PATTERN_RESOURCE calls=",row.calls," bytes=",row.bytes,
                " cap=",4096*row.calls+1048576)
        end
    finally
        A.finalize()
    end
end
end
