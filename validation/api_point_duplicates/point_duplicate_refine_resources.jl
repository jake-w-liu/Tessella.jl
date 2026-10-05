# Actual public refinement of distinct Point cells sharing one primary node.
using Tessella,Test,SHA,TOML
const PointAPI=Tessella.API

function point_duplicate_resource_fixture(count::Int)
    PointAPI.initialize()
    PointAPI.option("Mesh.Renumber",0)
    PointAPI.model.add_discrete_entity(0,91)
    mesh=Tessella.Elements.MixedMesh(reshape([3.,3.,0.],3,1),
        [Tessella.Elements.ElementBlock(15,ones(Int32,1,count))])
    owners=[(0,Int32(91))]
    class=PointAPI._mixed_classification(mesh,(0,Int32(91)),owners,owners,
        Dict{Tuple{Int,Int32},Vector{Int32}}(),Dict(15=>fill(Int32(91),count)))
    labels=UInt64.(10000 .+ 3collect(1:count))
    table=PointAPI._cache_public_tags(mesh;node_tags=[1],element_tags=labels)
    PointAPI._replace_mesh_cache_locked!(mesh,PointAPI._classification_with_public_tags(class,table))
    return labels
end

function point_duplicate_resource_rows()
    rows=Dict{String,Any}[]
    try
        for count in (1000,2000,4000)
            labels=point_duplicate_resource_fixture(count)
            PointAPI.mesh.refine();PointAPI.mesh.refine()
            samples=[@timed(PointAPI.mesh.refine()) for _ in 1:3]
            selected=samples[argmin([sample.bytes for sample in samples])]
            actual,connections=PointAPI.mesh.get_elements_by_type(15,91)
            @test actual==labels && connections==ones(UInt64,count)
            @test PointAPI.mesh.get_nodes(0,91)==(UInt64[1],[3.,3.,0.],Float64[])
            @test PointAPI.mesh.get_max_element_tag()==last(labels)
            @test PointAPI.mesh.get_element(first(labels))==(15,UInt64[1],0,91)
            @test PointAPI.mesh.get_element(last(labels))==(15,UInt64[1],0,91)
            public=PointAPI.LAST_MESH_CLASS[].public_tags
            @test public.element_tags==labels && length(public.element_indices)==count
            push!(rows,Dict("cells"=>count,"nodes"=>1,"allocated"=>selected.bytes,
                "allocation_samples"=>[sample.bytes for sample in samples],
                "allocation_count"=>Base.gc_alloc_count(selected.gcstats),"seconds"=>selected.time,
                "retained_bytes"=>Base.summarysize(public),
                "connectivity_sha256"=>bytes2hex(sha256(reinterpret(UInt8,connections))),
                "element_labels_sha256"=>bytes2hex(sha256(reinterpret(UInt8,actual)))))
        end
        for index in 2:length(rows)
            @test rows[index]["allocated"]<=2.15*rows[index-1]["allocated"]+65536
            @test rows[index]["retained_bytes"]<=2.15*rows[index-1]["retained_bytes"]+65536
        end
    finally
        PointAPI.finalize()
    end
    return rows
end

function point_duplicate_boxes(value)
    value isa Core.Box && return 1
    value isa GlobalRef && return value==GlobalRef(Core,:Box) ? 1 : 0
    value isa Core.CodeInfo && return sum(point_duplicate_boxes,value.code;init=0)
    value isa Expr && return sum(point_duplicate_boxes,value.args;init=0)
    value isa Core.NewvarNode && return point_duplicate_boxes(value.slot)
    value isa QuoteNode && return point_duplicate_boxes(value.value)
    return value===Core.Box ? 1 : 0
end

function point_duplicate_box_control()
    captured=0;read=()->captured;captured=1
    return read
end

function point_duplicate_keyword_box_control(;initial=0)
    captured=initial;read=()->captured;captured=1
    return read
end

function point_duplicate_keyword_bodies(method::Method)
    result=Set([method])
    prefix="#"*String(method.name)*"#"
    for name in names(method.module;all=true)
        startswith(String(name),prefix) || continue
        value=getfield(method.module,name)
        if value isa Function
            union!(result,methods(value))
        elseif value isa Type && value<:Function
            matches=Base._methods_by_ftype(Tuple{value,Vararg{Any}},-1,Base.get_world_counter())
            matches===nothing || foreach(match->push!(result,match.method),matches)
        end
    end
    return result
end

rows=Dict{String,Any}[];ast=Dict{String,Any}[]
@testset "Point cell refinement resources and actual implementation bodies" begin
    @test point_duplicate_boxes(Base.uncompressed_ast(first(methods(point_duplicate_box_control))))>0
    @test sum(method->point_duplicate_boxes(Base.uncompressed_ast(method)),
        point_duplicate_keyword_bodies(first(methods(point_duplicate_keyword_box_control)));init=0)>0
    inventory=Set{Method}()
    for function_value in (Tessella.Elements.validate,Tessella.Elements.write_mixed_msh,
                           PointAPI._mixed_refine_locked!,
                           point_duplicate_resource_fixture,point_duplicate_resource_rows)
        for method in methods(function_value)
            method.module in (Tessella.Elements,Tessella.API,@__MODULE__) || continue
            union!(inventory,point_duplicate_keyword_bodies(method))
        end
    end
    for method in sort!(collect(inventory);by=string)
        boxes=point_duplicate_boxes(Base.uncompressed_ast(method))
        @test boxes==0
        push!(ast,Dict("method"=>string(method),"boxes"=>boxes))
    end
    @test any(method->method.module===Tessella.Elements &&
        startswith(String(method.name),"#validate#"),inventory)
    @test any(method->method.name===:_mixed_refine_locked!,inventory)
    append!(rows,point_duplicate_resource_rows())
end
if !isempty(ARGS)
    @assert !ispath(only(ARGS))
    open(only(ARGS),"w") do stream
        TOML.print(stream,Dict("julia"=>string(VERSION),"rows"=>rows,"ast"=>ast,
                              "growth_factor"=>2.15,"growth_slack_bytes"=>65536))
    end
end
println("POINT_DUPLICATE_REFINE_RESOURCE_COMPLETE rows=",length(rows)," actual_bodies=",length(ast))
