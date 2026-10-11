module APINativeSelector4ViewTests
using Test,Tessella,TOML,SHA,Serialization
if Sys.iswindows() && Sys.ARCH===:x86_64
const A=Tessella.API
const B=getfield(parentmodule(@__MODULE__),:APINativeBipartiteTests)
const C=getfield(parentmodule(@__MODULE__),:APINativeSelector4CompleteTests)
call(f,args...;kwargs...)=Base.invokelatest(f,args...;kwargs...)
bits(values)=reinterpret(UInt64,Float64.(values))
fingerprint(value)=begin io=IOBuffer();serialize(io,value);take!(io);end
@testset "public selector4 view rows retain consumed UIDs and number complete rows" begin
    cases=Any[]
    for (name,count) in (("native_selector4_views_primary.toml",224),
                         ("native_selector_view_algorithm_neighbors.toml",24))
        corpus=TOML.parsefile(joinpath(@__DIR__,"../artifacts",name))
        @test corpus["cases_count"]==length(corpus["cases"])==count
        for (raw,hash) in corpus["raw_bindings"]
            @test bytes2hex(sha256(read(joinpath(@__DIR__,"../artifacts",raw))))==hash
        end
        append!(cases,corpus["cases"])
    end
    @test length(cases)==248
    for case in cases
        @testset "view native $(case["fixture_id"]) $(case["view_kind"])" begin
        call(A.initialize);call(A.mesh.generate,0)
        view_tag=nothing
        try
            B.fixture(case)
            call(A.option,"Mesh.RecombinationAlgorithm",get(case,"algorithm",4))
            view_tag=call(A.view.add,"consumed_selector4_view")
            before=case["views_before"]
            rows=Dict(UInt64(tag)=>values for (tag,values) in zip(before["tags"],before["values"]))
            birth=UInt64.(case["birth_tags"])
            call(A.view.add_model_data,view_tag,0,call(A.model.get_current),case["view_kind"],
                birth,[rows[tag] for tag in birth],before["time"],before["components"])
            view=A.SESSION_VIEWS[Int(view_tag)]
            key=(0,case["view_kind"])
            @test view.model_contexts[key].extent==UInt64(case["logical_extent"])
            original=A.CURRENT[];cached=A.LAST_MESH[]
            class=cached===nothing ? nothing : call(A._cached_classification_locked,cached)
            snapshot=fingerprint((original,cached,view))
            detached=call(A._api_native_recombine_plan,original,cached,class,"native view detached";
                options=A.OPTIONS,initial_max_node_tag=A.NODE_TAG_MAX[],initial_max_element_tag=A.ELEMENT_TAG_MAX[])
            if haskey(case,"native_view_element_map")
                @test detached.native_view_element_tag_map==Dict(UInt64(a)=>UInt64(b) for (a,b) in case["native_view_element_map"])
                @test detached.native_view_node_tag_map==Dict(UInt64(a)=>UInt64(b) for (a,b) in case["native_view_node_map"])
                @test isempty(detached.element_tag_map)
            end
            call(A._native_recombine_views_plan,detached;renumber=!iszero(case["renumber"]))
            @test fingerprint((original,cached,view))==snapshot
            call(A.logger.start)
            result,out,err=C.capture(()->call(A.mesh.recombine))
            @test result===nothing
            logs=call(A.logger.get)
            if haskey(case,"algorithm")
                # Adjacent selector controls bind view numbering events and
                # tables. Their pre-existing mesher summary diagnostics have
                # separate coverage; retain full native logs in the artifact.
                view_event(message)=occursin("model data",message) || startswith(message,"Warning: Wrong source index")
                @test filter(view_event,logs)==filter(view_event,case["logger"])
            else
                @test C.timing_projection.(logs)==C.timing_projection.(case["logger"])
            end
            @test out==case["stdout"]
            @test err==case["stderr"]
            expected=case["views_after"]
            kinds,tags,values=call(A.view.get_model_data,view_tag,0)
            @test kinds==[expected["kind"]]
            @test tags==UInt64.(expected["tags"])
            @test bits(values)==bits([value for row in expected["values"] for value in row])
            @test view.model_contexts[key].counts==Int32[length(row) for row in expected["values"]]
            @test bits([view.model_contexts[key].time])==bits([expected["time"]])
            tags.=UInt64(99999);values.=99.
            @test call(A.view.get_model_data,view_tag,0)[2]==UInt64.(expected["tags"])
            B.state(case["after"])
            @test call(A.logger.get_last_error)==case["after"]["last_error"]
            @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==0
            call(A.logger.stop)
        finally
            A.NATIVE_LOGGER_ACTIVE[] && call(A.logger.stop)
            view_tag===nothing || call(A.view.remove,view_tag)
            call(A.finalize)
        end
        end
    end
end
end
end
