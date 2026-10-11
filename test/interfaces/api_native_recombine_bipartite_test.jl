module APINativeBipartiteTests
using Test,Tessella,TOML,SHA,Serialization
if Sys.iswindows() && Sys.ARCH===:x86_64
const A=Tessella.API
call(f,args...;kwargs...)=Base.invokelatest(f,args...;kwargs...)
bits(values)=reinterpret(UInt64,Float64.(values))
fingerprint(value)=begin io=IOBuffer();serialize(io,value);take!(io);end
function state(expected)
    @test call(A.mesh.get_elements)==(expected["types"],
        [UInt64.(row) for row in expected["tags"]],[UInt64.(row) for row in expected["cells"]])
    nodes=get(expected,"nodes",Any[])
    @test Set(first(call(A.mesh.get_nodes)))==Set(UInt64(row["tag"]) for row in nodes)
    for row in nodes
        xyz,uv,dimension,owner=call(A.mesh.get_node,row["tag"])
        @test bits(xyz)==bits(row["xyz"])
        @test bits(uv)==bits(row["parameters"])
        @test (dimension,owner)==Tuple(row["owner"])
    end
    for (tags,owners) in zip(expected["tags"],expected["element_owners"]), (tag,owner) in zip(tags,owners)
        record=call(A.mesh.get_element,tag)
        @test (record[3],record[4])==Tuple(owner)
    end
    @test call(A.mesh.get_max_node_tag)==expected["max_node"]
    @test call(A.mesh.get_max_element_tag)==expected["max_element"]
end
function fixture(case)
    shapes=get(case,"faces",["pair"]);face=get(case,"face_tag",7)
    if get(case,"cad",false)
        points=((0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.))
        for (tag,point) in enumerate(points);call(A.model.add_point,point...;tag,meshSize=10.);end
        for tag in 1:4;call(A.model.add_line,tag,mod1(tag+1,4);tag);end
        call(A.model.add_curve_loop,[1,2,3,4];tag=face)
        call(A.model.add_plane_surface,[face];tag=face)
        get(case,"volume",false) && call(A.model.add_discrete_entity,3,47)
    else
        entities=Set{Tuple{Int,Int}}((2,face+10(index-1)) for index in eachindex(shapes))
        for row in get(case["before"],"nodes",Any[]);push!(entities,Tuple(row["owner"]));end
        for owners in case["before"]["element_owners"],owner in owners;push!(entities,Tuple(owner));end
        for entity in sort!(collect(entities));call(A.model.add_discrete_entity,entity...);end
    end
    before=case["before"]
    for row in get(before,"nodes",Any[])
        call(A.mesh.add_nodes,row["owner"]...,[row["tag"]],row["xyz"],row["parameters"])
    end
    for (kind,tags,cells,owners) in zip(before["types"],before["tags"],before["cells"],before["element_owners"])
        width=div(length(cells),length(tags))
        for index in eachindex(tags)
            call(A.mesh.add_elements_by_type,owners[index][2],kind,[tags[index]],cells[(index-1)*width+1:index*width])
        end
    end
    call(A.option,"General.Terminal",0)
    call(A.option,"General.AbortOnError",get(case,"abort",2))
    call(A.option,"Mesh.Renumber",get(case,"renumber",0))
    call(A.option,"Mesh.RecombinationAlgorithm",4)
end
@testset "native public experimental selector ordered transfer and lifecycle" begin
    raw=joinpath(@__DIR__,"../artifacts/native_recombine_bipartite_primary.json")
    projection=joinpath(@__DIR__,"../artifacts/native_recombine_bipartite_primary.toml")
    @test bytes2hex(sha256(read(raw)))=="874ceb7012b0ad6887c05b6fc341321bde96d1a9f11d17d234123d6aa019cc52"
    @test bytes2hex(sha256(read(projection)))=="140f59fdee7bb95a03880e6964cf027fb4d5714e31204c0b103037ed11cfbdf0"
    @test bytes2hex(sha256(read(joinpath(@__DIR__,"../artifacts/native_selector4_prior_lookup_primary.json"))))==
        "45a60ccde8ea9ef34667d7783887875d99ec6a05cff9586be22630e7a29f47db"
    corpus=TOML.parsefile(projection)
    @test corpus["raw_sha256"]==bytes2hex(sha256(read(raw)))
    @test length(corpus["cases"])==40
    for (index,case) in enumerate(corpus["cases"])
        @testset "native experimental control $index $(get(case,"faces",Any[]))" begin
        call(A.initialize);call(A.mesh.generate,0)
        try
            fixture(case);state(case["before"])
            if get(case,"prior_error",false)
                if iszero(get(case,"abort",2))
                    @test call(A.mesh.get_node,999999)==(Float64[],Float64[],-1,-1)
                else
                    @test_throws ErrorException call(A.mesh.get_node,999999)
                end
            end
            current=A.CURRENT[];cached=A.LAST_MESH[]
            snapshot=fingerprint((current,cached,A.NODE_TAG_MAX[],A.ELEMENT_TAG_MAX[],
                A.NATIVE_LOOKUP_QUERY_ERRORS[],A.NATIVE_MSG_LAST_ERROR[]))
            class=cached===nothing ? nothing : call(A._cached_classification_locked,cached)
            plan=call(A._api_native_recombine_plan,current,cached,class,"native experimental detached";
                options=A.OPTIONS,initial_max_node_tag=A.NODE_TAG_MAX[],initial_max_element_tag=A.ELEMENT_TAG_MAX[],
                skip_recombination=!iszero(A.OPTIONS["General.AbortOnError"]) && A.NATIVE_LOOKUP_QUERY_ERRORS[]>0)
            @test fingerprint((current,cached,A.NODE_TAG_MAX[],A.ELEMENT_TAG_MAX[],
                A.NATIVE_LOOKUP_QUERY_ERRORS[],A.NATIVE_MSG_LAST_ERROR[]))==snapshot
            if get(case,"repeat",1)==1
                diagnostics=plan===nothing ? String[] : [string(kind===:warning ? "Warning: " : "Error: ",message)
                    for (kind,message) in get(plan,:native_recombine_events,()) if kind in (:warning,:error)]
                @test diagnostics==case["diagnostics"]
            end
            word=Tessella.GmshLibm._win_x87_control_word()
            failures=String[]
            for step in 1:get(case,"repeat",1)
                failure=try call(A.mesh.recombine);nothing catch error;error;end
                expected=case["exceptions"][step]
                if isempty(expected)
                    @test failure===nothing
                else
                    @test failure isa ErrorException
                    @test failure!==nothing && sprint(showerror,failure)==expected
                end
                push!(failures,failure===nothing ? "" : sprint(showerror,failure))
            end
            @test failures==case["exceptions"]
            state(case["after"])
            @test call(A.logger.get_last_error)==case["after"]["last_error"]
            @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==Int(get(case,"prior_error",false))+
                count(startswith("Error: "),case["diagnostics"])
            @test Tessella.GmshLibm._win_x87_control_word()==word
        finally
            call(A.finalize)
        end
        end
    end
end
end
end
