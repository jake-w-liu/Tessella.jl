module APINativeSelector4CompleteTests
using Test,Tessella,TOML,SHA,Serialization
if Sys.iswindows() && Sys.ARCH===:x86_64
const A=Tessella.API
const B=getfield(parentmodule(@__MODULE__),:APINativeBipartiteTests)
const F=getfield(parentmodule(@__MODULE__),:APINativeFullQuadGenerationTests)
call(f,args...;kwargs...)=Base.invokelatest(f,args...;kwargs...)
bits(values)=reinterpret(UInt64,Float64.(values))
fingerprint(value)=begin io=IOBuffer();serialize(io,value);take!(io);end
function capture(f)
    return mktemp() do outpath,outstream
        mktemp() do errpath,errstream
            result=redirect_stdout(outstream) do
                redirect_stderr(errstream) do
                    try f();nothing catch error;error end
                end
            end
            flush(outstream);flush(errstream)
            (result,read(outpath,String),read(errpath,String))
        end
    end
end
function timing_projection(text)
    expression=r"\(Wall ([^ ]+)s, CPU ([^ ]+)s\)"
    for match in eachmatch(expression,text)
        @test all(value->isfinite(value) && value>=0,
            parse.(Float64,match.captures))
    end
    replace(text,expression=>"(Wall <time>s, CPU <time>s)")
end
@testset "native selector4 generation deletion and diagnostic publication" begin
@testset "native generation4 preserves lower identities and consumed factory clocks" begin
    raw=joinpath(@__DIR__,"../artifacts/native_bipartite_generation_primary.json")
    projection=joinpath(@__DIR__,"../artifacts/native_bipartite_generation_primary.toml")
    @test bytes2hex(sha256(read(raw)))=="c5db8e5e7d17dfaaed8c0da112fa38863deaffa442e08a66cb0a257c87571bcc"
    @test bytes2hex(sha256(read(projection)))=="da2b118df53c413fdecd01636c99e7f4328dcd9a5939187043b5433ff255c8d0"
    corpus=TOML.parsefile(projection)
    @test corpus["raw_sha256"]==bytes2hex(sha256(read(raw)))
    @test length(corpus["cases"])==24
    for (index,case) in enumerate(corpus["cases"])
        @testset "generation4 native case $index" begin
        call(A.initialize);call(A.mesh.generate,0)
        try
            F.fixture(case);call(A.option,"General.Terminal",0)
            original=A.CURRENT[];cache=A.LAST_MESH[]
            snapshot=fingerprint((original,cache))
            word=Tessella.GmshLibm._win_x87_control_word()
            result,out,err=capture(()->call(A.mesh.generate,2))
            @test result===nothing
            @test fingerprint((original,cache))==snapshot
            @test call(A.mesh.get_elements)==(case["types"],
                [UInt64.(row) for row in case["tags"]],[UInt64.(row) for row in case["cells"]])
            @test Set(first(call(A.mesh.get_nodes)))==Set(UInt64(row["tag"]) for row in case["nodes"])
            for row in case["nodes"]
                xyz,uv,dimension,owner=call(A.mesh.get_node,row["tag"])
                @test bits(xyz)==bits(row["xyz"])
                @test bits(uv)==bits(row["parameters"])
                @test (dimension,owner)==Tuple(row["owner"])
                @test dimension!=2
            end
            @test call(A.mesh.get_max_node_tag)==case["max_node"]
            @test call(A.mesh.get_max_element_tag)==case["max_element"]
            @test call(A.logger.get_last_error)==case["last_error"]
            @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==0
            @test out==join(["1  "*label for label in ("Quadrangulation of face",
                "Transfer to Half Edge done","Initial labelling done","internal labelling done",
                "optimum labelling done","steiner points added","quads created","2 gmsh done")],"\n")*"\n"
            @test err==""
            @test Tessella.GmshLibm._win_x87_control_word()==word
        finally
            call(A.finalize)
        end
        end
    end
end
@testset "native logger and terminal admit selector4 events at their own verbosity gates" begin
    raw=joinpath(@__DIR__,"../artifacts/native_selector4_diagnostics_primary.json")
    projection=joinpath(@__DIR__,"../artifacts/native_selector4_diagnostics_primary.toml")
    @test bytes2hex(sha256(read(raw)))=="95435b668d54b7498200923ff0c04ed7a006b73c4dff42c7a03d287aa38b55d2"
    @test bytes2hex(sha256(read(projection)))=="2ae5ce433d94c3a03964ed65bb06ee7aa19eac5b5ed8c3c0d4a8df7b0369620e"
    corpus=TOML.parsefile(projection)
    @test corpus["raw_sha256"]==bytes2hex(sha256(read(raw)))
    @test length(corpus["cases"])==48
    for (index,case) in enumerate(corpus["cases"])
        @testset "diagnostic native case $index" begin
        call(A.initialize);call(A.mesh.generate,0)
        try
            B.fixture(case)
            call(A.option,"General.Terminal",case["terminal"])
            call(A.option,"General.Verbosity",case["verbosity"])
            @test !A.NATIVE_LOGGER_ACTIVE[]
            call(A.logger.start)
            result,out,err=capture(()->call(A.mesh.recombine))
            expected=only(case["exceptions"])
            if isempty(expected)
                @test result===nothing
            else
                @test result isa ErrorException
                @test sprint(showerror,result)==expected
            end
            messages=call(A.logger.get)
            @test timing_projection.(messages)==case["logger"]
            @test timing_projection(out)==case["stdout"]
            @test err==case["stderr"]
            B.state(case["after"])
            @test call(A.logger.get_last_error)==case["after"]["last_error"]
            @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==Int(first(case["faces"])=="highorder")
            call(A.logger.stop)
            @test isempty(call(A.logger.get))
        finally
            A.NATIVE_LOGGER_ACTIVE[] && call(A.logger.stop)
            call(A.finalize)
        end
        end
    end
end
@testset "native logger process lifetime and copied buffer" begin
    raw=joinpath(@__DIR__,"../artifacts/native_logger_lifecycle_primary.json")
    @test bytes2hex(sha256(read(raw)))=="f1294673377d6423d4e172da452cf7e36350e50d0adff99b02599672772c19b4"
    call(A.initialize);call(A.option,"General.Terminal",0)
    call(A.mesh.generate,0)
    try
        call(A.logger.stop);call(A.logger.start);call(A.logger.start)
        @test call(A.logger.get)==["Warning: Logger already started - ignoring"]
        detached=call(A.logger.get);empty!(detached)
        @test call(A.logger.get)==["Warning: Logger already started - ignoring"]
        call(A.finalize);call(A.initialize);call(A.option,"General.Terminal",0)
        @test call(A.logger.get)==["Warning: Logger already started - ignoring"]
        @test call(A.logger.get_last_error)==""
        @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==0
        call(A.logger.stop);@test isempty(call(A.logger.get))
    finally
        A.NATIVE_LOGGER_ACTIVE[] && call(A.logger.stop)
        call(A.finalize)
    end
end
end
end
end
