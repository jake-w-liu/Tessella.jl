module APINativeSelector4GeoTests
using Test,Tessella,TOML,SHA,Serialization
if Sys.iswindows() && Sys.ARCH===:x86_64
const A=Tessella.API
const B=getfield(parentmodule(@__MODULE__),:APINativeBipartiteTests)
const C=getfield(parentmodule(@__MODULE__),:APINativeSelector4CompleteTests)
call(f,args...;kwargs...)=Base.invokelatest(f,args...;kwargs...)
fingerprint(value)=begin io=IOBuffer();serialize(io,value);take!(io);end
@testset "actual GEO selector4 parser/importer preserves native records and completed prefixes" begin
    directory=joinpath(@__DIR__,"../artifacts/native_selector4_geo_primary_v5")
    raw=joinpath(directory,"primary.json")
    projection=joinpath(directory,"primary.toml")
    @test bytes2hex(sha256(read(raw)))=="db7503299987f3a255515b15fab8ae4b169e85e7bcad2d1064503d0b9a92f4ba"
    corpus=TOML.parsefile(projection)
    @test corpus["raw_sha256"]==bytes2hex(sha256(read(raw)))
    @test length(corpus["cases"])==24
    for (index,case) in enumerate(corpus["cases"])
        @testset "GEO native case $index $(case["faces"])" begin
        call(A.initialize);call(A.mesh.generate,0)
        try
            source=abspath(joinpath(directory,case["source"]))
            @test bytes2hex(sha256(read(source)))==case["source_sha256"]
            mktempdir() do temporary
                prefix="General.Terminal=0;\nGeneral.Verbosity=5;\nGeneral.AbortOnError=$(case["abort"]);\n"*
                    "Mesh.RecombinationAlgorithm=4;\nMesh.Renumber=$(case["renumber"]);\nMerge \"$(replace(source,'\\'=>'/'))\";\n"
                baseline=joinpath(temporary,"source.geo")
                write(baseline,prefix)
                call(A.open_geo!,baseline)
                B.state(case["before"])
                if get(case,"prior_error",false)
                    if case["abort"]==0
                        @test call(A.mesh.get_node,999999)==(Float64[],Float64[],-1,-1)
                    else
                        @test_throws ErrorException call(A.mesh.get_node,999999)
                    end
                end
                original=A.CURRENT[];cached=A.LAST_MESH[]
                before=fingerprint((original,cached))
                path=joinpath(temporary,"recombine.geo")
                write(path,prefix*"RecombineMesh;\n")
                call(A.logger.start)
                failure,out,err=C.capture(()->call(A.open_geo!,path))
                expected=only(case["exceptions"])
                if isempty(expected)
                    @test failure===nothing
                else
                    @test failure isa ErrorException
                    @test failure!==nothing && sprint(showerror,failure)==expected
                end
                @test fingerprint((original,cached))==before
                @test C.timing_projection.(call(A.logger.get))==C.timing_projection.(case["logger"])
                @test out==case["stdout"]
                @test err==case["stderr"]
                B.state(case["after"])
                @test call(A.logger.get_last_error)==case["after"]["last_error"]
                @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==Int(get(case,"prior_error",false))+
                    count(startswith("Error: "),case["logger"])
                @test call(A.option,"Mesh.RecombinationAlgorithm")==4
                @test call(A.option,"Mesh.Renumber")==case["renumber"]
                call(A.logger.stop)
            end
        finally
            A.NATIVE_LOGGER_ACTIVE[] && call(A.logger.stop)
            call(A.finalize)
        end
        end
    end
end
end
end
