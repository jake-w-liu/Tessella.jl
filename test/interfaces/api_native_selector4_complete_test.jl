module APINativeSelector4CompleteTests
using Test,Tessella
const A=Tessella.API
call(f,args...;kwargs...)=Base.invokelatest(f,args...;kwargs...)
@testset "native logger process lifetime and copied buffer" begin
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
end # module
