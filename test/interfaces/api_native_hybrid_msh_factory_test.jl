module NativeHybridMshFactoryChecks
using Test,Tessella
const E=Tessella.Elements

@testset "MSH factory clocks reject invalid inputs" begin
    @test_throws ArgumentError E._mixed_msh_factory_context(true,0)
    @test_throws ArgumentError E._mixed_msh_factory_context(0,false)
    @test_throws ArgumentError E._mixed_msh_factory_context(-1,0)
    @test_throws ArgumentError E._mixed_msh_factory_context(0,big(typemax(UInt64))+1)
end
end
