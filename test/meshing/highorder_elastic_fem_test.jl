module HighOrderElasticTests
using Test,Tessella
const A=Tessella.API

@testset "Public elastic Poisson controls and boundary preservation" begin
    for (poisson,expected) in (
        (0.0,(
            Float64[0.20175685171910943, -0.5320220614384539, 0.0],
            Float64[0.14053908224064962, -0.04310945992554645, 0.0],
            Float64[0.2829488895169657, 0.1670132678276144, 0.0],
            Float64[-0.2787775764037904, 0.26382097950255884, 0.0],
            Float64[-0.25347445425899795, -0.6682642207031829, 0.0],
        )),
        (0.2,(
            Float64[0.20173595854019555, -0.5319348385167016, 0.0],
            Float64[0.1404843591192392, -0.043062207180750094, 0.0],
            Float64[0.2828303447155906, 0.16699882634023475, 0.0],
            Float64[-0.2787887160623528, 0.2637955620822512, 0.0],
            Float64[-0.2534039407533041, -0.6680901902845942, 0.0],
        )),
        (0.33,(
            Float64[0.20171200258373592, -0.53180403073948, 0.0],
            Float64[0.1404272155128989, -0.04298153194539482, 0.0],
            Float64[0.28268691854538697, 0.16698082048202523, 0.0],
            Float64[-0.2788906213464496, 0.26380289232468535, 0.0],
            Float64[-0.2533495779728241, -0.667994074251216, 0.0],
        )),
        (0.49,(
            Float64[0.2018250247328992, -0.5303706162729214, 0.0],
            Float64[0.14040321505746742, -0.04284194324930699, 0.0],
            Float64[0.2820644564622929, 0.16677276140885366, 0.0],
            Float64[-0.28081428179996115, 0.26415257398328157, 0.0],
            Float64[-0.25354843539618843, -0.6696500399334403, 0.0],
        )),
    )
        A.initialize()
        try
            A.option("General.Terminal",0)
            A.open(joinpath(@__DIR__,"..","fixtures","msh_high_order_elastic_disk.msh"))
            A.option("Mesh.HighOrderPoissonRatio",poisson)
            A.option("Mesh.HighOrderThresholdMin",.99)
            before=A.mesh.get_nodes()
            boundary=A.mesh.get_nodes(1,1,true)[1]
            A.mesh.optimize("HighOrderElastic",true,1)
            tags,coords,_=A.mesh.get_nodes()
            @test tags==before[1]
            @test any(coords .!= before[2])
            for tag in boundary
                k=findfirst(==(tag),tags)
                @test coords[3k-2:3k]==before[2][3k-2:3k]
            end
            for (indices,point) in zip((79:81,103:105,148:150,223:225,286:288),expected)
                @test coords[indices]≈point atol=1e-8 rtol=0
            end
        finally
            A.finalize()
        end
    end
end
end
