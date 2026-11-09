module OCCFuseImageChecks
using Test,Tessella,TOML
const A=Tessella.API
const CONTROL=TOML.parsefile(joinpath(@__DIR__,"..","fixtures","native_occ_fuse_images.toml"))
dimtags(rows)=Tuple{Int,Int}[(row[1],row[2]) for row in rows]

@testset "Native fuse images match public geometry controls" begin
    @test CONTROL["version"]=="4.15.2"
    @test length(CONTROL["cases"])==44
    for case in CONTROL["cases"]
        @testset "$(case["name"])" begin
            A.initialize()
            try
                A.option("General.Terminal",0)
                occ=A.model.occ
                occ.add_box(0.,0.,0.,1.,1.,1.;tag=1)
                occ.add_box(case["box"]...;tag=2)
                out,images=occ.fuse([(3,1)],[(3,2)];tag=case["tag"],
                    removeObject=case["remove_object"],removeTool=case["remove_tool"])
                @test out==dimtags(case["out"])
                @test images==dimtags.(case["images"])
                occ.synchronize()
                @test all(dt->dt in A.model.get_entities(3),Iterators.flatten(images))
            finally
                A.finalize()
            end
        end
    end
end
end
