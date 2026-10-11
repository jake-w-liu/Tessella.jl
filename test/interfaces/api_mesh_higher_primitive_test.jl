module HigherPrimitiveCompatibilityTests
using Test,Tessella
using Tessella.MeshTypes: mesh_crc,ntets,validate
const API=Tessella.API

@testset "Higher generation preserves supported analytic native Volume encodings" begin
    try
        for kind in (:sphere,:cylinder,:cone)
            API.finalize();API.initialize()
            if kind===:sphere
                API.model.add_sphere(0,0,0,1;tag=1)
            elseif kind===:cylinder
                API.model.add_cylinder(0,0,0,0,0,2,1;tag=1)
            else
                API.model.add_cone(0,0,0,0,0,2,1,.5;tag=1)
            end
            mesh=API.mesh.generate(3)
            @test validate(mesh).ok
            @test ntets(mesh)==(kind===:sphere ? 152 : 96)
            @test mesh_crc(mesh).sha==(kind===:sphere ?
                "86c41ffe7a27274dadcee65064f657a15f8d839807db1785bdd3c4c3c6d944b3" :
                "735491540012b556296701a5f4ab11b7b12552441fbb732b385eb67ae02286c4")
            types,tags,connections=API.mesh.get_elements()
            nodes=API.mesh.get_nodes()[1]
            @test allunique(nodes) && allunique(vcat(tags...))
            @test Set(vcat(connections...)) ⊆ Set(nodes)
        end
    finally
        API.finalize()
    end
end
end
