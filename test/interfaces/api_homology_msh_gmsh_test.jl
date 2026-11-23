module ApiHomologyMshGmshTest

using Test, Tessella

const API = Tessella.API
const FIXTURE = joinpath(@__DIR__, "..", "fixtures", "gmsh_square_physical80.msh")

# Homology on a mesh read from MSH (written by Gmsh 4.15.2: a square with
# physical surface 7 and physical curve 80). Gmsh numbers a chain's physical
# group after the largest physical number of any dimension (81 here);
# Tessella started dimension 0 at 1. The basis point is Gmsh's element 17.
@testset "homology chains on an MSH mesh take Gmsh's physical numbers" begin
    try
        API.initialize(); API.option("General.Terminal", 0)
        API.open(FIXTURE)
        API.mesh.add_homology_request("Homology", [7], [80])
        API.mesh.add_homology_request("Homology", [7], Int[])
        API.mesh.compute_homology()
        @test API.model.get_physical_groups() == [(0, 81), (1, 80), (2, 7)]
        @test API.model.get_physical_name(0, 81) == "H_0{7}1"
        entities = API.model.get_entities_for_physical_group(0, 81)
        @test length(entities) == 1
        @test API.mesh.get_elements(0, entities[1])[2] == [[17]]
    finally
        API.finalize()
    end
end

# Gmsh writes a $Nodes total above the blocks' sum when Mesh.SaveAll is off
# and entities are left out (its own reader only sizes an array with it):
# such a file reads, with the nodes it holds.
@testset "an MSH4 node total above the nodes written is read" begin
    path = joinpath(mktempdir(), "declared.msh")
    text = read(FIXTURE, String)
    write(path, replace(text, "9 12 1 12" => "9 14 1 12"))
    try
        API.initialize(); API.option("General.Terminal", 0)
        API.open(path)
        @test length(API.mesh.get_nodes()[1]) == 12
        @test sum(length, API.mesh.get_elements()[2]) == 16
    finally
        API.finalize()
    end
end

end
