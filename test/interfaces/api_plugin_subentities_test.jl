# `plugin` MeshSubEntities: vertex/edge sub-entity enumeration on a
# physically-grouped surface mesh, upstream option spellings (`Ouput*`),
# silent no-op dimensions, and the empty-group diagnostic.
module PluginSubEntitiesChecks
using Test, Tessella
const A=Tessella.API

@testset "MeshSubEntities plugin" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.model.add("m")
        A.model.add_box(0,0,0,1,1,1)
        A.option("Mesh.MeshSizeMin",0.3)
        A.option("Mesh.MeshSizeMax",0.3)
        surface_tags=[t for (_,t) in A.model.get_entities(2)]
        curve_tags=[t for (_,t) in A.model.get_entities(1)]
        A.model.add_physical_group(2,surface_tags;tag=1)
        A.mesh.generate(2)

        points0=length(A.model.get_entities(0))
        curves0=length(A.model.get_entities(1))
        faces0=length(A.model.get_entities(2))

        # OuputDimension 0: unique mesh vertices on the grouped surfaces
        # spawn discrete point entities; corner vertices already realized
        # by the box's point entities reuse those owners. Either way each
        # owning entity joins the output physical group.
        A.plugin.set_number("MeshSubEntities","InputDimension",2)
        A.plugin.set_number("MeshSubEntities","InputPhysicalGroup",1)
        A.plugin.set_number("MeshSubEntities","OuputDimension",0)
        A.plugin.set_number("MeshSubEntities","OuputPhysicalGroup",2)
        A.plugin.run("MeshSubEntities")
        @test length(A.model.get_entities(0))>points0
        @test !isempty(A.model.get_entities_for_physical_group(0,2))

        # OuputDimension 1: unique element edges spawn discrete curve
        # entities (boundary edges may already be realized on the box's
        # curves; interior edges create new ones).
        A.plugin.set_number("MeshSubEntities","OuputDimension",1)
        A.plugin.set_number("MeshSubEntities","OuputPhysicalGroup",3)
        A.plugin.run("MeshSubEntities")
        @test length(A.model.get_entities(1))>curves0
        @test !isempty(A.model.get_entities_for_physical_group(1,3))

        # OuputDimension 2 is a silent no-op upstream: no entities are
        # created and the output group never materializes.
        A.plugin.set_number("MeshSubEntities","OuputDimension",2)
        A.plugin.set_number("MeshSubEntities","OuputPhysicalGroup",4)
        A.plugin.run("MeshSubEntities")
        @test length(A.model.get_entities(2))==faces0
        @test_throws ArgumentError A.model.get_entities_for_physical_group(2,4)

        # Edge requests on 1-D input also produce nothing: line elements
        # have no proper 1-D sub-entities.
        A.model.add_physical_group(1,curve_tags;tag=5)
        A.plugin.set_number("MeshSubEntities","InputDimension",1)
        A.plugin.set_number("MeshSubEntities","InputPhysicalGroup",5)
        A.plugin.set_number("MeshSubEntities","OuputDimension",1)
        A.plugin.set_number("MeshSubEntities","OuputPhysicalGroup",6)
        curves_before=length(A.model.get_entities(1))
        A.plugin.run("MeshSubEntities")
        @test length(A.model.get_entities(1))==curves_before
        @test_throws ArgumentError A.model.get_entities_for_physical_group(1,6)

        # An empty input physical group reports the upstream diagnostic
        # and produces no output.
        A.logger.start()
        A.plugin.set_number("MeshSubEntities","InputDimension",2)
        A.plugin.set_number("MeshSubEntities","InputPhysicalGroup",99)
        A.plugin.set_number("MeshSubEntities","OuputDimension",0)
        A.plugin.set_number("MeshSubEntities","OuputPhysicalGroup",7)
        A.plugin.run("MeshSubEntities")
        msgs=A.logger.get()
        A.logger.stop()
        @test any(m->occursin("Physical group 99 (dimension 2) is empty",m),
                  msgs)
        @test_throws ArgumentError A.model.get_entities_for_physical_group(0,7)
    finally
        A.finalize()
    end
end

end # module PluginSubEntitiesChecks
