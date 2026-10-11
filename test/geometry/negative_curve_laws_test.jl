module NegativeCurveLawsTests
using Tessella,Test
const A=Tessella.API
include("../../validation/native_curve_events/negative_curve_laws_fixture_v1.jl")
float_bits(value)=string(reinterpret(UInt64,Float64(value));base=16,pad=16)
function check_owners(expected)
    @test Set(A.model.get_entities())==Set((row["dimension"],row["owner"]) for row in expected)
    for row in expected
        dim,owner=row["dimension"],row["owner"]
        tags,coordinates,parameters=A.mesh.get_nodes(dim,owner,false,true)
        @test tags==UInt64.(row["node_tags"])
        @test float_bits.(coordinates)==row["coordinate_bits"]
        @test float_bits.(parameters)==row["parameter_bits"]
        types,cell_tags,connectivity=A.mesh.get_elements(dim,owner)
        @test types==row["element_types"]
        @test cell_tags==[UInt64.(values) for values in row["element_tags"]]
        @test connectivity==[UInt64.(values) for values in row["connectivity"]]
        for (slot,tag) in enumerate(tags)
            point,parameter,dimension,classification=A.mesh.get_node(tag)
            @test (dimension,classification)==(dim,owner)
            @test float_bits.(point)==row["coordinate_bits"][3slot-2:3slot]
            @test float_bits.(parameter)==(dim==1 ? [row["parameter_bits"][slot]] : String[])
        end
    end
end
@testset "Native signed Curve laws and preserved positive/API controls" begin
    for (index,row) in enumerate(NEGATIVE_CURVE_LAW_ROWS)
        @testset "row=$index $(row["recipe"]) $(row["scope"]) $(row["law"])" begin
            A.initialize()
            try
                A.option("Mesh.Renumber",0)
                execution=mktemp() do path,stream
                    write(stream,row["geo"]);close(stream);A.open_geo!(path)
                end
                if row["scope"]=="api_negative_control"
                    A.mesh.set_transfinite_curve(11,row["count"],row["law"],row["coefficient"])
                elseif row["scope"]=="raw_geo_negative" && row["recipe"] in ("line","circle")
                    # The captured carrier remains independent of actual native births.
                    model=execution.model
                    spec=model.meshing.transfinite_curves[11]
                    params=Tessella.Model._model_curve_transfinite_params(
                        model,11,0.0,1.0,spec,"negative Curve carrier check")
                    @test params==[k/(row["count"]-1) for k in 0:row["count"]-1]
                end
                for expected in row["states"]
                    A.mesh.generate(1)
                    check_owners(expected)
                end
            finally
                A.finalize()
            end
        end
    end
end
end
