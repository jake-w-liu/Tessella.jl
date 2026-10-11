module NativeLineFallbackTests
using Tessella,Test
const A=Tessella.API
include("../../validation/native_geo_kernel/native_line_fallback_fixture_v1.jl")
include("../../validation/native_geo_kernel/native_line_maxima_fixture_v2.jl")
float_bits(value)=string(reinterpret(UInt64,Float64(value));base=16,pad=16)
function check_state(expected,maxima)
    @test A.mesh.get_max_node_tag()==maxima[1]
    @test A.mesh.get_max_element_tag()==maxima[2]
    @test float_bits(A.option("General.BoundingBoxSize"))==expected["lc_bits"]
    @test Set(A.model.get_entities())==Set((row["dimension"],row["owner"]) for row in expected["entities"])
    for row in expected["entities"]
        dim,owner=row["dimension"],row["owner"]
        tags,coordinates,parameters=A.mesh.get_nodes(dim,owner,false,true)
        @test tags==UInt64.(row["node_tags"])
        @test float_bits.(coordinates)==row["coordinate_bits"]
        @test float_bits.(parameters)==row["parameter_bits"]
        types,tags,nodes=A.mesh.get_elements(dim,owner)
        @test types==row["element_types"]
        @test tags==[UInt64.(values) for values in row["element_tags"]]
        @test nodes==[UInt64.(values) for values in row["connectivity"]]
    end
end
@testset "Native ordinary Line unsized and sized endpoint fallbacks" begin
    for (index,row) in enumerate(NATIVE_LINE_FALLBACK_ROWS)
        @testset "row=$index $(row["recipe"]) R=$(row["renumber"]) SOL=$(row["linear"])" begin
            A.initialize()
            try
                A.option("Mesh.Renumber",row["renumber"])
                A.option("Mesh.SecondOrderLinear",row["linear"])
                mktemp() do path,stream
                    write(stream,row["geo"]);close(stream);A.open_geo!(path)
                end
                A.mesh.generate(1);check_state(row["states"][1],NATIVE_LINE_FALLBACK_MAXIMA[index][1])
                A.mesh.set_order(3);check_state(row["states"][2],NATIVE_LINE_FALLBACK_MAXIMA[index][2])
                A.mesh.generate(1);check_state(row["states"][3],NATIVE_LINE_FALLBACK_MAXIMA[index][3])
            finally
                A.finalize()
            end
        end
    end
end
end
