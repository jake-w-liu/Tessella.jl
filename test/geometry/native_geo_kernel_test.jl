module NativeGeoKernelTests
using Tessella,Test
const A=Tessella.API
include("../../validation/native_geo_kernel/native_geo_kernel_fixture_v1.jl")
float_bits(value)=string(reinterpret(UInt64,Float64(value));base=16,pad=16)
function check_state(expected;dangling_word=false)
    @test float_bits(A.option("General.BoundingBoxSize"))==float_bits(expected["lc"])
    @test Set(A.model.get_entities())==Set((row["dimension"],row["owner"]) for row in expected["entities"])
    for row in expected["entities"]
        dim,owner=row["dimension"],row["owner"]
        tags,coordinates,parameters=A.mesh.get_nodes(dim,owner,false,true)
        @test tags==UInt64.(row["node_tags"])
        @test float_bits.(coordinates)==row["coordinate_bits"]
        @test float_bits.(parameters)==row["parameter_bits"]
        types,cell_tags,connectivity=A.mesh.get_elements(dim,owner)
        @test types==row["element_types"]
        @test cell_tags==[UInt64.(values) for values in row["element_tags"]]
        if dangling_word && dim==1 && owner==11
            # Preserved native repeats demonstrate one freed MVertex word.
            # Keep all other retained element incidence exact.
            @test length(connectivity)==length(row["connectivity"])
            @test connectivity[1][2:end]==UInt64.(row["connectivity"][1][2:end])
        else
            @test connectivity==[UInt64.(values) for values in row["connectivity"]]
        end
    end
end
@testset "Native built-in GEO versus live Model Curve authority" begin
    for (index,row) in enumerate(NATIVE_GEO_KERNEL_ROWS)
        @testset "row=$index $(row["origin"]) $(row["mutation"])" begin
            A.initialize()
            try
                A.option("Mesh.Renumber",0)
                mktemp() do path,stream
                    write(stream,row["geo"]);close(stream);A.open_geo!(path)
                end
                origin,mutation=row["origin"],row["mutation"]
                origin=="api" && A.mesh.set_transfinite_curve(11,4)
                A.mesh.generate(1);check_state(row["states"][1])
                if origin=="geo_internal_setter"
                    A.model.geo.mesh.set_transfinite_curve(11,7)
                    A.mesh.generate(1);check_state(row["states"][2])
                    A.mesh.set_transfinite_curve(11,9)
                    A.mesh.generate(1);check_state(row["states"][3])
                    A.model.geo.synchronize();A.mesh.generate(1)
                    check_state(row["states"][4])
                    continue
                end
                state=2
                if mutation=="add_line_sync"
                    A.model.geo.add_point(2.,.3,0.;tag=5)
                    A.model.geo.add_line(4,5;tag=21)
                elseif mutation=="curve_retag_sync"
                    A.model.set_tag(1,11,91)
                    check_state(row["states"][state]);state+=1
                elseif mutation=="point_retag_sync"
                    A.model.set_tag(0,1,101)
                    check_state(row["states"][state]);state+=1
                elseif mutation=="api_override_sync"
                    A.mesh.set_transfinite_curve(11,7)
                    A.mesh.generate(1)
                    check_state(row["states"][state]);state+=1
                end
                A.model.geo.synchronize()
                check_state(row["states"][state];dangling_word=
                    origin=="api" && mutation=="point_retag_sync")
                state+=1
                A.mesh.generate(1);check_state(row["states"][state])
                if mutation=="repeat_sync"
                    A.model.geo.synchronize();A.mesh.generate(1)
                    check_state(row["states"][state+1])
                end
            finally
                A.finalize()
            end
        end
    end
end
end
