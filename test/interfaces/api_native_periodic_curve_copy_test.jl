module NativePeriodicCurveCopyV1
using Tessella,Test
const A=Tessella.API
include(joinpath(@__DIR__,"..","..","validation","native_geometry","native_periodic_curve_copy_literal_v1.jl"))

@testset "Native affine periodic Curve copies project actual master vertices" begin
    for row in NATIVE_PERIODIC_CURVE_COPY_V1
        A.initialize()
        try
            A.option("Mesh.Renumber",row["renumber"])
            A.option("Mesh.ElementOrder",row["order"])
            for (tag,point) in enumerate(row["points"])
                A.model.add_point(point...;meshSize=1.,tag)
            end
            A.model.add_line(1,2;tag=1)
            A.model.add_line(3,4;tag=2)
            for curve in (1,2)
                A.mesh.set_transfinite_curve(curve,row["count"])
            end
            A.mesh.set_periodic(1,[2],[1],row["affine"])
            A.mesh.generate(1)
            tags,xyz,uv=A.mesh.get_nodes()
            @test tags==row["node_tags"]
            @test reinterpret(UInt64,xyz)==reinterpret(UInt64,Float64.(row["coordinates"]))
            @test uv==row["parameters"]
            types,etags,enodes=A.mesh.get_elements()
            @test types==row["element_types"]
            @test etags==row["element_tags"]
            @test enodes==row["element_nodes"]
            for curve in (1,2)
                expected=row["owned"][string(curve)]
                owned_tags,owned_xyz,owned_uv=A.mesh.get_nodes(1,curve)
                @test owned_tags==expected["node_tags"]
                @test reinterpret(UInt64,owned_xyz)==reinterpret(UInt64,Float64.(expected["coordinates"]))
                @test reinterpret(UInt64,owned_uv)==reinterpret(UInt64,Float64.(expected["parameters"]))
            end
            master,types,mtypes,keys,mkeys,coords,mcoords=A.mesh.get_periodic_keys(row["msh"],"Lagrange",2)
            @test master==row["master"]
            @test types==row["types"]
            @test mtypes==row["master_types"]
            @test keys==row["keys"]
            @test mkeys==row["master_keys"]
            @test reinterpret(UInt64,coords)==reinterpret(UInt64,Float64.(row["periodic_coordinates"]))
            @test reinterpret(UInt64,mcoords)==reinterpret(UInt64,Float64.(row["periodic_master_coordinates"]))
        finally
            A.finalize()
        end
    end
end
end
