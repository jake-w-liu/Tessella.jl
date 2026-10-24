module NativeHybridPublicChecks
using Test,TOML,SHA,Tessella
const A=Tessella.API
const M=Tessella.Model
const artifacts=joinpath(@__DIR__,"..","artifacts")
const fixture=TOML.parsefile(joinpath(artifacts,"native_hybrid_public_4152_v1.toml"))
hexpoint(xs)=Float64[Core.eval(@__MODULE__,Meta.parse(x)) for x in xs]
function opened(f,source)
    A.initialize()
    try
        mktempdir() do directory
            path=joinpath(directory,"case.geo");write(path,source)
            A.open_geo!(path;mesh_dim=0);f()
        end
    finally
        A.finalize()
    end
end
@testset "public native hybrid original Box and current lower constructor history" begin
    for (key,filename) in (("direct_primary_sha256","tessella_ordinary_layers_hybrid_direct_native_primary_20261010_v1.json"),
            ("staged_primary_sha256","tessella_ordinary_layers_hybrid_current_constructor_primary_20261010_v2.json"))
        path=joinpath(@__DIR__,"..","..","validation","native_geometry",filename)
        @test bytes2hex(sha256(read(path)))==fixture[key]
    end
    for case in fixture["case"]
        @testset "$(case["name"]) lower phase $(case["lower_phase"])" begin
            println("HYBRID_PUBLIC_CASE ",case["name"]," lower_phase=",case["lower_phase"]);flush(stdout)
            opened(case["source"]) do
                case["lower_phase"]==2 && A.mesh.generate(2)
                A.mesh.generate(3)
                @test A.mesh.get_max_node_tag()==case["node_high_water"]
                @test A.mesh.get_max_element_tag()==case["element_high_water"]
                for product in case["product"]
                    dim,tag=product["owner"]
                    ids,xyz,uv=A.mesh.get_nodes(dim,tag,false,true)
                    @test ids==UInt64.(product["node_tags"])
                    expected_xyz=Float64[value for point in product["xyz_hex"] for value in hexpoint(point)]
                    @test isequal(xyz,expected_xyz)
                    @test isequal(uv,hexpoint(product["uv_hex"]))
                    families,labels,nodes=A.mesh.get_elements(dim,tag)
                    @test families==Int32[cell["msh"] for cell in get(product,"cell",[])]
                    @test labels==[UInt64.(cell["labels"]) for cell in get(product,"cell",[])]
                    @test nodes==[UInt64[uid for cell_nodes in cell["nodes"] for uid in cell_nodes] for cell in get(product,"cell",[])]
                end
                actual=A._with_model() do model;deepcopy(model.meshing.mesh_authority.mesh);end
                @test Tessella.MeshTypes.validate(actual).ok
                @test Set(block.msh for block in actual.blocks if Tessella.Elements.msh_dimension(block.msh)==3)==Set([4,7])
                saved=A.mesh.get_nodes()
                saved_cells=A.mesh.get_elements()
                A._with_model() do model
                    authority=model.meshing.mesh_authority
                    @test length(authority.surface_constructor_inputs)==6
                    for (owner,constructor) in authority.surface_constructor_inputs
                        current=M._model_mesh_product(authority,2,Int(owner[2]);caller="public hybrid constructor")
                        @test M._model_surface_constructor_binding(current,constructor,"public hybrid constructor")===nothing
                    end
                end
                @test isequal(A.mesh.get_nodes(),saved)
                @test A.mesh.get_elements()==saved_cells
            end
        end
    end
end
end
