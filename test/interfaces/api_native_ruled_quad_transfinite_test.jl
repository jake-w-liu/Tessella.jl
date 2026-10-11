module NativeRuledQuadTransfiniteTests
using Test,TOML,SHA,Tessella
const API=Tessella.API
const FIXTURE=joinpath(@__DIR__,"..","artifacts","native_ruled_quad_transfinite_4152_v1.toml")
bits(values)=[string(reinterpret(UInt64,Float64(x));base=16,pad=16) for x in values]

@testset "Native ruled quad TF publication preserves coordinates and births" begin
    native=TOML.parsefile(FIXTURE)
    @test native["gmsh_version"]=="4.15.2"
    @test length(native["case"])==64
    @test bytes2hex(sha256(read(joinpath(@__DIR__,"..","artifacts","native_ruled_quad_transfinite_primary_4152_v1.json"))))==native["primary_sha256"]
    for case in native["case"]
        @testset "$(case["name"])" begin
            @test bytes2hex(sha256(codeunits(case["source"])))==case["source_sha256"]
            mktempdir() do directory
                path=joinpath(directory,"ruled_quad.geo")
                deferred=replace(case["source"],r"Mesh\s+2\s*;"=>"")
                @test deferred!=case["source"]
                write(path,deferred)
                API.initialize()
                try
                    API.open_geo!(path;mesh_dim=0)
                    API.mesh.generate(2)
                    ids,xyz,_=API.mesh.get_nodes()
                    @test ids==UInt64.(case["node_tags"])
                    @test bits(xyz)==reduce(vcat,case["xyz_bits"])
                    @test API.mesh.get_max_node_tag()==case["max_node_tag"]
                    @test API.mesh.get_max_element_tag()==case["max_element_tag"]
                    for (row,uid) in enumerate(ids)
                        _,parameters,dim,tag=API.mesh.get_node(uid)
                        @test bits(parameters)==case["parameter_bits"][row]
                        @test [dim,tag]==case["owners"][row]
                    end
                    for owner in unique(cell["owner"] for cell in case["cell"])
                        records=filter(cell->cell["owner"]==owner,case["cell"])
                        kinds,labels,nodes=API.mesh.get_elements(owner...)
                        @test kinds==[cell["msh"] for cell in records]
                        @test labels==[UInt64.(cell["labels"]) for cell in records]
                        @test nodes==[UInt64.(reduce(vcat,cell["nodes"])) for cell in records]
                    end
                    before=deepcopy((API.mesh.get_nodes(),API.mesh.get_elements(),
                        API.mesh.get_max_node_tag(),API.mesh.get_max_element_tag()))
                    API.mesh.set_size([(0,1)],0.13)
                    @test isequal(before,(API.mesh.get_nodes(),API.mesh.get_elements(),
                        API.mesh.get_max_node_tag(),API.mesh.get_max_element_tag()))
                finally
                    API.finalize()
                end
            end
        end
    end
end
end
