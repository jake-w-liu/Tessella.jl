module TransfiniteNativeOrphanBirthTests
using Test,TOML,SHA,Tessella
using Tessella.MeshTypes:nnodes
const API=Tessella.API
const FIXTURE=joinpath(@__DIR__,"..","artifacts","native_compact_transfinite_orphan_births_4152_v1.toml")
bits(values)=[string(reinterpret(UInt64,Float64(value));base=16,pad=16) for value in values]

@testset "Recorded compact TF orphan births survive actual publication" begin
    native=TOML.parsefile(FIXTURE)
    @test native["gmsh_version"]=="4.15.2"
    @test length(native["case"])==5
    for case in native["case"]
        @testset "$(case["name"])" begin
            @test bytes2hex(sha256(codeunits(case["source"])))==case["source_sha256"]
            mktempdir() do directory
                path=joinpath(directory,"native_compact.geo")
                write(path,case["source"])
                # The original .geo geometry projection also carries unused
                # native logical rows, separately from the complete cache.
                execution=Tessella.execute_geo(path)
                @test nnodes(execution.mesh)==length(case["node_tags"])
                @test Tessella.validate(execution.mesh).ok
                authority=execution.model.meshing.mesh_authority
                actual=authority.mesh
                @test actual.entity_data.external_node_tags==UInt64.(case["node_tags"])
                @test bits(vec(actual.coords))==case["xyz_bits"]
                @test [parameter===nothing ? String[] : bits(parameter) for parameter in actual.entity_data.node_parametric]==case["parameter_bits"]
                @test actual.entity_data.node_entities==Tuple{Int,Int32}[(owner[1],Int32(owner[2])) for owner in case["node_owners"]]
                for uid in UInt64.(case["orphan_tags"])
                    row=only(findall(==(uid),actual.entity_data.external_node_tags))
                    @test actual.entity_data.node_entities[row]==(3,Int32(1))
                    @test all(block->!(Int32(row) in block.nodes),actual.blocks)
                    @test row in authority.carrier_rows
                end
                # A separate public generation starts from the exact same CAD
                # and options, with the final Mesh command deferred to the API.
                deferred=replace(case["source"],r"Mesh\s+3\s*;"=>"")
                @test deferred!=case["source"]
                write(path,deferred)
                API.initialize()
                try
                    API.open_geo!(path;mesh_dim=0)
                    API.mesh.generate(3)
                    ids,xyz,_=API.mesh.get_nodes()
                    @test ids==UInt64.(case["node_tags"])
                    @test bits(xyz)==case["xyz_bits"]
                    @test [bits(API.mesh.get_node(uid)[2]) for uid in ids]==case["parameter_bits"]
                    @test API.mesh.get_nodes(3,1,false,true)[1]==UInt64.(case["volume_node_tags"])
                    @test API.mesh.get_max_node_tag()==case["max_node_tag"]
                    @test API.mesh.get_max_element_tag()==case["max_element_tag"]
                    for uid in UInt64.(case["orphan_tags"])
                        point,parameter,dimension,owner=API.mesh.get_node(uid)
                        native_row=only(findall(==(uid),case["node_tags"]))
                        @test bits(point)==case["xyz_bits"][3native_row-2:3native_row]
                        @test isempty(parameter)
                        @test (dimension,owner)==(3,1)
                    end
                    for owner in unique(cell["owner"] for cell in case["native_cells"])
                        records=filter(cell->cell["owner"]==owner,case["native_cells"])
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
