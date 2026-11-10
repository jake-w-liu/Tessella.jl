module NativeTri3MappedPointProbeV1
using Tessella,Test,SHA,TOML,Serialization
const A=Tessella.API
include(joinpath(@__DIR__,"native_tri3_mapped_point_fixtures_v1.jl"))
bits(xs)=[string(reinterpret(UInt64,Float64(x));base=16,pad=16) for x in xs]
digest(x)=let io=IOBuffer();serialize(io,x);bytes2hex(sha256(take!(io)));end
function run(output)
    records=Dict{String,Any}[]
    @testset "Native actual raw Tri3 original-coordinate mapped sums" begin
        for row in NATIVE_TRI3_MAPPED60
            @testset "$(row["recipe"])/$(row["repeat"])" begin
                A.initialize();A.model.add("raw literal Tri3")
                try
                    A.model.add_discrete_entity(2,50)
                    A.mesh.add_nodes(2,50,[2001,2002,2003],Float64.(row["input_coordinates"]))
                    A.mesh.add_elements_by_type(50,2,[9001],[2001,2002,2003])
                    model=A.CURRENT[];meshing=model.meshing;authority=meshing.mesh_authority
                    snapshot=digest((model,A.OPTIONS,A.SESSION_VIEWS,A.NODE_TAG_MAX[],A.ELEMENT_TAG_MAX[]))
                    @test A.LAST_MESH[]===nothing
                    local_coordinates=Float64.(row["local_coordinates"])
                    scalar=A.mesh.get_jacobian(9001,local_coordinates)
                    bulk=A.mesh.get_jacobians(2,local_coordinates,50)
                    @test all(bits(a)==bits(b) for (a,b) in zip(scalar,bulk))
                    @test bits(scalar[1])==row["jacobian_bits"]
                    @test bits(scalar[2])==row["determinant_bits"]
                    @test bits(scalar[3])==row["mapped_point_bits"]
                    components,values,orientations=A.mesh.get_basis_functions(2,local_coordinates,"Lagrange")
                    @test components==1 && orientations==1
                    @test bits(values)==row["values_bits"]
                    tags,coordinates,parameters=A.mesh.get_nodes()
                    @test Int.(tags)==row["node_tags"]
                    @test bits(coordinates)==row["node_coordinate_bits"]
                    @test bits(parameters)==row["node_parameter_bits"]
                    types,labels,nodes=A.mesh.get_elements()
                    @test Int.(types)==row["element_types"]
                    @test [Int.(r) for r in labels]==row["element_tags"]
                    @test [Int.(r) for r in nodes]==row["element_nodes"]
                    @test A.mesh.get_max_node_tag()==row["max_node_tag"]
                    @test A.mesh.get_max_element_tag()==row["max_element_tag"]
                    @test A.CURRENT[]===model && model.meshing===meshing && meshing.mesh_authority===authority
                    @test A.LAST_MESH[]===nothing && A.LAST_MESH_CLASS[]===nothing
                    @test digest((model,A.OPTIONS,A.SESSION_VIEWS,A.NODE_TAG_MAX[],A.ELEMENT_TAG_MAX[]))==snapshot
                    push!(records,Dict("recipe"=>row["recipe"],"repeat"=>row["repeat"],"mapped_sha256"=>digest(bits(scalar[3]))))
                finally;A.finalize();end
            end
        end
    end
    open(output,"w") do io;TOML.print(io,Dict("julia_version"=>string(VERSION),"rows"=>records,
        "native_bits"=>true,"callbacks_added"=>0));end
end
end
