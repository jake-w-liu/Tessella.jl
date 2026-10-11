module GeoMeshCarrierCellTests
using Tessella,Test,TOML,SHA
using Tessella.MeshTypes: Mesh
using Tessella.Elements: MixedMesh,ElementBlock
const API=Tessella.API
const GEO=Tessella.GeoExec
const Model=Tessella.Model

function square_source(order,renumber)
    return """
    Mesh.Renumber=$renumber; Mesh.ElementOrder=$order; Mesh.Smoothing=0;
    Point(1)={0,0,0,1}; Point(2)={1,0,0,1}; Point(3)={1,1,0,1}; Point(4)={0,1,0,1};
    Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
    Curve Loop(1)={1,2,3,4}; Plane Surface(1)={1};
    Transfinite Curve{1:4}=2; Transfinite Surface{1};
    """
end

function check_literal_carrier(result,row)
    authority=result.model.meshing.mesh_authority
    geometry=result.mesh;data=authority.mesh.entity_data
    expected=Dict{UInt64,Tuple{Int,Vector{UInt64}}}()
    for entity in row["entities"],(msh,tags,flat) in
            zip(entity["types"],entity["elements"],entity["connections"])
        width=Tessella.Elements.msh_num_nodes(msh)
        for (column,tag) in enumerate(tags)
            expected[UInt64(tag)]=(msh,UInt64.(flat[width*(column-1)+1:width*column]))
        end
    end
    blocks=Model._model_mesh_carrier_blocks(geometry)
    @test length(authority.carrier_cells)==length(blocks)
    @test authority.carrier===geometry
    @test authority.mesh.coords[:,authority.carrier_rows]==geometry.coords
    visible=UInt64[]
    for ((msh,cells),tags) in zip(blocks,authority.carrier_cells)
        @test length(tags)==size(cells,2)
        for (column,tag) in enumerate(tags)
            @test haskey(expected,tag)
            haskey(expected,tag) || continue
            actual_msh,connections=expected[tag]
            @test Tessella.Elements.msh_family(actual_msh)==Tessella.Elements.msh_family(msh)
            public=data.external_node_tags[authority.carrier_rows[cells[:,column]]]
            @test public==connections[1:length(public)]
            push!(visible,tag)
        end
    end
    expected_tags=sort!(UInt64[tag for (tag,(msh,_)) in expected if msh!=15])
    @test sort!(visible)==expected_tags
    return authority
end

@testset "Native carrier cells bind literal primary element identities" begin
    primary=TOML.parsefile(joinpath(@__DIR__,"../artifacts/onlyempty_carrier_cells_primary.toml"))
    @test primary["gmsh_version"]=="4.15.2"
    @test primary["primary_json_sha256"]==
        "9bc22cd39d526b6af80c4722f1b8896475c66624bb2aa06adc73d03233934dfb"
    for row in primary["rows"]
        mktempdir() do directory
            path=joinpath(directory,"cells.geo")
            tail=row["mode"]=="quad" ? "Recombine Surface{1}; Mesh 2;\n" :
                "Mesh 2; RefineMesh;\n"
            write(path,square_source(row["order"],row["renumber"])*tail)
            result=GEO.execute_geo(path)
            check_literal_carrier(result,row)
        end
    end
end

@testset "Cell traces follow duplicate and owner bucket merge events" begin
    line=Mesh([0. 1.;0. 0.;0. 0.];segs=reshape(Int32[1,2],2,1))
    identities=[[UInt64[0],UInt64[],UInt64[]],[UInt64[201],UInt64[],UInt64[]]]
    _,_,keeps,cells=API._merge_entity_meshes([(1,line),(1,line)],"cell provenance";
        node_identities=[Int32[1,2],Int32[1,2]],cell_identities=identities)
    @test keeps==[BitVector([true]),BitVector([false])]
    @test cells==[UInt64[201],UInt64[],UInt64[]]
    conflict=[[UInt64[201],UInt64[],UInt64[]],[UInt64[202],UInt64[],UInt64[]]]
    @test_throws ArgumentError API._merge_entity_meshes([(1,line),(1,line)],"cell provenance";
        node_identities=[Int32[1,2],Int32[1,2]],cell_identities=conflict)
    mixed=MixedMesh(line.coords,[ElementBlock(1,line.segs,Int32[0])])
    _,_,keeps,cells=API._merge_entity_meshes([(1,mixed),(1,line)],"cell provenance";
        node_identities=[Int32[1,2],Int32[1,2]],
        cell_identities=[[UInt64[0]],[UInt64[201],UInt64[],UInt64[]]])
    @test cells==[UInt64[201]]
    @test keeps[2][1]==BitVector([false])
    merged,_,_,cells=GEO._geo_merge_entity_meshes([(1,1,mixed),(1,2,line)];
        node_identities=[Int32[1,2],Int32[1,2]],
        cell_identities=[[UInt64[201]],[UInt64[202],UInt64[],UInt64[]]],return_remaps=true)
    @test length(merged.blocks)==2
    @test cells==[UInt64[201],UInt64[202]]
end

@testset "Compatibility and malformed cell binding admission" begin
    mktempdir() do directory
        path=joinpath(directory,"cells.geo")
        write(path,square_source(1,0)*"Mesh 2;\n")
        result=GEO.execute_geo(path);a=result.model.meshing.mesh_authority
        legacy=Model._ModelMeshAuthority(a.mesh,a.products,a.done,a.private_rank,
            a.native_inputs,a.carrier,a.carrier_rows,a.producer_rows)
        @test isempty(legacy.carrier_cells)
        @test_throws ArgumentError Model._model_mesh_authority(a.mesh,a.done;
            carrier=a.carrier,carrier_rows=a.carrier_rows,carrier_cells=[UInt64[1]])
        malformed=deepcopy(a.carrier_cells);malformed[2][1]=UInt64(1) # A Point is not a Surface.
        @test_throws ArgumentError Model._model_mesh_authority(a.mesh,a.done;
            carrier=a.carrier,carrier_rows=a.carrier_rows,carrier_cells=malformed)
        @test a.carrier_cells[2]==UInt64[9,10]
    end
end
end
