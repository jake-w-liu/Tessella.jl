using Test, Tessella
using Tessella.Elements: ElementBlock, MixedMesh, MixedEntityData, msh_dimension
using Tessella.MeshTypes: nnodes
include(joinpath(@__DIR__,"quadtri_nonew_certificates.jl"))
include(joinpath(@__DIR__,"quadtri_nonew_b4_free_strip_certificates.jl"))
include(joinpath(@__DIR__,"..","..","validation","native_geometry","native_occ_cube_generation_primary_v1.jl"))

const _MAP_Q=QuadTriNoNewCertificates
const _MAP_F=QuadTriNoNewB4FreeStripCertificates
const _MAP_M=Tessella.Model

function _map_cells(mesh;dimension=nothing)
    data=mesh.entity_data
    Dict(data.external_element_tags[bi][column] =>
        (block.msh,data.block_entities[bi][column],block.tags[column],
         Tuple(data.external_node_tags[block.nodes[:,column]]))
        for (bi,block) in enumerate(mesh.blocks) for column in axes(block.nodes,2)
        if dimension===nothing || msh_dimension(block.msh)==dimension)
end

function _map_reorder(mesh,permutation)
    inverse=invperm(permutation);data=mesh.entity_data
    blocks=[ElementBlock(block.msh,Int32.(inverse[block.nodes]),block.tags)
            for block in mesh.blocks]
    metadata=MixedEntityData(data.entities;
        node_entities=data.node_entities[permutation],
        node_parametric=data.node_parametric[permutation],
        external_node_tags=data.external_node_tags[permutation],
        block_entities=data.block_entities,external_element_tags=data.external_element_tags)
    MixedMesh(mesh.coords[:,permutation],blocks;physical_names=mesh.physical_names,
        entity_data=metadata,elementary_entities=mesh.elementary_entities,
        periodic_links=mesh.periodic_links)
end

@testset "actual projection boundaries" begin
@testset "current actual Volume product projection" begin
    for motion in (:rotation,:twist)
        fixture=_MAP_Q.fixture("actual_volume_$(motion)","Layers{3}",
            [0.,1/3,2/3,1.],false;motion)
        source=fixture.source*"""
        Physical Surface("source",101)={1};
        Physical Surface("top",102)={sweep[0]};
        Physical Volume("body",103)={sweep[1]};
        """
        execution=_MAP_Q.execute(source)
        volume=geo_entity_mesh(execution,3,1)
        certificate=_MAP_Q.certify_fixture(volume,fixture)
        @test validate(volume).ok
        @test volume.entity_data!==nothing
        projected=model_to_mixed(execution.model,volume,3,1)
        @test validate(projected).ok
        @test projected.physical_names==Dict((2,101)=>"source",(2,102)=>"top",(3,103)=>"body")
        @test _map_cells(projected;dimension=3)==_map_cells(volume;dimension=3)
        @test _MAP_Q.surface_faces([(2,0,projected)],volume)==Set(keys(certificate.boundary))
        reordered=_map_reorder(volume,collect(nnodes(volume):-1:1))
        again=model_to_mixed(execution.model,reordered,3,1)
        @test again.entity_data.external_node_tags[1:nnodes(reordered)]==reordered.entity_data.external_node_tags
        @test _map_cells(again)==_map_cells(projected)
        @test again.coords[:,1:nnodes(reordered)]==reordered.coords
    end
end

@testset "NoNew projection consumes actual retained Source and cap" begin
    for witness in _MAP_F.RETAINED_WITNESSES
        retained=_MAP_F.retained_product(witness)
        model=retained.geometry.model
        model.meshing.mesh_authority=_MAP_M._model_mesh_authority(retained.projected)
        region=Int(retained.geometry.lists["sweep"][2])
        @test _MAP_F.certify(retained.volume,retained.f,retained.source).ncell==retained.certificate.ncell
        # This naked geometry path still has the original strict projection
        # contract. Its completed plan must consume the model's actual products.
        projected=model_to_mixed(model,retained.volume,3,region)
        @test validate(projected).ok
        @test _MAP_Q.surface_faces([(2,0,projected)],retained.volume)==
              Set(keys(retained.certificate.boundary))
        source=_MAP_M._model_mesh_product(model.meshing.mesh_authority,2,retained.f.surface)
        source_projected=model_to_mixed(model,source,retained.f.surface)
        @test _map_cells(source_projected;dimension=2)==_map_cells(source;dimension=2)
    end
end
end

function _map_repack(mesh;coords=mesh.coords,blocks=mesh.blocks,
        entities=mesh.entity_data.entities,node_entities=mesh.entity_data.node_entities,
        parameters=mesh.entity_data.node_parametric,node_tags=mesh.entity_data.external_node_tags,
        owners=mesh.entity_data.block_entities,cell_tags=mesh.entity_data.external_element_tags)
    data=MixedEntityData(entities;node_entities=node_entities,node_parametric=parameters,
        external_node_tags=node_tags,block_entities=owners,external_element_tags=cell_tags)
    MixedMesh(coords,blocks;physical_names=mesh.physical_names,entity_data=data,
        elementary_entities=owners,periodic_links=mesh.periodic_links)
end

@testset "malformed actual projection is atomic" begin
    fixture=_MAP_Q.fixture("actual_projection_guards","Layers{3}",[0.,1/3,2/3,1.],false)
    execution=_MAP_Q.execute(fixture.source*"Physical Volume(103)={sweep[1]};")
    model=execution.model;volume=geo_entity_mesh(execution,3,1)
    authority=model.meshing.mesh_authority
    original=deepcopy(model.curve_params)
    node_tags=copy(volume.entity_data.external_node_tags);node_tags[1]=maximum(node_tags)+1000
    element_tags=deepcopy(volume.entity_data.external_element_tags)
    element_tags[1][1]=maximum(vcat(element_tags...))+1000
    coordinates=copy(volume.coords);coordinates[1,1]+=.01
    parameters=deepcopy(volume.entity_data.node_parametric)
    row=findfirst(owner->owner[1]==1,volume.entity_data.node_entities)
    parameters[row]=Float64[.123456]
    blocks=copy(volume.blocks)
    slot=findfirst(block->size(block.nodes,2)>=2,blocks)
    nodes=copy(blocks[slot].nodes);nodes[:,[1,2]]=nodes[:,[2,1]]
    blocks[slot]=ElementBlock(blocks[slot].msh,nodes,blocks[slot].tags)
    malformed=(_map_repack(volume;node_tags),_map_repack(volume;cell_tags=element_tags),
        _map_repack(volume;coords=coordinates),_map_repack(volume;parameters),
        _map_repack(volume;blocks))
    for input in malformed
        @test validate(input).ok # These are well-formed meshes, with forged current identities.
        @test_throws ArgumentError model_to_mixed(model,input,3,1)
        @test model.curve_params==original
        @test model.meshing.mesh_authority===authority
    end
    # A naked tagged mesh cannot use current authority to bypass the original
    # unclassified input contract.
    naked=MixedMesh(volume.coords,volume.blocks)
    @test_throws r"tags must be zero" model_to_mixed(model,naked,3,1)
    @test model.curve_params==original
    @test model.meshing.mesh_authority===authority
end

@testset "actual quadratic Source, cap and Volume identities" begin
    fixture=_MAP_F.fixture("actual_quadratic_projection";strip_length=5,layers=:one)
    geometry=_MAP_F.execute(fixture.source;dim=2)
    linear_source=_MAP_M._model_mesh_product(geometry.model.meshing.mesh_authority,
        2,fixture.surface;producer_order=true)
    execution=_MAP_F.execute(fixture.source*"Mesh.ElementOrder=2;")
    region=Int(execution.lists["sweep"][2])
    volume=geo_entity_mesh(execution,3,region)
    certificate=_MAP_F.certify_quadratic(volume,fixture,linear_source)
    @test certificate.nsupport==nnodes(volume)
    authority=execution.model.meshing.mesh_authority
    projected=model_to_mixed(execution.model,volume,3,region)
    @test validate(projected).ok
    @test _map_cells(projected;dimension=3)==_map_cells(volume;dimension=3)
    @test projected.entity_data.external_node_tags[1:nnodes(volume)]==volume.entity_data.external_node_tags
    for surface in (fixture.surface,Int(execution.lists["sweep"][1]))
        actual=_MAP_M._model_mesh_product(authority,2,surface)
        owned=model_to_mixed(execution.model,actual,surface)
        @test validate(owned).ok
        @test _map_cells(owned;dimension=2)==_map_cells(actual;dimension=2)
        @test owned.entity_data.node_parametric[1:nnodes(actual)]==actual.entity_data.node_parametric
        @test owned.entity_data.external_node_tags[1:nnodes(actual)]==actual.entity_data.external_node_tags
    end
end


@testset "ordinary analytic Cube actual products" begin
    source="""
    SetFactory("OpenCASCADE");
    Box(1)={0,0,0,1,1,1};
    Transfinite Curve{:}=2;
    Transfinite Surface{:}; Recombine Surface{:};
    Transfinite Volume{:};
    Physical Volume("solid",103)={1};
    """
    for order in (1,2)
        execution=_MAP_Q.execute(source*"Mesh.ElementOrder=$order;")
        volume=geo_entity_mesh(execution,3,1)
        @test length(volume.blocks)==1
        @test only(volume.blocks).msh==(order==1 ? 5 : 12)
        @test size(only(volume.blocks).nodes,2)==1
        native=only(row for row in NATIVE_OCC_CUBE_GENERATION_PRIMARY_V1 if row.order==order)
        authority=execution.model.meshing.mesh_authority
        actual=authority.mesh;data=actual.entity_data
        positions=Dict(uid=>row for (row,uid) in enumerate(data.external_node_tags))
        parameters=_MAP_M._model_mesh_runtime_parameters(execution.model,actual)
        @test Set(data.external_node_tags)==Set(node.uid for node in native.nodes)
        @test nnodes(volume)==length(native.nodes)
        for node in native.nodes
            row=positions[node.uid]
            @test data.node_entities[row]==node.owner
            @test reinterpret(UInt64,actual.coords[:,row])==reinterpret(UInt64,collect(node.xyz))
            stored=parameters[row]===nothing ? Float64[] : parameters[row]
            @test reinterpret(UInt64,stored)==reinterpret(UInt64,node.parameters)
        end
        native_cells=Dict(tag=>(cell.msh,cell.owner[2],
                Int32(cell.owner[1]==3 ? 103 : 0),Tuple(nodes[:,column]))
            for cell in native.cells
            for nodes in (reshape(cell.nodes,Tessella.Elements.msh_num_nodes(cell.msh),:),)
            for (column,tag) in enumerate(cell.tags))
        @test _map_cells(actual)==native_cells
        projected=model_to_mixed(execution.model,volume,3,1)
        @test validate(projected).ok
        @test projected.physical_names==Dict((3,103)=>"solid")
        @test _map_cells(projected;dimension=3)==_map_cells(volume;dimension=3)
        @test Dict(dimension=>sum(size(block.nodes,2) for block in projected.blocks
            if msh_dimension(block.msh)==dimension;init=0) for dimension in 0:3)==
              Dict(0=>8,1=>12,2=>6,3=>1)
    end
end

