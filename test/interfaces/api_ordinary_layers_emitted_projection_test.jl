module OrdinaryLayersEmittedProjectionChecks
using Tessella,Test
using Tessella.MeshTypes: nnodes
const M=Tessella.Model
const SOURCE=raw"""
Point(1)={0,0,0,0.5};Point(2)={1,0,0,0.5};Point(3)={0,1,0,0.5};
Line(1)={1,2};Line(2)={2,3};Line(3)={3,1};
Curve Loop(1)={1,2,3};Plane Surface(1)={1};
Extrude{0,0,2}{Surface{1};Layers{3};Recombine;}
"""
counts(mesh)=Dict(kind=>sum(size(block.nodes,2) for block in mesh.blocks if block.msh==kind)
    for kind in unique(block.msh for block in mesh.blocks))
function execute()
    mktemp() do path,io
        write(io,SOURCE);close(io);execute_geo(path;mesh_dim=3)
    end
end
function signed_volume(a,b,c,d)
    q=Rational{BigInt}
    u=ntuple(i->q(b[i])-q(a[i]),3);v=ntuple(i->q(c[i])-q(a[i]),3)
    w=ntuple(i->q(d[i])-q(a[i]),3)
    sum(u[i]*(v[mod1(i+1,3)]*w[mod1(i+2,3)]-v[mod1(i+2,3)]*w[mod1(i+1,3)])
        for i in 1:3)/6
end
function prism_volumes(mesh)
    [sum(signed_volume((Tuple(mesh.coords[:,cell[i]]) for i in slots)...)
         for slots in ((1,2,3,6),(1,2,6,5),(1,5,6,4)))
     for block in mesh.blocks if block.msh==6 for cell in eachcol(block.nodes)]
end
function receipt_parts(receipt)
    owners=sort!(collect(keys(receipt.parts)))
    parts=[(owner[1],Int(owner[2]),receipt.parts[owner][1]) for owner in owners]
    rows=[copy(receipt.parts[owner][2]) for owner in owners]
    parts,rows
end
@testset "emitted ordinary prism projection retains its owned constructor receipt" begin
    execution=execute()
    @test counts(execution.mesh)==Dict(1=>23,2=>14,3=>21,6=>21)
    @test nnodes(execution.mesh)==32
    volume=geo_entity_mesh(execution,3,1;geometry=true)
    native=geo_entity_mesh(execution,3,1)
    authority=execution.model.meshing.mesh_authority
    before=deepcopy(authority.mesh)
    @test authority.emitted_geometry!==nothing
    projected=model_to_mixed(execution.model,volume,3,1)
    @test validate(projected).ok
    @test counts(projected)[6]==21 && counts(projected)[3]==21
    @test projected.coords[:,1:nnodes(volume)]==volume.coords
    @test all(isnothing,projected.entity_data.node_parametric)
    @test !projected.entity_data.native_mesh_tags
    volumes=prism_volumes(volume)
    @test length(volumes)==21 && all(>(0),volumes)
    @test Float64(sum(volumes))≈1.0 atol=2e-11
    @test authority.mesh.coords==before.coords
    @test authority.mesh.entity_data.node_parametric==before.entity_data.node_parametric
    # Receipt publication checks ownership and row identity before replacing
    # the authority. Missing, foreign and repeated constructor entries fail.
    parts,rows=receipt_parts(authority.emitted_geometry)
    duplicate_parts=vcat(parts,parts[1:1]);duplicate_rows=vcat(rows,rows[1:1])
    @test_throws ArgumentError M._model_mesh_with_emitted_receipt(
        authority,duplicate_parts,duplicate_rows;caller="duplicate emitted owner")
    foreign=copy(parts);foreign[1]=(3,999,foreign[1][3])
    @test_throws ArgumentError M._model_mesh_with_emitted_receipt(
        authority,foreign,rows;caller="foreign emitted owner")
    missing=deepcopy(rows);pop!(missing[1])
    @test_throws ArgumentError M._model_mesh_with_emitted_receipt(
        authority,parts,missing;caller="missing emitted row")
    @test execution.model.meshing.mesh_authority===authority
    @test M._model_mesh_emitted_state_matches(authority.mesh,authority.emitted_geometry.native_state)
    # Typed actual event receipts must remain literal and current.
    event_model=deepcopy(execution.model)
    event_authority=event_model.meshing.mesh_authority
    block=findfirst(block->block.msh==6,event_authority.carrier.blocks)
    event_authority.carrier_cells[block][1]=nothing
    @test_throws ArgumentError model_to_mixed(event_model,volume,3,1)
    @test M._model_mesh_emitted_state_matches(event_authority.mesh,event_authority.emitted_geometry.native_state)
    state_model=deepcopy(execution.model)
    state_authority=state_model.meshing.mesh_authority
    state_authority.mesh.coords[1,1]=nextfloat(state_authority.mesh.coords[1,1])
    @test_throws ArgumentError model_to_mixed(state_model,volume,3,1)
    point_model=deepcopy(execution.model)
    M.set_point_coordinates!(point_model,1,0.125,0.0,0.0)
    @test_throws ArgumentError model_to_mixed(point_model,volume,3,1)
    @test point_model.points[1]==(0.125,0.0,0.0)
    @test authority.mesh.entity_data.external_node_tags==before.entity_data.external_node_tags
    @test authority.mesh.entity_data.external_element_tags==before.entity_data.external_element_tags
    @test all(left.nodes==right.nodes && left.tags==right.tags for
        (left,right) in zip(authority.mesh.blocks,before.blocks))
    # A one-ULP edit keeps constructor cells/row order, then fails exact XYZ
    # verification after the owned receipt was selected.
    changed=deepcopy(volume)
    changed.coords[1,1]=nextfloat(changed.coords[1,1])
    @test_throws ArgumentError model_to_mixed(execution.model,changed,3,1)
    @test authority.mesh.coords==before.coords
    @test authority.mesh.entity_data.node_parametric==before.entity_data.node_parametric
    # Rebuilding the authority cannot preserve an old emitted admission.
    updated=M._model_mesh_authority(authority.mesh,authority.done;
        private_rank=authority.private_rank,carrier=authority.carrier,
        carrier_rows=authority.carrier_rows,carrier_cells=authority.carrier_cells)
    @test updated.emitted_geometry===nothing
    @test M._model_mesh_emitted_projection(execution.model,volume,3,999,
        authority.emitted_geometry,"unmatched emitted owner")===nothing
    @test M._model_mesh_emitted_projection(execution.model,native,3,1,
        authority.emitted_geometry,"tagged actual input")===nothing
    @test_throws ArgumentError M._model_mesh_projection_input(
        execution.model,projected,3,1,"geometry role test")
    actual_projected=model_to_mixed(execution.model,native,3,1)
    @test validate(actual_projected;reject_duplicate_point_cells=false).ok
    @test actual_projected.entity_data.node_parametric[1:nnodes(native)]==native.entity_data.node_parametric
    @test actual_projected.entity_data.external_node_tags[1:nnodes(native)]==native.entity_data.external_node_tags
end
end
