# The suite includes every test file into Main; this file keeps its aliases to itself.
module EmittedReceiptStampTests
using Test
using Tessella
import Tessella.Model as M
const A = Tessella.API

# Every generation publishes a receipt of its emitted parts. The receipt
# authenticates the actual mesh by a fingerprint of its native state (not a
# deep copy of it) and keeps the owner of every carrier cell in tables with
# concrete keys. Both must tell apart everything the copies told apart.

@testset "Carrier cell owners are kept per arity and found again" begin
    owners = M._EmittedCellOwners()
    rows = Int32[11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    for (msh, arity) in ((1, 2), (2, 3), (4, 4), (7, 5), (6, 6), (5, 8), (11, 10))
        cells = reshape(Int32.(collect(1:arity)), arity, 1)
        carrier = reshape(rows[1:arity], arity, 1)
        @test M._emitted_owner(owners, msh, carrier, 1) === nothing
        @test M._emitted_owner_insert!(owners, msh, cells, 1, rows, (3, Int32(msh)))
        @test !M._emitted_owner_insert!(owners, msh, cells, 1, rows, (3, Int32(99)))   # already owned
        @test M._emitted_owner(owners, msh, carrier, 1) == (3, Int32(msh))
        # The element type is part of the key; so is the order of the rows.
        @test M._emitted_owner(owners, msh + 100, carrier, 1) === nothing
        @test M._emitted_owner(owners, msh, reverse(carrier; dims=1), 1) === nothing
    end
    @test length(owners) == 7
    @test (4, (Int32(11), Int32(12), Int32(13), Int32(14))) in M._emitted_owner_keys(owners)
    @test length(M._emitted_owner_keys(owners)) == 7
    # Rows given as Int index the same cell.
    @test !M._emitted_owner_insert!(owners, 4, reshape([1, 2, 3, 4], 4, 1), 1, collect(11:14), (3, Int32(4)))
end

@testset "The state fingerprint follows every word of the native state" begin
    A.initialize()
    try
        A.option("General.Terminal", 0)
        A.model.occ.add_box(0, 0, 0, 1, 1, 1); A.model.occ.synchronize()
        A.option("Mesh.MeshSizeMin", 0.5); A.option("Mesh.MeshSizeMax", 0.5)
        A.mesh.generate(3)
        authority = A.CURRENT[].meshing.mesh_authority
        receipt = authority.emitted_geometry
        @test receipt !== nothing
        mesh = authority.mesh; data = mesh.entity_data
        stamp = receipt.native_state
        @test stamp isa M._EmittedStateStamp
        @test stamp.nodes == size(mesh.coords, 2)
        @test stamp.cells == sum(block -> size(block.nodes, 2), mesh.blocks)
        @test M._model_mesh_emitted_state_matches(mesh, stamp)
        @test M._model_mesh_emitted_state(deepcopy(mesh)) == stamp
        # The owner of every carrier cell is on record.
        blocks = collect(M._emitted_projection_blocks(authority.carrier))
        @test length(receipt.cell_owners) == length(blocks)
        @test all(length(owners) == size(cells, 2) && all(!=(M._EMITTED_NO_OWNER), owners)
                  for (owners, (msh, cells, tags)) in zip(receipt.cell_owners, blocks))

        # One edit at a time, undone after it is seen.
        function seen(edit!, undo!)
            edit!()
            changed = !M._model_mesh_emitted_state_matches(mesh, stamp)
            undo!()
            return changed && M._model_mesh_emitted_state_matches(mesh, stamp)
        end
        x = mesh.coords[2, 3]
        @test seen(() -> (mesh.coords[2, 3] = nextfloat(x)), () -> (mesh.coords[2, 3] = x))
        @test seen(() -> (mesh.coords[2, 3] = x + 1.0), () -> (mesh.coords[2, 3] = x))
        tag = data.external_node_tags[end]
        @test seen(() -> (data.external_node_tags[end] = tag + 1), () -> (data.external_node_tags[end] = tag))
        owner = data.node_entities[1]
        @test seen(() -> (data.node_entities[1] = (owner[1], owner[2] + Int32(1))), () -> (data.node_entities[1] = owner))
        parameter = findfirst(p -> p !== nothing && !isempty(p), data.node_parametric)
        if parameter !== nothing
            saved = data.node_parametric[parameter]
            @test seen(() -> (data.node_parametric[parameter] = nothing), () -> (data.node_parametric[parameter] = saved))
            value = saved[1]
            @test seen(() -> (saved[1] = prevfloat(value)), () -> (saved[1] = value))
        end
        block = mesh.blocks[end]
        a, b = block.nodes[1, 1], block.nodes[2, 1]
        @test seen(() -> (block.nodes[1, 1] = b; block.nodes[2, 1] = a), () -> (block.nodes[1, 1] = a; block.nodes[2, 1] = b))
        physical = block.tags[1]
        @test seen(() -> (block.tags[1] = physical + Int32(1)), () -> (block.tags[1] = physical))
        entity_of = data.block_entities[end][1]
        @test seen(() -> (data.block_entities[end][1] = entity_of + Int32(1)), () -> (data.block_entities[end][1] = entity_of))
        label = data.external_element_tags[1][1]
        @test seen(() -> (data.external_element_tags[1][1] = label + 1), () -> (data.external_element_tags[1][1] = label))
        key = first(sort!(collect(keys(data.entities))))
        entity = data.entities[key]
        @test seen(() -> push!(entity.physical_tags, Int32(7)), () -> pop!(entity.physical_tags))
        @test seen(() -> delete!(data.entities, key), () -> (data.entities[key] = entity))

        # The order in which the entity table was filled is not part of the state.
        refilled = typeof(data.entities)()
        for entry in reverse(sort!(collect(keys(data.entities))))
            refilled[entry] = data.entities[entry]
        end
        @test M._emitted_mix(M._EMITTED_MIX_SEED, refilled) == M._emitted_mix(M._EMITTED_MIX_SEED, data.entities)
    finally
        A.finalize()
    end
end

@testset "Values are mixed as isequal compares them" begin
    seed = M._EMITTED_MIX_SEED
    mix(x) = M._emitted_mix(seed, x)
    @test mix(Int32(5)) == mix(5) == mix(UInt64(5))
    @test mix(Int32(-1)) == mix(-1)
    @test mix(0.0) != mix(-0.0)
    @test mix(nothing) != mix(0)
    @test mix([1, 2]) != mix([2, 1])
    @test mix(Int32[1, 2]) == mix([1, 2])
    @test mix((1, 2.5)) != mix((1, 2.5, nothing))
    @test mix(Float64[]) != mix([0.0])
    @test mix(Set([3, 1, 2])) == mix(Set([2, 3, 1]))
    @test mix(Dict(1 => [1.0], 2 => [2.0])) != mix(Dict(1 => [2.0], 2 => [1.0]))
    @test mix("ab") != mix("ba")
    @test mix(:name) == mix("name")
end
end
