module NativeHybridActualModelChecks
using Test,TOML,SHA,Tessella
using Tessella.MeshTypes: Mesh,nnodes,validate
using Tessella.Elements: ElementBlock,MixedMesh,MixedEntityData,MixedEntity
const M=Tessella.Model
const artifacts=joinpath(@__DIR__,"..","artifacts")
const raw=TOML.parsefile(joinpath(artifacts,"native_hybrid_current_lower_4152_v2.toml"))
hexpoint(xs)=Float64[Core.eval(@__MODULE__,Meta.parse(x)) for x in xs]
function box_model()
    m=M.GeoModel()
    points=((0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1))
    for (tag,p) in enumerate(points);M.add_point!(m,p...;tag,mesh_size=0.4);end
    curves=((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),(7,8),(8,5),(1,5),(2,6),(3,7),(4,8))
    for (tag,ends) in enumerate(curves);M.add_line!(m,ends...;tag);end
    loops=((1,2,3,4),(5,6,7,8),(1,10,-5,-9),(2,11,-6,-10),(3,12,-7,-11),(4,9,-8,-12))
    for (tag,curves) in enumerate(loops)
        M.add_curve_loop!(m,collect(curves);tag)
        M.add_plane_surface!(m,[tag];tag)
    end
    M.add_surface_loop!(m,collect(1:6);tag=1)
    M.add_volume!(m,[1];tag=1)
    return m
end
function actual_lower(m,case)
    products=filter(p->p["owner"][1]<3,case["product"])
    records=Dict{UInt64,NamedTuple}()
    for product in products
        dim,tag=product["owner"];ids=product["node_tags"]
        for (ordinal,uid) in enumerate(ids)
            @assert !haskey(records,UInt64(uid))
            uv=dim==0 ? nothing : hexpoint(product["uv_hex"][dim*(ordinal-1)+1:dim*ordinal])
            records[UInt64(uid)]=(owner=(Int(dim),Int32(tag)),xyz=hexpoint(product["xyz_hex"][ordinal]),uv=uv)
        end
    end
    ids=sort!(collect(keys(records)));positions=Dict(uid=>Int32(row) for (row,uid) in enumerate(ids))
    coords=reduce(hcat,records[uid].xyz for uid in ids)
    entities=Dict{Tuple{Int,Int},MixedEntity}()
    blocks=ElementBlock[];owners=Vector{Int32}[];labels=Vector{UInt64}[]
    for product in products
        dim,tag=Int.(product["owner"])
        supports=UInt64[uid for cell in get(product,"cell",[]) for nodes in cell["nodes"] for uid in nodes]
        append!(supports,UInt64.(product["node_tags"]));unique!(supports)
        bounds=isempty(supports) ? (0.,0.,0.,0.,0.,0.) :
            M._model_projection_bbox(Mesh(coords),Int32[positions[uid] for uid in supports],"native lower")
        boundary=dim==0 ? Int32[] : dim==1 ? Int32[m.curves[tag]...] :
            M._model_projection_surface_boundaries(m,tag)
        entities[(dim,tag)]=M._model_mesh_entity_metadata(m,dim,tag,bounds;boundaries=boundary)
        for cell in get(product,"cell",[])
            nodes=hcat((Int32[positions[UInt64(uid)] for uid in vertices] for vertices in cell["nodes"])...)
            push!(blocks,ElementBlock(cell["msh"],nodes))
            push!(owners,fill(Int32(tag),length(cell["labels"])))
            push!(labels,UInt64.(cell["labels"]))
        end
    end
    data=MixedEntityData(entities;node_entities=[records[uid].owner for uid in ids],
        node_parametric=[records[uid].uv for uid in ids],external_node_tags=ids,
        block_entities=owners,external_element_tags=labels)
    lower=MixedMesh(coords,blocks;entity_data=data)
    authority=M._model_mesh_authority(lower,keys(entities);caller="native owned lower")
    constructors=Dict((2,Int32(tag))=>M._model_mesh_product(authority,2,tag;caller="native owned lower") for tag in 1:6)
    return M._model_mesh_with_surface_constructor_inputs(authority,constructors)
end
@testset "actual Model native hybrid Volume production" begin
    for case in raw["case"]
        println("HYBRID_ACTUAL_MODEL_CASE ",case["name"]);flush(stdout)
        m=box_model();input=actual_lower(m,case);saved=deepcopy(input.mesh)
        constructor=M._model_native_hybrid_volume_input(m,1,input,"scoped hybrid")
        @test M._model_native_hybrid_volume_scope(m,1,constructor,"scoped hybrid")===nothing
        multiple=deepcopy(m);multiple.volumes[2]=copy(multiple.volumes[1])
        @test_throws ArgumentError M._model_native_hybrid_volume_scope(multiple,1,constructor,"multi Region")
        foreign=deepcopy(input)
        source=foreign.products[(2,Int32(1))];quad=first(source.blocks)
        foreign.products[(2,Int32(99))]=M._ModelMeshProduct((2,Int32(99)),copy(source.rows),
            [ElementBlock(2,quad.nodes[1:3,1:1])],[UInt64[999999]])
        @test_throws ArgumentError M._model_native_hybrid_volume_scope(m,1,merge(constructor,(input=foreign,)),"foreign sizing Tri")
        options=M._model_default_mesh1d_options(m,"native owned lower")
        # The bound native lower RNG continuation records no draws for these
        # four cases. The production consumes this state through its options.
        options.random_state[]=UInt32(1)
        options=endswith(case["name"],"no_tet_optimizer") ?
            M._ModelMesh1DOptions((name==:optimize ? 0 : getfield(options,name)
                for name in fieldnames(M._ModelMesh1DOptions))...) : options
        production=M.mesh_model_volume(m,1;_mesh_input=input,_mesh1d_options=options,
            _return_volume_production=true,_mesh_node_high_water=UInt64(case["node_clock"]),
            _mesh_element_high_water=UInt64(case["current_lower_element_high_water"]))
        @test production isa M._ModelNativeHybridVolumeProduction
        @test validate(production.mesh).ok
        @test production.input===input
        @test production.node_high_water==UInt64(case["node_high_water"])
        @test production.element_high_water==UInt64(case["element_high_water"])
        expected=only(filter(p->p["owner"]==[3,1],case["product"]))
        wanted_ids=UInt64[uid for cell in expected["cell"] for uid in cell["labels"]]
        wanted_nodes=[Tuple(UInt64.(nodes)) for cell in expected["cell"] for nodes in cell["nodes"]]
        @test production.element_tags==wanted_ids
        @test [Tuple(production.boundary_node_tags[nodes]) for block in production.mesh.blocks for nodes in eachcol(block.nodes)]==wanted_nodes
        actual_rows=Dict(uid=>row for (row,uid) in enumerate(production.boundary_node_tags))
        wanted_coords=reduce(hcat,hexpoint.(expected["xyz_hex"]))
        @test isequal(production.mesh.coords[:,[actual_rows[UInt64(uid)] for uid in expected["node_tags"]]],wanted_coords)
        @test isequal(input.mesh.coords,saved.coords)
        @test input.mesh.entity_data.external_node_tags==saved.entity_data.external_node_tags
        @test [block.nodes for block in input.mesh.blocks]==[block.nodes for block in saved.blocks]
        @test isequal(input.mesh.entity_data.node_parametric,saved.entity_data.node_parametric)
    end
end
end
