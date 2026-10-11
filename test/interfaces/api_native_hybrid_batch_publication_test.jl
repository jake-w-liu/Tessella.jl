module NativeHybridBatchPublicationChecks
using Test,TOML,SHA,Tessella
using Tessella.MeshTypes: Mesh,nnodes,validate
using Tessella.Elements: ElementBlock,MixedMesh,MixedEntityData,MixedEntity,msh_dimension
const M=Tessella.Model
const A=Tessella.API
const artifacts=joinpath(@__DIR__,"..","artifacts")
const native=TOML.parsefile(joinpath(artifacts,"native_hybrid_batch_scheduling_4152_v1.toml"))
function helper(expr,name)
    expr isa Expr || return nothing
    if expr.head===:function && expr.args[1] isa Expr && expr.args[1].head===:call && expr.args[1].args[1]===name
        return expr
    end
    for child in expr.args
        found=helper(child,name);found===nothing || return found
    end
    return nothing
end
for (file,names) in (("api_native_hybrid_actual_model_test.jl",(:box_model,:actual_lower)),
        ("api_native_hybrid_batch_scheduling_test.jl",(:disjoint_model,)),
        ("api_native_hybrid_connected_components_test.jl",(:adjacent_model,)))
    path=joinpath(@__DIR__,file)
    body=replace(read(path,String),"for tag in 1:6"=>"for tag in sort!(collect(keys(m.surfaces)))")
    for name in names
        selected=helper(Meta.parseall(body;filename=path),name)
        @assert selected!==nothing;Core.eval(@__MODULE__,selected)
    end
end
hexpoint(xs)=Float64[Core.eval(@__MODULE__,Meta.parse(x)) for x in xs]
function assembly(input,case)
    actual=input.mesh;data=actual.entity_data
    nodes=A._Dim01Node[A._Dim01Node(uid,ntuple(axis->actual.coords[axis,row],3),data.node_entities[row],
        data.node_parametric[row]===nothing ? nothing : copy(data.node_parametric[row]))
        for (row,uid) in enumerate(data.external_node_tags)]
    cells=A._Dim01Cell[A._Dim01Cell(Int(block.msh),(msh_dimension(block.msh),owners[column]),
        copy(block.nodes[:,column]),uid,uid,block.tags[column])
        for (block,labels,owners) in zip(actual.blocks,data.external_element_tags,data.block_entities)
        for (column,uid) in enumerate(labels)]
    return A._Dim01Assembly(nodes,cells,Dict(node.tag=>Int32(row) for (row,node) in enumerate(nodes)),
        UInt64(case["node_clock"]),UInt64(case["current_lower_element_high_water"]),
        1_000_000,1_000_000,"native batch publication",false,false)
end
function classification(production,tag,input)
    data=input.mesh.entity_data;rows=Dict(uid=>row for (row,uid) in enumerate(data.external_node_tags))
    owners=Tuple{Int,Int32}[];parameters=Union{Nothing,Vector{Float64}}[]
    for uid in production.boundary_node_tags
        row=get(rows,uid,0)
        push!(owners,row==0 ? (3,Int32(tag)) : data.node_entities[row])
        parameter=row==0 ? nothing : data.node_parametric[row]
        push!(parameters,parameter===nothing ? nothing : copy(parameter))
    end
    # This controlled classification seam carries literal native lower owners
    # and raw parameters. The tested append/allocator/receipt bodies are actual
    # Source; public lifecycle classification is qualified separately.
    return (mesh=production.mesh,entity=(3,Int32(tag)),node_entities=owners,node_parametric=parameters,
        cell_entities=Dict(Int(block.msh)=>fill(Int32(tag),size(block.nodes,2)) for block in production.mesh.blocks))
end
function state(part)
    return (nodes=Tuple((node.tag,node.point,node.owner,node.parameter===nothing ? nothing : Tuple(node.parameter)) for node in part.nodes),
        cells=Tuple((cell.tag,cell.msh,cell.owner,Tuple(cell.nodes),cell.physical,cell.source_tag) for cell in part.cells),
        positions=copy(part.node_positions),clocks=(part.node_max,part.cell_max))
end
function publish!(part,m,production,class,tag,receipt)
    return A._higher_append_part!(part,m,production.mesh,class,3,tag,Set{Tuple{Int,Int32}}();
        source_row_tags=production.boundary_node_tags,source_node_base=receipt.node_base,
        source_node_high_water=receipt.batch.node_high_water,source_cell_tags=production.element_tags,
        source_cell_high_water=receipt.batch.element_high_water,source_batch_receipt=receipt,
        source_row_order=production.registration_rows,source_volume_production=production)
end
@testset "operation-owned native hybrid publication preserves global births and strict collision guards" begin
    for case in native["case"]
        println("HYBRID_PUBLICATION_CASE ",case["name"]);flush(stdout)
        m=startswith(case["name"],"connected") ? adjacent_model() : disjoint_model()
        input=actual_lower(m,case);saved=deepcopy(input.mesh)
        options=M._model_default_mesh1d_options(m,"native batch publication")
        options=M._ModelMesh1DOptions((name===:optimize ? Int(case["optimize"]) : getfield(options,name)
            for name in fieldnames(M._ModelMesh1DOptions))...)
        options.random_state[]=UInt32(1)
        births=M._ModelNativeVolumeBirths(UInt64(case["node_clock"]),UInt64(case["current_lower_element_high_water"]))
        batch=M._model_native_hybrid_batch_production(m,[1,2],input,options,births,"native batch publication")
        initial=assembly(input,case);part=deepcopy(initial)
        receipt=A._higher_native_hybrid_batch_receipt(part,m,batch,input,[2,1],"native batch publication")
        @test receipt.batch===batch
        @test state(part)==state(initial)
        @test receipt.node_base==UInt64(case["node_clock"])
        @test receipt.cell_base==UInt64(case["current_lower_element_high_water"])
        @test A._higher_native_hybrid_batch_current(receipt,part,m,"native batch publication")===nothing
        production=batch.products[1];class=classification(production,1,input)
        before=state(part);random=options.random_state[]
        @test_throws ArgumentError A._higher_native_hybrid_batch_receipt(part,m,batch,deepcopy(input),[1,2],"foreign authority")
        @test_throws ArgumentError A._higher_native_hybrid_batch_receipt(part,m,batch,input,[1],"incomplete Region batch")
        duplicate=deepcopy(batch)
        first_owned=findfirst(uid->uid>receipt.node_base,duplicate.products[1].boundary_node_tags)
        second_owned=findfirst(uid->uid>receipt.node_base,duplicate.products[2].boundary_node_tags)
        duplicate.products[2].boundary_node_tags[second_owned]=duplicate.products[1].boundary_node_tags[first_owned]
        @test_throws ArgumentError A._higher_native_hybrid_batch_receipt(part,m,duplicate,duplicate.input,[1,2],"aliased Region node birth")
        duplicate=deepcopy(batch);duplicate.products[2].element_tags[1]=duplicate.products[1].element_tags[1]
        @test_throws ArgumentError A._higher_native_hybrid_batch_receipt(part,m,duplicate,duplicate.input,[1,2],"aliased Region cell birth")
        @test_throws ArgumentError A._higher_native_hybrid_batch_binding(receipt,part,m,deepcopy(production.mesh),class,3,1,
            production.boundary_node_tags,receipt.node_base,batch.node_high_water,production.element_tags,batch.element_high_water,"naked copy")
        @test_throws ArgumentError A._higher_native_hybrid_batch_binding(receipt,part,m,production.mesh,class,3,1,
            nothing,receipt.node_base,batch.node_high_water,production.element_tags,batch.element_high_water,"missing rows")
        @test_throws ArgumentError A._higher_native_hybrid_batch_cell_births(receipt,part,1,reverse(production.element_tags),batch.element_high_water,"reordered cells")
        @test_throws ArgumentError A._higher_native_hybrid_batch_cell_births(receipt,part,1,fill(first(production.element_tags),length(production.element_tags)),batch.element_high_water,"duplicate cells")
        foreign_class=merge(class,(node_entities=copy(class.node_entities),))
        row=findfirst(uid->uid>receipt.node_base,production.boundary_node_tags);@assert row!==nothing
        foreign_class.node_entities[row]=(3,Int32(2))
        @test_throws ArgumentError publish!(part,m,production,foreign_class,1,receipt)
        @test state(part)==before
        @test options.random_state[]==random
        changed=deepcopy(receipt);changed.batch.input.mesh.coords[1,1]=nextfloat(changed.batch.input.mesh.coords[1,1])
        @test_throws ArgumentError A._higher_native_hybrid_batch_current(changed,part,m,"one ULP actual Source")
        changed=deepcopy(receipt);owned=findfirst(uid->uid>receipt.node_base,changed.batch.products[1].boundary_node_tags)
        changed.batch.products[1].mesh.coords[1,owned]=nextfloat(changed.batch.products[1].mesh.coords[1,owned])
        @test_throws ArgumentError A._higher_native_hybrid_batch_current(changed,part,m,"one ULP product")
        changed=deepcopy(receipt);reverse!(changed.batch.products[1].element_tags)
        @test_throws ArgumentError A._higher_native_hybrid_batch_current(changed,part,m,"changed factory order")
        changed=deepcopy(receipt);reverse!(changed.batch.products[1].registration_rows)
        @test_throws ArgumentError A._higher_native_hybrid_batch_current(changed,part,m,"changed registration receipt")
        changed=deepcopy(receipt);delete!(changed.batch.products,2)
        @test_throws ArgumentError A._higher_native_hybrid_batch_current(changed,part,m,"missing producer")
        changed=deepcopy(receipt)
        constructor=first(values(changed.batch.input.surface_constructor_inputs))
        block_index=findfirst(block->size(block.nodes,2)>1,constructor.blocks);@assert block_index!==nothing
        block=constructor.blocks[block_index];block.nodes[:,:]=block.nodes[:,end:-1:1]
        reverse!(block.tags);reverse!(constructor.entity_data.external_element_tags[block_index])
        @test_throws ArgumentError A._higher_native_hybrid_batch_current(changed,part,m,"changed owned constructor order")
        changed=deepcopy(receipt);parameter=findfirst(node->node.parameter!==nothing,changed.input_nodes)
        @assert parameter!==nothing
        changed.batch.input.mesh.entity_data.node_parametric[parameter][1]=nextfloat(changed.batch.input.mesh.entity_data.node_parametric[parameter][1])
        @test_throws ArgumentError A._higher_native_hybrid_batch_current(changed,part,m,"changed raw Source UV")
        @test state(part)==before
        @test publish!(part,m,production,class,1,receipt) isa Vector{Int32}
        @test part.node_max==batch.node_high_water
        @test part.cell_max==batch.element_high_water
        @test A._higher_native_hybrid_batch_current(receipt,part,m,"first Region published")===nothing
        before=state(part)
        @test_throws ArgumentError publish!(part,m,production,class,1,receipt)
        @test state(part)==before
        second=batch.products[2];second_class=classification(second,2,input)
        generic=deepcopy(part)
        error=try
            A._higher_append_part!(generic,m,second.mesh,second_class,3,2,Set{Tuple{Int,Int32}}();
                source_row_tags=second.boundary_node_tags,source_node_base=receipt.node_base,
                source_node_high_water=batch.node_high_water,source_cell_tags=second.element_tags,
                source_cell_high_water=batch.element_high_water)
            nothing
        catch exception
            exception
        end
        @test error isa ArgumentError
        @test occursin("native cell births collide with retained identities",sprint(showerror,error))
        @test publish!(part,m,second,second_class,2,receipt) isa Vector{Int32}
        @test A._higher_native_hybrid_batch_current(receipt,part,m,"complete Regions published")===nothing
        @test part.node_max==UInt64(case["node_high_water"])
        @test part.cell_max==UInt64(case["element_high_water"])
        for expected in case["expected"]
            tag=Int32(expected["owner"][2]);cells=filter(cell->cell.owner==(3,tag),part.cells)
            @test [cell.tag for cell in cells]==UInt64[uid for block in expected["cell"] for uid in block["labels"]]
            @test [Tuple(part.nodes[row].tag for row in cell.nodes) for cell in cells]==Tuple.(UInt64.(nodes) for block in expected["cell"] for nodes in block["nodes"])
            @test [cell.msh for cell in cells]==[Int(block["msh"]) for block in expected["cell"] for _ in block["labels"]]
            @test UInt64[node.tag for node in part.nodes if node.owner==(3,tag)]==UInt64.(expected["node_tags"])
            positions=part.node_positions
            observed=reduce(hcat,collect(part.nodes[positions[UInt64(uid)]].point) for uid in expected["node_tags"])
            @test isequal(observed,reduce(hcat,hexpoint.(expected["xyz_hex"])))
        end
        @test allunique(node.tag for node in part.nodes)
        @test allunique(cell.tag for cell in part.cells)
        @test isequal(input.mesh.coords,saved.coords)
        @test isequal(input.mesh.entity_data.node_parametric,saved.entity_data.node_parametric)
        @test input.mesh.entity_data.external_node_tags==saved.entity_data.external_node_tags
        @test input.mesh.entity_data.external_element_tags==saved.entity_data.external_element_tags
        @test [block.nodes for block in input.mesh.blocks]==[block.nodes for block in saved.blocks]
        println("HYBRID_PUBLICATION_FINAL ",case["name"]," clocks=",(part.node_max,part.cell_max));flush(stdout)
    end
end
end
