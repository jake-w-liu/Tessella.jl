module APINativeTransfiniteGridTests
using Test,Tessella,TOML,SHA
const L=Tessella.Model
const E=Tessella.Elements
const TRANSPORT=joinpath(@__DIR__,"..","artifacts","native_transfinite_grid_boundary_4152_v1.toml")
const PRIMARY=joinpath(@__DIR__,"..","artifacts","native_transfinite_grid_boundary_4152_v1.json")
# The cells below are literal native products. Coordinates are immaterial to
# this incidence proof; the public classifier separately certifies geometry.
function grid_ids(panels)
    ids=zeros(UInt64,16)
    for (index,panel) in enumerate(panels)
        i=(index-1)÷3;j=(index-1)%3
        for (slot,(x,y)) in enumerate(((i,j),(i+1,j),(i+1,j+1),(i,j+1)))
            target=x+1+4*y
            @assert iszero(ids[target]) || ids[target]==panel[slot]
            ids[target]=panel[slot]
        end
    end
    @assert all(!iszero,ids)
    return ids
end
function fixture(cases)
    grids=grid_ids.([case["panels"] for case in cases])
    source_tags=sort!(unique(reduce(vcat,grids)))
    lookup=Dict(uid=>Int32(row) for (row,uid) in enumerate(source_tags))
    products=Dict{Tuple{Int,Int32},L._ModelMeshProduct}()
    blocks=L.ElementBlock[];groups=Vector{Int32}[];tags=Vector{UInt64}[]
    next=UInt64(1)
    for (owner,case) in enumerate(cases)
        ids=grids[owner]
        rows=Int32[lookup[uid] for uid in ids]
        local_lookup=Dict(uid=>Int32(row) for (row,uid) in enumerate(ids))
        owned=L.ElementBlock[];owned_tags=Vector{UInt64}[]
        for (msh,n) in ((2,3),(3,4))
            cells=[cell for cell in case["surface"] if length(cell)==n]
            isempty(cells) && continue
            local_cells=reduce(hcat,[Int32[local_lookup[uid] for uid in cell] for cell in cells])
            cells_global=reduce(hcat,[Int32[lookup[uid] for uid in cell] for cell in cells])
            cell_tags=collect(next:next+UInt64(length(cells))-1)
            next+=UInt64(length(cells))
            push!(owned,L.ElementBlock(msh,local_cells));push!(owned_tags,cell_tags)
            push!(blocks,L.ElementBlock(msh,cells_global));push!(groups,fill(Int32(owner),length(cells)))
            push!(tags,copy(cell_tags))
        end
        key=(2,Int32(owner))
        products[key]=L._ModelMeshProduct(key,rows,owned,owned_tags)
    end
    entities=Dict((2,i)=>E.MixedEntity(2,i,(0.,0.,0.,0.,0.,0.)) for i in 1:6)
    data=E._native_mixed_entity_data(entities;
        node_entities=fill((2,Int32(1)),length(source_tags)),
        node_parametric=Union{Nothing,Vector{Float64}}[nothing for _ in source_tags],external_node_tags=source_tags,
        block_entities=groups,external_element_tags=tags)
    mesh=L.MixedMesh(zeros(3,length(source_tags)),blocks;entity_data=data,elementary_entities=groups)
    authority=L._ModelMeshAuthority(mesh,products,Set(keys(products)),UInt64.(1:length(source_tags)),
        Dict{Tuple{Int,Int32},Union{L.Mesh,L.MixedMesh}}(),mesh,Int32.(1:length(source_tags)),
        Dict{Tuple{Int,Int32},Vector{Int32}}(),Vector{Union{Nothing,UInt64}}[],
        Dict{Tuple{Int,Int32},L.MixedMesh}())
    exterior=Set(Tuple(sort(Int32[lookup[uid] for uid in cell]))
        for case in cases for cell in case["volume"])
    return authority,grids,source_tags,exterior
end
@testset "Actual product receipt certifies independent native TF diagonals" begin
    data=TOML.parsefile(TRANSPORT)
    @test data["primary_sha256"]==bytes2hex(sha256(read(PRIMARY)))
    names=unique([case["name"] for case in data["case"]])
    @test length(names)==6
    for name in names
        cases=sort([case for case in data["case"] if case["name"]==name];by=case->case["owner"])
        input,grids,source_tags,exterior=fixture(cases)
        receipt=L._model_transfinite_grid_receipt(input,grids,Int32.(1:6),(3,3,3),"native TF")
        result=L._model_transfinite_grid_boundary(receipt,input,source_tags,exterior,"native TF")
        @test Set(keys(result))==exterior
        @test Set(values(result))==Set((2,Int32(i)) for i in 1:6)
        @test length(result)==sum(length(case["volume"]) for case in cases)
        @test all(Set(face.ids)==Set(grids[i]) for (i,face) in enumerate(receipt.faces))
        @test_throws ArgumentError L._model_transfinite_grid_boundary(receipt,deepcopy(input),source_tags,exterior,"stale")
        missing=setdiff(exterior,Set([first(exterior)]))
        @test_throws ArgumentError L._model_transfinite_grid_boundary(receipt,input,source_tags,missing,"missing")
        source_zero=copy(source_tags);source_zero[first(first(exterior))]=0
        @test_throws ArgumentError L._model_transfinite_grid_boundary(receipt,input,source_zero,exterior,"unbound")
        product=input.products[(2,Int32(1))]
        product.element_tags[1][1]+=1
        @test_throws ArgumentError L._model_transfinite_grid_boundary(receipt,input,source_tags,exterior,"retagged")
        product.element_tags[1][1]-=1
        old=product.blocks[1].nodes[1,1]
        product.blocks[1].nodes[1,1]=product.blocks[1].nodes[2,1]
        @test_throws ArgumentError L._model_transfinite_grid_boundary(receipt,input,source_tags,exterior,"edited")
        product.blocks[1].nodes[1,1]=old
        @test L._model_transfinite_grid_boundary(receipt,input,source_tags,exterior,"restored")==result
        grids[1][1]=grids[1][2]
        @test_throws ArgumentError L._model_transfinite_grid_receipt(input,grids,Int32.(1:6),(3,3,3),"foreign grid")
    end
end
end
