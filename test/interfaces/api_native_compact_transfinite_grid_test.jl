using Test,TOML,Tessella
const TCGA=Tessella.Model
const COMPACT_TRANSPORT=joinpath(@__DIR__,"..","artifacts","native_compact_transfinite_grid_4152_v1.toml")
# Native compact prisms retain separate Source Surface and newly evaluated
# Region diagonal records. The receipt authenticates logical grid events;
# it preserves their literal UIDs, owners, coordinates, cells and clocks.
@testset "Actual compact row bindings and recorded native exterior events" begin
    cases=TOML.parsefile(COMPACT_TRANSPORT)["case"]
    @test length(cases)==5
    for case in cases
        n=case["n"];counts=(n,n,n)
        grids=[UInt64.(ids) for ids in case["ids"]]
        face_ids=Any[grids[1],grids[2],grids[3],nothing,grids[4],grids[5]]
        count=length(case["native_node_ids"])
        dummy=TCGA.Mesh(zeros(3,count))
        tags,retained=TCGA._model_transfinite_volume_bindings(dummy,face_ids,counts,"actual compact";prism=true,compact=true)
        index=TCGA._model_transfinite_volume_index(face_ids,counts,"actual compact";prism=true,compact=true)
        actual=copy(tags);seen=falses(count);births=Int32[]
        for i in 0:n,j in 0:n,k in 0:n
            row=TCGA._model_transfinite_volume_row(index,i,j,k)
            seen[row] && continue
            seen[row]=true
            retained[row] || push!(births,row)
        end
        @test all(seen)
        @test length(births)==length(case["fresh_uids"])
        actual[births].=UInt64.(case["fresh_uids"])
        @test allunique(actual) && Set(actual)==Set(UInt64.(case["native_node_ids"]))
        native_rows=Dict(UInt64(uid)=>row for (row,uid) in enumerate(case["native_node_ids"]))
        coords=Matrix{Float64}(undef,3,count)
        for row in 1:count,axis in 1:3
            bits=case["native_coord_bits"][native_rows[actual[row]]][axis]
            coords[axis,row]=reinterpret(Float64,parse(UInt64,bits;base=16))
        end
        row_by_uid=Dict(uid=>Int32(row) for (row,uid) in enumerate(actual))
        tets=reduce(hcat,[Int32[row_by_uid[UInt64(uid)] for uid in cell] for cell in case["native_tets"]])
        for column in axes(tets,2)
            volume=TCGA.tet_signed_volume((coords[:,tets[i,column]] for i in 1:4)...)
            volume<0 && ((tets[1,column],tets[2,column])=(tets[2,column],tets[1,column]))
        end
        mesh=TCGA.Mesh(coords;tets=tets)
        @test TCGA.validate(mesh).ok
        @test length(case["native_tets"])==3n^3
        faces=TCGA._ModelTransfiniteGridFace[]
        panels=Vector{UInt64}[];source_cells=Tuple[]
        for (part,owner) in enumerate(case["owners"])
            cells=Tuple{Int,Tuple}[(2,Tuple(UInt64.(cell))) for cell in case["source_cells"][part]]
            labels=[UInt64.(case["source_labels"][part])]
            face=TCGA._ModelTransfiniteGridFace(Int32(owner),grids[part],(n,n),cells,labels)
            push!(faces,face)
            local_panels=TCGA._model_transfinite_grid_panels(face.ids,face.extents,"native compact")
            @test length(TCGA._model_transfinite_grid_partition(local_panels,last.(cells),"native compact"))==length(cells)
            append!(panels,local_panels);append!(source_cells,last.(cells))
        end
        TCGA._model_transfinite_grid_partition(panels,source_cells,"native compact")
        receipt=TCGA._model_transfinite_compact_alias_events(mesh,grids[3],counts,5,"owned compact")
        aliases=TCGA._model_transfinite_compact_alias_binding(receipt,faces,mesh,tags,"owned compact")
        @test length(aliases)==(n-1)^2
        @test Set(actual[row] for row in keys(aliases))==Set(UInt64.(case["volume_owned_exterior_uid"]))
        boundary=[Tuple(row_by_uid[UInt64(uid)] for uid in cell) for cell in case["native_boundary"]]
        volume_cells=[Tuple(iszero(tags[row]) ? aliases[row] : tags[row] for row in cell) for cell in boundary]
        @test length(TCGA._model_transfinite_grid_partition(panels,volume_cells,"native compact"))==length(boundary)
        # Authenticate the same aliases through the actual Source grid receipt
        # and actual native lower products, including ordered cells and labels.
        lower=[row for (row,owner) in enumerate(case["native_node_owners"]) if owner[1]<3]
        lower_tags=UInt64.(case["native_node_ids"][lower])
        lower_lookup=Dict(uid=>Int32(row) for (row,uid) in enumerate(lower_tags))
        lower_owners=Tuple{Int,Int32}[(case["native_node_owners"][row][1],Int32(case["native_node_owners"][row][2])) for row in lower]
        entity_owners=unique(vcat(lower_owners,Tuple{Int,Int32}[(2,Int32(owner)) for owner in case["owners"]]))
        entities=Dict(owner=>TCGA.Elements.MixedEntity(owner[1],owner[2],(0.,0.,0.,0.,0.,0.)) for owner in entity_owners)
        input_coords=Matrix{Float64}(undef,3,length(lower))
        for (target,source) in enumerate(lower),axis in 1:3
            input_coords[axis,target]=reinterpret(Float64,parse(UInt64,case["native_coord_bits"][source][axis];base=16))
        end
        parameters=Union{Nothing,Vector{Float64}}[parse.(Float64,case["native_node_parameter_hex"][row]) for row in lower]
        products=Dict{Tuple{Int,Int32},TCGA._ModelMeshProduct}()
        blocks=TCGA.ElementBlock[];groups=Vector{Int32}[];labels=Vector{UInt64}[]
        for (part,owner) in enumerate(case["owners"])
            uid_rows=unique(grids[part]);local_lookup=Dict(uid=>Int32(row) for (row,uid) in enumerate(uid_rows))
            cells=case["source_cells"][part]
            global_cells=reduce(hcat,[Int32[lower_lookup[UInt64(uid)] for uid in cell] for cell in cells])
            local_cells=reduce(hcat,[Int32[local_lookup[UInt64(uid)] for uid in cell] for cell in cells])
            cell_tags=UInt64.(case["source_labels"][part]);key=(2,Int32(owner))
            products[key]=TCGA._ModelMeshProduct(key,Int32[lower_lookup[uid] for uid in uid_rows],
                [TCGA.ElementBlock(2,local_cells)],[copy(cell_tags)])
            push!(blocks,TCGA.ElementBlock(2,global_cells));push!(groups,fill(Int32(owner),length(cell_tags)));push!(labels,cell_tags)
        end
        data=TCGA.Elements._native_mixed_entity_data(entities;node_entities=lower_owners,node_parametric=parameters,
            external_node_tags=lower_tags,block_entities=groups,external_element_tags=labels)
        lower_mesh=TCGA.MixedMesh(input_coords,blocks;entity_data=data,elementary_entities=groups)
        authority=TCGA._ModelMeshAuthority(lower_mesh,products,Set(keys(products)),UInt64.(1:length(lower)),
            Dict{Tuple{Int,Int32},Union{TCGA.Mesh,TCGA.MixedMesh}}(),lower_mesh,Int32.(1:length(lower)),
            Dict{Tuple{Int,Int32},Vector{Int32}}(),Vector{Union{Nothing,UInt64}}[],
            Dict{Tuple{Int,Int32},TCGA.MixedMesh}(),nothing)
        grid_receipt=TCGA._model_transfinite_grid_receipt(authority,face_ids,[3,4,5,0,1,2],counts,"actual receipt";mesh,compact=true)
        supports=TCGA._model_transfinite_grid_boundary(grid_receipt,authority,tags,boundary,"actual receipt";mesh)
        @test length(supports)==length(boundary)
        @test Set(values(supports))==Set((2,Int32(owner)) for owner in case["owners"])
        @test_throws ArgumentError TCGA._model_transfinite_grid_boundary(grid_receipt,authority,tags,boundary,"strict")
        @test_throws ArgumentError TCGA._model_transfinite_grid_boundary(grid_receipt,authority,tags,boundary,"strict";mesh=TCGA.Mesh(nextfloat.(coords);tets=tets))
        # Source cells and labels are rechecked by the owning grid receipt.
        first_product=products[(2,Int32(3))];saved=first_product.element_tags[1][1]
        first_product.element_tags[1][1]=saved+1
        @test_throws ArgumentError TCGA._model_transfinite_grid_boundary(grid_receipt,authority,tags,boundary,"strict";mesh)
        first_product.element_tags[1][1]=saved
        @test any(row->iszero(tags[row]),Iterators.flatten(boundary))
        @test TCGA._model_transfinite_compact_alias_binding(nothing,faces,mesh,tags,"strict")==Dict{Int32,UInt64}()
        @test_throws ArgumentError TCGA._model_transfinite_compact_alias_binding(receipt,vcat(faces,[faces[1]]),mesh,tags,"strict")
        @test_throws ArgumentError TCGA._model_transfinite_compact_alias_binding(receipt,faces,nothing,tags,"strict")
        bad=TCGA._ModelTransfiniteGridAliasReceipt(counts,Int32(99),receipt.events)
        @test_throws ArgumentError TCGA._model_transfinite_compact_alias_binding(bad,faces,mesh,tags,"strict")
        bad=TCGA._ModelTransfiniteGridAliasReceipt(counts,Int32(5),receipt.events[1:end-1])
        @test_throws ArgumentError TCGA._model_transfinite_compact_alias_binding(bad,faces,mesh,tags,"strict")
        event=first(receipt.events)
        wrong=TCGA._ModelTransfiniteGridAlias(event.row,event.source_uid+1,event.point)
        bad=TCGA._ModelTransfiniteGridAliasReceipt(counts,Int32(5),vcat([wrong],receipt.events[2:end]))
        @test_throws ArgumentError TCGA._model_transfinite_compact_alias_binding(bad,faces,mesh,tags,"strict")
        changed=copy(coords);changed[1,event.row]=nextfloat(changed[1,event.row])
        @test_throws ArgumentError TCGA._model_transfinite_compact_alias_binding(receipt,faces,TCGA.Mesh(changed;tets=tets),tags,"strict")
        changed_tags=copy(tags);changed_tags[event.row]=event.source_uid
        @test_throws ArgumentError TCGA._model_transfinite_compact_alias_binding(receipt,faces,mesh,changed_tags,"strict")
        @test_throws ArgumentError TCGA._model_transfinite_compact_alias_events(mesh,grids[3],(n,n+1,n),5,"strict")
        @test_throws ArgumentError TCGA._model_transfinite_compact_alias_events(mesh,grids[3],(typemax(Int),typemax(Int),n),5,"strict")
    end
end
