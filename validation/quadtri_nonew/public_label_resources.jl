# Actual native label plans over disjoint analytical unit-right triangles.
# The public rename is measured after warming; geometric proofs are uncharged.
function public_label_resource_fixture(count::Int)
    count>0 && count<=typemax(Int32)÷3 || throw(ArgumentError("invalid label resource cell count"))
    api=Tessella.API
    api.initialize()
    api.model.add_discrete_entity(2,1)
    coordinates=zeros(3,3count)
    for column in 1:count
        coordinates[:,3column-2:3column]=Float64[3column 3column+1 3column;0 0 1;0 0 0]
    end
    mesh=Mesh(coordinates;tris=reshape(Int32.(1:3count),3,count))
    class=api._MeshClassification(mesh,(2,Int32(1)),[(2,Int32(1))],
        fill((2,Int32(1)),3count),Dict{Tuple{Int,Int32},Vector{Int32}}(),
        Int32[],fill(Int32(1),count),Int32[])
    api._replace_mesh_cache_locked!(mesh,class)
    return mesh
end

function public_label_resource_rows()
    api=Tessella.API
    rows=Dict{String,Any}[]
    try
        for count in (1000,2000,4000)
            mesh=public_label_resource_fixture(count)
            coordinates=copy(mesh.coords);connectivity=copy(mesh.tris)
            api.mesh.renumber_elements()
            api.mesh.renumber_elements([1],[10001])
            api.mesh.renumber_elements()
            samples=[@timed(api.mesh.renumber_elements()) for _ in 1:3]
            selected=samples[argmin([sample.bytes for sample in samples])]
            public=api.LAST_MESH_CLASS[].public_tags
            @test api.LAST_MESH[]===mesh
            @test mesh.coords==coordinates && mesh.tris==connectivity
            @test public.node_tags==UInt64.(1:3count)
            @test public.element_tags==UInt64.(1:count)
            @test api.mesh.get_element_qualities([1,count],"volume")==[.5,.5]
            @test api.mesh.get_max_element_tag()==10000+count
            push!(rows,Dict("cells"=>count,"nodes"=>3count,"path"=>"native_element_labels",
                "allocated"=>selected.bytes,"allocation_count"=>Base.gc_alloc_count(selected.gcstats),
                "allocation_samples"=>[sample.bytes for sample in samples],
                "allocation_count_samples"=>[Base.gc_alloc_count(sample.gcstats) for sample in samples],
                "seconds"=>selected.time,"retained_bytes"=>Base.summarysize(public),
                "first_area"=>.5,"last_area"=>.5,"max_element_tag"=>10000+count,
                "coordinate_sha256"=>bytes2hex(sha256(reinterpret(UInt8,vec(mesh.coords)))),
                "connectivity_sha256"=>bytes2hex(sha256(reinterpret(UInt8,vec(mesh.tris))))))
        end
        for index in 2:length(rows)
            @test rows[index]["allocated"]<=2.15*rows[index-1]["allocated"]+65536
            @test rows[index]["retained_bytes"]<=2.15*rows[index-1]["retained_bytes"]+65536
        end
    finally
        api.finalize()
    end
    return rows
end

# Actual public label and lifecycle implementation bodies, including keyword
# wrappers and closure-call methods, are checked in the shared resource gate.
const PUBLIC_LABEL_RESOURCE_API_TARGETS=(
    :_api_rebind_class,:_append_native_boundary_nodes!,:_apply_mesh_order,
    :_cache_public_tags,:_commit_renumber_labels_locked!,:_dim01_actual_cache,
    :_dim01_add_node!,:_dim01_ingest!,:_dim01_ingest_record_cells!,
    :_dim01_quadratic!,:_dim01_reconcile_records!,:_edit_mesh_records_locked!,
    :_entity_tri_elements,:_generate,:_generation_bind_point_records!,
    :_generation_point_sources,:_generation_public_labels!,:_get_barycenters,
    :_get_basis_functions_orientation,:_get_basis_functions_orientation_for_element,:_get_element_by_coordinates,
    :_get_elements_by_coordinates,:_get_jacobian,:_get_jacobians,
    :_get_keys,:_get_keys_for_element,:_get_local_coordinates_in_element,
    :_get_nodes,:_get_periodic_keys,:_get_periodic_nodes,
    :_mesh_node_coords,:_mesh_public_node_coords,:_mixed_linear_cache,
    :_mixed_quadratic_cache,:_mixed_rebind_class,:_mixed_rebuild_metadata,
    :_mixed_refine_locked!,:_model_fresh_node_tag,:_native_dedup_node_identity_check,
    :_native_mixed_set_order_labels_plan,:_native_public_keys,:_native_public_node_remap,
    :_native_rebind_public_tags,:_native_set_order_labels_plan,:_public_point_locations,
    :_public_point_query_mesh,:_public_point_record_mesh,:_recombine,
    :_recombine_excluded_edges,:_reconcile_generated_records!,:_record_dedup_nodes!,
    :_record_edit_cache_nodes,:_record_nodes_of_type,:_refine,
    :_remove_duplicate_nodes,:_remove_public_elements_locked!,:_renumber_elements,
    :_renumber_nodes,:_set_order,:_split_quadrangles,
    :_tagged_mesh_max_tags,:_tagged_record_renumber!,:_tagged_renumber_lists,
    :_tagged_renumber_mapping,:_tagged_renumber_order,:_tagged_renumber_plan,
    :_validate_generated_cache_references,:_validate_record_cache_labels,:_validate_retained_cache_references,
)
