# GModel's default optimizer visits complete regions. niter and dimTags do
# not select work for these three method aliases; force admits discrete CAD.
function _native_optimize_selection(dim_tags,caller)
    (dim_tags isa AbstractVector || dim_tags isa Tuple) || throw(ArgumentError(
        "$caller: dim_tags must be a vector or tuple"))
    result=Tuple{Int,Int}[]
    for pair in dim_tags
        (pair isa Tuple || pair isa AbstractVector) && length(pair)==2 ||
            throw(ArgumentError("$caller: each dim_tags entry must contain two integers"))
        push!(result,(_mesh_query_integer(pair[1],caller,"dimension"),
                     _mesh_query_integer(pair[2],caller,"entity tag")))
    end
    return result
end

function _native_optimize_fully_discrete(m,region)
    record=get(m.discrete,(3,region),nothing)
    record!==nothing && !record.geometry_parametrized || return false
    for (dimension,face) in record.boundary
        dimension==2 || continue
        surface=get(m.discrete,(2,abs(face)),nothing)
        surface!==nothing && !surface.geometry_parametrized || return false
        for (curve_dimension,curve) in surface.boundary
            curve_dimension==1 || continue
            edge=get(m.discrete,(1,abs(curve)),nothing)
            edge!==nothing && !Model._model_discrete_curve_parametrized(edge) || return false
        end
    end
    return true
end

function _native_optimize_region_skipped(m,region,force)
    !force && _native_optimize_fully_discrete(m,region) && return true
    haskey(m.meshing.transfinite_volumes,region) && return true
    sweep=get(m.meshing.extrude,(3,region),nothing)
    source=get(m.meshing.extrude_sources,(3,region),nothing)
    return sweep!==nothing && !isempty(sweep.layers) && source!==nothing && source[1]==2
end

function _native_optimize_determinant(a,b,c,d)
    first=ntuple(i->b[i]-a[i],3)
    second=ntuple(i->c[i]-a[i],3)
    third=ntuple(i->d[i]-a[i],3)
    return first[1]*(second[2]*third[3]-third[2]*second[3])-
        second[1]*(first[2]*third[3]-third[2]*first[3])+
        third[1]*(first[2]*second[3]-second[2]*first[3])
end

function _native_optimize_volumes_positive!(assembly)
    for cell in assembly.cells
        cell.owner[1]==3 || continue
        family=msh_spec(cell.msh).family
        slots=family===:tet || family===:pri ? (1,2,3,4) :
              family===:hex || family===:pyr ? (1,2,4,5) : nothing
        slots===nothing && continue
        points=ntuple(i->assembly.nodes[cell.nodes[slots[i]]].point,4)
        _native_optimize_determinant(points...)<0 &&
            _record_reverse_connectivity!(Int32(cell.msh),cell.nodes)
    end
    return nothing
end

function _native_optimize_embedded(m,region,cells,positions)
    faces=Set{NTuple{3,Int32}}()
    edges=Set{Tuple{Int32,Int32}}()
    embedded=Set((dimension,Int32(tag))
        for (dimension,tag) in get(m.embeds,(3,region),NTuple{2,Int}[]))
    isempty(embedded) && return faces,edges
    for cell in cells
        cell.owner in embedded || continue
        family=msh_spec(cell.msh).family
        if family===:lin
            first=get(positions,cell.nodes[1],Int32(0))
            last=get(positions,cell.nodes[2],Int32(0))
            first!=0 && last!=0 && push!(edges,first<last ? (first,last) : (last,first))
        elseif family===:tri
            nodes=ntuple(i->get(positions,cell.nodes[i],Int32(0)),3)
            all(!iszero,nodes) && push!(faces,Tuple(sort(collect(nodes))))
        end
    end
    return faces,edges
end

function _native_optimize_region!(assembly,m,region,cells;netgen=false,relocate=true)
    rows=Int32[];positions=Dict{Int32,Int32}()
    for cell in cells,slot in 1:4
        node=cell.nodes[slot]
        if !haskey(positions,node)
            push!(rows,node);positions[node]=Int32(length(rows))
        end
    end
    coords=Matrix{Float64}(undef,3,length(rows))
    movement=falses(length(rows))
    for (row,node) in enumerate(rows)
        coords[:,row].=assembly.nodes[node].point
        movement[row]=assembly.nodes[node].owner[1]==3
    end
    connectivity=Matrix{Int32}(undef,4,length(cells))
    numbers=Vector{UInt64}(undef,length(cells))
    physical=Vector{Int32}(undef,length(cells))
    original=Dict{UInt64,_Dim01Cell}()
    for (column,cell) in enumerate(cells)
        connectivity[:,column]=Int32[positions[cell.nodes[slot]] for slot in 1:4]
        numbers[column]=cell.tag;physical[column]=cell.physical
        original[cell.tag]=cell
    end
    faces,edges=_native_optimize_embedded(m,region,assembly.cells,positions)
    product=if netgen
        Model._netgen_optimize_tetrahedra(Mesh(coords;tets=connectivity,tet_tag=physical);
            element_numbers=numbers,element_high_water=assembly.cell_max,movable=movement,
            embedded_faces=faces,embedded_edges=edges,positive=false)
    else
        Model._gmsh_optimize_tetrahedra(Mesh(coords;tets=connectivity,tet_tag=physical);
            element_numbers=numbers,element_high_water=assembly.cell_max,movable=movement,
            threshold=OPTIONS["Mesh.OptimizeThreshold"],embedded_faces=faces,
            embedded_edges=edges,positive=false,relocate=relocate)
    end
    assembly.cell_max=product.element_high_water
    for (row,node) in enumerate(rows)
        previous=assembly.nodes[node]
        assembly.nodes[node]=_Dim01Node(previous.tag,
            ntuple(i->product.mesh.coords[i,row],3),previous.owner,previous.parameter)
    end
    replacement=_Dim01Cell[]
    for column in axes(product.mesh.tets,2)
        tag=product.element_numbers[column]
        nodes=Int32[rows[product.mesh.tets[slot,column]] for slot in 1:4]
        survivor=get(original,tag,nothing)
        if survivor===nothing
            push!(replacement,_Dim01Cell(4,(3,Int32(region)),nodes,tag,tag,
                                        product.mesh.tet_tag[column]))
        else
            # Existing curved interpolation vertices belong to the surviving
            # MElement. Newly allocated optimizer tetrahedra are linear.
            survivor.nodes[1:4]==nodes || throw(ArgumentError(
                "$(assembly.caller): optimizer changed a surviving tetrahedron's primary slots"))
            push!(replacement,survivor)
        end
    end
    return replacement
end

function _native_default_optimize_locked!(force,pairs,caller;netgen=false)
    model=_model_locked();cached=LAST_MESH[]
    class=cached===nothing ? nothing : _cached_classification_locked(cached)
    # The legacy simplex carrier stores its actual quadratic cells and edited
    # interpolation vertices in an overlay. Native optimization preserves
    # surviving MElements and still renumbers after a surface-only no-op.
    actual,actual_class=_dim01_actual_cache(cached,class;physical_names=model.physical_names)
    isempty(pairs) || @warn "Entity selection in the default mesh optimizer is ignored"
    source=Model._model_mesh_transaction_copy(model)
    assembly=_Dim01Assembly(_Dim01Node[],_Dim01Cell[],Dict{UInt64,Int32}(),
        NODE_TAG_MAX[],ELEMENT_TAG_MAX[],typemax(Int32),typemax(Int32),caller)
    mapping=_dim01_ingest!(assembly,source,actual,actual_class;preserve_source_order=true)
    # OptimizeMesh checks Msg errors before work. GModel's final renumbering
    # still executes when that error guard skips the optimizer.
    admitted=iszero(OPTIONS["General.AbortOnError"]) || NATIVE_LOOKUP_QUERY_ERRORS[]==0
    if admitted
        regions=Dict{Int,Vector{_Dim01Cell}}()
        # Upstream's hybrid gates inspect only hexahedra, prisms and pyramids:
        # optimizeMeshGRegion still swaps hybrid edges but skips its relocate
        # pass, while optimizeMeshGRegionNetgen skips the region outright.
        hybrid=Set{Int}()
        for cell in assembly.cells
            cell.owner[1]==3 || continue
            family=msh_spec(cell.msh).family
            if family===:tet
                push!(get!(()->_Dim01Cell[],regions,Int(cell.owner[2])),cell)
            elseif family in (:hex,:pri,:pyr)
                push!(hybrid,Int(cell.owner[2]))
            end
        end
        replacements=Dict{Int,Vector{_Dim01Cell}}()
        for (_,region) in model_entities(source,3)
            _native_optimize_region_skipped(source,region,force) && continue
            if netgen
                if region in hybrid
                    _native_diagnostic_emit_locked(:info,
                        "Skipping Netgen optimizer for hybrid mesh")
                    continue
                end
                _native_diagnostic_emit_locked(:info,"Optimizing volume $region")
            end
            cells=get(regions,region,nothing)
            cells===nothing || (replacements[region]=
                _native_optimize_region!(assembly,source,region,cells;
                    netgen=netgen,relocate=!(region in hybrid)))
        end
        result=_Dim01Cell[];emitted=Set{Int}()
        for cell in assembly.cells
            region=Int(cell.owner[2])
            if cell.owner[1]==3 && msh_spec(cell.msh).family===:tet && haskey(replacements,region)
                if !(region in emitted)
                    append!(result,replacements[region]);push!(emitted,region)
                end
            else
                push!(result,cell)
            end
        end
        assembly.cells=result
        _native_optimize_volumes_positive!(assembly)
    end
    isempty(assembly.nodes) && isempty(assembly.cells) && return nothing
    renumber=!iszero(OPTIONS["Mesh.Renumber"])
    plan=_dim01_finish!(assembly,source,model_entities(source),Set{Tuple{Int,Int32}}(),
        falses(length(assembly.nodes)),renumber,false;
        save_all=!iszero(OPTIONS["Mesh.SaveAll"]),source_mesh=actual,
        source_node_positions=mapping)
    if cached isa MixedMesh && (!isempty(cached.ancillary_sections) ||
        !isempty(cached.data_sections) || cached.partition_data!==nothing)
        cached.entity_data!==nothing ||
            (isempty(cached.data_sections) && cached.partition_data===nothing) ||
            throw(ArgumentError("$caller: retained file rows require explicit public-ID metadata"))
        sections=_reclassify_data_sections(cached,plan.mesh,plan,caller)
        partition=_reclassify_partition_data(cached,plan.mesh,plan,caller)
        preserved=MixedMesh(plan.mesh.coords,plan.mesh.blocks;
            physical_names=plan.mesh.physical_names,entity_data=plan.mesh.entity_data,
            elementary_entities=plan.mesh.elementary_entities,periodic_links=plan.mesh.periodic_links,
            ancillary_sections=cached.ancillary_sections,data_sections=sections,partition_data=partition)
        plan=merge(plan,(mesh=preserved,class=_mixed_rebind_class(plan.class,preserved),))
    end
    _commit_tagged_plan_locked!(plan;preserve_visibility=true,
        element_tag_map=plan.visibility_tag_map)
    # With renumbering disabled native lookup memo retains misses for births.
    renumber && _native_lookup_invalidate_locked!()
    return nothing
end
