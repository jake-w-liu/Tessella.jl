# Actual original-source frames and read-only malformed source preflight.

function b4_source_frames(f)
    M=Tessella.Model
    execution=B4.execute(f.source;dim=0);m=execution.model
    out=Int.(execution.lists["sweep"]);t=out[2]
    caller="independent B4 original source frames"
    completed=M._extrude_nonew_plan(deepcopy(m),t,caller)
    original=completed.source_mesh
    nodes=only(B4.blocks(original)).nodes
    params=M._extrude_gate(m,3,t);spec=m.meshing.extrude_specs[(3,t)]
    ncells=size(nodes,2)
    orders=[collect(1:ncells),collect(ncells:-1:1),circshift(collect(1:ncells),1)]
    for order in orders,rotation in 0:3
        cyclic=hcat((circshift(nodes[:,order[column]],mod(rotation+column-1,4))
            for column in 1:ncells)...)
        source=MixedMesh(original.coords,[ElementBlock(3,cyclic)])
        independent=B4.source_complex(source,f)
        # Proposed Fresh source mode is explicit; rectangular defaults remain
        # untouched. Actual source IDs and cyclic frames must be retained.
        data=M._extrude_nonew_rect_grid_source(m,f.surface,source,spec,caller;mode=:b4_strip)
        @test data.source_cells==[Tuple(Int32.(cell)) for cell in eachcol(cyclic)]
        @test data.source_coordinates==[Tuple(original.coords[:,node]) for node in 1:nnodes(original)]
        @test all(data.boundary_vertices) && length(independent.interior)==0
        catalog=M._extrude_nonew_b4_strip_plan(data,completed.catalog.levels,
            completed.catalog.layer_refs,true,caller)
        @test catalog.cell_counts==(2ncells,ncells*(f.intervals-1),0,5ncells)
        @test catalog.face_capacity==(4ncells+1)*f.intervals+15ncells
        @test catalog.diagonal_count==ncells
        M._extrude_nonew_rect_grid_product_certify(completed.sweep.cols,catalog,spec,caller)
        plan=M._extrude_nonew_b4_strip_finish(m,t,params,f.surface,source,catalog,
            completed.sweep.cols,out[1],Tuple(out[3:end]),caller)
        certificate=B4.certify(plan.sweep.volume,f,source)
        @test length(certificate.domains)==ncells*f.intervals
        @test length(certificate.centers)==ncells
        @test certificate.source.area==independent.area
    end
end

function b4_malformed_sources(f)
    M=Tessella.Model
    execution=B4.execute(f.source;dim=0);m=execution.model
    source=mesh_model_surface(m,f.surface)
    t=Int(execution.lists["sweep"][2]);spec=m.meshing.extrude_specs[(3,t)]
    nodes=only(B4.blocks(source)).nodes
    before=(repr(m),copy(source.coords),copy(nodes))
    missing=MixedMesh(source.coords,[ElementBlock(3,nodes[:,1:end-1])])
    duplicate=deepcopy(source);only(B4.blocks(duplicate)).nodes[:,2]=nodes[:,1]
    reversed=deepcopy(source);only(B4.blocks(reversed)).nodes[:,1]=reverse(nodes[:,1])
    coincident=deepcopy(source);coincident.coords[:,2]=coincident.coords[:,1]
    # Pure source preflight remains read-only on every unsuccessful request.
    for malformed in (missing,duplicate,reversed,coincident)
        @test_throws ArgumentError M._extrude_nonew_rect_grid_source(
            m,f.surface,malformed,spec,"malformed B4 source";mode=:b4_strip)
        @test repr(m)==before[1] && source.coords==before[2] && nodes==before[3]
    end
end

function b4_node_frames(f)
    M=Tessella.Model;g=B4.execute(f.source;dim=0);m=g.model
    out=Int.(g.lists["sweep"]);t=out[2];caller="independent B4 node permutations"
    completed=M._extrude_nonew_plan(deepcopy(m),t,caller)
    original=completed.source_mesh;nodes=only(B4.blocks(original)).nodes;n=nnodes(original)
    params=M._extrude_gate(m,3,t);spec=m.meshing.extrude_specs[(3,t)]
    orders=[collect(1:n),collect(n:-1:1),circshift(collect(1:n),1),
        vcat(collect(1:2:n),collect(2:2:n)),circshift(collect(1:n),n÷2),
        vcat(collect(n÷2:-1:1),collect(n:-1:n÷2+1))]
    @test length(Set(Tuple(order) for order in orders))==6
    for order in orders
        inverse=invperm(order)
        source=MixedMesh(original.coords[:,order],[ElementBlock(3,Int32.(inverse[nodes]))])
        columns=completed.sweep.cols[:,order]
        data=M._extrude_nonew_rect_grid_source(m,f.surface,source,spec,caller;mode=:b4_strip)
        @test data.source_coordinates==[Tuple(source.coords[:,node]) for node in 1:n]
        @test data.source_cells==[Tuple(Int32.(cell)) for cell in eachcol(only(B4.blocks(source)).nodes)]
        catalog=M._extrude_nonew_b4_strip_plan(data,completed.catalog.levels,
            completed.catalog.layer_refs,true,caller)
        M._extrude_nonew_rect_grid_product_certify(columns,catalog,spec,caller)
        plan=M._extrude_nonew_b4_strip_finish(m,t,params,f.surface,source,catalog,
            columns,out[1],Tuple(out[3:end]),caller)
        certificate=B4.certify(plan.sweep.volume,f,source)
        @test length(certificate.centers)==f.strip_length
        @test certificate.total==certificate.source.area*abs(B4.Q(f.height))
    end
end
