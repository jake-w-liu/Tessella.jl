module HigherIncrementalAuthorityChecks
using Test,Tessella
const A=Tessella.API
const M=Tessella.Model

# The Surface pass publishes each merged Surface to the next producer. It used
# to rebuild the authority of the whole mesh for every Surface, which made the
# pass quadratic in the number of Surfaces (400 unit faces: 9.8 s against
# Gmsh 4.15.2's 0.2 s). It now seeds one mesh log and appends each merge to it.
# The counters pin that contract; the verification switch rebuilds every
# appended authority in full and fails the generation on any difference.

function plane_grid(n)
    geo=A.model.geo
    point(i,j)=i*(n+1)+j+1
    for i in 0:n,j in 0:n
        geo.add_point(Float64(i),Float64(j),0.0;tag=point(i,j))
    end
    along(i,j)=i*n+j+1
    across(i,j)=(n+1)*n+i*(n+1)+j+1
    for i in 0:n,j in 0:n-1
        geo.add_line(point(i,j),point(i,j+1);tag=along(i,j))
    end
    for i in 0:n-1,j in 0:n
        geo.add_line(point(i,j),point(i+1,j);tag=across(i,j))
    end
    face=0
    for i in 0:n-1,j in 0:n-1
        face+=1
        geo.add_curve_loop([along(i,j),across(i,j+1),-along(i+1,j),-across(i,j)];tag=face)
        geo.add_plane_surface([face];tag=face)
    end
    geo.synchronize()
    return nothing
end

function generate_counted(n,size)
    A.initialize()
    try
        A.option("General.Terminal",0)
        plane_grid(n)
        A.option("Mesh.MeshSizeMin",size);A.option("Mesh.MeshSizeMax",size)
        seeds=M._MODEL_MESH_LOG_SEEDS[];appends=M._MODEL_MESH_LOG_APPENDS[]
        A.mesh.generate(2)
        nodes,_,_=A.mesh.get_nodes()
        _,tags,_=A.mesh.get_elements(2)
        return (seeds=M._MODEL_MESH_LOG_SEEDS[]-seeds,
            appends=M._MODEL_MESH_LOG_APPENDS[]-appends,
            nodes=length(nodes),triangles=sum(length,tags;init=0))
    finally
        A.finalize()
    end
end

@testset "the Surface pass extends one authority" begin
    plain=generate_counted(3,0.5)
    @test plain.seeds==1
    @test plain.appends==9
    M._MODEL_MESH_LOG_VERIFY[]=true
    try
        # Every append is compared with the full rebuild of the same state.
        verified=generate_counted(3,0.5)
        @test verified.seeds==1
        @test verified.appends==9
        @test (verified.nodes,verified.triangles)==(plain.nodes,plain.triangles)
    finally
        M._MODEL_MESH_LOG_VERIFY[]=false
    end
end

@testset "mesh log appends" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        plane_grid(1)
        A.option("Mesh.MeshSizeMin",0.5);A.option("Mesh.MeshSizeMax",0.5)
        A.mesh.generate(2)
        authority=A.CURRENT[].meshing.mesh_authority
        log=M._model_mesh_log(M._model_mesh_authority(authority.mesh,authority.done;
            private_rank=authority.private_rank))
        @test log!==nothing
        first_authority=M._model_mesh_log_authority!(log,authority.done)
        count=size(first_authority.mesh.coords,2)
        rank=maximum(authority.private_rank)
        # A row the full constructors would reject invalidates the log rather
        # than entering it: a non-finite coordinate, then (on a fresh log) a
        # repeated node label.
        @test !M._model_mesh_log_node!(log,(NaN,0.0,0.0),(2,Int32(1)),[0.0,0.0],
            typemax(UInt64),rank+1)
        @test !log.valid
        @test M._model_mesh_log_authority!(log,authority.done)===nothing
        log=M._model_mesh_log(M._model_mesh_authority(authority.mesh,authority.done;
            private_rank=authority.private_rank))
        taken=first_authority.mesh.entity_data.external_node_tags[1]
        @test !M._model_mesh_log_node!(log,(0.25,0.25,0.0),(2,Int32(1)),[0.25,0.25],taken,rank+1)
        log=M._model_mesh_log(M._model_mesh_authority(authority.mesh,authority.done;
            private_rank=authority.private_rank))
        before=M._model_mesh_log_authority!(log,authority.done)
        @test M._model_mesh_log_node!(log,(0.25,0.25,0.0),(2,Int32(1)),[0.25,0.25,9.0],
            typemax(UInt64),rank+1)
        grown=M._model_mesh_log_authority!(log,authority.done)
        @test size(grown.mesh.coords,2)==count+1
        @test grown.mesh.entity_data.node_parametric[end]==[0.25,0.25]
        @test count+1 in grown.products[(2,Int32(1))].rows
        # An authority taken before the append keeps the rows it was taken with.
        @test size(before.mesh.coords,2)==count
        @test length(before.private_rank)==count
        @test before.mesh.coords==grown.mesh.coords[:,1:count]==first_authority.mesh.coords
    finally
        A.finalize()
    end
end

# A plan works on a copy of the model that shares the previous mesh authority
# with the session instead of cloning it. A plan that fails must still leave
# the session's authority exactly as it was.
const QUADTRI_WITHOUT_RECOMBINE="""
lc=0.5;
Point(1)={0,0,0,lc};Point(2)={1,0,0,lc};Point(3)={1,1,0,lc};Point(4)={0,1,0,lc};
Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
Transfinite Curve{1,2,3,4}=3;Transfinite Surface{1};Recombine Surface{1};
out[]=Extrude {0,0,1} { Surface{1}; Layers{2}; QuadTriNoNewVerts RecombLaterals; };
"""

@testset "a failed plan leaves the shared authority untouched" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        mktempdir() do directory
            path=joinpath(directory,"quadtri.geo")
            write(path,QUADTRI_WITHOUT_RECOMBINE)
            A.open_geo!(path;mesh_dim=0)
        end
        A.mesh.generate(2)
        authority=A.CURRENT[].meshing.mesh_authority
        saved=deepcopy(authority)
        nodes=A.mesh.get_nodes();elements=A.mesh.get_elements()
        @test_throws ArgumentError A.mesh.generate(3)
        @test A.CURRENT[].meshing.mesh_authority===authority
        @test M._model_mesh_same(saved,authority)
        @test A.mesh.get_nodes()==nodes
        @test A.mesh.get_elements()==elements
    finally
        A.finalize()
    end
end

# A Surface's constructor input is checked against its actual product. When
# both list the same nodes and cells in the same order the check compares
# them in place; any other valid input takes the general check, and an input
# that changes a cell is rejected by either.
function rebuilt_input(input,columns,labels)
    data=input.entity_data
    block=only(input.blocks)
    metadata=Tessella.Elements._native_mixed_entity_data(data.entities;
        node_entities=data.node_entities,node_parametric=data.node_parametric,
        external_node_tags=data.external_node_tags,
        block_entities=[data.block_entities[1][columns]],
        external_element_tags=[labels])
    return Tessella.MixedMesh(input.coords,
        [Tessella.ElementBlock(block.msh,block.nodes[:,columns],block.tags[columns])];
        physical_names=input.physical_names,entity_data=metadata,
        elementary_entities=[input.elementary_entities[1][columns]])
end

@testset "aligned and general Surface input checks agree" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        plane_grid(2)
        A.option("Mesh.MeshSizeMin",0.4);A.option("Mesh.MeshSizeMax",0.4)
        A.mesh.generate(2)
        authority=A.CURRENT[].meshing.mesh_authority
        actual=M._model_mesh_product(authority,2,1;producer_order=true,first_order=true)
        input=authority.surface_constructor_inputs[(2,Int32(1))]
        labels=only(input.entity_data.external_element_tags)
        cells=length(labels)
        @test cells>=4
        @test M._model_surface_nodes_aligned(actual,input)
        @test M._model_surface_cells_aligned(actual,input)
        @test M._model_surface_constructor_binding(actual,input,"test")===nothing
        # The same cells listed backwards: valid, but not aligned.
        backwards=rebuilt_input(input,collect(cells:-1:1),reverse(labels))
        @test !M._model_surface_cells_aligned(actual,backwards)
        @test M._model_surface_constructor_binding(actual,backwards,"test")===nothing
        # Two cells exchange their labels: the same label set, another mesh.
        exchanged=copy(labels);exchanged[1],exchanged[2]=labels[2],labels[1]
        changed=rebuilt_input(input,collect(1:cells),exchanged)
        @test !M._model_surface_cells_aligned(actual,changed)
        @test_throws ArgumentError M._model_surface_constructor_binding(actual,changed,"test")
        # A repeated label is never accepted.
        repeated=copy(labels);repeated[2]=labels[1]
        @test_throws ArgumentError M._model_surface_constructor_binding(actual,
            rebuilt_input(input,collect(1:cells),repeated),"test")
    finally
        A.finalize()
    end
end

# A producer reads a face's boundary Curves row by row from the authority
# instead of cutting a validated mesh per Curve per face. The rows must be the
# cut's rows, and reading them must not allocate with the size of the Curve.
@testset "product rows are the product cut, read in place" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        plane_grid(2)
        A.option("Mesh.MeshSizeMin",0.25);A.option("Mesh.MeshSizeMax",0.25)
        A.mesh.generate(2)
        authority=A.CURRENT[].meshing.mesh_authority
        @test M._model_mesh_product_rows(authority,1,10_000)===nothing
        for (dim,tag) in ((1,1),(1,7),(2,3),(0,5))
            cut=M._model_mesh_product(authority,dim,tag)
            view=M._model_mesh_product_rows(authority,dim,tag)
            @test (cut===nothing)==(view===nothing)
            cut===nothing && continue
            data=view.mesh.entity_data
            @test length(view.rows)==size(cut.coords,2)
            @test view.mesh.coords[:,view.rows]==cut.coords
            @test data.external_node_tags[view.rows]==cut.entity_data.external_node_tags
            @test data.node_entities[view.rows]==cut.entity_data.node_entities
            @test isequal(data.node_parametric[view.rows],cut.entity_data.node_parametric)
            @test [(block.msh,block.nodes,block.tags) for block in view.blocks]==
                [(block.msh,block.nodes,block.tags) for block in cut.blocks]
            @test view.element_tags==cut.entity_data.external_element_tags
        end
        if Base.JLOptions().compile_enabled==1
            read_rows(authority,tag)=M._model_mesh_product_rows(authority,1,tag)
            read_rows(authority,1)
            @test (@allocated read_rows(authority,1))<=256
        end
    finally
        A.finalize()
    end
end

# A generation plan reconciles the model's entity records with its mesh, and
# the commit of that very mesh skips deriving them a second time. That is
# sound because the second derivation changes nothing: here it is run by hand
# on the committed state and every record must come out as it was.
same_record(a,b)=all(isequal(getfield(a,name),getfield(b,name)) for name in fieldnames(typeof(a)))

@testset "a generation's entity records need no second reconciliation" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        plane_grid(2)
        A.model.add_physical_group(2,[1,4];tag=7)
        A.option("Mesh.MeshSizeMin",0.4);A.option("Mesh.MeshSizeMax",0.4)
        for dimension in (1,2)
            A.mesh.generate(dimension)
            model=A.CURRENT[];mesh=A.LAST_MESH[];class=A.LAST_MESH_CLASS[]
            before=deepcopy(model.meshing.attached)
            @test !isempty(before)
            @test any(record->!isempty(record.element_tags),values(before))
            identity=Dict(tag=>tag for tag in class.public_tags.node_tags)
            A._dim01_reconcile_records!(model,collect(class.public_tags.authority),mesh,identity;
                runtime_parameters=class.node_parametric)
            after=model.meshing.attached
            @test Set(keys(after))==Set(keys(before))
            @test all(same_record(before[key],after[key]) for key in keys(before))
        end
    finally
        A.finalize()
    end
end

# `generate` returns a detached copy of the mesh with its blocks grouped by
# element type. For a mesh with entity metadata and nothing else attached the
# copy is made once per array; it must be the mesh the public constructors
# build from the grouped arrays, also when a type occurs in separated blocks.
function same_content(a,b)
    typeof(a)==typeof(b) || return false
    a isa AbstractDict && return Set(keys(a))==Set(keys(b)) &&
        all(same_content(a[key],b[key]) for key in keys(a))
    a isa AbstractArray && return size(a)==size(b) &&
        all(same_content(x,y) for (x,y) in zip(a,b))
    a isa Tuple && return length(a)==length(b) &&
        all(same_content(x,y) for (x,y) in zip(a,b))
    (isstructtype(typeof(a)) && !(a isa Union{Number,AbstractString,Symbol})) || return isequal(a,b)
    return all(same_content(getfield(a,name),getfield(b,name)) for name in fieldnames(typeof(a)))
end

@testset "the returned mesh copy is the grouped mesh" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        plane_grid(2)
        A.option("Mesh.MeshSizeMin",0.4);A.option("Mesh.MeshSizeMax",0.4)
        A.mesh.generate(2)
        source=A.CURRENT[].meshing.mesh_authority.mesh
        plain=A._copy_mesh(source)
        @test same_content(plain,source)
        @test plain.coords!==source.coords
        @test plain.entity_data.external_node_tags!==source.entity_data.external_node_tags
        # Split the triangles into two blocks around the other blocks.
        data=source.entity_data
        triangles=findfirst(block->block.msh==2,source.blocks)
        block=source.blocks[triangles];half=size(block.nodes,2)÷2
        @test half>=1
        others=[index for index in eachindex(source.blocks) if index!=triangles]
        part(values,range)=values[range]
        first_half=1:half;second_half=half+1:size(block.nodes,2)
        blocks=vcat([Tessella.ElementBlock(2,block.nodes[:,first_half],block.tags[first_half])],
            [source.blocks[index] for index in others],
            [Tessella.ElementBlock(2,block.nodes[:,second_half],block.tags[second_half])])
        split(values)=vcat([part(values[triangles],first_half)],[values[index] for index in others],
            [part(values[triangles],second_half)])
        separated=Tessella.MixedMesh(source.coords,blocks;physical_names=source.physical_names,
            entity_data=Tessella.Elements._native_mixed_entity_data(data.entities;
                node_entities=data.node_entities,node_parametric=data.node_parametric,
                external_node_tags=data.external_node_tags,
                block_entities=split(data.block_entities),
                external_element_tags=split(data.external_element_tags)),
            elementary_entities=split(source.elementary_entities))
        copied=A._canonical_mixed_copy(separated)
        @test [block.msh for block in copied.blocks]==vcat([2],[source.blocks[index].msh for index in others])
        grouped(values)=vcat([values[triangles]],[values[index] for index in others])
        expected=Tessella.MixedMesh(source.coords,
            vcat([block],[source.blocks[index] for index in others]);
            physical_names=source.physical_names,
            entity_data=Tessella.Elements._native_mixed_entity_data(data.entities;
                node_entities=data.node_entities,node_parametric=data.node_parametric,
                external_node_tags=data.external_node_tags,
                block_entities=grouped(data.block_entities),
                external_element_tags=grouped(data.external_element_tags)),
            elementary_entities=grouped(source.elementary_entities))
        @test same_content(copied,expected)
    finally
        A.finalize()
    end
end

staged_records(plan)=plan.model.meshing.attached

# `_dim01_finish!` orders, renumbers and groups cells on columns. A caller
# that holds whole blocks passes them with `_dim01_add_cell_block!` instead
# of one object and one node vector per cell; the finished mesh and every
# label map must be those of the object route. The blocks are given with the
# triangles first, so neither route receives its cells in publication order.
function finished(mesh,model,blockwise,renumber,p2_pass)
    data=mesh.entity_data
    staged=M._model_mesh_transaction_copy(model)
    assembly=A._Dim01Assembly(A._Dim01Node[],A._Dim01Cell[],Dict{UInt64,Int32}(),
        A.NODE_TAG_MAX[],A.ELEMENT_TAG_MAX[],Int(typemax(Int32)),Int(typemax(Int32)),"test")
    for row in axes(mesh.coords,2)
        A._dim01_add_node!(assembly,data.external_node_tags[row],
            ntuple(axis->mesh.coords[axis,row],3),data.node_entities[row],data.node_parametric[row])
    end
    order=sortperm([block.msh for block in mesh.blocks];rev=true)
    columns=blockwise ? A._Dim01CellColumns() : nothing
    for index in order
        block=mesh.blocks[index]
        owners=data.block_entities[index];labels=data.external_element_tags[index]
        if blockwise
            A._dim01_add_cell_block!(columns,assembly,block.msh,block.nodes,owners,labels,block.tags)
        else
            dimension=Int(Tessella.Elements.msh_dimension(block.msh))
            for column in axes(block.nodes,2)
                A._dim01_add_cell!(assembly,block.msh,(dimension,owners[column]),
                    block.nodes[:,column],labels[column];physical=block.tags[column])
            end
        end
    end
    plan=A._dim01_finish!(assembly,staged,Tessella.Model.model_entities(staged),
        Set{Tuple{Int,Int32}}(),falses(length(assembly.nodes)),renumber,p2_pass;cells=columns)
    return plan,assembly
end

@testset "a finish from cell blocks is the finish from cell objects" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        plane_grid(2)
        A.model.add_physical_group(2,[2,3];tag=5)
        A.option("Mesh.MeshSizeMin",0.4);A.option("Mesh.MeshSizeMax",0.4)
        A.mesh.generate(2)
        model=A.CURRENT[];mesh=model.meshing.mesh_authority.mesh
        @test length(mesh.blocks)>=2
        for renumber in (true,false),p2_pass in (false,true)
            objects,assembly=finished(mesh,model,false,renumber,p2_pass)
            blocks,empty_assembly=finished(mesh,model,true,renumber,p2_pass)
            @test isempty(empty_assembly.cells)
            @test same_content(objects.mesh,blocks.mesh)
            @test objects.node_tag_map==blocks.node_tag_map
            @test objects.element_tag_map==blocks.element_tag_map
            @test objects.visibility_tag_map==blocks.visibility_tag_map
            @test (objects.max_node_tag,objects.max_element_tag)==(blocks.max_node_tag,blocks.max_element_tag)
            @test same_content(staged_records(objects),staged_records(blocks))
            # The object route leaves its cells in publication order with
            # their published labels.
            @test issorted(assembly.cells;by=A._dim01_cell_order)
            @test sort!([cell.tag for cell in assembly.cells])==
                sort!(reduce(vcat,objects.mesh.entity_data.external_element_tags))
        end
        # A block is checked as each of its cells would be, and a finish
        # takes objects or columns, not both.
        columns=A._dim01_add_cell_block!(A._Dim01CellColumns(),
            A._Dim01Assembly(A._Dim01Node[],A._Dim01Cell[],Dict{UInt64,Int32}(),UInt64(0),UInt64(0),
                Int(typemax(Int32)),Int(typemax(Int32)),"test"),
            15,reshape(Int32[1],1,1),Int32[1],UInt64[1],Int32[0])
        @test length(columns)==1
        @test_throws ArgumentError A._dim01_add_cell_block!(A._Dim01CellColumns(),
            A._Dim01Assembly(A._Dim01Node[],A._Dim01Cell[],Dict{UInt64,Int32}(),UInt64(0),UInt64(0),
                Int(typemax(Int32)),Int(typemax(Int32)),"test"),
            2,reshape(Int32[1,2],2,1),Int32[1],UInt64[1],Int32[0])
        _,assembly=finished(mesh,model,false,true,false)
        @test !isempty(assembly.cells)
        @test_throws ArgumentError A._dim01_finish!(assembly,M._model_mesh_transaction_copy(model),
            Tessella.Model.model_entities(model),Set{Tuple{Int,Int32}}(),
            falses(length(assembly.nodes)),true,false;cells=columns)
    finally
        A.finalize()
    end
end

# A refinement labels the cells of the refined carrier parent by parent, after
# checking that each parent is the actual cell of its label on the same rows.
# The check used a catalog of the whole previous mesh; it now looks each
# parent up by position. `catalog_cell_tags` is the catalog route, kept here
# as the reference.
function catalog_cell_tags(plan,previous)
    catalog,_=M._model_mesh_primary_catalog(previous.mesh)
    output=Vector{Union{Nothing,UInt64}}[]
    for (bi,(msh,cells)) in enumerate(M._model_mesh_carrier_blocks(previous.carrier))
        children=msh==1 ? 2 : msh==2 ? 4 : 8
        tags=Union{Nothing,UInt64}[]
        for parent in axes(cells,2)
            source=previous.carrier_cells[bi][parent]
            source===nothing && return Vector{Union{Nothing,UInt64}}[]
            actual=get(catalog,source,nothing)
            actual!==nothing && actual[2]==previous.carrier_rows[cells[:,parent]] ||
                return Vector{Union{Nothing,UInt64}}[]
            append!(tags,plan.refine_cell_tags[source])
        end
        push!(output,tags)
    end
    return output
end

function with_carrier(previous,carrier;cells=previous.carrier_cells)
    return M._ModelMeshAuthority(previous.mesh,previous.products,previous.done,
        previous.private_rank,previous.native_inputs,carrier,previous.carrier_rows,
        previous.producer_rows,cells,previous.surface_constructor_inputs,
        previous.emitted_geometry,previous.derived)
end

@testset "refined carrier cell labels match the catalog route" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        plane_grid(2)
        A.option("Mesh.MeshSizeMin",0.4);A.option("Mesh.MeshSizeMax",0.4)
        A.mesh.generate(2)
        previous=A.CURRENT[].meshing.mesh_authority
        carrier=previous.carrier
        @test carrier isa Tessella.Mesh
        @test !isempty(previous.carrier_cells)
        children=Dict{UInt64,Vector{UInt64}}()
        for (bi,(msh,cells)) in enumerate(M._model_mesh_carrier_blocks(carrier)),parent in axes(cells,2)
            label=previous.carrier_cells[bi][parent]
            label===nothing && continue
            children[label]=UInt64[1000*label+child for child in 1:(msh==1 ? 2 : msh==2 ? 4 : 8)]
        end
        plan=(refine_cell_tags=children,)
        expected=catalog_cell_tags(plan,previous)
        @test !isempty(expected) && sum(length,expected)>0
        @test M._model_mesh_refined_cell_tags(plan,previous)==expected
        # A parent listed on the same rows in another order is not the actual
        # cell: both routes give up.
        triangles=copy(carrier.tris);triangles[1,1],triangles[2,1]=triangles[2,1],triangles[1,1]
        permuted=with_carrier(previous,Tessella.Mesh(carrier.coords;segs=carrier.segs,tris=triangles,
            tets=carrier.tets,seg_tag=carrier.seg_tag,tri_tag=carrier.tri_tag,tet_tag=carrier.tet_tag))
        @test isempty(catalog_cell_tags(plan,permuted))
        @test isempty(M._model_mesh_refined_cell_tags(plan,permuted))
        # So do they for a parent whose label names no actual cell.
        unknown=[copy(tags) for tags in previous.carrier_cells]
        block=findfirst(!isempty,unknown);unknown[block][1]=typemax(UInt64)
        @test isempty(M._model_mesh_refined_cell_tags(plan,with_carrier(previous,carrier;cells=unknown)))
        @test isempty(catalog_cell_tags(plan,with_carrier(previous,carrier;cells=unknown)))
    finally
        A.finalize()
    end
end

# The refined carrier of an authority is the uniform refinement of the old
# carrier, traced for its node and cell order. Its caller certifies every
# child against the refined mesh, so the trace runs `trusted`: without
# validating its input or its result. The trace itself must not change.
@testset "a trusted refinement trace is the validated trace" begin
    R=Tessella.Refine
    coords=[0.0 1.0 0.0 0.0 1.0; 0.0 0.0 1.0 0.0 1.0; 0.0 0.0 0.0 1.0 1.0]
    mesh=Tessella.Mesh(coords;segs=Int32[1; 2;;],tris=Int32[1 2; 2 5; 3 3],
        tets=Int32[1 2; 2 5; 3 3; 4 4],seg_tag=Int32[7],tri_tag=Int32[8,9],tet_tag=Int32[10,11])
    checked=R._refine_uniform_traced(mesh;require_positive_tets=false)
    trusted=R._refine_uniform_traced(mesh;require_positive_tets=false,trusted=true)
    for name in (:coords,:segs,:tris,:tets,:seg_tag,:tri_tag,:tet_tag)
        @test getfield(trusted.mesh,name)==getfield(checked.mesh,name)
    end
    @test trusted.primary_map==checked.primary_map
    @test trusted.midpoint_ids==checked.midpoint_ids
    @test size(trusted.mesh.tets,2)==16
    # Only the validating trace refuses a negatively oriented tetrahedron when
    # positive ones are required.
    inverted=Tessella.Mesh(coords;tets=Int32[1; 3; 2; 4;;])
    @test_throws ArgumentError R._refine_uniform_traced(inverted)
    @test size(R._refine_uniform_traced(inverted;trusted=true).mesh.tets,2)==8
end

# Entity records are built and copied through `MixedEntity`'s checks many
# times per generation. Values that already have the stored types (`Int32`
# lists, six `Float64` bounds) take a route without conversions; it must
# accept, reject and return exactly what the general route does.
function int32_outcome(values;kwargs...)
    try
        Tessella.Elements._mixed_int32_vector(values,"test";kwargs...)
    catch err
        err isa ArgumentError || rethrow()
        err.msg
    end
end

@testset "typed entity record values are checked like any others" begin
    E=Tessella.Elements
    long=Int32.(1:40)
    for values in (Int32[],Int32[3,1,2],Int32[1,0,2],Int32[2,-1],Int32[4,4],Int32[typemin(Int32)],
            long,vcat(long,Int32[7])),
            options in ((;),(;positive=true),(;nonzero=true),(;absolute_fits=true),
            (;unique_values=true),(;positive=true,unique_values=true),
            (;nonzero=true,absolute_fits=true))
        typed=int32_outcome(values;options...)
        general=int32_outcome(Int.(values);options...)
        @test typed==general
        typed isa Vector && @test typed!==values && typed isa Vector{Int32}
    end
    bounds=(0.0,1.0,2.0,3.0,4.0,5.0)
    tuple_record=E.MixedEntity(2,7,bounds;physical_tags=Int32[3,1],boundaries=Int32[1,-2],
        embedded_curves=Int32[9,5])
    vector_record=E.MixedEntity(2,7,collect(bounds);physical_tags=[3,1],boundaries=[1,-2],
        embedded_curves=[9,5])
    @test same_content(tuple_record,vector_record)
    @test tuple_record.embedded_curves==Int32[5,9]
    @test_throws ArgumentError E.MixedEntity(2,7,(0.0,1.0,2.0,3.0,Inf,5.0))
    @test_throws ArgumentError E.MixedEntity(2,7,(4.0,1.0,2.0,3.0,4.0,5.0))
    @test_throws ArgumentError E.MixedEntity(2,7,bounds;boundaries=Int32[5],embedded_curves=Int32[5])
    @test_throws ArgumentError E.MixedEntity(2,7,bounds;embedded_curves=Int32[5,5])
    @test_throws ArgumentError E.MixedEntity(1,7,bounds;boundaries=Int32[1,0])
    @test same_content(E.MixedEntity(0,3,(1.0,2.0,3.0)),E.MixedEntity(0,3,(1.0,2.0,3.0,1.0,2.0,3.0)))
end

# `generate(1)` keeps a geometry carrier beside the mesh: the same nodes and
# cells with the public (unsigned) Curve descriptors. It is copied from the
# validated mesh once per array. The copy must be a valid mesh in normal
# form (the public constructors rebuild it unchanged), detached from the
# mesh it was made from.
@testset "the dimension-one geometry carrier is one detached copy" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        plane_grid(2)
        A.option("Mesh.MeshSizeMin",0.4);A.option("Mesh.MeshSizeMax",0.4)
        A.mesh.generate(1)
        authority=A.CURRENT[].meshing.mesh_authority
        mesh=authority.mesh;carrier=authority.carrier
        @test carrier isa Tessella.MixedMesh && carrier!==mesh
        @test Tessella.validate(carrier).ok
        data=carrier.entity_data
        @test carrier.coords!==mesh.coords && data.external_node_tags!==mesh.entity_data.external_node_tags
        @test all(first(pair)!==last(pair) for pair in zip(carrier.blocks,mesh.blocks))
        @test data.external_node_tags==mesh.entity_data.external_node_tags
        @test [block.nodes for block in carrier.blocks]==[block.nodes for block in mesh.blocks]
        curves=[record for record in values(data.entities) if record.dim==1]
        @test !isempty(curves) && all(record->all(>(0),record.boundaries),curves)
        @test any(record->any(<(0),record.boundaries),values(mesh.entity_data.entities))
        rebuilt=Tessella.MixedMesh(carrier.coords,carrier.blocks;physical_names=carrier.physical_names,
            entity_data=Tessella.Elements._native_mixed_entity_data(data.entities;
                node_entities=data.node_entities,node_parametric=data.node_parametric,
                external_node_tags=data.external_node_tags,block_entities=data.block_entities,
                external_element_tags=data.external_element_tags),
            elementary_entities=carrier.elementary_entities)
        @test same_content(rebuilt,carrier)
    finally
        A.finalize()
    end
end

# A caller that only reads an entity's product takes a shared cut: the rows
# sliced into arrays of its own, the entity records, parameter vectors, cell
# labels and names left in the authority's storage, nothing checked again.
# It must read exactly like the detached cut, in every order a cut can have.
@testset "a shared product cut reads like the detached cut" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        plane_grid(2)
        A.option("Mesh.MeshSizeMin",0.3);A.option("Mesh.MeshSizeMax",0.3)
        A.mesh.generate(2)
        authority=A.CURRENT[].meshing.mesh_authority
        source=authority.mesh
        for (dim,tag) in ((0,5),(1,1),(1,7),(2,1),(2,4)),
                flags in ((;),(;producer_order=true),(;producer_order=true,first_order=true))
            detached=M._model_mesh_product(authority,dim,tag;flags...)
            shared=M._model_mesh_product(authority,dim,tag;flags...,shared=true)
            @test (detached===nothing)==(shared===nothing)
            detached===nothing && continue
            @test same_content(detached,shared)
            @test shared.coords!==source.coords
            @test shared.entity_data.external_element_tags===authority.products[(dim,Int32(tag))].element_tags
            @test shared.physical_names===source.physical_names
            @test detached.physical_names!==source.physical_names
            @test all(record->any(other->other===record,values(source.entity_data.entities)),
                values(shared.entity_data.entities))
        end
        @test M._model_mesh_product(authority,1,10_000;shared=true)===nothing
        if Base.JLOptions().compile_enabled==1
            cut(authority,shared)=M._model_mesh_product(authority,1,1;shared)
            cut(authority,true);cut(authority,false)
            @test 2*(@allocated cut(authority,true))<(@allocated cut(authority,false))
        end
        # The input of a merged Surface adopts the shared cut made for it. It
        # must be the input the copying route builds from a detached cut.
        stored=authority.surface_constructor_inputs[(2,Int32(1))]
        data=stored.entity_data
        production=(node_tags=data.external_node_tags,node_owners=data.node_entities,
            node_parameters=data.node_parametric,
            element_tags=reduce(vcat,data.external_element_tags),constructor_mesh=stored)
        copying=M._model_surface_constructor_input(production,
            M._model_mesh_product(authority,2,1;producer_order=true),"test")
        cut=M._model_mesh_product(authority,2,1;producer_order=true,shared=true)
        adopting=M._model_surface_constructor_input(production,cut,"test";adopt=true)
        @test same_content(copying,adopting)
        @test adopting.coords===cut.coords
        @test adopting.entity_data.external_node_tags===cut.entity_data.external_node_tags
        # What an input must own although the cut shares it with the authority.
        @test all(first(pair)===nothing || first(pair)!==last(pair)
            for pair in zip(adopting.entity_data.node_parametric,cut.entity_data.node_parametric))
        @test adopting.physical_names!==source.physical_names
        @test M._model_surface_constructor_binding(
            M._model_mesh_product(authority,2,1;producer_order=true,first_order=true),adopting,"test")===nothing
    finally
        A.finalize()
    end
end

# A transfinite block of hexahedra is a lattice of rows. Its carrier cells are
# matched to lattice cells by position, and its boundary faces and edges are
# read from the lattice, instead of sorting and hashing the corners, faces
# and edges of every cell. Each must agree with the generic table it replaces.
@testset "a lattice block's cells, boundary and edges come from the lattice" begin
    counts=(2,3,2);spanj=counts[2]+1;spank=counts[3]+1
    total=(counts[1]+1)*spanj*spank
    # Lattice point (i,j,k) is row `rows[(i*spanj+j)*spank+k+1]`, in a scrambled numbering.
    rows=Int32[mod(7*position,total)+1 for position in 1:total]
    @test sort(rows)==1:total
    place=zeros(Int32,total)
    for (position,row) in enumerate(rows)
        place[row]=position
    end
    at(i,j,k)=rows[(i*spanj+j)*spank+k+1]
    coords=zeros(3,total)
    for i in 0:counts[1],j in 0:counts[2],k in 0:counts[3]
        coords[:,at(i,j,k)]=[i,j,k]
    end
    lattice_cells=[(i,j,k) for i in 0:counts[1]-1 for j in 0:counts[2]-1 for k in 0:counts[3]-1]
    corners(i,j,k)=Int32[at(i,j,k),at(i+1,j,k),at(i+1,j+1,k),at(i,j+1,k),
        at(i,j,k+1),at(i+1,j,k+1),at(i+1,j+1,k+1),at(i,j+1,k+1)]
    # Carrier cells in another order than the lattice's.
    order=reverse(eachindex(lattice_cells))
    cells=reduce(hcat,[corners(lattice_cells[index]...) for index in order])
    mesh=Tessella.MixedMesh(coords,[Tessella.ElementBlock(5,cells)])
    @test Tessella.validate(mesh).ok
    columns=A._higher_transfinite_hex_columns!(zeros(Int,counts...),cells,rows,place,spanj,spank)
    @test all(columns[(lattice_cells[index].+1)...]==column for (column,index) in enumerate(order))
    # The same cell with its corners listed in another order is still found;
    # a cell that is not a lattice cell is not.
    scrambled=copy(cells);scrambled[:,1]=scrambled[[5,6,7,8,1,2,3,4],1]
    @test A._higher_transfinite_hex_columns!(zeros(Int,counts...),scrambled,rows,place,spanj,spank)==columns
    foreign=copy(cells);foreign[1,1]=foreign[2,1]
    missing_one=A._higher_transfinite_hex_columns!(zeros(Int,counts...),foreign,rows,place,spanj,spank)
    @test count(iszero,missing_one)==1
    @test A._higher_lattice_boundary_faces(rows,counts)==first(M._volume_cell_boundary_faces(mesh))
    edges=M._volume_edge_set(mesh)
    @test all(A._higher_lattice_edge(place,counts,a,b)==(minmax(Int32(a),Int32(b)) in edges)
        for a in 1:total for b in 1:total if a!=b)
    @test !A._higher_lattice_edge(place,counts,0,1) && !A._higher_lattice_edge(place,counts,1,total+1)
end

# The geometry carrier asks, for every traced cell, whether the assembly
# still holds that cell object. With distinct dense labels the answer is read
# from a vector by label; otherwise from a table of the objects. Both must
# name the same cells live.
@testset "live traced cells are found by label or by object" begin
    cell(tag)=A._Dim01Cell(4,(3,Int32(1)),Int32[1,2,3,4],UInt64(tag),nothing,Int32(0))
    cells=[cell(tag) for tag in (3,1,2,5)]
    stale=cell(2)
    by_label=A._higher_live_cells(cells)
    @test by_label isa Vector{Int32}
    block=Union{Nothing,A._Dim01Cell}[cells[2],nothing,cells[4]]
    @test A._higher_live_cell_tags(by_label,cells,block,"test")==[UInt64(1),nothing,UInt64(5)]
    @test_throws ArgumentError A._higher_live_cell_tags(by_label,cells,Union{Nothing,A._Dim01Cell}[stale],"test")
    repeated=[cell(1),cell(1),cell(2)]
    by_object=A._higher_live_cells(repeated)
    @test by_object isa IdDict
    @test A._higher_live_cell_tags(by_object,repeated,Union{Nothing,A._Dim01Cell}[repeated[2]],"test")==[UInt64(1)]
    @test_throws ArgumentError A._higher_live_cell_tags(by_object,repeated,Union{Nothing,A._Dim01Cell}[stale],"test")
    @test A._higher_live_cells([cell(0)]) isa IdDict
    @test A._higher_live_cells([cell(10_000_000)]) isa IdDict
end

const TRANSFINITE_BLOCK="""
Point(1)={0,0,0}; Point(2)={1,0,0}; Point(3)={1,1,0}; Point(4)={0,1,0};
Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
Curve Loop(1)={1,2,3,4}; Plane Surface(1)={1};
Extrude {0,0,1} { Surface{1}; }
Transfinite Curve{:} = 4;
Transfinite Surface{:};
Recombine Surface{:};
Transfinite Volume{1};
"""

@testset "a transfinite hexahedral block meshes through the lattice route" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        mktempdir() do directory
            path=joinpath(directory,"block.geo")
            write(path,TRANSFINITE_BLOCK)
            A.open_geo!(path;mesh_dim=0)
        end
        A.mesh.generate(3)
        types,tags,_=A.mesh.get_elements(3)
        @test types==[5] && length(only(tags))==27
        @test length(A.mesh.get_nodes()[1])==64
        types2,tags2,_=A.mesh.get_elements(2)
        @test types2==[3] && length(only(tags2))==54
    finally
        A.finalize()
    end
end
end
