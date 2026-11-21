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
end
