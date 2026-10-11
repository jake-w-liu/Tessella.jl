module OrdinaryLayersSubdivisionHistoryChecks
using Tessella,Test,TOML,SHA
const A=Tessella.API
const M=Tessella.Model
const PRIMARY=TOML.parsefile(joinpath(@__DIR__,"..","artifacts","ordinary_layers_subdivision_gmsh4152_v1.toml"))
function source(levels;axis="-1")
    """
    Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
    Point(3)={1,1,0,1}; Point(4)={0,1,0,1};
    Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
    Curve Loop(1)={1,2,3,4}; Plane Surface(1)={1};
    Transfinite Curve{:}=2; Transfinite Surface{1};
    Extrude {{0,1,0},{$axis,0,0},Pi/3} { Surface{1}; Layers{$levels}; }
    """
end
function opened(f,text)
    A.initialize()
    try
        mktemp() do path,io
            write(io,text);close(io)
            A.open_geo!(path;mesh_dim=2)
            f()
        end
    finally
        A.finalize()
    end
end
cells(values,width)=Set(Tuple(sort(values[i:i+width-1])) for i in 1:width:length(values))
function signed_volume(a,b,c,d)
    # Independent world-coordinate determinant, evaluated exactly on the
    # binary64 inputs supplied by the actual constructor.
    q=Rational{BigInt}
    u=ntuple(i->q(b[i])-q(a[i]),3);v=ntuple(i->q(c[i])-q(a[i]),3)
    w=ntuple(i->q(d[i])-q(a[i]),3)
    sum(u[i]*(v[mod1(i+1,3)]*w[mod1(i+2,3)]-v[mod1(i+2,3)]*w[mod1(i+1,3)]) for i in 1:3)/6
end
function world_volumes(production,input,levels)
    mesh=production.mesh
    cells=mesh isa Tessella.Mesh ? mesh.tets : only(block.nodes for block in mesh.blocks if block.msh==4)
    tetrahedra=[signed_volume((Tuple(mesh.coords[:,cell[i]]) for i in 1:4)...) for cell in eachcol(cells)]
    source=M._model_mesh_product(input,2,1;caller="ordinary world-coordinate source volume")
    columns=Dict(uid=>column for (column,uid) in enumerate(production.source_node_tags))
    prisms=Rational{BigInt}[]
    for block in source.blocks,cell in eachcol(block.nodes),level in 1:levels
        @test block.msh==2
        rows=Int[columns[source.entity_data.external_node_tags[row]] for row in cell]
        points=ntuple(k->Tuple(mesh.coords[:,production.column_rows[k<=3 ? level : level+1,rows[mod1(k,3)]]]),6)
        volume=sum(signed_volume((points[i] for i in slots)...) for slots in ((1,2,3,6),(1,2,6,5),(1,5,6,4)))
        push!(prisms,abs(volume))
    end
    tetrahedra,prisms
end
function certify(primary)
    model=A._with_model() do current;deepcopy(current);end
    input=model.meshing.mesh_authority
    pass=M._extrude_volume_pass(model,"ordinary late-subdivision certificate";_mesh_input=input)
    production=pass.productions[1]
    assembly=A._Dim01Assembly(A._Dim01Node[],A._Dim01Cell[],Dict{UInt64,Int32}(),
        UInt64(0),UInt64(0),typemax(Int32),typemax(Int32),"ordinary late-subdivision certificate")
    parameters=M._model_mesh_runtime_parameters(model,input.mesh)
    for (row,uid) in enumerate(input.mesh.entity_data.external_node_tags)
        A._dim01_add_node!(assembly,uid,Tuple(input.mesh.coords[:,row]),
            input.mesh.entity_data.node_entities[row],parameters[row])
    end
    result=A._higher_ordinary_sweep_classification(assembly,model,production,input,3,1,
        "ordinary late-subdivision certificate")
    @test validate(production.mesh).ok
    @test M._volume_projection_contract(production.mesh,"ordinary late-subdivision certificate")==primary["tetrahedra"]
    tetrahedra,prisms=world_volumes(production,input,primary["levels"])
    @test all(>(0),tetrahedra) && all(>(0),prisms)
    @test Float64(sum(tetrahedra))≈primary["volume"] atol=2e-11
    @test Float64(sum(prisms))≈primary["volume"] atol=2e-11
    @test sum(tetrahedra)≈sum(prisms) atol=2e-11
    @test length(M._volume_cell_boundary_faces(production.mesh)[1])==primary["outer_faces"]
    @test allunique(result.node_tags[result.retained])
    @test pass.temporary_elements[1]==2primary["levels"]
    @test length(production.cell_events)==primary["tetrahedra"]
    model,input,production,assembly
end
@testset "ordinary late subdivision certifies actual remeshed incidence and native cell history" begin
    @test PRIMARY["native"]=="4.15.2"
    for key in ("capture","science","orphan_capture")
        path=joinpath(@__DIR__,"..","..","validation","native_geometry",PRIMARY[key])
        @test bytes2hex(sha256(read(path)))==PRIMARY[key*"_sha256"]
    end
    for (dimension,tag,point,parameters) in ((1,1,[.25,0.,0.],[.25]),
                                          (2,13,[.75,0.,-.2],[.75,.2]))
        opened(source(3)) do
            A.mesh.add_nodes(dimension,tag,[999],point,parameters)
            A.mesh.generate(3)
            ids,xyz,uv=A.mesh.get_nodes(dimension,tag)
            @test xyz==point && uv==parameters
            @test isempty(intersect(Set(ids),Set(vcat(A.mesh.get_elements(dimension,tag)[3]...))))
            @test length(A.mesh.get_nodes()[1])==19
            @test A.mesh.get_max_node_tag()==999
            @test A.mesh.get_max_element_tag()==154
            @test length(A.mesh.get_elements_by_type(4,1)[1])==18
            certify(PRIMARY["layers_3"])
            # Constrained remeshing preserves the actual Surface vertex vector
            # and never creates new children of an unused Source-Curve row.
            if dimension==1
                @test isempty(A.mesh.get_nodes(2,13)[1])
            end
        end
    end
    for levels in (1,3)
        primary=PRIMARY["layers_$levels"]
        opened(source(levels)) do
            for iteration in 1:3
                A.mesh.generate(3)
                @test length(A.mesh.get_nodes()[1])==primary["nodes"]
                @test A.mesh.get_max_node_tag()==primary["nodes"]
                @test A.mesh.get_max_element_tag()==primary["max_element_tags"][iteration]
                tags,connectivity=A.mesh.get_elements_by_type(4,1)
                @test length(tags)==primary["tetrahedra"]
                @test cells(connectivity,4)==cells(primary["tetrahedron_nodes"],4)
                certify(primary)
            end
        end
    end
    # The coarse source on the rotation-fixed axis also fails in the native
    # SDK. Admission does not turn a failed subdivision into fabricated cells.
    opened(source(3;axis="0")) do
        before=(nodes=deepcopy(A.mesh.get_nodes()),cells=deepcopy(A.mesh.get_elements()),
            node_clock=A.mesh.get_max_node_tag(),element_clock=A.mesh.get_max_element_tag(),
            model=A._with_model() do current;current;end)
        failure=nothing
        try
            A.mesh.generate(3)
        catch error
            failure=error
        end
        # Native reports "Unable to subdivide extruded mesh". The tagged
        # Julia path rejects earlier at its strict actual-shell certificate;
        # the native spelling and exception class are preserved separately.
        @test failure isa ArgumentError
        @test occursin("lower Surface cells do not exactly cover its emitted Volume boundary",sprint(showerror,failure))
        @test A.mesh.get_nodes()==before.nodes
        @test A.mesh.get_elements()==before.cells
        @test A.mesh.get_max_node_tag()==before.node_clock
        @test A.mesh.get_max_element_tag()==before.element_clock
        @test A._with_model(current->current===before.model)
    end
end
end
