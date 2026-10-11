module OrdinaryLayersActualSourceChecks
using Tessella,Test,TOML,SHA
const A=Tessella.API
const M=Tessella.Model
const PRIMARY=TOML.parsefile(joinpath(@__DIR__,"..","artifacts","ordinary_layers_gmsh4152_v1.toml"))
const SQUARE="""
    Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
    Point(3)={1,1,0,1}; Point(4)={0,1,0,1};
    Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
    Curve Loop(1)={1,2,3,4}; Plane Surface(1)={1};
    """
function fixture(count,angle="Pi/3")
    mktemp() do path,io
        write(io,SQUARE*"""
            Transfinite Curve{:}=$count;
            Transfinite Surface{1}; Recombine Surface{1};
            Extrude {{0,1,0},{0,0,0},$angle} {
                Surface{1}; Layers{3}; Recombine;
            }
            """)
        close(io)
        execute_geo(path;mesh_dim=2)
    end
end
function assembled(model,input)
    assembly=A._Dim01Assembly(A._Dim01Node[],A._Dim01Cell[],Dict{UInt64,Int32}(),
        UInt64(0),UInt64(0),typemax(Int32),typemax(Int32),"ordinary Layers actual Source")
    parameters=M._model_mesh_runtime_parameters(model,input.mesh)
    data=input.mesh.entity_data
    for (row,uid) in enumerate(data.external_node_tags)
        A._dim01_add_node!(assembly,uid,Tuple(input.mesh.coords[:,row]),data.node_entities[row],parameters[row])
    end
    assembly
end
function classify(model,input,production)
    A._higher_ordinary_sweep_classification(assembled(model,input),model,production,input,3,1,
        "ordinary Layers actual Source")
end
function certify(model,input,production;faces)
    product=production.mesh
    result=classify(model,input,production)
    @test validate(product).ok
    @test M._volume_projection_contract(product,"ordinary actual Source certificate")>0
    @test length(M._volume_cell_boundary_faces(product)[1])==faces
    @test all(row->!result.retained[row] || result.node_tags[row]>0,eachindex(result.retained))
    @test allunique(result.node_tags[result.retained])
    # Every retained UID carries its actual raw owner and parameter, with no
    # CAD projection of the ruled side or the edited Source.
    data=input.mesh.entity_data
    rows=Dict(uid=>row for (row,uid) in enumerate(data.external_node_tags))
    runtime=M._model_mesh_runtime_parameters(model,input.mesh)
    for row in findall(result.retained)
        source=rows[result.node_tags[row]]
        @test result.class.node_entities[row]==data.node_entities[source]
        @test isequal(result.class.node_parametric[row],runtime[source])
        @test isequal(product.coords[:,row],input.mesh.coords[:,source])
    end
    result
end
@testset "ordinary curved Layers consumes the actual Source and lower UID events" begin
    @test PRIMARY["native"]=="4.15.2"
    for (key,name) in (("fixed_axis","ordinary_layers_fixed_axis_gmsh4152_v1.json"),
                       ("edited_source","ordinary_layers_edited_source_gmsh4152_v1.json"))
        path=joinpath(@__DIR__,"..","..","validation","native_geometry",name)
        @test bytes2hex(sha256(read(path)))==PRIMARY[key]["capture_sha256"]
    end
    for angle in ("Pi/3","-Pi/3")
        model=fixture(2,angle).model
        input=model.meshing.mesh_authority
        pass=M._extrude_volume_pass(model,"ordinary fixed-axis";_mesh_input=input)
        production=pass.productions[1]
        @test size(production.mesh.coords,2)==PRIMARY["fixed_axis"]["nodes"]
        @test Dict(block.msh=>size(block.nodes,2) for block in production.mesh.blocks)==
            Dict(6=>PRIMARY["fixed_axis"]["prisms"])
        certify(model,input,production;faces=PRIMARY["fixed_axis"]["boundary_faces"])
    end

    model=fixture(3).model
    input=model.meshing.mesh_authority
    data=input.mesh.entity_data
    source=only(findall(==((2,Int32(1))),data.node_entities))
    cap=only(findall(==((2,Int32(21))),data.node_entities))
    input.mesh.coords[:,source].=(.5,.5,.025)
    params=M._extrude_gate(model,3,1)
    spec=model.meshing.extrude_specs[(3,1)]
    us=M._extrude_level_us(params)
    transforms=M._extrude_level_transforms(spec,us,"edited actual Source")
    input.mesh.coords[:,cap].=M._extrude_level_at(spec,(.5,.5,.025),us[end],transforms,length(us))
    pass=M._extrude_volume_pass(model,"edited actual Source";_mesh_input=input)
    production=pass.productions[1]
    @test size(production.mesh.coords,2)==PRIMARY["edited_source"]["nodes"]
    @test Dict(block.msh=>size(block.nodes,2) for block in production.mesh.blocks)==
        Dict(5=>PRIMARY["edited_source"]["hexahedra"],6=>PRIMARY["edited_source"]["prisms"])
    result=certify(model,input,production;faces=PRIMARY["edited_source"]["boundary_faces"])
    column=only(findall(==(data.external_node_tags[source]),production.source_node_tags))
    @test production.mesh.coords[:,production.column_rows[1,column]]==[.5,.5,.025]
    for level in eachindex(us)
        @test Tuple(production.mesh.coords[:,production.column_rows[level,column]])==
            (level==1 ? (.5,.5,.025) : M._extrude_level_at(spec,(.5,.5,.025),us[level],transforms,level))
    end

    # Renderer receipts, actual UIDs and literal lower incidence are mandatory.
    # Every failure leaves the detached assembly and original production intact.
    assembly=assembled(model,input)
    state=([(node.tag,node.point,node.owner,node.parameter===nothing ? nothing : copy(node.parameter)) for node in assembly.nodes],
        length(assembly.cells),copy(assembly.node_positions),assembly.node_max,assembly.cell_max)
    check(p,authority=input)=A._higher_ordinary_sweep_classification(assembly,model,p,authority,3,1,"ordinary rejection")
    duplicated=deepcopy(production);duplicated.source_node_tags[end]=duplicated.source_node_tags[1]
    @test_throws ArgumentError check(duplicated)
    missing=deepcopy(production);missing.source_node_tags[end]=typemax(UInt64)
    @test_throws ArgumentError check(missing)
    changed=deepcopy(production);changed.mesh.coords[1,changed.column_rows[1,column]]=nextfloat(.5)
    @test_throws ArgumentError check(changed)
    foreign=deepcopy(input);foreign.mesh.entity_data.node_entities[source]=(3,Int32(999))
    @test_throws ArgumentError check(production,foreign)
    unsupported=deepcopy(input)
    owner=(1,Int32(1));lower=unsupported.products[owner]
    unsupported.products[owner]=M._ModelMeshProduct(owner,lower.rows,
        Tessella.Elements.ElementBlock[],Vector{UInt64}[])
    @test_throws ArgumentError check(production,unsupported)
    repeated=deepcopy(input);lower=repeated.products[(1,Int32(6))]
    repeated.products[(1,Int32(6))]=M._ModelMeshProduct((1,Int32(6)),lower.rows,
        vcat(lower.blocks,lower.blocks),vcat(lower.element_tags,lower.element_tags))
    @test_throws ArgumentError check(production,repeated)
    overlap=deepcopy(input)
    first=overlap.products[(2,Int32(1))];second=overlap.products[(2,Int32(21))]
    overlap.products[(2,Int32(21))]=M._ModelMeshProduct((2,Int32(21)),first.rows,first.blocks,first.element_tags)
    @test_throws ArgumentError check(production,overlap)
    @test [(node.tag,node.point,node.owner,node.parameter) for node in assembly.nodes]==state[1]
    @test length(assembly.cells)==state[2] && assembly.node_positions==state[3]
    @test (assembly.node_max,assembly.cell_max)==state[4:5]
end
end
