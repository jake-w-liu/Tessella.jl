module OrdinaryLayersOrphanHistoryChecks
using Tessella,Test,SHA,TOML
const A=Tessella.API
const M=Tessella.Model
const PRIMARY=TOML.parsefile(joinpath(@__DIR__,"..","artifacts","ordinary_layers_orphan_gmsh4152_v1.toml"))
const CURVE="""
    Point(1)={0,0,0,1}; Point(2)={1,0,0,1}; Line(1)={1,2};
    Transfinite Curve{1}=3;
    Extrude {0,0,1} { Curve{1}; Layers{3}; Recombine; }
    """
const SURFACE="""
    Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
    Point(3)={1,1,0,1}; Point(4)={0,1,0,1};
    Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
    Curve Loop(1)={1,2,3,4}; Plane Surface(1)={1};
    Transfinite Curve{:}=2; Transfinite Surface{1}; Recombine Surface{1};
    Extrude {{0,1,0},{0,0,0},Pi/3} { Surface{1}; Layers{3}; Recombine; }
    """
function opened(f,source,dimension)
    A.initialize()
    try
        mktemp() do path,io
            write(io,source);close(io)
            A.open_geo!(path;mesh_dim=dimension==1 ? 0 : dimension)
            dimension==1 && A.mesh.generate(1)
            f()
        end
    finally
        A.finalize()
    end
end
used(dimension,tag)=Set(vcat(A.mesh.get_elements(dimension,tag)[3]...))
function actual_certificate()
    model=A._with_model() do current;deepcopy(current);end
    authority=model.meshing.mesh_authority
    production=M._extrude_volume_pass(model,"ordinary orphan certificate";_mesh_input=authority).productions[1]
    assembly=A._Dim01Assembly(A._Dim01Node[],A._Dim01Cell[],Dict{UInt64,Int32}(),
        UInt64(0),UInt64(0),typemax(Int32),typemax(Int32),"ordinary orphan certificate")
    parameters=M._model_mesh_runtime_parameters(model,authority.mesh)
    for (row,uid) in enumerate(authority.mesh.entity_data.external_node_tags)
        A._dim01_add_node!(assembly,uid,Tuple(authority.mesh.coords[:,row]),
            authority.mesh.entity_data.node_entities[row],parameters[row])
    end
    certificate=A._higher_ordinary_sweep_classification(assembly,model,production,authority,3,1,
        "ordinary orphan certificate")
    @test M._volume_projection_contract(production.mesh,"ordinary orphan certificate")>0
    @test length(M._volume_cell_boundary_faces(production.mesh)[1])==11
    @test allunique(certificate.node_tags[certificate.retained])
    model,authority,production
end
@testset "ordinary Layers preserves owned orphan events and allocator history" begin
    @test PRIMARY["native"]=="4.15.2"
    for name in ("curve","surface","region")
        path=joinpath(@__DIR__,"..","..","validation","native_geometry",PRIMARY[name]["capture"])
        @test bytes2hex(sha256(read(path)))==PRIMARY[name]["capture_sha256"]
    end
    @testset "Curve Source orphan creates only interior lateral events" begin
    opened(CURVE,1) do
        A.mesh.add_nodes(1,1,[999],[.25,0.,0.],[.25])
        A.mesh.generate(2)
        @test length(A.mesh.get_nodes()[1])==PRIMARY["curve"]["nodes"]
        @test A.mesh.get_max_node_tag()==PRIMARY["curve"]["max_node_tag"]
        @test A.mesh.get_max_element_tag()==PRIMARY["curve"]["max_element_tag"]
        @test length(A.mesh.get_elements_by_type(3,5)[1])==6
        ids,coordinates,parameters=A.mesh.get_nodes(1,1)
        orphan=only(ids[findall(==(.25),coordinates[1:3:end])])
        @test !(orphan in used(1,1))
        @test parameters[only(findall(==(orphan),ids))]==.25
        ids,coordinates,parameters=A.mesh.get_nodes(2,5)
        rows=findall(==(.25),coordinates[1:3:end])
        @test length(rows)==2
        @test all(row->!(ids[row] in used(2,5)),rows)
        @test coordinates[3rows]==[1/3,2/3]
        @test parameters[2rows.-1]==[.25,.25]
        @test parameters[2rows]==[1/3,2/3]
        # An unused source row creates neither a copied top-Curve vertex nor
        # a new boundary Line. Actual lower cells still cover the whole strip.
        @test !(.25 in A.mesh.get_nodes(1,2)[2][1:3:end])
        @test length(A.mesh.get_elements_by_type(1,1)[1])==2
        @test length(A.mesh.get_elements_by_type(1,2)[1])==2
    end
    end
    @testset "Surface Source orphan creates only interior Region events" begin
    opened(SURFACE,2) do
        A.mesh.add_nodes(2,1,[999],[.25,.25,0.],[.25,.25])
        A.mesh.generate(3)
        @test length(A.mesh.get_nodes()[1])==PRIMARY["surface"]["nodes"]
        @test A.mesh.get_max_node_tag()==PRIMARY["surface"]["max_node_tag"]
        @test A.mesh.get_max_element_tag()==PRIMARY["surface"]["max_element_tag"]
        source_ids,source_xyz,source_uv=A.mesh.get_nodes(2,1)
        @test source_xyz==[.25,.25,0.] && source_uv==[.25,.25]
        @test isempty(intersect(Set(source_ids),used(2,1)))
        ids,xyz,uv=A.mesh.get_nodes(3,1)
        @test length(ids)==2 && isempty(uv)
        @test isempty(intersect(Set(ids),used(3,1)))
        @test xyz[2:3:end]==[.25,.25]
        @test length(A.mesh.get_elements_by_type(6,1)[1])==3
        @test isempty(A.mesh.get_nodes(2,21)[1])
        model,authority,production=actual_certificate()
        column=only(findall(==(only(source_ids)),production.source_node_tags))
        @test production.column_rows[[1,end],column]==[0,0]
        @test all(>(0),production.column_rows[2:end-1,column])
        @test size(production.mesh.coords,2)==12
        # Regeneration removes the old Region orphans, consumes two new
        # interior source events, and retains the original Source orphan.
        old=copy(ids);clock=A.mesh.get_max_node_tag()
        A.mesh.generate(3)
        @test A.mesh.get_nodes(2,1)==(source_ids,source_xyz,source_uv)
        @test length(A.mesh.get_nodes()[1])==13
        @test A.mesh.get_max_node_tag()==clock+2
        @test isempty(intersect(Set(A.mesh.get_nodes(3,1)[1]),used(3,1)))
        actual_certificate()
    end
    end
    @testset "Region orphan is removed with its consumed UID retained" begin
    opened(SURFACE,3) do
        A.mesh.add_nodes(3,1,[999],[.25,.25,-.05])
        @test 999 in A.mesh.get_nodes(3,1)[1]
        for iteration in 1:3
            A.mesh.generate(3)
            @test isempty(A.mesh.get_nodes(3,1)[1])
            @test length(A.mesh.get_nodes()[1])==10
            @test A.mesh.get_max_node_tag()==999
            @test A.mesh.get_max_element_tag()==PRIMARY["region"]["max_element_tags"][iteration]
            @test length(A.mesh.get_elements_by_type(6,1)[1])==3
            actual_certificate()
        end
    end
    end
end
end
