using Test
using Tessella
using Tessella.MeshTypes: mesh_crc, nnodes, node, ntets, tet_volume
using Tessella.Elements: mixed_crc

const _GEO_GEOMETRY_EXPRESSION_FIXTURE=normpath(joinpath(
    @__DIR__,"..","fixtures","geo_geometry_expressions.geo"))

function _execute_geometry_expression_source(source::AbstractString;mesh_dim=0)
    return mktemp() do path,io
        write(io,source)
        close(io)
        execute_geo(path;mesh_dim=mesh_dim)
    end
end

function _geometry_expression_error(source::AbstractString)
    try
        _execute_geometry_expression_source(source)
        return nothing
    catch err
        return err
    end
end

@testset "bounded .geo geometry expressions" begin
    parsed=execute_geo(_GEO_GEOMETRY_EXPRESSION_FIXTURE)
    model=parsed.model
    @test sort!(collect(keys(model.points)))==[1,2,3,4,5,6,7,8,100]
    @test sort!(collect(keys(model.curves)))==collect(10:21)
    @test sort!(collect(keys(model.loops)))==collect(30:35)
    @test sort!(collect(keys(model.surfaces)))==collect(40:45)
    @test sort!(collect(keys(model.surface_loops)))==[50]
    @test sort!(collect(keys(model.volumes)))==[60]
    @test model.points[1]==(0.0,0.0,0.0)
    @test model.points[3]==(1.0,1.0,0.0)
    @test model.points[100]==(0.5,0.5,0.5)
    @test all(==(0.45),values(model.point_size))
    @test model.loops[30]==collect(10:13)
    @test model.loops[32]==[10,19,-14,-18]
    @test model.surface_loops[50]==collect(40:45)
    @test model.volumes[60]==[50]
    @test get(model.embeds,(3,60),NTuple{2,Int}[])==[(0,100)]
    @test model.physical==Dict(
        (0,70)=>collect(1:8),(0,71)=>[100],(1,72)=>collect(10:21),
        (2,73)=>collect(40:45),(3,74)=>[60])
    @test model.physical_names==Dict(
        (0,70)=>"corners",(0,71)=>"probe",(1,72)=>"edges",
        (2,73)=>"boundary",(3,74)=>"domain")

    meshed=execute_geo(_GEO_GEOMETRY_EXPRESSION_FIXTURE;mesh_dim=3)
    @test validate(meshed.mesh).ok
    @test nnodes(meshed.mesh)==125
    @test ntets(meshed.mesh)==384
    @test mesh_crc(meshed.mesh).sha==
          "e59cb4912df62d2927553d4cd5cffbdba2e5d6616db319558ed7a0c24d031cc3"
    volume=sum(tet_volume(
        node(meshed.mesh,meshed.mesh.tets[1,cell]),
        node(meshed.mesh,meshed.mesh.tets[2,cell]),
        node(meshed.mesh,meshed.mesh.tets[3,cell]),
        node(meshed.mesh,meshed.mesh.tets[4,cell]))
        for cell in 1:ntets(meshed.mesh))
    @test volume≈1.0 atol=1e-12
    projected=model_to_mixed(meshed.model,meshed.mesh,3,60)
    @test validate(projected).ok
    @test projected.physical_names==model.physical_names
    @test mixed_crc(projected).sha==
          "fb81cf837080eb18dbdffcd092f3e4efb4c344f5ec16c0287d6356899064f3bd"

    optional_size=_execute_geometry_expression_source(
        "Point(1 + 0.9) = {0, 1, 2};")
    @test optional_size.model.points[1]==(0.0,1.0,2.0)
    # `VExpr_Single`'s 4th component defaults to 0 — the "unset" mesh size.
    @test optional_size.model.point_size[1]==0.0

    transformed=_execute_geometry_expression_source(raw"""
        SetFactory("OpenCASCADE");
        tag = 1.9;
        Box(Max(1, tag)) = {0, 0, 0, Sqrt(1), 2, 3};
        Translate {1 / 2, -1 / 2, 1} { Volume{tag}; };
        Dilate {{0, 0, 0}, 2 / 1} { Volume{tag}; };
        Rotate {{0, 0, 2}, {0, 0, 0}, Pi / 2} { Volume{tag}; };
        """)
    @test collect(transformed.model.box_extents[1])≈[-3.0,1.0,2.0,4.0,2.0,6.0]

    primitives=_execute_geometry_expression_source(raw"""
        SetFactory("OpenCASCADE");
        base = 1.9;
        Cylinder(base) = {0, 0, 0, 0, 0, 2, Sqrt(1) / 2};
        Sphere(base + 1) = {1 / 2, 1 / 2, 1 / 2, 1 / 4};
        Cone(base + 2) = {0, 0, 0, 0, 2, 0, 1 / 2, 1 / 4};
        """)
    @test primitives.model.cylinders[1].radius==0.5
    @test primitives.model.spheres[2].center==(0.5,0.5,0.5)
    @test primitives.model.cones[3].axis==(0.0,2.0,0.0)

    boolean=_execute_geometry_expression_source(raw"""
        SetFactory("OpenCASCADE");
        first = 1.9;
        Box(first) = {0, 0, 0, 2, 1, 1};
        Box(first + 1) = {0, 0, 0, 1, 1, 1};
        BooleanDifference(first + 2) =
          { Volume{first}; Delete; }{ Volume{first + 1}; Delete; };
        """)
    @test sort!(collect(keys(boolean.model.volumes)))==[3]
    @test boolean.model.booleans[3]==(op=:difference,operands=[1,2])

    # `Point(0.9)` truncates to `Point(0)` under Gmsh's `(int)` tag cast — a
    # legal entity, not an error.
    @test _execute_geometry_expression_source(
        "Point(0.9) = {0, 0, 0, 1};").model.points[0]==(0.0,0.0,0.0)
    invalid_sources=(
        "Point(missing) = {0, 0, 0, 1};",
        "Point(1) = {0, 0};",
        "Point(1) = {0, 0, 0, 1, 2, 3};",
        "Point(1) = {0, 0, 0, 1}; Line(1) = {1:65537};",
        "Point(1) = {0, 0, 0, 1}; Line(1) = {1, 1, 1};",
        "Curve Loop(1) = {0, 1, 2};",
        "Box(1) = {0, 0, 0, 1, 1};",
        "SetFactory(\"OpenCASCADE\"); Box(1) = {0, 0, 0, 1, 1, 1}; " *
            "Translate {1, 2} { Volume{1}; };",
        "SetFactory(\"OpenCASCADE\"); Box(1) = {0, 0, 0, 1, 1, 1}; " *
            "Box(2) = {2, 0, 0, 1, 1, 1}; " *
            "BooleanUnion(3) = {Volume{1, 2};}{Volume{2};};",
        "Physical Point(1) = {missing};",
    )
    for source in invalid_sources
        @test _geometry_expression_error(source) isa ArgumentError
    end
    point_count_error=_geometry_expression_error("Point(1) = {0, 0};")
    # Upstream's `VExpr_Single` needs 3–5 components — `{a,b}` is a syntax
    # error at the closing `}`.
    @test occursin("syntax error",sprint(showerror,point_count_error))
    # multi-operand lists are legal; a tag appearing in both operand groups
    # is the rejected case
    boolean_count_error=_geometry_expression_error(
        "SetFactory(\"OpenCASCADE\"); Box(1)={0,0,0,1,1,1}; Box(2)={2,0,0,1,1,1}; " *
        "BooleanUnion(3)={Volume{1,2};}{Volume{2};};")
    @test occursin("Boolean operands must be distinct volumes",
                   sprint(showerror,boolean_count_error))
    range_error=_geometry_expression_error(
        "Point(1)={0,0,0,1}; Line(1)={1:65537};")
    @test occursin("expanded list exceeds 65536 entries",
                   sprint(showerror,range_error))

    @test isempty(Docs.undocumented_names(Tessella.GeoExec;private=false))
    @test isempty(Test.detect_ambiguities(Tessella.GeoExec;recursive=true))
end

@testset "bounded .geo Mesh.TransfiniteTri option" begin
    unset=_execute_geometry_expression_source("Point(1) = {0,0,0};")
    @test unset.model.meshing.transfinite_tri==0
    @test unset.transfinite_tri===nothing
    on=_execute_geometry_expression_source(
        "Mesh.TransfiniteTri = 1; Point(1) = {0,0,0};")
    @test on.model.meshing.transfinite_tri==1
    @test on.transfinite_tri==1
    expr=_execute_geometry_expression_source(
        "k = 0; Mesh.TransfiniteTri = k + 0;")
    @test expr.model.meshing.transfinite_tri==0
    @test expr.transfinite_tri==0
    for source in ("Mesh.TransfiniteTri = 2;",
                   "Mesh.TransfiniteTri = 0.5;",
                   "Mesh.TransfiniteTri = -1;")
        @test _geometry_expression_error(source) isa ArgumentError
    end
end
