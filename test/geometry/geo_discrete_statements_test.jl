using Test
using Tessella

function _execute_geo_source(source::AbstractString;mesh_dim=0)
    return mktemp() do path,io
        write(io,source)
        close(io)
        execute_geo(path;mesh_dim=mesh_dim)
    end
end

function _execute_geo_error(source::AbstractString;mesh_dim=0)
    return try
        _execute_geo_source(source;mesh_dim=mesh_dim)
        nothing
    catch err
        err isa InterruptException && rethrow()
        err
    end
end

const _GEO_SQUARE = """
    Point(1) = {0, 0, 0, 0.4};
    Point(2) = {1, 0, 0, 0.4};
    Point(3) = {1, 1, 0, 0.4};
    Point(4) = {0, 1, 0, 0.4};
    Line(1) = {1, 2};
    Line(2) = {2, 3};
    Line(3) = {3, 4};
    Line(4) = {4, 1};
    Curve Loop(1) = {1, 2, 3, 4};
    Plane Surface(1) = {1};
    Physical Surface(7) = {1};
    """

const _GEO_ANNULUS = """
    Point(1) = {0, 0, 0, 0.5};
    Point(2) = {4, 0, 0, 0.5};
    Point(3) = {4, 4, 0, 0.5};
    Point(4) = {0, 4, 0, 0.5};
    Point(5) = {1.5, 1.5, 0, 0.5};
    Point(6) = {2.5, 1.5, 0, 0.5};
    Point(7) = {2.5, 2.5, 0, 0.5};
    Point(8) = {1.5, 2.5, 0, 0.5};
    Line(1) = {1, 2};
    Line(2) = {2, 3};
    Line(3) = {3, 4};
    Line(4) = {4, 1};
    Line(5) = {5, 6};
    Line(6) = {6, 7};
    Line(7) = {7, 8};
    Line(8) = {8, 5};
    Curve Loop(1) = {1, 2, 3, 4};
    Curve Loop(2) = {5, 6, 7, 8};
    Plane Surface(1) = {1, 2};
    Physical Surface(7) = {1};
    """

const _MSH_TWO_SURFACES = """
    \$MeshFormat
    2.2 0 8
    \$EndMeshFormat
    \$Nodes
    5
    1 0 0 0
    2 1 0 0
    3 1 1 0
    4 0 1 0
    5 2 0 0
    \$EndNodes
    \$Elements
    3
    1 2 2 9 5 1 2 3
    2 2 2 9 5 1 3 4
    3 2 2 8 7 2 5 3
    \$EndElements
    """

function _physical_names_of(m)
    out=Dict{Tuple{Int,Int},String}()
    for (key,name) in m.physical_names
        out[key]=name
    end
    return out
end

_ntris_nonempty(ex)=ex.mesh!==nothing && size(ex.mesh.tris,2)>0

@testset ".geo Homology/Cohomology/Betti request parsing" begin
    function requests(src)
        ex=_execute_geo_source(_GEO_SQUARE*src)
        return ex.model.meshing.homology_requests
    end
    # Bare forms queue the whole-model default request.
    for (word,kind) in (("Homology","Homology"),("Cohomology","Cohomology"),
                        ("Betti","Betti"))
        r=requests("$word;")
        @test length(r)==1
        @test r[1].kind==kind
        @test isempty(r[1].domain) && isempty(r[1].subdomain)
        @test r[1].dims==Int[0,1,2,3]
    end
    r=requests("Homology{7};")
    @test r[1].domain==[7] && isempty(r[1].subdomain)
    r=requests("Homology{7,8};")
    @test r[1].domain==[7,8] && isempty(r[1].subdomain)
    # `'{' ListOfDouble ',' ListOfDouble '}'` — each braced group is a list.
    r=requests("Homology{{7},{8}};")
    @test r[1].domain==[7] && r[1].subdomain==[8]
    r=requests("Homology{{7,9},{8,10}};")
    @test r[1].domain==[7,9] && r[1].subdomain==[8,10]
    r=requests("Cohomology{{7},{8}};")
    @test r[1].kind=="Cohomology" && r[1].subdomain==[8]
    r=requests("Betti{{7},{8}};")
    @test r[1].kind=="Betti"
    # Explicit dimensions require both lists.
    r=requests("Homology(0,1){{7},{8}};")
    @test r[1].dims==Int[0,1] && r[1].subdomain==[8]
    r=requests("Homology(2){{7},{}};")
    @test r[1].dims==Int[2]
    # Flat-list grammar: these are syntax errors upstream too.
    for src in ("Homology{7},{8};","Homology{7,{8}};","Homology{{7},{8},{9}};",
                "Homology(0,1);","Homology(0,1){7};","Homology foo;",
                "Homology = 5;","Homology{7}{8};")
        err=_execute_geo_error(_GEO_SQUARE*src)
        @test err isa ArgumentError
        @test occursin("syntax error",sprint(showerror,err))
    end
    # Negative list values are accepted by the grammar but rejected as tags.
    err=_execute_geo_error(_GEO_SQUARE*"Homology{-{7},{8}};")
    @test err isa ArgumentError
    @test occursin("non-negative",sprint(showerror,err))
    # Multiple requests accumulate in order.
    r=requests("Homology{7};\nCohomology{7};\nBetti;")
    @test [q.kind for q in r]==["Homology","Cohomology","Betti"]
end

@testset ".geo discrete-model statement parsing" begin
    # Bare and braced forms parse; runtime errors are fine (no discrete
    # entities exist yet), syntax errors are not.
    for src in ("CreateTopology;","CreateTopology{1,0};",
                "ClassifySurfaces{40*Pi/180,1,0};",
                "ClassifySurfaces{40*Pi/180,1,0,1};",
                "CreateGeometry;","CreateGeometry{Curve{1};};")
        err=_execute_geo_error(_GEO_SQUARE*src)
        @test err===nothing ||
            !occursin("syntax error",sprint(showerror,err))
    end
    for src in ("CreateTopology{1};","CreateTopology{1,0,0};",
                "ClassifySurfaces{1,0};","ClassifySurfaces{1,0,0,0,0};",
                "CreateGeometry{Curve{1};Surface{1}};",
                "CreateTopology foo;")
        err=_execute_geo_error(_GEO_SQUARE*src)
        @test err isa ArgumentError
    end
    # `CreateTopology{1}` is a syntax error upstream (arity 2 required).
    err=_execute_geo_error(_GEO_SQUARE*"CreateTopology{1};")
    @test occursin("syntax error",sprint(showerror,err))
end

@testset ".geo Homology generator lifecycle" begin
    # H_0 = 1 component, H_1 = 1 hole loop for the annulus.
    ex=_execute_geo_source(_GEO_ANNULUS*"Homology{7};\nMesh 2;\n")
    m=ex.model
    names=_physical_names_of(m)
    @test count(==( "H_0{7}1"),values(names))==1
    @test count(==( "H_1{7}1"),values(names))==1
    h0=findfirst(==( "H_0{7}1"),names)
    h1=findfirst(==( "H_1{7}1"),names)
    @test h0[1]==0 && h1[1]==1
    @test length(m.discrete)==2
    @test length(m.physical[h0])==1 && length(m.physical[h1])==1
    @test (0,m.physical[h0][1]) in keys(m.discrete)
    @test (1,m.physical[h1][1]) in keys(m.discrete)
    # Cohomology stores H^ names.
    ex=_execute_geo_source(_GEO_ANNULUS*"Cohomology{7};\nMesh 2;\n")
    names=_physical_names_of(ex.model)
    @test "H^0{7}1" in values(names) && "H^1{7}1" in values(names)
    # Betti computes ranks only — no entities or groups are added.
    ex=_execute_geo_source(_GEO_ANNULUS*"Betti{7};\nMesh 2;\n")
    @test isempty(ex.model.discrete)
    names=_physical_names_of(ex.model)
    @test !any(startswith(n,"H_")||startswith(n,"H^") for n in values(names))
    # Empty domain reports {0} and defaults to top-dimensional entities.
    ex=_execute_geo_source(_GEO_ANNULUS*"Homology;\nMesh 2;\n")
    names=_physical_names_of(ex.model)
    @test "H_0{0}1" in values(names) && "H_1{0}1" in values(names)
    # A repeated request stores nothing extra (upstream isHomologyComputed).
    ex=_execute_geo_source(_GEO_ANNULUS*"Homology{7};\nHomology{7};\nMesh 2;\n")
    names=_physical_names_of(ex.model)
    @test count(n->startswith(n,"H_"),values(names))==2
    @test length(ex.model.discrete)==2
    # A partially computed re-request re-stores chains for every requested
    # dimension — upstream recomputes all of `dim` when any dim is missing,
    # and the re-stored generator repeats the "H_0{7}1" name (pinned Gmsh
    # writes `H_0{7}1` at two physical tags plus `H_1{7}1`).
    ex=_execute_geo_source(_GEO_ANNULUS*
        "Homology(0){{7},{}};\nHomology{{7},{}};\nMesh 2;\n")
    names=_physical_names_of(ex.model)
    @test count(==("H_0{7}1"),values(names))==2
    @test "H_1{7}1" in values(names)
    @test length(ex.model.discrete)==3
    # Explicit dims and subdomain lists work.
    ex=_execute_geo_source(_GEO_ANNULUS*"Homology(0,1){{7},{}};\nMesh 2;\n")
    names=_physical_names_of(ex.model)
    @test "H_0{7}1" in values(names) && "H_1{7}1" in values(names)
    # Requests queued after `Mesh` never execute (upstream computes inside
    # `GModel::mesh`), and the generated mesh itself is preserved.
    ex=_execute_geo_source(_GEO_ANNULUS*"Mesh 2;\nHomology{7};\n")
    @test isempty(ex.model.discrete)
    @test _ntris_nonempty(ex)
    # A domain resolving to no entities is upstream's "Domain is empty".
    err=_execute_geo_error(_GEO_ANNULUS*"Homology{99};\nMesh 2;\n")
    @test err isa ArgumentError
    @test occursin("domain is empty",sprint(showerror,err))
end

@testset ".geo Merge imports .msh as discrete entities" begin
    msh_path=mktempdir()*"/two_surfaces.msh"
    write(msh_path,_MSH_TWO_SURFACES)
    # Elementary tags become discrete entities; physical tags become groups.
    ex=_execute_geo_source("Merge \"$msh_path\";\n")
    m=ex.model
    @test sort!(collect(keys(m.discrete)))==[(2,5),(2,7)]
    @test m.physical[(2,9)]==[5]
    @test m.physical[(2,8)]==[7]
    r5=m.discrete[(2,5)]; r7=m.discrete[(2,7)]
    @test length(r5.element_tags)==2 && length(r7.element_tags)==1
    @test all(t==2 for t in r5.element_types)
    @test length(r5.node_tags)>=3
    # The simplex cells also fold into the mid-file mesh.
    @test ex.mesh!==nothing && size(ex.mesh.tris,2)==3
    @test sort!(collect(ex.mesh.tri_tag))==[8,9,9]
    # `CreateTopology` derives boundary entities like upstream.
    ex=_execute_geo_source("Merge \"$msh_path\";\nCreateTopology;\n")
    m=ex.model
    dims=sort!(collect(keys(m.discrete)))
    @test any(d==2 for (d,_) in dims)
    @test any(d==1 for (d,_) in dims)
    @test any(d==0 for (d,_) in dims)
    # `CreateGeometry` parametrizes the merged surfaces in place.
    ex=_execute_geo_source("Merge \"$msh_path\";\nCreateGeometry;\n")
    @test (2,5) in keys(ex.model.discrete)
    # `CreateGeometry{...}` accepts an explicit shape list.
    ex=_execute_geo_source("Merge \"$msh_path\";\nCreateGeometry{Surface{5};};\n")
    @test (2,5) in keys(ex.model.discrete)
    # `ClassifySurfaces` regroups the faces into classified surfaces.
    ex=_execute_geo_source("Merge \"$msh_path\";\nClassifySurfaces{40*Pi/180,1,0};\n")
    @test any(d==2 for (d,_) in keys(ex.model.discrete))
end

@testset ".geo Merge element-tag and tag-0 edge cases" begin
    dir=mktempdir()
    # Elements without an elementary tag land on entity 0 (upstream rule).
    path=dir*"/noelem.msh"
    write(path,"""
        \$MeshFormat
        2.2 0 8
        \$EndMeshFormat
        \$Nodes
        4
        1 0 0 0
        2 1 0 0
        3 1 1 0
        4 0 1 0
        \$EndNodes
        \$Elements
        2
        1 2 1 9 1 2 3
        2 2 0 1 3 4
        \$EndElements
        """)
    ex=_execute_geo_source("Merge \"$path\";\n")
    m=ex.model
    @test haskey(m.discrete,(2,0))
    @test m.physical[(2,9)]==[0]
    # Physical names merge into the model.
    path2=dir*"/named.msh"
    write(path2,"""
        \$MeshFormat
        2.2 0 8
        \$EndMeshFormat
        \$PhysicalNames
        1
        2 9 "side wall"
        \$EndPhysicalNames
        \$Nodes
        4
        1 0 0 0
        2 1 0 0
        3 1 1 0
        4 0 1 0
        \$EndNodes
        \$Elements
        2
        1 2 2 9 5 1 2 3
        2 2 2 9 5 1 3 4
        \$EndElements
        """)
    ex=_execute_geo_source("Merge \"$path2\";\n")
    @test ex.model.physical_names[(2,9)]=="side wall"
    # Discrete-only models have no meshable entities (upstream remeshes
    # discrete entities — not yet implemented).
    err=_execute_geo_error("Merge \"$path\";\nMesh 2;\n")
    @test err isa ArgumentError
end
