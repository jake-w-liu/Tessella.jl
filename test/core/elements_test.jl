using Test
using SHA
using Tessella

const ElementsUnderTest = Tessella.Elements

const ExpectedSpec = NamedTuple{
    (:family, :dim, :order, :nnodes, :serendipity),
    Tuple{Symbol,Int,Int,Int,Bool},
}
const EXPECTED_SPECS = Dict{Int,ExpectedSpec}()

function add_expected!(family, dim, orders, tags, node_counts, serendipity=false)
    length(orders) == length(tags) == length(node_counts) || error("bad oracle table")
    for (order, tag, nnodes) in zip(orders, tags, node_counts)
        haskey(EXPECTED_SPECS, tag) && error("duplicate oracle tag $tag")
        EXPECTED_SPECS[tag] = (; family, dim, order, nnodes, serendipity)
    end
end

# Independent transcription of Gmsh 4.15.2's fixed-node type map.
add_expected!(:pnt, 0, [0], [15], [1])
add_expected!(:lin, 1, 0:10, [84,1,8,26,27,28,62,63,64,65,66], collect(1:11))
add_expected!(:tri, 2, 0:10, [85,2,9,21,23,25,42,43,44,45,46],
              [1,3,6,10,15,21,28,36,45,55,66])
add_expected!(:qua, 2, 0:10, [86,3,10,36,37,38,47,48,49,50,51],
              [1,4,9,16,25,36,49,64,81,100,121])
add_expected!(:tet, 3, 0:10, [87,4,11,29,30,31,71,72,73,74,75],
              [1,4,10,20,35,56,84,120,165,220,286])
add_expected!(:hex, 3, 0:9, [88,5,12,92,93,94,95,96,97,98],
              [1,8,27,64,125,216,343,512,729,1000])
add_expected!(:pri, 3, 0:9, [89,6,13,90,91,106,107,108,109,110],
              [1,6,18,40,75,126,196,288,405,550])
add_expected!(:pyr, 3, 0:9, [132,7,14,118,119,120,121,122,123,124],
              [1,5,14,30,55,91,140,204,285,385])

add_expected!(:tri, 2, 3:10, [20,22,24,52,53,54,55,56],
              [9,12,15,18,21,24,27,30], true)
add_expected!(:qua, 2, 2:10, [16,39,40,41,57,58,59,60,61],
              [8,12,16,20,24,28,32,36,40], true)
add_expected!(:tet, 3, 3:10, [137,32,33,79,80,81,82,83],
              [16,22,28,34,40,46,52,58], true)
add_expected!(:hex, 3, 2:9, [17,99,100,101,102,103,104,105],
              [20,32,44,56,68,80,92,104], true)
add_expected!(:pri, 3, 2:9, [18,111,112,113,114,115,116,117],
              [15,24,33,42,51,60,69,78], true)
add_expected!(:pyr, 3, 2:9, [19,125,126,127,128,129,130,131],
              [13,21,29,37,45,53,61,69], true)
EXPECTED_SPECS[140] = (family=:trih, dim=3, order=1, nnodes=4,
                       serendipity=false)

const EXPECTED_SPECIAL = Dict(
    34 => (family=:polygon, dim=2, order=1, nnodes=nothing, kind=:decomposed),
    35 => (family=:polyhedron, dim=3, order=1, nnodes=nothing, kind=:decomposed),
    67 => (family=:line_border, dim=1, order=1, nnodes=2, kind=:border),
    68 => (family=:triangle_border, dim=2, order=1, nnodes=3, kind=:border),
    69 => (family=:polygon_border, dim=2, order=1, nnodes=nothing, kind=:border),
    70 => (family=:line_child, dim=1, order=1, nnodes=2, kind=:child),
    133 => (family=:point_xfem, dim=0, order=1, nnodes=1, kind=:subelement),
    134 => (family=:line_xfem, dim=1, order=1, nnodes=2, kind=:subelement),
    135 => (family=:triangle_xfem, dim=2, order=1, nnodes=3, kind=:subelement),
    136 => (family=:tetrahedron_xfem, dim=3, order=1, nnodes=4, kind=:subelement),
    138 => (family=:triangle_mini, dim=2, order=3, nnodes=4, kind=:basis_only),
    139 => (family=:tetrahedron_mini, dim=3, order=3, nnodes=5, kind=:basis_only),
)

const GMSH_FAMILY_NAME = Dict(
    :pnt => "Point", :lin => "Line", :tri => "Triangle",
    :qua => "Quadrangle", :tet => "Tetrahedron", :hex => "Hexahedron",
    :pri => "Prism", :pyr => "Pyramid",
)
const GMSH_NAME_PREFIX = Dict(
    :pnt => "Point", :lin => "Line", :tri => "Triangle",
    :qua => "Quadrilateral", :tet => "Tetrahedron", :hex => "Hexahedron",
    :pri => "Prism", :pyr => "Pyramid",
)
const PRIMARY_NODES = Dict(
    :pnt => 1, :lin => 2, :tri => 3, :qua => 4, :tet => 4,
    :hex => 8, :pri => 6, :pyr => 5,
)
const PROPERTY_API_GAPS = Set([90,91,106,107,108,109,110,
                               111,112,113,114,115,116,117,140])
const PRISM_PROPERTY_GAPS = setdiff(PROPERTY_API_GAPS, Set([140]))

function gmsh_python_environment()
    python_name = get(ENV, "TESSELLA_PYTHON", "python3")
    python = Sys.which(python_name)
    python === nothing && error("Python executable '$python_name' was not found")

    paths = String[]
    explicit = get(ENV, "TESSELLA_GMSH_PYTHONPATH", "")
    if !isempty(explicit)
        push!(paths, explicit)
    else
        gmsh = Sys.which("gmsh")
        if gmsh !== nothing
            candidate = normpath(joinpath(dirname(realpath(gmsh)), "..", "lib"))
            isfile(joinpath(candidate, "gmsh.py")) && push!(paths, candidate)
        end
    end
    inherited = get(ENV, "PYTHONPATH", "")
    isempty(inherited) || push!(paths, inherited)
    return python, join(paths, Sys.iswindows() ? ';' : ':')
end

function run_gmsh_oracle()
    tags = join(sort!(collect(keys(EXPECTED_SPECS))), ",")
    shapes = String[]
    for tag in sort!(collect(keys(EXPECTED_SPECS)))
        spec = EXPECTED_SPECS[tag]
        spec.family === :trih && continue
        name = repr(GMSH_FAMILY_NAME[spec.family])
        serendipity = spec.serendipity ? "True" : "False"
        push!(shapes, "($tag,$name,$(spec.order),$serendipity)")
    end
    shape_table = join(shapes, ",")
    special_tags = join(sort!(collect(keys(EXPECTED_SPECIAL))), ",")

    script = """
import gmsh

TAGS = [$tags]
SHAPES = [$shape_table]
SPECIAL_TAGS = [$special_tags]

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
try:
    print("VERSION\\t" + gmsh.__version__)
    for element_type in TAGS:
        try:
            name, dim, order, nnodes, coords, nprimary = \\
                gmsh.model.mesh.getElementProperties(element_type)
            flat = ",".join(format(float(x), ".17g") for x in coords)
            print("PROP\\t{}\\tOK\\t{}\\t{}\\t{}\\t{}\\t{}\\t{}".format(
                element_type, name.replace("\\t", " "), dim, order,
                nnodes, nprimary, flat))
        except Exception:
            print("PROP\\t{}\\tERROR".format(element_type))

    for element_type in SPECIAL_TAGS:
        try:
            name, dim, order, nnodes, coords, nprimary = \\
                gmsh.model.mesh.getElementProperties(element_type)
            print("SPECIAL\\t{}\\tOK\\t{}\\t{}\\t{}\\t{}\\t{}".format(
                element_type, name.replace("\\t", " "), dim, order,
                nnodes, nprimary))
        except Exception:
            print("SPECIAL\\t{}\\tERROR".format(element_type))

    for expected, family, order, serendipity in SHAPES:
        actual = gmsh.model.mesh.getElementType(family, order, serendipity)
        print("TYPE\\t{}\\t{}".format(expected, actual))

    # getElementProperties currently throws while constructing high-order prism
    # closures. Generate a straight one-prism mesh instead: its global (x,y,z)
    # equal local (u,v,(w+1)/2), independently exposing every node in order.
    for order in range(3, 10):
        for incomplete in (False, True):
            expected = gmsh.model.mesh.getElementType("Prism", order, incomplete)
            gmsh.clear()
            gmsh.model.add("reference_prism")
            geo = gmsh.model.geo
            p = [geo.addPoint(0, 0, 0), geo.addPoint(1, 0, 0),
                 geo.addPoint(0, 1, 0)]
            lines = [geo.addLine(p[0], p[1]), geo.addLine(p[1], p[2]),
                     geo.addLine(p[2], p[0])]
            surface = geo.addPlaneSurface([geo.addCurveLoop(lines)])
            geo.extrude([(2, surface)], 0, 0, 1, numElements=[1],
                        recombine=True)
            geo.synchronize()
            for _, curve in gmsh.model.getEntities(1):
                gmsh.model.mesh.setTransfiniteCurve(curve, 2)
            gmsh.model.mesh.setTransfiniteSurface(surface)
            gmsh.model.mesh.generate(3)
            gmsh.option.setNumber("Mesh.SecondOrderIncomplete", int(incomplete))
            gmsh.model.mesh.setOrder(order)

            types, element_tags, connectivities = gmsh.model.mesh.getElements(3)
            volume = [(int(t), ts, ns) for t, ts, ns in
                      zip(types, element_tags, connectivities) if len(ts)]
            if len(volume) != 1 or volume[0][0] != expected or len(volume[0][1]) != 1:
                raise RuntimeError("did not generate the expected single prism")
            connectivity = volume[0][2]
            node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
            node_index = {int(tag): i for i, tag in enumerate(node_tags)}
            flat = []
            for tag in connectivity:
                i = node_index[int(tag)]
                flat.extend(float(x) for x in node_coords[3*i:3*i+3])
            text = ",".join(format(x, ".17g") for x in flat)
            print("PRISM\\t{}\\t{}\\t{}".format(expected, len(connectivity), text))
finally:
    gmsh.finalize()
"""

    python, pythonpath = gmsh_python_environment()
    command = Cmd([python, "-c", script])
    isempty(pythonpath) || (command = addenv(command, "PYTHONPATH" => pythonpath))
    output = read(command, String)

    version = ""
    properties = Dict{Int,Any}()
    special_properties = Dict{Int,Any}()
    type_lookups = Dict{Int,Int}()
    prism_nodes = Dict{Int,Tuple{Int,Vector{Float64}}}()
    for line in eachline(IOBuffer(output))
        fields = split(line, '\t'; keepempty=true)
        if fields[1] == "VERSION"
            version = fields[2]
        elseif fields[1] == "PROP"
            tag = parse(Int, fields[2])
            if fields[3] == "ERROR"
                properties[tag] = nothing
            else
                coordinates = isempty(fields[9]) ? Float64[] :
                    parse.(Float64, split(fields[9], ','))
                properties[tag] = (
                    name=fields[4], dim=parse(Int, fields[5]),
                    order=parse(Int, fields[6]), nnodes=parse(Int, fields[7]),
                    nprimary=parse(Int, fields[8]), coordinates=coordinates,
                )
            end
        elseif fields[1] == "SPECIAL"
            tag = parse(Int, fields[2])
            if fields[3] == "ERROR"
                special_properties[tag] = nothing
            else
                special_properties[tag] = (
                    name=fields[4], dim=parse(Int, fields[5]),
                    order=parse(Int, fields[6]), nnodes=parse(Int, fields[7]),
                    nprimary=parse(Int, fields[8]),
                )
            end
        elseif fields[1] == "TYPE"
            type_lookups[parse(Int, fields[2])] = parse(Int, fields[3])
        elseif fields[1] == "PRISM"
            prism_nodes[parse(Int, fields[2])] = (
                parse(Int, fields[3]), parse.(Float64, split(fields[4], ',')))
        else
            error("unexpected Gmsh oracle output: $line")
        end
    end
    return (; version, properties, special_properties, type_lookups, prism_nodes)
end

@testset "Gmsh fixed-node catalog" begin
    @test length(EXPECTED_SPECS) == 125
    @test Set(keys(ElementsUnderTest.MSH_CATALOG)) == Set(keys(EXPECTED_SPECS))
    for tag in sort!(collect(keys(EXPECTED_SPECS)))
        expected = EXPECTED_SPECS[tag]
        actual = ElementsUnderTest.msh_spec(tag)
        @test (actual.family, actual.dim, actual.order, actual.nnodes,
               actual.serendipity) ==
              (expected.family, expected.dim, expected.order, expected.nnodes,
               expected.serendipity)
        @test ElementsUnderTest.msh_num_nodes(tag) == expected.nnodes
        @test ElementsUnderTest.msh_dimension(tag) == expected.dim
        @test ElementsUnderTest.msh_order(tag) == expected.order
        @test ElementsUnderTest.msh_family(tag) == expected.family

        coordinates = ElementsUnderTest.lagrange_nodes(tag)
        @test size(coordinates) == (3, expected.nnodes)
        @test all(isfinite, coordinates)
        @test length(Set(Tuple(column) for column in eachcol(coordinates))) ==
              expected.nnodes
        @test ElementsUnderTest.lagrange_nodes(
            expected.family, expected.order;
            serendipity=expected.serendipity) == coordinates
    end

    # MTrihedron.h::getNode is the sole local-node authority for type 140.
    @test ElementsUnderTest.lagrange_nodes(140) ==
          Float64[-1 1 1 -1; -1 -1 1 1; 0 0 0 0]
end

@testset "special tags and integer bounds" begin
    @test ElementsUnderTest.MSH_SPECIAL_TYPES == EXPECTED_SPECIAL
    for tag in keys(EXPECTED_SPECIAL)
        @test_throws ArgumentError ElementsUnderTest.msh_spec(tag)
        @test_throws ArgumentError ElementsUnderTest.lagrange_nodes(tag)
    end
    @test_throws ArgumentError ElementsUnderTest.msh_spec(0)
    @test_throws ArgumentError ElementsUnderTest.msh_spec(141)
    @test_throws ArgumentError ElementsUnderTest.msh_spec(big(typemax(Int)) + 1)
    @test_throws ArgumentError ElementsUnderTest.lagrange_nodes(:lin, -1)
    @test_throws ArgumentError ElementsUnderTest.lagrange_nodes(
        :lin, big(typemax(Int)) + 1)
    @test_throws ArgumentError ElementsUnderTest.lagrange_nodes(:tri, 11)
    @test_throws ArgumentError ElementsUnderTest.lagrange_nodes(
        :tri, 2; serendipity=true)
end

let source_root = get(ENV, "TESSELLA_GMSH_SOURCE", "")
    # The installed 4.15.2 API differential below is mandatory. A full source
    # checkout is an additional provenance gate when explicitly supplied; package
    # tests must not depend on an ephemeral machine-local `/tmp` checkout.
    if isempty(source_root)
        local_checkout = "/tmp/tessella-gmsh-U09p6u"
        isdir(local_checkout) && (source_root = local_checkout)
    end
    if isempty(source_root) || !isdir(source_root)
        @info "TESSELLA_GMSH_SOURCE not set; skipping optional pinned-source hashes"
    else
        @testset "pinned Gmsh source oracle" begin
            expected_hashes = Dict(
                "src/numeric/pointsGenerators.cpp" =>
                    "d24ced74ba29242507d676fe614b6294b909ddf2ff729397519aff1c7cd76699",
                "src/numeric/ElementType.cpp" =>
                    "8e470f6c2dfed971e464fb098823e6609fab5d7875b1cf01532ed099577b9489",
                "src/common/GmshDefines.h" =>
                    "5e1d894ee09b43e644fa3449945378c569f2dca27292127fe643371a7f4cc1a4",
                "src/geo/MTrihedron.h" =>
                    "afddb334f64d8a42e08e7a01f6d1c20fb407c40debc87cb4542dc9409c77c7c4",
                "src/geo/MElement.cpp" =>
                    "969b99203046f3a6f8bc74c8228f1bc62cb2388946c7ab463c2b90297b0d16b6",
                "src/geo/MElementCut.h" =>
                    "c13f6393b7ffa2baf4524072223bad9fa9054fd41270d749199ef07ad7fb2e74",
                "src/geo/MSubElement.h" =>
                    "67642d2b7d4179c232f5573c6f0ebc82d3c5b6b524c010fe651dc41eaca050e9",
                "src/geo/GModelIO_MSH2.cpp" =>
                    "0f90244e233849fe2620e4a53a34a725d7678359dd2bd600b8f0a05b4819bd31",
                "src/geo/GModelIO_MSH4.cpp" =>
                    "183993dca10f546cac5f3304a73beb5f3e83c4b12aeb76929ee16a7f98d6d7cd",
                "src/numeric/BasisFactory.cpp" =>
                    "aaf7ab81db649473ecc9351dcba04c30300ff5e43bd6b7d62e6087296295eb58",
                "src/numeric/miniBasis.cpp" =>
                    "90166cbf75c0f93f7e7412627843dd8247563b5e38cc2708704336cf70e111db",
            )
            for (relative_path, expected_hash) in expected_hashes
                path = joinpath(source_root, relative_path)
                @test isfile(path)
                isfile(path) && @test bytes2hex(sha256(read(path))) == expected_hash
            end
            element_source=read(joinpath(source_root,"src/geo/MElement.cpp"),String)
            for (symbol,width) in (
                ("MSH_LIN_B",2),("MSH_LIN_C",2),("MSH_TRI_B",3),
                ("MSH_POLYG_",0),("MSH_POLYG_B",0),("MSH_POLYH_",0),
                ("MSH_PNT_SUB",1),("MSH_LIN_SUB",2),
                ("MSH_TRI_SUB",3),("MSH_TET_SUB",4),
            )
                pattern=Regex("case "*symbol*":\\s+if\\(name\\) \\*name = [^;]+;\\s+return "*
                              string(width)*";")
                @test occursin(pattern,element_source)
            end
            @test !occursin("MSH_TRI_MINI",element_source)
            @test !occursin("MSH_TET_MINI",element_source)
            cut_source=read(joinpath(source_root,"src/geo/MElementCut.h"),String)
            @test occursin("return _parts.size() * 3",cut_source)
            @test occursin("return _parts.size() * 4",cut_source)
            basis_source=read(
                joinpath(source_root,"src/numeric/BasisFactory.cpp"),String)
            @test occursin("tag == MSH_TRI_MINI",basis_source)
            @test occursin("tag == MSH_TET_MINI",basis_source)
        end
    end
end

@testset "installed Gmsh 4.15.2 differential" begin
    oracle = run_gmsh_oracle()
    @test oracle.version == "4.15.2"
    @test Set(keys(oracle.properties)) == Set(keys(EXPECTED_SPECS))
    @test Set(tag for (tag, value) in oracle.properties if value === nothing) ==
          PROPERTY_API_GAPS

    @test Set(keys(oracle.special_properties)) == Set(keys(EXPECTED_SPECIAL))
    @test Set(tag for (tag, value) in oracle.special_properties if value === nothing) ==
          Set([67,68,70,138,139])
    expected_special_properties = Dict(
        34 => ("Polygon", 2),
        35 => ("Polyhedron", 3),
        69 => ("Polygon Border", 2),
        133 => ("Point Xfem", 0),
        134 => ("Line Xfem", 1),
        135 => ("Triangle Xfem", 2),
        136 => ("Tetrahedron Xfem", 3),
    )
    for (tag, (name, dim)) in expected_special_properties
        actual = oracle.special_properties[tag]
        @test (actual.name, actual.dim, actual.order,
               actual.nnodes, actual.nprimary) == (name, dim, 1, 0, 0)
    end

    @test Set(keys(oracle.type_lookups)) == setdiff(Set(keys(EXPECTED_SPECS)), Set([140]))
    for (expected, actual) in oracle.type_lookups
        @test actual == expected
    end

    for tag in sort!(collect(setdiff(Set(keys(EXPECTED_SPECS)),
                                     PROPERTY_API_GAPS)))
        expected = EXPECTED_SPECS[tag]
        actual = oracle.properties[tag]
        @test startswith(actual.name, GMSH_NAME_PREFIX[expected.family])
        @test actual.dim == expected.dim
        @test actual.order == expected.order
        @test actual.nnodes == expected.nnodes
        @test actual.nprimary == PRIMARY_NODES[expected.family]

        coordinate_dimension = max(expected.dim, 1)
        coordinates = ElementsUnderTest.lagrange_nodes(tag)
        wanted = vec(coordinates[1:coordinate_dimension, :])
        @test length(actual.coordinates) == length(wanted)
        @test isapprox(actual.coordinates, wanted; atol=8eps(Float64),
                       rtol=8eps(Float64))
    end

    @test Set(keys(oracle.prism_nodes)) == PRISM_PROPERTY_GAPS
    for tag in sort!(collect(PRISM_PROPERTY_GAPS))
        expected = EXPECTED_SPECS[tag]
        nnodes, global_coordinates = oracle.prism_nodes[tag]
        @test nnodes == expected.nnodes
        coordinates = copy(ElementsUnderTest.lagrange_nodes(tag))
        coordinates[3, :] .= (coordinates[3, :] .+ 1) ./ 2
        @test length(global_coordinates) == length(coordinates)
        @test isapprox(global_coordinates, vec(coordinates); atol=2e-8, rtol=2e-8)
    end
end

function mixed_io_fixture(; escaped_name::Bool=true)
    coordinates = Float64[
        -2 0 1 2 3 5 6 5 8 9 8 8;
         0 0 0 0 0 0 0 1 0 0 1 0;
         0 0 0 0 0 0 0 0 0 0 0 1
    ]
    blocks = ElementsUnderTest.ElementBlock[
        ElementsUnderTest.ElementBlock(15, reshape(Int32[1], 1, 1), Int32[9]),
        ElementsUnderTest.ElementBlock(1, Int32[2 4; 3 5], Int32[6, 8]),
        ElementsUnderTest.ElementBlock(2, reshape(Int32[6,7,8], 3, 1), Int32[5]),
        ElementsUnderTest.ElementBlock(4, reshape(Int32[9,10,11,12], 4, 1), Int32[7]),
    ]
    names = Dict(
        (0,9) => "probe",
        (1,6) => "wire A",
        (1,8) => "wire B",
        (2,5) => escaped_name ? "surface \"quoted\" \\ path\nline\tend" : "surface",
        (3,7) => "volume",
    )
    return ElementsUnderTest.MixedMesh(coordinates, blocks; physical_names=names)
end

function special_v2_fixture()
    coordinates=Float64[
        0 1 1 0 0 1 1 0;
        0 0 1 1 0 0 1 1;
        0 0 0 0 1 1 1 1
    ]
    domains(a,b)=reshape(Any[a,b],2,1)
    blocks=Any[
        ElementsUnderTest.ElementBlock(1,Int32[1 3;2 4],Int32[1,2]),
        ElementsUnderTest.ElementBlock(2,Int32[1 2;2 4;3 3],Int32[3,4]),
        ElementsUnderTest.ElementBlock(4,Int32[1 2;2 4;3 3;5 8],Int32[5,6]),
        ElementsUnderTest.SpecialElementBlock(
            34,[Int32[1,2,3,1,3,4]],Int32[7];parent_refs=[(2,1)]),
        ElementsUnderTest.SpecialElementBlock(
            35,[Int32[1,2,3,5,2,3,5,7]],Int32[8];parent_refs=[(3,1)]),
        ElementsUnderTest.SpecialElementBlock(
            67,reshape(Int32[1,2],2,1),Int32[9];
            domain_refs=domains((2,1),(2,2))),
        ElementsUnderTest.SpecialElementBlock(
            68,reshape(Int32[1,2,3],3,1),Int32[10];
            domain_refs=domains((3,1),(3,2))),
        ElementsUnderTest.SpecialElementBlock(
            69,[Int32[1,2,3,1,3,4]],Int32[11];
            domain_refs=domains((3,1),(3,2))),
        ElementsUnderTest.SpecialElementBlock(
            70,reshape(Int32[3,4],2,1),Int32[12];parent_refs=[(1,1)]),
        ElementsUnderTest.SpecialElementBlock(
            133,reshape(Int32[1],1,1),Int32[13];parent_refs=[(1,1)]),
        ElementsUnderTest.SpecialElementBlock(
            134,reshape(Int32[1,2],2,1),Int32[14];parent_refs=[(1,1)]),
        ElementsUnderTest.SpecialElementBlock(
            135,reshape(Int32[1,2,3],3,1),Int32[15];parent_refs=[(2,1)]),
        ElementsUnderTest.SpecialElementBlock(
            136,reshape(Int32[1,2,3,5],4,1),Int32[16];parent_refs=[(3,1)]),
    ]
    return ElementsUnderTest.MixedMesh(coordinates,blocks)
end

function special_v2_binary_fixture()
    full=special_v2_fixture()
    blocks=Any[full.blocks[1],full.blocks[2],full.blocks[3],
        ElementsUnderTest.SpecialElementBlock(
            67,reshape(Int32[1,2],2,1),Int32[9]),
        ElementsUnderTest.SpecialElementBlock(
            68,reshape(Int32[1,2,3],3,1),Int32[10]),
        full.blocks[9],full.blocks[10],full.blocks[11],full.blocks[12],full.blocks[13]]
    return ElementsUnderTest.MixedMesh(full.coords,blocks)
end

function special_v4_fixture()
    coordinates=Float64[0 1 0 0;0 0 1 0;0 0 0 1]
    blocks=Any[
        ElementsUnderTest.SpecialElementBlock(
            67,reshape(Int32[1,2],2,1),Int32[1]),
        ElementsUnderTest.SpecialElementBlock(
            68,reshape(Int32[1,2,3],3,1),Int32[2]),
        ElementsUnderTest.SpecialElementBlock(
            70,reshape(Int32[1,2],2,1),Int32[3]),
        ElementsUnderTest.SpecialElementBlock(
            133,reshape(Int32[1],1,1),Int32[4]),
        ElementsUnderTest.SpecialElementBlock(
            134,reshape(Int32[1,2],2,1),Int32[5]),
        ElementsUnderTest.SpecialElementBlock(
            135,reshape(Int32[1,2,3],3,1),Int32[6]),
        ElementsUnderTest.SpecialElementBlock(
            136,reshape(Int32[1,2,3,4],4,1),Int32[7]),
    ]
    return ElementsUnderTest.MixedMesh(coordinates,blocks)
end

function special_v4_metadata_fixture()
    records=[
        (msh=133,tag=4,coordinates=reshape(Float64[9,0,0],3,1)),
        (msh=67,tag=1,coordinates=Float64[0 1;0 0;0 0]),
        (msh=70,tag=3,coordinates=Float64[6 7;0 0;0 0]),
        (msh=134,tag=5,coordinates=Float64[11 12;0 0;0 0]),
        (msh=68,tag=2,coordinates=Float64[3 4 3;0 0 1;0 0 0]),
        (msh=135,tag=6,coordinates=Float64[14 15 14;0 0 1;0 0 0]),
        (msh=136,tag=7,coordinates=Float64[17 18 17 17;0 0 1 0;0 0 0 1]),
    ]
    coordinates=hcat((record.coordinates for record in records)...)
    blocks=Any[]
    entities=Dict{Tuple{Int,Int},ElementsUnderTest.MixedEntity}()
    node_entities=Tuple{Int,Int32}[]
    physical_names=Dict{Tuple{Int,Int},String}()
    position=1
    for record in records
        dim=EXPECTED_SPECIAL[record.msh].dim
        count=size(record.coordinates,2)
        connectivity=reshape(Int32.(position:position+count-1),count,1)
        push!(blocks,ElementsUnderTest.SpecialElementBlock(
            record.msh,connectivity,Int32[record.tag]))
        lower=ntuple(d->minimum(record.coordinates[d,:]),3)
        upper=ntuple(d->maximum(record.coordinates[d,:]),3)
        entities[(dim,record.tag)]=ElementsUnderTest.MixedEntity(
            dim,record.tag,(lower...,upper...);physical_tags=[record.tag])
        append!(node_entities,fill((dim,Int32(record.tag)),count))
        physical_names[(dim,record.tag)]="special $(record.msh)"
        position+=count
    end
    data=ElementsUnderTest.MixedEntityData(entities;
        node_entities=node_entities,
        node_parametric=fill(nothing,size(coordinates,2)),
        external_node_tags=UInt64.(101:100+size(coordinates,2)),
        block_entities=[[Int32(record.tag)] for record in records],
        external_element_tags=[[UInt64(200+i)] for i in eachindex(records)])
    return ElementsUnderTest.MixedMesh(
        coordinates,blocks;physical_names=physical_names,entity_data=data)
end

function legacy_crc(mesh)
    projection=ElementsUnderTest.MixedMesh(
        mesh.coords,mesh.blocks;physical_names=mesh.physical_names)
    return ElementsUnderTest.mixed_crc(projection)
end

function exhaustive_catalog_mesh(tags=keys(ElementsUnderTest.MSH_CATALOG))
    tags = sort!(collect(tags))
    count = sum(ElementsUnderTest.msh_num_nodes(tag) for tag in tags)
    coordinates = Matrix{Float64}(undef, 3, count)
    blocks = ElementsUnderTest.ElementBlock[]
    first_node = 1
    for (index, tag) in pairs(tags)
        local_coordinates = ElementsUnderTest.lagrange_nodes(tag)
        nlocal = size(local_coordinates, 2)
        range = first_node:first_node+nlocal-1
        coordinates[:, range] = local_coordinates
        coordinates[1, range] .+= 4index
        connectivity = reshape(Int32.(range), nlocal, 1)
        push!(blocks, ElementsUnderTest.ElementBlock(
            tag, connectivity, Int32[mod(index, 13)]))
        first_node += nlocal
    end
    return ElementsUnderTest.MixedMesh(coordinates, blocks)
end

function gmsh_check(path)
    executable = Sys.which("gmsh")
    executable === nothing && error("gmsh executable was not found")
    output = IOBuffer()
    process = run(pipeline(ignorestatus(
        `$executable $path -check -parse_and_exit -v 5`),
        stdout=output, stderr=output))
    return success(process), String(take!(output))
end

function gmsh_rewrite(input,output;version::Float64,binary::Bool)
    script = raw"""
import gmsh
import sys

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
try:
    gmsh.open(sys.argv[1])
    gmsh.option.setNumber("Mesh.MshFileVersion", float(sys.argv[3]))
    gmsh.option.setNumber("Mesh.Binary", int(sys.argv[4]))
    gmsh.write(sys.argv[2])
    print(gmsh.__version__)
finally:
    gmsh.finalize()
"""
    python,pythonpath=gmsh_python_environment()
    command=Cmd([python,"-c",script,input,output,string(version),binary ? "1" : "0"])
    isempty(pythonpath) || (command=addenv(command,"PYTHONPATH"=>pythonpath))
    return strip(read(command,String))
end

function gmsh_element_summary_without_finalize(path;probe_tags=Int[])
    # Gmsh 4.15.2 constructs type 69 correctly but crashes while destroying it.
    # Exit the subprocess directly after the independent API has parsed and
    # enumerated the mesh. This is a parser probe only: normal Gmsh-compatible
    # output is rejected because a normal process lifecycle is unsafe.
    script = raw"""
import gmsh
import os
import sys

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.logger.start()
gmsh.open(sys.argv[1])
node_tags, _, _ = gmsh.model.mesh.getNodes()
types, element_blocks, _ = gmsh.model.mesh.getElements()
counts = {int(t): len(tags) for t, tags in zip(types, element_blocks)}
print("VERSION\t" + gmsh.__version__)
print("NODES\t{}".format(len(node_tags)))
print("ELEMENTS\t{}".format(sum(counts.values())))
for element_type, count in sorted(counts.items()):
    print("TYPE\t{}\t{}".format(element_type, count))
for tag in (int(x) for x in sys.argv[2].split(",") if x):
    element_type, nodes, _, _ = gmsh.model.mesh.getElement(tag)
    print("PROBE\t{}\t{}\t{}".format(tag, element_type, len(nodes)))
for message in gmsh.logger.get():
    print("LOG\t" + message.replace("\t", " "))
sys.stdout.flush()
os._exit(0)
"""
    python,pythonpath=gmsh_python_environment()
    probes=join(probe_tags,",")
    command=Cmd([python,"-c",script,path,probes])
    isempty(pythonpath) || (command=addenv(command,"PYTHONPATH"=>pythonpath))
    version=""; nodes=0; elements=0; counts=Dict{Int,Int}()
    probed=Dict{Int,Tuple{Int,Int}}(); messages=String[]
    for line in eachline(IOBuffer(read(command,String)))
        fields=split(line,'\t')
        if fields[1]=="VERSION"
            version=fields[2]
        elseif fields[1]=="NODES"
            nodes=parse(Int,fields[2])
        elseif fields[1]=="ELEMENTS"
            elements=parse(Int,fields[2])
        elseif fields[1]=="TYPE"
            counts[parse(Int,fields[2])]=parse(Int,fields[3])
        elseif fields[1]=="PROBE"
            probed[parse(Int,fields[2])]=(
                parse(Int,fields[3]),parse(Int,fields[4]))
        elseif fields[1]=="LOG"
            push!(messages,fields[2])
        else
            error("unexpected Gmsh summary output: $line")
        end
    end
    return (;version,nodes,elements,counts,probed,messages)
end

function gmsh_physical_element_counts(path)
    script = raw"""
import gmsh
import sys

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
try:
    gmsh.open(sys.argv[1])
    print("VERSION\t" + gmsh.__version__)
    for dim, physical in sorted(gmsh.model.getPhysicalGroups()):
        counts = {}
        for entity in gmsh.model.getEntitiesForPhysicalGroup(dim, physical):
            types, tags, _ = gmsh.model.mesh.getElements(dim, entity)
            for element_type, element_tags in zip(types, tags):
                counts[int(element_type)] = counts.get(int(element_type), 0) + len(element_tags)
        for element_type, count in sorted(counts.items()):
            print("COUNT\t{}\t{}\t{}\t{}".format(
                dim, physical, element_type, count))
finally:
    gmsh.finalize()
"""
    python, pythonpath = gmsh_python_environment()
    command = Cmd([python, "-c", script, path])
    isempty(pythonpath) || (command = addenv(command, "PYTHONPATH" => pythonpath))
    version = ""
    counts = Dict{NTuple{3,Int},Int}()
    for line in eachline(IOBuffer(read(command, String)))
        fields = split(line, '\t')
        if fields[1] == "VERSION"
            version = fields[2]
        elseif fields[1] == "COUNT"
            counts[(parse(Int, fields[2]), parse(Int, fields[3]),
                    parse(Int, fields[4]))] = parse(Int, fields[5])
        else
            error("unexpected Gmsh MSH oracle output: $line")
        end
    end
    return version, counts
end

function gmsh_physical_name(path, dimension, tag)
    script = raw"""
import gmsh
import sys

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
try:
    gmsh.open(sys.argv[1])
    print(gmsh.__version__)
    print(gmsh.model.getPhysicalName(int(sys.argv[2]), int(sys.argv[3])))
finally:
    gmsh.finalize()
"""
    python,pythonpath=gmsh_python_environment()
    command=Cmd([python,"-c",script,path,string(dimension),string(tag)])
    isempty(pythonpath) || (command=addenv(command,"PYTHONPATH"=>pythonpath))
    lines=split(chomp(read(command,String)),'\n';keepempty=true)
    length(lines)==2 || error("unexpected Gmsh physical-name oracle output")
    return lines[1],lines[2]
end

function gmsh_physical_entities(path)
    script = raw"""
import gmsh
import sys

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
try:
    gmsh.open(sys.argv[1])
    print("VERSION\t" + gmsh.__version__)
    for dim, physical in sorted(gmsh.model.getPhysicalGroups()):
        entities = gmsh.model.getEntitiesForPhysicalGroup(dim, physical)
        print("ENTITY\t{}\t{}\t{}".format(
            dim, physical, ",".join(str(int(x)) for x in entities)))
finally:
    gmsh.finalize()
"""
    python,pythonpath=gmsh_python_environment()
    command=Cmd([python,"-c",script,path])
    isempty(pythonpath) || (command=addenv(command,"PYTHONPATH"=>pythonpath))
    version=""; groups=Dict{Tuple{Int,Int},Vector{Int}}()
    for line in eachline(IOBuffer(read(command,String)))
        fields=split(line,'\t';keepempty=true)
        if fields[1]=="VERSION"
            version=fields[2]
        elseif fields[1]=="ENTITY"
            entities=isempty(fields[4]) ? Int[] : parse.(Int,split(fields[4],','))
            groups[(parse(Int,fields[2]),parse(Int,fields[3]))]=entities
        else
            error("unexpected Gmsh physical-entity oracle output: $line")
        end
    end
    return version,groups
end

function gmsh_external_tags(path)
    script = raw"""
import gmsh
import sys

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
try:
    gmsh.open(sys.argv[1])
    nodes, _, _ = gmsh.model.mesh.getNodes()
    _, element_blocks, _ = gmsh.model.mesh.getElements()
    elements = [int(x) for block in element_blocks for x in block]
    print("VERSION\t" + gmsh.__version__)
    print("NODES\t" + ",".join(str(int(x)) for x in sorted(nodes)))
    print("ELEMENTS\t" + ",".join(str(x) for x in sorted(elements)))
finally:
    gmsh.finalize()
"""
    python,pythonpath=gmsh_python_environment()
    command=Cmd([python,"-c",script,path])
    isempty(pythonpath) || (command=addenv(command,"PYTHONPATH"=>pythonpath))
    lines=collect(eachline(IOBuffer(read(command,String))))
    length(lines)==3 || error("unexpected Gmsh external-tag oracle output")
    parse_tags(line)=begin
        fields=split(line,'\t';keepempty=true)
        isempty(fields[2]) ? Int[] : parse.(Int,split(fields[2],','))
    end
    return split(lines[1],'\t')[2],parse_tags(lines[2]),parse_tags(lines[3])
end

function write_gmsh_reference_mesh(directory)
    script = raw"""
import gmsh
import os
import sys

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
try:
    gmsh.model.add("mixed_reference")
    cases = [
        (1, 11, 6, 1, [101, 305],
         [0, 0, 0, 1, 0, 0]),
        (2, 22, 5, 2, [701, 702, 999],
         [3, 0, 0, 4, 0, 0, 3, 1, 0]),
        (3, 33, 7, 4, [1201, 1202, 1203, 1400],
         [6, 0, 0, 7, 0, 0, 6, 1, 0, 6, 0, 1]),
    ]
    for dim, entity, physical, element_type, nodes, coordinates in cases:
        gmsh.model.addDiscreteEntity(dim, entity)
        gmsh.model.mesh.addNodes(dim, entity, nodes, coordinates)
        gmsh.model.mesh.addElementsByType(
            entity, element_type, [2000 + entity], nodes)
        gmsh.model.addPhysicalGroup(dim, [entity], physical)
        gmsh.model.setPhysicalName(dim, physical, "group {}".format(physical))
    for binary in (0, 1):
        gmsh.option.setNumber("Mesh.Binary", binary)
        for version in (2.2, 4.1):
            gmsh.option.setNumber("Mesh.MshFileVersion", version)
            suffix = "-binary" if binary else ""
            gmsh.write(os.path.join(
                sys.argv[1], "gmsh-reference-{}{}.msh".format(version, suffix)))
    gmsh.clear()
    gmsh.model.add("empty_reference")
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    for binary in (0, 1):
        gmsh.option.setNumber("Mesh.Binary", binary)
        suffix = "-binary" if binary else ""
        gmsh.write(os.path.join(
            sys.argv[1], "gmsh-empty-4.1{}.msh".format(suffix)))
    print(gmsh.__version__)
finally:
    gmsh.finalize()
"""
    python,pythonpath=gmsh_python_environment()
    command=Cmd([python,"-c",script,directory])
    isempty(pythonpath) || (command=addenv(command,"PYTHONPATH"=>pythonpath))
    return strip(read(command,String))
end

@testset "MixedMesh structure, CRC, and simplex conversion" begin
    mesh = mixed_io_fixture()
    diagnostic = ElementsUnderTest.validate(mesh)
    @test diagnostic.ok
    crc = ElementsUnderTest.mixed_crc(mesh)
    @test crc.n_nodes == 12
    @test crc.n_blocks == 4
    @test crc.n_cells == 5
    @test crc.bbox == ((-2.0,0.0,0.0), (9.0,1.0,1.0))
    @test crc.sha == "b219f5afde8b589ce8c31c0fb174ebd2373811ab0f4e37a564f18934499e00c0"

    # Canonical cell ordering is insensitive to blocks/cell iteration, but the
    # local node order is deliberately orientation-sensitive.
    reordered_blocks = ElementsUnderTest.ElementBlock[]
    for block in reverse(mesh.blocks)
        columns = reverse(axes(block.nodes, 2))
        push!(reordered_blocks, ElementsUnderTest.ElementBlock(
            block.msh, block.nodes[:, columns], block.tags[columns]))
    end
    reordered = ElementsUnderTest.MixedMesh(
        mesh.coords, reordered_blocks; physical_names=mesh.physical_names)
    @test ElementsUnderTest.mixed_crc(reordered).sha == crc.sha

    changed_coordinates = copy(mesh.coords)
    changed_coordinates[1,1] = nextfloat(changed_coordinates[1,1])
    @test ElementsUnderTest.mixed_crc(ElementsUnderTest.MixedMesh(
        changed_coordinates, mesh.blocks; physical_names=mesh.physical_names)).sha != crc.sha

    changed_orientation = copy(mesh.blocks)
    line_nodes = copy(changed_orientation[2].nodes)
    line_nodes[1,1], line_nodes[2,1] = line_nodes[2,1], line_nodes[1,1]
    changed_orientation[2] = ElementsUnderTest.ElementBlock(
        1, line_nodes, changed_orientation[2].tags)
    @test ElementsUnderTest.mixed_crc(ElementsUnderTest.MixedMesh(
        mesh.coords, changed_orientation; physical_names=mesh.physical_names)).sha != crc.sha

    changed_tags = copy(mesh.blocks)
    tags = copy(changed_tags[3].tags); tags[1] += 1
    changed_tags[3] = ElementsUnderTest.ElementBlock(2, changed_tags[3].nodes, tags)
    @test ElementsUnderTest.mixed_crc(ElementsUnderTest.MixedMesh(
        mesh.coords, changed_tags; physical_names=mesh.physical_names)).sha != crc.sha

    changed_names = copy(mesh.physical_names); changed_names[(3,7)] = "other"
    @test ElementsUnderTest.mixed_crc(ElementsUnderTest.MixedMesh(
        mesh.coords, mesh.blocks; physical_names=changed_names)).sha != crc.sha

    empty = ElementsUnderTest.MixedMesh(zeros(3,0), ElementsUnderTest.ElementBlock[])
    @test ElementsUnderTest.validate(empty).ok
    @test ElementsUnderTest.mixed_crc(empty).bbox ==
          ((0.0,0.0,0.0), (0.0,0.0,0.0))

    simplex = ElementsUnderTest.mixed_to_simplex(mesh)
    @test size(simplex.segs) == (2,2)
    @test size(simplex.tris) == (3,1)
    @test size(simplex.tets) == (4,1)
    @test simplex.seg_tag == Int32[6,8]
    @test simplex.tri_tag == Int32[5]
    @test simplex.tet_tag == Int32[7]
    without_point = ElementsUnderTest.MixedMesh(
        mesh.coords, mesh.blocks[2:end]; physical_names=mesh.physical_names)
    @test ElementsUnderTest.mixed_crc(
        ElementsUnderTest.simplex_to_mixed(simplex;
            physical_names=mesh.physical_names)).sha ==
          ElementsUnderTest.mixed_crc(without_point).sha
    high_order = ElementsUnderTest.MixedMesh(
        Float64[0 0.5 1; 0 0 0; 0 0 0],
        [ElementsUnderTest.ElementBlock(8,
            reshape(Int32[1,3,2],3,1),Int32[1])])
    @test_throws ArgumentError ElementsUnderTest.mixed_to_simplex(high_order)
end

@testset "MixedMesh structural validation and linear allocation growth" begin
    @test_throws ArgumentError ElementsUnderTest.ElementBlock(
        1, reshape(BigInt[1, big(typemax(Int32))+1],2,1))
    @test_throws ArgumentError ElementsUnderTest.ElementBlock(
        1, reshape(Int32[1,2],2,1), [big(typemax(Int32))+1])
    @test_throws ArgumentError ElementsUnderTest.MixedMesh(
        reshape(Float64[Inf,0,0],3,1), ElementsUnderTest.ElementBlock[])
    @test_throws ArgumentError ElementsUnderTest.MixedMesh(
        zeros(3,1), [ElementsUnderTest.ElementBlock(
            1,reshape(Int32[1,2],2,1))])
    @test_throws ArgumentError ElementsUnderTest.MixedMesh(
        zeros(3,0), [ElementsUnderTest.ElementBlock(
            1,Matrix{Int32}(undef,2,0))])
    @test_throws ArgumentError ElementsUnderTest.MixedMesh(
        zeros(3,0), ElementsUnderTest.ElementBlock[];
        physical_names=Dict((4,1)=>"bad"))
    @test_throws ArgumentError ElementsUnderTest.MixedMesh(
        zeros(3,0), ElementsUnderTest.ElementBlock[];
        physical_names=Dict((1,0)=>"bad"))
    @test_throws ArgumentError ElementsUnderTest.MixedMesh(
        zeros(3,0), ElementsUnderTest.ElementBlock[];
        physical_names=Dict((1,1)=>"bad\0name"))

    appended=ElementsUnderTest.MixedMesh(Float64[0 1;0 0;0 0],
                                          ElementsUnderTest.ElementBlock[])
    block=ElementsUnderTest.ElementBlock(
        1,reshape(Int32[1,2],2,1),Int32[3])
    @test ElementsUnderTest.add_block!(appended,block)===appended
    @test appended.blocks==[block]
    @test_throws ArgumentError ElementsUnderTest.add_block!(
        appended,ElementsUnderTest.ElementBlock(
            1,reshape(Int32[2,3],2,1),Int32[3]))
    @test_throws ArgumentError ElementsUnderTest.add_block!(
        appended,ElementsUnderTest.ElementBlock(
            1,Matrix{Int32}(undef,2,0)))
    mutated=ElementsUnderTest.ElementBlock(
        1,reshape(Int32[1,2],2,1),Int32[3])
    mutated.tags[1]=Int32(-1)
    @test_throws ArgumentError ElementsUnderTest.add_block!(appended,mutated)
    empty!(mutated.tags)
    @test_throws ArgumentError ElementsUnderTest.add_block!(appended,mutated)

    repeated = ElementsUnderTest.MixedMesh(
        Float64[0 1;0 0;0 0],
        [ElementsUnderTest.ElementBlock(1,
            reshape(Int32[1,1],2,1),Int32[1])])
    @test !ElementsUnderTest.validate(repeated).ok
    duplicate = ElementsUnderTest.MixedMesh(
        Float64[0 1;0 0;0 0],
        [ElementsUnderTest.ElementBlock(1,
            Int32[1 2;2 1],Int32[1,2])])
    @test !ElementsUnderTest.validate(duplicate).ok
    @test ElementsUnderTest.validate(
        duplicate; reject_duplicate_cells=false).ok

    function fragmented_lines(count)
        coordinates = Float64[0 1;0 0;0 0]
        blocks = [ElementsUnderTest.ElementBlock(
            1,reshape(Int32[1,2],2,1),Int32[mod(i,7)]) for i in 1:count]
        return ElementsUnderTest.MixedMesh(coordinates, blocks)
    end
    small=fragmented_lines(1_000); large=fragmented_lines(2_000)
    ElementsUnderTest.mixed_to_simplex(small) # compile before measurement
    ElementsUnderTest.mixed_crc(small)
    GC.gc(); allocated_small=@allocated ElementsUnderTest.mixed_to_simplex(small)
    GC.gc(); allocated_large=@allocated ElementsUnderTest.mixed_to_simplex(large)
    @test allocated_large <= 2.25allocated_small + 16_384
    GC.gc(); crc_allocated_small=@allocated ElementsUnderTest.mixed_crc(small)
    GC.gc(); crc_allocated_large=@allocated ElementsUnderTest.mixed_crc(large)
    @test crc_allocated_large <= 2.25crc_allocated_small + 16_384
end

@testset "special-element structure, links, validation, and CRC" begin
    @test Tessella.ElementRef === ElementsUnderTest.ElementRef
    @test Tessella.SpecialElementBlock === ElementsUnderTest.SpecialElementBlock
    missing=ElementsUnderTest.ElementRef()
    @test (missing.block,missing.cell)==(Int32(0),Int32(0))
    @test ElementsUnderTest.ElementRef(2,3)==ElementsUnderTest.ElementRef(2,3)
    @test_throws ArgumentError ElementsUnderTest.ElementRef(0,1)
    @test_throws ArgumentError ElementsUnderTest.ElementRef(-1,-1)
    @test_throws ArgumentError ElementsUnderTest.ElementRef(big(typemax(Int32))+1,1)

    polygon=ElementsUnderTest.SpecialElementBlock(
        34,[Int32[1,2,3,1,3,4]],Int32[7];parent_refs=[(1,1)])
    @test polygon.connectivity==Int32[1,2,3,1,3,4]
    @test polygon.offsets==Int32[1,7]
    @test polygon.parent_refs==[ElementsUnderTest.ElementRef(1,1)]
    @test size(polygon.domain_refs)==(2,1)
    line_border=ElementsUnderTest.SpecialElementBlock(
        67,reshape(Int32[1,2],2,1),Int32[9])
    @test line_border.offsets==Int32[1,3]
    @test_throws ArgumentError ElementsUnderTest.SpecialElementBlock(
        138,reshape(Int32[1,2,3,4],4,1))
    @test_throws ArgumentError ElementsUnderTest.SpecialElementBlock(
        139,reshape(Int32[1,2,3,4,5],5,1))
    @test_throws ArgumentError ElementsUnderTest.SpecialElementBlock(
        34,[Int32[1,2,3,4]])
    @test_throws ArgumentError ElementsUnderTest.SpecialElementBlock(
        67,[Int32[1,2,3]])
    @test_throws ArgumentError ElementsUnderTest.SpecialElementBlock(
        34,Int32[1,2,3],Int32[0,4],Int32[1])
    @test_throws ArgumentError ElementsUnderTest.SpecialElementBlock(
        34,Int32[1,2,3],Int32[1,5],Int32[1])
    @test_throws ArgumentError ElementsUnderTest.SpecialElementBlock(
        34,[Int32[1,2,3]],Int32[-1])
    @test_throws ArgumentError ElementsUnderTest.SpecialElementBlock(
        34,[Int32[1,2,3]];domain_refs=reshape(Any[(1,1),(1,1)],2,1))
    @test_throws ArgumentError ElementsUnderTest.SpecialElementBlock(
        67,reshape(Int32[1,2],2,1);parent_refs=[(1,1)])
    @test_throws ArgumentError ElementsUnderTest.SpecialElementBlock(
        67,reshape(Int32[1,2],2,1);
        domain_refs=reshape(Any[(),(1,1)],2,1))

    valid=special_v2_fixture()
    diagnostic=ElementsUnderTest.validate(valid)
    @test diagnostic.ok
    crc=ElementsUnderTest.mixed_crc(valid)
    @test crc==(n_nodes=8,n_blocks=13,n_cells=16,
        bbox=((0.0,0.0,0.0),(1.0,1.0,1.0)),
        sha="944dbe417200e316bf265d2d6f8a154dfb4c6bec3eb853401f289990a9be4089")
    @test_throws ArgumentError ElementsUnderTest.mixed_to_simplex(valid)

    changed=special_v2_fixture()
    changed.blocks[11].parent_refs[1]=ElementsUnderTest.ElementRef(1,2)
    @test ElementsUnderTest.validate(changed).ok
    @test ElementsUnderTest.mixed_crc(changed).sha!=crc.sha

    cycle=special_v2_fixture()
    cycle.blocks[11].parent_refs[1]=ElementsUnderTest.ElementRef(11,1)
    @test !ElementsUnderTest.validate(cycle).ok
    @test_throws ArgumentError ElementsUnderTest.mixed_crc(cycle)

    @test_throws ArgumentError ElementsUnderTest.MixedMesh(
        valid.coords,[ElementsUnderTest.SpecialElementBlock(
            134,reshape(Int32[1,2],2,1);parent_refs=[(2,1)])])
    @test_throws ArgumentError ElementsUnderTest.MixedMesh(
        valid.coords,[ElementsUnderTest.SpecialElementBlock(
            134,reshape(Int32[1,2],2,1);parent_refs=[(1,2)])])

    shared_nodes=ElementsUnderTest.MixedMesh(
        valid.coords,[ElementsUnderTest.SpecialElementBlock(
            34,[Int32[1,2,3,1,3,4]])])
    @test ElementsUnderTest.validate(shared_nodes).ok
    repeated=ElementsUnderTest.MixedMesh(
        valid.coords,[ElementsUnderTest.SpecialElementBlock(
            34,[Int32[1,1,3]])])
    @test !ElementsUnderTest.validate(repeated).ok
    repeated_part=ElementsUnderTest.MixedMesh(
        valid.coords,[ElementsUnderTest.SpecialElementBlock(
            34,[Int32[1,2,3,3,2,1]])])
    @test !ElementsUnderTest.validate(repeated_part).ok
    duplicate=ElementsUnderTest.MixedMesh(
        valid.coords,[ElementsUnderTest.SpecialElementBlock(34,
            [Int32[1,2,3,1,3,4],Int32[3,2,1,4,3,1]])])
    @test !ElementsUnderTest.validate(duplicate).ok
    @test ElementsUnderTest.validate(
        duplicate;reject_duplicate_cells=false).ok
    fixed_repeated=ElementsUnderTest.MixedMesh(
        valid.coords,[ElementsUnderTest.SpecialElementBlock(
            67,reshape(Int32[1,1],2,1))])
    @test !ElementsUnderTest.validate(fixed_repeated).ok

    coordinates=Float64[0 1 2;0 0 0;0 0 0]
    ordinary=ElementsUnderTest.ElementBlock(
        1,Int32[1 2;2 3],Int32[1,2])
    ordered=ElementsUnderTest.MixedMesh(coordinates,Any[
        ordinary,ElementsUnderTest.SpecialElementBlock(
            134,reshape(Int32[2,3],2,1),Int32[3];parent_refs=[(1,2)])])
    reordered=ElementsUnderTest.MixedMesh(coordinates,Any[
        ElementsUnderTest.SpecialElementBlock(
            134,reshape(Int32[2,3],2,1),Int32[3];parent_refs=[(2,2)]),ordinary])
    @test ElementsUnderTest.mixed_crc(reordered).sha==
          ElementsUnderTest.mixed_crc(ordered).sha

    appended=ElementsUnderTest.MixedMesh(coordinates,Any[ordinary])
    invalid_link=ElementsUnderTest.SpecialElementBlock(
        134,reshape(Int32[1,2],2,1);parent_refs=[(3,1)])
    @test_throws ArgumentError ElementsUnderTest.add_block!(appended,invalid_link)
    @test length(appended.blocks)==1
end

@testset "mixed MSH v2.2/v4.1 round-trip and Gmsh acceptance" begin
    directory=mktempdir()
    mesh=mixed_io_fixture(escaped_name=false)
    expected_crc=ElementsUnderTest.mixed_crc(mesh).sha
    expected_counts=Dict(
        (0,9,15)=>1,
        (1,6,1)=>1,
        (1,8,1)=>1,
        (2,5,2)=>1,
        (3,7,4)=>1,
    )
    entities,groups,node_entity=ElementsUnderTest._mixed_v4_layout(mesh)
    entity_bounds=Dict((entity.dim,Int(entity.physical))=>entity.bounds
                       for entity in entities)
    @test entity_bounds[(0,9)][1:3]==(-2.0,0.0,0.0)
    @test entity_bounds[(1,6)]==(0.0,0.0,0.0,1.0,0.0,0.0)
    @test entity_bounds[(1,8)]==(2.0,0.0,0.0,3.0,0.0,0.0)
    @test entity_bounds[(2,5)]==(5.0,0.0,0.0,6.0,1.0,0.0)
    @test entity_bounds[(3,7)]==(8.0,0.0,0.0,9.0,1.0,1.0)
    @test entity_bounds[(3,0)]==(-2.0,0.0,0.0,9.0,1.0,1.0)
    @test node_entity>0
    @test length(groups)==5
    for version in (2.2,4.1), binary in (false,true)
        mode=binary ? "binary" : "ascii"
        path=joinpath(directory,"mixed-$version-$mode.msh")
        @test ElementsUnderTest.write_mixed_msh(
            path,mesh;version=version,binary=binary)==path
        open(path,"r") do io
            @test readline(io)=="\$MeshFormat"
            fields=split(readline(io))
            @test fields==[string(version),binary ? "1" : "0","8"]
            if binary
                @test read(io,Int32)==1
            end
        end
        back=ElementsUnderTest.read_mixed_msh(path)
        @test back.coords==mesh.coords
        @test back.physical_names==mesh.physical_names
        @test ElementsUnderTest.validate(back).ok
        @test legacy_crc(back).sha==expected_crc
        ok,output=gmsh_check(path)
        @test ok
        @test !occursin("Error",output)
        @test occursin("12 nodes",output)
        @test occursin("5 elements",output)
        oracle_version,counts=gmsh_physical_element_counts(path)
        @test oracle_version=="4.15.2"
        @test counts==expected_counts
    end

    escaped=mixed_io_fixture()
    escaped_crc=ElementsUnderTest.mixed_crc(escaped).sha
    for version in (2.2,4.1), binary in (false,true)
        path=joinpath(directory,"escaped-$version-$binary.msh")
        @test_throws ArgumentError ElementsUnderTest.write_mixed_msh(
            path,escaped;version=version,binary=binary)
        ElementsUnderTest.write_mixed_msh(
            path,escaped;version=version,binary=binary,gmsh_compatible=false)
        @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(path)
        back=ElementsUnderTest.read_mixed_msh(path;tessella_extensions=true)
        @test back.physical_names==escaped.physical_names
        @test legacy_crc(back).sha==escaped_crc
    end

    # Cross-format conversion does not depend on the writer's block grouping.
    v2=joinpath(directory,"cross-v2.msh")
    v4=joinpath(directory,"cross-v4.msh")
    v2_again=joinpath(directory,"cross-v2-again.msh")
    ElementsUnderTest.write_mixed_msh(v2,mesh;version=2.2,binary=true)
    ElementsUnderTest.write_mixed_msh(
        v4,ElementsUnderTest.read_mixed_msh(v2);version=4.1,binary=true)
    ElementsUnderTest.write_mixed_msh(
        v2_again,ElementsUnderTest.read_mixed_msh(v4);version=2.2)
    @test ElementsUnderTest.mixed_crc(
        ElementsUnderTest.read_mixed_msh(v2_again)).sha==expected_crc

    empty=ElementsUnderTest.MixedMesh(zeros(3,0),ElementsUnderTest.ElementBlock[];
                                      physical_names=Dict((3,1)=>"unused"))
    for version in (2.2,4.1), binary in (false,true)
        path=joinpath(directory,"empty-$version-$binary.msh")
        ElementsUnderTest.write_mixed_msh(
            path,empty;version=version,binary=binary)
        back=ElementsUnderTest.read_mixed_msh(path)
        @test legacy_crc(back).sha==
              ElementsUnderTest.mixed_crc(empty).sha
        @test back.physical_names==empty.physical_names
        ok,output=gmsh_check(path)
        @test ok
        @test !occursin("Error",output)
    end

    @test write_gmsh_reference_mesh(directory)=="4.15.2"
    gmsh_reference_crc=nothing
    gmsh_reference_v4_crc=nothing
    for version in (2.2,4.1), binary in (false,true)
        suffix=binary ? "-binary" : ""
        path=joinpath(directory,"gmsh-reference-$version$suffix.msh")
        reference=ElementsUnderTest.read_mixed_msh(path)
        @test [(block.msh,size(block.nodes,2),block.tags) for block in reference.blocks]==[
            (1,1,Int32[6]),(2,1,Int32[5]),(4,1,Int32[7])]
        @test reference.physical_names==Dict(
            (1,6)=>"group 6",(2,5)=>"group 5",(3,7)=>"group 7")
        crc=legacy_crc(reference).sha
        if gmsh_reference_crc===nothing
            gmsh_reference_crc=crc
        else
            @test crc==gmsh_reference_crc
        end
        if version==4.1
            full_crc=ElementsUnderTest.mixed_crc(reference).sha
            if gmsh_reference_v4_crc===nothing
                gmsh_reference_v4_crc=full_crc
            else
                @test full_crc==gmsh_reference_v4_crc
            end
            @test reference.entity_data.external_node_tags==
                  UInt64[101,305,701,702,999,1201,1202,1203,1400]
            @test reference.entity_data.external_element_tags==
                  [UInt64[2011],UInt64[2022],UInt64[2033]]
        end
        ok,output=gmsh_check(path)
        @test ok
        @test !occursin("Error",output)
    end
    @test gmsh_reference_crc==
          "bce5e5f31440ecf9b0765e116210b6584080663061be00482505d8b658341cf9"
    @test gmsh_reference_v4_crc!==nothing
    for binary in (false,true)
        suffix=binary ? "-binary" : ""
        path=joinpath(directory,"gmsh-empty-4.1$suffix.msh")
        reference=ElementsUnderTest.read_mixed_msh(path)
        @test size(reference.coords)==(3,0)
        @test isempty(reference.blocks)
        @test reference.entity_data!==nothing
        @test isempty(reference.entity_data.entities)
        @test ElementsUnderTest.validate(reference).ok
        ok,output=gmsh_check(path)
        @test ok
        @test !occursin("Error",output)
    end
end

@testset "special-element MSH round-trip and Gmsh 4.15.2 acceptance" begin
    directory=mktempdir()
    full=special_v2_fixture()
    expected_crc="944dbe417200e316bf265d2d6f8a154dfb4c6bec3eb853401f289990a9be4089"
    ascii=joinpath(directory,"special-full-v2-ascii.msh")
    @test_throws ArgumentError ElementsUnderTest.write_mixed_msh(
        ascii,full;version=2.2,binary=false)
    @test !isfile(ascii)
    @test ElementsUnderTest.write_mixed_msh(
        ascii,full;version=2.2,binary=false,gmsh_compatible=false)==ascii
    back=ElementsUnderTest.read_mixed_msh(ascii)
    @test ElementsUnderTest.validate(back).ok
    @test ElementsUnderTest.mixed_crc(back).sha==expected_crc
    @test Set(block.msh for block in back.blocks)==
          Set([1,2,4,34,35,67,68,69,70,133,134,135,136])
    summary=gmsh_element_summary_without_finalize(
        ascii;probe_tags=[10,13,16])
    @test summary.version=="4.15.2"
    @test summary.nodes==7
    @test summary.elements==10
    @test summary.counts==Dict(
        1=>1,2=>1,4=>1,67=>1,68=>1,
        70=>1,133=>1,134=>1,135=>1,136=>1)
    @test summary.probed==Dict(10=>(34,4),13=>(35,5),16=>(69,4))
    @test "Info: 8 nodes" in summary.messages
    @test "Info: 16 elements" in summary.messages
    @test any(message->endswith(message,"Done reading '$ascii'"),summary.messages)

    binary_mesh=special_v2_binary_fixture()
    binary_crc="4ef8de3b81b4449aa1b3efc3905f4418e00b5ac96ffd609e6eef2957b1156b68"
    binary=joinpath(directory,"special-fixed-v2-binary.msh")
    ElementsUnderTest.write_mixed_msh(binary,binary_mesh;version=2.2,binary=true)
    binary_back=ElementsUnderTest.read_mixed_msh(binary)
    @test ElementsUnderTest.validate(binary_back).ok
    @test ElementsUnderTest.mixed_crc(binary_back).sha==binary_crc
    @test binary_back.blocks[7].parent_refs==
          [ElementsUnderTest.ElementRef(1,1)]
    ok,output=gmsh_check(binary)
    @test ok
    @test !occursin("Error",output)
    @test occursin("8 nodes",output)
    @test occursin("13 elements",output)

    # Gmsh itself independently reopens and rewrites representative variable,
    # parent-linked and domain-linked MSH2 records. Polygon-border type 69 has
    # only the isolated parser probe above: Gmsh 4.15.2 crashes during a normal
    # close or rewrite, so it is Tessella-only output.
    polygon=ElementsUnderTest.MixedMesh(
        Float64[0 1 1 0;0 0 1 1;0 0 0 0],Any[
            ElementsUnderTest.ElementBlock(
                2,reshape(Int32[1,2,3],3,1),Int32[2]),
            ElementsUnderTest.SpecialElementBlock(
                34,[Int32[1,2,3,1,3,4]],Int32[2];parent_refs=[(1,1)])])
    polyhedron=ElementsUnderTest.MixedMesh(
        Float64[0 1 0 0 0;0 0 1 0 0;0 0 0 1 -1],Any[
            ElementsUnderTest.ElementBlock(
                4,reshape(Int32[1,2,3,4],4,1),Int32[2]),
            ElementsUnderTest.SpecialElementBlock(
                35,[Int32[1,2,3,4,2,1,3,5]],Int32[2];parent_refs=[(1,1)])])
    border=ElementsUnderTest.MixedMesh(
        Float64[1 0 0 1;0 1 0 1;0 0 0 0],Any[
            ElementsUnderTest.ElementBlock(
                2,Int32[3 1;1 4;2 2],Int32[2,2]),
            ElementsUnderTest.SpecialElementBlock(
                67,reshape(Int32[1,2],2,1),Int32[2];
                domain_refs=reshape(Any[(1,1),(1,2)],2,1))])
    subelement=ElementsUnderTest.MixedMesh(
        Float64[0 1;0 0;0 0],Any[
            ElementsUnderTest.ElementBlock(
                1,reshape(Int32[1,2],2,1),Int32[2]),
            ElementsUnderTest.SpecialElementBlock(
                134,reshape(Int32[1,2],2,1),Int32[2];parent_refs=[(1,1)])])
    for (name,mesh,special_type,arity) in (
        ("polygon",polygon,34,6),
        ("polyhedron",polyhedron,35,8),
        ("border",border,67,2),
        ("subelement",subelement,134,2),
    )
        source=joinpath(directory,"gmsh-$name-source.msh")
        rewritten=joinpath(directory,"gmsh-$name-rewritten.msh")
        ElementsUnderTest.write_mixed_msh(source,mesh;version=2.2)
        @test gmsh_rewrite(
            source,rewritten;version=2.2,binary=false)=="4.15.2"
        ok,output=gmsh_check(rewritten)
        @test ok
        @test !occursin("Error",output)
        rewritten_mesh=ElementsUnderTest.read_mixed_msh(rewritten)
        @test ElementsUnderTest.validate(rewritten_mesh).ok
        special=only(filter(block->block.msh==special_type,
                            rewritten_mesh.blocks))
        @test length(special.connectivity)==arity
        if special_type in (34,35,134)
            @test special.parent_refs==[ElementsUnderTest.ElementRef(1,1)]
        else
            @test special.domain_refs[:,1]==[
                ElementsUnderTest.ElementRef(1,1),
                ElementsUnderTest.ElementRef(1,2)]
        end
        @test Set(Tuple(column) for column in eachcol(rewritten_mesh.coords))==
              Set(Tuple(column) for column in eachcol(mesh.coords))
        @test ElementsUnderTest.mixed_crc(rewritten_mesh).sha==
              ElementsUnderTest.mixed_crc(mesh).sha
    end
    for binary_output in (false,true)
        source=joinpath(directory,"gmsh-sub-source-$binary_output.msh")
        rewritten=joinpath(directory,"gmsh-sub-rewritten-$binary_output.msh")
        ElementsUnderTest.write_mixed_msh(source,subelement;version=2.2)
        @test gmsh_rewrite(source,rewritten;
                           version=2.2,binary=binary_output)=="4.15.2"
        @test ElementsUnderTest.mixed_crc(
            ElementsUnderTest.read_mixed_msh(rewritten)).sha==
              ElementsUnderTest.mixed_crc(subelement).sha
    end

    # Pinned Gmsh 4.15.2 preserves distinct parent links when reading Tessella's
    # binary MSH2 and writing ASCII, but its own binary writer leaves the parent
    # in an unwritten blob slot and emits a constant third tag of 1
    # (`MElement.cpp:1447-1454`). Keep this
    # upstream corruption explicit instead of claiming a binary Gmsh rewrite
    # round trip that the pinned implementation cannot provide.
    distinct_parents=ElementsUnderTest.MixedMesh(
        Float64[0 1 2;0 0 0;0 0 0],Any[
            ElementsUnderTest.ElementBlock(
                1,Int32[1 2;2 3],Int32[2,2]),
            ElementsUnderTest.SpecialElementBlock(
                134,Int32[1 2;2 3],Int32[2,2];
                parent_refs=[(1,1),(1,2)])])
    parent_source=joinpath(directory,"distinct-parent-source.msh")
    ElementsUnderTest.write_mixed_msh(
        parent_source,distinct_parents;version=2.2,binary=true)
    parent_source_back=ElementsUnderTest.read_mixed_msh(parent_source)
    @test parent_source_back.blocks[2].parent_refs==[
        ElementsUnderTest.ElementRef(1,1),
        ElementsUnderTest.ElementRef(1,2)]
    parent_ascii=joinpath(directory,"distinct-parent-gmsh-ascii.msh")
    @test gmsh_rewrite(parent_source,parent_ascii;
                       version=2.2,binary=false)=="4.15.2"
    @test ElementsUnderTest.read_mixed_msh(parent_ascii).blocks[2].parent_refs==[
        ElementsUnderTest.ElementRef(1,1),
        ElementsUnderTest.ElementRef(1,2)]
    parent_binary=joinpath(directory,"distinct-parent-gmsh-binary.msh")
    @test gmsh_rewrite(parent_source,parent_binary;
                       version=2.2,binary=true)=="4.15.2"
    @test ElementsUnderTest.read_mixed_msh(parent_binary).blocks[2].parent_refs==[
        ElementsUnderTest.ElementRef(1,1),
        ElementsUnderTest.ElementRef(1,1)]

    v4_raw=special_v4_fixture()
    @test ElementsUnderTest.mixed_crc(v4_raw).sha==
          "c43a87878c4f94bb352de2b16b97b86c274685f6b594d99af1ee5239493e680e"
    v4_names=Dict(
        (1,1)=>"line border",(2,2)=>"triangle border",(1,3)=>"line child",
        (0,4)=>"point xfem",(1,5)=>"line xfem",(2,6)=>"triangle xfem",
        (3,7)=>"tetrahedron xfem")
    v4=ElementsUnderTest.MixedMesh(
        v4_raw.coords,v4_raw.blocks;physical_names=v4_names)
    v4_crc=ElementsUnderTest.mixed_crc(v4).sha
    for binary_output in (false,true)
        mode=binary_output ? "binary" : "ascii"
        path=joinpath(directory,"special-v4-$mode.msh")
        @test_throws ArgumentError ElementsUnderTest.write_mixed_msh(
            path,v4;version=4.1,binary=binary_output)
        @test !isfile(path)
        ElementsUnderTest.write_mixed_msh(
            path,v4;version=4.1,binary=binary_output,gmsh_compatible=false)
        v4_back=ElementsUnderTest.read_mixed_msh(path)
        @test ElementsUnderTest.validate(v4_back).ok
        @test legacy_crc(v4_back).sha==v4_crc
        @test v4_back.physical_names==v4_names
        @test v4_back.entity_data!==nothing
        @test v4_back.entity_data.external_node_tags==UInt64[1,2,3,4]
        @test sort!(reduce(vcat,v4_back.entity_data.external_element_tags))==
              UInt64[1,2,3,4,5,6,7]
        metadata_path=joinpath(directory,"special-v4-metadata-$mode.msh")
        @test_throws ArgumentError ElementsUnderTest.write_mixed_msh(
            metadata_path,v4_back;version=4.1,binary=!binary_output)
        @test !isfile(metadata_path)
        ElementsUnderTest.write_mixed_msh(
            metadata_path,v4_back;version=4.1,binary=!binary_output,
            gmsh_compatible=false)
        metadata_back=ElementsUnderTest.read_mixed_msh(metadata_path)
        @test ElementsUnderTest.mixed_crc(metadata_back).sha==
              ElementsUnderTest.mixed_crc(v4_back).sha
    end


    # With one neutral physical projection, Gmsh rewrites these fixed-width
    # special records in both v4 encodings. Its writer emits empty node blocks
    # for entities that own elements but no classified nodes; those blocks are
    # legal, bounded by max_blocks, and must not be mistaken for truncation.
    neutral_blocks=Any[
        ElementsUnderTest.SpecialElementBlock(
            block.msh,block.connectivity,block.offsets,
            zeros(Int32,length(block.tags))) for block in v4_raw.blocks]
    neutral_v4=ElementsUnderTest.MixedMesh(v4_raw.coords,neutral_blocks)
    neutral_crc=ElementsUnderTest.mixed_crc(neutral_v4).sha
    neutral_source=joinpath(directory,"gmsh-neutral-v4-source.msh")
    ElementsUnderTest.write_mixed_msh(neutral_source,neutral_v4;version=4.1)
    for binary_output in (false,true)
        mode=binary_output ? "binary" : "ascii"
        rewritten=joinpath(directory,"gmsh-neutral-v4-$mode.msh")
        @test gmsh_rewrite(neutral_source,rewritten;
                           version=4.1,binary=binary_output)=="4.15.2"
        rewritten_mesh=ElementsUnderTest.read_mixed_msh(rewritten)
        @test ElementsUnderTest.validate(rewritten_mesh).ok
        @test legacy_crc(rewritten_mesh).sha==neutral_crc
        ok,output=gmsh_check(rewritten)
        @test ok
        @test !occursin("Error",output)
    end


    classified_v4=special_v4_metadata_fixture()
    classified_crc="4948d49bd52c04b3fb755ba8baf3dc7569cf8d5d140ce960dc9c1886434be7a3"
    @test ElementsUnderTest.validate(classified_v4).ok
    @test ElementsUnderTest.mixed_crc(classified_v4).sha==classified_crc
    classified_source=joinpath(directory,"gmsh-classified-v4-source.msh")
    ElementsUnderTest.write_mixed_msh(
        classified_source,classified_v4;version=4.1)
    for binary_output in (false,true)
        mode=binary_output ? "binary" : "ascii"
        rewritten=joinpath(directory,"gmsh-classified-v4-$mode.msh")
        @test gmsh_rewrite(classified_source,rewritten;
                           version=4.1,binary=binary_output)=="4.15.2"
        rewritten_mesh=ElementsUnderTest.read_mixed_msh(rewritten)
        @test ElementsUnderTest.validate(rewritten_mesh).ok
        @test ElementsUnderTest.mixed_crc(rewritten_mesh).sha==classified_crc
        @test rewritten_mesh.physical_names==classified_v4.physical_names
        @test rewritten_mesh.entity_data.external_node_tags==UInt64.(101:117)
        @test sort!(reduce(vcat,rewritten_mesh.entity_data.external_element_tags))==
              UInt64.(201:207)
        ok,output=gmsh_check(rewritten)
        @test ok
        @test !occursin("Error",output)
    end

    atomic=joinpath(directory,"special-atomic.msh")
    write(atomic,"sentinel")
    @test_throws ArgumentError ElementsUnderTest.write_mixed_msh(
        atomic,full;version=2.2,binary=true)
    @test read(atomic,String)=="sentinel"
    @test_throws ArgumentError ElementsUnderTest.write_mixed_msh(
        atomic,full;version=4.1,binary=false)
    @test read(atomic,String)=="sentinel"
    @test_throws ArgumentError ElementsUnderTest.write_mixed_msh(
        atomic,binary_mesh;version=4.1,binary=true)
    @test read(atomic,String)=="sentinel"
    @test_throws ArgumentError ElementsUnderTest.write_mixed_msh(
        atomic,border;version=2.2,binary=true)
    @test read(atomic,String)=="sentinel"
end

@testset "lossless v4 entity metadata and multiple physical memberships" begin
    directory=mktempdir()
    source="""
    \$MeshFormat
    4.1 0 8
    \$EndMeshFormat
    \$PhysicalNames
    2
    1 7 "shared"
    1 8 "selected"
    \$EndPhysicalNames
    \$Entities
    3 2 0 0
    1 0 0 0 0
    2 1 0 0 0
    3 2 0 0 0
    11 0 0 0 1 0 0 2 7 8 2 1 -2
    22 1 0 0 2 0 0 1 7 2 2 -3
    \$EndEntities
    \$Nodes
    3 3 10 30
    0 1 0 1
    10
    0 0 0
    0 2 0 1
    20
    1 0 0
    0 3 0 1
    30
    2 0 0
    \$EndNodes
    \$Elements
    2 2 101 202
    1 11 1 1
    101 10 20
    1 22 1 1
    202 20 30
    \$EndElements
    """
    input=joinpath(directory,"entity-input.msh"); write(input,source)
    mesh=ElementsUnderTest.read_mixed_msh(input)
    @test ElementsUnderTest.validate(mesh).ok
    @test mesh.blocks[1].tags==Int32[7,7]
    data=mesh.entity_data
    @test data!==nothing
    @test Set(keys(data.entities))==Set([(0,1),(0,2),(0,3),(1,11),(1,22)])
    @test data.entities[(1,11)].bbox==(0.0,0.0,0.0,1.0,0.0,0.0)
    @test data.entities[(1,11)].physical_tags==Int32[7,8]
    @test data.entities[(1,11)].boundaries==Int32[1,-2]
    @test data.entities[(1,22)].physical_tags==Int32[7]
    @test data.entities[(1,22)].boundaries==Int32[2,-3]
    @test data.external_node_tags==[10,20,30]
    @test data.node_entities==[(0,Int32(1)),(0,Int32(2)),(0,Int32(3))]
    @test data.node_parametric==[nothing,nothing,nothing]
    @test data.block_entities==[Int32[11,22]]
    @test data.external_element_tags==[[101,202]]

    input_ok,input_output=gmsh_check(input)
    @test input_ok
    @test !occursin("Error",input_output)
    output=joinpath(directory,"entity-output.msh")
    ElementsUnderTest.write_mixed_msh(output,mesh;version=4.1)
    output_ok,output_text=gmsh_check(output)
    @test output_ok
    @test !occursin("Error",output_text)
    version,physical_entities=gmsh_physical_entities(output)
    @test version=="4.15.2"
    @test physical_entities==Dict((1,7)=>[11,22],(1,8)=>[11])
    oracle_version,physical_counts=gmsh_physical_element_counts(output)
    @test oracle_version=="4.15.2"
    @test physical_counts==Dict((1,7,1)=>2,(1,8,1)=>1)
    external_version,external_nodes,external_elements=gmsh_external_tags(output)
    @test external_version=="4.15.2"
    @test external_nodes==[10,20,30]
    @test external_elements==[101,202]

    roundtrip=ElementsUnderTest.read_mixed_msh(output)
    @test roundtrip.entity_data.external_node_tags==[10,20,30]
    @test roundtrip.entity_data.external_element_tags==[[101,202]]
    @test roundtrip.entity_data.block_entities==[Int32[11,22]]
    @test roundtrip.entity_data.entities[(1,11)].physical_tags==Int32[7,8]
    @test roundtrip.entity_data.entities[(1,11)].boundaries==Int32[1,-2]
    crc=ElementsUnderTest.mixed_crc(mesh).sha
    @test crc=="c08800cbabb0eb2ffc048481320b9d12e7f3ad5ffce0ffd25beb75f6bd3d608a"
    @test ElementsUnderTest.mixed_crc(roundtrip).sha==crc

    binary_output=joinpath(directory,"entity-output-binary.msh")
    ElementsUnderTest.write_mixed_msh(
        binary_output,mesh;version=4.1,binary=true)
    binary_roundtrip=ElementsUnderTest.read_mixed_msh(binary_output)
    @test ElementsUnderTest.validate(binary_roundtrip).ok
    @test ElementsUnderTest.mixed_crc(binary_roundtrip).sha==crc
    @test binary_roundtrip.entity_data.external_node_tags==[10,20,30]
    @test binary_roundtrip.entity_data.external_element_tags==[[101,202]]
    @test binary_roundtrip.entity_data.node_entities==data.node_entities
    @test binary_roundtrip.entity_data.block_entities==[Int32[11,22]]
    @test binary_roundtrip.entity_data.entities[(1,11)].physical_tags==Int32[7,8]
    @test binary_roundtrip.entity_data.entities[(1,11)].boundaries==Int32[1,-2]
    binary_ok,binary_text=gmsh_check(binary_output)
    @test binary_ok
    @test !occursin("Error",binary_text)
    binary_version,binary_entities=gmsh_physical_entities(binary_output)
    @test binary_version=="4.15.2"
    @test binary_entities==Dict((1,7)=>[11,22],(1,8)=>[11])
    binary_external_version,binary_nodes,binary_elements=
        gmsh_external_tags(binary_output)
    @test binary_external_version=="4.15.2"
    @test binary_nodes==[10,20,30]
    @test binary_elements==[101,202]

    function rebuild(;entities=data.entities,node_entities=data.node_entities,
                     node_parametric=data.node_parametric,
                     external_node_tags=data.external_node_tags,
                     block_entities=data.block_entities,
                     external_element_tags=data.external_element_tags)
        return ElementsUnderTest.MixedEntityData(entities;
            node_entities=node_entities,node_parametric=node_parametric,
            external_node_tags=external_node_tags,block_entities=block_entities,
            external_element_tags=external_element_tags)
    end
    function with_data(replacement;blocks=mesh.blocks)
        return ElementsUnderTest.MixedMesh(mesh.coords,blocks;
            physical_names=mesh.physical_names,entity_data=replacement)
    end
    function changed_entity(key;bounds=data.entities[key].bbox,
                            physical_tags=data.entities[key].physical_tags,
                            boundaries=data.entities[key].boundaries)
        entities=copy(data.entities); old=data.entities[key]
        entities[key]=ElementsUnderTest.MixedEntity(old.dim,old.tag,bounds;
            physical_tags=physical_tags,boundaries=boundaries)
        return rebuild(entities=entities)
    end
    variants=ElementsUnderTest.MixedMesh[]
    push!(variants,with_data(changed_entity((1,11);physical_tags=Int32[7,8,9])))
    push!(variants,with_data(changed_entity(
        (1,11);bounds=(0.0,0.0,0.0,1.25,0.0,0.0))))
    push!(variants,with_data(changed_entity((1,11);boundaries=Int32[-1,-2])))
    changed_nodes=copy(data.external_node_tags); changed_nodes[1]=11
    push!(variants,with_data(rebuild(external_node_tags=changed_nodes)))
    classifications=copy(data.node_entities); classifications[1]=(0,Int32(2))
    push!(variants,with_data(rebuild(node_entities=classifications)))
    cell_entities=deepcopy(data.block_entities); cell_entities[1][1]=Int32(22)
    push!(variants,with_data(rebuild(block_entities=cell_entities)))
    element_tags=deepcopy(data.external_element_tags); element_tags[1][1]=303
    push!(variants,with_data(rebuild(external_element_tags=element_tags)))
    parameters=deepcopy(data.node_parametric); parameters[1]=Float64[]
    push!(variants,with_data(rebuild(node_parametric=parameters)))
    @test all(variant->ElementsUnderTest.validate(variant).ok,variants)
    @test all(variant->ElementsUnderTest.mixed_crc(variant).sha!=crc,variants)

    reordered_block=ElementsUnderTest.ElementBlock(
        mesh.blocks[1].msh,mesh.blocks[1].nodes[:,end:-1:1],
        mesh.blocks[1].tags[end:-1:1])
    reordered_data=rebuild(
        block_entities=[data.block_entities[1][end:-1:1]],
        external_element_tags=[data.external_element_tags[1][end:-1:1]])
    reordered=with_data(reordered_data;blocks=[reordered_block])
    @test ElementsUnderTest.mixed_crc(reordered).sha==crc

    duplicate_membership=with_data(changed_entity(
        (1,11);physical_tags=Int32[7,8,7]))
    duplicate_path=joinpath(directory,"duplicate-membership.msh")
    ElementsUnderTest.write_mixed_msh(duplicate_path,duplicate_membership;version=4.1)
    duplicate_back=ElementsUnderTest.read_mixed_msh(duplicate_path)
    @test duplicate_back.entity_data.entities[(1,11)].physical_tags==Int32[7,8,7]
    duplicate_ok,duplicate_output=gmsh_check(duplicate_path)
    @test duplicate_ok
    @test !occursin("Error",duplicate_output)

    point=data.entities[(0,1)]
    @test_throws ArgumentError ElementsUnderTest.MixedEntity(
        0,1,(-0.0,0.0,0.0,0.0,0.0,0.0))
    @test_throws ArgumentError ElementsUnderTest.MixedEntityData(
        Dict((0.0,1)=>point))
    @test_throws ArgumentError ElementsUnderTest.MixedEntityData(
        Dict((0,1)=>point);node_entities=[(0.0,1)])

    missing_cell=with_data(data)
    missing_cell.entity_data.block_entities[1][1]=Int32(99)
    @test !ElementsUnderTest.validate(missing_cell).ok
    duplicate_node=with_data(data)
    duplicate_node.entity_data.external_node_tags[2]=10
    @test !ElementsUnderTest.validate(duplicate_node).ok
    wrong_parameters=with_data(data)
    wrong_parameters.entity_data.node_parametric[1]=Float64[0]
    @test !ElementsUnderTest.validate(wrong_parameters).ok
    bad_boundary=with_data(data)
    push!(bad_boundary.entity_data.entities[(1,11)].boundaries,Int32(99))
    @test !ElementsUnderTest.validate(bad_boundary).ok
    copied_blocks=[ElementsUnderTest.ElementBlock(
        block.msh,block.nodes,block.tags) for block in mesh.blocks]
    wrong_projection=with_data(data;blocks=copied_blocks)
    wrong_projection.blocks[1].tags[1]=Int32(8)
    @test !ElementsUnderTest.validate(wrong_projection).ok
    @test_throws ArgumentError ElementsUnderTest.add_block!(
        mesh,ElementsUnderTest.ElementBlock(
            1,reshape(Int32[1,2],2,1),Int32[7]))
end

@testset "repeated v4 entity/element sections and full-width external tags" begin
    directory=mktempdir()
    repeated="""
    \$MeshFormat
    4.1 0 8
    \$EndMeshFormat
    \$PhysicalNames
    2
    1 7 "shared"
    1 8 "selected"
    \$EndPhysicalNames
    \$Entities
    3 1 0 0
    1 0 0 0 0
    2 1 0 0 0
    3 2 0 0 0
    11 0 0 0 1 0 0 2 7 8 2 1 -2
    \$EndEntities
    \$Entities
    0 1 0 0
    22 1 0 0 2 0 0 1 7 2 2 -3
    \$EndEntities
    \$Nodes
    3 3 10 30
    0 1 0 1
    10
    0 0 0
    0 2 0 1
    20
    1 0 0
    0 3 0 1
    30
    2 0 0
    \$EndNodes
    \$Elements
    1 1 101 101
    1 11 1 1 101 10 20
    \$EndElements
    \$Elements
    1 1 202 202
    1 22 1 1 202 20 30
    \$EndElements
    """
    repeated_path=joinpath(directory,"repeated-v4-sections.msh")
    write(repeated_path,repeated)
    mesh=ElementsUnderTest.read_mixed_msh(repeated_path)
    @test ElementsUnderTest.validate(mesh).ok
    @test mesh.entity_data.entities[(1,11)].physical_tags==Int32[7,8]
    @test mesh.entity_data.entities[(1,11)].boundaries==Int32[1,-2]
    @test mesh.entity_data.entities[(1,22)].physical_tags==Int32[7]
    @test mesh.entity_data.entities[(1,22)].boundaries==Int32[2,-3]
    @test mesh.entity_data.external_node_tags==UInt64[10,20,30]
    @test mesh.entity_data.block_entities==[Int32[11,22]]
    @test mesh.entity_data.external_element_tags==[UInt64[101,202]]
    @test ElementsUnderTest.mixed_crc(mesh).sha==
          "c08800cbabb0eb2ffc048481320b9d12e7f3ad5ffce0ffd25beb75f6bd3d608a"
    @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(
        repeated_path;max_entities=4)
    @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(
        repeated_path;max_nodes=2)
    @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(
        repeated_path;max_elements=1)
    @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(
        repeated_path;max_blocks=2)
    ok,output=gmsh_check(repeated_path)
    @test ok
    @test !occursin("Error",output)
    version,nodes,elements=gmsh_external_tags(repeated_path)
    @test version=="4.15.2"
    @test nodes==[10,20,30]
    @test elements==[101,202]

    canonical_path=joinpath(directory,"canonical-v4-sections.msh")
    ElementsUnderTest.write_mixed_msh(canonical_path,mesh;version=4.1)
    canonical_lines=readlines(canonical_path)
    @test count(==("\$Entities"),canonical_lines)==1
    @test count(==("\$Nodes"),canonical_lines)==1
    @test count(==("\$Elements"),canonical_lines)==1
    @test ElementsUnderTest.mixed_crc(
        ElementsUnderTest.read_mixed_msh(canonical_path)).sha==
        ElementsUnderTest.mixed_crc(mesh).sha
    ok,output=gmsh_check(canonical_path)
    @test ok
    @test !occursin("Error",output)

    packed_entities="""
    \$MeshFormat
    4.1 0 8
    \$EndMeshFormat
    \$PhysicalNames
    2
    0 7 "shared points"
    0 8 "selected point"
    \$EndPhysicalNames
    \$Entities
    2 0 0 0
    1 0 0 0 1 7 2 1 0 0 2 7 8
    \$EndEntities
    \$Nodes
    2 2 10 20
    0 1 0 1
    10
    0 0 0
    0 2 0 1
    20
    1 0 0
    \$EndNodes
    \$Elements
    0 0 0 0
    \$EndElements
    """
    packed_entities_path=joinpath(directory,"packed-v4-entities.msh")
    write(packed_entities_path,packed_entities)
    packed_mesh=ElementsUnderTest.read_mixed_msh(packed_entities_path)
    @test packed_mesh.entity_data.entities[(0,1)].physical_tags==Int32[7]
    @test packed_mesh.entity_data.entities[(0,2)].physical_tags==Int32[7,8]
    @test packed_mesh.entity_data.external_node_tags==UInt64[10,20]
    @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(
        packed_entities_path;max_entities=1)
    ok,output=gmsh_check(packed_entities_path)
    @test ok
    @test !occursin("Error",output)
    version,nodes,elements=gmsh_external_tags(packed_entities_path)
    @test version=="4.15.2"
    @test nodes==[10,20]
    @test isempty(elements)

    high_node=UInt64(typemax(Int))+one(UInt64)
    high_element=high_node+one(UInt64)
    full_width="""
    \$MeshFormat
    4.1 0 8
    \$EndMeshFormat
    \$Entities
    1 0 0 0
    1 -0 0 0 1 9
    \$EndEntities
    \$Nodes
    1 1 $high_node $high_node
    0 1 0 1
    $high_node
    -0 0 0
    \$EndNodes
    \$Elements
    1 1 $high_element $high_element
    0 1 15 1
    $high_element $high_node
    \$EndElements
    """
    full_width_path=joinpath(directory,"full-width-tags.msh")
    write(full_width_path,full_width)
    full_width_mesh=ElementsUnderTest.read_mixed_msh(full_width_path)
    @test signbit(full_width_mesh.entity_data.entities[(0,1)].bbox[1])
    @test signbit(full_width_mesh.entity_data.entities[(0,1)].bbox[4])
    @test full_width_mesh.entity_data.external_node_tags==UInt64[high_node]
    @test full_width_mesh.entity_data.external_element_tags==[UInt64[high_element]]
    @test full_width_mesh.blocks[1].tags==Int32[9]
    full_width_outputs=String[]
    for binary in (false,true)
        output_path=joinpath(directory,"full-width-tags-output-$binary.msh")
        ElementsUnderTest.write_mixed_msh(
            output_path,full_width_mesh;version=4.1,binary=binary)
        output_mesh=ElementsUnderTest.read_mixed_msh(output_path)
        @test output_mesh.entity_data.external_node_tags==UInt64[high_node]
        @test output_mesh.entity_data.external_element_tags==[UInt64[high_element]]
        @test ElementsUnderTest.mixed_crc(output_mesh).sha==
              ElementsUnderTest.mixed_crc(full_width_mesh).sha
        push!(full_width_outputs,output_path)
    end
    for path in (full_width_path,full_width_outputs...)
        ok,output=gmsh_check(path)
        @test ok
        @test !occursin("Error",output)
        @test occursin("1 node",output) && occursin("1 element",output)
    end
    @test_throws ArgumentError ElementsUnderTest.MixedEntityData(
        Dict{Tuple{Int,Int},ElementsUnderTest.MixedEntity}();
        external_node_tags=[big(typemax(UInt64))+1])
end

@testset "all 125 fixed-node types survive explicit Tessella-only MSH round-trip" begin
    mesh=exhaustive_catalog_mesh()
    @test ElementsUnderTest.validate(mesh).ok
    expected=ElementsUnderTest.mixed_crc(mesh).sha
    @test expected=="17f74d0e0185e18bf87eb4d987fe820af43d0ca28c66d29ce7989892aecd3933"
    @test Set(ElementsUnderTest.GMSH_4_15_2_MSH_READER_GAPS_V4)==Set([
        84,85,86,87,88,100,101,102,103,104,105,
        125,126,127,128,129,130,131,132,
    ])
    @test Set(ElementsUnderTest.GMSH_4_15_2_MSH_READER_GAPS_V2)==union(
        Set(ElementsUnderTest.GMSH_4_15_2_MSH_READER_GAPS_V4),Set([69,89,140]))
    directory=mktempdir()
    for version in (2.2,4.1), binary in (false,true)
        path=joinpath(directory,"all-types-$version-$binary.msh")
        @test_throws ArgumentError ElementsUnderTest.write_mixed_msh(
            path,mesh;version=version,binary=binary)
        ElementsUnderTest.write_mixed_msh(
            path,mesh;version=version,binary=binary,gmsh_compatible=false)
        back=ElementsUnderTest.read_mixed_msh(path)
        @test Set(block.msh for block in back.blocks)==Set(keys(EXPECTED_SPECS))
        @test sum(size(block.nodes,2) for block in back.blocks)==125
        @test ElementsUnderTest.validate(back).ok
        @test legacy_crc(back).sha==expected
    end
    for version in (2.2,4.1), binary in (false,true)
        gaps=version==2.2 ? ElementsUnderTest.GMSH_4_15_2_MSH_READER_GAPS_V2 :
                           ElementsUnderTest.GMSH_4_15_2_MSH_READER_GAPS_V4
        compatible=exhaustive_catalog_mesh(setdiff(Set(keys(EXPECTED_SPECS)),gaps))
        expected_count=length(setdiff(Set(keys(EXPECTED_SPECS)),Set(gaps)))
        @test sum(size(block.nodes,2) for block in compatible.blocks)==expected_count
        path=joinpath(directory,"gmsh-compatible-types-$version-$binary.msh")
        ElementsUnderTest.write_mixed_msh(
            path,compatible;version=version,binary=binary)
        ok,output=gmsh_check(path)
        @test ok
        @test !occursin("Error",output)
        @test occursin("$expected_count elements",output)
    end
end

@testset "sparse external tags and v4 entity physical tags" begin
    directory=mktempdir()
    v2="""
    \$MeshFormat
    2.2 0 8
    \$EndMeshFormat
    \$PhysicalNames
    1
    2 5 "surface"
    \$EndPhysicalNames
    \$Nodes
    3
    100 0 0 0
    7 1 0 0
    999 0 1 0
    \$EndNodes
    \$Elements
    1
    77 2 2 5 17 999 100 7
    \$EndElements
    """
    v4="""
    \$MeshFormat
    4.1 0 8
    \$EndMeshFormat
    \$PhysicalNames
    1
    2 5 "surface"
    \$EndPhysicalNames
    \$Entities
    0 0 1 1
    17 0 0 0 1 1 0 1 5 0
    42 0 0 0 1 1 0 0 0
    \$EndEntities
    \$Nodes
    1 3 7 999
    3 42 0 3
    100
    7
    999
    0 0 0
    1 0 0
    0 1 0
    \$EndNodes
    \$Elements
    1 1 1234 1234
    2 17 2 1
    1234 999 100 7
    \$EndElements
    """
    for (name,text) in (("sparse-v2.msh",v2),("sparse-v4.msh",v4))
        path=joinpath(directory,name); write(path,text)
        mesh=ElementsUnderTest.read_mixed_msh(path)
        @test mesh.coords==Float64[0 1 0;0 0 1;0 0 0]
        @test length(mesh.blocks)==1
        @test mesh.blocks[1].msh==2
        @test mesh.blocks[1].nodes[:,1]==Int32[3,1,2]
        @test mesh.blocks[1].tags==Int32[5]
        @test mesh.physical_names==Dict((2,5)=>"surface")
    end
    packed_v4="""
    \$MeshFormat
    4.1 0 8
    \$EndMeshFormat
    \$PhysicalNames
    1
    2 5 "surface"
    \$EndPhysicalNames
    \$Entities
    0 0 1 0
    17 0 0 0 1 1 0 1 5 0
    \$EndEntities
    \$Nodes
    1 3 7 999 2 17 0 3 100 7 999 0 0 0 1 0 0 0 1 0
    \$EndNodes
    \$Elements
    1 1 1234 1234
    2 17 2 1
    1234 999 100 7
    \$EndElements
    """
    path=joinpath(directory,"packed-v4-nodes.msh"); write(path,packed_v4)
    packed=ElementsUnderTest.read_mixed_msh(path)
    @test packed.coords==Float64[0 1 0;0 0 1;0 0 0]
    @test packed.blocks[1].nodes[:,1]==Int32[3,1,2]
    @test packed.blocks[1].tags==Int32[5]
    @test packed.physical_names==Dict((2,5)=>"surface")
    ok,output=gmsh_check(path)
    @test ok
    @test !occursin("Error",output)

    repeated_names="""
    \$MeshFormat
    2.2 0 8
    \$EndMeshFormat
    \$PhysicalNames
    1
    1 3 "curve"
    \$EndPhysicalNames
    \$PhysicalNames
    2
    1 3 "curve"
    2 5 "surface"
    \$EndPhysicalNames
    \$Nodes
    0
    \$EndNodes
    """
    path=joinpath(directory,"repeated-physical-names.msh"); write(path,repeated_names)
    @test ElementsUnderTest.read_mixed_msh(path).physical_names==
          Dict((1,3)=>"curve",(2,5)=>"surface")
    @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(
        path;max_physical_names=2)
    ok,output=gmsh_check(path)
    @test ok
    @test !occursin("Error",output)
    @test gmsh_physical_name(path,1,3)==("4.15.2","curve")
    @test gmsh_physical_name(path,2,5)==("4.15.2","surface")

    conflicting_names=replace(repeated_names,"1 3 \"curve\"\n2 5"=>
                                              "1 3 \"other\"\n2 5")
    path=joinpath(directory,"conflicting-physical-names.msh")
    write(path,conflicting_names)
    @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(path)

    literal_backslash="""
    \$MeshFormat
    2.2 0 8
    \$EndMeshFormat
    \$PhysicalNames
    1
    1 3 "C:\\new\\tab\\mesh"
    \$EndPhysicalNames
    \$Nodes
    0
    \$EndNodes
    """
    path=joinpath(directory,"literal-backslash.msh"); write(path,literal_backslash)
    @test ElementsUnderTest.read_mixed_msh(path).physical_names[(1,3)]==
          raw"C:\new\tab\mesh"
    @test gmsh_physical_name(path,1,3)==("4.15.2",raw"C:\new\tab\mesh")
    @test ElementsUnderTest.read_mixed_msh(
        path;tessella_extensions=true).physical_names[(1,3)]==
          "C:\n"*"ew\t"*"ab\\mesh"
    implicit_entity="""
    \$MeshFormat
    4.1 0 8
    \$EndMeshFormat
    \$Nodes
    1 4 1 4
    3 42 0 4
    1
    2
    3
    4
    0 0 0
    1 0 0
    0 1 0
    0 0 1
    \$EndNodes
    \$Elements
    1 1 1 1
    3 42 4 1
    1 1 2 3 4
    \$EndElements
    """
    path=joinpath(directory,"implicit-entity.msh"); write(path,implicit_entity)
    implicit_mesh=ElementsUnderTest.read_mixed_msh(path)
    @test implicit_mesh.blocks[1].msh==4
    @test implicit_mesh.blocks[1].tags==Int32[0]
    @test isempty(implicit_mesh.entity_data.entities)
    @test implicit_mesh.entity_data.node_entities==fill((3,Int32(42)),4)
    @test implicit_mesh.entity_data.block_entities==[Int32[42]]
    ok,output=gmsh_check(path)
    @test ok
    @test !occursin("Error",output)
    implicit_output=joinpath(directory,"implicit-entity-output.msh")
    ElementsUnderTest.write_mixed_msh(implicit_output,implicit_mesh;version=4.1)
    implicit_roundtrip=ElementsUnderTest.read_mixed_msh(implicit_output)
    @test ElementsUnderTest.mixed_crc(implicit_roundtrip).sha==
          ElementsUnderTest.mixed_crc(implicit_mesh).sha
    ok,output=gmsh_check(implicit_output)
    @test ok
    @test !occursin("Error",output)
    parametric_nodes="""
    \$MeshFormat
    4.1 0 8
    \$EndMeshFormat
    \$Entities
    0 1 0 0
    5 0 0 0 1 0 0 0 0
    \$EndEntities
    \$Nodes
    1 2 10 20
    1 5 1 2
    10
    20
    0 0 0 0
    1 0 0 1
    \$EndNodes
    \$Elements
    1 1 1 1
    1 5 1 1
    1 10 20
    \$EndElements
    """
    path=joinpath(directory,"parametric-nodes.msh"); write(path,parametric_nodes)
    parametric_mesh=ElementsUnderTest.read_mixed_msh(path)
    @test parametric_mesh.coords==Float64[0 1;0 0;0 0]
    @test parametric_mesh.blocks[1].nodes==reshape(Int32[1,2],2,1)
    @test parametric_mesh.entity_data.external_node_tags==[10,20]
    @test parametric_mesh.entity_data.node_entities==fill((1,Int32(5)),2)
    @test parametric_mesh.entity_data.node_parametric==[Float64[0],Float64[1]]
    @test parametric_mesh.entity_data.external_element_tags==[[1]]
    parametric_output=joinpath(directory,"parametric-nodes-output.msh")
    ElementsUnderTest.write_mixed_msh(parametric_output,parametric_mesh;version=4.1)
    parametric_roundtrip=ElementsUnderTest.read_mixed_msh(parametric_output)
    @test parametric_roundtrip.entity_data.node_parametric==[Float64[0],Float64[1]]
    @test ElementsUnderTest.mixed_crc(parametric_roundtrip).sha==
          ElementsUnderTest.mixed_crc(parametric_mesh).sha
    ok,output=gmsh_check(parametric_output)
    @test ok
    @test !occursin("Error",output)
    parametric_binary=joinpath(directory,"parametric-nodes-output-binary.msh")
    ElementsUnderTest.write_mixed_msh(
        parametric_binary,parametric_mesh;version=4.1,binary=true)
    parametric_binary_roundtrip=ElementsUnderTest.read_mixed_msh(parametric_binary)
    @test parametric_binary_roundtrip.entity_data.node_parametric==
          [Float64[0],Float64[1]]
    @test ElementsUnderTest.mixed_crc(parametric_binary_roundtrip).sha==
          ElementsUnderTest.mixed_crc(parametric_mesh).sha
    ok,output=gmsh_check(parametric_binary)
    @test ok
    @test !occursin("Error",output)
end

_write_swapped(io,value::Int32)=write(io,bswap(value))
_write_swapped(io,value::UInt64)=write(io,bswap(value))
_write_swapped(io,value::Float64)=write(io,bswap(reinterpret(UInt64,value)))

function write_swapped_binary_v2(path)
    open(path,"w") do io
        println(io,"\$MeshFormat"); println(io,"2.2 1 8")
        _write_swapped(io,Int32(1)); write(io,UInt8('\n'))
        println(io,"\$EndMeshFormat")
        println(io,"\$PhysicalNames\n1\n1 6 \"curve\"\n\$EndPhysicalNames")
        println(io,"\$Nodes\n2")
        _write_swapped(io,Int32(10))
        for value in (0.0,0.0,0.0); _write_swapped(io,value); end
        _write_swapped(io,Int32(20))
        for value in (1.0,0.0,0.0); _write_swapped(io,value); end
        write(io,UInt8('\n')); println(io,"\$EndNodes")
        println(io,"\$Elements\n1")
        for value in Int32[1,1,2,77,6,11,10,20]
            _write_swapped(io,value)
        end
        write(io,UInt8('\n')); println(io,"\$EndElements")
    end
    return path
end

function write_swapped_binary_v4(path)
    open(path,"w") do io
        println(io,"\$MeshFormat"); println(io,"4.1 1 8")
        _write_swapped(io,Int32(1)); write(io,UInt8('\n'))
        println(io,"\$EndMeshFormat")
        println(io,"\$PhysicalNames\n1\n1 6 \"curve\"\n\$EndPhysicalNames")
        println(io,"\$Entities")
        for count in UInt64[2,1,0,0]; _write_swapped(io,count); end
        _write_swapped(io,Int32(1))
        for value in (0.0,0.0,0.0); _write_swapped(io,value); end
        _write_swapped(io,UInt64(0))
        _write_swapped(io,Int32(2))
        for value in (1.0,0.0,0.0); _write_swapped(io,value); end
        _write_swapped(io,UInt64(0))
        _write_swapped(io,Int32(11))
        for value in (0.0,0.0,0.0,1.0,0.0,0.0); _write_swapped(io,value); end
        _write_swapped(io,UInt64(1)); _write_swapped(io,Int32(6))
        _write_swapped(io,UInt64(2)); _write_swapped(io,Int32(1))
        _write_swapped(io,Int32(-2))
        write(io,UInt8('\n')); println(io,"\$EndEntities")
        println(io,"\$Nodes")
        for value in UInt64[1,2,10,20]; _write_swapped(io,value); end
        for value in Int32[1,11,1]; _write_swapped(io,value); end
        _write_swapped(io,UInt64(2))
        for value in UInt64[10,20]; _write_swapped(io,value); end
        for value in (0.0,0.0,0.0,0.0,1.0,0.0,0.0,1.0)
            _write_swapped(io,value)
        end
        write(io,UInt8('\n')); println(io,"\$EndNodes")
        println(io,"\$Elements")
        for value in UInt64[1,1,77,77]; _write_swapped(io,value); end
        for value in Int32[1,11,1]; _write_swapped(io,value); end
        _write_swapped(io,UInt64(1))
        for value in UInt64[77,10,20]; _write_swapped(io,value); end
        write(io,UInt8('\n')); println(io,"\$EndElements")
    end
    return path
end

function write_swapped_binary_v2_special(path)
    open(path,"w") do io
        println(io,"\$MeshFormat"); println(io,"2.2 1 8")
        _write_swapped(io,Int32(1)); write(io,UInt8('\n'))
        println(io,"\$EndMeshFormat")
        println(io,"\$PhysicalNames\n1\n1 6 \"curve\"\n\$EndPhysicalNames")
        println(io,"\$Nodes\n2")
        _write_swapped(io,Int32(10))
        for value in (0.0,0.0,0.0); _write_swapped(io,value); end
        _write_swapped(io,Int32(20))
        for value in (1.0,0.0,0.0); _write_swapped(io,value); end
        write(io,UInt8('\n')); println(io,"\$EndNodes")
        println(io,"\$Elements\n2")
        for value in Int32[1,1,2,77,6,11,10,20]
            _write_swapped(io,value)
        end
        for value in Int32[134,1,3,88,6,11,77,10,20]
            _write_swapped(io,value)
        end
        write(io,UInt8('\n')); println(io,"\$EndElements")
    end
    return path
end

function write_swapped_binary_v4_special(path)
    open(path,"w") do io
        println(io,"\$MeshFormat"); println(io,"4.1 1 8")
        _write_swapped(io,Int32(1)); write(io,UInt8('\n'))
        println(io,"\$EndMeshFormat")
        println(io,"\$PhysicalNames\n1\n2 6 \"surface\"\n\$EndPhysicalNames")
        println(io,"\$Entities")
        for count in UInt64[0,0,1,0]; _write_swapped(io,count); end
        _write_swapped(io,Int32(11))
        for value in (0.0,0.0,0.0,1.0,1.0,0.0); _write_swapped(io,value); end
        _write_swapped(io,UInt64(1)); _write_swapped(io,Int32(6))
        _write_swapped(io,UInt64(0))
        write(io,UInt8('\n')); println(io,"\$EndEntities")
        println(io,"\$Nodes")
        for value in UInt64[1,3,10,30]; _write_swapped(io,value); end
        for value in Int32[2,11,0]; _write_swapped(io,value); end
        _write_swapped(io,UInt64(3))
        for value in UInt64[10,20,30]; _write_swapped(io,value); end
        for value in (0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0)
            _write_swapped(io,value)
        end
        write(io,UInt8('\n')); println(io,"\$EndNodes")
        println(io,"\$Elements")
        for value in UInt64[1,1,77,77]; _write_swapped(io,value); end
        for value in Int32[2,11,135]; _write_swapped(io,value); end
        _write_swapped(io,UInt64(1))
        for value in UInt64[77,10,20,30]; _write_swapped(io,value); end
        write(io,UInt8('\n')); println(io,"\$EndElements")
    end
    return path
end

@testset "binary MSH opposite-endian decoding" begin
    directory=mktempdir()
    for (version,writer) in ((2.2,write_swapped_binary_v2),
                             (4.1,write_swapped_binary_v4))
        path=writer(joinpath(directory,"swapped-$version.msh"))
        open(path,"r") do io
            @test readline(io)=="\$MeshFormat"
            @test split(readline(io))==[string(version),"1","8"]
            marker=read(io,UInt32)
            @test marker!=UInt32(1)
            @test bswap(marker)==UInt32(1)
        end
        mesh=ElementsUnderTest.read_mixed_msh(path)
        @test mesh.coords==Float64[0 1;0 0;0 0]
        @test mesh.blocks[1].msh==1
        @test mesh.blocks[1].nodes==reshape(Int32[1,2],2,1)
        @test mesh.blocks[1].tags==Int32[6]
        @test mesh.physical_names==Dict((1,6)=>"curve")
        if version==4.1
            @test mesh.entity_data.external_node_tags==UInt64[10,20]
            @test mesh.entity_data.external_element_tags==[UInt64[77]]
            @test mesh.entity_data.node_parametric==[Float64[0],Float64[1]]
            @test mesh.entity_data.entities[(1,11)].boundaries==Int32[1,-2]
        end
        ok,output=gmsh_check(path)
        @test ok
        @test !occursin("Error",output)
    end

    v2_special=write_swapped_binary_v2_special(
        joinpath(directory,"swapped-v2-special.msh"))
    mesh=ElementsUnderTest.read_mixed_msh(v2_special)
    @test [block.msh for block in mesh.blocks]==[1,134]
    @test mesh.blocks[2].connectivity==Int32[1,2]
    @test mesh.blocks[2].parent_refs==[ElementsUnderTest.ElementRef(1,1)]
    @test mesh.physical_names==Dict((1,6)=>"curve")
    @test ElementsUnderTest.validate(mesh).ok
    ok,output=gmsh_check(v2_special)
    @test ok
    @test !occursin("Error",output)

    v4_special=write_swapped_binary_v4_special(
        joinpath(directory,"swapped-v4-special.msh"))
    mesh=ElementsUnderTest.read_mixed_msh(v4_special)
    @test length(mesh.blocks)==1
    @test mesh.blocks[1].msh==135
    @test mesh.blocks[1].connectivity==Int32[1,2,3]
    @test mesh.blocks[1].parent_refs==[ElementsUnderTest.ElementRef()]
    @test mesh.entity_data.external_node_tags==UInt64[10,20,30]
    @test mesh.entity_data.external_element_tags==[UInt64[77]]
    @test mesh.entity_data.block_entities==[Int32[11]]
    @test ElementsUnderTest.validate(mesh).ok
    ok,output=gmsh_check(v4_special)
    @test ok
    @test !occursin("Error",output)
end

@testset "mixed MSH malformed input, resource gates, and atomic output" begin
    directory=mktempdir()
    function malformed(name,text)
        path=joinpath(directory,name); write(path,text); return path
    end
    cases=Dict(
        "missing-format.msh" => "\$Nodes\n0\n\$EndNodes\n",
        "bad-version.msh" =>
            "\$MeshFormat\n3.0 0 8\n\$EndMeshFormat\n",
        "bad-data-size.msh" =>
            "\$MeshFormat\n4.1 0 4\n\$EndMeshFormat\n",
        "bad-file-type.msh" =>
            "\$MeshFormat\n4.1 2 8\n\$EndMeshFormat\n",
        "binary.msh" =>
            "\$MeshFormat\n4.1 1 8\n\$EndMeshFormat\n",
        "duplicate-node.msh" =>
            "\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n\$Nodes\n2\n1 0 0 0\n1 1 0 0\n\$EndNodes\n",
        "nonfinite-node.msh" =>
            "\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n\$Nodes\n1\n1 NaN 0 0\n\$EndNodes\n",
        "unknown-node.msh" =>
            "\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n\$Nodes\n2\n1 0 0 0\n2 1 0 0\n\$EndNodes\n\$Elements\n1\n1 1 0 1 9\n\$EndElements\n",
        "duplicate-element.msh" =>
            "\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n\$Nodes\n2\n1 0 0 0\n2 1 0 0\n\$EndNodes\n\$Elements\n2\n1 1 0 1 2\n1 1 0 1 2\n\$EndElements\n",
        "special-element.msh" =>
            "\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n\$Nodes\n1\n1 0 0 0\n\$EndNodes\n\$Elements\n1\n1 34 0 1\n\$EndElements\n",
        "truncated-variable-element.msh" =>
            "\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n\$Nodes\n4\n1 0 0 0\n2 1 0 0\n3 0 1 0\n4 1 1 0\n\$EndNodes\n\$Elements\n1\n1 34 0 6 1 2 3 1 3\n\$EndElements\n",
        "extra-variable-element-token.msh" =>
            "\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n\$Nodes\n3\n1 0 0 0\n2 1 0 0\n3 0 1 0\n\$EndNodes\n\$Elements\n1\n1 34 0 3 1 2 3 99\n\$EndElements\n",
        "huge-variable-element.msh" =>
            "\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n\$Nodes\n0\n\$EndNodes\n\$Elements\n1\n1 34 0 2147483646\n\$EndElements\n",
        "dangling-special-parent.msh" =>
            "\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n\$Nodes\n2\n1 0 0 0\n2 1 0 0\n\$EndNodes\n\$Elements\n1\n1 134 5 0 1 1 0 99 1 2\n\$EndElements\n",
        "cyclic-special-parents.msh" =>
            "\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n\$Nodes\n2\n1 0 0 0\n2 1 0 0\n\$EndNodes\n\$Elements\n2\n1 134 5 0 1 1 0 2 1 2\n2 134 5 0 1 1 0 1 1 2\n\$EndElements\n",
        "ambiguous-special-metadata.msh" =>
            "\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n\$Nodes\n2\n1 0 0 0\n2 1 0 0\n\$EndNodes\n\$Elements\n1\n1 134 3 0 1 7 1 2\n\$EndElements\n",
        "truncated-special-partitions.msh" =>
            "\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n\$Nodes\n2\n1 0 0 0\n2 1 0 0\n\$EndNodes\n\$Elements\n1\n1 134 4 0 1 2 0 1 2\n\$EndElements\n",
        "one-domain-special-metadata.msh" =>
            "\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n\$Nodes\n2\n1 0 0 0\n2 1 0 0\n\$EndNodes\n\$Elements\n1\n1 67 5 0 1 1 0 9 1 2\n\$EndElements\n",
        "basis-only-element.msh" =>
            "\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n\$Nodes\n4\n1 0 0 0\n2 1 0 0\n3 0 1 0\n4 0.25 0.25 0\n\$EndNodes\n\$Elements\n1\n1 138 0 1 2 3 4\n\$EndElements\n",
        "v4-variable-element.msh" =>
            "\$MeshFormat\n4.1 0 8\n\$EndMeshFormat\n\$Entities\n0 0 1 0\n1 0 0 0 1 1 0 0 0\n\$EndEntities\n\$Nodes\n1 3 1 3\n2 1 0 3\n1 2 3\n0 0 0\n1 0 0\n0 1 0\n\$EndNodes\n\$Elements\n1 1 1 1\n2 1 34 1\n1 1 2 3\n\$EndElements\n",
        "truncated-section.msh" =>
            "\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n\$Comments\nmissing end\n",
        "huge-truncated-count.msh" =>
            "\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n\$Nodes\n2147483647\n",
        "repeated-node-sections.msh" =>
            "\$MeshFormat\n4.1 0 8\n\$EndMeshFormat\n" *
            "\$Nodes\n0 0 0 0\n\$EndNodes\n" *
            "\$Nodes\n0 0 0 0\n\$EndNodes\n",
        "bad-boundary.msh" =>
            "\$MeshFormat\n4.1 0 8\n\$EndMeshFormat\n\$Entities\n0 1 0 0\n1 0 0 0 1 0 0 0 1 999\n\$EndEntities\n",
        "bad-v4-range.msh" =>
            "\$MeshFormat\n4.1 0 8\n\$EndMeshFormat\n\$Entities\n0 0 0 1\n1 0 0 0 0 0 0 0 0\n\$EndEntities\n\$Nodes\n1 1 1 2\n3 1 0 1\n2\n0 0 0\n\$EndNodes\n",
        "bad-v4-dimension.msh" =>
            "\$MeshFormat\n4.1 0 8\n\$EndMeshFormat\n\$Entities\n0 1 0 1\n1 0 0 0 1 0 0 0 0\n1 0 0 0 1 0 0 0 0\n\$EndEntities\n\$Nodes\n1 2 1 2\n3 1 0 2\n1\n2\n0 0 0\n1 0 0\n\$EndNodes\n\$Elements\n1 1 1 1\n1 1 2 1\n1 1 2 1\n\$EndElements\n",
    )
    for (name,text) in cases
        @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(
            malformed(name,text))
    end

    function malformed_binary(writer,name)
        path=joinpath(directory,name)
        open(path,"w") do io
            writer(io)
        end
        return path
    end
    binary_format(io,version="4.1")=begin
        println(io,"\$MeshFormat"); println(io,version," 1 8")
        write(io,Int32(1)); write(io,UInt8('\n')); println(io,"\$EndMeshFormat")
    end
    invalid_marker=malformed_binary("invalid-marker.msh") do io
        println(io,"\$MeshFormat\n4.1 1 8")
        write(io,Int32(2)); write(io,UInt8('\n')); println(io,"\$EndMeshFormat")
    end
    unknown_binary=malformed_binary("unknown-binary-section.msh") do io
        binary_format(io)
        println(io,"\$Comments\nnot safely skippable\n\$EndComments")
    end
    truncated_v2_nodes=malformed_binary("truncated-v2-nodes.msh") do io
        binary_format(io,"2.2")
        println(io,"\$Nodes\n1")
        write(io,Int32(1))
    end
    huge_v2_element=malformed_binary("huge-v2-element-width.msh") do io
        binary_format(io,"2.2")
        println(io,"\$Nodes\n0"); write(io,UInt8('\n')); println(io,"\$EndNodes")
        println(io,"\$Elements\n1")
        write(io,Int32(1)); write(io,Int32(1)); write(io,typemax(Int32))
    end
    binary_v2_variable=malformed_binary("binary-v2-variable.msh") do io
        binary_format(io,"2.2")
        println(io,"\$Nodes\n0"); write(io,UInt8('\n')); println(io,"\$EndNodes")
        println(io,"\$Elements\n1")
        write(io,Int32(34)); write(io,Int32(1)); write(io,Int32(0))
    end
    truncated_v2_special=malformed_binary("truncated-v2-special.msh") do io
        binary_format(io,"2.2")
        println(io,"\$Nodes\n2")
        write(io,Int32(1)); for value in (0.0,0.0,0.0); write(io,value); end
        write(io,Int32(2)); for value in (1.0,0.0,0.0); write(io,value); end
        write(io,UInt8('\n')); println(io,"\$EndNodes")
        println(io,"\$Elements\n1")
        write(io,Int32(134)); write(io,Int32(1)); write(io,Int32(2))
        for value in Int32[1,0,1,1]; write(io,value); end
    end
    binary_v2_domain_link=malformed_binary("binary-v2-domain-link.msh") do io
        binary_format(io,"2.2")
        println(io,"\$Nodes\n2")
        write(io,Int32(1)); for value in (0.0,0.0,0.0); write(io,value); end
        write(io,Int32(2)); for value in (1.0,0.0,0.0); write(io,value); end
        write(io,UInt8('\n')); println(io,"\$EndNodes")
        println(io,"\$Elements\n1")
        for value in Int32[67,1,3,1,0,1,9,1,2]; write(io,value); end
    end
    huge_v4_entities=malformed_binary("huge-v4-entity-membership.msh") do io
        binary_format(io)
        println(io,"\$Entities")
        for count in UInt64[1,0,0,0]; write(io,count); end
        write(io,Int32(1)); for value in (0.0,0.0,0.0); write(io,value); end
        write(io,typemax(UInt64))
    end
    huge_v4_nodes=malformed_binary("huge-v4-node-count.msh") do io
        binary_format(io)
        println(io,"\$Nodes")
        for value in UInt64[1,typemax(UInt64),1,typemax(UInt64)]
            write(io,value)
        end
    end
    bad_binary_newline=malformed_binary("bad-binary-newline.msh") do io
        binary_format(io,"2.2")
        println(io,"\$Nodes\n0"); write(io,UInt8('x')); println(io,"\$EndNodes")
    end
    for path in (invalid_marker,unknown_binary,truncated_v2_nodes,
                 huge_v2_element,binary_v2_variable,truncated_v2_special,
                 binary_v2_domain_link,huge_v4_entities,huge_v4_nodes,
                 bad_binary_newline)
        @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(path)
    end

    valid=joinpath(directory,"valid.msh")
    ElementsUnderTest.write_mixed_msh(
        valid,mixed_io_fixture(escaped_name=false);version=4.1)
    @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(valid;max_nodes=11)
    @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(valid;max_elements=4)
    @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(valid;max_blocks=0)
    @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(valid;max_entities=5)
    @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(valid;max_physical_names=4)
    @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(valid;max_name_bytes=3)
    @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(
        valid;max_file_bytes=filesize(valid)-1)
    @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(valid;max_nodes=-1)
    @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(valid;max_nodes=1.0)
    @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(
        valid;max_connectivity=-1)
    @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(
        valid;max_connectivity=1.0)
    @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(
        valid;tessella_extensions=1)
    valid_binary=joinpath(directory,"valid-binary.msh")
    ElementsUnderTest.write_mixed_msh(
        valid_binary,mixed_io_fixture(escaped_name=false);
        version=4.1,binary=true)
    @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(
        valid_binary;max_nodes=11)
    @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(
        valid_binary;max_elements=4)
    @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(
        valid_binary;max_blocks=0)
    @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(
        valid_binary;max_entities=5)
    @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(
        valid_binary;max_physical_names=4)
    @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(
        valid_binary;max_name_bytes=3)
    @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(
        valid_binary;max_file_bytes=filesize(valid_binary)-1)

    special_valid=joinpath(directory,"special-valid-ascii.msh")
    ElementsUnderTest.write_mixed_msh(
        special_valid,special_v2_fixture();version=2.2,
        gmsh_compatible=false)
    @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(
        special_valid;max_connectivity=54)
    @test ElementsUnderTest.mixed_crc(
        ElementsUnderTest.read_mixed_msh(
            special_valid;max_connectivity=55)).sha==
          "944dbe417200e316bf265d2d6f8a154dfb4c6bec3eb853401f289990a9be4089"
    special_binary=joinpath(directory,"special-valid-binary.msh")
    ElementsUnderTest.write_mixed_msh(
        special_binary,special_v2_binary_fixture();version=2.2,binary=true)
    @test_throws ArgumentError ElementsUnderTest.read_mixed_msh(
        special_binary;max_connectivity=34)
    @test ElementsUnderTest.mixed_crc(
        ElementsUnderTest.read_mixed_msh(
            special_binary;max_connectivity=35)).sha==
          "4ef8de3b81b4449aa1b3efc3905f4418e00b5ac96ffd609e6eef2957b1156b68"

    function write_rejected_variable_record(path,arity)
        open(path,"w") do io
            println(io,"\$MeshFormat\n2.2 0 8\n\$EndMeshFormat")
            println(io,"\$Nodes\n3")
            println(io,"1 0 0 0\n2 1 0 0\n3 0 1 0")
            println(io,"\$EndNodes\n\$Elements\n1")
            print(io,"1 34 0 ",arity)
            for _ in 1:arity
                print(io," 1")
            end
            println(io,"\n\$EndElements")
        end
        return path
    end
    function reject_zero_connectivity(path)
        try
            ElementsUnderTest.read_mixed_msh(path;max_connectivity=0)
        catch err
            err isa ArgumentError || rethrow()
            return nothing
        end
        error("oversized variable record unexpectedly passed max_connectivity=0")
    end
    rejected_small=write_rejected_variable_record(
        joinpath(directory,"rejected-variable-small.msh"),30_000)
    rejected_large=write_rejected_variable_record(
        joinpath(directory,"rejected-variable-large.msh"),300_000)
    reject_zero_connectivity(rejected_small)
    GC.gc(); rejected_small_bytes=@allocated reject_zero_connectivity(rejected_small)
    GC.gc(); rejected_large_bytes=@allocated reject_zero_connectivity(rejected_large)
    @test rejected_large_bytes<=rejected_small_bytes+32_768

    target=joinpath(directory,"atomic.msh"); write(target,"sentinel")
    invalid=mixed_io_fixture(escaped_name=false); invalid.blocks[1].tags[1]=Int32(-1)
    @test_throws ArgumentError ElementsUnderTest.write_mixed_msh(target,invalid)
    @test read(target,String)=="sentinel"
    @test_throws ArgumentError ElementsUnderTest.write_mixed_msh(
        target,mixed_io_fixture(escaped_name=false);version=3.0)
    @test read(target,String)=="sentinel"
    @test_throws ArgumentError ElementsUnderTest.write_mixed_msh(
        target,mixed_io_fixture(escaped_name=false);gmsh_compatible=1)
    @test read(target,String)=="sentinel"
    @test_throws ArgumentError ElementsUnderTest.write_mixed_msh(
        target,mixed_io_fixture(escaped_name=false);binary=1)
    @test read(target,String)=="sentinel"
    @test_throws ArgumentError ElementsUnderTest.write_mixed_msh(
        target,invalid;binary=true)
    @test read(target,String)=="sentinel"
    long_name=mixed_io_fixture(escaped_name=false)
    long_name.physical_names[(3,7)]=repeat("x",129)
    for binary in (false,true)
        @test_throws ArgumentError ElementsUnderTest.write_mixed_msh(
            target,long_name;binary=binary)
        @test read(target,String)=="sentinel"
    end
end

@testset "binary MSH read allocation grows linearly" begin
    function line_chain(n)
        coordinates=zeros(3,n); coordinates[1,:].=0:n-1
        connectivity=Matrix{Int32}(undef,2,n-1)
        @inbounds for j in 1:n-1
            connectivity[1,j]=j; connectivity[2,j]=j+1
        end
        block=ElementsUnderTest.ElementBlock(
            1,connectivity,fill(Int32(3),n-1))
        return ElementsUnderTest.MixedMesh(coordinates,[block])
    end
    directory=mktempdir()
    small=joinpath(directory,"small-binary.msh")
    large=joinpath(directory,"large-binary.msh")
    ElementsUnderTest.write_mixed_msh(small,line_chain(2_000);binary=true)
    ElementsUnderTest.write_mixed_msh(large,line_chain(4_000);binary=true)
    small_mesh=ElementsUnderTest.read_mixed_msh(small)
    large_mesh=ElementsUnderTest.read_mixed_msh(large)
    @test ElementsUnderTest.validate(small_mesh).ok
    @test ElementsUnderTest.validate(large_mesh).ok
    @test size(small_mesh.coords,2)==2_000
    @test size(large_mesh.coords,2)==4_000
    @test filesize(large)<=2.1filesize(small)
    GC.gc(); allocated_small=@allocated ElementsUnderTest.read_mixed_msh(small)
    GC.gc(); allocated_large=@allocated ElementsUnderTest.read_mixed_msh(large)
    @test allocated_large<=2.8allocated_small+131_072
end

@testset "variable-connectivity ASCII read allocation grows linearly" begin
    function write_triangle_decompositions(path,count)
        open(path,"w") do io
            println(io,"\$MeshFormat\n2.2 0 8\n\$EndMeshFormat")
            println(io,"\$Nodes\n",count+2)
            for i in 1:count+2
                println(io,i," ",i-1," 0 0")
            end
            println(io,"\$EndNodes\n\$Elements\n",count)
            for i in 1:count
                println(io,i," 34 0 3 ",i," ",i+1," ",i+2)
            end
            println(io,"\$EndElements")
        end
        return path
    end
    directory=mktempdir()
    small=write_triangle_decompositions(
        joinpath(directory,"small-variable.msh"),2_000)
    large=write_triangle_decompositions(
        joinpath(directory,"large-variable.msh"),4_000)
    small_mesh=ElementsUnderTest.read_mixed_msh(small)
    large_mesh=ElementsUnderTest.read_mixed_msh(large)
    @test ElementsUnderTest.validate(small_mesh).ok
    @test ElementsUnderTest.validate(large_mesh).ok
    @test length(small_mesh.blocks[1].offsets)==2_001
    @test length(large_mesh.blocks[1].offsets)==4_001
    @test filesize(large)<=2.2filesize(small)
    GC.gc(); allocated_small=@allocated ElementsUnderTest.read_mixed_msh(small)
    GC.gc(); allocated_large=@allocated ElementsUnderTest.read_mixed_msh(large)
    @test allocated_large<=2.8allocated_small+131_072
end
