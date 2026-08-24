#!/usr/bin/env julia
# Differential validation of Tessella's size-field catalog against Gmsh 4.15.2.
#
# Gmsh's public field API does not expose a coordinate evaluator. The official
# MeshSizeFieldView plugin does: it writes a field onto model-backed NodeData.
# The direct checks below use that supported path. Mesh grading remains an
# intentionally approximate check, and model/algorithm-dependent gaps are
# reported as CONTEXT_SKIP entries instead of being presented as parity.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella
using Tessella.Mesh1D: mesh_segment
using Printf: @sprintf
using Statistics: median

const TARGET_GMSH_VERSION = "4.15.2"
const DIRECT_RESULTS = NamedTuple[]
const MESH_RESULTS = NamedTuple[]
const PLUGIN_CALLS = Ref(0)
const PROBE_SERIAL = Ref(1_000_000)

const CONTEXT_SKIPS = [
    (name="anisotropic metric tensors",
     reason="the public Gmsh API/plugin exposes each field's scalar operator, not its SMetric3; MathEvalAniso, MinAniso, IntersectAniso, and AttractorAnisoCurve metric entries are therefore not claimed pointwise"),
    (name="anisotropic field-driven meshing",
     reason="algorithm- and dimension-dependent metric consumption is not a scalar oracle (Gmsh 4.15.2 scalarizes the tested 1-D MathEvalAniso background)"),
    (name="boundary-layer element construction",
     reason="the scalar BoundaryLayer law is checked, but fan/quads/extrusion topology belongs to the boundary-layer meshing track"),
    (name="PostView differential scope",
     reason="closest scalar-point views are checked directly; first-order scalar/vector standard list elements, multiple-time-step selection and the tensor scalar operator are separately tested but are not exercised here; high-order/custom-interpolation, mixed-component and tensor-metric views remain unsupported"),
    (name="AutomaticMeshSizeField",
     reason="Tessella exposes a documented discrete sphere-fit analogue while Gmsh uses a global HXT/P4EST model pipeline; no equivalent input/state oracle exists"),
]

function find_gmsh_executable()
    explicit = get(ENV, "GMSH_EXECUTABLE", "")
    if !isempty(explicit)
        isfile(explicit) || error("GMSH_EXECUTABLE does not name a file: $explicit")
        return realpath(explicit)
    end
    executable = Sys.which("gmsh")
    executable !== nothing && return realpath(executable)
    fallback = "/opt/homebrew/bin/gmsh"
    isfile(fallback) && return realpath(fallback)
    error("Gmsh executable is required; install 4.15.2 or set GMSH_EXECUTABLE")
end

function find_gmsh_api(executable)
    explicit = get(ENV, "GMSH_JULIA_API", "")
    if !isempty(explicit)
        isfile(explicit) || error("GMSH_JULIA_API does not name a file: $explicit")
        return realpath(explicit)
    end
    prefix = dirname(dirname(executable))
    candidates = (joinpath(prefix, "lib", "gmsh.jl"),
                  joinpath(prefix, "lib64", "gmsh.jl"),
                  "/opt/homebrew/lib/gmsh.jl",
                  "/opt/homebrew/opt/gmsh/lib/gmsh.jl",
                  "/usr/local/opt/gmsh/lib/gmsh.jl")
    for path in candidates
        isfile(path) && return realpath(path)
    end
    error("could not locate gmsh.jl for $executable; set GMSH_JULIA_API")
end

const GMSH_EXECUTABLE = find_gmsh_executable()
const GMSH_API_FILE = find_gmsh_api(GMSH_EXECUTABLE)
const GMSH_CLI_VERSION = strip(read(`$GMSH_EXECUTABLE --version`, String))
(GMSH_CLI_VERSION == TARGET_GMSH_VERSION ||
 startswith(GMSH_CLI_VERSION, TARGET_GMSH_VERSION * "-")) || error(
    "expected Gmsh $TARGET_GMSH_VERSION, got CLI $GMSH_CLI_VERSION from $GMSH_EXECUTABLE")
include(GMSH_API_FILE)

function setup_model(name)
    gmsh.clear()
    gmsh.model.add(name)
end

function add_field(kind, tag; numbers=(), strings=(), lists=())
    actual = gmsh.model.mesh.field.add(kind, tag)
    actual == tag || error("Gmsh allocated Field[$actual], expected Field[$tag]")
    for (option, value) in numbers
        gmsh.model.mesh.field.setNumber(tag, option, value)
    end
    for (option, value) in strings
        gmsh.model.mesh.field.setString(tag, option, value)
    end
    for (option, value) in lists
        gmsh.model.mesh.field.setNumbers(tag, option, value)
    end
    return tag
end

function add_probe_nodes(points; entity_tags=nothing)
    isempty(points) && error("cannot create an empty Gmsh field probe")
    entity_tags === nothing || length(entity_tags) == length(points) || error(
        "entity tag count does not match probe point count")
    tags = Vector{UInt64}(undef, length(points))
    for (i, point) in enumerate(points)
        length(point) == 3 || error("probe point $i is not three-dimensional")
        PROBE_SERIAL[] += 1
        entity = entity_tags === nothing ? PROBE_SERIAL[] : Int(entity_tags[i])
        node = UInt64(PROBE_SERIAL[] + 1_000_000_000)
        gmsh.model.addDiscreteEntity(0, entity)
        gmsh.model.mesh.addNodes(0, entity, UInt64[node], Float64[point...])
        gmsh.model.mesh.addElementsByType(entity, 15, UInt64[node], UInt64[node])
        tags[i] = node
    end
    return tags
end

function gmsh_probe_nodes(field_tag, node_tags)
    isempty(node_tags) && error("cannot evaluate a Gmsh field on zero nodes")
    view = gmsh.view.add("field_probe_$(PLUGIN_CALLS[] + 1)")
    gmsh.view.addModelData(view, 0, gmsh.model.getCurrent(), "NodeData", node_tags,
                           [[0.0] for _ in node_tags], 0.0, 1)
    view_tags = gmsh.view.getTags()
    position = findfirst(==(view), view_tags)
    position === nothing && error("Gmsh did not register probe View[$view]")
    gmsh.plugin.setNumber("MeshSizeFieldView", "MeshSizeField", field_tag)
    gmsh.plugin.setNumber("MeshSizeFieldView", "View", position - 1) # plugin uses index
    gmsh.plugin.setNumber("MeshSizeFieldView", "Component", 0)
    returned = gmsh.plugin.run("MeshSizeFieldView")
    returned == view || error("MeshSizeFieldView returned View[$returned], expected View[$view]")
    data_type, output_tags, data, _, components = gmsh.view.getModelData(view, 0)
    data_type == "NodeData" || error("field probe returned $data_type, expected NodeData")
    components == 1 || error("field probe returned $components components, expected 1")
    length(output_tags) == length(data) || error("field probe tag/data length mismatch")
    values = Dict{UInt64,Float64}()
    for (tag, datum) in zip(output_tags, data)
        length(datum) == 1 || error("field probe node $tag has $(length(datum)) values")
        values[UInt64(tag)] = only(datum)
    end
    result = Float64[get(values, tag) do
        error("field probe omitted node $tag")
    end for tag in node_tags]
    all(isfinite, result) || error("Gmsh field probe returned non-finite values: $result")
    PLUGIN_CALLS[] += 1
    return result
end

gmsh_probe(field_tag, points; entity_tags=nothing) =
    gmsh_probe_nodes(field_tag, add_probe_nodes(points; entity_tags=entity_tags))

function tessella_values(field, points; entities=nothing)
    if entities === nothing
        return Float64[field_value(field, point...) for point in points]
    end
    length(entities) == length(points) || error("entity count does not match point count")
    return Float64[field_value(field, point..., entity)
                   for (point, entity) in zip(points, entities)]
end

function check_direct(name, gmsh_values, tessella_output; atol=2e-12, rtol=2e-12,
                      classification="direct")
    length(gmsh_values) == length(tessella_output) || error(
        "$name: result lengths differ (Gmsh=$(length(gmsh_values)), Tessella=$(length(tessella_output)))")
    isempty(gmsh_values) && error("$name: no samples")
    max_abs = 0.0
    max_rel = 0.0
    for i in eachindex(gmsh_values, tessella_output)
        g = gmsh_values[i]
        t = tessella_output[i]
        all(isfinite, (g, t)) || error("$name sample $i is non-finite: Gmsh=$g Tessella=$t")
        delta = abs(g - t)
        allowed = atol + rtol * max(abs(g), abs(t))
        delta <= allowed || error(
            "$name sample $i mismatch: Gmsh=$g Tessella=$t abs=$delta allowed=$allowed")
        max_abs = max(max_abs, delta)
        max_rel = max(max_rel, delta / max(abs(g), abs(t), floatmin(Float64)))
    end
    result = (name=name, classification=classification, samples=length(gmsh_values),
              max_abs=max_abs, max_rel=max_rel)
    push!(DIRECT_RESULTS, result)
    println(@sprintf("  DIRECT %-29s samples=%2d max_abs=%9.3g max_rel=%9.3g",
                     name, result.samples, max_abs, max_rel))
    return result
end

function check_matheval()
    setup_model("matheval")
    expression = "0.31 + 0.07*x - 0.04*y + 0.03*z + Sin(x)/100"
    add_field("MathEval", 1; strings=(("F", expression),))
    # VERIFIED probe constraint: a MathEval -> MathEval F<n> view probe did not
    # return. A Box dependency exercises the F-reference path without relying
    # on that probe combination.
    add_field("Box", 2; numbers=(("VIn", 0.2), ("VOut", 0.4),
                                 ("XMin", -0.5), ("XMax", 0.5),
                                 ("YMin", -1.0), ("YMax", 1.0),
                                 ("ZMin", -1.0), ("ZMax", 1.0)))
    add_field("MathEval", 3; strings=(("F", "F2 + 0.05"),))
    points = [(-0.7, 0.2, 0.4), (0.0, 0.0, 0.0), (0.8, -0.3, 1.1)]
    base = MathEvalField(expression)
    box = BoxField(-0.5, 0.5, -1.0, 1.0, -1.0, 1.0; vin=0.2, vout=0.4)
    composed = MathEvalField("F2 + 0.05"; fields=Dict(2 => box))
    check_direct("MathEval", gmsh_probe(1, points), tessella_values(base, points))
    check_direct("MathEval F-reference", gmsh_probe(3, points),
                 tessella_values(composed, points))
end

function check_differential_operators()
    setup_model("differential_operators")
    expression = "x^2+2*y^2+3*z^2"
    add_field("MathEval", 1; strings=(("F", expression),))
    # A 1e-2 step keeps Gmsh's Hessian eigensolve away from cancellation while
    # the quadratic input makes the centered finite differences exact in theory.
    delta = 1e-2
    specs = (("Gradient", 2, (("InField", 1.0), ("Kind", 3.0), ("Delta", delta))),
             ("Laplacian", 3, (("InField", 1.0), ("Delta", delta))),
             ("Mean", 4, (("InField", 1.0), ("Delta", delta))),
             ("Curvature", 5, (("InField", 1.0), ("Delta", delta))),
             ("MaxEigenHessian", 6, (("InField", 1.0), ("Delta", delta))))
    for (kind, tag, options) in specs
        add_field(kind, tag; numbers=options)
    end
    input = MathEvalField(expression)
    fields = (GradientField(input; kind=3, delta=delta),
              LaplacianField(input; delta=delta), MeanField(input; delta=delta),
              CurvatureField(input; delta=delta),
              MaxEigenHessianField(input; delta=delta))
    points = [(1.0, 2.0, 3.0), (-0.6, 0.8, 1.2)]
    for ((kind, tag, _), field) in zip(specs, fields)
        check_direct(kind, gmsh_probe(tag, points), tessella_values(field, points);
                     atol=8e-7, rtol=2e-8)
    end
end

function check_coordinate_maps()
    setup_model("coordinate_maps")
    expression = "0.4 + 0.1*x + 0.2*y + 0.03*z"
    add_field("MathEval", 1; strings=(("F", expression),))
    add_field("LonLat", 2; numbers=(("InField", 1.0), ("RadiusStereo", 2.0)))
    add_field("LonLat", 3; numbers=(("InField", 1.0), ("FromStereo", 1.0),
                                             ("RadiusStereo", 2.0)))
    add_field("Param", 4; numbers=(("InField", 1.0),),
              strings=(("FX", "2*x"), ("FY", "y+1"), ("FZ", "z/2")))
    input = MathEvalField(expression)
    spherical = [(2.0, 0.0, 0.0), (0.0, sqrt(3.0), 1.0)]
    stereo = [(0.0, 0.0, 0.0), (0.8, -0.4, 0.0)]
    parametric = [(0.1, 0.2, 0.4), (0.3, -0.1, 0.8)]
    check_direct("LonLat", gmsh_probe(2, spherical),
                 tessella_values(LonLatField(input; radius=2.0), spherical))
    check_direct("LonLat stereographic", gmsh_probe(3, stereo),
                 tessella_values(LonLatField(input; from_stereo=true, radius=2.0), stereo))
    check_direct("Param", gmsh_probe(4, parametric),
                 tessella_values(ParametricField(input; fx="2*x", fy="y+1", fz="z/2"),
                                  parametric))
end

function check_structured()
    mktempdir() do directory
        path = joinpath(directory, "structured.txt")
        data = Array{Float64}(undef, 2, 2, 2)
        for i in 1:2, j in 1:2, k in 1:2
            data[i, j, k] = 0.1 + 0.2 * (i - 1) + 0.3 * (j - 1) + 0.4 * (k - 1)
        end
        open(path, "w") do io
            println(io, "0 0 0")
            println(io, "1 1 1")
            println(io, "2 2 2")
            println(io, join((data[i, j, k] for i in 1:2, j in 1:2, k in 1:2), ' '))
        end
        setup_model("structured")
        add_field("Structured", 1;
                  numbers=(("TextFormat", 1.0), ("SetOutsideValue", 1.0),
                           ("OutsideValue", 0.91)),
                  strings=(("FileName", path),))
        field = StructuredField(path; text=true, outside=0.91)
        points = [(0.25, 0.5, 0.75), (0.8, 0.2, 0.4), (-0.1, 0.5, 0.5)]
        check_direct("Structured", gmsh_probe(1, points), tessella_values(field, points))
    end
end

function check_entity_fields()
    setup_model("entity_fields")
    points = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]
    entity_tags = [101, 102]
    nodes = add_probe_nodes(points; entity_tags=entity_tags)
    add_field("MathEval", 1; strings=(("F", "0.37"),))
    add_field("Restrict", 2; numbers=(("InField", 1.0), ("IncludeBoundary", 0.0)),
              lists=(("PointsList", [101]),))
    add_field("Constant", 3; numbers=(("VIn", 0.23), ("VOut", 0.89),
                                               ("IncludeBoundary", 0.0)),
              lists=(("PointsList", [101]),))
    entities = [(0, 101), (0, 102)]
    restrict = RestrictField(ConstantSize(0.37); points=[101], include_boundary=false)
    constant = ConstantField(; vin=0.23, vout=0.89, points=[101], include_boundary=false)
    check_direct("Restrict entity dispatch", gmsh_probe_nodes(2, nodes),
                 tessella_values(restrict, points; entities=entities))
    check_direct("Constant entity dispatch", gmsh_probe_nodes(3, nodes),
                 tessella_values(constant, points; entities=entities))
end

function check_postview()
    setup_model("postview")
    input_view = gmsh.view.add("scalar_points")
    gmsh.view.addListData(input_view, "SP", 3,
        [0.0, 0.0, 0.0, 0.2, 1.0, 0.0, 0.0, 0.6, 2.0, 0.0, 0.0, -0.4])
    add_field("PostView", 1;
              numbers=(("ViewTag", Float64(input_view)), ("UseClosest", 1.0),
                       ("CropNegativeValues", 1.0)))
    coordinates = [0.0 1.0 2.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    field = PostViewField(coordinates, [0.2, 0.6, -0.4];
                          crop_negative=true, use_closest=true)
    points = [(0.1, 0.0, 0.0), (0.8, 0.0, 0.0), (1.8, 0.0, 0.0)]
    check_direct("PostView scalar points", gmsh_probe(1, points),
                 tessella_values(field, points); classification="direct subset")
end

function check_anisotropic_scalar_operators()
    setup_model("anisotropic_scalars")
    diagonal_a = (("M11", "1"), ("M22", "4"), ("M33", "9"),
                  ("M12", "0"), ("M13", "0"), ("M23", "0"))
    diagonal_b = (("M11", "4"), ("M22", "1"), ("M33", "16"),
                  ("M12", "1.2"), ("M13", "0.5"), ("M23", "-0.7"))
    add_field("MathEvalAniso", 1; strings=diagonal_a)
    add_field("MathEvalAniso", 2; strings=diagonal_b)
    add_field("MinAniso", 3; lists=(("FieldsList", [1, 2]),))
    add_field("IntersectAniso", 4; lists=(("FieldsList", [1, 2]),))
    a = MathEvalAnisoField(; m11="1", m22="4", m33="9")
    b = MathEvalAnisoField(; m11="4", m22="1", m33="16",
                           m12="1.2", m13="0.5", m23="-0.7")
    fields = (a, b, MinAnisoField((a, b)), IntersectAnisoField((a, b)))
    names = ("MathEvalAniso M11 (A)", "MathEvalAniso M11 (B)",
             "MinAniso scalar", "IntersectAniso scalar")
    points = [(0.1, 0.2, 0.3), (-0.4, 0.7, 0.2)]
    for (tag, name, field) in zip(1:4, names, fields)
        check_direct(name, gmsh_probe(tag, points), tessella_values(field, points);
                     atol=2e-10, rtol=2e-10, classification="scalar operator")
    end
end

function check_boundary_layer_scalar()
    setup_model("boundary_layer_scalar")
    gmsh.model.geo.addPoint(0.0, 0.0, 0.0, 1.0, 1)
    gmsh.model.geo.synchronize()
    add_field("BoundaryLayer", 1;
              numbers=(("Size", 0.05), ("Ratio", 1.2), ("Thickness", 0.4),
                       ("SizeFar", 0.5)), lists=(("PointsList", [1]),))
    field = BoundaryLayerField(DistanceField(points=[(0.0, 0.0, 0.0)]);
                               hwall=0.05, ratio=1.2, thickness=0.4, hfar=0.5)
    points = [(0.0, 0.0, 0.0), (0.1, 0.0, 0.0),
              (0.45, 0.0, 0.0), (0.6, 0.0, 0.0)]
    check_direct("BoundaryLayer scalar law", gmsh_probe(1, points),
                 tessella_values(field, points); classification="direct scalar")
end

function check_attractor_aniso_scalar()
    setup_model("attractor_aniso_scalar")
    gmsh.model.geo.addPoint(0.0, 0.0, 0.0, 1.0, 1)
    gmsh.model.geo.addPoint(1.0, 0.0, 0.0, 1.0, 2)
    gmsh.model.geo.addLine(1, 2, 1)
    gmsh.model.geo.synchronize()
    sampling = 5
    add_field("AttractorAnisoCurve", 1;
              numbers=(("Sampling", Float64(sampling)), ("DistMin", 0.1),
                       ("DistMax", 0.5), ("SizeMinTangent", 0.2),
                       ("SizeMaxTangent", 0.5), ("SizeMinNormal", 0.05),
                       ("SizeMaxNormal", 0.5)), lists=(("CurvesList", [1]),))
    samples = [(i / (sampling - 1), 0.0, 0.0) for i in 0:sampling-1]
    tangents = fill((1.0, 0.0, 0.0), sampling)
    field = AttractorAnisoCurveField(samples, tangents; dist_min=0.1, dist_max=0.5,
        size_min_tangent=0.2, size_max_tangent=0.5,
        size_min_normal=0.05, size_max_normal=0.5)
    points = [(0.1, 0.1, 0.0), (0.5, 0.2, 0.0), (1.1, 0.0, 0.0)]
    check_direct("AttractorAnisoCurve scalar", gmsh_probe(1, points),
                 tessella_values(field, points); classification="sampled scalar")
end

function check_octree()
    setup_model("octree")
    gmsh.model.occ.addBox(0.0, 0.0, 0.0, 1.0, 1.0, 1.0)
    gmsh.model.occ.synchronize()
    expression = "0.2+0.1*x+0.05*y"
    add_field("MathEval", 1; strings=(("F", expression),))
    add_field("Octree", 2; numbers=(("InField", 1.0),))
    raw_bbox = gmsh.model.getBoundingBox(-1, -1)
    bbox = (raw_bbox[1], raw_bbox[4], raw_bbox[2], raw_bbox[5],
            raw_bbox[3], raw_bbox[6])
    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    # VERIFIED probe constraint: before any mesh pass, Gmsh 4.15.2 crashed in
    # OctreeField::Cell::evaluate. A real mesh pass initializes this fixture.
    gmsh.model.mesh.generate(1)
    field = OctreeField(MathEvalField(expression), bbox...; max_level=4)
    points = [(0.0, 0.0, 0.0), (1.0, 1.0, 1.0),
              (0.23, 0.41, 0.62), (0.77, 0.12, 0.5)]
    check_direct("Octree sampled field", gmsh_probe(2, points),
                 tessella_values(field, points); atol=2e-12, rtol=2e-12,
                 classification="direct approximation")
end

function curve_samples_and_sizes(curve_tag)
    types, _, connectivities = gmsh.model.mesh.getElements(1, curve_tag)
    position = findfirst(==(Int32(1)), types) # first-order line
    position === nothing && error("Extend fixture curve $curve_tag has no type-1 elements")
    connectivity = connectivities[position]
    iseven(length(connectivity)) || error("Extend fixture has malformed line connectivity")
    samples = NTuple{3,Float64}[]
    sizes = Float64[]
    for i in 1:2:length(connectivity)
        a, _, _, _ = gmsh.model.mesh.getNode(connectivity[i])
        b, _, _, _ = gmsh.model.mesh.getNode(connectivity[i + 1])
        push!(samples, ((a[1] + b[1]) / 2, (a[2] + b[2]) / 2, (a[3] + b[3]) / 2))
        push!(sizes, hypot(b[1] - a[1], b[2] - a[2], b[3] - a[3]))
    end
    return samples, sizes
end

function check_extend()
    setup_model("extend")
    for (tag, (x, y)) in enumerate(((0.0, 0.0), (1.0, 0.0),
                                     (1.0, 1.0), (0.0, 1.0)))
        gmsh.model.geo.addPoint(x, y, 0.0, 1.0, tag)
    end
    for (tag, (a, b)) in enumerate(((1, 2), (2, 3), (3, 4), (4, 1)))
        gmsh.model.geo.addLine(a, b, tag)
    end
    gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)
    gmsh.model.geo.addPlaneSurface([1], 1)
    gmsh.model.geo.synchronize()
    for tag in 1:4
        gmsh.model.mesh.setTransfiniteCurve(tag, 5)
    end
    gmsh.model.mesh.setTransfiniteSurface(1, "Left", [1, 2, 3, 4])
    gmsh.model.mesh.generate(2)
    add_field("Extend", 1; numbers=(("DistMax", 1.0), ("SizeMax", 1.0),
                                             ("Power", 1.0)),
              lists=(("CurvesList", [1]),))
    samples, sizes = curve_samples_and_sizes(1)
    field = ExtendField(samples, sizes; dist_max=1.0, size_max=1.0, power=1.0)
    node_tags, coordinates, _ = gmsh.model.mesh.getNodes(2, 1, false, false)
    isempty(node_tags) && error("Extend fixture has no surface-owned nodes")
    points = [(coordinates[3i - 2], coordinates[3i - 1], coordinates[3i])
              for i in eachindex(node_tags)]
    entities = fill((2, 1), length(points))
    check_direct("Extend sampled boundary", gmsh_probe_nodes(1, node_tags),
                 tessella_values(field, points; entities=entities);
                 atol=2e-12, rtol=2e-12, classification="direct sampled context")
end

function check_external_process()
    mktempdir() do directory
        helper = joinpath(directory, "external_field.jl")
        open(helper, "w") do io
            write(io, """
                while true
                    eof(stdin) && break
                    x = read(stdin, Float64)
                    y = read(stdin, Float64)
                    z = read(stdin, Float64)
                    any(isnan, (x, y, z)) && break
                    write(stdout, hypot(x, y, z) + 0.07)
                    flush(stdout)
                end
                """)
        end
        julia = Sys.which("julia")
        julia === nothing && error("Julia executable is not on PATH for ExternalProcess")
        command = Base.shell_escape(julia, "--startup-file=no", helper)
        setup_model("external_process")
        add_field("ExternalProcess", 1; strings=(("CommandLine", command),))
        points = [(0.0, 0.0, 0.0), (3.0, 4.0, 0.0), (-1.0, 2.0, 2.0)]
        gmsh_values = try
            gmsh_probe(1, points)
        finally
            gmsh.model.mesh.field.remove(1) # closes the process launched by this field
        end
        field = ExternalProcessField(command)
        tessella_result = try
            tessella_values(field, points)
        finally
            close(field)
        end
        check_direct("ExternalProcess protocol", gmsh_values, tessella_result;
                     classification="direct protocol")
    end
end

# Existing mesh-observed checks. Exact coordinates are intentionally not
# compared: the two curve discretizers solve different placement problems.
function common_gmsh_line(name)
    setup_model(name)
    gmsh.model.geo.addPoint(-2.0, 0.0, 0.0, 1.0, 1)
    gmsh.model.geo.addPoint(2.0, 0.0, 0.0, 1.0, 2)
    gmsh.model.geo.addLine(1, 2, 1)
end

function finish_gmsh_line()
    gmsh.model.geo.synchronize()
    for (option, value) in (("Mesh.MeshSizeExtendFromBoundary", 0.0),
                            ("Mesh.MeshSizeFromPoints", 0.0),
                            ("Mesh.MeshSizeFromCurvature", 0.0),
                            ("Mesh.MeshSizeMin", 0.01),
                            ("Mesh.MeshSizeMax", 1.0))
        gmsh.option.setNumber(option, value)
    end
    gmsh.model.mesh.generate(1)
    _, coordinates, _ = gmsh.model.mesh.getNodes(1, 1, true, false)
    nodes = sort!(unique(coordinates[1:3:end]))
    length(nodes) >= 3 || error("Gmsh line mesh returned only $(length(nodes)) nodes")
    return nodes
end

function gmsh_threshold_line()
    common_gmsh_line("threshold_line")
    gmsh.model.geo.addPoint(0.0, 0.0, 0.0, 1.0, 3)
    gmsh.model.geo.synchronize()
    add_field("Distance", 1; lists=(("PointsList", [3]),))
    add_field("Threshold", 2;
              numbers=(("InField", 1.0), ("DistMin", 0.0), ("DistMax", 1.0),
                       ("SizeMin", 0.1), ("SizeMax", 0.5)))
    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    return finish_gmsh_line()
end

function gmsh_box_line()
    common_gmsh_line("box_line")
    gmsh.model.geo.synchronize()
    add_field("Box", 1;
              numbers=(("VIn", 0.1), ("VOut", 0.5), ("XMin", -0.5),
                       ("XMax", 0.5), ("YMin", -1.0), ("YMax", 1.0),
                       ("ZMin", -1.0), ("ZMax", 1.0), ("Thickness", 0.5)))
    gmsh.model.mesh.field.setAsBackgroundMesh(1)
    return finish_gmsh_line()
end

function gmsh_ball_line()
    common_gmsh_line("ball_line")
    gmsh.model.geo.synchronize()
    add_field("Ball", 1;
              numbers=(("VIn", 0.1), ("VOut", 0.5), ("XCenter", 0.0),
                       ("YCenter", 0.0), ("ZCenter", 0.0), ("Radius", 0.5),
                       ("Thickness", 0.5)))
    gmsh.model.mesh.field.setAsBackgroundMesh(1)
    return finish_gmsh_line()
end

function gmsh_cylinder_line()
    common_gmsh_line("cylinder_line")
    gmsh.model.geo.synchronize()
    add_field("Cylinder", 1;
              numbers=(("VIn", 0.1), ("VOut", 0.5), ("XCenter", 0.0),
                       ("YCenter", 0.0), ("ZCenter", 0.0), ("XAxis", 0.5),
                       ("YAxis", 0.0), ("ZAxis", 0.0), ("Radius", 1.0)))
    gmsh.model.mesh.field.setAsBackgroundMesh(1)
    return finish_gmsh_line()
end

function gmsh_frustum_line()
    common_gmsh_line("frustum_line")
    gmsh.model.geo.synchronize()
    add_field("Frustum", 1;
              numbers=(("X1", 0.0), ("Y1", 0.0), ("Z1", -1.0),
                       ("X2", 0.0), ("Y2", 0.0), ("Z2", 1.0),
                       ("InnerR1", 0.0), ("OuterR1", 1.0),
                       ("InnerR2", 0.0), ("OuterR2", 1.0),
                       ("InnerV1", 0.1), ("OuterV1", 0.5),
                       ("InnerV2", 0.1), ("OuterV2", 0.5)))
    add_field("Box", 2; numbers=(("VIn", 0.5), ("VOut", 0.5)))
    add_field("Min", 3; lists=(("FieldsList", [1, 2]),))
    gmsh.model.mesh.field.setAsBackgroundMesh(3)
    return finish_gmsh_line()
end

tessella_line(field) = first.(mesh_segment((-2.0, 0.0, 0.0),
                                            (2.0, 0.0, 0.0), field)[1])

function spacing_metrics(nodes)
    gaps = diff(nodes)
    midpoints = (nodes[1:end-1] + nodes[2:end]) ./ 2
    near = gaps[abs.(midpoints) .< 0.4]
    far = gaps[abs.(midpoints) .> 1.2]
    isempty(near) && error("line mesh has no near-region gaps")
    isempty(far) && error("line mesh has no far-region gaps")
    return (count=length(nodes), minimum=minimum(gaps), maximum=maximum(gaps),
            near=median(near), far=median(far))
end

function mesh_result(name, gmsh_metrics, tessella_metrics)
    push!(MESH_RESULTS, (name=name, gmsh=gmsh_metrics, tessella=tessella_metrics))
    println("  MESH_APPROX ", rpad(name, 20), " Gmsh=", gmsh_metrics,
            " Tessella=", tessella_metrics)
end

function check_mesh_grading()
    gt = gmsh_threshold_line()
    threshold = ThresholdField(DistanceField(points=[(0.0, 0.0, 0.0)]);
        dist_min=0.0, dist_max=1.0, size_min=0.1, size_max=0.5)
    gm = spacing_metrics(gt); tm = spacing_metrics(tessella_line(threshold))
    abs(gm.count - tm.count) <= 1 || error("Threshold node-count mismatch: Gmsh=$gm Tessella=$tm")
    0.75 <= tm.minimum / gm.minimum <= 1.25 || error("Threshold fine-spacing mismatch: Gmsh=$gm Tessella=$tm")
    0.90 <= tm.maximum / gm.maximum <= 1.10 || error("Threshold coarse-spacing mismatch: Gmsh=$gm Tessella=$tm")
    mesh_result("Distance→Threshold", gm, tm)

    gb = spacing_metrics(gmsh_box_line())
    tb = spacing_metrics(tessella_line(BoxField(-0.5, 0.5, -1.0, 1.0, -1.0, 1.0;
                                               vin=0.1, vout=0.5, thickness=0.5)))
    abs(gb.count - tb.count) <= 4 || error("Box node-count mismatch: Gmsh=$gb Tessella=$tb")
    0.95 <= tb.near / gb.near <= 1.05 || error("Box interior-spacing mismatch: Gmsh=$gb Tessella=$tb")
    0.80 <= tb.far / gb.far <= 1.25 || error("Box exterior-spacing mismatch: Gmsh=$gb Tessella=$tb")
    tb.near < tb.far && gb.near < gb.far || error("Box did not refine its interior")
    mesh_result("Box", gb, tb)

    gba = spacing_metrics(gmsh_ball_line())
    tba = spacing_metrics(tessella_line(BallField((0.0, 0.0, 0.0), 0.5;
                                                   vin=0.1, vout=0.5, thickness=0.5)))
    abs(gba.count - tba.count) <= 4 || error("Ball node-count mismatch: Gmsh=$gba Tessella=$tba")
    0.90 <= tba.near / gba.near <= 1.10 || error("Ball interior-spacing mismatch: Gmsh=$gba Tessella=$tba")
    0.80 <= tba.far / gba.far <= 1.25 || error("Ball exterior-spacing mismatch: Gmsh=$gba Tessella=$tba")
    tba.near < tba.far && gba.near < gba.far || error("Ball did not refine its interior")
    mesh_result("Ball", gba, tba)

    gc = spacing_metrics(gmsh_cylinder_line())
    tc = spacing_metrics(tessella_line(CylinderField((0.0, 0.0, 0.0),
        (0.5, 0.0, 0.0), 1.0; vin=0.1, vout=0.5)))
    abs(gc.count - tc.count) <= 4 || error("Cylinder node-count mismatch: Gmsh=$gc Tessella=$tc")
    0.75 <= tc.near / gc.near <= 1.25 || error("Cylinder interior-spacing mismatch: Gmsh=$gc Tessella=$tc")
    0.75 <= tc.far / gc.far <= 1.25 || error("Cylinder exterior-spacing mismatch: Gmsh=$gc Tessella=$tc")
    tc.near < tc.far && gc.near < gc.far || error("Cylinder did not refine its interior")
    mesh_result("Cylinder", gc, tc)

    gf = spacing_metrics(gmsh_frustum_line())
    frustum = MinSize((FrustumField((0.0, 0.0, -1.0), (0.0, 0.0, 1.0);
        inner_r1=0.0, outer_r1=1.0, inner_r2=0.0, outer_r2=1.0,
        inner_v1=0.1, outer_v1=0.5, inner_v2=0.1, outer_v2=0.5),
        ConstantSize(0.5)))
    tf = spacing_metrics(tessella_line(frustum))
    abs(gf.count - tf.count) <= 4 || error("Frustum node-count mismatch: Gmsh=$gf Tessella=$tf")
    0.75 <= tf.near / gf.near <= 1.25 || error("Frustum inner-spacing mismatch: Gmsh=$gf Tessella=$tf")
    0.80 <= tf.far / gf.far <= 1.25 || error("Frustum outer-spacing mismatch: Gmsh=$gf Tessella=$tf")
    tf.near < tf.far && gf.near < gf.far || error("Frustum did not grade radially")
    mesh_result("Frustum", gf, tf)
end

gmsh.initialize([GMSH_EXECUTABLE, "-v", "0"])
try
    string(gmsh.GMSH_API_VERSION) == TARGET_GMSH_VERSION || error(
        "expected Gmsh API $TARGET_GMSH_VERSION, got $(gmsh.GMSH_API_VERSION)")
    println("GMSH_RUNTIME_OK cli=$GMSH_CLI_VERSION api=$(gmsh.GMSH_API_VERSION)")
    println("  executable=$GMSH_EXECUTABLE")
    println("  julia_api=$GMSH_API_FILE")
    println("  library=$(gmsh.lib)")
    println("POINTWISE_ORACLE MeshSizeFieldView")

    check_matheval()
    check_differential_operators()
    check_coordinate_maps()
    check_structured()
    check_entity_fields()
    check_postview()
    check_anisotropic_scalar_operators()
    check_boundary_layer_scalar()
    check_attractor_aniso_scalar()
    check_octree()
    check_extend()
    check_external_process()
    check_mesh_grading()

    println("CONTEXT_DEPENDENT_SKIPS count=$(length(CONTEXT_SKIPS))")
    for skip in CONTEXT_SKIPS
        println("  CONTEXT_SKIP $(skip.name): $(skip.reason)")
    end
    direct_samples = sum(result.samples for result in DIRECT_RESULTS; init=0)
    println("SIZE_FIELD_DIFFERENTIAL_OK gmsh=$(gmsh.GMSH_API_VERSION) " *
            "plugin_calls=$(PLUGIN_CALLS[]) direct_cases=$(length(DIRECT_RESULTS)) " *
            "direct_samples=$direct_samples mesh_cases=$(length(MESH_RESULTS)) " *
            "context_skips=$(length(CONTEXT_SKIPS))")
finally
    gmsh.finalize()
end
