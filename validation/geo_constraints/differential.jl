#!/usr/bin/env julia
# Differential for Gmsh 4.15.2 `.geo` meshing-constraint and lifecycle
# statements against `Tessella.GeoExec.execute_geo`. Each case runs the same
# source through `gmsh.open` (+ `gmsh.model.mesh.generate` for mesh cases) and
# through `execute_geo`, then compares the observable state: entity inventory,
# physical groups and names, and — for the deterministic structured cases —
# the full node-coordinate multiset and per-entity element counts.
#
# Intentional nonclaims: recombine/algorithm/smoother effects inside the
# `.geo` mesh path (Tessella records the constraints; Gmsh auto-applies
# recombine at generate), compound-entity meshing (an API `mesh.generate`
# unit-grouping feature), and unstructured-mesh topology equivalence — only
# targeted probes (embedded-node presence, element counts) compare there.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella
using Tessella.Model: model_physical_groups, model_entities_for_physical_group,
                      model_physical_name

const TARGET_GMSH_VERSION = "4.15.2"

function find_gmsh_api()
    configured = get(ENV, "GMSH_JULIA_API", "")
    if !isempty(configured)
        isfile(configured) || error(
            "GMSH_JULIA_API does not name a file: $configured")
        return realpath(configured)
    end
    candidates = String[]
    executable = Sys.which("gmsh")
    if executable !== nothing
        prefix = dirname(dirname(realpath(executable)))
        append!(candidates, (joinpath(prefix, "lib", "gmsh.jl"),
                             joinpath(prefix, "lib64", "gmsh.jl")))
    end
    append!(candidates, ("/opt/homebrew/lib/gmsh.jl",
                         "/opt/homebrew/opt/gmsh/lib/gmsh.jl",
                         "/usr/local/opt/gmsh/lib/gmsh.jl"))
    for candidate in unique(candidates)
        isfile(candidate) && return realpath(candidate)
    end
    error("Gmsh 4.15.2 Julia API not found; set GMSH_JULIA_API")
end

include(find_gmsh_api())
gmsh.GMSH_API_VERSION == TARGET_GMSH_VERSION || error(
    "expected Gmsh API $TARGET_GMSH_VERSION, got $(gmsh.GMSH_API_VERSION)")

const SQUARE = """
Point(1) = {0,0,0,0.25};
Point(2) = {1,0,0,0.25};
Point(3) = {1,1,0,0.25};
Point(4) = {0,1,0,0.25};
Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
"""

const BOX = """
Point(1) = {0,0,0,0.5}; Point(2) = {1,0,0,0.5};
Point(3) = {1,1,0,0.5}; Point(4) = {0,1,0,0.5};
Point(5) = {0,0,1,0.5}; Point(6) = {1,0,1,0.5};
Point(7) = {1,1,1,0.5}; Point(8) = {0,1,1,0.5};
Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
Line(5) = {5,6}; Line(6) = {6,7}; Line(7) = {7,8}; Line(8) = {8,5};
Line(9) = {1,5}; Line(10) = {2,6}; Line(11) = {3,7}; Line(12) = {4,8};
Curve Loop(1) = {1,2,3,4}; Curve Loop(2) = {5,6,7,8};
Curve Loop(3) = {1,10,-5,-9}; Curve Loop(4) = {2,11,-6,-10};
Curve Loop(5) = {3,12,-7,-11}; Curve Loop(6) = {4,9,-8,-12};
Plane Surface(1) = {1}; Plane Surface(2) = {2};
Plane Surface(3) = {3}; Plane Surface(4) = {4};
Plane Surface(5) = {5}; Plane Surface(6) = {6};
Surface Loop(1) = {1,2,3,4,5,6};
Volume(1) = {1};
"""

# Two coplanar squares sharing edge Line(2) (surface 2 references it as -2).
const TWO_SQUARES = """
Point(1) = {0,0,0,0.25}; Point(2) = {1,0,0,0.25};
Point(3) = {1,1,0,0.25}; Point(4) = {0,1,0,0.25};
Point(5) = {2,0,0,0.25}; Point(6) = {2,1,0,0.25};
Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
Line(5) = {5,6}; Line(6) = {6,3}; Line(7) = {2,5};
Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
Curve Loop(2) = {7,5,6,-2};
Plane Surface(2) = {2};
"""

# A standalone dangling curve — negative delete tags hit the mirror/abs lookup.
const SQUARE_PLUS_CURVE = SQUARE * "Line(9) = {1,4};\n"

# mode :state  — compare entities + physicals only
# mode :mesh   — additionally compare mesh output (see compare_mesh)
# mode :error  — both sides must reject the source
const CASES = (
    (name=:transfinite_quad_surface, mode=:mesh, dim=2,
     source=SQUARE * """
     Transfinite Curve{1,2,3,4} = 5;
     Transfinite Surface{1} = {1,2,3,4} Right;
     """),
    (name=:transfinite_signed_curve, mode=:mesh, dim=2,
     # The negative tag flips the progression direction on curve 1; the
     # boundary node set (and therefore the whole structured grid) differs
     # observably from the positive form. Opposite edges share node counts.
     source=SQUARE * """
     Transfinite Curve{-1} = 6 Using Progression 2;
     Transfinite Curve{2,3,4} = 6 Using Progression 2;
     Transfinite Surface{1};
     """),
    (name=:transfinite_volume, mode=:mesh, dim=3,
     source=BOX * """
     Transfinite Curve{:} = 4;
     Transfinite Surface{:};
     Transfinite Volume{1} = {1,2,3,4,5,6,7,8};
     """),
    (name=:transfquadtri_blocker, mode=:error, dim=3,
     # Native kernel emits tetrahedra only — the QuadTri flag is an explicit
     # blocker (Gmsh succeeds with its HAVE_QUADTRI path).
     source=BOX * """
     Transfinite Curve{:} = 3;
     Transfinite Surface{:};
     Transfinite Volume{1};
     TransfQuadTri{1};
     """),
    (name=:delete_owned_refused, mode=:state, dim=0,
     source=SQUARE * "Delete{Curve{1};}\n"),
    (name=:delete_surface_then_curves, mode=:state, dim=0,
     source=SQUARE * "Delete{Surface{1};}\nDelete{Curve{1,2};}\n"),
    (name=:delete_signed_tags, mode=:state, dim=0,
     # Surface{-1} is a signed no-op; Point{-2}/Curve{-9} delete via the
     # absolute/mirrored lookup (curve 9 is dangling, point 2 is owned).
     source=SQUARE_PLUS_CURVE * "Delete{Surface{-1};}\nDelete{Curve{-9};}\nDelete{Point{-2};}\n"),
    (name=:recursive_delete_shared, mode=:state, dim=0,
     # Curve 2 and its endpoints are still owned by surviving surface 2.
     source=TWO_SQUARES * "Recursive Delete{Surface{1};}\n"),
    (name=:recursive_delete_owned, mode=:state, dim=0,
     # A box face is owned by the volume — Recursive Delete refuses it.
     source=BOX * "Recursive Delete{Surface{1};}\n"),
    (name=:delete_wildcard, mode=:state, dim=0,
     source=BOX * "Delete{Volume{1};}\nDelete{Surface{:};}\n"),
    (name=:delete_physical_selector, mode=:state, dim=0,
     source=SQUARE * """
     Physical Curve(5) = {1,2};
     Delete{Surface{1};}
     Delete{Physical Curve{5};}
     """),
    (name=:physical_stale_filter, mode=:state, dim=0,
     source=SQUARE * """
     Physical Curve("border",5) = {1,2};
     Delete{Surface{1};}
     Delete{Curve{1,2,3,4};}
     """),
    (name=:physical_resurrect, mode=:state, dim=0,
     source=SQUARE * """
     Physical Point(5) = {1};
     Delete{Surface{1};}
     Delete{Curve{1,2,3,4};}
     Delete{Point{1,2,3,4};}
     Point(1) = {0,0,0,0.1};
     """),
    (name=:delete_physicals, mode=:state, dim=0,
     source=SQUARE * """
     Physical Curve("border",5) = {1,2};
     Delete Physicals;
     Physical Point("later") = {1};
     """),
    (name=:delete_model, mode=:state, dim=0,
     source=SQUARE * """
     x = 42;
     Physical Point("kept",7) = {1};
     Delete Model;
     Point(3) = {0,0,0,1};
     Physical Point("after",8) = {3};
     """),
    (name=:delete_all, mode=:state, dim=0,
     source=SQUARE * """
     x = 42;
     Physical Point("gone",7) = {1};
     Delete All;
     Point(3) = {0,0,0,1};
     """),
    (name=:setmaxtag_assign, mode=:state, dim=0,
     source="""
     Point(10) = {0,0,0,1};
     SetMaxTag Point(3);
     Point(newp) = {1,0,0,1};
     Line(newl) = {4,10};
     """),
    (name=:compound_curve_state, mode=:state, dim=0,
     # Missing members are skipped; the later spec wins duplicate membership.
     source=SQUARE * """
     Compound Curve{1,2};
     Compound Curve{2,3};
     Compound Curve{4,99,-6};
     """),
    (name=:recombine_state, mode=:recombine_gap, dim=2,
     source=SQUARE * "Recombine Surface{1} = 30;\n"),
    (name=:smoother_algorithm_state, mode=:state, dim=0,
     source=SQUARE * """
     Smoother Surface{1} = 4;
     MeshAlgorithm Surface{1} = 6;
     MeshSizeFromBoundary Surface{1} = 1;
     Degenerated Curve{1};
     ReverseMesh Curve{1};
     RelocateMesh Point{1};
     """),
    (name=:embedded_delete, mode=:embed, dim=2,
     # Off-grid coordinates so the embed node cannot coincide with a natural
     # triangulation vertex.
     source=SQUARE * """
     Point(9) = {0.37,0.51,0,0.1};
     Point{9} In Surface{1};
     Delete Embedded{Surface{1};}
     """),
    (name=:embedded_kept, mode=:embed, dim=2,
     source=SQUARE * """
     Point(9) = {0.37,0.51,0,0.1};
     Point{9} In Surface{1};
     """),
    (name=:err_bad_surface_corners, mode=:error, dim=0,
     source=SQUARE * "Transfinite Surface{1} = {1,2};\n",
     expect="corner"),
    (name=:err_unknown_corner, mode=:error, dim=0,
     source=SQUARE * "Transfinite Surface{1} = {1,2,3,99};\n"),
    (name=:err_settag, mode=:error, dim=0,
     source=SQUARE * "SetTag Surface(1, 9);\n"),
    (name=:err_recursive_nonlist, mode=:error, dim=0,
     source=SQUARE * "Recursive Delete Model;\n"),
    (name=:err_delete_unknown_word, mode=:error, dim=0,
     source=SQUARE * "Delete Magic{Surface{1};}\n"),
    (name=:err_delete_unknown_name, mode=:error, dim=0,
     source=SQUARE * "Delete not_a_variable;\n"),
)

function tessella_run(source, dim)
    return mktemp() do path, io
        write(io, source)
        close(io)
        try
            return (execution=execute_geo(path; mesh_dim=dim), error=nothing)
        catch err
            err isa InterruptException && rethrow()
            return (execution=nothing, error=err)
        end
    end
end

function gmsh_run(path, dim)
    gmsh.clear()
    gmsh.parser.clear()
    # `gmsh.merge` rather than `gmsh.open`: OpenProject holds CTX::lock from
    # entry until normal exit, so a parse error thrown mid-merge leaves the
    # lock stuck and every later `gmsh.open` silently no-ops ("I'm busy!").
    # `gmsh.merge` calls MergeFile directly; `gmsh.clear`+`parser.clear`
    # reproduce the model/symbol reset OpenProject performs.
    try
        gmsh.merge(path)
    catch err
        return (error=sprint(showerror, err),)
    end
    if dim > 0
        try
            gmsh.model.mesh.generate(dim)
        catch err
            return (error=sprint(showerror, err),)
        end
    end
    return (error=nothing,)
end

# Node-coordinate sets compared by tolerance matching at 1e-7 — the same
# tolerance the transfinite-curve differential uses. Gmsh obtains node
# positions through recursive trapezoidal integration (≈1e-8 jitter that
# even reorders nodes within a grid column) while Tessella evaluates the
# distribution law analytically, so exact multisets and sorted-positional
# pairing are not meaningful parity criteria.
const COORD_TOLERANCE = 1e-7

function gmsh_coord_list(dim)
    tags, coords, _ = gmsh.model.mesh.getNodes(dim, -1, true, false)
    return [Tuple(coords[3i-2:3i]) for i in eachindex(tags)]
end

function tessella_coord_list(mesh)
    return [Tuple(col) for col in eachcol(mesh.coords)]
end

coord_close(l, r) = all(i -> abs(l[i] - r[i]) <= COORD_TOLERANCE,
                        eachindex(l))

function coords_equal(left, right)
    length(left) == length(right) || return false
    unmatched = collect(right)
    for l in left
        index = findfirst(r -> coord_close(l, r), unmatched)
        index === nothing && return false
        deleteat!(unmatched, index)
    end
    return true
end

function gmsh_elements(dim, tag)
    types, _, nodes = gmsh.model.mesh.getElements(dim, tag)
    counts = Dict{Int,Int}()
    for (type, list) in zip(types, nodes)
        _, _, _, num_nodes, _, _ = gmsh.model.mesh.getElementProperties(
            Int(type))
        counts[Int(type)] = get(counts, Int(type), 0) +
                            length(list) ÷ Int(num_nodes)
    end
    return counts
end

function compare_entities(tmodel, name)
    for (dim, dict) in ((0, tmodel.points), (1, tmodel.curves),
                        (2, tmodel.surfaces), (3, tmodel.volumes))
        tessella_tags = sort!(collect(keys(dict)))
        gmsh_tags = sort!(Int.(last.(gmsh.model.getEntities(dim))))
        extra = setdiff(gmsh_tags, tessella_tags)
        missing = setdiff(tessella_tags, gmsh_tags)
        isempty(missing) || error(
            "$name: Tessella dim-$dim entities $missing absent in Gmsh " *
            "(gmsh=$gmsh_tags tessella=$tessella_tags)")
        for tag in extra
            # Compounds materialize as discrete entities at sync — the only
            # Gmsh-side extras a `.geo` constraint script can add.
            kind = gmsh.model.getType(dim, tag)
            occursin("Discrete", kind) || error(
                "$name: Gmsh dim-$dim entity $tag ($kind) absent in Tessella " *
                "(gmsh=$gmsh_tags tessella=$tessella_tags)")
        end
    end
    return nothing
end

function compare_physicals(tmodel, name)
    gmsh_groups = sort!(Tuple{Int,Int}[
        (Int(d), Int(t)) for (d, t) in gmsh.model.getPhysicalGroups()])
    tessella_groups = sort!(collect(model_physical_groups(tmodel)))
    gmsh_groups == tessella_groups || error(
        "$name: physical groups differ (gmsh=$gmsh_groups " *
        "tessella=$tessella_groups)")
    for (dim, tag) in gmsh_groups
        gmsh_members = sort!(Int.(last.(
            gmsh.model.getEntitiesForPhysicalGroup(dim, tag))))
        tessella_members = sort!(
            model_entities_for_physical_group(tmodel, dim, tag))
        gmsh_members == tessella_members || error(
            "$name: physical ($dim,$tag) members differ " *
            "(gmsh=$gmsh_members tessella=$tessella_members)")
        gmsh_name = try
            gmsh.model.getPhysicalName(dim, tag)
        catch
            ""
        end
        tessella_name = model_physical_name(tmodel, dim, tag)
        gmsh_name == tessella_name || error(
            "$name: physical ($dim,$tag) name differs " *
            "(gmsh=$(repr(gmsh_name)) tessella=$(repr(tessella_name)))")
    end
    # Names outlive their groups (and `Delete Model` keeps them while
    # `Delete All` drops them) — sweep every recorded name, not just the
    # live groups.
    for ((dim, tag), tessella_name) in tmodel.physical_names
        gmsh_name = try
            gmsh.model.getPhysicalName(dim, tag)
        catch
            ""
        end
        gmsh_name == tessella_name || error(
            "$name: physical ($dim,$tag) name differs " *
            "(gmsh=$(repr(gmsh_name)) tessella=$(repr(tessella_name)))")
    end
    return nothing
end

function compare_mesh_counts(tmesh, dim, name)
    gmsh_tags, _, _ = gmsh.model.mesh.getNodes()
    total_gmsh = length(gmsh_tags)
    total_tessella = size(tmesh.coords, 2)
    total_gmsh == total_tessella || error(
        "$name: node counts differ (gmsh=$total_gmsh " *
        "tessella=$total_tessella)")
    coords_equal(gmsh_coord_list(dim), tessella_coord_list(tmesh)) || error(
        "$name: node coordinate sets differ")
    if dim == 2
        gmsh_tris = get(gmsh_elements(2, -1), 2, 0)
        gmsh_tris == size(tmesh.tris, 2) || error(
            "$name: triangle counts differ (gmsh=$gmsh_tris " *
            "tessella=$(size(tmesh.tris, 2)))")
    end
    # dim == 3: Gmsh emits structured hexahedra while the native kernel
    # decomposes the same node lattice into tetrahedra — the node multiset
    # above is the parity check, element types differ by construction.
    return nothing
end

function check_embed_point(tmesh, present, name)
    # An embedded point is observable as a vertex of the surface's elements —
    # the point entity itself keeps its own (dim-0) mesh vertex either way,
    # so plain node presence does not distinguish embedded from free.
    target = (0.37, 0.51, 0.0)
    _, _, gmsh_enodes = gmsh.model.mesh.getElements(2, -1)
    referenced = Set{Int}(Int(t) for list in gmsh_enodes for t in list)
    gmsh_tags, gmsh_coords, _ = gmsh.model.mesh.getNodes()
    gmsh_present = any(
        i -> gmsh_tags[i] in referenced &&
             coord_close(Tuple(gmsh_coords[3i-2:3i]), target),
        eachindex(gmsh_tags))
    gmsh_present == present || error(
        "$name: Gmsh embed surface-node presence $gmsh_present, " *
        "expected $present")
    tessella_present = any(
        ((i, t),) -> coord_close(Tuple(tmesh.coords[:, tmesh.tris[i, t]]),
                                 target),
        Iterators.product(axes(tmesh.tris, 1), axes(tmesh.tris, 2)))
    tessella_present == present || error(
        "$name: Tessella embed surface-node presence $tessella_present, " *
        "expected $present")
    return nothing
end

results = String[]
gaps = String[]

gmsh.initialize(String[], false, false)
try
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("General.NumThreads", 1)
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    gmsh.option.setNumber("Mesh.FlexibleTransfinite", 0)
    runtime_version = gmsh.option.getString("General.Version")
    (runtime_version == TARGET_GMSH_VERSION ||
     startswith(runtime_version, TARGET_GMSH_VERSION * "-")) || error(
        "expected Gmsh runtime $TARGET_GMSH_VERSION, got $runtime_version")

    for case in CASES
        t = tessella_run(case.source, case.dim)
        g = mktemp() do path, io
            write(io, case.source)
            close(io)
            gmsh_run(path, case.dim)
        end

        if case.mode == :error
            # TransfQuadTri is a Tessella-only blocker — Gmsh succeeds.
            if case.name == :transfquadtri_blocker
                t.error !== nothing || error(
                    "$(case.name): Tessella did not reject TransfQuadTri")
                g.error === nothing || error(
                    "$(case.name): Gmsh rejected TransfQuadTri: $(g.error)")
                push!(gaps, "transfquadtri: Tessella native-kernel blocker " *
                            "(Gmsh HAVE_QUADTRI path only)")
            else
                t.error === nothing && error(
                    "$(case.name): Tessella accepted an invalid source")
                g.error === nothing && error(
                    "$(case.name): Gmsh accepted an invalid source")
                expect = get(case, :expect, "")
                isempty(expect) ||
                    occursin(expect, sprint(showerror, t.error)) ||
                    error("$(case.name): Tessella error " *
                          "$(sprint(showerror, t.error)) lacks " *
                          repr(expect))
            end
            push!(results, string(case.name))
            continue
        end

        t.error === nothing || error(
            "$(case.name): Tessella threw $(sprint(showerror, t.error))")
        g.error === nothing || error(
            "$(case.name): Gmsh threw $(g.error)")
        tm = t.execution.model

        compare_entities(tm, case.name)
        compare_physicals(tm, case.name)

        if case.mode == :mesh
            compare_mesh_counts(t.execution.mesh, case.dim, case.name)
        elseif case.mode == :embed
            check_embed_point(t.execution.mesh,
                              case.name == :embedded_kept, case.name)
        elseif case.mode == :recombine_gap
            tm.meshing.recombine[(2, 1)] == 30.0 || error(
                "$(case.name): Tessella lost the Recombine 30 constraint")
            gmsh_quads = get(gmsh_elements(2, -1), 3, 0)
            gmsh_quads > 0 || error(
                "$(case.name): Gmsh produced no quads under Recombine")
            push!(gaps,
                  "recombine: constraint stored; `.geo` mesh path keeps " *
                  "simplices (recombine is the API post-pass) — gmsh_quads=" *
                  "$gmsh_quads tessella_tris=$(size(t.execution.mesh.tris, 2))")
        elseif case.name == :compound_curve_state
            # The last spec wins curve 2; missing/negative members are
            # skipped. Gmsh additionally materializes compound entities.
            specs = tm.meshing.compounds
            any(spec -> spec.second == [2, 3], specs) || error(
                "$(case.name): last-wins compound membership broken: $specs")
        elseif case.name == :smoother_algorithm_state
            tm.meshing.smoothing[(2, 1)] == 4 || error("smoothing record lost")
            tm.meshing.algorithm[(2, 1)] == 6 || error("algorithm record lost")
            tm.meshing.size_from_boundary[(2, 1)] ||
                error("size-from-boundary record lost")
            1 in tm.meshing.degenerated || error("degenerated flag lost")
            get(tm.meshing.reverse, (1, 1), false) ||
                error("reverse flag lost")
        elseif case.name == :delete_physicals
            tm.physical_names[(1, 5)] == "border" || error(
                "$(case.name): physical name did not survive Delete Physicals")
        elseif case.name == :delete_all
            isempty(tm.physical_names) || error(
                "$(case.name): Delete All kept physical names")
        end
        push!(results, string(case.name))
    end

    println("GEO_CONSTRAINTS_DIFFERENTIAL_OK gmsh=$runtime_version " *
            "cases=$(length(results)) " *
            "checks=entities+physicals+mesh+errors " *
            "documented_gaps=$(length(gaps))")
    for gap in gaps
        println("  gap: $gap")
    end
finally
    gmsh.isInitialized() != 0 && gmsh.finalize()
end
