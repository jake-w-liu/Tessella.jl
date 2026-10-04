using Test
using Tessella

@testset "Tessella" begin
    @testset "repository Julia-file layout" begin
        repository_root = normpath(joinpath(@__DIR__, ".."))
        function julia_files(relative_root)
            root = joinpath(repository_root, relative_root)
            paths = String[
                relpath(joinpath(directory, filename), root)
                for (directory, _, filenames) in walkdir(root)
                for filename in filenames if endswith(filename, ".jl")
            ]
            return sort!(paths)
        end

        source_files = julia_files("src")
        @test filter(path -> length(splitpath(path)) == 1, source_files) ==
              ["Tessella.jl"]
        source_domains = Set(("core", "fields", "geometry", "interfaces",
                              "meshing", "structured"))
        source_members = filter(!=("Tessella.jl"), source_files)
        @test all(length(splitpath(path)) == 2 for path in source_members)
        @test Set(first(splitpath(path)) for path in source_members) == source_domains

        test_files = julia_files("test")
        @test filter(path -> length(splitpath(path)) == 1, test_files) ==
              ["runtests.jl"]
        test_domains = Set(("core", "fields", "geometry", "integration",
                            "interfaces", "meshing", "structured"))
        test_members = filter(!=("runtests.jl"), test_files)
        @test all(length(splitpath(path)) == 2 for path in test_members)
        discovered_test_domains = Set(
            first(splitpath(path)) for path in test_members if
            first(splitpath(path)) != "tmp")
        @test discovered_test_domains == test_domains

        validation_files = julia_files("validation")
        @test filter(path -> length(splitpath(path)) == 1, validation_files) ==
              ["run_all.jl"]
        validation_members = filter(!=("run_all.jl"), validation_files)
        @test all(length(splitpath(path)) == 2 for path in validation_members)
        @test "support" in Set(first(splitpath(path)) for path in validation_members)
    end

    # Stage 0 — Foundations (CRC-gated, DEVELOPMENT.md discipline).
    include("core/predicates_test.jl")   # exact predicates vs exact-rational oracle
    include("core/allocation_audit_test.jl") # closure-boxing scan + hot-kernel allocation bounds
    include("core/meshtypes_test.jl")    # mesh container, topology, quality, checksum
    include("core/mesh_entity_topology_test.jl") # automatic/manual global edge/face ids
    include("core/mesh_point_location_test.jl") # robust simplex inversion and AABB lookup
    include("core/mesh_reference_geometry_test.jl") # robust simplex forward maps and Jacobians
    include("core/mesh_quadrature_test.jl") # bounded fixed-family rules and analytic moments
    include("core/mesh_function_spaces_test.jl") # arbitrary-order nodal plus H1/simplex-Whitney bases
    include("core/mesh_element_quality_test.jl") # Gmsh-shaped simplex quality measures
    include("core/transform_test.jl")    # validated affine transforms + orientation preservation
    include("interfaces/io_test.jl")     # .msh v2/v4 round-trip, STL, .geo scan
    include("interfaces/stream_cleanup_test.jl") # reader failure closes OS handles
    include("core/elements_test.jl")     # fixed/special Gmsh 4.15.2 records + mixed entity I/O
    include("meshing/recombine_test.jl") # deterministic triangle-to-quad recombination
    include("meshing/refine_test.jl")    # deterministic one-level uniform simplex refinement
    include("structured/transfinite_test.jl") # validated four-sided planar structured patches
    include("structured/transfinite_curve_test.jl") # straight Progression/Bump/Beta/HWall laws
    include("structured/transfinite_triangle_test.jl") # specific three-sided triangle/quad patches
    include("geometry/transfinite_triangle_orientation_test.jl") # public CAD-relative triangle winding
    include("geometry/transfinite_quad_orientation_test.jl") # public four-sided frames and CAD winding
    include("structured/transfinite_quad_test.jl") # recombined four-sided quadrangle patches
    include("structured/transfinite_volume_test.jl") # affine six-face structured volumes
    include("structured/structured_quadtri_test.jl") # warped/folded transition certificates
    include("structured/transfinite_prism_test.jl") # affine five-face transfinite prisms
    include("structured/transfinite_hex_test.jl") # affine six-face recombined hexahedra
    include("geometry/model_test.jl")    # entity kernel + .geo execution
    include("geometry/model_discrete_test.jl") # record atomicity and classification
    include("geometry/model_boolean_snapshot_test.jl") # owned Boolean operands and Delete cleanup
    include("geometry/model_boolean_multi_test.jl") # N-way OCC-cell Boolean operands
    include("geometry/model_topology_query_test.jl") # explicit topology query API
    include("geometry/model_entity_identity_test.jl") # entity names and atomic retagging
    include("geometry/model_entity_removal_test.jl") # dependency-safe recursive removal
    include("geometry/model_spatial_query_test.jl") # analytical bounding boxes
    include("geometry/model_entity_metadata_test.jl") # native type/property/partition metadata
    include("geometry/model_entity_evaluation_test.jl") # native geometry evaluation
    include("geometry/model_ruled_parametrization_test.jl") # ruled inverse and actual P2 node parameters
    include("geometry/model_entity_state_test.jl") # visibility/color/coordinates/attributes
    include("geometry/geo_geometry_expression_test.jl") # expression-backed geometry execution
    include("geometry/geo_control_flow_test.jl") # If/For/While/Function control flow + conditional operators
    include("geometry/geo_struct_test.jl")     # Struct/NameSpace/NameStruct namespaces + EOF rules
    include("geometry/geo_transform_test.jl") # Translate/Dilate/Rotate/Symmetry + Duplicata + coherence
    include("geometry/model_coherence_precision_test.jl") # local tolerance and rounded spatial bins
    include("geometry/geo_extrude_test.jl")  # translational Extrude + params + lateral merge
    include("geometry/geo_curved_test.jl")   # circle/ellipse arcs, ruled surfaces, .geo curved statements
    include("geometry/geo_spline_test.jl")   # Spline/BSpline/Bezier/Nurbs records + .geo statements
    include("geometry/occ_primitives_test.jl") # materialized Cylinder/Sphere/Cone/Torus OCC boundaries
    include("geometry/geo_recoverable_add_test.jl") # lenient built-in/OCC entity adds
    include("geometry/geo_list_variable_test.jl") # bounded numeric list variables
    include("geometry/geo_dynamic_tag_test.jl") # geometry/Physical allocation and lifecycle
    include("geometry/geo_set_max_tag_test.jl") # factory-aware max-tag counters
    include("geometry/geo_constraints_test.jl") # .geo Transfinite/Recombine/Delete meshing constraints
    include("geometry/geo_transfquadtri_test.jl") # boundary-diagonal mixed-volume transitions
    include("geometry/geo_quadtri_extrude_test.jl") # QuadTriAddVerts structured extrusion
    include("geometry/quadtri_nonew_templates_test.jl") # exact full-hex face relation
    include("geometry/quadtri_nonew_prism_templates_test.jl") # true triangular-prism relation
    include("geometry/quadtri_nonew_jacobian_test.jl") # complete P1 cell map certificates
    include("geometry/quadtri_nonew_nonhex_jacobian_test.jl") # exact pyramid/prism domain checks
    include("geometry/quadtri_nonew_prism_global_test.jl") # bounded six-corner hull separation
    include("geometry/geo_quadtri_nonew_test.jl") # isolated NoNew source-quad sweeps
    include("geometry/geo_quadtri_nonew_triangle_test.jl") # isolated triangular-source sweeps
    include("geometry/geo_quadtri_nonew_two_tri_test.jl") # jointly certified two-triangle source grids
    include("geometry/geo_quadtri_nonew_quad_patch_test.jl") # existing-pivot 2-by-2 quadrangle grids
    include("geometry/geo_mesh_identity_test.jl") # coincident orphan point identities
    include("geometry/model_mesh_identity_helpers_test.jl") # discrete and closed curve mesh-node ownership
    include("geometry/geo_mesh_size_test.jl") # Point sizing and topology-derived Physical groups
    include("geometry/geo_periodic_test.jl") # expression-backed periodic curves
    include("geometry/model_periodic_io_test.jl") # classified periodic/embedded MSH projection
    include("geometry/model_volume_io_test.jl") # classified volume-embedding MSH projection
    include("geometry/nurbs_test.jl")    # De Boor vs Bernstein/circle oracles
    include("meshing/boundarylayer_test.jl") # prismatic layer extrusion
    include("meshing/boundary_layer_fan_test.jl") # multi-region fan layers + core fill
    include("meshing/periodic_test.jl")  # periodic identification
    include("interfaces/api_test.jl") # synchronized session, detached mesh cache
    include("interfaces/api_mixed_cache_test.jl") # native mixed cache ownership and mutations
    include("interfaces/api_p2_certification_test.jl") # full P2 maps and atomic publication
    include("interfaces/api_p2_boundary_ownership_test.jl") # actual primary support and boundary closure
    include("interfaces/api_nonew_two_tri_boundary_test.jl") # two-triangle grid P2 carriers and lifecycle
    include("interfaces/api_nonew_quad_patch_boundary_test.jl") # quad-patch P2 carriers and lifecycle
    include("interfaces/api_generate01_test.jl") # native 0D/1D generation and sparse identities
    include("interfaces/api_mixed_queries_test.jl") # native mixed reference and function-space queries
    include("interfaces/api_mixed_advanced_test.jl") # native duplicate removal and partitioning
    include("interfaces/api_mixed_refine_test.jl") # native family-preserving mixed refinement
    include("interfaces/api_mixed_refine_support_test.jl") # inherited actual edge and face carriers
    include("interfaces/api_mesh_lifecycle_test.jl") # cached refinement and clearing
    include("interfaces/api_record_mutation_test.jl") # atomic sparse/dense tag mutation
    include("interfaces/api_refinement_classification_test.jl") # refinement identity and ownership
    include("interfaces/api_mesh_transform_test.jl") # atomic whole-cache affine transforms
    include("interfaces/api_mesh_data_test.jl") # detached node/element block queries
    include("interfaces/api_mesh_entity_topology_test.jl") # cached edge/face topology lifecycle
    include("interfaces/api_mesh_point_location_test.jl") # cached simplex point location
    include("interfaces/api_mesh_jacobian_test.jl") # cached simplex Jacobian maps
    include("interfaces/api_legacy_p2_jacobian_test.jl") # actual retained quadratic maps
    include("interfaces/api_legacy_p2_query_test.jl") # published quadratic types and node data
    include("interfaces/api_legacy_p2_quality_test.jl") # actual quadratic quality definitions
    include("interfaces/api_legacy_p2_quality_bounds_test.jl") # interior extrema and isotropy bounds
    include("interfaces/api_mesh_quadrature_test.jl") # session-independent reference rules
    include("interfaces/api_mesh_function_spaces_test.jl") # nodal/hierarchical basis-key API
    include("interfaces/api_mesh_element_quality_test.jl") # cached simplex quality queries
    include("interfaces/api_element_catalog_test.jl") # fixed-element type and property queries
    include("interfaces/api_boolean_lifecycle_test.jl") # synchronized Boolean ownership
    include("interfaces/api_topology_query_test.jl") # explicit topology queries
    include("interfaces/api_entity_identity_test.jl") # synchronized entity identity
    include("interfaces/api_entity_removal_test.jl") # synchronized entity removal
    include("interfaces/api_spatial_query_test.jl") # synchronized bounding boxes
    include("interfaces/api_entity_metadata_test.jl") # synchronized entity metadata
    include("interfaces/api_entity_evaluation_test.jl") # synchronized geometry evaluation
    include("interfaces/api_entity_state_test.jl") # synchronized entity presentation state
    include("interfaces/cli_test.jl") # bounded parser + non-destructive output
    include("interfaces/gui_test.jl") # validated headless command/state machine
    include("interfaces/post_test.jl") # owned scalar views + synchronized plugins
    include("interfaces/post_view_io_test.jl") # Gmsh .pos list-format view I/O

    # Stage 1 — 2-D meshing (CRC-gated).
    include("meshing/mesh2d_test.jl") # Delaunay: exact empty-circumcircle oracle

    # Stage 2 — Gmsh-compatible size fields + 1-D/surface meshing (CRC-gated).
    include("fields/sizefield_test.jl") # Distance/Threshold/Box/Min fields + local 3-D sizing
    include("fields/postview_field_test.jl") # PostView scalar/vector/tensor size fields
    include("fields/postview_highorder_test.jl") # high-order .pos records + PostViewField schemes
    include("meshing/mesh1d_test.jl") # size fields + graded edge meshing vs arc length
    include("meshing/meshsurface_test.jl") # planar/cylinder/parametric surface meshing

    # Stage 3 — 3-D meshing and certified boundary recovery (CRC-gated).
    include("meshing/mesh3d_test.jl") # 3-D Delaunay + volume filling + coax junction

    # Stage 4 — optimization, quality, and sliver removal (CRC-gated).
    include("meshing/optimize_test.jl") # tet quality report + Laplacian & ODT smoothing

    # Stage 5 — heal + native primitives + top-level pipeline.
    include("geometry/heal_test.jl") # surface diagnostics (open/non-manifold/degenerate…)
    include("geometry/geometry_test.jl") # box / cylinder / box-tunnel primitives
    include("geometry/cad_test.jl") # native analytical geometry (surfaces + exact imprints), no OCC
    include("geometry/brep_test.jl") # ISO-10303-21 STEP / IGES classified-solid import
    include("integration/pipeline_test.jl") # mesh_volume: validated-or-explicit-blocker

    # Stage 6 — high-order elements.
    include("meshing/highorder_test.jl") # quadratic (P2) tet generation + type-11 I/O
    include("meshing/highorder_jacobian_test.jl") # exact P2 Jacobian bounds for quad/hex/prism
    include("meshing/high_order_pyramid_test.jl") # full rational Pyr14 map bounds

    # Application: native meshes for all 22 HFSS User Guide case geometries (no gmsh/OCC).
    include("integration/hfss_cases_test.jl") # STATUS #12 meshing half — valid+watertight+conforming

    @testset "stage banner" begin
        @test Tessella.stage() isa Int
        @test Tessella.stage() == 6     # all package stages through P2/I/O are green
    end
end
