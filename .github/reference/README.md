# Feature reference

- [Model views as mesh-size fields](model_view_fields.md)

These notes document specific implemented domains and their certificates.
Some describe an earlier bounded increment; consult the
[current checklist](../development/CHECKLIST.md) and current source before
interpreting a historical restriction as the complete present API surface.
Examples use native Julia production.

| Area | References |
| --- | --- |
| Command line | [Controls](cli_controls.md), [pyramid option scope](geo_pyramid_option_scope.md) |
| Model files | [Import](model_import.md), [merge](model_merge.md), [export](model_export.md), [CAD imports](cad_mesh_imports.md) |
| Tensor sizing | [Explicit Surface API](api_metric_surfaces.md), [ordinary Surfaces](metric_model_surfaces.md), [constant-metric boxes](constant_metric_boxes.md) |
| Hybrid constraints | [Points](embedded_hybrid_points.md), [Curves](hybrid_embedded_curves.md), [Surfaces](hybrid_embedded_surfaces.md), [collections](hybrid_constraint_collections.md), [polygon sheets](hybrid_polygon_sheets.md) |
| Hybrid generation | [Prismatic fill](prismatic_fill.md), [Region batches](hybrid_batch_pyramids.md), [normal apices](hybrid_normal_apices.md), [relocation](hybrid_relocate_modes.md), [Winslow](hybrid_winslow.md) |
| Periodicity | [Hybrid boundaries](periodic_hybrid_boundaries.md), [shared hybrids](shared_periodic_hybrids.md) |
| Partitions | [Tetrahedra](partitioned_tets.md), [BRep](partitioned_brep.md), [declared topology](partitioned_declared_brep.md), [open seams](partitioned_open_seams.md) |
| STEP/IGES | [Mapped occurrences](step_mapped_occurrences.md), [assemblies](step_relationship_assemblies.md), [units](step_relationship_si_units.md), [Cartesian targets](step_cartesian_rigid_targets.md), [plane charts](iges_plane_charts.md) |
| Views | [POS import](session_pos_import.md), [time](view_time.md), [representations](view_representations.md), [list export](list_view_exports.md), [model export](model_view_exports.md) |
