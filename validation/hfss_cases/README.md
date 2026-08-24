# HFSS 22-case native meshing (STATUS #12, meshing half)

Every geometry in the HFSS v10 User Guide worked examples (ch. 5–10) meshed **from
scratch with Tessella's own primitives / raw triangulated surfaces — no gmsh, no
OpenCASCADE**. Each `build_case(id)` returns a valid, watertight, conforming tet mesh
of a representative geometry of that case's class (topology / shape; exact guide
dimensions are not required for the meshing demonstration).

Run the whole set:

```
julia --project=<Tessella.jl> validation/hfss_cases/hfss_case_meshes.jl
```

It prints, per case, `tets / valid / watertight / volume`. The set is also
regression-pinned in `test/integration/hfss_cases_test.jl` (all 22 assert valid + watertight;
box-assembly cases assert exact volume; multi-region cases assert conforming
interfaces + every region filled).

## Coverage (22/22 valid + watertight)

| case | geometry | Tessella route |
|---|---|---|
| 5.1 | UHF monopole probe over ground | `cylinder_surface` + `box_surface`, `tetrahedralize_conforming_exact` |
| 5.2 | conical horn | raw frustum surface → `tetrahedralize` |
| 5.3 | probe-fed patch (ground/substrate/patch/via) | `mesh_box_regions` |
| 5.4 | slot patch | `mesh_box_regions` |
| 5.5 | SAR spherical bowl | raw lat/long sphere surface → `tetrahedralize` |
| 5.6 | CPW bowtie | raw triangular-prism plate → `tetrahedralize` |
| 5.7 | endfire array unit cell (periodic) | `mesh_box` |
| 6.1 | magic tee | `mesh_box_regions` (fused waveguide arms) |
| 6.2 | coax bend | `mesh_cylinder` (per leg) |
| 6.3 | rat-race ring hybrid | raw annulus (washer) surface → `tetrahedralize` (genus-1) |
| 6.4 | coax stub | `mesh_cylinder` |
| 6.5 | microstrip wave port | `mesh_box_regions` |
| 6.6 | dielectric resonator | `cylinder_surface` puck + `box_surface`, `tetrahedralize_conforming_exact` |
| 7.1 | coupled-line bandpass filter | `mesh_box_regions` |
| 7.2 | stub bandstop filter | `mesh_box_regions` |
| 8.1 | LVDS differential pair | `mesh_box_regions` |
| 8.2 | segmented return plane | `mesh_box_regions` |
| 8.3 | non-ideal finite planes | `mesh_box_regions` |
| 8.4 | return-path discontinuity + via | `mesh_box_regions` |
| 9.1 | finned heat sink | `mesh_box_regions` (base + fins union) |
| 9.2 | enclosure coax feed-through | `cylinder_surface` pin + `box_surface`, `tetrahedralize_conforming_exact` |
| 10.1 | square-spiral inductor | `mesh_box_regions` (winding + substrate) |

## Scope note

This is the **meshing** half of the 22-case regression. The full-wave **solve** of
each case is the ASCENT project's external campaign. Case **9.2 (the enclosure) is
the only gmsh-impossible case** — the one that motivates Tessella's existence — and
it is meshed at its literal `.geo` dimensions **and solved in ASCENT**
(`validation/enclosure_literal/`, `CASE_9_2_OK`). The compact 9.2 stand-in here
(pin-in-cavity) is the same class; the literal fixture lives in `enclosure_literal/`.
The other 21 cases are not gmsh failures — standard tools mesh them — so meshing them
with Tessella demonstrates coverage rather than filling a capability gap.
