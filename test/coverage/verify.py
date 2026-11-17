"""Require executed function-body counters from the selected coverage process."""
from pathlib import Path
import sys

root = Path(__file__).resolve().parents[2]
mode = sys.argv[1] if len(sys.argv) == 2 else "kernel"
if mode not in ("kernel", "api") or len(sys.argv) > 2:
    raise SystemExit("Usage: verify.py [kernel|api]")
pid = int((root / "coverage" / ("api-runner-pid" if mode == "api" else "runner-pid")).read_text())
signatures = (
    ("src/geometry/Model.jl", "function add_point!(m::GeoModel, x, y, z; tag::Integer=0, mesh_size::Real=1.0,"),
    ("src/geometry/GeoExec.jl", "function _geo_emit_statement!(statements,locs,buf,buf_nl,buf_line0,cur_line,line_lfs)"),
    ("src/interfaces/APIViewProbe.jl", "function _view_probe_precise_coordinates!(location,geometry)"),
) if mode == "api" else (
    ("src/core/Predicates.jl", "function orient2(a, b, c)::Int"),
    ("src/core/MeshTypes.jl", "function triangle_area(a, b, c)"),
    ("src/geometry/NURBS.jl", "function nurbs_eval(c::NURBSCurve, u)"),
)
for source, signature in signatures:
    path = root / source
    coverage = path.with_name(f"{path.name}.{pid}.cov")
    lines = path.read_text().splitlines()
    start = lines.index(signature)
    counters = coverage.read_text().splitlines()
    # A top-level include/definition counter is insufficient: inspect body lines.
    executed = sum(
        counter.split(maxsplit=1)[0].isdigit()
        and int(counter.split(maxsplit=1)[0]) > 0
        for counter in counters[start + 1 : start + 12]
    )
    if executed == 0:
        raise RuntimeError(f"No executed function-body coverage for {signature}")
    print(f"Executed function-body counters: {source}: {executed}")
