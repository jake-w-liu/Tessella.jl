"""Require executed function-body counters from the selected coverage run.

    verify.py [kernel|api]           .cov files written by the coverage process
    verify.py api-lcov FILE...       LCOV reports of the API coverage parts

The API suite runs in parts on separate machines (TESSELLA_TEST_SHARD); a
function body counts as executed when any part executed it.
"""
from pathlib import Path
import sys

root = Path(__file__).resolve().parents[2]
modes = ("kernel", "api", "api-lcov")
mode = sys.argv[1] if len(sys.argv) >= 2 else "kernel"
if mode not in modes or (mode != "api-lcov" and len(sys.argv) > 2) or (mode == "api-lcov" and len(sys.argv) < 3):
    raise SystemExit("Usage: verify.py [kernel|api] | verify.py api-lcov FILE...")

signatures = (
    ("src/geometry/Model.jl", "function add_point!(m::GeoModel, x, y, z; tag::Integer=0, mesh_size::Real=1.0,"),
    ("src/geometry/GeoExec.jl", "function _geo_emit_statement!(statements,locs,buf,buf_nl,buf_line0,cur_line,line_lfs)"),
    ("src/interfaces/APIViewProbe.jl", "function _view_probe_precise_coordinates!(location,geometry)"),
) if mode != "kernel" else (
    ("src/core/Predicates.jl", "function orient2(a, b, c)::Int"),
    ("src/core/MeshTypes.jl", "function triangle_area(a, b, c)"),
    ("src/geometry/NURBS.jl", "function nurbs_eval(c::NURBSCurve, u)"),
)


def lcov_hits(paths):
    """Executed line numbers per source file (relative to the repository), over all reports."""
    hits = {}
    for report in paths:
        current = None
        for line in Path(report).read_text(encoding="utf-8").splitlines():
            if line.startswith("SF:"):
                name = Path(line[3:].strip())
                try:
                    name = name.resolve().relative_to(root)
                except ValueError:
                    pass
                current = hits.setdefault(name.as_posix().split("Tessella.jl/")[-1], set())
            elif line.startswith("DA:") and current is not None:
                number, count = line[3:].split(",")[:2]
                if int(count) > 0:
                    current.add(int(number))
    return hits


hits = lcov_hits(sys.argv[2:]) if mode == "api-lcov" else None
if mode != "api-lcov":
    pid = int((root / "coverage" / ("api-runner-pid" if mode == "api" else "runner-pid")).read_text(encoding="utf-8"))

for source, signature in signatures:
    path = root / source
    lines = path.read_text(encoding="utf-8").splitlines()
    start = lines.index(signature)
    # A top-level include/definition counter is insufficient: inspect body lines.
    if hits is not None:
        executed = sum(number in hits.get(source, set()) for number in range(start + 2, start + 13))
    else:
        counters = path.with_name(f"{path.name}.{pid}.cov").read_text(encoding="utf-8").splitlines()
        executed = sum(
            counter.split(maxsplit=1)[0].isdigit()
            and int(counter.split(maxsplit=1)[0]) > 0
            for counter in counters[start + 1 : start + 12]
        )
    if executed == 0:
        raise RuntimeError(f"No executed function-body coverage for {signature}")
    print(f"Executed function-body counters: {source}: {executed}")
