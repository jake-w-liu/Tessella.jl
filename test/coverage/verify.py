"""Require executed function-body counters from the compiled coverage process."""
from pathlib import Path

root = Path(__file__).resolve().parents[2]
pid = int((root / "coverage" / "runner-pid").read_text())
for source, signature in (
    ("src/core/Predicates.jl", "function orient2(a, b, c)::Int"),
    ("src/core/MeshTypes.jl", "function triangle_area(a, b, c)"),
    ("src/geometry/NURBS.jl", "function nurbs_eval(c::NURBSCurve, u)"),
):
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
