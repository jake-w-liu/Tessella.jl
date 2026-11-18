# Solid meshing benchmark, Gmsh side: the cases of solid_meshing.jl through
# Gmsh's Python API on one thread.
#
#     python benchmark/solid_meshing_gmsh.py [results.toml]
#     python benchmark/solid_meshing_gmsh.py --faces [results.toml]
#
# Every case is meshed twice and the second run is reported, like the
# Tessella driver. Requires the `gmsh` Python package (pip install gmsh).

import sys
import time

import gmsh
import numpy as np

occ = gmsh.model.occ


def drilled_block():
    block = occ.addBox(0, 0, 0, 2, 2, 2)
    drill = occ.addCylinder(1, 1, -1, 0, 0, 4, 0.5)
    occ.cut([(3, block)], [(3, drill)])


CASES = [
    ("box", lambda: occ.addBox(0, 0, 0, 1, 1, 1), (0.1, 0.05, 0.03, 0.02)),
    ("sphere", lambda: occ.addSphere(0, 0, 0, 1), (0.2, 0.1, 0.06, 0.04)),
    ("cylinder", lambda: occ.addCylinder(0, 0, 0, 0, 0, 2, 0.5), (0.1, 0.05, 0.03)),
    ("torus", lambda: occ.addTorus(0, 0, 0, 2, 0.5), (0.15, 0.1, 0.07)),
    ("drilled_block", drilled_block, (0.15, 0.08, 0.05)),
]


def measure(build, size):
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("General.NumThreads", 1)
        build()
        occ.synchronize()
        gmsh.option.setNumber("Mesh.MeshSizeMin", size)
        gmsh.option.setNumber("Mesh.MeshSizeMax", size)
        start = time.perf_counter()
        gmsh.model.mesh.generate(2)
        surface = time.perf_counter() - start
        start = time.perf_counter()
        gmsh.model.mesh.generate(3)
        volume = time.perf_counter() - start
        node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
        points = coordinates.reshape(-1, 3)
        position = np.zeros(int(node_tags.max()) + 1, dtype=np.int64)
        position[node_tags.astype(np.int64)] = np.arange(len(node_tags))
        triangles = 0
        types, tags, _ = gmsh.model.mesh.getElements(2)
        for kind, block in zip(types, tags):
            if kind == 2:
                triangles += len(block)
        types, tags, nodes = gmsh.model.mesh.getElements(3)
        tet_tags = np.array([], dtype=np.uint64)
        total = 0.0
        for kind, block, connectivity in zip(types, tags, nodes):
            if kind == 4:
                tet_tags = block
                p = points[position[connectivity.astype(np.int64)]].reshape(-1, 4, 3)
                total = float(np.einsum("ij,ij->i", p[:, 1] - p[:, 0],
                                        np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0])).sum() / 6)
        quality = gmsh.model.mesh.getElementQualities(tet_tags, "gamma") if len(tet_tags) else np.array([np.nan])
        return dict(nodes=len(node_tags), triangles=triangles, tetrahedra=len(tet_tags), volume=total,
                    gamma_min=float(np.min(quality)), gamma_mean=float(np.mean(quality)),
                    surface_seconds=surface, volume_seconds=volume)
    finally:
        gmsh.finalize()


def plane_grid(n):
    geo = gmsh.model.geo
    point = lambda i, j: i * (n + 1) + j + 1
    for i in range(n + 1):
        for j in range(n + 1):
            geo.addPoint(i, j, 0, tag=point(i, j))
    along = lambda i, j: i * n + j + 1
    across = lambda i, j: (n + 1) * n + i * (n + 1) + j + 1
    for i in range(n + 1):
        for j in range(n):
            geo.addLine(point(i, j), point(i, j + 1), tag=along(i, j))
    for i in range(n):
        for j in range(n + 1):
            geo.addLine(point(i, j), point(i + 1, j), tag=across(i, j))
    face = 0
    for i in range(n):
        for j in range(n):
            face += 1
            geo.addCurveLoop([along(i, j), across(i, j + 1), -along(i + 1, j), -across(i, j)], tag=face)
            geo.addPlaneSurface([face], tag=face)
    geo.synchronize()


def measure_faces(n):
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("General.NumThreads", 1)
        plane_grid(n)
        gmsh.option.setNumber("Mesh.MeshSizeMin", 0.2)
        gmsh.option.setNumber("Mesh.MeshSizeMax", 0.2)
        start = time.perf_counter()
        gmsh.model.mesh.generate(1)
        curves = time.perf_counter() - start
        start = time.perf_counter()
        gmsh.model.mesh.generate(2)
        surfaces = time.perf_counter() - start
        nodes = len(gmsh.model.mesh.getNodes()[0])
        triangles = sum(len(block) for block in gmsh.model.mesh.getElements(2)[1])
        return dict(nodes=nodes, triangles=triangles, curve_seconds=curves, surface_seconds=surfaces)
    finally:
        gmsh.finalize()


def main(arguments):
    faces = bool(arguments) and arguments[0] == "--faces"
    if faces:
        arguments = arguments[1:]
    out = open(arguments[0], "w") if arguments else sys.stdout
    try:
        gmsh.initialize()
        version = gmsh.option.getString("General.Version")
        gmsh.finalize()
        print(f'implementation = "Gmsh {version}"', file=out)
        print("threads = 1", file=out)
        if faces:
            for n in (5, 10, 20):
                measure_faces(n)
                result = measure_faces(n)
                print("\n[[grid]]", file=out)
                print(f"faces = {n * n}", file=out)
                print(f"nodes = {result['nodes']}", file=out)
                print(f"triangles = {result['triangles']}", file=out)
                print(f"curve_seconds = {result['curve_seconds']:.4f}", file=out)
                print(f"surface_seconds = {result['surface_seconds']:.4f}", file=out)
                out.flush()
            return
        for name, build, sizes in CASES:
            for size in sizes:
                measure(build, size)
                result = measure(build, size)
                print("\n[[case]]", file=out)
                print(f'solid = "{name}"', file=out)
                print(f"size = {size}", file=out)
                print(f"nodes = {result['nodes']}", file=out)
                print(f"triangles = {result['triangles']}", file=out)
                print(f"tetrahedra = {result['tetrahedra']}", file=out)
                print(f"volume = {result['volume']:.6f}", file=out)
                print(f"gamma_min = {result['gamma_min']:.4f}", file=out)
                print(f"gamma_mean = {result['gamma_mean']:.4f}", file=out)
                print(f"surface_seconds = {result['surface_seconds']:.4f}", file=out)
                print(f"volume_seconds = {result['volume_seconds']:.4f}", file=out)
                out.flush()
    finally:
        if out is not sys.stdout:
            out.close()


if __name__ == "__main__":
    main(sys.argv[1:])
