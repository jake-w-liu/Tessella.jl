# Gmsh's own .geo files, Gmsh side: the files of gmsh_examples.jl through
# Gmsh's Python API on one thread.
#
#     python benchmark/gmsh_examples_gmsh.py <gmsh source> [results.toml]
#
# Every file is run twice and the second run is reported, like the Tessella
# driver. Requires the `gmsh` Python package (pip install gmsh).

import json
import os
import sys
import time

import gmsh

DIRECTORIES = ("tutorials", "examples/simple_geo")

# Files that need what neither driver provides: a post-processing view to
# merge (t8, t9), a surface mesh to remesh (t13) or a background mesh (t17).
SKIPPED = ("t8.geo", "t9.geo", "t13.geo", "t17.geo")


def count_cells(dim):
    return sum(len(tags) for tags in gmsh.model.mesh.getElements(dim)[1])


def measure(path):
    opened = None
    seconds, nodes, surface, volume = [], [], [], []
    failure = ""
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("General.NumThreads", 1)
        start = time.perf_counter()
        gmsh.open(path)
        opened = time.perf_counter() - start
        for dim in (1, 2, 3):
            start = time.perf_counter()
            gmsh.model.mesh.generate(dim)
            seconds.append(time.perf_counter() - start)
            nodes.append(len(gmsh.model.mesh.getNodes()[0]))
            surface.append(count_cells(2))
            volume.append(count_cells(3))
    except Exception as err:
        step = "open" if opened is None else "generate(%d)" % (len(seconds) + 1)
        failure = step + ": " + " ".join(str(err).split())[:120]
    finally:
        gmsh.finalize()
    return opened, seconds, nodes, surface, volume, failure


def main(arguments):
    if not arguments:
        sys.exit("usage: gmsh_examples_gmsh.py <gmsh source directory> [results.toml]")
    root = os.path.abspath(arguments[0])
    out = open(arguments[1], "w") if len(arguments) > 1 else sys.stdout
    start_directory = os.getcwd()
    try:
        gmsh.initialize()
        version = gmsh.option.getString("General.Version")
        gmsh.finalize()
        print("# Gmsh %s, Python API, 1 thread" % version, file=out)
        for directory in DIRECTORIES:
            for name in sorted(os.listdir(os.path.join(root, directory))):
                if not name.endswith(".geo") or name in SKIPPED:
                    continue
                path = os.path.join(root, directory, name)
                os.chdir(os.path.dirname(path))
                measure(path)
                opened, seconds, nodes, surface, volume, failure = measure(path)
                os.chdir(start_directory)
                print("\n[[file]]", file=out)
                print('name = "%s/%s"' % (directory, name), file=out)
                print("open_seconds = %.4f" % (-1.0 if opened is None else opened), file=out)
                print("generate_seconds = [%s]" % ", ".join("%.4f" % s for s in seconds), file=out)
                print("nodes = %s" % json.dumps(nodes), file=out)
                print("surface_cells = %s" % json.dumps(surface), file=out)
                print("volume_cells = %s" % json.dumps(volume), file=out)
                print("error = %s" % json.dumps(failure), file=out)
                out.flush()
    finally:
        os.chdir(start_directory)
        if out is not sys.stdout:
            out.close()


main(sys.argv[1:])
