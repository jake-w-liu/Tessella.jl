"""Gmsh 4.15.2 oracle for plugin batch-3 differential testing.

For each (fixture, plugin, options) case: load the fixture .pos, set the
options, run the plugin, and dump every resulting view's list data as JSON.
Run with the Python that has the gmsh 4.15.2 SDK installed.
"""
import gmsh, json, sys, os

FX = os.path.dirname(os.path.abspath(__file__))


def list_data(tag):
    """Return [(recordType, count, flatPayload), ...] for a list view."""
    out = []
    try:
        dtypes, numeles, datas = gmsh.view.getListData(tag)
    except Exception:
        return out
    for dtype, numele, payload in zip(dtypes, numeles, datas):
        out.append([dtype, int(numele), [float(v) for v in payload]])
    return out


def run_case(fixture, plugin, numbers, strings=None, pre=None):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    views_before = set(gmsh.view.getTags())
    if fixture:
        gmsh.merge(os.path.join(FX, fixture))
    if pre:
        for fn, args in pre:
            getattr(gmsh.plugin, fn)(*args)
    for k, v in numbers.items():
        gmsh.plugin.setNumber(plugin, k, v)
    for k, v in (strings or {}).items():
        gmsh.plugin.setString(plugin, k, v)
    try:
        gmsh.plugin.run(plugin)
        err = ""
    except Exception as e:
        err = str(e)
    views_after = gmsh.view.getTags()
    new_views = [t for t in views_after if t not in views_before]
    result = {"error": err, "new_views": {}}
    for t in new_views:
        name = ""
        try:
            name = gmsh.view.option.getString(t, "Name") or ""
        except Exception:
            pass
        result["new_views"][str(t)] = {
            "name": name, "records": list_data(t)}
    # In-place plugins mutate the input view — also dump view 0.
    if not new_views and fixture:
        allv = gmsh.view.getTags()
        if allv:
            result["new_views"][str(allv[0])] = {
                "name": "(inplace)", "records": list_data(allv[0])}
    gmsh.finalize()
    return result


CASES = [
    ("scalar_tri.pos", "CutGrid",
     {"X0": 0, "Y0": 0, "Z0": 0, "X1": 1, "Y1": 0, "Z1": 0,
      "X2": 0, "Y2": 1, "Z2": 0, "NumPointsU": 3, "NumPointsV": 2,
      "ConnectPoints": 1, "View": 0}, None),
    ("scalar_tri.pos", "CutGrid",
     {"NumPointsU": 3, "NumPointsV": 2, "ConnectPoints": 0, "View": 0}, None),
    ("scalar_tri.pos", "CutBox",
     {"X0": 0, "Y0": 0, "Z0": 0, "X1": 1, "Y1": 0, "Z1": 0,
      "X2": 0, "Y2": 1, "Z2": 0, "X3": 0, "Y3": 0, "Z3": 1,
      "NumPointsU": 3, "NumPointsV": 2, "NumPointsW": 2,
      "ConnectPoints": 1, "Boundary": 0, "View": 0}, None),
    ("scalar_tri.pos", "CutBox",
     {"NumPointsU": 3, "NumPointsV": 2, "NumPointsW": 2,
      "ConnectPoints": 0, "Boundary": 1, "View": 0}, None),
    ("scalar_tri.pos", "CutBox",
     {"X0": 0, "Y0": 0, "Z0": 0, "X1": 1, "Y1": 0, "Z1": 0,
      "X2": 0, "Y2": 1, "Z2": 0, "X3": 0, "Y3": 0, "Z3": 1,
      "NumPointsU": 3, "NumPointsV": 2, "NumPointsW": 2,
      "ConnectPoints": 1, "Boundary": 1, "View": 0}, None),
    ("scalar_tri.pos", "CutGrid",
     {"X0": 0, "Y0": 0, "Z0": 0, "X1": 1, "Y1": 0, "Z1": 0,
      "X2": 0, "Y2": 1, "Z2": 0, "NumPointsU": 3, "NumPointsV": 1,
      "ConnectPoints": 1, "View": 0}, None),
    ("scalar_tri.pos", "CutParametric",
     {"MinU": 0, "MaxU": 1, "NumPointsU": 4, "MinV": 0, "MaxV": 1,
      "NumPointsV": 3, "ConnectPoints": 1, "View": 0},
     {"X": "u", "Y": "v", "Z": "0"}),
    ("vector_tet.pos", "StreamLines",
     {"X0": 0.25, "Y0": 0.25, "Z0": 0.1, "X1": 0.5, "Y1": 0.25, "Z1": 0.1,
      "X2": 0.25, "Y2": 0.5, "Z2": 0.1, "NumPointsU": 2, "NumPointsV": 1,
      "DT": 0.1, "MaxIter": 5, "TimeStep": 0, "View": 0}, None),
    # StreamLines with OtherView: SL records valued by the scalar source
    # (view index 1 = the scalar tet view in the combined fixture).
    ("vector_scalar.pos", "StreamLines",
     {"X0": 0.25, "Y0": 0.25, "Z0": 0.1, "X1": 0.5, "Y1": 0.25, "Z1": 0.1,
      "X2": 0.25, "Y2": 0.5, "Z2": 0.1, "NumPointsU": 2, "NumPointsV": 1,
      "DT": 0.1, "MaxIter": 5, "TimeStep": 0, "View": 0,
      "OtherView": 1}, None),
    ("vector_tet.pos", "Particles",
     {"X0": 0.25, "Y0": 0.25, "Z0": 0.1, "NumPointsU": 2, "NumPointsV": 1,
      "A2": 1, "A1": 0, "A0": 0, "DT": 0.1, "MaxIter": 5, "TimeStep": 0,
      "View": 0}, None),
    # points5 is non-cospherical: the Delaunay tessellation is unique, so the
    # Tessella differential is meaningful.
    ("points5.pos", "Triangulate", {"View": 0}, None),
    ("points5.pos", "Tetrahedralize", {"View": 0}, None),
    ("points4.pos", "NearestNeighbor", {"View": 0}, None),
    ("scalar_tri.pos", "Scal2Tens",
     {"NumberOfComponents": 3, "View0": 0, "View1": 0, "View2": 0}, None),
    ("tensor_tet.pos", "Eigenvectors",
     {"ScaleByEigenvalues": 1, "View": 0}, None),
    ("lonlat.pos", "LongitudeLatitude", {"View": 0}, None),
]


def main():
    results = []
    for fixture, plugin, numbers, strings in CASES:
        r = run_case(fixture, plugin, numbers, strings)
        r["plugin"] = plugin
        r["fixture"] = fixture
        r["numbers"] = numbers
        r["strings"] = strings or {}
        results.append(r)
    # Summation needs a second fixture merged — run it separately.
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.merge(os.path.join(FX, "scalar_tri.pos"))
    gmsh.merge(os.path.join(FX, "scalar_b.pos"))
    gmsh.plugin.setNumber("Summation", "View 0", 0)
    gmsh.plugin.setNumber("Summation", "View 1", 1)
    gmsh.plugin.setString("Summation", "Resuling View Name", "sum")
    gmsh.plugin.run("Summation")
    tags = gmsh.view.getTags()
    for t in tags:
        if t > 1:
            results.append({
                "plugin": "Summation2", "fixture": "scalar_tri+scalar_b",
                "records": list_data(t)})
    gmsh.finalize()
    with open(os.path.join(FX, "oracle_out.json"), "w") as f:
        json.dump(results, f, indent=1)
    print("oracle_out.json written:", len(results), "cases")


if __name__ == "__main__":
    main()
