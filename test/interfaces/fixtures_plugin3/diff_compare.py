# Compare `tessella_out.txt` (Tessella) against `oracle_out.json`
# (Gmsh 4.15.2) for the batch-3 plugin cases. Cases are paired in order;
# records are canonicalized per node (sorted (coords, values)) so record and
# vertex ordering do not produce false diffs.
import json, os

FX = os.path.dirname(os.path.abspath(__file__))
WIDTH = {"P": 1, "L": 2, "T": 3, "Q": 4, "S": 4, "H": 8, "I": 6, "Y": 5}
FIXTURE_NAMES = {"scalar", "scalarB", "vector", "tensor", "pts", "pts2",
                 "lonlat"}


def canon_records(records):
    out = []
    for dtype, count, payload in records:
        arity = WIDTH[dtype[1]]
        w = len(payload) // count
        nval = w - 3 * arity
        for r in range(count):
            seg = payload[r * w:(r + 1) * w]
            nodes = []
            for n in range(arity):
                xyz = (round(seg[n], 6), round(seg[arity + n], 6),
                       round(seg[2 * arity + n], 6))
                vals = tuple(round(seg[3 * arity + n + k * arity], 6)
                             for k in range(nval // arity))
                nodes.append((xyz, vals))
            out.append((dtype, tuple(sorted(nodes))))
    return sorted(out)


def tessella_cases():
    cases = []
    for line in open(os.path.join(FX, "tessella_out.txt")):
        line = line.rstrip("\n")
        if line.startswith("CASE "):
            _, plugin, fixture = line.split(" ", 2)
            cases.append({"plugin": plugin, "fixture": fixture,
                          "views": [], "error": ""})
        elif line.startswith("ERROR "):
            cases[-1]["error"] = line[6:]
        elif line.startswith("VIEW "):
            head, name = line.split(" name=", 1)
            cases[-1]["views"].append({"name": name.strip('"'),
                                       "records": []})
        elif line.startswith("REC "):
            ty, cnt, nums = line[4:].split(" ", 2)
            payload = [float(x) for x in nums.split(",")]
            cases[-1]["views"][-1]["records"].append(
                [ty, int(cnt), payload])
    return cases


def main():
    oracle = json.load(open(os.path.join(FX, "oracle_out.json")))
    tess = tessella_cases()
    tidx = 0
    for case in oracle:
        plugin = case["plugin"]
        t = tess[tidx]; tidx += 1
        if plugin == "Summation2":
            # Oracle stores the new view records directly under "records".
            assert t["plugin"] == "Summation", (t["plugin"], plugin)
            orecs = canon_records(
                [(e[0], e[1], e[2]) for e in case["records"]])
            trecs = []
            for v in t["views"]:
                trecs.extend(canon_records(
                    [(e[0], e[1], e[2]) for e in v["records"]]))
            status = "MATCH" if orecs == trecs else "DIFF"
            print(f"{'Summation':20s} oracle={len(orecs)} "
                  f"tessella={len(trecs)} {status}")
            if status == "DIFF":
                for r in orecs[:4]:
                    print("   oracle:", r)
                for r in trecs[:4]:
                    print("   tess:  ", r)
            continue
        assert t["plugin"] == plugin, (t["plugin"], plugin)
        orecs, trecs = [], []
        onames, tnames = [], []
        views = list(case.get("new_views", {}).values())
        for v in views:
            nm = v.get("name", "")
            # The fixture view only appears alongside a fresh output view;
            # for in-place plugins the (renamed) fixture view IS the output.
            if nm in FIXTURE_NAMES and len(views) > 1:
                continue
            onames.append(nm)
            orecs.extend(canon_records(
                [(e[0], e[1], e[2]) for e in v["records"]]))
        for v in t["views"]:
            tnames.append(v["name"])
            trecs.extend(canon_records(
                [(e[0], e[1], e[2]) for e in v["records"]]))
        status = "MATCH" if orecs == trecs else "DIFF"
        print(f"{plugin:20s} names={onames}/{tnames} "
              f"oracle={len(orecs)} tessella={len(trecs)} {status}")
        if status == "DIFF":
            os_, ts_ = set(map(str, orecs)), set(map(str, trecs))
            for r in [r for r in orecs if str(r) not in ts_][:4]:
                print("   oracle only:", r)
            for r in [r for r in trecs if str(r) not in os_][:4]:
                print("   tess only:", r)
        if t["error"]:
            print("   tessella error:", t["error"])


main()
