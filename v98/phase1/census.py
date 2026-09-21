#!/usr/bin/env python3
"""Enumerate cycle-closed contact graphs from the NetworkX atlas (n ≤ 7).

Writes a TSV of isomorphism types and a short summary to stdout.
No kernel. Companion to CLOSURE.md.
"""
from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path

import networkx as nx
from networkx.generators.atlas import graph_atlas_g


OUT = Path(__file__).resolve().parent / "runs"


def invariants(g: nx.Graph) -> dict:
    n = g.number_of_nodes()
    e = g.number_of_edges()
    b1 = e - n + 1
    bridges = list(nx.bridges(g))
    cycle_closed = n >= 1 and nx.is_connected(g) and len(bridges) == 0 and b1 >= 1
    # connectivity: isolated / small graphs
    if n <= 1:
        kappa = 0
    elif not nx.is_connected(g):
        kappa = 0
    else:
        kappa = nx.node_connectivity(g)
    try:
        girth = nx.girth(g)
        if girth is None or girth > n + e:
            girth = 0
    except (nx.NetworkXError, ValueError):
        girth = 0
    odd = False
    if nx.is_connected(g) and n >= 3:
        try:
            odd = not nx.is_bipartite(g)
        except nx.NetworkXError:
            odd = True
    # Laplacian spectrum, rounded for grouping (exact polys live in FriCAS)
    spec = []
    if n >= 1:
        ev = nx.laplacian_spectrum(g)
        spec = [round(float(x), 6) for x in ev]
    return {
        "n": n,
        "e": e,
        "b1": b1,
        "bridges": len(bridges),
        "cycle_closed": int(cycle_closed),
        "kappa": kappa,
        "f_cut": int(kappa >= 2),
        "f_3": int(kappa >= 3),
        "girth": girth if girth != 0 else "",
        "f_odd": int(odd),
        "bipartite": int(nx.is_bipartite(g) if n else 0),
        "spec": spec,
        "name": g.name or "",
    }


def atlas_name(g: nx.Graph, idx: int) -> str:
    # Prefer small well-known names when the atlas graph matches.
    named = {
        (3, 3): "C3=K3",
        (4, 4): "C4",
        (4, 5): "K4-e",
        (4, 6): "K4",
        (5, 5): "C5",
        (6, 6): "C6",
        (5, 10): "K5",
        (6, 15): "K6",
        (7, 7): "C7",
        (7, 21): "K7",
    }
    key = (g.number_of_nodes(), g.number_of_edges())
    # Only unique (n,e) pairs among cycle-closed get the short name;
    # collisions stay as atlas index.
    return named.get(key, f"G{idx}")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    closed = []
    for idx, g in enumerate(graph_atlas_g()):
        if g.number_of_nodes() < 1 or g.number_of_nodes() > 7:
            continue
        inv = invariants(g)
        inv["atlas"] = idx
        inv["label"] = atlas_name(g, idx)
        rows.append(inv)
        if inv["cycle_closed"]:
            closed.append(inv)

    tsv = OUT / "census_atlas.tsv"
    fields = [
        "atlas",
        "label",
        "n",
        "e",
        "b1",
        "bridges",
        "cycle_closed",
        "kappa",
        "f_cut",
        "f_3",
        "girth",
        "f_odd",
        "bipartite",
        "spec",
    ]
    with tsv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in closed:
            rec = dict(r)
            rec["spec"] = ",".join(str(x) for x in r["spec"])
            w.writerow(rec)

    def count(pred) -> Counter:
        c: Counter = Counter()
        for r in closed:
            if pred(r):
                c[r["n"]] += 1
        return c

    all_n = count(lambda r: True)
    cut_n = count(lambda r: r["f_cut"])
    k3_n = count(lambda r: r["f_3"])
    even_n = count(lambda r: not r["f_odd"])
    even_cut = count(lambda r: not r["f_odd"] and r["f_cut"])
    even_k3 = count(lambda r: not r["f_odd"] and r["f_3"])

    def smallest(pred):
        cand = [r for r in closed if pred(r)]
        if not cand:
            return None
        cand.sort(key=lambda r: (r["n"], r["e"], r["atlas"]))
        first_n = cand[0]["n"]
        at_n = [r for r in cand if r["n"] == first_n]
        return cand[0], at_n

    summary = []
    summary.append(f"atlas graphs scanned: {len(rows)}")
    summary.append(f"cycle-closed (C1-C3): {len(closed)}")
    summary.append(f"by n (raw): {dict(sorted(all_n.items()))}")
    summary.append(f"by n (F-cut): {dict(sorted(cut_n.items()))}")
    summary.append(f"by n (F-3): {dict(sorted(k3_n.items()))}")
    summary.append(f"by n (F-odd, i.e. bipartite): {dict(sorted(even_n.items()))}")
    summary.append(f"by n (F-odd + F-cut): {dict(sorted(even_cut.items()))}")
    summary.append(f"by n (F-odd + F-3): {dict(sorted(even_k3.items()))}")

    classes = [
        ("C1-C3", lambda r: True),
        ("F-cut", lambda r: r["f_cut"]),
        ("F-3", lambda r: r["f_3"]),
        ("F-odd", lambda r: not r["f_odd"]),
        ("F-odd+F-cut", lambda r: not r["f_odd"] and r["f_cut"]),
        ("F-odd+F-3", lambda r: not r["f_odd"] and r["f_3"]),
    ]
    for name, pred in classes:
        got = smallest(pred)
        if not got:
            summary.append(f"smallest {name}: NONE")
            continue
        first, at_n = got
        labels = [r["label"] for r in at_n]
        summary.append(
            f"smallest {name}: n={first['n']} e={first['e']} "
            f"count_at_n={len(at_n)} labels={labels} "
            f"kappa={first['kappa']} girth={first['girth']} b1={first['b1']}"
        )

    # zoo vs line: how many nonisomorphic closed graphs per n
    zoo = defaultdict(list)
    for r in closed:
        zoo[r["n"]].append((r["e"], r["label"], r["kappa"], r["f_odd"]))
    summary.append("zoo (n: list of (e,label,kappa,odd)):")
    for n in sorted(zoo):
        summary.append(f"  n={n} k={len(zoo[n])}")

    text = "\n".join(summary) + "\n"
    (OUT / "census_summary.txt").write_text(text)
    print(text, end="")
    print(f"wrote {tsv}")


if __name__ == "__main__":
    main()
