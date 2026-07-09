#!/usr/bin/env python3
"""Drift summary from scp_sim diag.tsv (complex_phi columns)."""
import sys

def load(path):
    rows = []
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) < len(hdr):
                continue
            rows.append({k: float(v) for k, v in zip(hdr, p)})
    return rows

def pick(rows, t):
    return min(rows, key=lambda r: abs(r["t"] - t))

for path in sys.argv[1:]:
    r = load(path)
    a, b, c = pick(r, 0), pick(r, 10), pick(r, 60)
    name = path.split("/")[-1].replace("_diag.tsv", "")
    dE_full = 100 * (c["E_total"] - a["E_total"]) / a["E_total"]
    dE_set  = 100 * (c["E_total"] - b["E_total"]) / b["E_total"]
    dQ = 100 * (c["Q_total"] - a["Q_total"]) / a["Q_total"] if a.get("Q_total") else float("nan")
    print(f"{name:12s} E0={a['E_total']:9.3f}  drift(0->60)={dE_full:+7.3f}%  "
          f"drift(10->60)={dE_set:+7.3f}%  dQ={dQ:+7.4f}%  "
          f"s_max {a['s_max']:.4f}->{c['s_max']:.4f}  Qth0={a['Q_theta']:.2f}->{c['Q_theta']:.2f}")
