#!/usr/bin/env python3
"""Fit a_rel from mf_pair_track / sfa_qball_track TSVs; score F/R/O gates."""
from __future__ import annotations
import csv, math, sys
from pathlib import Path

def load_mf_track(path: Path):
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    t, D = [], []
    for i, r in enumerate(rows):
        # mf_pair_track may not have physical t; use frame index * snap if needed
        tt = float(r.get("t", i))
        dd = float(r["D"])
        t.append(tt); D.append(dd)
    return t, D

def load_qball_track(path: Path):
    """sfa_qball_track --tsv: expect columns that allow two-cluster D."""
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    # group by frame
    from collections import defaultdict
    byf = defaultdict(list)
    for r in rows:
        fr = int(float(r.get("frame", r.get("f", 0))))
        byf[fr].append(r)
    t, D = [], []
    for fr in sorted(byf):
        cl = byf[fr]
        if len(cl) < 2:
            continue
        # pick two heaviest
        cl = sorted(cl, key=lambda r: -float(r.get("mass", r.get("M", 0))))
        def cen(r):
            return (float(r.get("cx", r.get("x", 0))),
                    float(r.get("cy", r.get("y", 0))),
                    float(r.get("cz", r.get("z", 0))))
        a, b = cen(cl[0]), cen(cl[1])
        d = math.sqrt(sum((a[i]-b[i])**2 for i in range(3)))
        t.append(float(cl[0].get("t", fr)))
        D.append(d)
    return t, D

def fit_arel(t, D, tmin_frac=0.15, tmax_frac=0.85):
    """Quadratic D = D0 + v0 t + 0.5 a t^2 on middle window; return a, D0, n."""
    if len(t) < 5:
        return float("nan"), float("nan"), 0
    t0, t1 = t[0], t[-1]
    lo = t0 + tmin_frac * (t1 - t0)
    hi = t0 + tmax_frac * (t1 - t0)
    ts, Ds = [], []
    for ti, di in zip(t, D):
        if lo <= ti <= hi and di > 0.5:
            ts.append(ti - t0)  # relative
            Ds.append(di)
    n = len(ts)
    if n < 5:
        return float("nan"), float("nan"), n
    # least squares for D = a0 + a1 s + a2 s^2  => a_rel = 2 a2
    # normal equations
    S0 = n
    S1 = sum(ts); S2 = sum(s*s for s in ts); S3 = sum(s**3 for s in ts); S4 = sum(s**4 for s in ts)
    T0 = sum(Ds); T1 = sum(s*d for s,d in zip(ts,Ds)); T2 = sum(s*s*d for s,d in zip(ts,Ds))
    # solve 3x3
    import numpy as np
    A = np.array([[S0,S1,S2],[S1,S2,S3],[S2,S3,S4]], dtype=float)
    b = np.array([T0,T1,T2], dtype=float)
    try:
        coef = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        return float("nan"), float("nan"), n
    a0, a1, a2 = coef
    a_rel = 2.0 * a2
    return float(a_rel), float(a0), n

def main(resdir: Path):
    resdir = Path(resdir)
    lines = ["# Track analysis", ""]
    force_rows = []
    for D in [8, 12, 16, 20, 24]:
        p = resdir / f"mf_force_D{D}_track.tsv"
        if not p.exists():
            lines.append(f"- mf_force_D{D}: **MISSING track**")
            continue
        t, Dd = load_mf_track(p)
        # physical time: if t is 0..n-1 * dt_wrong, use frame index * snap_dt=2.5
        if len(t) > 1 and (t[-1] - t[0]) < 1e-6:
            t = [i * 2.5 for i in range(len(Dd))]
        elif len(t) > 1 and t[-1] < 5 and len(t) > 10:
            t = [i * 2.5 for i in range(len(Dd))]
        a, D0, n = fit_arel(t, Dd)
        force_rows.append((D, a, D0, n, Dd[0] if Dd else float("nan"), Dd[-1] if Dd else float("nan")))
        sign = "attract" if a < 0 else ("repel" if a > 0 else "flat")
        # convention: D decreasing => attract => a_rel for separation should be negative if D'' < 0
        lines.append(f"- D={D}: a_rel={a:.6e} ({sign}) D0_fit={D0:.3f} n={n} D_start={Dd[0]:.3f} D_end={Dd[-1]:.3f}")

    lines.append("")
    lines.append("## Same-sign (elite)")
    same_rows = []
    for D in [12, 16, 20]:
        p = resdir / f"elite_same_D{D}_track.tsv"
        if not p.exists():
            lines.append(f"- elite_same_D{D}: **MISSING track**")
            continue
        t, Dd = load_qball_track(p)
        if len(t) > 1 and t[-1] < 5:
            t = [i * 2.5 for i in range(len(Dd))]
        a, D0, n = fit_arel(t, Dd)
        same_rows.append((D, a, n, Dd[0] if Dd else float("nan"), Dd[-1] if Dd else float("nan")))
        sign = "attract" if a < 0 else ("repel" if a > 0 else "flat")
        lines.append(f"- D={D}: a_rel={a:.6e} ({sign}) n={n} D_start={Dd[0]:.3f} D_end={Dd[-1]:.3f}")

    lines.append("")
    lines.append("## Orbit")
    for tag, vt in [("0p03", 0.03), ("0p05", 0.05)]:
        p = resdir / f"mf_orbit_R16_vt{tag}_track.tsv"
        if not p.exists():
            lines.append(f"- vt={vt}: **MISSING track**")
            continue
        t, Dd = load_mf_track(p)
        if len(Dd) < 3:
            lines.append(f"- vt={vt}: **INCOMPLETE track** (frames={len(Dd)}; re-run after orbit SFA ready)")
            continue
        if len(t) > 1 and t[-1] < 5:
            t = [i * 5.0 for i in range(len(Dd))]
        # crude: count local minima of D as half-revs if oscillatory
        dmin, dmax = min(Dd), max(Dd)
        # angle from in-plane if we only have D: use d'' oscillation count
        n_cross = 0
        mean = sum(Dd)/len(Dd)
        for i in range(1, len(Dd)):
            if (Dd[i-1]-mean)*(Dd[i]-mean) < 0:
                n_cross += 1
        lines.append(f"- vt={vt}: frames={len(Dd)} D_min={dmin:.3f} D_max={dmax:.3f} mean_crossings={n_cross} D0={Dd[0]:.3f} Dend={Dd[-1]:.3f}")

    # Gates
    lines.append("")
    lines.append("## Gates (AND scorecard)")
    # F: opposite a_rel < 0 for most D (D decreasing / attract)
    if force_rows:
        n_att = sum(1 for D,a,_,_,ds,de in force_rows if a < 0 or de < ds - 0.05)
        n = len(force_rows)
        if n_att >= max(2, n-1) and all(not math.isnan(a) for _,a,_,_,_,_ in force_rows):
            f_status = "PASS" if n_att == n else "SOFT"
        elif n_att == 0:
            f_status = "FAIL"
        else:
            f_status = "SOFT"
        lines.append(f"- **F Force:** {f_status} ({n_att}/{n} separations show attraction signature)")
    else:
        lines.append("- **F Force:** FAIL (no tracks)")
        f_status = "FAIL"

    if same_rows:
        n_rep = sum(1 for D,a,n,ds,de in same_rows if a > 0 or de > ds + 0.05)
        if n_rep >= 1:
            r_status = "PASS" if n_rep == len(same_rows) else "SOFT"
        else:
            r_status = "FAIL" if all(a < 0 for _,a,_,_,_ in same_rows if not math.isnan(a)) else "SOFT"
        lines.append(f"- **R Same-sign:** {r_status}")
    else:
        lines.append("- **R Same-sign:** SOFT (missing)")
        r_status = "SOFT"

    lines.append("- **L1b Instrumentation:** PASS if this file lists tracks for planned jobs")
    lines.append("- **O Orbit:** see orbit D(t) notes above — need ≥2 revs for PASS")

    out = resdir / "TRACK_ANALYSIS.md"
    out.write_text("\n".join(lines) + "\n")
    print(out.read_text())
    return 0

if __name__ == "__main__":
    import math
    d = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
    sys.exit(main(d))
