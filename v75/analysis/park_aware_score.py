#!/usr/bin/env python3
"""Park-aware atom packaging score from mf_pair_track + diag.

c_Q_seed  = (Qc0 - Qc_end)/Qc0          # legacy (penalizes nuclear park)
c_Q_park  = (Qc_mid - Qc_end)/Qc_mid    # post-park drift only
c_L       = (mL0 - mL_end)/mL0
c_E       = |E_end - E0| / E0 / 0.20

PASS_park: c_Q_park<=0.15, c_L<=0.15, c_G=0, massL_end>0.5*mL0
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path

def parse_diag(path: Path):
    rows = []
    with open(path) as f:
        for line in f:
            p = line.split()
            if not p or p[0] == "t" or not p[0][0].isdigit() and p[0][0] != "-":
                try:
                    float(p[0])
                except Exception:
                    continue
            try:
                vals = [float(x) for x in p]
            except ValueError:
                continue
            if len(vals) < 26:
                continue
            rows.append({"t": vals[0], "E": vals[9], "Q": vals[14], "s_max": vals[17], "gauss": vals[25]})
    return rows

def parse_track(path: Path):
    rows = []
    with open(path) as f:
        for line in f:
            p = line.split()
            if not p or p[0] == "frame":
                continue
            try:
                vals = [float(x) for x in p]
            except ValueError:
                continue
            if len(vals) < 7:
                continue
            rows.append({
                "t": vals[1], "massC": vals[2], "massL": vals[3],
                "Qc": vals[4], "Ql": vals[5],
                "D": vals[12] if len(vals) > 12 else float("nan"),
            })
    return rows

def score(diag_path: Path, track_path: Path) -> dict:
    d = parse_diag(diag_path)
    tr = parse_track(track_path)
    if not d or not tr:
        return {"ok": False, "error": "missing diag/track"}
    mid = d[len(d)//2]
    end = d[-1]
    t0, te = tr[0], tr[-1]
    Qc0, Qc_end = t0["Qc"], te["Qc"]
    Qc_mid = mid["Q"]  # diag Q_phi ≈ nuclear
    mL0, mL = t0["massL"], te["massL"]
    E0, E = d[0]["E"], end["E"]
    gmax = max(abs(r["gauss"]) for r in d)

    def clip01(x):
        return max(0.0, min(1.0, x))

    c_Q_seed = clip01((Qc0 - Qc_end) / max(abs(Qc0), 1e-12))
    c_Q_park = clip01((abs(Qc_mid) - abs(Qc_end)) / max(abs(Qc_mid), 1e-12))
    c_L = clip01((mL0 - mL) / max(mL0, 1e-12))
    c_E = clip01(abs(E - E0) / max(abs(E0), 1e-12) / 0.20)
    c_G = 1.0 if gmax > 1e-10 else 0.0
    cost_seed = 0.35*c_Q_seed + 0.35*c_L + 0.20*c_E + 0.10*c_G
    cost_park = 0.35*c_Q_park + 0.35*c_L + 0.20*c_E + 0.10*c_G
    pass_park = (
        c_Q_park <= 0.15 and c_L <= 0.15 and c_G == 0
        and mL > 0.5 * mL0 and abs(Qc_end) > 0.5 * abs(Qc_mid)
    )
    return {
        "ok": True,
        "Qc0": Qc0, "Qc_mid": Qc_mid, "Qc_end": Qc_end,
        "massL0": mL0, "massL": mL, "E0": E0, "E": E, "gauss_max": gmax,
        "c_Q_seed": c_Q_seed, "c_Q_park": c_Q_park, "c_L": c_L, "c_E": c_E,
        "cost_seed": cost_seed, "cost_park": cost_park,
        "PASS_seed": cost_seed <= 0.15 and c_Q_seed <= 0.15 and c_L <= 0.15 and c_G == 0,
        "PASS_park": pass_park,
        "D0": t0["D"], "Dend": te["D"],
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("diag")
    ap.add_argument("track")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()
    s = score(Path(args.diag), Path(args.track))
    if args.json:
        print(json.dumps(s, indent=2))
    else:
        print(f"PASS_park={s.get('PASS_park')} PASS_seed={s.get('PASS_seed')}")
        print(f"  cost_park={s.get('cost_park'):.4f} cost_seed={s.get('cost_seed'):.4f}")
        print(f"  c_Q_park={s.get('c_Q_park'):.4f} c_Q_seed={s.get('c_Q_seed'):.4f} c_L={s.get('c_L'):.4f}")
        print(f"  Qc {s.get('Qc0'):.1f} → mid {s.get('Qc_mid'):.1f} → end {s.get('Qc_end'):.1f}")
        print(f"  massL {s.get('massL0'):.1f} → {s.get('massL'):.1f}")
    return 0 if s.get("ok") else 1

if __name__ == "__main__":
    raise SystemExit(main())
