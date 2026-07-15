#!/usr/bin/env python3
"""Park-aware + Z/N score from multi-fabric diag.tsv (GPU or CPU).

Nuclear-only (no L): uses Q_phi (bag) and Q_flux (EM proxy).
Atom (with L): also uses mf_pair_track massL if track provided.

PASS_nuc:
  c_Q_park  = |Q_mid - Q_end| / |Q_mid|  <= 0.15
  c_Qem_park = |Qf_mid - Qf_end| / max(|Qf_mid|,1) <= 0.20
  gauss floor, Q_end > 0.5 |Q_mid|

PASS_atom (needs track): PASS_nuc + c_L<=0.15 + massL_end>0.5*mL0
"""
from __future__ import annotations
import argparse, json, csv
from pathlib import Path


def load_diag_tsv(path: Path):
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        rows = []
        for row in r:
            try:
                rows.append({k: float(v) for k, v in row.items() if v not in (None, "")})
            except ValueError:
                continue
    return rows


def load_track(path: Path):
    rows = []
    with open(path) as f:
        for line in f:
            p = line.split()
            if not p or p[0] in ("frame", "t"):
                continue
            try:
                vals = [float(x) for x in p]
            except ValueError:
                continue
            if len(vals) < 6:
                continue
            # frame t nC nL massC massL Qc Ql ...
            rows.append({
                "t": vals[1] if len(vals) > 1 else vals[0],
                "massC": vals[4] if len(vals) > 4 else float("nan"),
                "massL": vals[5] if len(vals) > 5 else float("nan"),
                "Qc": vals[6] if len(vals) > 6 else float("nan"),
                "Ql": vals[7] if len(vals) > 7 else float("nan"),
            })
    return rows


def clip01(x: float) -> float:
    return max(0.0, min(1.0, x))


def score(diag_path: Path, track_path: Path | None = None,
          Z: int | None = None, N: int | None = None) -> dict:
    d = load_diag_tsv(diag_path)
    if not d:
        return {"ok": False, "error": f"empty diag {diag_path}"}
    mid = d[len(d) // 2]
    end = d[-1]
    t0 = d[0]

    def g(row, *keys, default=float("nan")):
        for k in keys:
            if k in row:
                return row[k]
        return default

    Q0 = g(t0, "Q_phi", "Q")
    Qm = g(mid, "Q_phi", "Q")
    Qe = g(end, "Q_phi", "Q")
    Qf0 = g(t0, "Q_flux")
    Qfm = g(mid, "Q_flux")
    Qfe = g(end, "Q_flux")
    E0 = g(t0, "E_total", "E")
    Ee = g(end, "E_total", "E")
    gmax = max(abs(g(r, "gauss_max", "gauss", default=0.0)) for r in d)
    s0 = g(t0, "s_max")
    se = g(end, "s_max")

    c_Q_seed = clip01((abs(Q0) - abs(Qe)) / max(abs(Q0), 1e-12))
    c_Q_park = clip01((abs(Qm) - abs(Qe)) / max(abs(Qm), 1e-12))
    c_Qem_park = clip01(abs(abs(Qfm) - abs(Qfe)) / max(abs(Qfm), 1.0))
    c_E = clip01(abs(Ee - E0) / max(abs(E0), 1e-12) / 0.20)
    c_G = 1.0 if gmax > 1e-10 else 0.0

    pass_nuc = (
        c_Q_park <= 0.15
        and c_Qem_park <= 0.20
        and c_G == 0
        and abs(Qe) > 0.5 * abs(Qm)
    )

    out = {
        "ok": True,
        "Z": Z, "N": N,
        "Q0": Q0, "Q_mid": Qm, "Q_end": Qe,
        "Qf0": Qf0, "Qf_mid": Qfm, "Qf_end": Qfe,
        "E0": E0, "E_end": Ee,
        "s0": s0, "s_end": se,
        "gauss_max": gmax,
        "c_Q_seed": c_Q_seed, "c_Q_park": c_Q_park,
        "c_Qem_park": c_Qem_park, "c_E": c_E, "c_G": c_G,
        "PASS_nuc": pass_nuc,
        "Qem_per_Z": (abs(Qfe) / Z) if Z and Z > 0 else float("nan"),
        "Q_per_A": (abs(Qe) / (Z + N)) if Z is not None and N is not None and (Z + N) > 0 else float("nan"),
    }

    if track_path and track_path.exists():
        tr = load_track(track_path)
        if tr:
            mL0, mL = tr[0]["massL"], tr[-1]["massL"]
            c_L = clip01((mL0 - mL) / max(mL0, 1e-12)) if mL0 > 1e-6 else 0.0
            pass_atom = pass_nuc and c_L <= 0.15 and (mL0 < 1e-6 or mL > 0.5 * mL0)
            out.update({
                "massL0": mL0, "massL_end": mL, "c_L": c_L,
                "PASS_atom": pass_atom,
                "PASS_park": pass_atom,  # alias for F16 scorecard name
            })
        else:
            out["PASS_atom"] = False
            out["error_track"] = "empty track"
    else:
        out["PASS_atom"] = None  # N/A nuclear-only

    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("diag")
    ap.add_argument("track", nargs="?", default=None)
    ap.add_argument("--Z", type=int, default=None)
    ap.add_argument("--N", type=int, default=None)
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()
    s = score(Path(args.diag), Path(args.track) if args.track else None, args.Z, args.N)
    if args.json:
        print(json.dumps(s, indent=2))
    else:
        print(f"PASS_nuc={s.get('PASS_nuc')} PASS_atom={s.get('PASS_atom')}")
        print(f"  c_Q_park={s.get('c_Q_park'):.4f} c_Qem_park={s.get('c_Qem_park'):.4f} "
              f"c_E={s.get('c_E'):.4f} c_G={s.get('c_G')}")
        print(f"  Q_phi {s.get('Q0'):.1f} → mid {s.get('Q_mid'):.1f} → end {s.get('Q_end'):.1f}")
        print(f"  Q_flux {s.get('Qf0'):.1f} → mid {s.get('Qf_mid'):.1f} → end {s.get('Qf_end'):.1f}")
        if s.get("Z"):
            print(f"  Qem/Z end = {s.get('Qem_per_Z'):.1f}")
        if "massL0" in s:
            print(f"  massL {s.get('massL0'):.1f} → {s.get('massL_end'):.1f} c_L={s.get('c_L'):.4f}")
    return 0 if s.get("ok") else 1


if __name__ == "__main__":
    raise SystemExit(main())
