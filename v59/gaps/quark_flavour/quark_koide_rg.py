#!/usr/bin/env python3
"""
quark_koide_rg.py  --  v59 gap attack, quark-flavour cluster (G5/G6).

HONEST RG / scheme analysis of the quark Koide ratios
    Q_d = sum m_i / (sum sqrt(m_i))^2     (down-type d,s,b)
    Q_u = ...                              (up-type u,c,t)

v59 targets (from t_N^2 = 1 - 14/D_N, Q = (1+2 t^2)/3):
    Q_d -> 11/15 = 0.733333...
    Q_u -> 23/27 = 0.851851...

The point of this script is NOT to celebrate a 0.3% match. It is to quantify
how much Q_d, Q_u MOVE when you change:
    (a) the renormalization scale mu (2 GeV vs M_Z vs m_t vs 1 TeV),
    (b) the input scheme (MS-bar running mass vs pole-ish / on-shell-like inputs),
    (c) the PDG uncertainty band on the masses.
and to compare that spread to the |target - central| gap. If the RG/scheme/PDG
spread is LARGER than the gap, then "0.3%" is not a prediction -- it is a number
that happens to land near a structural ratio at one particular (scale,scheme).

All mass numbers are documented MS-bar(mu) values compiled from PDG 2024 and
standard RG-evolution tables (e.g. Xing-Zhang-Zhou, Antusch-Maria style runs).
They are entered explicitly so the provenance and the scale are auditable.

NO claim of new precision is made; this is a falsification/honesty exercise.
"""

import numpy as np
import json
import itertools

# --------------------------------------------------------------------------
# v59 structural targets
# --------------------------------------------------------------------------
Q_d_target = 11/15          # 0.7333...
Q_u_target = 23/27          # 0.85185...
Q_lep_target = 2/3          # for reference

def koide_Q(masses):
    m = np.asarray(masses, dtype=float)
    return float(np.sum(m) / np.sum(np.sqrt(m))**2)

# --------------------------------------------------------------------------
# Mass datasets at several scales / schemes.
#
# Each entry: label -> dict with up=(u,c,t), down=(d,s,b), and a note.
# Masses in GeV. MS-bar masses m_q(mu) unless noted.
#
# Sources / construction:
#  * PDG2024 "current quark masses": light quarks quoted at mu=2 GeV;
#    m_c(m_c), m_b(m_b); m_t from direct (pole-like) measurement.
#    -> "PDG_native" mixes scales (the WORST but most-quoted convention).
#  * RG-evolved common-scale tables (MS-bar all run to one mu):
#    values below at 2 GeV, M_Z (91.19 GeV), m_t (~173 GeV), 1 TeV
#    are representative of the standard literature (Xing et al. 2008/2020,
#    Huang et al.); entered to ~few % accuracy which is the physical point.
# --------------------------------------------------------------------------

datasets = {
    # The most commonly quoted "textbook" numbers -- but each quark at a
    # DIFFERENT scale. This is the convention 11_quark_sector.py used.
    "PDG_native_mixed_scales": {
        "up":   [2.16e-3, 1.27, 172.69],   # u,c at their own ref; t pole
        "down": [4.67e-3, 93.4e-3, 4.18],  # d,s @2GeV; b @ m_b
        "note": "PDG quoted values, MIXED reference scales (u,d,s@2GeV; c@mc; b@mb; t pole). Not a single-scale set.",
    },
    # All MS-bar run to mu = 2 GeV
    "MSbar_2GeV": {
        "up":   [2.16e-3, 1.24, 384.0],    # m_u(2),  m_c(2)~1.24, m_t(2)~384 (large running up)
        "down": [4.67e-3, 93.4e-3, 7.4],   # m_d(2), m_s(2), m_b(2)~7.4 (b runs up from 4.18)
        "note": "MS-bar all evolved to 2 GeV (heavy-quark masses run UP at low scale).",
    },
    # All MS-bar run to mu = M_Z = 91.1876 GeV
    "MSbar_MZ": {
        "up":   [1.29e-3, 0.627, 171.7],   # m_u(MZ),m_c(MZ),m_t(MZ)
        "down": [2.75e-3, 0.0535, 2.89],   # m_d(MZ),m_s(MZ),m_b(MZ)
        "note": "MS-bar all evolved to M_Z=91.19 GeV (Xing-Zhang-Zhou style).",
    },
    # All MS-bar run to mu = m_t ~ 173.3 GeV
    "MSbar_mt": {
        "up":   [1.10e-3, 0.532, 162.5],   # m_u(mt),m_c(mt),m_t(mt)~162.5
        "down": [2.34e-3, 0.0455, 2.43],   # m_d(mt),m_s(mt),m_b(mt)
        "note": "MS-bar all evolved to m_t ~ 173 GeV.",
    },
    # All MS-bar run to mu = 1 TeV
    "MSbar_1TeV": {
        "up":   [0.96e-3, 0.466, 150.7],
        "down": [2.04e-3, 0.0398, 2.13],
        "note": "MS-bar all evolved to 1 TeV (GUT-direction trend).",
    },
}

# PDG uncertainty (1-sigma, GeV) for the *native* set, for a Monte-Carlo band.
# Light-quark errors are large fractionally; this dominates Q_d, Q_u spread.
pdg_unc = {
    "up":   [(0.07e-3), (0.02), (0.30)],      # u,c,t
    "down": [(0.07e-3), (0.86e-2), (0.025)],  # d,s,b  (s ~ +8.6/-?  ~ 8.6 MeV)
}

# --------------------------------------------------------------------------
# 1. Q at each scale / scheme
# --------------------------------------------------------------------------
print("=" * 78)
print("Quark Koide Q at several scales / schemes")
print("=" * 78)
print(f"  v59 targets:  Q_d = 11/15 = {Q_d_target:.6f}   Q_u = 23/27 = {Q_u_target:.6f}")
print()
hdr = f"  {'dataset':<26} {'Q_d':>9} {'gap_d%':>8} {'Q_u':>9} {'gap_u%':>8}"
print(hdr)
print("  " + "-" * (len(hdr) - 2))

table = {}
Qd_vals = []
Qu_vals = []
for label, d in datasets.items():
    Qd = koide_Q(d["down"])
    Qu = koide_Q(d["up"])
    gd = (Q_d_target - Qd) / Qd * 100
    gu = (Q_u_target - Qu) / Qu * 100
    Qd_vals.append(Qd)
    Qu_vals.append(Qu)
    table[label] = {"Q_d": Qd, "Q_u": Qu, "gap_d_pct": gd, "gap_u_pct": gu,
                    "note": d["note"]}
    print(f"  {label:<26} {Qd:9.5f} {gd:+8.3f} {Qu:9.5f} {gu:+8.3f}")

Qd_vals = np.array(Qd_vals)
Qu_vals = np.array(Qu_vals)

print()
print("  Spread across scales/schemes:")
print(f"    Q_d:  min {Qd_vals.min():.5f}  max {Qd_vals.max():.5f}  "
      f"range {Qd_vals.max()-Qd_vals.min():.5f}  "
      f"({(Qd_vals.max()-Qd_vals.min())/Qd_vals.mean()*100:.2f}% of mean)")
print(f"    Q_u:  min {Qu_vals.min():.5f}  max {Qu_vals.max():.5f}  "
      f"range {Qu_vals.max()-Qu_vals.min():.5f}  "
      f"({(Qu_vals.max()-Qu_vals.min())/Qu_vals.mean()*100:.2f}% of mean)")
print()
print(f"    |11/15 - Q_d| at PDG_native = {abs(Q_d_target-table['PDG_native_mixed_scales']['Q_d']):.5f}")
print(f"    |11/15 - Q_d| at M_Z        = {abs(Q_d_target-table['MSbar_MZ']['Q_d']):.5f}")
print(f"    => RG spread of Q_d ({Qd_vals.max()-Qd_vals.min():.4f}) "
      f"is {(Qd_vals.max()-Qd_vals.min())/abs(Q_d_target-table['PDG_native_mixed_scales']['Q_d']):.1f}x "
      f"the native gap.")
print(f"    => RG spread of Q_u ({Qu_vals.max()-Qu_vals.min():.4f}) "
      f"is {(Qu_vals.max()-Qu_vals.min())/abs(Q_u_target-table['PDG_native_mixed_scales']['Q_u']):.1f}x "
      f"the native gap.")

# --------------------------------------------------------------------------
# 2. PDG-uncertainty Monte Carlo on the native set
# --------------------------------------------------------------------------
print()
print("=" * 78)
print("PDG-uncertainty Monte Carlo (native mixed-scale set)")
print("=" * 78)
rng = np.random.default_rng(59)
N = 200000

def mc_band(central, unc):
    c = np.array(central); u = np.array(unc)
    # truncated-normal-ish: resample, clip to positive
    samples = rng.normal(c, u, size=(N, 3))
    samples = np.clip(samples, 1e-6, None)
    Qs = samples.sum(axis=1) / (np.sqrt(samples).sum(axis=1))**2
    return Qs

Qd_mc = mc_band(datasets["PDG_native_mixed_scales"]["down"], pdg_unc["down"])
Qu_mc = mc_band(datasets["PDG_native_mixed_scales"]["up"], pdg_unc["up"])

for nm, mc, tgt in [("Q_d", Qd_mc, Q_d_target), ("Q_u", Qu_mc, Q_u_target)]:
    lo, med, hi = np.percentile(mc, [16, 50, 84])
    inside = (tgt > lo) and (tgt < hi)
    nsig = abs(tgt - med) / ((hi - lo) / 2)
    print(f"  {nm}: median {med:.5f}  68%CI [{lo:.5f}, {hi:.5f}]  "
          f"target {tgt:.5f}  -> {nsig:.2f} sigma, "
          f"{'INSIDE' if inside else 'OUTSIDE'} 68% band")

# --------------------------------------------------------------------------
# 3. The cross-sector invariant (1 - t^2) * D = 14, restated at each scale.
#    This is the *cleaner* statement: it is just an algebraic recast of Q,
#    so it inherits exactly the same RG dependence. Show it explicitly.
# --------------------------------------------------------------------------
print()
print("=" * 78)
print("Cross-sector invariant  (1 - t^2)*D = dimG2 = 14  --  is it scale stable?")
print("=" * 78)
print("  t^2 = (3Q - 1)/2 ;  with D_d=35, D_u=63 the v59 claim is (1-t^2)*D = 14.")
print(f"  {'dataset':<26} {'(1-t2_d)*35':>12} {'(1-t2_u)*63':>12}")
print("  " + "-" * 52)
inv_d = []; inv_u = []
for label, d in datasets.items():
    Qd = koide_Q(d["down"]); Qu = koide_Q(d["up"])
    t2d = (3*Qd - 1)/2; t2u = (3*Qu - 1)/2
    vd = (1 - t2d) * 35; vu = (1 - t2u) * 63
    inv_d.append(vd); inv_u.append(vu)
    print(f"  {label:<26} {vd:12.3f} {vu:12.3f}")
print(f"  target = 14.000 for both.")
print(f"  Q_d invariant spans [{min(inv_d):.2f}, {max(inv_d):.2f}],  "
      f"Q_u invariant spans [{min(inv_u):.2f}, {max(inv_u):.2f}].")

# --------------------------------------------------------------------------
# 4. Overfitting control: how many simple rationals p/q (q<=27) land within
#    the RG band of Q_d, Q_u? If many do, "11/15" is not special.
# --------------------------------------------------------------------------
print()
print("=" * 78)
print("Overfitting control: simple rationals inside the RG band")
print("=" * 78)
for nm, vals, tgt in [("Q_d", Qd_vals, Q_d_target), ("Q_u", Qu_vals, Q_u_target)]:
    lo, hi = vals.min(), vals.max()
    hits = []
    for q in range(2, 28):
        for p in range(1, q):
            r = p / q
            if lo <= r <= hi:
                hits.append((p, q, r))
    # unique by value
    seen = set(); uniq = []
    for p, q, r in sorted(hits, key=lambda z: z[2]):
        key = round(r, 6)
        if key not in seen:
            seen.add(key); uniq.append((p, q, r))
    print(f"  {nm}: RG band [{lo:.4f}, {hi:.4f}] (width {hi-lo:.4f})")
    print(f"       {len(uniq)} distinct rationals p/q (q<=27) fall inside:")
    line = "       " + ", ".join(f"{p}/{q}" for p, q, r in uniq)
    print(line)
    print(f"       v59 target {('11/15' if nm=='Q_d' else '23/27')} is "
          f"{'among them' if any((p,q)==((11,15) if nm=='Q_d' else (23,27)) for p,q,_ in uniq) else 'NOT in band'}.")

# --------------------------------------------------------------------------
# Save
# --------------------------------------------------------------------------
out = {
    "targets": {"Q_d": Q_d_target, "Q_u": Q_u_target},
    "per_dataset": table,
    "Q_d_spread": {"min": float(Qd_vals.min()), "max": float(Qd_vals.max()),
                   "range": float(Qd_vals.max()-Qd_vals.min())},
    "Q_u_spread": {"min": float(Qu_vals.min()), "max": float(Qu_vals.max()),
                   "range": float(Qu_vals.max()-Qu_vals.min())},
    "mc_native": {
        "Q_d": {"median": float(np.median(Qd_mc)),
                "ci68": [float(np.percentile(Qd_mc,16)), float(np.percentile(Qd_mc,84))]},
        "Q_u": {"median": float(np.median(Qu_mc)),
                "ci68": [float(np.percentile(Qu_mc,16)), float(np.percentile(Qu_mc,84))]},
    },
}
with open("quark_koide_rg.json", "w") as f:
    json.dump(out, f, indent=2)
print()
print("Saved quark_koide_rg.json")
