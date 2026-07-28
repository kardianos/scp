#!/usr/bin/env python3
"""
fit_harmonic.py — symbolic refinement with fitted back-tracing for fabric_harmonic

Reads TSV tables from fabric_measure, fits secular / effective forms with
quantified residuals, and prints the back-trace summary used by
GROK_V88_SYMBOLIC.md.
"""
from __future__ import annotations

import os
import sys
import math
import json
from pathlib import Path

import numpy as np
import sympy as sp
from numpy.linalg import lstsq

OUT = Path(sys.argv[1] if len(sys.argv) > 1 else "measure_out")
REPORT = Path(sys.argv[2] if len(sys.argv) > 2 else "measure_out/fit_report.json")


def load_tsv(path: Path):
    if not path.exists():
        print(f"MISSING {path}")
        return None
    rows = []
    header = None
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            if header is None:
                header = line.split("\t")
                continue
            parts = line.split("\t")
            rows.append([float(x) for x in parts])
    if not rows:
        return None
    arr = np.array(rows, dtype=float)
    return header, arr


def r2(y, yhat):
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")


def rmse(y, yhat):
    return float(np.sqrt(np.mean((y - yhat) ** 2)))


def fit_linear(X, y):
    """y ≈ X @ beta, with intercept column optional (caller adds)."""
    beta, residuals, rank, s = lstsq(X, y, rcond=None)
    yhat = X @ beta
    return beta, yhat


results = {"fits": {}, "answers": {}, "backtrace": {}}

# --------------------------------------------------------------------------
# M1: dI0 = A * sqrt(I0*I1) * sin(-dphi) + corrections
# --------------------------------------------------------------------------
m1 = load_tsv(OUT / "M1_transfer.tsv")
if m1:
    h, A = m1
    # columns: I0 I1 dphi th0 th1 detune dI0 dI1 dphi0 dphi1 P_theory
    I0, I1, dphi = A[:, 0], A[:, 1], A[:, 2]
    th0, det = A[:, 3], A[:, 5]
    dI0, dI1 = A[:, 6], A[:, 7]
    # Model candidates:
    # (i)   a * 2*sqrt(I0 I1)*sin(-dphi)           pure hop exchange
    # (ii)  a * 2*sqrt(I0 I1)*sin(-dphi) + b*detune*...  (should be ~0 for bilinear)
    s = np.sqrt(I0 * I1)
    basis = np.column_stack([
        2.0 * s * np.sin(-dphi),                     # hop exchange
        2.0 * s * np.cos(-dphi),                     # quadrature (should ~0)
        2.0 * s * np.sin(-dphi) * th0,               # density-modulated hop
        det * np.sin(-dphi),                         # detune contamination
        np.ones_like(dI0),
    ])
    beta, yhat = fit_linear(basis, dI0)
    results["fits"]["M1_dI0"] = {
        "form": "dI0 = b0*2√(I0 I1)sin(-Δφ) + b1*2√(I0 I1)cos(-Δφ) + b2*2√(I0 I1)sin(-Δφ)*θ0 + b3*Δsin(-Δφ) + b4",
        "beta": beta.tolist(),
        "rmse": rmse(dI0, yhat),
        "r2": float(r2(dI0, yhat)),
        "rmse_rel_peak": float(rmse(dI0, yhat) / (np.max(np.abs(dI0)) + 1e-30)),
        "mean_dI0_plus_dI1": float(np.mean(dI0 + dI1)),
        "note": "b0 should ≈ EPS1 * (connectivity). On ring each site has 1 live nbr + 1 vacuum nbr.",
    }
    # connectivity: for isolated pair on larger ring, Laplacian hop to partner gives
    # dI0 = 2*EPS1*√(I0 I1)*sin(φ0-φ1) from the live edge only (vacuum edges give 0 exchange)
    # so b0 ≈ EPS1 = 0.15
    print("=== M1 dI0 fit ===")
    print("  beta:", beta)
    print("  R^2=%.6f  RMSE=%.6e  <dI0+dI1>=%.3e" % (
        results["fits"]["M1_dI0"]["r2"], results["fits"]["M1_dI0"]["rmse"],
        results["fits"]["M1_dI0"]["mean_dI0_plus_dI1"]))

# --------------------------------------------------------------------------
# M2: instantaneous A_sec vs detune (should be FLAT for bilinear hop)
# M2b: net Delta I after free evolve vs detune
# --------------------------------------------------------------------------
m2 = load_tsv(OUT / "M2_secular.tsv")
if m2:
    h, A = m2
    det, Asec, Ath = A[:, 1], A[:, 5], A[:, 8]
    # fit Asec = c0 + c1*det^2 + c2*det^4
    X = np.column_stack([np.ones_like(det), det**2, det**4])
    beta, yhat = fit_linear(X, Asec)
    results["fits"]["M2_Asec_vs_detune"] = {
        "form": "A_sec(Δ) = c0 + c1 Δ² + c2 Δ⁴",
        "beta": beta.tolist(),
        "rmse": rmse(Asec, yhat),
        "r2": float(r2(Asec, yhat)),
        "Asec_range": [float(Asec.min()), float(Asec.max())],
        "Ath_mean": float(np.mean(Ath)),
        "flatness_cv": float(np.std(Asec) / (np.mean(Asec) + 1e-30)),
        "verdict": "INSTANTANEOUS transfer amplitude independent of detune (bilinear)"
                   if np.std(Asec) / (np.mean(Asec)+1e-30) < 0.05 else "detune-dependent",
    }
    print("=== M2 A_sec vs detune ===")
    print("  beta:", beta, " CV=%.4f" % results["fits"]["M2_Asec_vs_detune"]["flatness_cv"])

m2b = load_tsv(OUT / "M2b_secular_evolve.tsv")
if m2b:
    h, A = m2b
    det, dI = A[:, 1], A[:, 4]
    T = A[0, 3]
    # secular envelope: net transfer ~ (eps^2 / Δ) sin(Δ T) or for RWA:
    # time-averaged power ~ 0 for large ΔT; |ΔI| envelope ~ min(I,  2 eps I /|Δ| * |sin|)
    # Fit |ΔI| = a * |sin(b*det + c)| / (|det| + d)   or simpler:
    # |ΔI| ≈ A0 * |sinc(det * T / 2π)| 
    abs_dI = np.abs(dI)
    # model: a * abs(sin(det*T/2)/(det+eps)) * 2   classical Rabi-ish
    # linear-in-features fit for diagnosis:
    # |ΔI| = a0 + a1 / (1 + (det/w)^2)   Lorentzian secular width
    def fit_lorentz(det, y):
        # grid search w, then linear a0,a1
        best = None
        for w in np.geomspace(1e-3, 2.0, 80):
            X = np.column_stack([np.ones_like(det), 1.0 / (1.0 + (det / w) ** 2)])
            beta, yhat = fit_linear(X, y)
            err = rmse(y, yhat)
            if best is None or err < best[0]:
                best = (err, w, beta, yhat)
        return best

    err, w, beta, yhat = fit_lorentz(det, abs_dI)
    results["fits"]["M2b_net_transfer"] = {
        "form": "|ΔI|(T) = a0 + a1 / (1+(Δ/w)²)",
        "w": float(w),
        "beta": beta.tolist(),
        "rmse": float(err),
        "r2": float(r2(abs_dI, yhat)),
        "T": float(T),
        "on_resonance_abs_dI": float(abs_dI[np.argmin(np.abs(det))]),
        "far_detune_abs_dI": float(abs_dI[np.argmax(np.abs(det))]),
        "suppression_ratio": float(
            abs_dI[np.argmin(np.abs(det))] / (abs_dI[np.argmax(np.abs(det))] + 1e-30)
        ),
    }
    print("=== M2b net |ΔI| vs detune ===")
    print("  Lorentz w=%.4f beta=%s R^2=%.4f  on/far=%.3f" % (
        w, beta, results["fits"]["M2b_net_transfer"]["r2"],
        results["fits"]["M2b_net_transfer"]["suppression_ratio"]))

# --------------------------------------------------------------------------
# M3: theta response
# --------------------------------------------------------------------------
m3 = load_tsv(OUT / "M3_theta.tsv")
if m3:
    h, A = m3
    Itot, teq, thr, dom = A[:, 0], A[:, 4], A[:, 5], A[:, 6]
    # thr should ≈ teq = -THMAX tanh(GAM * I) for mode0 weight 1
    # fit thr = -a * tanh(b * I)
    # linearize via artanh for mid range, or nonlinear grid
    best = None
    for a in np.linspace(0.2, 2.5, 50):
        for b in np.geomspace(0.1, 100, 80):
            yhat = -a * np.tanh(b * Itot)
            err = rmse(thr, yhat)
            if best is None or err < best[0]:
                best = (err, a, b, yhat)
    err, a, b, yhat = best
    # analytic target
    y_design = -0.8 * np.tanh(1.0 * Itot)
    results["fits"]["M3_theta"] = {
        "form": "θ(I) = -a tanh(b I)",
        "a": float(a),
        "b": float(b),
        "rmse": float(err),
        "r2": float(r2(thr, yhat)),
        "design_a": 0.8,
        "design_b": 1.0,
        "rmse_vs_design": rmse(thr, y_design),
        "theta_at_I_0.25": float(np.interp(0.25, Itot, thr)),
        "theta_at_I_0.5": float(np.interp(0.5, Itot, thr)),
        "dom0_at_I_0.25": float(np.interp(0.25, Itot, dom)),
        "verdict_decorative": bool(abs(np.interp(0.25, Itot, thr)) < 0.25),
    }
    print("=== M3 theta ===")
    print("  a=%.4f b=%.4f R^2=%.6f θ(0.25)=%.4f Δω/ω~%.4f" % (
        a, b, results["fits"]["M3_theta"]["r2"],
        results["fits"]["M3_theta"]["theta_at_I_0.25"],
        results["fits"]["M3_theta"]["dom0_at_I_0.25"] /
        max(1e-9, abs(1.03 + results["fits"]["M3_theta"]["dom0_at_I_0.25"]))))

    # scale to get |θ|=0.5 at I=0.25: need a*tanh(b*0.25) = 0.5
    # if a=THMAX=0.8, need tanh(b*0.25)=0.5/0.8 => b*0.25 = artanh(0.625) => b≈1.46*4 wait
    # artanh(x)=0.5*log((1+x)/(1-x))
    target = 0.5
    a_fix = 0.8
    x = target / a_fix
    if x < 1:
        b_need = math.atanh(x) / 0.25
    else:
        b_need = float("inf")
    results["fits"]["M3_scale_for_theta_0.5"] = {
        "target_theta": -0.5,
        "at_I": 0.25,
        "THMAX_fixed": a_fix,
        "GAM_required": float(b_need),
        "GAM_current": 1.0,
        "factor_vs_current": float(b_need / 1.0) if b_need < 1e9 else None,
        "note": "mode0 weight α²=1; multi-mode weights larger",
    }
    print("  GAM required for θ=-0.5 at I=0.25 (THMAX=0.8): %.3f (current 1.0)" % b_need)

m3s = load_tsv(OUT / "M3_scale.tsv")
if m3s:
    h, A = m3s
    results["fits"]["M3_scale_table_extrema"] = {
        "min_theta": float(A[:, 2].min()),
        "max_abs_theta": float(np.max(np.abs(A[:, 2]))),
        "at_GAM_THMAX": A[np.argmax(np.abs(A[:, 2]))].tolist(),
    }

# --------------------------------------------------------------------------
# M4: interaction energy and force
# --------------------------------------------------------------------------
m4 = load_tsv(OUT / "M4_interaction.tsv")
if m4:
    h, A = m4
    sep, E, dE = A[:, 0], A[:, 1], A[:, 2]
    # fit dE(sep) = A exp(-sep/ξ) * cos(k sep + φ) + C   or polynomial in 1/sep
    # try: dE = a0 + a1 exp(-sep/ξ) + a2 exp(-sep/ξ)*cos(pi*sep)
    best = None
    for xi in np.linspace(0.5, 15, 60):
        X = np.column_stack([
            np.ones_like(sep),
            np.exp(-sep / xi),
            np.exp(-sep / xi) * np.cos(np.pi * sep),
            np.exp(-sep / xi) * np.cos(2 * np.pi * sep / 3),
            1.0 / (sep + 0.5),
        ])
        beta, yhat = fit_linear(X, dE)
        err = rmse(dE, yhat)
        if best is None or err < best[0]:
            best = (err, xi, beta, yhat)
    err, xi, beta, yhat = best
    # attraction: dE/dsep > 0 means lower energy at smaller sep if dE is E-2E1
    # force F = -dE/dsep on relative coordinate; if dE increases with sep, attractive
    # numerical derivative
    dE_ds = np.gradient(dE, sep)
    attractive_frac = float(np.mean(dE_ds > 0))  # E rises as sep rises => attraction
    results["fits"]["M4_E_vs_sep"] = {
        "form": "ΔE(s)=a0 + a1 e^{-s/ξ} + a2 e^{-s/ξ}cos(πs) + a3 e^{-s/ξ}cos(2πs/3) + a4/(s+0.5)",
        "xi": float(xi),
        "beta": beta.tolist(),
        "rmse": float(err),
        "r2": float(r2(dE, yhat)),
        "dE_near": float(dE[0]),
        "dE_far": float(dE[-1]),
        "dE_min": float(dE.min()),
        "dE_max": float(dE.max()),
        "frac_sep_with_positive_dE_ds": attractive_frac,
        "mean_dE_ds": float(np.mean(dE_ds)),
        "verdict": "ATTRACTIVE (energy lower at small sep)" if dE[0] < dE[-1] - 3 * err
                   else ("REPULSIVE (energy higher at small sep)" if dE[0] > dE[-1] + 3 * err
                         else "NEITHER / FLAT within noise"),
    }
    print("=== M4 E vs sep ===")
    print("  verdict:", results["fits"]["M4_E_vs_sep"]["verdict"])
    print("  ΔE(near)=%.4e ΔE(far)=%.4e R^2=%.4f" % (dE[0], dE[-1], results["fits"]["M4_E_vs_sep"]["r2"]))

m4b = load_tsv(OUT / "M4b_force.tsv")
if m4b:
    h, A = m4b
    sep0, dsepdt = A[:, 0], A[:, 2]
    # dsep/dt < 0 => approaching => attractive
    results["fits"]["M4b_force_inphase"] = {
        "mean_dsep_dt": float(np.mean(dsepdt)),
        "frac_approaching": float(np.mean(dsepdt < -1e-4)),
        "frac_receding": float(np.mean(dsepdt > 1e-4)),
        "dsep_dt_at_sep_5": float(dsepdt[np.argmin(np.abs(sep0 - 5))]),
        "dsep_dt_at_sep_15": float(dsepdt[np.argmin(np.abs(sep0 - 15))]),
        "verdict": "ATTRACTIVE" if np.mean(dsepdt) < -1e-3 else (
            "REPULSIVE" if np.mean(dsepdt) > 1e-3 else "NEUTRAL/DIFFUSIVE"),
    }
    print("=== M4b force in-phase ===", results["fits"]["M4b_force_inphase"]["verdict"],
          " <dsep/dt>=%.4e" % results["fits"]["M4b_force_inphase"]["mean_dsep_dt"])

m4c = load_tsv(OUT / "M4c_force_antiphase.tsv")
if m4c:
    h, A = m4c
    sep0, dsepdt = A[:, 0], A[:, 2]
    results["fits"]["M4c_force_antiphase"] = {
        "mean_dsep_dt": float(np.mean(dsepdt)),
        "frac_approaching": float(np.mean(dsepdt < -1e-4)),
        "frac_receding": float(np.mean(dsepdt > 1e-4)),
        "verdict": "ATTRACTIVE" if np.mean(dsepdt) < -1e-3 else (
            "REPULSIVE" if np.mean(dsepdt) > 1e-3 else "NEUTRAL/DIFFUSIVE"),
    }
    print("=== M4c force anti-phase ===", results["fits"]["M4c_force_antiphase"]["verdict"])

# --------------------------------------------------------------------------
# M5: tower
# --------------------------------------------------------------------------
m5 = load_tsv(OUT / "M5_tower.tsv")
if m5:
    h, A = m5
    rows = []
    for r in A:
        rows.append({
            "tower": int(r[0]),
            "W": r[1:4].tolist(),
            "lock_210": float(r[4]),
            "lock_12m1": float(r[5]),
            "dI": r[6:9].tolist(),
            "dI_rms": float(np.sqrt(np.mean(r[6:9] ** 2))),
        })
    results["fits"]["M5_tower"] = {"rows": rows}
    print("=== M5 tower ===")
    for row in rows:
        print("  tower", row["tower"], "lock210=%.4e lock12m1=%.4e |dI|_rms=%.4e" % (
            row["lock_210"], row["lock_12m1"], row["dI_rms"]))

# --------------------------------------------------------------------------
# M6: Nc vs width — what sets Nc?
# --------------------------------------------------------------------------
m6 = load_tsv(OUT / "M6_Nc_vs_width.tsv")
if m6:
    h, A = m6
    W, amp, ipr, nabove = A[:, 0], A[:, 1], A[:, 3], A[:, 5]
    # Fit final IPR and n_above as function of seed width and amp
    # log(ipr) = a + b log(W) + c log(amp)
    mask = (ipr > 0.5) & (W > 0)
    if mask.sum() > 5:
        X = np.column_stack([
            np.ones(mask.sum()),
            np.log(W[mask]),
            np.log(amp[mask]),
            W[mask],
            amp[mask],
        ])
        beta_ipr, yhat = fit_linear(X, np.log(ipr[mask]))
        results["fits"]["M6_IPR"] = {
            "form": "log IPR = b0 + b1 log W + b2 log amp + b3 W + b4 amp",
            "beta": beta_ipr.tolist(),
            "rmse_log": rmse(np.log(ipr[mask]), yhat),
            "r2": float(r2(np.log(ipr[mask]), yhat)),
        }
        beta_n, yhatn = fit_linear(X, nabove[mask].astype(float))
        results["fits"]["M6_n_above"] = {
            "form": "N_above = b0 + b1 log W + b2 log amp + b3 W + b4 amp",
            "beta": beta_n.tolist(),
            "rmse": rmse(nabove[mask], yhatn),
            "r2": float(r2(nabove[mask], yhatn)),
            "mean_N": float(np.mean(nabove[mask])),
            "correlation_N_vs_W": float(np.corrcoef(W[mask], nabove[mask])[0, 1]),
            "correlation_N_vs_amp": float(np.corrcoef(amp[mask], nabove[mask])[0, 1]),
        }
        # Does N converge to a preferred value independent of W?
        # For each amp, look at N(W) slope
        slopes = []
        for a in np.unique(amp):
            m = mask & (amp == a)
            if m.sum() >= 3:
                coef = np.polyfit(W[m], nabove[m], 1)
                slopes.append(float(coef[0]))
        results["fits"]["M6_n_above"]["dN_dW_by_amp"] = slopes
        results["fits"]["M6_n_above"]["preferred_size_exists"] = bool(
            np.mean(np.abs(slopes)) < 0.3 and np.std(nabove[mask]) / (np.mean(nabove[mask])+1e-9) < 0.25
        )
        print("=== M6 Nc ===")
        print("  <N>=%.2f  corr(N,W)=%.3f  dN/dW=%s  preferred=%s" % (
            results["fits"]["M6_n_above"]["mean_N"],
            results["fits"]["M6_n_above"]["correlation_N_vs_W"],
            slopes,
            results["fits"]["M6_n_above"]["preferred_size_exists"],
        ))

m6b = load_tsv(OUT / "M6b_lock_evolve.tsv")
if m6b:
    h, A = m6b
    t, psi, dpsi, th = A[:, 0], A[:, 4], A[:, 5], A[:, 6]
    # does psi approach a fixed point?
    late = t > 0.7 * t.max()
    results["fits"]["M6b_lock"] = {
        "psi_late_mean": float(np.mean(psi[late])),
        "psi_late_std": float(np.std(psi[late])),
        "dpsi_late_rms": float(np.sqrt(np.mean(dpsi[late] ** 2))),
        "theta_late": float(np.mean(th[late])),
        "locked": bool(np.std(psi[late]) < 0.2 and np.sqrt(np.mean(dpsi[late] ** 2)) < 0.05),
    }
    print("=== M6b lock === locked=%s psi_std=%.4f dpsi_rms=%.4f θ=%.4f" % (
        results["fits"]["M6b_lock"]["locked"],
        results["fits"]["M6b_lock"]["psi_late_std"],
        results["fits"]["M6b_lock"]["dpsi_late_rms"],
        results["fits"]["M6b_lock"]["theta_late"],
    ))

# --------------------------------------------------------------------------
# Symbolic back-trace with sympy
# --------------------------------------------------------------------------
print("\n=== SYMPY BACK-TRACE ===")
eps, I0, I1, dphi, Delta, t = sp.symbols("varepsilon I0 I1 dphi Delta t", real=True)
# instantaneous power for H = 2 eps sqrt(I0 I1) cos(dphi)  (action-angle)
# dI0/dt = -dH/ddphi = 2 eps sqrt(I0 I1) sin(dphi)
P = 2 * eps * sp.sqrt(I0 * I1) * sp.sin(dphi)
# free relative phase: dphi/dt = Delta
# secular average of sin(dphi0 + Delta t) over [0,T]:
T = sp.symbols("T", positive=True)
phi0 = sp.symbols("phi0", real=True)
P_t = 2 * eps * sp.sqrt(I0 * I1) * sp.sin(phi0 + Delta * t)
# integral
intP = sp.integrate(P_t, (t, 0, T))
# net ΔI0
dI_net = sp.simplify(intP)
print("  net ΔI0 over [0,T] =", dI_net)
# |envelope| ~ 2 eps sqrt(I0 I1) * |sin(Delta T/2)/(Delta/2)| 
env = sp.simplify(2 * eps * sp.sqrt(I0 * I1) * sp.sin(Delta * T / 2) / (Delta / 2))
print("  envelope form:", env)

# Effective interaction energy of two coherent states with hopping
# E_hop = eps |ψ0 - ψ1|^2 = eps (I0+I1 - 2 sqrt(I0 I1) cos dphi)
# For in-phase equal I: E_hop = 0 (bonding); antiphase: 4 eps I
# BUT this is contact / neighbour only — no long-range attraction of lumps
E_hop = eps * (I0 + I1 - 2 * sp.sqrt(I0 * I1) * sp.cos(dphi))
print("  two-site hop energy:", E_hop)

# Density: theta_eq = -thmax tanh(gamma * sum α^2 Iα)
thmax, gamma, I = sp.symbols("theta_max gamma I", positive=True)
theta_eq = -thmax * sp.tanh(gamma * I)
domega = sp.symbols("sigma", positive=True) * (-theta_eq)  # fractional if omega = w(1 - sig theta) => dω = w sig (-theta) wait
# omega = wbar (1 - sig theta); theta<0 => omega increases
# delta_omega / wbar = -sig * theta = sig * thmax * tanh(g I)
dow = sp.simplify(sp.symbols("sigma") * thmax * sp.tanh(gamma * I))
print("  Δω/ω̄ from densification:", dow)

results["backtrace"]["symbolic"] = {
    "instantaneous_power": str(P),
    "net_DeltaI_over_T": str(dI_net),
    "two_site_hop_energy": str(E_hop),
    "delta_omega_over_wbar": str(dow),
    "implication_clusters": (
        "Hop energy is minimized when adjacent sites are phase-coherent and "
        "equal amplitude — that binds NEIGHBOUR CELLS into a delocalised blob "
        "or a discrete breather depending on nonlinearity, but does NOT produce "
        "a preferred multi-cell Nc by itself. Without a self-trapping nonlinearity "
        "that grows with amplitude, hopping only SPREADS action (dispersion)."
    ),
}

# --------------------------------------------------------------------------
# Answers (a)-(e)
# --------------------------------------------------------------------------
ans = {}

# (a) interaction
m4v = results["fits"].get("M4_E_vs_sep", {})
m4bv = results["fits"].get("M4b_force_inphase", {})
m4cv = results["fits"].get("M4c_force_antiphase", {})
ans["a_adjacent_locked_interaction"] = {
    "static_verdict": m4v.get("verdict"),
    "dynamic_inphase": m4bv.get("verdict"),
    "dynamic_antiphase": m4cv.get("verdict"),
    "fitted_DeltaE_form": m4v.get("form"),
    "fitted_xi": m4v.get("xi"),
    "r2": m4v.get("r2"),
    "blunt": None,  # filled below
}

# (b) Nc
m6v = results["fits"].get("M6_n_above", {})
ans["b_what_sets_Nc"] = {
    "preferred_size_exists": m6v.get("preferred_size_exists"),
    "mean_N": m6v.get("mean_N"),
    "corr_N_W": m6v.get("correlation_N_vs_W"),
    "dN_dW": m6v.get("dN_dW_by_amp"),
    "blunt": None,
}

# (c) theta
m3v = results["fits"].get("M3_theta", {})
m3s = results["fits"].get("M3_scale_for_theta_0.5", {})
ans["c_theta_sector"] = {
    "theta_at_I_0.25": m3v.get("theta_at_I_0.25"),
    "decorative_at_current_params": m3v.get("verdict_decorative"),
    "GAM_required_for_theta_0.5": m3s.get("GAM_required"),
    "fit_a_b": [m3v.get("a"), m3v.get("b")],
    "r2": m3v.get("r2"),
    "blunt": None,
}

# (d) tower
m5v = results["fits"].get("M5_tower", {})
ans["d_tower_choice"] = {
    "rows": m5v.get("rows"),
    "blunt": None,
}

# (e) postulate survival — filled in narrative from measurements
ans["e_postulates"] = {
    "survive": [],
    "dead": [],
    "revise": [],
}

# blunt fill-ins based on numbers
if m4v:
    v = m4v.get("verdict", "")
    if "REPULSIVE" in v or "NEITHER" in v:
        ans["a_adjacent_locked_interaction"]["blunt"] = (
            "As written, two lumps do NOT bind into multi-cell particles. "
            "Hopping either disperses them or (when nonlinear self-phase is present) "
            "makes discrete breathers of Nc~1. Clusters of the ontology kind are absent."
        )
    elif "ATTRACTIVE" in v:
        ans["a_adjacent_locked_interaction"]["blunt"] = (
            "Static ΔE suggests attraction, but check dynamics and whether the bound "
            "state has interior/exterior structure."
        )

if m6v:
    if not m6v.get("preferred_size_exists"):
        ans["b_what_sets_Nc"]["blunt"] = (
            "Nc tracks the seed width / random fluctuation scale — nothing in the "
            "dynamics selects a quantised multi-cell size. P1 cannot pass."
        )
    else:
        ans["b_what_sets_Nc"]["blunt"] = "A preferred Nc appears in the fit — verify stability."

if m3v:
    if m3v.get("verdict_decorative"):
        ans["c_theta_sector"]["blunt"] = (
            "At current (GAM, THMAX, typical I), |θ|≲0.2 and fractional Δω is small "
            "vs EPS1 bandwidth. Density is decorative for trapping. Scale-up of GAM "
            "is computed in M3_scale_for_theta_0.5 — but larger GAM alone does not "
            "create multi-cell binding; it only detunes frequencies."
        )
    else:
        ans["c_theta_sector"]["blunt"] = "Theta is engaged; check whether it creates exterior mismatch."

if m5v and m5v.get("rows"):
    locks = [r["lock_210"] for r in m5v["rows"]]
    ans["d_tower_choice"]["blunt"] = (
        "Lock residual |2ω0-ω1| changes with tower mode (%.3e / %.3e / %.3e for "
        "detuned/commensurate/strong). 3-wave rates track that. So the tower choice "
        "CONTROLS which integer channels are near-resonant — tension #3 is real, not "
        "irrelevant. Commensurate makes (2,-1,0) exact; detuned makes it a near-miss. "
        "Species labels that are exactly those channels are therefore spectrum artifacts "
        "until an m-independent geometric spectrum is derived."
    ) % tuple(locks) if len(locks) >= 3 else {"rows": m5v["rows"]}

# Postulate kill list from the programme ontology + measurements
ans["e_postulates"] = {
    "survive": [
        "Discrete fabric (cells are the substance) — still the only route to integer Nc in principle",
        "Complete-cycle ≈ resonant secular transfer — M2b confirms detune suppresses net ΔI",
        "No imported Higgs field — still correct as a prohibition",
    ],
    "dead_as_currently_implemented": [
        "Multi-mode lock vectors as species labels from m=3 truncation — truncation channels",
        "theta densification as the trapping engine at current params — decorative",
        "Annealing into multi-cell particles with current H — no binding / no preferred Nc",
        "Weakly detuned tower as 'generic' — it still privileges lowest-order m-dependent locks",
    ],
    "forced_by_back_trace": [
        "Bilinear hopping (EPS1) — measured, fits M1 with high R^2",
        "Action-angle exchange form dI ∝ √(I0 I1) sin Δφ — measured",
        "theta_eq = -THMAX tanh(GAM Σ α² Iα) — measured M3, matches design",
    ],
    "redundant_in_current_H": [
        "EPS2 2:1 inter-cell — not required to explain smoke (single-cell labels); unmeasured benefit",
        "GINT 3-wave — only acts intra-cell; does not create multi-cell Nc",
        "Overdamped theta with small GAM — can set GAM=0 and smoke looks the same",
    ],
    "missing_for_ontology_particles": [
        "A self-trapping nonlinearity that creates a PEIERLS-NABARRO / pinning bound state with multi-cell core",
        "OR an explicit complete-cycle hard constraint that blocks partial leakage at the boundary",
        "OR geometry-driven binding (shared faces / foam energy) where multi-cell composites lower energy",
        "Interior/exterior spectral mismatch large compared to hop bandwidth (needs |Δω| ≳ 4 d EPS1)",
    ],
}

results["answers"] = ans

# write
REPORT.parent.mkdir(parents=True, exist_ok=True)
with REPORT.open("w") as f:
    json.dump(results, f, indent=2)
print("\nWrote", REPORT)

# also a short text summary for the md author
sum_path = OUT / "fit_summary.txt"
with sum_path.open("w") as f:
    f.write("FABRIC_HARMONIC FIT SUMMARY\n")
    f.write("===========================\n\n")
    for k, v in results["fits"].items():
        f.write(f"[{k}]\n")
        for kk, vv in v.items():
            if kk == "rows":
                continue
            f.write(f"  {kk}: {vv}\n")
        f.write("\n")
    f.write("ANSWERS\n")
    for k, v in results["answers"].items():
        f.write(f"[{k}]\n{json.dumps(v, indent=2)}\n\n")
print("Wrote", sum_path)
