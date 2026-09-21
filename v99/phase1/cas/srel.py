#!/usr/bin/env python3
"""v99 Phase 1 — relative entropy of a W=1 wrap vs vacuum on C_n.

Paper instrument, not a kernel. Free real scalar on the cycle graph,
ground-state covariance, coherent displacement = classical winding.

S_rel on a proper subsystem A (the door). Globally two pures give
S_rel = 0 or ∞; that is why A is a proper subset (Dorau: wedge algebra).

Run: python3 v99/phase1/cas/srel.py
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent / "out"
OUT.mkdir(exist_ok=True)

NS = list(range(3, 25))
MS = (0.0, 0.1, 1.0)


def omega_k(n: int, k: int, m: float) -> float:
    lam = 4.0 * math.sin(math.pi * k / n) ** 2
    return math.sqrt(lam + m * m)


def correlators(n: int, m: float) -> tuple[np.ndarray, np.ndarray]:
    """⟨φ_i φ_j⟩ and ⟨π_i π_j⟩ of the Gaussian vacuum on C_n.

    k=0, m=0 is the IR zero mode and is omitted (quotient by constants).
    A W=1 wrap is already orthogonal to the constant mode on the full cycle.
    """
    X = np.zeros((n, n))
    P = np.zeros((n, n))
    for k in range(n):
        om = omega_k(n, k, m)
        if om < 1e-14:
            continue
        ks = np.arange(n)
        # real cosine transform of a circulant mode
        # X_ij += cos(2π k (i-j)/n) / (2 ω) / n
        phases = 2.0 * math.pi * k / n * ks
        c = np.cos(phases)
        s = np.sin(phases)
        # outer circulant via real/imag of e^{ikx}
        vec = c + 1j * s
        outer = np.real(np.outer(np.conj(vec), vec))
        X += outer / (2.0 * om * n)
        P += outer * (om / (2.0 * n))
    return X, P


def wrap_mean(n: int, A: float, traveling: bool, m: float):
    """Classical W=1 wrap. Two chiral planes: φ=A cos θ, χ=A sin θ.

    Static: π=0. Traveling / on-shell: π_φ = A ω sin θ, π_χ = -A ω cos θ
    so the complex combination is a single k=+1 coherent state.
    """
    theta = 2.0 * math.pi * np.arange(n) / n
    phi = A * np.cos(theta)
    chi = A * np.sin(theta)
    if not traveling:
        return phi, chi, np.zeros(n), np.zeros(n)
    om = omega_k(n, 1, m)
    pi_phi = A * om * np.sin(theta)
    pi_chi = -A * om * np.cos(theta)
    return phi, chi, pi_phi, pi_chi


def srel_block(mean: np.ndarray, cov: np.ndarray) -> float:
    """(1/2) μ^T C^{-1} μ, C the restricted covariance of one quadrature set."""
    if mean.size == 0:
        return 0.0
    # Symmetric positive: solve rather than invert.
    # Tiny ridge if the IR skip left a near-kernel.
    C = cov + 1e-12 * np.eye(cov.shape[0])
    try:
        y = np.linalg.solve(C, mean)
    except np.linalg.LinAlgError:
        y = np.linalg.pinv(C) @ mean
    val = 0.5 * float(mean @ y)
    return val if val >= 0.0 else float("nan")


def srel_door(n: int, m: float, ell: int, A: float, traveling: bool) -> float:
    """S_rel(vacuum || wrap) on a contiguous arc of `ell` sites."""
    X, P = correlators(n, m)
    phi, chi, pi_phi, pi_chi = wrap_mean(n, A, traveling, m)
    Aidx = slice(0, ell)
    XA, PA = X[Aidx, Aidx], P[Aidx, Aidx]
    s = 0.0
    s += srel_block(phi[Aidx], XA)
    s += srel_block(chi[Aidx], XA)
    if traveling:
        s += srel_block(pi_phi[Aidx], PA)
        s += srel_block(pi_chi[Aidx], PA)
    return s


def trough_Ew(n: int, A: float = 1.0) -> float:
    """TROUGH.md E_w for a complex wrap of amplitude A per site."""
    return 4.0 * n * (math.sin(math.pi / n) ** 2) * (A * A)


def omega1(n: int, m: float = 0.0) -> float:
    return omega_k(n, 1, m)


def atom_eps(n: int, m: float = 0.0) -> float:
    """ε = ω/2π at the wrap's pitch, A₀=1 lattice units."""
    return omega1(n, m) / (2.0 * math.pi)


def casini_flux(n: int, ell: int, A: float = 1.0) -> float:
    """Modular-weighted energy flux of a uniform wrap on an interval.

    β(x) = 2π x (ℓ-x)/ℓ   (Casini interval kernel, lattice sites 1..ell-1).
    Site energy density of |ψ|=A wrap: λ_1 A²  (Dirichlet per site).
    """
    lam = 4.0 * math.sin(math.pi / n) ** 2 * (A * A)
    acc = 0.0
    for j in range(1, ell):  # endpoints have β=0
        beta = 2.0 * math.pi * j * (ell - j) / ell
        acc += beta * lam
    return acc


def amp_unit_L2(n: int) -> float:
    """|ψ_j|=A with Σ_j |ψ_j|² = 1  ⇒  A = 1/√n  (complex)."""
    return 1.0 / math.sqrt(n)


def amp_one_quantum(n: int, m: float) -> float:
    """On-shell energy = ω (one quantum, ħ=1).

    Complex wrap |ψ|=A: TROUGH E_w = n λ A². Identify on-shell E with E_w
    (both quadratures; virial). One quantum: n λ A² = ω, A² = ω/(n λ).
    For m=0, ω=√λ, A² = 1/(n √λ) * 1/√λ wait: ω/(n λ) = 1/(n ω).
    """
    om = omega1(n, m)
    lam = om * om if m == 0.0 else 4.0 * math.sin(math.pi / n) ** 2 + m * m
    # E_w = n * (4 sin²) A² ; for m>0 use n * λ_m A² with λ_m = ω²
    lam = om * om  # ω² = λ + m²; Dirichlet is λ A² n not ω².
    # Use Dirichlet of the spatial wrap: n * 4 sin²(π/n) A² = ω
    spatial = 4.0 * math.sin(math.pi / n) ** 2
    if spatial < 1e-18:
        return 0.0
    return math.sqrt(om / (n * spatial))


def argmin_int(ns, vals):
    pairs = [(n, v) for n, v in zip(ns, vals) if v == v and math.isfinite(v)]
    if not pairs:
        return None, None
    m = min(v for _, v in pairs)
    best = [n for n, v in pairs if abs(v - m) < 1e-9]
    return best, m


def fmt(x: float) -> str:
    if x != x or not math.isfinite(x):
        return "nan"
    return f"{x:.6g}"


def dump_table(path: Path, header: list[str], rows: list[list]):
    with path.open("w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(c) for c in r) + "\n")


def main() -> None:
    lines: list[str] = []

    def log(msg: str = "") -> None:
        print(msg)
        lines.append(msg)

    log("=== v99 srel.py — W=1 wrap vs vacuum on C_n ===")
    log("A proper subsystem is required (global pures: S_rel 0 or inf).")
    log("")

    # --- Dorau analog: S_rel vs Casini-weighted flux, half-cycle, A=1 ---
    log("=== D1: S_rel (half-cycle, A=1, static) vs Casini flux; m=0.1 ===")
    log(f"{'n':>4} {'ell':>4} {'S_rel':>12} {'K_casini':>12} {'S/K':>10} {'E_w':>10}")
    d1_rows = []
    for n in NS:
        if n % 2:
            continue
        ell = n // 2
        s = srel_door(n, 0.1, ell, 1.0, traveling=False)
        k = casini_flux(n, ell, 1.0)
        ew = trough_Ew(n)
        ratio = s / k if abs(k) > 1e-18 else float("nan")
        log(f"{n:4d} {ell:4d} {s:12.5g} {k:12.5g} {ratio:10.4g} {ew:10.5g}")
        d1_rows.append([n, ell, fmt(s), fmt(k), fmt(ratio), fmt(ew)])
    dump_table(OUT / "d1_half_casini.tsv",
               ["n", "ell", "S_rel", "K_casini", "S_over_K", "E_w"], d1_rows)

    log("")
    log("=== D2: door size (A=1, static, m=0.1)  1-site / 2-site / half ===")
    log(f"{'n':>4} {'S_1':>12} {'S_2':>12} {'S_half':>12} {'E_w':>10} {'S_1/E_w':>10}")
    d2_rows = []
    for n in NS:
        s1 = srel_door(n, 0.1, 1, 1.0, False)
        s2 = srel_door(n, 0.1, min(2, n - 1), 1.0, False)
        sh = srel_door(n, 0.1, max(1, n // 2), 1.0, False)
        ew = trough_Ew(n)
        log(f"{n:4d} {s1:12.5g} {s2:12.5g} {sh:12.5g} {ew:10.5g} {s1/ew:10.4g}")
        d2_rows.append([n, fmt(s1), fmt(s2), fmt(sh), fmt(ew), fmt(s1 / ew)])
    dump_table(OUT / "d2_doors_A1.tsv",
               ["n", "S_1", "S_2", "S_half", "E_w", "S1_over_Ew"], d2_rows)

    log("")
    log("=== D3: one-quantum amplitude (E_w = ω), S_rel on 1-site and half, m=0.1 ===")
    log(f"{'n':>4} {'A':>10} {'S_1':>12} {'S_half':>12} {'E_w':>10} {'ω':>10}")
    d3_rows = []
    for n in NS:
        A = amp_one_quantum(n, 0.1)
        s1 = srel_door(n, 0.1, 1, A, False)
        sh = srel_door(n, 0.1, max(1, n // 2), A, False)
        ew = trough_Ew(n, A)
        om = omega1(n, 0.1)
        log(f"{n:4d} {A:10.5g} {s1:12.5g} {sh:12.5g} {ew:10.5g} {om:10.5g}")
        d3_rows.append([n, fmt(A), fmt(s1), fmt(sh), fmt(ew), fmt(om)])
    dump_table(OUT / "d3_one_quantum.tsv",
               ["n", "A", "S_1", "S_half", "E_w", "omega"], d3_rows)

    log("")
    log("=== D4: unit-L2 wrap (A=1/√n), S_rel on 1-site, m=0.1 ===")
    log(f"{'n':>4} {'S_1':>12} {'S_half':>12}")
    d4_rows = []
    for n in NS:
        A = amp_unit_L2(n)
        s1 = srel_door(n, 0.1, 1, A, False)
        sh = srel_door(n, 0.1, max(1, n // 2), A, False)
        log(f"{n:4d} {s1:12.5g} {sh:12.5g}")
        d4_rows.append([n, fmt(s1), fmt(sh)])
    dump_table(OUT / "d4_unit_L2.tsv", ["n", "S_1", "S_half"], d4_rows)

    log("")
    log("=== D5: mass robustness, S_1 (A=1, static) ===")
    log(f"{'n':>4} {'m=0':>12} {'m=0.1':>12} {'m=1':>12}")
    d5_rows = []
    for n in NS:
        row = [n]
        logn = f"{n:4d}"
        for m in MS:
            s = srel_door(n, m, 1, 1.0, False)
            row.append(fmt(s))
            logn += f" {s:12.5g}"
        log(logn)
        d5_rows.append(row)
    dump_table(OUT / "d5_mass.tsv", ["n", "m0", "m01", "m1"], d5_rows)

    log("")
    log("=== D6: traveling vs static, half-cycle A=1 m=0.1 (even n) ===")
    log(f"{'n':>4} {'static':>12} {'travel':>12} {'ratio':>10}")
    d6_rows = []
    for n in NS:
        if n % 2:
            continue
        ell = n // 2
        ss = srel_door(n, 0.1, ell, 1.0, False)
        st = srel_door(n, 0.1, ell, 1.0, True)
        log(f"{n:4d} {ss:12.5g} {st:12.5g} {st/ss if ss else float('nan'):10.4g}")
        d6_rows.append([n, fmt(ss), fmt(st), fmt(st / ss if ss else float("nan"))])
    dump_table(OUT / "d6_travel.tsv", ["n", "static", "travel", "ratio"], d6_rows)

    # --- conservation selectors (no S=A/4) ---
    log("")
    log("=== S1: integer argmin of candidate totals, n=3..24 (and even) ===")

    def series():
        Ew = [trough_Ew(n) for n in NS]
        om = [omega1(n, 0.0) for n in NS]
        eps = [atom_eps(n, 0.0) for n in NS]
        # S_rel 1-site A=1 m=0.1 (numeric; m=0 IR on 1-site is noisier)
        S1 = [srel_door(n, 0.1, 1, 1.0, False) for n in NS]
        Sh = [srel_door(n, 0.1, max(1, n // 2), 1.0, False) for n in NS]
        return Ew, om, eps, S1, Sh

    Ew, om, eps, S1, Sh = series()

    candidates = {
        "E_w + n  (ρ=1, v98 1B)": [Ew[i] + NS[i] for i in range(len(NS))],
        "E_w + n ε  (space = n atoms at wrap pitch)": [
            Ew[i] + NS[i] * eps[i] for i in range(len(NS))
        ],
        "E_w/ε + n  (pattern in atoms + site count)": [
            Ew[i] / eps[i] + NS[i] for i in range(len(NS))
        ],
        "S_1 + n": [S1[i] + NS[i] for i in range(len(NS))],
        "S_half + n": [Sh[i] + NS[i] for i in range(len(NS))],
        "S_1  (no space leg)": S1,
        "S_half (no space leg)": Sh,
        "E_w  (no space leg)": Ew,
        "n ε + const  (space atoms only)": [NS[i] * eps[i] for i in range(len(NS))],
        "2π E_w + n  (import Dorau 2π, ρ=1)": [
            2 * math.pi * Ew[i] + NS[i] for i in range(len(NS))
        ],
    }

    even_idx = [i for i, n in enumerate(NS) if n % 2 == 0]
    for name, vals in candidates.items():
        best, mv = argmin_int(NS, vals)
        best_e, mv_e = argmin_int([NS[i] for i in even_idx], [vals[i] for i in even_idx])
        log(f"  {name}")
        log(f"    all n:  n*={best}  min={fmt(mv) if mv is not None else 'nan'}")
        log(f"    even:   n*={best_e}  min={fmt(mv_e) if mv_e is not None else 'nan'}")

    log("")
    log("=== S2: n-dependence, selected columns (m=0.1 A=1 static) ===")
    log(f"{'n':>4} {'E_w':>10} {'ω':>10} {'ε':>10} {'nε':>10} {'S_1':>12} {'S_h':>12} {'E_w/ε':>10}")
    s2_rows = []
    for i, n in enumerate(NS):
        log(
            f"{n:4d} {Ew[i]:10.5g} {om[i]:10.5g} {eps[i]:10.5g} "
            f"{n*eps[i]:10.5g} {S1[i]:12.5g} {Sh[i]:12.5g} {Ew[i]/eps[i]:10.5g}"
        )
        s2_rows.append(
            [n, fmt(Ew[i]), fmt(om[i]), fmt(eps[i]), fmt(n * eps[i]),
             fmt(S1[i]), fmt(Sh[i]), fmt(Ew[i] / eps[i])]
        )
    dump_table(
        OUT / "s2_columns.tsv",
        ["n", "E_w", "omega", "eps", "n_eps", "S_1", "S_half", "Ew_over_eps"],
        s2_rows,
    )

    log("")
    log("=== S3: KMS unit-choice table, ω_1(n)=2 sin(π/n)  [named readings, not 1C] ===")
    log(f"{'target':>12} {'meaning':<44} {'n nearest'}")
    targets = [
        (1.0, "ω = grain (√λ or λ both give n=6 at 1)"),
        (1.0 / (2 * math.pi), "ω = Unruh T=κ/2π at κ=1"),
        (2 * math.pi, "ω = 2π (modular 2π as a frequency)"),
        (2.0, "ω = bandwidth edge 2|sin|≤2"),
    ]
    for val, meaning in targets:
        # find n in reals: 2 sin(π/n) = val ⇒ n = π / arcsin(val/2)
        if 0 < val <= 2:
            n_real = math.pi / math.asin(val / 2.0)
        else:
            n_real = float("nan")
        # nearest integer ≥ 3
        if n_real == n_real:
            cand = max(3, int(round(n_real)))
            # also report even
            even = cand if cand % 2 == 0 else cand + 1
            log(f"  {val:12.5g}  {meaning:<44}  n_real={n_real:.4g}  nearest={cand}  even≥={even}")
        else:
            log(f"  {val:12.5g}  {meaning:<44}  (no real n: target outside bandwidth)")

    log("")
    log("=== S4: dE_w/dn + ρ, ρ forced as ε(n) (same-n, not a constant ρ) ===")
    log("  If space cost per added site is ε(n)=ω/2π, stationarity is")
    log("  E_w(n+1)-E_w(n) + ε(n) = 0.  Discrete Δ, not a derivative.")
    log(f"{'n':>4} {'ΔE_w':>12} {'ε(n)':>12} {'ΔE_w+ε':>12} {'ΔE_w+1':>12}")
    s4_rows = []
    for i, n in enumerate(NS[:-1]):
        dEw = trough_Ew(n + 1) - trough_Ew(n)
        e = atom_eps(n, 0.0)
        log(f"{n:4d} {dEw:12.5g} {e:12.5g} {dEw+e:12.5g} {dEw+1:12.5g}")
        s4_rows.append([n, fmt(dEw), fmt(e), fmt(dEw + e), fmt(dEw + 1)])
    dump_table(OUT / "s4_discrete_dE.tsv",
               ["n", "dEw", "eps", "dEw_plus_eps", "dEw_plus_1"], s4_rows)

    log("")
    log("DONE")
    (OUT / "srel.out").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
