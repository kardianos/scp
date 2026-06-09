#!/usr/bin/env python3
"""
v65 — Quantitative grounding for the "magic = feature-learning onset" claim.

The neural-net mapping (NN_PHYSICS_EXPLORATION.md §2.2c) identifies the lazy/kernel(NTK)
regime vs the feature-learning regime with the LINEAR vs NONLINEAR sector of the SCP
theory. This script makes that quantitative using the ACTUAL potential the kernel carries,
  V(P) = (mu/2) P^2 / (1 + kappa P^2),   V'(P) = mu P / (1 + kappa P^2)^2,
and pins WHERE the lazy->feature crossover sits in kappa, so the v64-Thread-C / v65-E3 run
knows what to look for (a threshold near kappa ~ 1/P^2, NOT a linear-from-zero response).

Order parameter (relative nonlinearity of the back-reaction force):
  chi(kappa,P) = 1 - V'(P;kappa) / V'(P;0) = 1 - 1/(1 + kappa P^2)^2
  - lazy   (kappa P^2 << 1):  chi ~ 2 kappa P^2      (LINEAR in kappa => "lazy/NTK")
  - feature(kappa P^2 >> 1):  chi -> 1               (SATURATED => "feature learning")
  - crossover at kappa P^2 ~ O(1)  =>  kappa_c ~ 1/P^2.

Falsifiable prediction for E3: scan kappa; the background-reorganization order parameter
must show this crossover near kappa_c ~ 1/P^2. If it is linear in kappa all the way from 0
(no knee), "magic" is just a coefficient, not a regime -> refutes the regime framing.

Embedded self-checks (exit 1 on failure). Analytic, no simulation.
Refs: scp_sim.c:408-410 (den=1+kappa*P^2), v64/magnitude_precheck.py, v64/CONCEPT.md.
"""
import math

KAPPA_STD = 50.0          # CLAUDE.md standard
P_CORE = 1.0              # representative core density (bracketed below)

def chi(kappa, P):
    """Relative nonlinearity of the back-reaction force = feature-learning order param."""
    x = kappa * P * P
    return 1.0 - 1.0 / (1.0 + x)**2

def chi_lazy(kappa, P):
    """Linear (lazy/NTK) leading-order prediction: chi ~ 2 kappa P^2."""
    return 2.0 * kappa * P * P

def crossover_kappa(P, level):
    """kappa where chi(kappa,P) = level. Solve 1/(1+x)^2 = 1-level => x = (1-level)^-0.5 - 1."""
    x = (1.0 - level)**-0.5 - 1.0
    return x / (P * P), x

def report():
    print("=" * 70)
    print("v65 E3 grounding — lazy(NTK) -> feature-learning crossover in kappa")
    print("=" * 70)
    print("Order parameter chi(kappa,P) = 1 - 1/(1+kappa P^2)^2  (relative nonlinearity")
    print("of the back-reaction force V'(P); kernel den = 1+kappa P^2).\n")

    # Where is the crossover, for a bracket of core densities?
    print("Crossover kappa_c (chi = 0.5, i.e. force half-quenched by nonlinearity):")
    print(f"  {'P_core':>7} {'kappa_c':>10} {'kappa P^2 at c':>14}   regime of standard kappa=50")
    for P in [0.3, 0.5, 1.0, 2.0, 3.4]:
        kc, xc = crossover_kappa(P, 0.5)
        xstd = KAPPA_STD * P * P
        regime = "FEATURE (x>>1)" if xstd > 3 else ("crossover" if xstd > 0.3 else "lazy")
        print(f"  {P:>7.2f} {kc:>10.3f} {xc:>14.3f}   x_std={xstd:>7.1f} -> {regime}")
    print()

    # chi vs kappa at P=1: show the lazy-linear region, the knee, and saturation
    print(f"chi(kappa) and the lazy-linear prediction at P={P_CORE} (knee near kappa~1):")
    print(f"  {'kappa':>8} {'chi(full)':>11} {'chi_lazy(2kP^2)':>16} {'chi/lazy':>10}")
    for k in [0.0, 0.1, 0.3, 1.0, 3.0, 10.0, 50.0]:
        c = chi(k, P_CORE); cl = chi_lazy(k, P_CORE)
        ratio = c/cl if cl > 0 else float('nan')
        print(f"  {k:>8.2f} {c:>11.4f} {cl:>16.4f} {ratio:>10.3f}")
    print()

    # Susceptibility chi'(kappa) — peaks at kappa=0 then decays; the "response" curve
    kc, xc = crossover_kappa(P_CORE, 0.5)
    print("VERDICT")
    print("-" * 7)
    print(f"  At standard kappa={KAPPA_STD:.0f}, P~1: kappa P^2 = {KAPPA_STD:.0f} >> 1")
    print(f"  => the theory sits DEEP in the feature-learning (saturated nonlinear)")
    print(f"     regime; chi = {chi(KAPPA_STD,1.0):.3f} (force {100*chi(KAPPA_STD,1.0):.0f}% quenched).")
    print(f"  Lazy/NTK regime is kappa << kappa_c = {kc:.2f} (at P=1); the crossover is")
    print(f"     near kappa ~ 1/P^2 ~ O(1).  E3 PREDICTION: a background-reorganization")
    print(f"     order parameter scanned in kappa must show a KNEE here, not a straight")
    print(f"     line from 0. Linear-all-the-way refutes 'magic = a regime'.")
    print("=" * 70)

def self_checks():
    ok = True
    def ck(name, cond):
        nonlocal ok
        if not cond: ok = False
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    print("\nSELF-CHECKS\n" + "-"*11)
    # 1. lazy limit: chi -> 2 kappa P^2 as kappa->0
    k=1e-4
    ck("chi ~ 2 kappa P^2 as kappa->0", abs(chi(k,1.0)/(2*k) - 1.0) < 1e-3)
    # 2. saturation: chi -> 1 as kappa->inf
    ck("chi -> 1 as kappa->inf", abs(chi(1e6,1.0) - 1.0) < 1e-5)
    # 3. monotone increasing in kappa
    ck("chi monotone in kappa", chi(1.0,1.0) > chi(0.1,1.0) and chi(10.0,1.0) > chi(1.0,1.0))
    # 4. crossover scales as 1/P^2
    kc1,_ = crossover_kappa(1.0,0.5); kc2,_ = crossover_kappa(2.0,0.5)
    ck("kappa_c scales as 1/P^2 (factor 4 for 2x P)", abs(kc1/kc2 - 4.0) < 1e-6)
    # 5. standard kappa is in the saturated (feature) regime at P=1
    ck("standard kappa=50,P=1 is feature regime (chi>0.9)", chi(50.0,1.0) > 0.9)
    return ok

if __name__ == "__main__":
    import sys
    report()
    if not self_checks():
        print("\nSELF-CHECKS FAILED"); sys.exit(1)
    print("\nAll self-checks passed.")
