#!/usr/bin/env python3
"""
v63/integrate_density.py  --  READING A: "integrate according to local density".

Weight the integration by the local field density rho(x).  Two arenas:

  A1 (generation-space density):  the Brannen masses m_k ARE a density over the
     Z3 generation cycle.  Integrate the generation angle theta_k against that
     density -> the density-weighted CIRCULAR MEAN.  This is a clean, closed-form
     phase-EXTRACTION (synthesis) operation.  We compute it and find it produces
     a transcendental phase observable (depends on cos(3phi)) -- the RIGHT
     number-type -- but it does NOT derive the value 2/9 (phi is input via the
     amplitudes), and at the physical point it is a nontrivial function of phi.

  A2 (spatial soliton density):  weight by the actual soliton profile rho(r)
     (v3 data).  Any dimensionless density-weighted observable inherits the
     soliton's parameter-dependence (the massless and massive profiles differ in
     SHAPE -- already shown by v_min/c = 0.670 vs 0.625).  So this route cannot
     give a universal constant -- same failure mode as the c_eff lead.

Net (Reading A): density-weighting is a valid phase-EXTRACTION and naturally lands
in the transcendental number-type, but it does not DERIVE a parameter-free phase.
"""
import os
import math
import cmath
import numpy as np

SEP = "=" * 78
V3 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "v3")


def banner(s):
    print(SEP); print(s); print(SEP)


# ---------------------------------------------------------------------------
# A1 -- generation-space density-weighted circular mean
# ---------------------------------------------------------------------------
def brannen_masses(t, phi, a=1.0):
    th = [2 * math.pi * k / 3 for k in range(3)]
    s = [a * (1 + 2 * t * math.cos(thk + phi)) for thk in th]
    return th, [sk * sk for sk in s]


def circular_mean_numeric(t, phi, a=1.0):
    th, m = brannen_masses(t, phi, a)
    z = sum(mk * cmath.exp(1j * thk) for mk, thk in zip(m, th))
    return z


def circular_mean_closed(t, phi, a=1.0):
    # Derived: sum_k m_k e^{i theta_k} = a^2 (6 t e^{-i phi} + 3 t^2 e^{2 i phi}).
    return a * a * (6 * t * cmath.exp(-1j * phi) + 3 * t * t * cmath.exp(2j * phi))


def A1_generation_density():
    banner("A1: density-weighted circular mean over the Z3 generation cycle")
    t = math.sqrt(0.5)      # lepton symmetric point t^2 = 1/2
    phi = 2.0 / 9.0         # the Brannen phase (INPUT via the amplitudes)

    z_num = circular_mean_numeric(t, phi)
    z_cf = circular_mean_closed(t, phi)
    print(f"  sum_k m_k e^(i theta_k):  numeric = {z_num:.6f}")
    print(f"                            closed  = {z_cf:.6f}   (match: "
          f"{abs(z_num - z_cf) < 1e-12})")

    Phi = cmath.phase(z_num)
    # closed-form phase: Phi = -phi + atan2(3 t^2 sin 3phi, 6 t + 3 t^2 cos 3phi)
    corr = math.atan2(3 * t * t * math.sin(3 * phi), 6 * t + 3 * t * t * math.cos(3 * phi))
    Phi_cf = -phi + corr
    print(f"\n  density-weighted circular-mean phase  Phi = arg(sum) = {Phi:+.6f} rad")
    print(f"  closed form  -phi + atan2(...) = {Phi_cf:+.6f} rad  (match: "
          f"{abs(Phi - Phi_cf) < 1e-9})")
    print(f"  (input phi = {phi:.6f};  Phi is a NONTRIVIAL function of phi, not -phi)")
    print(f"""
  The correction term carries cos(3 phi) = cos(2/3) = {math.cos(2/3):.6f} -- the
  TRANSCENDENTAL invariant -- so the density-weighted phase observable lives in
  the right number-type.  But phi entered via the amplitudes m_k(phi): this is
  a valid phase-EXTRACTION (the 'integrate to get the form/phase' step), NOT a
  derivation of the VALUE 2/9.
""")
    return {"Phi": Phi, "match_closed": abs(Phi - Phi_cf) < 1e-9}


# ---------------------------------------------------------------------------
# A2 -- spatial soliton density: dimensionless observables are parameter-dependent
# ---------------------------------------------------------------------------
def load_profile(path):
    r, rho = [], []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            c = line.split()
            r.append(float(c[0])); rho.append(float(c[3]))  # col 3 = baryon_density
    return np.array(r), np.array(rho)


def density_moment_ratio(r, rho):
    """Dimensionless, scale-free shape observable: <r^2>_rho / <r>_rho^2,
       with the 3D radial measure r^2 dr.  (Independent of the overall size.)"""
    w = rho * r * r
    w = np.maximum(w, 0)
    m1 = np.trapz(w * r, r)
    m2 = np.trapz(w * r * r, r)
    m0 = np.trapz(w, r)
    mean_r = m1 / m0
    mean_r2 = m2 / m0
    return mean_r2 / (mean_r * mean_r)


def A2_spatial_density():
    banner("A2: density-weighted observables on the actual soliton are parameter-dependent")
    profs = {
        "massless": os.path.join(V3, "data", "profile_sigma_e1.dat"),
        "massive":  os.path.join(V3, "data", "profile_massive_e1_mpi0.398.dat"),
    }
    vals = {}
    for label, p in profs.items():
        if os.path.isfile(p):
            r, rho = load_profile(p)
            vals[label] = density_moment_ratio(r, rho)
            print(f"  {label:9s}: <r^2>/<r>^2 (density-weighted) = {vals[label]:.6f}")
        else:
            print(f"  {label}: profile missing ({p})")
    if len(vals) == 2:
        d = abs(vals["massless"] - vals["massive"])
        print(f"\n  massless vs massive differ by {d:.4f} "
              f"({d/vals['massless']*100:.2f}%) -> parameter-DEPENDENT.")
    print(f"""
  This dimensionless, scale-free density-weighted observable differs between the
  two solitons (just like v_min/c = 0.670 vs 0.625 in ceff_test.py).  The profile
  SHAPE moves with the model parameters (e.g. m_pi), so any spatial
  density-weighted quantity inherits that dependence -- it cannot be a universal
  constant.  Same failure mode as the c_eff lead.
""")
    return vals


# ---------------------------------------------------------------------------
def run_checks(a1, vals):
    banner("EMBEDDED SELF-CHECKS (Reading A)")
    ok = True

    def check(name, cond):
        nonlocal ok
        ok = ok and bool(cond)
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")

    # A1 closed form matches numeric.
    t, phi = math.sqrt(0.5), 2 / 9
    zn, zc = circular_mean_numeric(t, phi), circular_mean_closed(t, phi)
    check("A1 closed form matches numeric circular mean", abs(zn - zc) < 1e-12)
    check("A1 phase matches -phi + atan2 closed form", a1["match_closed"])
    # A1 phase is NOT simply -phi (nontrivial cos(3phi) correction at physical point).
    check("A1 circular-mean phase != -phi (nontrivial transcendental correction)",
          abs(cmath.phase(zn) - (-phi)) > 1e-2)
    # A2 parameter dependence.
    if len(vals) == 2:
        check("A2 spatial density observable is parameter-dependent (massless != massive)",
              abs(vals["massless"] - vals["massive"]) > 1e-3)

    print(SEP)
    print("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
    print(SEP)
    return ok


if __name__ == "__main__":
    a1 = A1_generation_density()
    print()
    vals = A2_spatial_density()
    print()
    ok = run_checks(a1, vals)
    import sys
    sys.exit(0 if ok else 1)
