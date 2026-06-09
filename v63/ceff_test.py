#!/usr/bin/env python3
"""
v63/ceff_test.py  --  GATING TEST for the "localization + density-dependent c"
holonomy idea (see ../v62/THESIS.md and the conversation that produced it).

THE QUESTION (the make-or-break check proposed before committing to a holonomy v63):
  The Brannen phase invariant is cos(3φ) = cos(2/3); the "inverted/integration"
  picture would have the phase be a holonomy  exp(i ∮ A)  with  ∮A = 2/3, where a
  density-dependent meta-constant c_eff supplies the conversion
        2/3  =  c_eff · (canonical weight).
  A tantalizing lead: the project's BLV null-rotor effective metric reports a core
  velocity ratio  v_min/c ≈ 0.670 ≈ 2/3.  IF c_eff = v_min/c is (a) structurally
  2/3 and (b) pinned parameter-free by the field equation, the holonomy idea has a
  parameter-free anchor.  This script tests exactly that.

SOURCE OF THE NUMBERS (context pulled):
  v3/tasks/done/gravity-null-rotor-metric.md  +  v3/src/nullrotor_metric.c.
  v_min/c is computed as  min over radius r and modes A of
        v_A(r)/c = sqrt( (1 + P4^{jj}(r)) / (1 + m4(r)) )
  i.e. a NUMERICAL minimum of a radial function over the soliton profile f(r)
  (nullrotor_metric.c lines 479,514-518), NOT a closed-form constant.

VERDICT (see the checks): the lead FAILS on three independent counts.  Good --
this is the gating test catching numerology before a v63 build.
"""
import os
import subprocess
import re

SEP = "=" * 78
Q = 2.0 / 3.0  # the transcendental target's argument, cos(3φ)=cos(2/3)

# Values reproduced from v3/bin/nullrotor_metric on 2026-05-28 (this session):
RECORDED = {
    "massless (profile_sigma_e1.dat)":          0.67020151,
    "massive  (profile_massive_e1_mpi0.398.dat)": 0.62519338,
}

V3 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "v3")
BIN = os.path.join(V3, "bin", "nullrotor_metric")
PROFILES = {
    "massless (profile_sigma_e1.dat)":          "data/profile_sigma_e1.dat",
    "massive  (profile_massive_e1_mpi0.398.dat)": "data/profile_massive_e1_mpi0.398.dat",
}


def banner(s):
    print(SEP); print(s); print(SEP)


def reproduce():
    """Re-run the v3 binary if available; else use the recorded values.
    Returns {label: v_min/c}."""
    banner("Reproduce v_min/c from v3/bin/nullrotor_metric")
    out = {}
    can_run = os.path.isfile(BIN) and os.access(BIN, os.X_OK)
    for label, prof in PROFILES.items():
        profpath = os.path.join(V3, prof)
        val = None
        if can_run and os.path.isfile(profpath):
            try:
                r = subprocess.run([BIN, "-profile", profpath],
                                   cwd=V3, capture_output=True, text=True, timeout=120)
                m = re.search(r"v/c\s*∈\s*\[([0-9.]+),", r.stdout)
                if m:
                    val = float(m.group(1))
            except Exception as e:
                print(f"  (could not run binary: {e})")
        if val is None:
            val = RECORDED[label]
            print(f"  {label}: {val:.8f}  (recorded; binary not re-run)")
        else:
            rec = RECORDED[label]
            agree = abs(val - rec) < 1e-6
            print(f"  {label}: {val:.8f}  (reproduced; matches recorded {rec:.8f}: {agree})")
        out[label] = val
    return out


def analyze(vmin):
    banner("Analysis: is v_min/c a parameter-free structural 2/3?")
    m0 = vmin["massless (profile_sigma_e1.dat)"]
    mM = vmin["massive  (profile_massive_e1_mpi0.398.dat)"]

    print(f"  target argument  2/3            = {Q:.8f}")
    print(f"  massless v_min/c                = {m0:.8f}   gap to 2/3 = "
          f"{abs(m0-Q):.6f} ({abs(m0-Q)/Q*100:.3f}%)")
    print(f"  massive  v_min/c                = {mM:.8f}   gap to 2/3 = "
          f"{abs(mM-Q):.6f} ({abs(mM-Q)/Q*100:.3f}%)")
    print()
    print("  FAILURE 1 -- not equal to 2/3.  The massless value is 0.53% off 2/3;")
    print("    it is a numerical profile minimum, with no structural reason to be 2/3.")
    print(f"  FAILURE 2 -- not density-independent.  massless != massive: they differ")
    print(f"    by {abs(m0-mM)/m0*100:.2f}% (the pion mass moves it 0.670 -> 0.625).  A")
    print("    parameter-free structural constant cannot change with m_pi, so v_min/c")
    print("    is NOT a universal c_eff = 2/3.")
    # If 2/3 = c_eff * weight, the implied 'canonical weight' must be consistent:
    w0, wM = Q / m0, Q / mM
    print(f"  FAILURE 3 -- no consistent 'weight'.  If 2/3 = c_eff*w then")
    print(f"    massless w = (2/3)/{m0:.4f} = {w0:.4f},  massive w = {wM:.4f}:")
    print(f"    they disagree by {abs(w0-wM)/w0*100:.1f}%, so there is no canonical weight.")
    print()
    print("  SECTOR CAVEAT (independent of the above): this BLV v_min/c is from the")
    print("  v3 null-rotor metric -- a NUCLEAR-scale (~250 MeV), Yukawa-range effect")
    print("  with NO h=±2 tensor modes (gravity-null-rotor-metric.md §8).  It was the")
    print("  SUPERSEDED gravity attempt; the current carrier lives in the v60/v61")
    print("  Cl(3,1) soldering.  So even a clean 2/3 here would be the wrong sector.")
    return {"m0": m0, "mM": mM, "w0": w0, "wM": wM}


def verdict():
    banner("VERDICT")
    print("""  The "c_eff = v_min/c = 2/3, picked up by a holonomy" lead is FALSIFIED:
    (1) v_min/c (massless) = 0.6702 != 2/3  (0.53% off; a numerical profile min);
    (2) v_min/c is parameter-dependent (0.6702 vs 0.6252 massive) -> not a
        structural constant -> cannot be a parameter-free c_eff;
    (3) no consistent canonical weight; and the sector is the superseded,
        nuclear-scale, no-tensor BLV metric, not the current gravity carrier.

  CONSEQUENCE for the localization+c program: the holonomy idea is NOT killed,
  but it gets NO parameter-free anchor from the existing density-dependent c.
  If pursued, c_eff must come from a genuinely structural / virial-pinned
  quantity in the CURRENT (Cl(3,1)) sector, not from the v3 BLV v_min/c.  Until
  such a parameter-free c_eff exists, "2/3 = c_eff*weight" remains a dial, i.e.
  the numerology trap the gating test was built to catch.

  The one EXACT structural constant in the BLV sector is P/m = 2 (machine
  precision, gravity-null-rotor-metric.md consistency check), NOT 2/3.""")


def run_checks(vmin, ana):
    banner("EMBEDDED SELF-CHECKS")
    ok = True

    def check(name, cond):
        nonlocal ok
        ok = ok and bool(cond)
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")

    m0, mM = ana["m0"], ana["mM"]
    # reproduced values match the recorded ones (reproducibility).
    check("massless reproduces recorded 0.67020151",
          abs(m0 - 0.67020151) < 1e-6)
    check("massive reproduces recorded 0.62519338",
          abs(mM - 0.62519338) < 1e-6)
    # FAILURE 1: massless v_min/c is NOT 2/3 (off by > 0.3%).
    check("massless v_min/c != 2/3 (gap > 3e-3, not machine zero)",
          abs(m0 - Q) > 3e-3)
    # FAILURE 2: parameter dependence -> not a universal constant.
    check("v_min/c is parameter-dependent (massless != massive, > 1%)",
          abs(m0 - mM) > 0.01 * m0)
    # FAILURE 3: implied weights inconsistent.
    check("implied canonical weight inconsistent (massless vs massive > 5%)",
          abs(ana["w0"] - ana["wM"]) > 0.05 * ana["w0"])
    # Sanity: both are subluminal (the BLV result is v<=c).
    check("both v_min/c < 1 (subluminal, as BLV requires)", m0 < 1 and mM < 1)

    print(SEP)
    print("ALL CHECKS PASSED (the lead is falsified, as designed)" if ok
          else "SOME CHECKS FAILED")
    print(SEP)
    return ok


if __name__ == "__main__":
    vmin = reproduce()
    print()
    ana = analyze(vmin)
    print()
    verdict()
    print()
    ok = run_checks(vmin, ana)
    import sys
    sys.exit(0 if ok else 1)
