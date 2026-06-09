#!/usr/bin/env python3
"""
v60 GEN8 -- end-to-end regression harness for the whole dynamical-Lagrangian loop.

Re-runs every generation's verification (Python + C + Maxima) and reports a single
pass/fail table.  This confirms the GEN1-7 arc is internally consistent and
reproducible as a whole (the Lean modules are built separately via `lake env lean`).

Run:  python3 17_verify_all.py
"""
import subprocess, sys, os, shutil

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)

PY = [
    ("GEN1 first-order parent",        "10_firstorder_parent.py"),
    ("GEN2 covariant first-order",     "11_covariant_firstorder.py"),
    ("GEN3 matter sector",             "12_matter_sector.py"),
    ("GEN4 coupling + EP",             "13_coupling_ep.py"),
    ("GEN5 spectrum",                  "14_spectrum.py"),
    ("GEN6 selection + rank tension",  "15_selection_rank.py"),
    ("GEN7 dynamics cross-check",      "16_dynamics_check.py"),
]
MAC = [
    ("GEN2 grade projection (Maxima)", "11_grade_projection.mac"),
    ("GEN3 Koide invariants (Maxima)", "12_koide_invariants.mac"),
    ("GEN5 Hessian PSD (Maxima)",      "14_hessian_psd.mac"),
    ("GEN6 Koide universal (Maxima)",  "15_koide_universal.mac"),
]
C = [
    ("GEN4 Newton/EP (C)",  "13_newton_ep.c"),
    ("GEN7 dynamics (C)",   "16_dynamics.c"),
]

results = []

def run(label, cmd, ok_check):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        ok = ok_check(r)
    except Exception as e:
        ok = False
        r = None
    results.append((label, ok))
    print(f"  [{'OK ' if ok else 'FAIL'}] {label}")
    return ok

print("=" * 72)
print("v60 GEN8 -- full regression of the dynamical-Lagrangian loop (GEN1-7)")
print("=" * 72)

print("\nPython verifications:")
for label, f in PY:
    run(label, [sys.executable, f], lambda r: r.returncode == 0 and "ALL CHECKS PASSED" in r.stdout)

print("\nMaxima cross-checks:")
have_maxima = shutil.which("maxima") is not None
for label, f in MAC:
    if not have_maxima:
        results.append((label, None)); print(f"  [skip] {label} (no maxima)"); continue
    # Maxima batch echoes the source (incl. the literal error("...FAILED") string),
    # so grepping "FAILED" is unreliable. The robust failure marker is Maxima's own
    # error banner "To debug this try"; on success a "PASSED" output line appears.
    run(label, ["maxima", "--very-quiet", "-b", f],
        lambda r: "To debug this try" not in (r.stdout + (r.stderr or "")) and "PASSED" in r.stdout)

print("\nC programs (compile + run):")
for label, f in C:
    exe = "/tmp/" + os.path.splitext(f)[0]
    cc = subprocess.run(["gcc", "-O2", "-o", exe, f, "-lm"], capture_output=True, text=True)
    if cc.returncode != 0:
        results.append((label, False)); print(f"  [FAIL] {label} (compile)"); continue
    run(label, [exe], lambda r: r.returncode == 0 and "PASSED" in r.stdout)

print("\n" + "=" * 72)
passed = sum(1 for _, ok in results if ok is True)
failed = sum(1 for _, ok in results if ok is False)
skipped = sum(1 for _, ok in results if ok is None)
print(f"REGRESSION SUMMARY: {passed} passed, {failed} failed, {skipped} skipped")
print("=" * 72)
if failed:
    print("REGRESSION FAILED")
    sys.exit(1)
print("REGRESSION PASSED -- the full GEN1-7 arc is reproducible and consistent.")
