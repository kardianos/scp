#!/usr/bin/env python3
"""
v61 GEN6 -- end-to-end regression harness for the v61 loop (GEN1-5).

Re-runs every generation's SymPy script and every Maxima cross-check, reporting a
single pass/fail table. (Lean modules are built separately via `lake env lean`.)

Maxima note (learned in GEN5): batch mode ECHOES the source, and `error()` exits 0,
so neither "FAILED" grep nor exit code is reliable. We check the OUTPUT text:
PASS = a line begins with "MAXIMA CHECK PASSED" AND no error banner
("To debug this try") AND no interactive hang ("positive or negative").

Run:  python3 06_verify_all.py
"""
import subprocess, sys, os, shutil

HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)

PY = [
    ("GEN1 curved backreaction",  "01_curved_backreaction.py"),
    ("GEN2 matter backreaction",  "02_backreaction.py"),
    ("GEN3 EW-vev home (R1)",     "03_ew_vev_home.py"),
    ("GEN4 gravitational waves",  "04_gravitational_waves.py"),
    ("GEN5 perihelion precession","05_perihelion.py"),
]
MAC = [
    ("GEN1 Schwarzschild (Maxima ctensor)", "01_schwarzschild.mac"),
    ("GEN2 backreaction (Maxima ctensor)",  "02_backreaction.mac"),
    ("GEN3 Frobenius hat (Maxima)",         "03_frobenius_hat.mac"),
    ("GEN4 quadrupole (Maxima)",            "04_quadrupole.mac"),
    ("GEN5 perihelion (Maxima)",            "05_perihelion.mac"),
]

results = []

def record(label, ok):
    results.append((label, ok)); print(f"  [{'OK ' if ok else 'FAIL'}] {label}")

print("=" * 72)
print("v61 GEN6 -- full regression of the v61 loop (GEN1-5)")
print("=" * 72)

print("\nPython verifications:")
for label, f in PY:
    try:
        r = subprocess.run([sys.executable, f], capture_output=True, text=True, timeout=300)
        record(label, r.returncode == 0 and "ALL CHECKS PASSED" in r.stdout)
    except Exception:
        record(label, False)

print("\nMaxima cross-checks (text-based pass detection):")
have = shutil.which("maxima") is not None
for label, f in MAC:
    if not have:
        results.append((label, None)); print(f"  [skip] {label} (no maxima)"); continue
    try:
        r = subprocess.run(["maxima", "--very-quiet", "-b", f],
                           capture_output=True, text=True, timeout=240)
        out = r.stdout + (r.stderr or "")
        passed_line = any(ln.lstrip().startswith("MAXIMA CHECK PASSED")
                          for ln in out.splitlines())
        no_error = "To debug this try" not in out
        no_hang = "positive or negative" not in out
        record(label, passed_line and no_error and no_hang)
    except Exception:
        record(label, False)

print("\n" + "=" * 72)
passed = sum(1 for _, ok in results if ok is True)
failed = sum(1 for _, ok in results if ok is False)
skipped = sum(1 for _, ok in results if ok is None)
print(f"REGRESSION SUMMARY: {passed} passed, {failed} failed, {skipped} skipped")
print("=" * 72)
if failed:
    print("REGRESSION FAILED"); sys.exit(1)
print("REGRESSION PASSED -- the full v61 GEN1-5 arc is reproducible and consistent.")
