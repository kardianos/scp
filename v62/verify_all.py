#!/usr/bin/env python3
"""
v62/verify_all.py  -- regression harness for the v62 number-type program.

Runs every self-verifying script and aggregates the embedded-check results, in
the project convention of v61/06_verify_all.py and v60/lagrangian/17_verify_all.py.

  python3 v62/verify_all.py

Lean (separate, needs the v59 Mathlib lake env):
  cd v59/furey_construction/lean && lake env lean ../../../v62/lean/PhaseNoGo.lean
  cd v59/furey_construction/lean && lake env lean ../../../v62/lean/EwVevHome.lean
"""
import subprocess
import sys
import os

HERE = os.path.dirname(os.path.abspath(__file__))

SCRIPTS = [
    ("no_go/transcendence_nogo.py",            "Brannen-phase transcendentality no-go"),
    ("flat_direction/flat_direction_demo.py",  "flat-direction law (phase + democracy)"),
    ("residual_audit/democracy_selection.py",  "democracy selection (3 criteria)"),
    ("residual_audit/phase_dynamical_home.py", "phase P3 loop-ratio home"),
]

LEAN = [
    ("lean/PhaseNoGo.lean",  "phase no-go (1 cited sorry: Lindemann-Weierstrass)"),
    ("lean/EwVevHome.lean",  "EW democracy / S^783 (no sorry)"),
]


def main():
    print("=" * 78)
    print("v62 verify_all  --  number-type map of the SCP residuals")
    print("=" * 78)
    results = []
    for rel, desc in SCRIPTS:
        path = os.path.join(HERE, rel)
        print(f"\n>>> {rel}  ({desc})")
        r = subprocess.run([sys.executable, path], capture_output=True, text=True)
        passed = (r.returncode == 0) and ("ALL CHECKS PASSED" in r.stdout)
        # surface the embedded-check summary line(s)
        for line in r.stdout.splitlines():
            if line.strip().startswith("[PASS]") or line.strip().startswith("[FAIL]"):
                print("    " + line.strip())
        print(f"    -> {'PASS' if passed else 'FAIL'} (exit {r.returncode})")
        if not passed and r.stderr:
            print("    stderr:", r.stderr.strip().splitlines()[-1] if r.stderr.strip() else "")
        results.append((rel, passed))

    print("\n" + "=" * 78)
    npass = sum(1 for _, p in results if p)
    for rel, p in results:
        print(f"  [{'PASS' if p else 'FAIL'}] {rel}")
    print(f"\n  PYTHON REGRESSION: {npass}/{len(results)} pass")
    print("\n  LEAN (run separately in the v59 lake env):")
    for rel, desc in LEAN:
        print(f"    - {rel}: {desc}")
    print("=" * 78)
    return 0 if npass == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
