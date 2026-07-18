#!/usr/bin/env python3
"""
Offline Round-1 number generator for ND (pure stdlib).
Recomputes closed-form J5/form-factor and discrete free-wave dispersion + leapfrog.
Writes outputs/ — same as sandboxes when run online.
"""
from __future__ import annotations

import json
import math
import os
import time

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "outputs")
os.makedirs(OUT, exist_ok=True)

# Import sibling sandboxes when available
import importlib.util


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    t0 = time.time()
    j5 = load("j5", os.path.join(ROOT, "sandbox_j5_formfactor.py"))
    fw = load("fw", os.path.join(ROOT, "sandbox_free_wave_c.py"))

    # Run J5 analytic path + small grid
    import sys

    sys.argv = ["sandbox_j5_formfactor.py", "--N", "16", "--iters", "120"]
    j5.main()

    sys.argv = ["sandbox_free_wave_c.py", "--c", "1.0", "--Nx", "301"]
    fw.main()

    # Extra: discrete dispersion relation analytic check
    c = 1.0
    L = 40.0
    Nx = 401
    dx = L / (Nx - 1)
    k = 0.8
    # ω_num = (2c/dx) sin(k dx / 2) for standard second-order FD wave
    omega_num = (2.0 * c / dx) * math.sin(k * dx / 2.0)
    v_phase_num = omega_num / k
    rel = abs(v_phase_num / c - 1.0)
    disp = {
        "formula": "omega = (2c/dx) sin(k dx/2); v_phase=omega/k",
        "c": c,
        "k": k,
        "dx": dx,
        "omega_continuum": c * k,
        "omega_discrete_fd2": omega_num,
        "v_phase_discrete": v_phase_num,
        "rel_err": rel,
        "pass_5pct": rel < 0.05,
        "pass_0.1pct": rel < 0.001,
        "note": "Exact discrete dispersion for leapfrog/FD2 free wave; continuum limit dx→0 recovers v=c",
    }
    with open(os.path.join(OUT, "free_wave_dispersion.json"), "w") as f:
        json.dump(disp, f, indent=2)

    print("\n=== DISCRETE DISPERSION (analytic FD2) ===")
    print(json.dumps(disp, indent=2))
    print(f"\noffline total wall {time.time()-t0:.2f}s")


if __name__ == "__main__":
    main()
