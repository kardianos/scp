#!/usr/bin/env python3
"""v86 HC-6 -- emit flavoured radial profiles for the converse-decay test.

HC-6 was previously gated: HC-3 found n(D) = 1 on the one detuning ray it
scanned, so there were no GSS-unstable partitions to seed, and NEXT_PROGRAM.md
recorded the live possibility that the flavoured branch has no unstable region
at all -- in which case HC-6 would be unrunnable.

HC-3-VOLUME SETTLED IT: it does. 20 window-valid partitions carry n(D) = 0
rather than 1, all clustered near the TOP of the omega window, past the VK
turning point of the symmetric branch at omega ~ 1.481 (where dQ/dw changes
sign; measured directly). The original ray sat at wbar = 1.42, far below that
turning point, which is why it saw only n(D) = 1.

GSS requires n(H_omega) = n(D). HC-1-gauged establishes n(H_omega) = 1 on the
branch. So n(D) = 0 partitions FAIL the index match and are decay targets.

This script writes the 4-column "r f0 f1 f2" profiles gen_qball_flavored eats,
for a TARGET (n(D)=0) and a CONTROL (n(D)=1) at comparable total charge. Both
must be run; a target that decays while the control persists is the result. A
target that decays while the control ALSO decays says nothing.

Usage: hc6_profiles.py [--out DIR]
"""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import n4_hc3_flavored as H3            # noqa: E402

DEG = np.pi / 180.0


def delta(theta, rho):
    return rho * np.cos(theta - 120.0 * DEG * np.arange(3))


def solve_partition(wbar, theta, rho):
    F0 = H3.symmetric_seed(wbar)
    w0 = np.array([wbar] * 3)
    F, w = F0.copy(), w0.copy()
    steps = max(2, int(rho / 0.005) + 2)
    F, w = H3.walk_to(F, w, w0 + delta(theta * DEG, rho), nstep=steps)
    if F is None:
        return None, None, None
    return F, w, H3.observables(F, w)


def charge_response(F, w, hw=1e-3):
    D = np.zeros((3, 3))
    for b in range(3):
        Q = []
        for sgn in (+1, -1):
            wp = np.array(w, float)
            wp[b] += sgn * hw
            Fp, wpp = H3.walk_to(F, w, wp, nstep=1)
            if Fp is None:
                return None
            Q.append(H3.observables(Fp, wpp)["Q"])
        D[:, b] = (Q[0] - Q[1]) / (2 * hw)
    ev = np.linalg.eigvalsh(0.5 * (D + D.T))
    tol = 1e-6 * np.max(np.abs(ev))
    return ev, int((ev < -tol).sum())


def write_profile(path, F, w, o, ev, nD, label, rmax=24.0):
    keep = H3.r <= rmax
    r = H3.r[keep]
    with open(path, "w") as fh:
        fh.write("# v86 HC-6 %s partition\n" % label)
        fh.write("# w = (%.10g, %.10g, %.10g)\n" % tuple(w))
        fh.write("# Q = (%.6f, %.6f, %.6f)  Qtot = %.6f  E = %.6f\n"
                 % (o["Q"][0], o["Q"][1], o["Q"][2], o["Qtot"], o["E"]))
        fh.write("# eig(dQ_a/dw_b) = [%.4f %.4f %.4f]   n(D) = %d\n"
                 % (ev[0], ev[1], ev[2], nD))
        fh.write("# GSS: n(H_omega)=1 (HC-1-gauged) vs n(D)=%d -> %s\n"
                 % (nD, "MISMATCH, decay target" if nD != 1 else "match, control"))
        fh.write("# r f0 f1 f2\n")
        for i in range(len(r)):
            fh.write("%.6f %.10e %.10e %.10e\n"
                     % (r[i], F[0][keep][i], F[1][keep][i], F[2][keep][i]))
    return path


def main():
    out = HERE
    if "--out" in sys.argv:
        out = sys.argv[sys.argv.index("--out") + 1]
    os.makedirs(out, exist_ok=True)

    cases = [
        # label,          wbar, theta, rho,   expectation
        ("target_nD0", 1.480, 60.0, 0.030, 0),
        ("target_nD0b", 1.480, 0.0, 0.015, 0),
        ("control_nD1", 1.420, 60.0, 0.030, 1),
        ("control_nD1b", 1.450, 60.0, 0.030, 1),
    ]
    print("v86 HC-6 profile generation")
    print("%-14s %-26s %10s %10s %8s %8s" %
          ("label", "w = (w0,w1,w2)", "Qtot", "E", "n(D)", "expect"))
    made = []
    for label, wbar, th, rho, exp in cases:
        F, w, o = solve_partition(wbar, th, rho)
        if F is None:
            print("%-14s  SOLVE FAILED" % label)
            continue
        cr = charge_response(F, w)
        if cr is None:
            print("%-14s  charge response failed" % label)
            continue
        ev, nD = cr
        p = os.path.join(out, "hc6_%s.txt" % label)
        write_profile(p, F, w, o, ev, nD, label)
        flag = "" if nD == exp else "   *** n(D) NOT AS EXPECTED ***"
        print("%-14s (%.4f,%.4f,%.4f) %10.2f %10.2f %8d %8d%s"
              % (label, w[0], w[1], w[2], o["Qtot"], o["E"], nD, exp, flag))
        made.append((label, p, o, nD))

    print()
    print("wrote %d profiles to %s" % (len(made), out))
    print()
    print("Seeding (ungauged, complex kernel):")
    for label, p, o, nD in made:
        print("  gen_qball_flavored 96 16.0 seeds/hc6_%s.sfa \\" % label)
        print("      %s <w0> <w1> <w2> 0 0 0" % os.path.basename(p))
    print()
    print("PRE-REGISTERED READING (state it before the runs, not after):")
    print("  * target n(D)=0 decays / redistributes charge between components,")
    print("    control n(D)=1 persists  -> GSS converse CONFIRMED")
    print("  * both persist   -> n(D)=0 is not sufficient for decay on the")
    print("    timescale probed; report the bound, do not claim confirmation")
    print("  * both decay     -> the seeding or the box is the instability,")
    print("    not the partition; the test is void and must be redesigned")


if __name__ == "__main__":
    main()
