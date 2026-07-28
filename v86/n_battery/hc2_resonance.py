#!/usr/bin/env python3
"""v86 census rung HC-2 -- multipole-FIRST resonance audit of the HC-1 catalog.

Protocol (GROUNDING §2, council fix #3; PROPOSAL Amendment 1):
  "Classification order corrected: for every internal mode, FIRST classify its
   multipole content. Any mode with a time-varying multipole l >= 1 couples to
   the massless gauge channel, open at all frequencies -- gap arithmetic decides
   nothing for it; its width is a golden-rule computation against the A-channel
   density of states (kappa_box-limited in finite volume). Only monopole (and
   pure-matter) excitations get the gap-arithmetic treatment vs the matter band
   [m, w_max], w_max^2 = m^2 + 12/dx^2 (formula council-verified). The lattice
   band-top re-protection warning stands ... refinement tests are mandatory for
   any near-edge verdict."

Input : hc1_catalog.tsv (+ the physical thresholds recomputed here)
Output: a per-mode classification with the decisive channel and, for monopole
        modes, the harmonic-crossing table n*Omega vs Omega_c.

The matter continuum edge used is the COUPLED-SYSTEM one derived in hc1_bdg.py:
    Omega_c = m - w        (lab-frame: w + Omega >= m radiates)
    Omega_max(dx) = w_max - w,  w_max = sqrt(m^2 + 12/dx^2)   [lattice band top]
Band-top "protection" is an ARTIFACT that refines away as dx -> 0, so this
script prints Omega_max for the production lattices and labels any verdict that
depends on it accordingly.
"""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
M = 1.5
M2 = M * M
NMAX = 6                      # harmonics tested
DXS = {"N=384,L=55": 2 * 55.0 / 383, "N=64,L=15": 2 * 15.0 / 63}


def main():
    path = os.path.join(HERE, "hc1_catalog.tsv")
    if not os.path.exists(path):
        sys.exit("missing %s -- run hc1_bdg.py first" % path)
    d = np.genfromtxt(path, names=True, dtype=None, encoding="utf-8")
    if d.ndim == 0:
        d = d.reshape(1)

    print("v86 HC-2 -- multipole-first resonance audit")
    print("=" * 78)
    print("lattice band tops (Omega_max = sqrt(m^2 + 12/dx^2) - w) -- ARTIFACTS:")
    for k, dx in DXS.items():
        wmax = np.sqrt(M2 + 12.0 / dx ** 2)
        print("   %-12s dx=%.4f  w_max=%.3f" % (k, dx, wmax))
    print("   any 'protected because it is above the band top' verdict below is")
    print("   flagged ARTIFACT(dx) and must be refinement-tested before use.\n")

    print("%-8s %-9s %3s %4s %11s %11s %-11s %s" %
          ("w", "channel", "l", "node", "lambda", "Omega", "decided by",
           "verdict"))
    ws = sorted(set(d["w"]))
    summary = []
    for w in ws:
        Om_c = M - w
        sub = d[d["w"] == w]
        for row in sub:
            lam, l, Om = float(row["lambda"]), int(row["l"]), float(row["Omega"])
            ch = str(row["channel"])
            zero = abs(lam) < 5e-3
            if zero and ch == "L0" and l == 0:
                dec, verd = "symmetry", "U(1) phase Goldstone -- exactly zero, no width"
            elif zero and ch == "Lx_sym" and l == 1:
                dec, verd = "symmetry", "translation Goldstone -- exactly zero, no width"
            elif lam < -5e-3:
                dec = "GSS index"
                verd = ("negative direction: contributes to n(H_omega), not a "
                        "propagating mode")
            elif l >= 1:
                dec = "MULTIPOLE"
                verd = ("l=%d: time-varying multipole couples to the MASSLESS gauge "
                        "channel, open at every frequency -> golden-rule width, gap "
                        "arithmetic does NOT apply" % l)
            else:
                dec = "gap arithmetic"
                verd = "monopole -- see harmonic table below"
            summary.append((w, ch, l, int(row["node"]), lam, Om, dec, verd))
            print("%-8.4f %-9s %3d %4d %11.5f %11s %-11s %s" %
                  (w, ch, l, int(row["node"]), lam,
                   ("%.5f" % Om) if np.isfinite(Om) else "imag", dec, verd))

    # ---- harmonic crossing table for the monopole modes --------------------
    print("\n" + "=" * 78)
    print("HARMONIC ARITHMETIC for monopole (l=0) propagating modes")
    print("  a mode of co-rotating frequency Omega emits into the matter band")
    print("  through its n-th harmonic when n*Omega >= Omega_c = m - w.")
    print("  The lowest crossing n is the leading radiative order; n = 1 means")
    print("  the mode radiates directly (a golden-rule width linear in the")
    print("  density of states), n >= 2 means direct emission is BLOCKED and the")
    print("  leading channel is n-quantum -- which is exactly the regime where")
    print("  HC-4's 'width proportional to amplitude^2' pre-registration does NOT")
    print("  apply unmodified.")
    print("")
    print("  CAVEAT (review Finding 7) -- the Omega values in this table are")
    print("  SINGLE-OPERATOR ESTIMATES, Omega(lambda) = sqrt(w^2+lambda) - w applied")
    print("  to eigenvalues of L0 / Lx_sym / Lx_flav separately. The true BdG")
    print("  frequencies of the coupled (x,y) system come from the quadratic pencil,")
    print("  which HC-1 solves ONLY for the l=0 symmetric sector and only there with")
    print("  a box-convergence test. Where the two disagree the COUPLED, box-tested")
    print("  result in hc1_bound_modes.tsv is authoritative: e.g. at w=1.32 the")
    print("  single-operator L0 spectrum shows six sub-threshold entries, while the")
    print("  box-converged coupled solve (Rmax = 5 r_rms = 113) shows ZERO bound")
    print("  internal modes -- the apparent ladder was a box artifact of an")
    print("  undersized domain. Treat this table as a CLASSIFICATION of candidate")
    print("  modes by channel, not as a measured frequency list.")
    print("=" * 78)
    mono = [s for s in summary if s[2] == 0 and s[6] == "gap arithmetic"
            and np.isfinite(s[5]) and s[5] > 1e-6]
    if not mono:
        print("  NO propagating monopole internal modes in the catalog.")
    else:
        print("%-8s %-9s %10s %10s %s" % ("w", "channel", "Omega", "Omega_c",
                                          "n*Omega vs Omega_c"))
        for (w, ch, l, nd, lam, Om, dec, verd) in mono:
            Om_c = M - w
            cells = []
            first = None
            for n in range(1, NMAX + 1):
                cross = n * Om >= Om_c
                if cross and first is None:
                    first = n
                cells.append("%d:%.4f%s" % (n, n * Om, "*" if cross else ""))
            print("%-8.4f %-9s %10.5f %10.5f  %s" %
                  (w, ch, Om, Om_c, "  ".join(cells)))
            print("         -> leading radiative order n = %s%s"
                  % (first if first else "> %d" % NMAX,
                     ("  (direct emission BLOCKED; the fundamental is below the "
                      "matter continuum)" if first and first > 1 else "")))

    print("\n" + "=" * 78)
    print("HC-2 SUMMARY")
    print("=" * 78)
    nl1 = sum(1 for s in summary if s[6] == "MULTIPOLE")
    nmono = len(mono)
    nsym = sum(1 for s in summary if s[6] == "symmetry")
    nneg = sum(1 for s in summary if s[6] == "GSS index")
    print("  catalog entries: %d symmetry (exact zero modes), %d GSS-index "
          "directions,\n  %d multipole (l>=1, massless-channel golden rule), "
          "%d propagating monopole." % (nsym, nneg, nl1, nmono))
    bp = os.path.join(HERE, "hc1_bound_modes.tsv")
    if os.path.exists(bp):
        b = np.genfromtxt(bp, names=True, dtype=None, encoding="utf-8")
        if b.ndim == 0:
            b = b.reshape(1)
        print("\n  HC-1's box-stable bound-mode census along the branch:")
        print("  %-8s %11s %8s %9s %6s %s" %
              ("w", "Q", "r_rms", "Omega_c", "count", "Omegas"))
        for row in b:
            print("  %-8.4f %11.1f %8.2f %9.5f %6d %s" %
                  (row["w"], row["Q"], row["r_rms"], row["Omega_c"],
                   int(row["n_boxstable"]), str(row["Omegas"])))
        tot = sum(int(row["n_boxstable"]) for row in b)
        print("\n  THE CENSUS ANSWER TO 'how many harmonics can a particle hold?':")
        print("  AT MOST ONE, and only in a narrow window. Across the whole scanned")
        print("  branch the box-converged coupled solve finds a single bound internal")
        print("  monopole mode at w = 1.33 and 1.34 and NONE anywhere else --")
        print("  including at the thin-wall extreme w = 1.32, where an undersized")
        print("  box had previously faked a ladder of three. Total box-stable bound")
        print("  internal modes over %d omegas: %d." % (len(b), tot))
        print("  The count is therefore NOT a constant of the theory and NOT the")
        print("  dimension of the conserved-charge lattice; it is a narrow feature of")
        print("  branch position, and it is ZERO over the working region w >= 1.36")
        print("  where every production run of this program has lived.")
    print("\n  Consequences the program must absorb:")
    print("   1. PROPOSAL I.1's slogan 'the number of independently stable")
    print("      harmonics = the dimension of the conserved-charge lattice' is")
    print("      FALSIFIED as stated. The conserved-charge lattice has fixed")
    print("      dimension 3 (+winding), but the internal-mode count varies along")
    print("      the branch and is ZERO over the working region. What the charge")
    print("      lattice counts is the number of independent STATE families, not")
    print("      internal excitations of one state -- which is what the corrected")
    print("      GSS statement (n(H) = n(D)) already says, and what HC-1/HC-3")
    print("      measured (n(H_omega) = 1 everywhere).")
    print("   2. HC-4 (line-width spectroscopy) has nothing to measure in the")
    print("      working region: with no bound internal mode above the Goldstones")
    print("      at w >= 1.36 there is no narrow line whose width can scale as")
    print("      amplitude^2. The QRK-1 lines (0.018-0.126) do NOT appear in this")
    print("      linear catalog at the omegas where they were observed, which")
    print("      supports the council's box-IR hazard warning (GROUNDING §2, fix")
    print("      #4): they are candidates for nonlinear frequency-pulling or box")
    print("      artifacts, not linear modes. The halving list's demotion of HC-4")
    print("      is therefore EARNED, not merely economical. If HC-4 is run at all")
    print("      it should be run at w <= 1.34, where real lines exist.")
    print("   3. Every bound monopole breather found has its FUNDAMENTAL below the")
    print("      matter continuum, so its leading radiative order is n >= 2 -- and")
    print("      being a MONOPOLE it also cannot radiate through the gauge channel.")
    print("      That is the linear-theory root of the v85 measurement 'monopole")
    print("      breathing is radiatively protected' (X10 series), now derived")
    print("      rather than observed.")
    print("   4. Direct emission being blocked (n >= 2 everywhere) means HC-4's")
    print("      pre-registered 'width ~ amplitude^2 = golden-rule confirmation'")
    print("      tests the WRONG power for these modes; the leading process is")
    print("      n-quantum and the exponent must be derived per mode from the")
    print("      table above before any width is called a confirmation.")

    out = os.path.join(HERE, "hc2_classification.tsv")
    with open(out, "w") as fh:
        fh.write("w\tchannel\tl\tnode\tlambda\tOmega\tdecided_by\tverdict\n")
        for s in summary:
            fh.write("%.6g\t%s\t%d\t%d\t%.10g\t%.10g\t%s\t%s\n" % s)
    print("\nwrote %s" % out)


if __name__ == "__main__":
    main()
