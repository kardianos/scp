#!/usr/bin/env python3
"""v87 B1 — formal construction of an in-kernel Bell test from the fabric's own
phase structure, and its numerical demonstration.

THE IDEA (crank 4). For a stationary fabric object Phi_a = f(r) e^{i omega t},
the profile f(r) is static: the ONLY thing that evolves is the phase. Time IS
phase advance. Two consequences chain:

  (i)  For two IDENTICAL objects (same omega, hence same Q by the branch
       relation Q = f(omega)), the relative phase is FROZEN: both advance at
       omega, so Delta phi(t) = Delta phi(0) for all t. A conserved, transport-
       surviving shared hidden variable.

  (ii) A detector is itself a fabric object (O2). If it carries the SAME omega,
       then in the relative phase the omega*t CANCELS:
              theta_A = (lambda + omega t) - (a + omega t) = lambda - a
       so the measured angle is time-independent. The setting is the detector's
       phase OFFSET a, which is freely dialable at fixed Q. (A setting cannot be
       a phase RATE: changing omega changes Q via the branch, i.e. changes the
       detector into a different object.)

Together these supply all four ingredients a Bell test needs, which the crank-2
review listed as gap #1 "what is a setting in the kernel? undefined".

This file constructs the protocol explicitly and measures what it gives.
"""
import numpy as np

TAU = 2 * np.pi
rng = np.random.default_rng(20260726)
N = 4_000_000


# ---------------------------------------------------------------- the protocol
def outcomes(lam, a, b, Delta0=np.pi):
    """Local, deterministic, dichotomic readouts.

    Object 1 phase: lam + omega t          -> wing A, detector offset a
    Object 2 phase: lam + Delta0 + omega t -> wing B, detector offset b
    Detectors carry the same omega, so omega*t cancels in both relative phases.

    Readout is the A6 fuse/repel dichotomy: co-phase (|rel| < pi/2) -> +1,
    anti-phase -> -1, i.e. sgn(cos(relative phase)).
    """
    A = np.sign(np.cos(lam - a))
    B = np.sign(np.cos(lam + Delta0 - b))
    return A, B


def E_corr(lam, a, b, w=None):
    A, B = outcomes(lam, a, b)
    if w is None:
        return float(np.mean(A * B))
    return float(np.average(A * B, weights=w))


def chsh(lam_fn, coarse=97):
    """max |S| over settings, where lam_fn(b) returns the lambda sample that
    is realised when Bob's setting is b (this is where (MI) lives)."""
    best, arg = 0.0, None
    grid = np.linspace(0, TAU, coarse, endpoint=False)
    for ap in grid:
        for b in grid:
            for bp in grid:
                lb, wb = lam_fn(b)
                lbp, wbp = lam_fn(bp)
                S = (E_corr(lb, 0.0, b, wb) + E_corr(lbp, 0.0, bp, wbp)
                     + E_corr(lb, ap, b, wb) - E_corr(lbp, ap, bp, wbp))
                if abs(S) > best:
                    best, arg = abs(S), (0.0, ap, b, bp)
    return best, arg


# ---------------------------------------------------------------- checks
def check_time_cancellation():
    print("=" * 74)
    print("1.  The omega*t cancellation -- the measured angle is time-independent")
    print("=" * 74)
    omega = 1.4300                      # v86 branch clock
    lam0 = 0.7
    a, b = 0.3, 1.1
    print("   omega = %.4f, lambda(0) = %.3f, a = %.3f, b = %.3f" % (omega, lam0, a, b))
    print("   %10s %14s %14s %8s %8s" % ("t", "obj1 phase", "detA phase", "rel", "A"))
    for t in (0.0, 1.0, 7.3, 55.0, 1000.0):
        ph_obj = (lam0 + omega * t) % TAU
        ph_det = (a + omega * t) % TAU
        rel = (ph_obj - ph_det) % TAU
        A = np.sign(np.cos(rel))
        print("   %10.1f %14.5f %14.5f %8.4f %8.0f" % (t, ph_obj, ph_det, rel, A))
    print("   -> relative phase constant at lambda - a = %.4f for all t.\n"
          % ((lam0 - a) % TAU))
    print("   Same-omega detectors are REQUIRED. With a detuned detector the")
    print("   relative phase drifts at (omega_obj - omega_det) t, and the outcome")
    print("   depends on WHEN the interaction happens -- reintroducing timing as")
    print("   an uncontrolled variable. Q = f(omega) then forces the detector to")
    print("   carry the same CHARGE as the object: the analyser is the same")
    print("   species as the thing analysed, which is O2-consistent.\n")


def check_frozen_relative_phase():
    print("=" * 74)
    print("2.  Frozen relative phase for identical objects")
    print("=" * 74)
    omega = 1.4300
    d0 = 1.234
    ts = np.array([0.0, 10.0, 100.0, 2000.0])
    drift = ((d0 + omega * ts) - (omega * ts)) % TAU - d0
    print("   identical pair (same omega): max |Delta phi(t) - Delta phi(0)| = %.2e"
          % np.max(np.abs(drift)))
    for dw in (1e-4, 1e-3, 1e-2):
        print("   detuned by Delta omega = %.0e : Delta phi drifts %.3f rad over "
              "T=2000  -> %s" % (dw, dw * 2000,
                                 "FROZEN" if dw * 2000 < 0.1 else "SCRAMBLED"))
    print("""
   This is the part that answers G14 in part. GROK's obstruction was that the
   kernel's mixing destroys any correlation prepared in the common past. For
   IDENTICAL monochromatic objects there is nothing to mix: both clocks tick at
   the same rate and the relative phase is a constant of the pair. The fabric
   already relies on this -- A6's co-phase-fuses / anti-phase-repels distinction
   requires a persistent Delta phi, and the v85 pair runs are built on it.
   Tolerance: the pair must be matched to |Delta omega| << 1/T.
""")


def arm1():
    print("=" * 74)
    print("3.  ARM 1 -- (MI) respected. Independent streams for lambda and settings.")
    print("=" * 74)
    def lam_fn(b):
        return rng.uniform(0, TAU, N), None     # lambda independent of b
    S, arg = chsh(lam_fn, coarse=49)
    print("   creation phase uncontrolled, settings from a separate stream")
    print("   max|S| = %.5f      (classical bound 2, triangle rung)" % S)
    ths = np.linspace(0, np.pi, 7)
    lam = rng.uniform(0, TAU, N)
    meas = [E_corr(lam, 0.0, np.pi - t) for t in ths]     # Delta0=pi folds in
    tri = [-1 + 2 * t / np.pi for t in ths]
    print("   E(theta):  " + " ".join("%+.4f" % v for v in meas))
    print("   triangle:  " + " ".join("%+.4f" % v for v in tri))
    print("   max dev from the triangle wave = %.4f"
          % max(abs(m - t) for m, t in zip(meas, tri)))
    print("\n   PREDICTION for the kernel: S = 2.00, not merely S <= 2.\n")
    return S


def arm2():
    print("=" * 74)
    print("4.  ARM 2 -- (MI) violated BY DESIGN. Can the fabric HOST Door A?")
    print("=" * 74)
    def lam_fn(b):
        lam = rng.uniform(0, TAU, N)
        w = np.abs(np.cos(lam + np.pi - b))     # p = 1 tilt on the B-side angle
        return lam, w
    S, arg = chsh(lam_fn, coarse=49)
    print("   creation phase drawn with the p=1 tilt rho(lambda|b) ~ |cos(lambda-b)|")
    print("   max|S| = %.5f      (Tsirelson 2sqrt2 = %.5f)" % (S, 2 * np.sqrt(2)))

    # cost of the tilt, measured
    nb, nl = 2, 64
    bs = [0.0, np.pi / 2]
    joint = np.zeros((nb, nl))
    for i, b in enumerate(bs):
        lam = rng.uniform(0, TAU, N // 4)
        w = np.abs(np.cos(lam + np.pi - b))
        h, _ = np.histogram(lam, bins=nl, range=(0, TAU), weights=w)
        joint[i] = h / h.sum() * 0.5
    pl = joint.sum(axis=0)
    pb = joint.sum(axis=1)
    I = 0.0
    for i in range(nb):
        for j in range(nl):
            if joint[i, j] > 0:
                I += joint[i, j] * np.log2(joint[i, j] / (pb[i] * pl[j]))
    print("   measured cost I(lambda ; b) = %.4f bits" % I)
    print("""
   NOTE the honest comparison: GROK's information-OPTIMAL cost for CHSH at
   Tsirelson is 0.046274 bits. The |cos| tilt used here is the EXACT-singlet
   construction, not the information-optimal one, so it costs more. Both reach
   2sqrt2; they differ in efficiency, not in reachability.

   What arm 2 establishes: the fabric's phase structure can CARRY the Door A
   mechanism. What it does NOT establish: that fabric dynamics GENERATE the
   tilt. The tilt is imposed here by drawing the creation phase from a
   setting-dependent distribution -- exactly the thing that is postulated, not
   derived.
""")
    return S, I


def arm3_readout_fuzziness():
    print("=" * 74)
    print("5.  ARM 3 -- readout fuzziness is a DETECTION LOOPHOLE, not noise")
    print("=" * 74)
    print("""
The A6 fuse/repel dichotomy is sharp only away from |relative phase| = pi/2.
Near the boundary the dynamical outcome is marginal. If marginal events are
DISCARDED rather than assigned, the surviving sub-ensemble is biased toward
well-aligned configurations -- which is precisely crank-2 Reading 2, already
refuted on efficiency grounds. Quantifying the damage:
""")
    lam = rng.uniform(0, TAU, N)
    print("   %10s %10s %12s" % ("dead zone", "efficiency", "max|S| (MI respected)"))
    for dz in (0.0, 0.10, 0.25, 0.50):
        def lam_fn(b, dz=dz):
            keep = np.abs(np.cos(lam - b)) > np.sin(dz)
            return lam[keep], None
        # efficiency
        eff = float(np.mean(np.abs(np.cos(lam)) > np.sin(dz)))
        S, _ = chsh(lam_fn, coarse=25)
        print("   %10.2f %10.4f %12.5f" % (dz, eff, S))
    print("""
   Discarding marginal events inflates S above 2 with NO measurement dependence
   at all -- a pure post-selection artifact. Any in-kernel run MUST assign every
   event (e.g. by the sign of the residual interaction, or a deterministic
   tie-break), and MUST report the fraction of marginal events. An S > 2 obtained
   with discards is the closed detection loophole and must be rejected.
""")


if __name__ == "__main__":
    check_time_cancellation()
    check_frozen_relative_phase()
    s1 = arm1()
    s2, I = arm2()
    arm3_readout_fuzziness()
    print("=" * 74)
    print("SUMMARY")
    print("=" * 74)
    print("   arm 1 (MI respected, fabric's own uniform creation phase) : S = %.4f" % s1)
    print("   arm 2 (MI violated by design, p=1 tilt)                   : S = %.4f "
          "at %.3f bits" % (s2, I))
    print("   arm 3 (marginal events discarded)                         : see above")
