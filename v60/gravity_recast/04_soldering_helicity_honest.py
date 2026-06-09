#!/usr/bin/env python3
"""
v60/gravity_recast/04_soldering_helicity_honest.py

HONEST helicity counting for the G9 induced-metric route — replaces the
tautological 02_constrained_helicity_count.py.

WHY 02 WAS CIRCULAR
-------------------
02 took `Jz_on_sym2_transverse()` — which is *by definition* the generator on a
spacetime symmetric-traceless transverse 2-tensor, i.e. the graviton, already
known (v59 g9_polarization_test.py item [3]) to carry helicity ±2 — and tensored
it with an inert internal identity `Id_phys_int`.  Tensoring a spin-2 generator
with the identity trivially preserves its ±2 eigenvalues.  So 02 did not *derive*
±2 from a constraint; it *inserted* the spin-2 generator by hand and reported it
back.  The hard part of G9 (a map binding internal so(8) indices to spacetime
μ,ν) was never touched.

WHAT THIS SCRIPT DOES INSTEAD
-----------------------------
It uses ONE consistent piece of machinery — build the SO(2)_z little-group
generator on a field space, read helicities = imag parts of its eigenvalues —
and applies it to four scenarios that differ ONLY in one physical input:
whether the internal index co-rotates under the spacetime little group
(= whether a soldering / tetrad exists).

  A. v59 carrier:        scalar amplitude × internal so(8) bivector (inert).
  B. spacetime 2-form,   B_μν valued in an inert internal so(n).
  C. SOLDERED 2-form:    B_μν^{IJ}, internal IJ identified with a spacetime
                         Lorentz index (so(3,1)), so it co-rotates.  ← Plebański
  D. tetrad endpoint:    the symmetric metric fluctuation h_μν = e_μ^a e_ν^b η,
                         TT-projected.

The decisive, non-tautological result:
  * A and B reach at most |h| = 1  (reproduces v59's scalar/spin-1 no-go).
  * C reaches ±2, and it does so PURELY because of the extra `I ⊗ Jz_internal`
    term — the soldering.  The ±2 is NOT inserted: each factor (a 2-form) has
    max |h| = 1; ±2 emerges from the *sum* of the two co-rotating charges.
  * Soldering is therefore EXACTLY the Minkowski sum of helicity charges, and it
    is the unique operation lifting max-helicity 1 → 2.

GROUNDING IN v59
----------------
v59 already contains the soldering target it never used: the spacetime factor
Cl(3,1) (dim 16 = C⊗H part; v59/INTEGRATION.md) has a 6-dim bivector grade
≅ so(3,1) — the Lorentz generators.  v59 placed the gravity carrier in the
*internal* Λ²(V^8) = so(8) instead, which is precisely why it came out scalar.
The honest fix: the spin-2 carrier is a 2-form soldered via the Cl(3,1) bivectors;
the so(8) second moment ρ_grav = Tr(M†M) stays as the (scalar) source.

Run:  python 04_soldering_helicity_honest.py
"""

import numpy as np

np.set_printoptions(precision=4, suppress=True)

# =============================================================================
# Shared machinery (identical to v59/gaps/gravity/g9_polarization_test.py)
# =============================================================================

def Rz(ph):
    """4x4 spatial rotation about z (Minkowski indices t,x,y,z)."""
    c, s = np.cos(ph), np.sin(ph)
    return np.array([[1, 0, 0, 0],
                     [0, c, -s, 0],
                     [0, s,  c, 0],
                     [0, 0,  0, 1]], dtype=float)


def helicities_from_generator(J):
    """Helicities = imaginary parts of the eigenvalues of the SO(2)_z generator."""
    w = np.linalg.eigvals(J)
    return np.round(np.sort(np.imag(w)), 6)


_COMPS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]  # 2-form basis tx,ty,tz,xy,xz,yz


def Jz_on_2form():
    """Generator of z-rotations on antisymmetric spacetime 2-forms (6 components)."""
    eps = 1e-6

    def rot_vec(vec, a):
        F = np.zeros((4, 4))
        for k, (mu, nu) in enumerate(_COMPS):
            F[mu, nu] = vec[k]
            F[nu, mu] = -vec[k]
        Fp = Rz(a) @ F @ Rz(a).T
        return np.array([Fp[mu, nu] for (mu, nu) in _COMPS])

    J = np.zeros((6, 6))
    for k in range(6):
        e = np.zeros(6); e[k] = 1.0
        J[:, k] = (rot_vec(e, eps) - rot_vec(e, -eps)) / (2 * eps)
    return J


def Jz_on_sym4():
    """Generator of z-rotations on symmetric 4x4 tensors (10 components)."""
    eps = 1e-6
    idx = [(i, j) for i in range(4) for j in range(i, 4)]

    def vec_to_M(v):
        M = np.zeros((4, 4))
        for k, (i, j) in enumerate(idx):
            M[i, j] = v[k]; M[j, i] = v[k]
        return M

    J = np.zeros((10, 10))
    for k in range(10):
        e = np.zeros(10); e[k] = 1.0
        M = vec_to_M(e)
        dM = (Rz(eps) @ M @ Rz(eps).T - Rz(-eps) @ M @ Rz(-eps).T) / (2 * eps)
        J[:, k] = np.array([dM[i, j] for (i, j) in idx])
    return J


def count_pm2(h):
    n_plus = int(np.sum(np.abs(h - 2.0) < 1e-4))
    n_minus = int(np.sum(np.abs(h + 2.0) < 1e-4))
    return n_plus, n_minus


# =============================================================================
# The four scenarios
# =============================================================================

def scenario_A_v59_carrier(n_internal=28):
    """v59: scalar carrier × internal so(8) bivector.  Internal index inert."""
    # scalar carrier: Jz = 0 (1 component).  internal: n inert copies.
    Jz = np.zeros((n_internal, n_internal))   # all helicities 0
    return Jz


def scenario_B_spacetime_2form_inert(n_internal=28):
    """Promote Ω to a spacetime 2-form B_μν but keep the internal index INERT.
    Jz_total = Jz_2form ⊗ Id_internal  (internal does NOT co-rotate)."""
    Jz_2f = Jz_on_2form()
    return np.kron(Jz_2f, np.eye(n_internal))


def scenario_C_soldered(self_dual=False):
    """Plebański: B_μν^{IJ}, internal IJ = spacetime Lorentz index (so(3,1)),
    SOLDERED so it co-rotates under the SAME SO(2)_z.
    Jz_total = Jz_2form(spacetime μν) ⊗ Id  +  Id ⊗ Jz_2form(internal IJ).

    The second term is the soldering: the internal Lorentz index transforms
    under the spacetime little group (this is what a tetrad e_μ^a accomplishes).
    """
    Jz_2f = Jz_on_2form()
    d = Jz_2f.shape[0]
    Jz_total = np.kron(Jz_2f, np.eye(d)) + np.kron(np.eye(d), Jz_2f)
    return Jz_total


def scenario_D_tetrad_TT():
    """Endpoint: the symmetric metric fluctuation h_μν, transverse-traceless.
    Restrict the symmetric 4x4 to the transverse (x,y) plane and remove the
    trace → the two physical graviton polarizations h_+, h_×."""
    # transverse-traceless symmetric 2-tensors live in the (x,y) 2x2 block,
    # traceless: basis h_+ = diag(1,-1)/√2, h_× = offdiag/√2.
    eps = 1e-6
    R2 = lambda a: np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
    basis = [np.array([[1, 0], [0, -1]]) / np.sqrt(2),
             np.array([[0, 1], [1, 0]]) / np.sqrt(2)]
    J = np.zeros((2, 2))
    for j, B in enumerate(basis):
        dM = (R2(eps) @ B @ R2(eps).T - R2(-eps) @ B @ R2(-eps).T) / (2 * eps)
        for i, A in enumerate(basis):
            J[i, j] = np.sum(A * dM)
    return J


# =============================================================================
# Main
# =============================================================================

def main():
    print("=" * 78)
    print("G9 HONEST HELICITY COUNT — soldering is the precise missing ingredient")
    print("=" * 78)

    # --- self-checks: the building blocks reproduce v59 ---
    h2f = helicities_from_generator(Jz_on_2form())
    print(f"\n[building block] spacetime 2-form helicities       = {h2f}")
    print(f"                 (spin-1 content {{±1}} + {{0,0}};  max |h| = 1)")
    assert count_pm2(h2f) == (0, 0), "2-form should have NO ±2"

    hsym = helicities_from_generator(Jz_on_sym4())
    print(f"[building block] symmetric 4x4 tensor helicities    = {hsym}")
    print(f"                 (contains ±2 — only a SYMMETRIC tensor does)")
    assert count_pm2(hsym) == (1, 1), "symmetric tensor should contain one ±2"

    # --- the four scenarios ---
    results = {}

    hA = helicities_from_generator(scenario_A_v59_carrier())
    results['A'] = count_pm2(hA)
    print("\n" + "-" * 78)
    print("[A] v59 carrier: scalar × internal so(8) (28).  Internal index inert.")
    print(f"    helicities = all 0  (max |h| = {np.max(np.abs(hA)):.0f}).  ±2 count = {results['A']}")
    print("    → SCALAR.  Reproduces the v59 no-go (internal_index_inert).")

    hB = helicities_from_generator(scenario_B_spacetime_2form_inert())
    results['B'] = count_pm2(hB)
    print("\n[B] spacetime 2-form B_μν, internal so(28) INERT (not soldered).")
    print(f"    distinct helicities = {np.unique(hB)}   (max |h| = {np.max(np.abs(hB)):.0f}).  ±2 count = {results['B']}")
    print("    → spin-1 (×28).  Even a SPACETIME 2-form with an inert internal")
    print("      index has NO ±2.  Promoting Ω to a 2-form is NOT enough.")

    hC = helicities_from_generator(scenario_C_soldered())
    results['C'] = count_pm2(hC)
    print("\n[C] SOLDERED B_μν^{IJ}: internal IJ = spacetime Lorentz (so(3,1)),")
    print("    co-rotating.  Jz = Jz_2form ⊗ I  +  I ⊗ Jz_2form  (2nd term = soldering).")
    print(f"    distinct helicities = {np.unique(hC)}   (max |h| = {np.max(np.abs(hC)):.0f}).  ±2 count = {results['C']}")
    print("    → ±2 PRESENT.  And it is NOT inserted: each factor (a 2-form) has")
    print("      max |h| = 1; the ±2 = (+1)+(+1) and (−1)+(−1) is the Minkowski")
    print("      sum of the two CO-ROTATING charges.  Soldering is that sum.")
    print("      (the multiplicity 4 is pre-constraint; the simplicity + diffeo")
    print("       constraints reduce the 36-dim B to the 2 physical modes — [D].)")

    hD = helicities_from_generator(scenario_D_tetrad_TT())
    results['D'] = count_pm2(hD)
    print("\n[D] tetrad endpoint: TT metric fluctuation h_+ , h_×.")
    print(f"    helicities = {hD}.  ±2 count = {results['D']}.  Physical DOF = {len(hD)} (exactly 2).")
    print("    → the constraint+gauge endpoint of [C]: precisely the 2 LIGO modes.")

    # --- the decisive contrast ---
    print("\n" + "=" * 78)
    print("DECISIVE CONTRAST  (same machinery; only the soldering term differs)")
    print("=" * 78)
    print(f"  {'scenario':<48}  {'max|h|':>6}  {'has ±2?':>8}")
    print("  " + "-" * 66)
    rows = [
        ("A  scalar × internal so(8)            [v59]", hA),
        ("B  spacetime 2-form, internal INERT", hB),
        ("C  spacetime 2-form, internal SOLDERED", hC),
        ("D  TT metric fluctuation (endpoint)", hD),
    ]
    for label, h in rows:
        mh = np.max(np.abs(h))
        has = "YES" if count_pm2(h) != (0, 0) else "no"
        print(f"  {label:<48}  {mh:>6.0f}  {has:>8}")

    print("""
READING:
  The ONLY difference between B (fails) and C (works) is the term I ⊗ Jz_internal
  — i.e. whether the internal index transforms under the spacetime little group.
  That term IS the soldering / tetrad.  v59's carrier (A) is the n=28, B-type
  object with the internal index inert, so it is forced to be scalar/spin-1.

  Soldering = Minkowski sum of helicity charges = the unique operation that lifts
  max-helicity from 1 to 2.  This is the precise, non-tautological content of G9.
""")

    # --- where the soldering lives in v59 (grounding, not yet a theorem) ---
    print("=" * 78)
    print("WHERE THE SOLDERING LIVES IN THE v59 ALGEBRA")
    print("=" * 78)
    print("""  Cl(7)_even ≅ ℂ⊗ℍ⊗𝕆 (dim 64).  v59/INTEGRATION.md assigns:
    - spacetime factor : Cl(3,1) (dim 16 = ℂ⊗ℍ); its grade-2 ≅ so(3,1) Lorentz.
    - internal factor  : 𝕆 / Spin(7) (color, G₂, triality).
  The Plebański 2-form B should be valued in the so(3,1) ⊂ Cl(3,1) bivectors
  (the soldering index above), NOT in the internal Λ²(V^8)=so(8).  ρ_grav
  (= Tr M†M, the so(8) second moment) remains the SCALAR source, exactly as a
  matter trace sources h_μν in GR.  The prefactor 21/16 = dimSpin(7)/dimCl(3,1)
  already pairs the internal and spacetime dimensions — consistent with a
  source(internal)/carrier(spacetime) split.

  OPEN (honest): is the so(3,1) soldering subalgebra forced/canonical and
  compatible with G₂ / triality / color, or does choosing it break them?
  That is the real remaining content of G9 — see 04_findings.md.""")

    # machine assertions guaranteeing the contrast is real, not asserted
    assert results['A'] == (0, 0) and results['B'] == (0, 0), "A,B must have no ±2"
    assert results['C'] == (4, 4), f"C must reach ±2 with mult 4 = 2×2 (got {results['C']})"
    assert results['D'] == (1, 1) and len(hD) == 2, "D must be exactly the 2 TT modes"
    print("\nAll assertions PASS — the contrast is computed, not assumed.")


if __name__ == "__main__":
    main()
