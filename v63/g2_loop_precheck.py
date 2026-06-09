#!/usr/bin/env python3
"""
v63/g2_loop_precheck.py  --  CHEAP PRE-CHECK before building the full
G2-equivariant finite-sum loop (the v59 Step-6, the only route that could make
the phase loop-ratio cos(3phi)=-c3/(4c6) parameter-free; see ../LOOP_RATIO.md).

Two decisive, cheap questions:

  (Q1) REP-THEORY PREMISE.  v59/algebra/v_eff_loop.py asserts the d-quark loop is
       a finite sum over a "14-dim G2-equivariant subspace of Lambda^4(R^7)"
       (N_X = dim G2 = 14).  Is that true?  We build G2 = Der(O) (the stabilizer
       of the associative 3-form phi in so(7)) from the repo's own Fano table and
       decompose Lambda^2 and Lambda^4 via the quadratic Casimir.

  (Q2) CAN N_X EVEN HELP?  v59's loop formula is
       V_eff^X(phi) = -(N_X/64pi^2) Tr[M_X(phi)^2 log ...].  N_X is an OVERALL,
       phi-INDEPENDENT prefactor -> it scales c3 and c6 EQUALLY -> it CANCELS in
       the ratio c3/c6.  So no G2 counting can un-suppress c6.

FINDING (negative): (Q1) the 14 (adjoint) lives in Lambda^2, NOT Lambda^4
(Lambda^4 = 1+7+27) -- the v_eff_loop counting is misattributed; (Q2) N_X cancels
in c3/c6 anyway.  The c6-suppression is a GENERATION-Z3 effect, orthogonal to the
G2/internal structure.  => the finite-G2-sum escape hatch is illusory; the full
build is NOT warranted; the loop-ratio's c6 boost stays a free (2-loop) dial.
"""
import numpy as np
from itertools import combinations
import math

SEP = "=" * 78

# ---- repo's authoritative Fano table (SevenDAlgebra.lean) -------------------
octMultTable = [
    [(1,0),(1,1),(1,2),(1,3),(1,4),(1,5),(1,6),(1,7)],
    [(1,1),(-1,0),(1,3),(-1,2),(1,5),(-1,4),(-1,7),(1,6)],
    [(1,2),(-1,3),(-1,0),(1,1),(1,6),(1,7),(-1,4),(-1,5)],
    [(1,3),(1,2),(-1,1),(-1,0),(1,7),(-1,6),(1,5),(-1,4)],
    [(1,4),(-1,5),(-1,6),(-1,7),(-1,0),(1,1),(1,2),(1,3)],
    [(1,5),(1,4),(-1,7),(1,6),(-1,1),(-1,0),(-1,3),(1,2)],
    [(1,6),(1,7),(1,4),(-1,5),(-1,2),(1,3),(-1,0),(-1,1)],
    [(1,7),(-1,6),(1,5),(1,4),(-1,3),(-1,2),(1,1),(-1,0)],
]


def banner(s):
    print(SEP); print(s); print(SEP)


def assoc_3form():
    """phi_{abc} (0-indexed over e1..e7) from the Fano table: e_a e_b = +e_c => phi=+1."""
    phi = np.zeros((7, 7, 7))
    for a in range(1, 8):
        for b in range(1, 8):
            sgn, tgt = octMultTable[a][b]
            if a != b and tgt >= 1:
                phi[a - 1, b - 1, tgt - 1] = sgn
    return phi


def so7_basis():
    """21 antisymmetric generators E_{pq} (p<q): entry (p,q)=+1,(q,p)=-1."""
    B = []
    for p in range(7):
        for q in range(p + 1, 7):
            A = np.zeros((7, 7)); A[p, q] = 1; A[q, p] = -1
            B.append(A)
    return B


def act_on_3form(A, phi):
    """(A.phi)_{ijk} = -sum_m (A_{mi}phi_{mjk}+A_{mj}phi_{imk}+A_{mk}phi_{ijm})."""
    t1 = np.einsum('mi,mjk->ijk', A, phi)
    t2 = np.einsum('mj,imk->ijk', A, phi)
    t3 = np.einsum('mk,ijm->ijk', A, phi)
    return -(t1 + t2 + t3)


def build_g2(phi):
    """g2 = {A in so(7): A.phi = 0} = ker of the (21 -> 35) stabilizer map."""
    basis = so7_basis()
    rows = [act_on_3form(A, phi).reshape(-1) for A in basis]
    M = np.array(rows)                      # 21 x 343
    # null space of M^T acting on coefficient vectors: solve c.M = 0  -> ker of M^T?
    # We want combinations c_a A_a with sum_a c_a (A_a.phi)=0  => c in left-null of M.
    U, s, Vt = np.linalg.svd(M, full_matrices=True)
    tol = 1e-9
    rank = int(np.sum(s > tol))
    # left-null space of M (dim 21-rank): rows of U beyond rank
    null = U[:, rank:]                      # 21 x (21-rank)
    gens = []
    for col in range(null.shape[1]):
        c = null[:, col]
        A = sum(c[a] * basis[a] for a in range(21))
        gens.append(A)
    return gens, 21 - rank


def gram_orthonormalize(gens):
    """Orthonormalize w.r.t. <A,B> = -Tr(AB) (= Tr(A^T B), positive-definite on so)."""
    mats = [g.copy() for g in gens]
    out = []
    for A in mats:
        for B in out:
            A = A - (-np.trace(A @ B)) * B
        nrm = math.sqrt(-np.trace(A @ A))
        out.append(A / nrm)
    return out


def wedge_action(A, p):
    """Induced action of A in so(7) on Lambda^p(R^7): returns (C(7,p) x C(7,p)) matrix."""
    subs = list(combinations(range(7), p))
    idx = {s: i for i, s in enumerate(subs)}
    n = len(subs)
    M = np.zeros((n, n))
    for col, S in enumerate(subs):
        for pos, j in enumerate(S):
            for m in range(7):
                if abs(A[m, j]) < 1e-15:
                    continue
                if m in S and m != j:
                    continue            # repeated index -> 0
                newlist = list(S); newlist[pos] = m
                # sort, track sign
                arr = newlist
                sign = 1
                a = arr[:]
                for i in range(len(a)):
                    for k in range(len(a) - 1 - i):
                        if a[k] > a[k + 1]:
                            a[k], a[k + 1] = a[k + 1], a[k]; sign = -sign
                T = tuple(a)
                if len(set(T)) != p:
                    continue
                M[idx[T], col] += sign * A[m, j]
    return M


def casimir_spectrum(g2on, p):
    """Casimir = sum_a rho_p(T_a)^2 on Lambda^p; return sorted |eigenvalues| with mults."""
    subs_dim = math.comb(7, p)
    C = np.zeros((subs_dim, subs_dim))
    for A in g2on:
        Mp = wedge_action(A, p)
        C += Mp @ Mp
    ev = np.linalg.eigvalsh(C)             # <= 0 (sum of squares of antisym)
    ev = np.round(-ev, 4)                  # |Casimir| values, rounded
    vals, counts = np.unique(ev, return_counts=True)
    return list(zip(vals.tolist(), counts.tolist()))


# ---------------------------------------------------------------------------
def main():
    banner("CHEAP PRE-CHECK: G2-equivariant finite-sum loop premise")
    phi = assoc_3form()
    # antisymmetry sanity
    asym = np.max(np.abs(phi + np.transpose(phi, (1, 0, 2))))
    print(f"  associative 3-form phi built from Fano table; antisymmetry residual = {asym:.1e}")

    gens, dimg2 = build_g2(phi)
    print(f"  dim g2 = dim Stab_so(7)(phi) = {dimg2}   (expect 14)")
    g2on = gram_orthonormalize(gens)

    print("\n  Casimir spectra (|eigenvalue|: multiplicity) on Lambda^p(R^7):")
    spec = {}
    for p, dim in [(1, 7), (2, 21), (4, 35)]:
        sp = casimir_spectrum(g2on, p)
        spec[p] = sp
        mults = {v: c for v, c in sp}
        print(f"    Lambda^{p} (dim {dim}): { {round(v,3): c for v,c in sp} }")

    # interpret
    banner("Interpretation")
    # fundamental-7 casimir value (from Lambda^1, the single nonzero block)
    c7 = spec[1][0][0]  # the 7's Casimir (Lambda^1 is the irreducible 7)
    print(f"  Casimir(7, fundamental) = {c7:.3f}  (from Lambda^1, irreducible)")
    # 14 from Lambda^2: the multiplicity-14 eigenvalue
    c14 = next((v for v, c in spec[2] if c == 14), None)
    print(f"  Casimir(14, adjoint)    = {c14}  (the mult-14 block, lives in Lambda^2)")
    # Lambda^4 multiplicities
    l4 = {c: v for v, c in spec[4]}
    print(f"  Lambda^4 multiplicities = { {c: round(v,3) for v,c in spec[4]} }")
    has_singlet = any(abs(v) < 1e-6 for v, c in spec[4])
    has_14_in_l4 = any(c == 14 for v, c in spec[4])
    print(f"\n  Lambda^4 has a singlet (the coassociative 4-form)? {has_singlet}")
    print(f"  Lambda^4 contains a 14 (adjoint)?                    {has_14_in_l4}")
    print(f"""
  => Lambda^4(R^7) = 1 + 7 + 27 (mult {{1,7,27}}); the 14 (adjoint) is in
     Lambda^2 = 7 + 14, NOT in Lambda^4.  v59 v_eff_loop's "N_X = 14 from Lambda^4"
     is a MISATTRIBUTION -- there is no 14-dim G2 subspace of Lambda^4.
""")

    banner("(Q2) Can any N_X change the ratio c3/c6?  No.")
    print("""  v59's loop formula:  V_eff^X(phi) = -(N_X/64pi^2) Tr[M_X(phi)^2 log(...)].
  N_X is an OVERALL, phi-INDEPENDENT prefactor.  A constant times V scales c3 and
  c6 by the SAME factor, so it CANCELS in cos(3phi) = -c3/(4c6).  Therefore no
  G2-equivariant counting -- 14, 27, 7, whatever -- can un-suppress c6 relative
  to c3.  The c6 suppression is a GENERATION-Z3 harmonic effect (cos 6phi appears
  at higher order in t than cos 3phi), ORTHOGONAL to the G2/internal index.""")
    # demonstrate cancellation numerically on a toy V
    g = np.linspace(0, 2 * math.pi, 4000, endpoint=False)
    V = np.cos(3 * g) * 1.98 + np.cos(6 * g) * 0.0643      # toy c3,c6 (from loop_ratio)
    def coeff(Vv, n): d = g[1]-g[0]; return (1/math.pi)*np.sum(Vv*np.cos(n*g))*d
    r1 = -coeff(V, 3) / (4 * coeff(V, 6))
    r2 = -coeff(14 * V, 3) / (4 * coeff(14 * V, 6))   # scale by N_X=14
    print(f"\n  ratio with N_X=1 : {r1:.4f};  with N_X=14 : {r2:.4f}  (identical: {abs(r1-r2)<1e-9})")

    return {"dimg2": dimg2, "c7": c7, "c14": c14, "spec4": spec[4],
            "has_14_in_l4": has_14_in_l4, "has_singlet": has_singlet,
            "ratio_invariant": abs(r1 - r2) < 1e-9}


def run_checks(R):
    banner("EMBEDDED SELF-CHECKS")
    ok = True

    def check(name, cond):
        nonlocal ok
        ok = ok and bool(cond)
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")

    check("dim g2 = 14 (derivation algebra from Fano table)", R["dimg2"] == 14)
    mults4 = sorted(c for v, c in R["spec4"])
    check("Lambda^4 decomposes as {1,7,27} (mults)", mults4 == [1, 7, 27])
    check("Lambda^4 contains NO 14 (adjoint is in Lambda^2, not Lambda^4)",
          not R["has_14_in_l4"])
    check("Lambda^4 has exactly one singlet (coassociative 4-form)", R["has_singlet"])
    check("Casimir(adjoint 14) exists in Lambda^2", R["c14"] is not None)
    check("ratio c3/c6 is invariant under overall N_X factor (N_X cannot help)",
          R["ratio_invariant"])

    print(SEP)
    print("ALL CHECKS PASSED (premise misattributed AND N_X cannot help -> build NOT warranted)"
          if ok else "SOME CHECKS FAILED")
    print(SEP)
    return ok


if __name__ == "__main__":
    R = main()
    print()
    ok = run_checks(R)
    import sys
    sys.exit(0 if ok else 1)
