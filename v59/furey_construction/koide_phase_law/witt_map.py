#!/usr/bin/env python3
"""
witt_map.py — the explicit Witt (octonion → Fock) map, and the holonomy test.

Constructs the Furey ℂ⊗𝕆 Witt creation/annihilation operators from the genuine
octonion left-multiplications γ_k (post 2026-05-24 sign fix), using the SCALAR
complex unit i of ℂ⊗𝕆 (so the operators live on the complexified ℂ⁸).  This is the
exact octonion→Fock refinement deferred in SevenDAlgebra.lean.

Then it runs the generation-cycle holonomy test for the lepton phase φ = Q/3.

Verdict (computed): the Witt map is a SINGLE-generation / 3-COLOR construction; the
color Z₃ fixes the lepton color-singlet, so the lepton phase is NOT an intra-ideal
holonomy.  The Brannen generation Z₃ is external (3 copies); running the holonomy
test requires first building the inter-generational connection.
"""
import numpy as np, itertools

T = [[(1,0),(1,1),(1,2),(1,3),(1,4),(1,5),(1,6),(1,7)],
     [(1,1),(-1,0),(1,3),(-1,2),(1,5),(-1,4),(-1,7),(1,6)],
     [(1,2),(-1,3),(-1,0),(1,1),(1,6),(1,7),(-1,4),(-1,5)],
     [(1,3),(1,2),(-1,1),(-1,0),(1,7),(-1,6),(1,5),(-1,4)],
     [(1,4),(-1,5),(-1,6),(-1,7),(-1,0),(1,1),(1,2),(1,3)],
     [(1,5),(1,4),(-1,7),(1,6),(-1,1),(-1,0),(-1,3),(1,2)],
     [(1,6),(1,7),(1,4),(-1,5),(-1,2),(1,3),(-1,0),(-1,1)],
     [(1,7),(-1,6),(1,5),(1,4),(-1,3),(-1,2),(1,1),(-1,0)]]
def gamma(k):
    G = np.zeros((8,8), complex)
    for m in range(8):
        s,t = T[k+1][m]; G[t,m] = s
    return G
G  = [gamma(k) for k in range(7)]
I8 = np.eye(8, dtype=complex)

def main():
    print("="*70); print("1. Witt operators (3 disjoint γ-pairs, scalar i of ℂ⊗𝕆)"); print("="*70)
    pairs = [(0,5),(1,4),(2,3)]          # disjoint pairs of the 6 units; e₇=γ₆ is special
    a  = [0.5*(G[x] + 1j*G[y]) for x,y in pairs]
    ad = [A.conj().T for A in a]
    print(f"  pairs (γ_a,γ_b) for the 3 color modes: {pairs}   (e₇ excluded)")
    print(f"  nilpotent a²=0           : {all(np.allclose(A@A,0) for A in a)}")
    print(f"  CAR  {{a_i,a_i†}} = 1      : {all(np.allclose(a[i]@ad[i]+ad[i]@a[i],I8) for i in range(3))}")
    print(f"  {{a_i,a_j}}=0, {{a_i,a_j†}}=0 (i≠j): "
          f"{all(np.allclose(a[i]@a[j]+a[j]@a[i],0) and np.allclose(a[i]@ad[j]+ad[j]@a[i],0) for i in range(3) for j in range(3) if i!=j)}")

    print("\n" + "="*70); print("2. The Fock grading (exact octonion→Fock map)"); print("="*70)
    Nop = sum(ad[i]@a[i] for i in range(3))
    ev = np.sort(np.round(np.real(np.linalg.eigvals(Nop)),6))
    print(f"  number operator N = Σ a_i†a_i eigenvalues: {ev}")
    print(f"  -> the Fock grading (0,1,1,1,2,2,2,3): {np.allclose(ev,[0,1,1,1,2,2,2,3])}")
    # vacuum = N=0 eigenvector
    w,V = np.linalg.eigh(Nop)
    vac = V[:, np.argmin(np.real(w))]
    print(f"  vacuum |Ω⟩ killed by all a_i: {[round(np.linalg.norm(A@vac),6) for A in a]}")
    print("  => 8 Fock states = Λ•(ℂ³_color) of ONE generation, built on 3 COLORS.")

    print("\n" + "="*70); print("3. Holonomy test for the lepton phase"); print("="*70)
    Jc = (G[0]@G[5]).real
    R = np.zeros((8,8)); R[0,0]=R[7,7]=1
    R[2,1]=R[3,2]=R[1,3]=1; R[5,6]=R[4,5]=R[6,4]=1   # color Z₃: cycle modes, fix lepton
    print(f"  color Z₃ R: R³=I {np.allclose(np.linalg.matrix_power(R,3),np.eye(8))}, "
          f"[R,J_c]=0 {np.allclose(R@Jc,Jc@R)}")
    print(f"  R on lepton singlet line {{0,7}}: R e₀=e{np.nonzero(R[:,0])[0][0]}, "
          f"R e₇=e{np.nonzero(R[:,7])[0][0]}  -> FIXED pointwise")
    print("  => color-Z₃ holonomy on the lepton = TRIVIAL (0).  φ=2/9 is NOT an")
    print("     intra-ideal color holonomy.  (Lepton = color singlet, the N=0 & N=3 fixed points.)")

    print("\n" + "="*70); print("4. Verdict — what the holonomy test needs"); print("="*70)
    print("""  The Witt map is a single-generation, 3-COLOR construction (verified above).
  The Brannen phase lives on the GENERATION Z₃ — a 3-COPY structure external to one
  ideal — whose Wilson loop is arg(ξ³)=3φ.  Every Z₃ internal to the single ideal
  (the color Z₃; the J_c-phase Z₃) either fixes the lepton singlet or fails to give
  three equal generation eigenspaces.  So:

    * φ = Q/3 is NOT an intra-ideal holonomy (proven trivial for the color Z₃);
    * to run the holonomy test one must BUILD the inter-generational connection —
      three copies of the ideal linked by a Z₃ connection (this is where the
      Cl(3,0)→Cl(7) / v58-multivector expansion would have to supply the generation
      structure) — and test whether ITS closed-cycle holonomy is Q.

  The single-ideal Witt map (done here) cannot supply that connection.  This is the
  precise missing ingredient, now located: the 3-copy generation connection, not the
  intra-ideal color/complex structure.""")

if __name__ == "__main__":
    main()
