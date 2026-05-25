#!/usr/bin/env python3
"""
sedenion_s3.py — analysis of φ = Q/3 in the SEDENION S₃ factor (not G₂).

Per the literature (Gresnigt et al. 2024, arXiv:2407.01580; Eakin–Sathaye), the three
fermion generations are grounded in  Aut(𝕊) = Aut(𝕆) × S₃ = G₂ × S₃: the S₃ factor — the
"/3" — is NOT present in the octonions (Aut(𝕆)=G₂) but appears in the sedenions 𝕊.

This script builds 𝕊, the explicit S₃ generators (ε order-2, ψ order-3 with √3), verifies
they are genuine automorphisms forming S₃, and analyses what the S₃ does — and does NOT —
fix for the Brannen generation phase φ = Q/3 = 2/9.

Conclusion (computed): the generation Z₃ = ψ IS a genuine sedenion automorphism, so the
"/3" is structural/emergent.  But ψ's intrinsic phases are {0, ±2π/3} (cube roots, π-rational),
whereas the coupling phase φ=2/9 is NOT π-rational (lean/PhaseExclusions.koide_not_pi_rational).
So the S₃ fixes the generation STRUCTURE, not the phase VALUE 'Q': the mystery splits into a
now-explained "/3" and a residual "Q".
"""
import numpy as np

# ---- octonions (fixed table) -------------------------------------------------
_T = [[(1,0),(1,1),(1,2),(1,3),(1,4),(1,5),(1,6),(1,7)],[(1,1),(-1,0),(1,3),(-1,2),(1,5),(-1,4),(-1,7),(1,6)],
      [(1,2),(-1,3),(-1,0),(1,1),(1,6),(1,7),(-1,4),(-1,5)],[(1,3),(1,2),(-1,1),(-1,0),(1,7),(-1,6),(1,5),(-1,4)],
      [(1,4),(-1,5),(-1,6),(-1,7),(-1,0),(1,1),(1,2),(1,3)],[(1,5),(1,4),(-1,7),(1,6),(-1,1),(-1,0),(-1,3),(1,2)],
      [(1,6),(1,7),(1,4),(-1,5),(-1,2),(1,3),(-1,0),(-1,1)],[(1,7),(-1,6),(1,5),(1,4),(-1,3),(-1,2),(1,1),(-1,0)]]
def omul(x,y):
    z=np.zeros(8)
    for i in range(8):
        if x[i]:
            for j in range(8):
                if y[j]:
                    s,k=_T[i][j]; z[k]+=s*x[i]*y[j]
    return z
def oconj(x): r=-x.copy(); r[0]=x[0]; return r

# ---- sedenions via Cayley–Dickson (a,b)(c,d)=(ac − d̄b, da + bc̄) -------------
def smul(u,v):
    a,b,c,d=u[:8],u[8:],v[:8],v[8:]
    return np.concatenate([omul(a,c)-omul(oconj(d),b), omul(d,a)+omul(b,oconj(c))])
def e(i): z=np.zeros(16); z[i]=1; return z

# ---- the S₃ automorphisms (Gresnigt et al.) ----------------------------------
def eps_v(u): return np.concatenate([u[:8], -u[8:]])                     # ε: A+Bs8 -> A-Bs8
def psi_v(u):                                                            # ψ (order 3, the √3 map)
    A,B=u[:8],u[8:]; As,Bs=oconj(A),oconj(B)
    return np.concatenate([0.25*(A+3*As+np.sqrt(3)*(B-Bs)),
                           0.25*(B+3*Bs-np.sqrt(3)*(A-As))])
E=np.column_stack([eps_v(e(i)) for i in range(16)])
P=np.column_stack([psi_v(e(i)) for i in range(16)])

def is_auto(fmat,fv):
    return all(np.allclose(fv(smul(e(i),e(j))), smul(fmat@e(i),fmat@e(j)))
               for i in range(16) for j in range(16))

def main():
    print("="*70); print("1. The sedenion S₃ = ⟨ε, ψ⟩ is genuine (the home of the '/3')"); print("="*70)
    print(f"  ε² = I           : {np.allclose(E@E,np.eye(16))}")
    print(f"  ψ³ = I           : {np.allclose(np.linalg.matrix_power(P,3),np.eye(16))}")
    print(f"  εψ = ψ²ε  (S₃)   : {np.allclose(E@P, np.linalg.matrix_power(P,2)@E)}")
    print(f"  ε is an automorphism of 𝕊 : {is_auto(E,eps_v)}")
    print(f"  ψ is an automorphism of 𝕊 : {is_auto(P,psi_v)}")
    print("  => the generation Z₃ = ⟨ψ⟩ is a GENUINE sedenion automorphism (not external).")
    print("     The '/3' of φ=Q/3 is structural: it IS emergent from (extended) octo-space.")

    print("\n" + "="*70); print("2. ψ acts as a 120° rotation of the two 7-dim Im-octonion blocks"); print("="*70)
    print("  ψ on (Im A, Im B): [[-1/2, √3/2],[-√3/2,-1/2]]  (rot by -2π/3, the √3 = sin120°)")
    print("  real parts e₀, e₈ are fixed.  So ⟨ψ⟩ realises the Z₃ Fourier split 1, ω, ω².")
    ev = np.linalg.eigvals(P)
    print(f"  ψ eigenphases / π : {np.unique(np.sort(np.round(np.angle(ev)/np.pi,4)))}")
    print("  => the SYMMETRY's intrinsic phases are {0, ±2π/3}: cube roots, π-RATIONAL.")

    print("\n" + "="*70); print("3. What the S₃ does NOT fix: the phase VALUE 2/9"); print("="*70)
    print("  The Brannen generation phase is the COUPLING phase φ = arg(ξ) in M = a(I+ξψ+ξ̄ψ²),")
    print("  NOT an eigenphase of ψ.  ψ-covariance only makes M circulant — any ξ is allowed.")
    print(f"  φ = 2/9 = {2/9:.5f} rad;  ψ's phases are multiples of 2π/3 = {2*np.pi/3:.5f}.")
    print("  And 2/9 is NOT π-rational (lean: PhaseExclusions.koide_not_pi_rational), so it is")
    print("  not ANY symmetry/holonomy phase.  => the S₃ does NOT fix φ=2/9.")

    print("\n" + "="*70); print("4. The clean decomposition of the mystery"); print("="*70)
    print("""  φ = Q / 3 splits into:
    * '/3' = N_generations = the sedenion S₃ automorphism factor  -> EMERGENT / structural (✓).
    * 'Q' = 2/3 = the phase MAGNITUDE (3φ) = the coupling value     -> RESIDUAL, not symmetry-
            fixed, provably not π-rational, hence not a holonomy/symmetry phase.
  So proceeding in the S₃ factor CONFIRMS the '/3' is emergent and ISOLATES the open part to
  the single question: why is the coupling magnitude exactly 3φ = Q = dimG₂/dimSpin7 = 2/3?
  That residual is orthogonal to the generation symmetry — a dynamical/mass-sector input.""")

if __name__ == "__main__":
    main()
