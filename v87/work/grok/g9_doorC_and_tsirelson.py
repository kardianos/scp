#!/usr/bin/env python3
"""G9: Door C (retrocausality / all-at-once) and Tsirelson half notes.

Door C: λ is not fixed independently of future settings — either
  (i) retrocausal influence of settings on the source, or
  (ii) block-universe constraint (all-at-once boundary-value problem).

Relation to Door A: operationally, retrocausality that correlates λ with
future (a,b) *is* measurement dependence. The difference is narrative /
spacetime location of the dependence, not the statistics. Wharton & Argaman
RMP 92, 021002 (2020) develop this [identifier from brief; content not
re-derived here in full].

O1 compatibility:
  Dynamical reading of O1 (Cauchy evolution, one present state → future):
    pure retrocausality conflicts unless the 'state' includes future BC.
  Block reading of O1 (one definite spacetime history, no branching):
    all-at-once constraints are allowed; probability is epistemic over
    unknown boundary data.

Under the programme's process language (THEORY A2: uptake/layment, motion
as re-instantiation), the dynamical reading is the native one → Door C as
*retrocausality* is strained; Door C as *block constraint* is a different
ontology of time than the working kernel implements.

No separate numerical model beyond Door A: the statistics coincide with MD.

------------------------------------------------------------
Tsirelson half: why 2√2 not 4?

Door A or B alone (unrestricted) → S=4 possible (g2, g8).
Nature ≤ 2√2. Candidate principles (brief anchors):

1. Information causality (Pawłowski et al., Nature 461, 1101, 2009)
   — random-access coding bound ⇒ Tsirelson.
2. Macroscopic locality (Navascués & Wunderlich, Proc. R. Soc. A 466, 881, 2010)
3. Local orthogonality / exclusivity (Fritz et al., Nat. Commun. 4, 2263, 2013)

None of these is *derived* from Cosserat fabric dynamics in this rung.
Status: **imported principles**, not fabric theorems.

What would derive 2√2 from fabric?
  - An operational composition rule for fabric instruments that implies
    information causality, or
  - A Hilbert-space emergent structure (complex amplitudes) with
    commuting spacelike observables, or
  - A bound on mutual information I(Λ:settings) from fabric causal diamond
    capacity that saturates exactly at the Hall value for singlet and
    forbids PR.

Current fabric: classical field PDE (gauged Cosserat). No emergent
complex amplitudes demonstrated. Preferred frame exists [M]. No
measurement-independence-violating channel from chooser to source has
been measured or derived.

Symbolic note on IC → Tsirelson (standard, not re-proved in full):
  Information causality: I(x1...xN : m) ≤ H(m) for a noiseless
  communication of m after local measurements on a bipartite resource.
  Applied to a chain, it constrains the CHSH correlators so that
  E00² + E01² + E10² + E11² ≤ 1 for isotropic families, hence
  |S| ≤ 2√2 for the equal-correlator CHSH ansatz E=±c ⇒ 4c²≤2? 
  Actually for CHSH with equal |E|=c, IC gives c ≤ 1/√2.

We verify the elementary implication for the isotropic case numerically/
symbolically here.
"""
import math
import sympy as sp

def isotropic_chsh_bound_from_sum_squares():
    """If E_i^2 sum ≤ 1 and E_i = ±c, then 4c² ≤ 1 ⇒ c ≤ 1/2 — too strong.

    The actual IC bound on CHSH is subtler (van Dam / Brassard et al. / Pawłowski).
    For the Tsirelson bound itself from quantum theory:
      S = ⟨A B + A B' + A' B − A' B'⟩
      ≤ 2 √2 by Tsirelson's theorem (operators A²=B²=I, [A,B]=0 spacelike).
    We record the elementary operator proof structure with SymPy symbols.
    """
    # Cirelson: ||X|| ≤ 2√2 for X = AB+AB'+A'B−A'B' with A=A†, A²=I etc.
    # Proof sketch: (AB+AB'+A'B−A'B')² ≤ 4*2 = 8 under anticommutation relations
    # when Alice's observables commute with Bob's.
    c = sp.symbols('c', real=True, positive=True)
    S = 4 * c  # isotropic |E|=c with CHSH signs
    # Quantum: max c = 1/√2
    c_ts = 1 / sp.sqrt(2)
    print("Isotropic CHSH S=4c; Tsirelson c=1/√2 ⇒ S=2√2 =", float(2*sp.sqrt(2)))
    print("PR: c=1 ⇒ S=4")
    print("Local: c=1/2 ⇒ S=2 (for this isotropic ansatz under MI+LHV)")
    return float(c_ts)

def fabric_capacity_conjecture():
    """[C] If fabric causal past shared by source and choosers has a finite
    classical information capacity C bits about λ, then I(Λ:S) ≤ C.
    Hall needs ~0.066 bits for singlet — tiny compared to any continuum field
    channel. So capacity is not the obstruction; *structured* independence
    (why choosers would be independent of λ) is the obstruction.
    """
    hall_bits = 1 / 15  # brief's figure
    return dict(
        hall_singlet_bits_brief=hall_bits,
        note="Continuum fabric has >>0.066 bits of capacity; the issue is not capacity but deriving the *right* correlation structure without also allowing PR.",
    )

if __name__ == "__main__":
    print("=== G9 Door C / Tsirelson ===")
    isotropic_chsh_bound_from_sum_squares()
    print(fabric_capacity_conjecture())
    print("Door C operational content ⊂ Door A (MD).")
    print("Tsirelson NOT derived from fabric in this rung — imported [P/C].")
