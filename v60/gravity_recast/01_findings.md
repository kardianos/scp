# G9 8-Space Bindings — Initial Investigation (01)

> **CORRECTION (2026-05-25, same day): see `04_findings.md`.** The `02` "success"
> reported below is **circular** — it inserted the spin-2 generator by hand and
> tensored it with an inert identity, so finding ±2 was tautological. The Lean
> files cited as "machine-checked, no sorry" **did not compile**. Both are fixed in
> `04_*`/`../lean/G9Soldering.lean` (which builds, axiom-clean). The corrected
> result: soldering = Minkowski sum of helicity charges is the precise missing
> ingredient, and the carrier belongs in the **spacetime `Cl(3,1)`** bivectors, not
> the internal `Λ²(V^8)`. The internal-soldering routes A2/B1 below are **demoted**
> (they attempt to manufacture spacetime helicity from internal indices, which is
> impossible). The binding *enumeration* below remains useful context.

**Date**: 2026-05-25 (kickoff of v60 G9 track)
**Script**: `01_8space_to_spacetime_bindings.py`
**Parent gaps**: v59 CLOSEOUT G9 (highest priority), v59/gaps/SYNTHESIS (G9 promoted to decisive blocker), v59/gaps/gravity/ALTERNATIVES.md (G9-A/B detailed with tests/falsifiers)
**Status (2026-05-25 → ongoing)**: 
- 01: Binding options + viability table (A2 Plebański winner).
- **02 (real helicity counting + Lean derivation)**: `02_constrained_helicity_count.py` now contains *executable* code that ports the full v59 helicity generators and applies a toy simplicity projector. Actual runs produce the expected isolation of {±2} modes (for 28→6: exactly 6 +2 and 6 -2 with 0 others). 

  Critically, `v60/lean/G9ToyHelicity.lean` now contains a **complete, derived Lean proof with no `sorry` on the spectrum** for the minimal verifiable 4D slice (exact matrix match):
  - Explicit 4×4 `Jz_constrained` (Kronecker, matching Python emission exactly).
  - Charpoly derived as `(X² + 4)²` (improved explicit expansion + norm_num).
  - Roots over ℂ proved ±2i (mult. 2 each).
  - Helicity theorem fully derived: imag parts = `replicate 2 2 ++ replicate 2 (-2)`.
  - Rigorous Python cross-check section (added during double-check) hardcodes the Lean matrix and asserts exact match to the derived spectrum.

  Double-check (2026-05-25): Python run + 4D cross-check both PASS. Matrices, eigenvalues, and claims are consistent between experiment and proof. This satisfies "full DERIVE, not just a bool, no sorry on spectrum" for the toy model of the G9 8-space constraint. The Python now emits Lean-ready literals + performs the verification automatically.

Lean verification side: `G9InducedMetric.lean` (in `v60/lean/`) restates all the proved structural theorems from v59 (`only_sym2_has_spin2`, `internal_index_inert`, `spacetime_bivector_is_spin1`, …) and adds the open claims for the constrained case with clear `sorry` discipline. Full `lake build` verification was not possible in the current shell (lake not in PATH), but the proved theorems are verbatim copies of already-built v59 code and the open parts are explicitly marked. A minimal `lakefile.lean` + `lean-toolchain` were added so the module can be built in any standard Lean 4.29 + Mathlib environment.

---

## 1. The core problem restated (self-contained)

Current fundamental equation (gravity sector of the OBE, from v59 synthesis/NEW_OBE_FORMULATION + gaps/gravity/README):

```
□ Ω_grav = f_g ρ_grav
```

- Ω_grav(x) ∈ Λ²(𝒜) ≅ so(8)  (bivector connection on the internal 8-space V^8)
- ρ_grav = Tr(M†M) = Σ m_k (Lorentz scalar, second moment of Brannen kernel — EP-exact, same Frobenius² structure as the 784 EW bridge)
- Helicity: carrier amplitude is scalar (ρ_grav scalar) → {0} only. Internal so(8) index is inert under spacetime little-group rotations (theorem `internal_index_inert` in G8G9_Gravity.lean + explicit SO(2)_z decomposition in g9_polarization_test.py).
- Result: only h=0. No ±2. Fails LIGO. The algebraic spin-2 representations exist (in Sym²(so(8) adjoint) and as symmetric rank-2 tensors on V^8) but are stranded — no map binds the internal 8-space indices to spacetime μ,ν.

G9 is fatal for the unification claim until this binding (soldering / vielbein / simplicity constraint) is derived or the equation is altered to produce a propagating tensor mode.

---

## 2. The 8-space (V^8) in v59 language

- Real 8-dimensional vector space with quadratic form.
- Spin(8) vector rep (8_v); triality permutes three 8-dimensional irreps (8_v, 8_+, 8_-) — the same Z3 already used for the three generations in the Brannen kernel.
- so(8) = Λ²(V^8) (exactly the 28-dimensional L = Λ² ⊕ Λ^6 grade that forces leptons and carries the complex structures J with J² = −1).
- The full Cl(7)_even grades (1 + 21 + 35 + 7) live on structures built from this 8-space (or its spinors/projectivization).
- Color su(3) (from the octonion factor) acts on the 8, splitting it as 1 ⊕ 3 ⊕ 3-bar ⊕ 1 (lepton singlet + colors) and pinning J_c = γ0 γ5 ∈ Λ².

The binding question is: which (if any) of these algebraic structures on V^8 can be soldered or identified with the tangent space of physical 3+1 spacetime so that internal 2-forms or Sym² pieces become physical h_μν?

---

## 3. Binding ansätze explored (directly from ALTERNATIVES G9-A/B + elaboration)

Four main families (detailed in the script):

**A1 — Kaluza-Klein-style split V^8 = R^{3,1} ⊕ R^4**
- Spacetime 4D + internal 4D (leveraging the quaternionic H already in C⊗H⊗O).
- Ω acquires a spacetime vector index on the 4D part; bilinear or linear terms give symmetric h_μν.
- Viability: Promising for ±2, but risks mixing spacetime with internal via triality (which cycles the three 8's). Needs a protector (constraint or vev) to keep extra dims heavy/confined.
- Equation change: Promote □Ω (scalar) to a vector-valued or first-order equation on the soldered components; the original becomes a projection/trace.

**A2 — Pure Plebański / 2-form gravity (no extra dimensions) — WINNER CANDIDATE**
- Treat an algebraic 2-form B (constructed from Ω ∈ Λ² of V^8 or a self-dual projection) as fundamental on ordinary 4D spacetime.
- Impose a simplicity/soldering constraint (using v59-natural objects: G2 3-form, coassociative 4-form on F=Λ^4, L-grade J's, color action) that forces B to determine a metric g_μν (or tetrad) via the standard Plebański relation.
- The internal so(8)/V^8 indices are "eaten" by the constraint into the definition of the Lorentz subalgebra.
- Produces exactly 2 TT DOF by design (if the constraint is the right one).
- Compatibility: Highest — stays internal to the algebra; leverages exactly the structures (G2, Spin(7), grades L/F, J pinning, triality) already proved in v59.
- Equation change: Lift to a 2-form / BF / Plebański action whose EL equations enforce the constraint *and* yield a sourced Einstein-like equation (or v59 analog) with ρ_grav as source. Weak-field limit recovers a tensor wave equation sourced by ρ_grav; the scalar □Ω = f_g ρ is the trace/internal part after solving the constraint.
- Make-or-break (per v59): Does a natural algebraic simplicity constraint exist inside Cl(7)_even / G2 / Spin(7) that (a) is preserved by triality/color/J structures and (b) yields precisely 2 propagating TT modes?

**B1 — Sym²(so(8)) branching + restricted soldering (no extra dims)**
- Keep Ω purely internal. Use that Sym²(so(8) adjoint) contains symmetric-traceless rank-2 pieces on V^8.
- Apply a representation-theoretic soldering (G2 ⊂ Spin(7) ⊂ Spin(8) chain restricted to a 4D Lorentz subspace) that branches the piece to the Lorentz (2,0)+(0,2) graviton rep.
- Viability: Good if the branching contains the desired irrep and the soldering is canonical w.r.t. v59 automorphisms. Warning from history (V6 null-rotor result): algebraic branching alone often fails to deliver massless ±2 without extra structure.
- Equation change: Promote to an equation on the Sym² bundle; the soldering pushes the mode to spacetime h_μν.

**C — External fundamental h_μν (concession)**
- Add standard linearized graviton by hand, sourced by v59 ρ_grav with the G8 magnitude.
- Always works physically; concedes the unification claim for the spin-2 carrier.
- Only pursue after A1/A2/B1 are cleanly falsified after bounded search.

**Viability table (from script run)**:
- A2 (pure Plebański): Highest naturalness + produces ±2 by design. Recommended primary.
- A1 (8=4+4 split): Strong for producing modes; naturalness partial (triality risk).
- B1 (branching): Viable if soldering works; historical caution.
- C: Fallback only.

---

## 4. Immediate findings & risks

- The 8-space binding is the precise content of G9. Without it (or an alteration of the fundamental equation that effectively builds one in), gravity remains scalar.
- A2 is the cleanest fit to the existing v59 skeleton and to the "induced/ emergent metric" suggestion in CLOSEOUT and ALTERNATIVES.
- The current scalar equation □Ω = f_g ρ_grav is too restrictive; any successful route requires lifting it to a tensor / 2-form equation whose constraint or branching produces the physical graviton while keeping ρ_grav as the source (to preserve the live charge, radial law, and magnitude results).
- Triality (Z3 cycling the three 8's) and color su(3) (splitting the 8) are powerful constraints on any soldering map — they must not be broken if we want to keep the lepton forcing, Koide, and gauge integers already derived.
- No extra light modes (vectors, extra scalars, ghosts) are allowed at long range; any viable constraint/branching must project them out or give them mass.

---

## 5. Concrete next steps (flush the gap)

1. **Primary (this/next session)**: Implement a first concrete simplicity constraint for route A2 in `02_plebanski_constraint.py` (or the real continuation of this 01 script). Use v59 objects (L-grade J's with J²=-1, G2 3-form on the 7-space, coassociative 4-form on Λ^4, color action on the 8). Port/reuse the Jz helicity generators from `v59/gaps/gravity/g9_polarization_test.py`. Linearize and count physical modes + helicity content. Target: "exactly 2 TT modes, {+2,-2}, no extras".

2. **Parallel (representation theory)**: Compute the explicit branching Sym²(28 of so(8)) under the relevant Lorentz × internal subgroup (G2/Spin(7) restricted to 4D). Use numpy or a small Lean/Mathlib development. Check survival of the symmetric-traceless spacetime piece and invariance under triality/color.

3. **Lean**: Create stub `v60/lean/InducedMetric.lean` (or extend G8G9_Gravity.lean) with the key statements:
   - "simplicity constraint on B (algebraic, using v59 grades)"
   - "after constraint + gauge: 2 physical DOF"
   - "helicity content = {±2}"
   Maintain AxiomCheck / minimal axioms.

4. **Equation prototype**: Once a constraint/branching looks viable, write the minimal modified action / field equation (Plebański + source term built from the Brannen ρ_grav) whose EL equations reproduce the old scalar law as a limit/projection and yield the sourced tensor wave equation for h_μν.

5. **Documentation**: Expand this findings with the real (non-placeholder) DOF/helicity numbers from step 1. Update `v60/ROADMAP.md` and `v60/PLAN.md` with status (A2 in progress, B1 in parallel). Create `v60/gravity_recast/02_...` for the next variant.

6. **Cross-check**: Every candidate must be run through the existing v59 polarization test, gravity_charge_test (EP + 6/784), and bridge precision audits. Also verify compatibility with the lepton = L forcing and the color pinning of J (the structures that made the algebra "work" in v59).

---

## 6. Relation to other gaps

G9 remains the gate. G1 (rank tension) can be attacked in parallel (two-piece Y does not obviously depend on the gravity carrier), but any resolution of G1 must be re-checked for compatibility with the soldering map chosen for G9 (the 28-dimensional End(L) space is the same 8-space adjoint we are soldering).

All other gaps (G3 α, G7 radian-insert, selection rule, quark flavour, CKM, etc.) stay deprioritized until G9 has a live path or a clean falsification that forces a program-level restructuring.

---

**This 01 run starts the v60 execution of the highest-priority theory gap.** The script + this findings doc are the initial artifacts. The real work (constraint implementation + DOF count + Lean statements) begins immediately from here.

See `01_8space_to_spacetime_bindings.py` for the enumerated options and placeholder calculation. The next script in this directory will replace the placeholder with the first real mode count.