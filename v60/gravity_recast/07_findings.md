# C3 closed — the unified algebra is v59's own Cl(3,1) ⊗ Cl(7)_even

**Date**: 2026-05-25 (builds out the v59 residual flagged in `05_findings.md` / `06`)
**Artifacts**:
- `07_unified_algebra.py` (runs clean; all commutators exactly 0)
- `../lean/G9Unification.lean` (compiles clean; axiom-only — general theorem)
- `06_lorentz_commutant.py` (the Schur no-go that forces the factorization)

## The fork, resolved from v59's own documents

C3 turned on one question: is the spacetime `ℍ` a **separate tensor factor** of
`ℂ⊗ℍ⊗𝕆`, or a quaternion **subalgebra inside** the internal `𝕆`? If the latter,
the Lorentz would be part of the internal symmetry and the whole split would be
illusory. v59's own sources settle it:

- `furey_construction/02_findings.md`: **the `ℍ` factor in `ℂ⊗ℍ⊗𝕆` is `SU(2)_L`**
  (weak isospin) + the generation structure — it is *internal*, not spacetime.
- `synthesis/SYNTHESIS.md §8`: the full algebra is **`Cl(spacetime) ⊗ Cl(internal)`**;
  "v58 `M` lives in `Cl(7)` (internal/octonion); the EM/spacetime scale comes from
  `Cl(3,1)` (spacetime) … the two algebras are factorized."
- `synthesis/SYNTHESIS.md §9`: `G_e = (dim Spin(7)/dim Cl(3,1))·α²¹ = (21/16)α²¹`,
  i.e. the prefactor is the **internal/spacetime dimension ratio** of this split.

So spacetime Lorentz lives in a **separate `Cl(3,1)` factor**; the entire internal
`Cl(7)_even` (including its `ℍ = SU(2)_L`) is the other factor. The gap was that
v59 *stated* this factorization but never built it or checked the commutation.

## What was built (`07_unified_algebra.py`)

Explicit `Cl(3,1) ⊗ Cl(7)_even` on the module `Dirac(4) ⊗ octonion(8) = 32`:

1. **Cl(3,1)** — real Dirac `γ^μ` (metric `diag(+,−,−,−)` recovered from
   `{γ^μ,γ^ν}/2`), the 6 bivectors `σ^{μν} = ¼[γ^μ,γ^ν] = so(3,1)`. Verified: 6
   independent, closes under bracket, **Lorentzian** (boosts `σ^{0i}` non-compact /
   real eigenvalues, rotations `σ^{ij}` compact / imaginary). Identified with the
   v59 **6-field Cosserat**: `σ^{ij}=3φ`, `σ^{0i}=3θ` — the soldering 2-form's home.
2. **Cl(7)_even** — the internal `Spin(7)` (21 generators) on the octonion 8
   (same rep as v59; commutant in `End(ℝ⁸)=1` by Schur, from `06`).
3. **Commutation** — `max |[σ^{μν}⊗I , I⊗Spin(7)]| = 0` exactly.
4. **Internal untouched** — `max |[σ^{μν}⊗I , I⊗(any internal op)]| = 0`: color,
   `G₂`, lepton=L are internal operators, all spacetime-blind.
5. **Triality** — a generation (Z₃-type) rotation is an internal-only op; commutes
   with all Lorentz. The 3 generations are spacetime-blind.
6. **Prefactor** — `21/16 = dim Spin(7)/dim Cl(3,1)` is the factorization's
   internal/spacetime dimension ratio (matches v59 SYNTHESIS §9).

## What was proved (`G9Unification.lean`, axiom-clean)

- `spacetime_internal_commute` — **for ALL matrices `A, B`**: `(A⊗I)(I⊗B) =
  (I⊗B)(A⊗I)`. This is the general reason spacetime ⊥ internal, not just the
  numerical instance. (Via `mul_kronecker_mul`.)
- `spacetime_internal_bracket_zero` — the commutator is `0`.
- `dims_eq`, `module_dim`, `unified_algebra_dim`, `prefactor_is_dim_ratio` —
  `dim Cl(3,1)=16`, `dim Spin(7)=21`, module `4·8=32`, algebra `16·64=1024`,
  `21/16 = dim Spin(7)/dim Cl(3,1)`.

## Status of C3

**Closed** (derived + constructed + machine-checked), to the level v59's framework
supports:

- **Derived** (`06`): Schur ⇒ no `so(3,1)` fits inside `Cl(7)_even`; the carrier is
  forced into the spacetime factor. (`G9Soldering.no_internal_lorentz`.)
- **Resolved fork** (v59 docs): the spacetime factor is a separate `Cl(3,1)`, the
  `ℍ` inside the internal algebra is `SU(2)_L`.
- **Constructed + verified** (`07`): explicit `Cl(3,1)⊗Cl(7)_even`; spacetime
  Lorentz = Cosserat bivectors; commutes with the full internal symmetry, internal
  physics and triality untouched.
- **Proved in general** (`G9Unification.lean`): tensor-factor commutation for all
  operators; axiom-clean.

**Residual (now a minor refinement, not a blocker):** whether the unification is
the *ungraded* tensor product `Cl(3,1)⊗Cl(7)_even` (dim 1024) or the `Z₂`-graded
(super) tensor product. The even/Lie-algebra commutation proved here holds either
way; choosing the grading is a model decision for the full Lagrangian (`L`-track),
not a gap in the gravity-sector compatibility.

## Net effect on the G9 picture

The soldered-tetrad route now has: ±2 iff soldered (`04`), exactly 2 TT DOF (`05`),
v59 scalar law as the trace sector (`05`), the carrier *forced* into the spacetime
`Cl(3,1)` factor (`06`), and that factor *explicitly* commuting with the entire
internal sector in v59's own factorization (`07`). The one genuinely open piece is
now the **full Plebański action + EOM** (write `S[B]` for the `Cl(3,1)` 2-form with
simplicity constraint, sourced by `ρ_grav`, and show its weak-field limit
reproduces `05`).
