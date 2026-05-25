# Forward Proposal: the v58 field equation as the dynamical home of the v59 constraints

*Brings the v58 "unified multivector force law" and its constraints into v59 as a
forward proposal, and uses it to sharpen the current `φ = Q/3` search.*

## 1. The v58 base equation and its constraints

The v58 pre-geometric program (`v58/pregeometric/unified_multivector_force/`) locked a
single retarded multivector equation in **Cl(3,0)**:

>   ⟨ D Ω + λ Ω² + μ ⟨Ω,Ω⟩ ⟩_{0,2}  =  f_g(ρ) · J_ρ  +  f_em(ρ) · J_χ

with the *constraints* its Lean track identified as necessary for the classical limits:

- **(C-surface)** a density `ρ_M = ½(M M̃ − v²)` measuring deviation of the multivector
  norm from a **Higgs-like vacuum manifold** `M M̃ = v²`;
- **(C-grade)** the equation lives on the **even sector** `⟨·⟩_{0,2}` = grade-0 (scalar)
  ⊕ grade-2 (bivector) — which in Cl(3,0) *is the quaternions* `Cl⁺(3,0) ≅ ℍ`;
- **(C-split)** **scalar grade = density/gravity (real)**, **bivector grade = protected
  chirality/EM (complex)**, with **no cross-talk** between the channels;
- **(C-quad)** the quadratic self-interaction `λΩ²+μ⟨Ω,Ω⟩` must be orthogonal to / higher
  order in both channels at linear order (safe band `|λ|≤0.005, |μ|≤0.001`);
- **(C-conn)** density gradients induce a **bivector connection** `ω ∝ ∇(M M̃)/|M|²`
  that governs how every excitation (including the chiral/light modes) propagates.

## 2. The constraint correspondence (v58 ↔ v59)

These are not the same algebra (v58 = Cl(3,0); v59 = Cl(7)_even), so the bridge is at the
level of the **even/ℍ sector and the constraint structure**, made rigorous by one
signature-independent fact: **a grade-2 bivector squares to −1 in *any* signature**
(`(eᵢeⱼ)² = −eᵢ²eⱼ² = −1` whether `e²=+1` or `−1`) — so the grade↔(complex/real) split is
universal (`lean/BladeSquareSign.lean`).

| v58 constraint | v59 constraint | bridge |
|---|---|---|
| vacuum manifold `M M̃ = v²` (Higgs-like) | Koide surface `\|ξ\|² = 1/2` + ξ-Higgs potential | both = fixed-norm vacuum manifold (`lean/XiVacuum`, `SilentDirection`) |
| even sector `⟨⟩_{0,2}` = ℍ | Brannen ℍ-slice `(e₀₁,e₀₂,e₁₂)`, silent-direction ℍ | `Cl⁺(3,0) ≅ ℍ` = the lepton quaternion slice |
| scalar = density/gravity (real) | symmetric `F`+scalar = real / Koide amplitude `Re ξ` | grade-0 squares to `+1`: real channel |
| bivector = protected chirality/EM (complex) | `J_c ∈ Λ²` = complex structure / `Im ξ` / phase | grade-2 squares to `−1`: complex/chiral channel |
| **no cross-talk** (scalar ⊥ bivector) | **L(skew)/F(symmetric) dichotomy**, F-subspace has no complex structure | `lean/LeptonRealityForcing`, `BladeSquareSign` **derive** the split v58 imposed |
| connection `ω ∝ ∇(MM̃)/\|M\|²` | (no analog yet) | **the missing piece — see §3** |

**Net (rigorous part):** v59's algebraic-rigidity theorems *derive* the grade-separation
(C-split) that v58 had to *assume*: in any signature the bivector grade is forced to be the
complex/chiral carrier and the scalar grade the real/density carrier, and the symmetric
sector carries no complex structure (no cross-talk).  And the v58 Higgs-like vacuum manifold
is the v59 Koide surface, with v59's proven ξ-Higgs potential (`XiVacuum`) supplying its
potential.

## 3. The forward proposal (what each supplies the other)

- **v58 → v59 (the dynamics):** v59 is "Kepler" — exact algebraic patterns (Koide `Q=2/3`,
  phase `φ=Q/3`, the grade split) with no equation of motion.  v58 is "Newton" — a retarded
  field equation with forces and a **connection Ω**.  The v59 order parameter `ξ` is the v58
  multivector `M` restricted to the even/ℍ sector on the vacuum manifold; the v59 constraints
  are the constraint structure of the v58 dynamics there.

- **v59 → v58 (the missing angular constraint):** v58's modulation `f_g(ρ)` depends only on
  the **magnitude** (density); v58 has **no phase-selecting mechanism** on its vacuum manifold.
  v59 supplies exactly that: the Koide amplitude (`\|ξ\|²=1/2`) fixes the radial position and
  the phase law (`φ=Q/3`) fixes the **angular** position of the order parameter on the ℍ
  vacuum manifold.  So `φ=Q/3` is the angular boundary condition v58's vacuum was missing.

## 4. Using this in the current search (the `φ=Q/3` mechanism)

This **concretizes the leading mechanism candidate M2/M3** (geometric/Berry phase of the Z₃
generation cycle) by giving it the connection it needs:

> **Reformulated mechanism:** the phase `φ` is the **holonomy of the v58 bivector connection
> `ω ∝ ∇(M M̃)/|M|²` transported around the Z₃ generation cycle** on the Koide vacuum manifold
> `|ξ|²=1/2`.  The law `3φ = Q` is the conjecture that this closed-cycle holonomy equals the
> Koide ratio `Q = dimG₂/dimSpin7`.

Why this is the right move:
- v58's `ω` is *already a bivector (grade-2 = complex/chiral) connection* — exactly the object
  whose holonomy is a phase.  M2 needed "a connection whose Z₃ holonomy is `Q`"; v58 provides
  the connection from first dynamical principles, not by hand.
- The "/3 = #generations" is automatic for a phase accumulated over the 3-step Z₃ cycle.
- The lepton-specificity follows: the holonomy lives on the lepton color-singlet ℂ-line
  `{0,7}` (where `J_c` defines the complex phase); the quark triplet carries a different
  (non-abelian) holonomy → quarks need not satisfy `φ=Q/3` (matches `study.py`).

### What the Witt map + holonomy test actually returned (2026-05-24)

The Witt map was constructed and verified (`witt_map.py`): three **disjoint γ-pairs**
`(γ₀,γ₅),(γ₁,γ₄),(γ₂,γ₃)` with the *scalar* `ℂ⊗𝕆` unit `i` give three fermionic modes
(nilpotent, CAR `{aᵢ,aᵢ†}=1`, cross-anticommutators zero), and the number operator
`N = Σ aᵢ†aᵢ` has eigenvalues exactly `(0,1,1,1,2,2,2,3)` — the Fock grading.  So the
deferred octonion→Fock map is done: **8 Fock states = `Λ•(ℂ³_color)` of ONE generation,
built on the three COLORS.**

Running the holonomy test on this structure gives a **decisive negative** (machine-checked,
`lean/ColorSU3.lean`: `colorZ3_cubed`, `colorZ3_commutes_Jc`, `colorZ3_fixes_lepton`):

> the color `Z₃` `R` (cycling the three quark modes, `R³=I`, `[R,J_c]=0`) **fixes the lepton
> singlet line `{0,7}` pointwise**.  So a color-`Z₃` holonomy on the lepton is identically
> trivial — **`φ = 2/9` is NOT an intra-ideal color holonomy.**

**Consequence (the located missing ingredient):** the Brannen phase lives on the *generation*
`Z₃` (Wilson loop `arg(ξ³)=3φ`), which is a **3-copy structure external to one ideal**.  Every
`Z₃` internal to the single ideal either fixes the lepton singlet (the color `Z₃`) or fails to
produce three equal generation eigenspaces (the `J_c`-phase `Z₃`).  So the holonomy mechanism
**cannot** be read off the single-ideal Witt/color structure; it requires first **building the
inter-generational connection** — three copies of the ideal linked by a `Z₃` connection.

### Concrete next experiment (revised)
This is exactly where the **Cl(3,0)→Cl(7) / v58-multivector expansion** the user flagged must
do its work: supply the 3-copy generation structure and its connection (the v58 `ω` lifted to
the generation index), then test whether *that* closed-cycle holonomy is forced to
`Q = 14/21` (⇒ `φ=Q/3`).  Steps:
1. **Generation extension:** build `(3 generations) ⊗ (one ideal)` with an inter-generation
   connection (candidate: the v58 bivector connection `ω ∝ ∇(MM̃)/|M|²` lifted along the
   generation index, or the Brannen coupling `ξ` itself as the connection).
2. **Holonomy:** compute `arg(ξ³)` of the closed generation cycle on the `|ξ|²=1/2` manifold.
3. **Test `= Q`.**  Positive ⇒ first dynamical derivation of `φ=Q/3`; negative ⇒ rigorous
   exclusion of the geometric-phase mechanism, leaving `φ=Q/3` a sharp unexplained regularity.

**Net:** the Witt map is done and the *intra-ideal* holonomy route is now closed (negative,
machine-checked).  The live route is the inter-generational connection — a smaller, sharper
target than before.

## 5. Honest caveats
- The full algebras differ (Cl(3,0) vs Cl(7)); the bridge is the **even/ℍ sector + the
  signature-independent grade↔complex law**, not an identification of the whole algebras.
  (Cf. the 2026-05-24 retraction of the *octonionic* living-candidate crossover, which
  conflated the dynamics; this proposal connects only the **constraints** and the **connection**,
  which the universal grade law makes legitimate.)
- v58's connection `ω` is defined on its 2D retarded lattice; transporting it on the *generation*
  (internal Z₃) space is a hypothesis to be tested, not yet a theorem.
- `φ = Q/3` remains empirically precise but mechanism-open (`CONCEPT.md`); this proposal makes
  the mechanism search concrete, it does not yet close it.
