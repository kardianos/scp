# Gravity gaps — solution space

Opening the solution space for **G8** (magnitude mechanism for `α^(21/2)`) and
**G9** (a genuine spin-2 mode). Each candidate has a **test** and a **falsifier**.
Numerics: `g8_exponent_test.py`, `g9_polarization_test.py`.

---

## G8 — mechanisms for `f_g = α^(21/2)` / the `21`-generator product

The robust fact (data-forced): the exponent is **exactly 21 = dim Spin(7)**. The
question is *what coherent object over the 21 generators produces a power 21*, and
whether it is forced (mechanism) or merely consistent (numerology).

### G8-A. Determinant / top-form over the 21 generators **(LIVE, cleanest)**
The determinant of the 21×21 scalar coupling matrix `g = √α · I₂₁` is `(√α)²¹ =
α^(21/2)` — **exactly** `f_g` (verified `det_scalar_matrix_pow`, B2). Equivalently
the **top exterior form** `e₁ ∧ … ∧ e₂₁` of the Spin(7) Lie algebra (a 1-dim Λ²¹),
with each leg `√α`, has "volume" `α^(21/2)` (B3). This is the natural object that
turns "all 21 generators at once" into a power 21 — the **multiplicative** dual of
the **additive** gauge index (`g_W² = 5√α`, a *sum* over 5 generators).
- **Test**: does any coupling in a candidate Lagrangian appear as `det(coupling-
  matrix)` or as a Λ²¹-top-form/Haar-volume term on the 21-dim group manifold?
  A Wess-Zumino-like `∫_{Spin(7)} Tr(g⁻¹dg)^{21}` (Haar volume) is the prototype.
- **Falsifier**: if the only natural group-invariant couplings are **traces**
  `Tr(F^k)` (additive, even-degree characteristic classes), a determinant/odd
  21-form never appears — then α^(21/2) is not generated and G8-A is dead.
  Note dim Spin(7)=21 is **odd**, so `Tr(F^k)` (even forms) cannot directly give a
  21-form; only the volume/Haar top-form (odd-dim group) can — this is the crux to
  check.

### G8-B. Instanton / `e^{-S}` with S = 21·ln(1/α) **(LIVE, reframing)**
`α²¹ = exp(21 ln α) = exp(−S)` with `S = 21·ln(1/α) = 103.3` (D2). The *measured*
`−ln(α_G(e)) = 103.06`; the **0.27 difference is exactly `ln(prefactor) = ln 1.31`**.
So the entire prefactor question is *the sub-leading correction to a 21-instanton
action* — "21 instanton-units, one per Spin(7) generator." A single SU(2)_L
instanton `S = 8π²/g_W² ≈ 185` is ~1.8× too large; you need exactly the per-
generator (21-fold) version.
- **Test**: is there a configuration whose Euclidean action is `21·ln(1/α)` — e.g.
  a fractional/constrained instanton summed over the 21 generator directions, or a
  21-fold covering? The clean signature: the prefactor should emerge as the
  one-loop determinant around the 21-instanton, fixing 21/16 (or 17/13) non-freely.
- **Falsifier**: standard instanton actions are `8π²/g²` (∝ 1/α, **not** ∝ ln(1/α)).
  A `ln(1/α)` "action" is anomalous — it is really an `α²¹` *power*, i.e. a 21-loop
  / 21-leg perturbative object, not a semiclassical instanton. If no 21-leg diagram
  or 21-fold topological term is identifiable, this is numerology dressed as e^{-S}.

### G8-C. Dimensional transmutation **(DEAD as a derivation)**
`(m/M_Pl)² ~ exp(−2π/(bα))` requires `b = 8.35` (D1) — not a SM-natural integer,
and it produces a *scale* not a *ratio to M_Pl* without inputting M_Pl. It relabels
the hierarchy, does not derive it. **Falsifier already triggered** (non-integer b,
needs M_Pl input).

### G8-D. The prefactor as a one-loop determinant **(open, the real prize)**
The exponent is law-like; the prefactor (21/16 vs 17/13 vs 4/3) is the weak link.
Any genuine mechanism (G8-A or G8-B) must *also* fix the O(1) prefactor as the
fluctuation determinant / normalization of the 21-object. Until that pins 21/16
uniquely, G8 stays a **value conjecture with a striking exponent**.
- **Test**: compute the Gaussian fluctuation determinant around the
  top-form/instanton; does it equal 21/16 = dim Spin(7)/dim Cl(3,1)?
- **Falsifier**: if it equals 17/13 or an irrational, the "21/16" reading is wrong
  (and the exponent stands alone).

---

## G9 — routes to a spin-2 (`h = ±2`) mode

The internal `Λ² = so(8)` bivector does **not** carry spacetime helicity (proved:
`internal_index_inert`). A spacetime *antisymmetric* 2-form is spin-1
(`spacetime_bivector_is_spin1`). The **only** `h=±2` carrier is a **symmetric**
spacetime rank-2 tensor. So every viable route must produce a *symmetric spacetime
tensor* `h_μν`.

### G9-A. Induced / emergent metric from the bivector connection **(most promising)**
If the bivector connection `Ω_μ^{ab}` (a spacetime 1-form valued in so(8)) defines
a **vielbein** `e_μ^a` (so(8) acting on an 8-dim "internal frame" that is identified
with an 8-dim spacetime-extension, or a Plebański-style construction where the
2-form *is* the gravitational field), then `g_μν = e_μ^a e_ν^b δ_ab` is **symmetric**
and its TT fluctuation is `h = ±2`. This is the bilinear `T_μν = Ω_μ^a Ω_ν^a`
promoted from "source" to "the metric itself" (Plebański: GR *is* a constrained
BF/2-form theory; the 2-form is fundamental, the metric is derived).
- **Test**: can the OBE be rewritten as a constrained `BF`/Plebański theory where
  `Ω ∈ Λ²` is the self-dual 2-form `B`, and the metric/`h_μν` emerges with the
  correct 2 propagating TT DOF? Count: a massless spin-2 has 2 DOF; check the
  constrained 2-form delivers exactly 2 (not 0, 1, or 5).
- **Falsifier**: if the so(8) 2-form has no constraint reducing it to a metric (no
  simplicity/Plebański constraint available in the algebra), or if the DOF count is
  wrong, no symmetric tensor emerges and gravity stays scalar. The decisive check:
  does the bilinear `Ω_μ^a Ω_ν^a` propagate a **massless pole** (long-range) or only
  appear at `O(α)` short-range? Current analysis says the latter ⇒ falsified unless
  a fundamental (non-composite) symmetric mode is posited.

### G9-B. Bivector ⊗ bivector → symmetric tensor in the so(8)/Λ² product
Decompose `Λ² ⊗ Λ²` of so(8): it contains the **symmetric traceless** rep
(the "Weyl/Riemann-like" piece) alongside the adjoint and singlet. A spin-2 field
could be the component of `Ω ⊗ Ω` (or of a `Λ⁴` curvature `R = dΩ + Ω∧Ω`) living in
the symmetric-traceless sub-rep — *if* an internal-to-spacetime index identification
(a soldering/vielbein) maps it to a symmetric *spacetime* tensor.
- **Test**: in `so(8)`, `Sym²(adjoint) ⊃` a rep that, under a Spin(7)→ Lorentz
  soldering (the G₂⊂Spin(7)⊂Spin(8) chain restricted to a 4D subspace), maps to the
  symmetric-traceless `(2,0)+(0,2)` Lorentz rep (the Weyl tensor / graviton). Test
  the branching numerically.
- **Falsifier**: if the branching of `Sym²(so(8)-adjoint)` to the embedded Lorentz
  group contains **no** symmetric-traceless rank-2 spacetime piece, there is no
  graviton in the bivector sector. (The historical V6/null-rotor result — only
  `h=0,±1`, **no `h=±2`** — is the cautionary precedent: a different algebra, but
  the same failure mode.)

### G9-C. Add a fundamental spacetime graviton (honest fallback)
Accept that gravity's spin-2 sector is **not** inside the internal Furey algebra and
add a standard symmetric `h_μν` (linearized GR) sourced by `ρ_grav`, with coupling
`f_g = α^(21/2)` from G8. v59 then *predicts the magnitude* but does **not** unify
the graviton with the gauge bivector — gravity becomes an external spin-2 field with
a v59-fixed strength.
- **Test**: does sourcing a standard `h_μν` by `T_μν` built from `ρ_grav` reproduce
  EP + `1/r²` + `h=±2`? (Yes by construction — this is GR with a fixed `G_N`.)
- **Falsifier**: none physically, but it **concedes the unification claim** for
  gravity: the spin-2 mode is not emergent from the algebra. This is the
  fall-back, not a success.

### G9-D. Tensor from a `Λ⁴` (F-grade) curvature **(speculative)**
The F-grade `Λ⁴` (35-dim, the G₂ associative-form sector) squares to `+1` (real
involutions). A *symmetric* object naturally lives in the real (F) grade, not the
complex (L) grade. Perhaps the spin-2 graviton couples to the **F** current
`𝒥_F = ⟨Φ Φ̃⟩_F` (symmetric) rather than the **L** current (antisymmetric/spin-1).
- **Test**: is the F-grade current `𝒥_F` symmetric under the relevant index swap,
  and does its long-range mode carry `h=±2` after soldering?
- **Falsifier**: `𝒥_F` is a color/strong current (short-range, confined) in v59, so
  it has no massless long-range mode — likely falsified, but worth the branching
  check because "real grade ↔ symmetric tensor" is the right parity.

---

## Summary of the solution space

| gap | most promising | cleanest test | key falsifier |
|---|---|---|---|
| **G8** | top-form/Haar volume over 21 gens (A) + its fluctuation determinant fixing the prefactor (D) | does a Λ²¹ Haar-volume term appear, and does its one-loop det = 21/16? | only `Tr(F^k)` (additive, even-degree) couplings exist ⇒ no 21-form ⇒ no α^(21/2) |
| **G9** | induced/Plebański metric from `Ω∈Λ²` (A), checked via `Sym²(so(8))` branching (B) | does constrained-2-form/`Sym²(adjoint)→Lorentz` give a symmetric-traceless rank-2 with exactly 2 TT DOF? | no symmetric-traceless rank-2 in the branching, or only an `O(α)` short-range composite ⇒ no graviton |

**The decisive ordering**: G9 dominates. A correct magnitude on a scalar force is
still ruled out by LIGO. If G9-A/B fail the branching/DOF test, v59 gravity is
falsified as a unification (G9-C concedes it), regardless of G8.
