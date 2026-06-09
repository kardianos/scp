# G9 Geometric Origin — The Octo-Space 3-Form / 4-Form and the Fundamental OBE

> **CORRECTION (2026-05-25, same day): see `04_findings.md`.** This document
> proposes soldering an algebraic 2-form **valued in the internal `Λ²(V^8)=so(8)`**
> via `G₂` 3-form/4-form constraints. The machine-checked helicity no-go
> (`../lean/G9Soldering.lean`) shows that **no constraint on a purely internal 2-form
> can produce spacetime helicity ±2** — helicity is a spacetime little-group label.
> The fix is not a cleverer internal constraint; it is to put the 2-form in the
> **spacetime `Cl(3,1)` bivectors** (`so(3,1)`, already present as the 6-field
> Cosserat sector) and keep `ρ_grav` as the scalar source. The `G₂`/4-form ideas
> below may still matter for the *simplicity constraint's coefficients*, but the
> carrier's index placement is the decisive correction.

**Date**: 2026-05-25 (continuation of v60 G9 attack)
**Parents**:
- `v59/gaps/gravity/ALTERNATIVES.md` (G9-A Plebański route)
- `v59/synthesis/NEW_OBE_FORMULATION.md` (the current fundamental multivector equation)
- `v59/gaps/quark_flavour/README.md` (explicit statement of the coassociative 4-form *φ)
- `v59/gaps/gravity/FINDINGS.md` and `g9_spin2_route_test.py` (the "no soldering" diagnosis)
- v60 `01_8space_to_spacetime_bindings.py` + `02_constrained_helicity_count.py` + Lean `G9ToyHelicity.lean` (the current toy model of the constraint)

---

## 1. The Diagnosis from v59 (the missing geometric object)

v59 established (theorem-grade in Lean + numerical helicity decomposition):

- The propagating gravity mode in the OBE is carried by a connection / 2-form valued in the **L-grade** Λ² ⊕ Λ⁶ ≅ so(8) acting on the **8-dimensional octo-space** V^8 (the vector representation underlying Spin(8), triality 8_v / 8_+ / 8_-, Cl(7)_even).
- This gives a perfectly good **scalar** long-range force sourced by the second-moment charge ρ_grav = Tr(M†M) = 9Qa² (EP-exact, same Frobenius² structure as the 784 EW bridge).
- The algebra contains spin-2 representations (Sym²(so(8) adjoint) contains symmetric-traceless rank-2 pieces on V^8; also symmetric objects naturally live in the real F-grade).
- **The fatal gap (G9)**: there is currently **no soldering / simplicity constraint** that maps those internal spin-2 pieces (or the 2-form itself) to a **spacetime** symmetric traceless tensor h_μν carrying exactly two transverse-traceless degrees of freedom with helicities {+2, -2}.

The documents repeatedly say the same thing:

> "v59 has the G₂ ⊂ Spin(7) ⊂ Spin(8) chain, but **NO soldering form identifying 4 of the 8 internal directions with spacetime**."  
> (g9_spin2_route_test.py, FINDINGS.md, ALTERNATIVES.md)

The **geometric origin** we are looking for is precisely the object (or pair of objects) inside the octo-space structure that can serve as that soldering data.

---

## 2. The Octo-Space and Its Built-in Geometric Forms

The parent algebra is Cl(7)_even ≅ ℂ ⊗ 𝕆 (dim 64 real), acting on structures built from the **real 8-dimensional space V^8** (vector rep of Spin(8), with triality permuting three 8-dimensional irreps).

Its Clifford-grade decomposition on ℝ^7 (or the associated 8D picture) is:

```
Λ⁰ = 1
Λ² = 21 (= dim Spin(7))          ← so(7) content, complex structures J (L-grade)
Λ⁴ = 35                           ← contains the coassociative 4-form *φ
Λ⁶ = 7 (= dim Im 𝕆 = dim S^7)
```

**The L ⊕ F bisection** (involution μ on grades, already theorem-grade):

- **L = Λ² ⊕ Λ⁶** (μ = −1, dim 28 = dim Spin(8))  
  "Lie-algebra content".  
  Carries the **complex structures** J ∈ Λ² with J² = −1 (used to force lepton = L, pin the Brannen phase via the sedenion S₃ automorphism, and define the skew part of the connection).

- **F = Λ⁴** (μ = +1, dim 35)  
  "**G₂-form content**".  
  Contains the **coassociative 4-form *φ** (the unique even-grade G₂ singlet).  
  The associative 3-form φ (G₂ structure on the 7D imaginary octonions) lives in the odd-grade picture and calibrates associative 3-planes; its dual *φ calibrates coassociative 4-planes.

These are **exactly** the calibrated forms that appear in G₂-geometry and in 7D/8D generalizations of Plebański / BF / 2-form gravity.

The same objects that already do the heavy lifting for:
- Lepton forcing (B² = −1 in L-grade)
- Color splitting (F-grade supplies the Cartan)
- Triality / three generations (Z₃ action on the three 8's of the octo-space)
- The structural origin of Q = 14/21 and φ = 2/9

are the natural candidates for the **simplicity / soldering constraint** that turns the algebraic 2-form B (valued in Λ² or the full algebra) into a spacetime metric + connection with the correct tensor polarizations.

---

## 3. How the Forms Provide the Soldering (Geometric Origin)

In standard Plebański theory (and its G₂ / octonionic extensions), a 2-form B (possibly with internal indices) is subject to a constraint of the schematic form:

```
B ∧ B  ~  vol₄ ⊗ (structure form on the internal space)
```

or

```
B^{ab} ∧ B^{cd} ⋅ φ_{abcd...}   =   metric volume form
```

Concretely, using the v59 objects:

1. **Self-duality from L-grade J**  
   The complex structures J ∈ Λ² (already central to the Brannen phase and lepton = L) can be used to define a notion of self-dual / anti-self-dual 2-forms on the internal indices, exactly as in 4D Plebański gravity (where one uses the spacetime Hodge * to split into self-dual 2-forms).

2. **Metric extraction / soldering from F-grade *φ**  
   The coassociative 4-form *φ (G₂ singlet in Λ⁴) provides a natural 4-form on the 8-space (or its 7D reduction). Contracting two algebraic 2-forms B and B with *φ (or the associative 3-form φ) yields a symmetric bilinear form on spacetime indices:

   ```
   g_μν ~ (B_μ^{ab} B_ν^{cd} *φ_{abcd})   or   ε_{μνρσ} B^{ρ ab} B^{σ cd} φ_{abcd}
   ```

   This is the direct analog of the Plebański constraint that forces the 2-form to determine a metric (or tetrad) while automatically enforcing the correct number of degrees of freedom (2 TT for a massless graviton in 4D).

3. **Compatibility with existing structure**  
   - The L/F bisection already separates "complex / Lie-algebra" (L, good for self-duality and so(8) connection) from "real G₂-form" (F, good for the calibrating 4-form).
   - Triality (Z₃ on the three 8's) acts on both the mass kernels (Brannen) and on the geometric forms; any constraint must be compatible with it or explain why gravity sees a different combination.
   - Color su(3) (from the octonion factor, acting on the 8) must remain unbroken or be broken in a controlled way — the constraint should not disturb the color splitting already proved in the quark sector.

This is why the "octo-space" is not an arbitrary internal manifold — its built-in G₂-calibrated geometry (3-form + 4-form + complex structures on the grades) supplies the soldering data "for free" once we decide to treat the algebraic 2-form as fundamental.

---

## 4. Tying Back to the Fundamental Physics Equations (the OBE)

The current v59 fundamental equation (NEW_OBE_FORMULATION.md) is written in integrated form:

```
Ω(x) = ∫ K(x,x') [ f_g ∇' (∑ ρ_N) + f_W J_L + f_C J_F ] d³x'
```

with the gravity piece reducing (for massless K and nonzero monopole) to the explicit wave equation

```
□ Ω_grav = f_g ρ_grav     (ρ_grav = Tr(M†M) = second moment on the generation space)
```

Here Ω_grav lives in the L-grade Λ² (so(8) on the 8-space).

**Geometric reinterpretation once the binding exists:**

Treat the fundamental object as a **spacetime 2-form B_μν** whose values are in the algebraic grades (primarily L-grade for the connection part, with F-grade forms entering the constraint):

- The **simplicity constraint** (using J from L + *φ from F) is imposed algebraically or dynamically. Solving it defines both the metric g_μν (or tetrad) **and** projects the algebraic connection onto a spacetime spin connection + the correct internal so(8) part.

- The **dynamics** become a 2-form field equation (Plebański / BF-type action possibly with cubic terms coming from the octonion associator or the original multivector quadratic terms):

  ```
  dB + B ∧ B  +  constraint terms involving φ and *φ   =   source current built from ρ_grav
  ```

- The original integrated OBE expression for Ω_grav is recovered as the **on-shell / integrated solution** for the connection part of B after the constraint is solved. The source ρ_grav (second moment of the Brannen kernel) appears naturally because the mass kernels themselves live on triality-selected 3-planes inside the same 8-space whose G₂ structure is providing the soldering.

- The quadratic and higher terms in older multivector formulations (Ω², Ω³, associator terms) arise from the non-linearities of the 2-form curvature F = dB + B∧B once the constraint is solved for the metric and the internal components.

- The L/F currents J_L and J_F become the projections of the full algebraic current onto the grades that couple to the respective fermion sectors (exactly as in the current OBE).

- The scale bridges (v_Higgs = 28² a_ℓ², G_e ~ α^{21/2}, etc.) remain as before; they are statements about dimensions and Frobenius norms on the 8-space and its grades, now interpreted geometrically as volumes or instanton actions calibrated by the same G₂ forms.

In short: the current OBE is the **effective, integrated, scalar-sourced** description. The geometric origin (G₂ 3-form + coassociative 4-form + L-grade J's on the octo-space) supplies the missing **local, constrained, tensorial** formulation whose EL equations + simplicity constraint reproduce the old equation while automatically giving LIGO-compatible gravity.

---

## 5. Concrete Next Steps (toward v60 deliverable)

1. **Formalize the constraint in the toy model** (already begun in 02_*.py + G9ToyHelicity.lean). Replace the schematic projector with an explicit contraction involving a model 3-form / 4-form on the 8-space (even if the numbers are still toy).

2. **Write the geometric action** in v60/gravity_recast/ (Plebański-style S[B] = ∫ B ∧ B ∧ φ + ... + source term from the Brannen ρ on the triality 3-space). Derive the EL equation and show it reduces to the current □Ω_grav = f_g ρ_grav when the constraint is solved.

3. **Lean formalization**: Extend G9ToyHelicity.lean (or a new G9G2Constraint.lean) with statements such as:
   - "The coassociative 4-form *φ (G₂ singlet in Λ⁴) defines a simplicity projector on 2-forms valued in Λ²."
   - "The resulting physical helicity content is hel_sym2TT."
   - Prove invariance under the structures already theorem-grade in v59 (triality, color splitting, J pinning lepton = L).

4. **Tie to the full OBE**: Show that the L-grade and F-grade currents in the current integrated equation arise as the projections of the full algebraic current once the 2-form B is decomposed according to the L/F bisection and the constraint is solved.

5. **Falsification / alternative**: If no natural constraint using φ / *φ / J produces exactly 2 TT DOF while preserving the existing bridges and selection rule, then G9 is cleanly falsified inside the octonionic framework (and we fall back to G9-C or a program-level restructuring).

---

This document closes the loop: the geometric origin of the required soldering is not an external addition — it is the **G₂-calibrated 3-form and 4-form** (plus the already-used complex structures) living in the graded octo-space that already explains the lepton forcing, Koide ratio, color structure, and triality. Once we promote the algebraic 2-form to the fundamental carrier and impose the constraint built from those same forms, the current OBE becomes the effective description of a geometrically consistent tensor theory of gravity sourced by the Brannen second-moment densities on the 8-space.

The v60 task is now to make this explicit (action + constraint + Lean + numerical test) rather than schematic.

Next artifact suggestion: a Python sketch of the constraint using a model G₂ 3-form on an 8D space + the resulting metric extraction, followed by the helicity count (extending 02_*.py). Then the corresponding Lean statements.