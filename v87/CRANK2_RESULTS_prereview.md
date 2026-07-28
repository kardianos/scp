# v87 B0-G, crank 2 — why p = 1, and how it fits the fabric
**Date:** 2026-07-26 · Follow-on to the three-seat rung (`FABLE_WORK.md`,
`GROK_WORK.md`, `GEOM_WORK.md`). Worksheets: `work/geom/nd_shape_limit.py`,
`work/geom/why_p1.py`, both run clean end to end.

Tags: **[D]** derived · **[M]** measured · **[P]** postulated · **[C]** conjecture

---

## 0. Where crank 1 left it

The three seats established: the ontology relaxes **(MI)**; explicit local
deterministic models reach exactly 2√2 at a cost of **0.046 bits**; and the
Tsirelson half was **not** explained, because measurement dependence gives a
smooth dial from 2 through 2√2 to 4 with nothing distinguishing the quantum
point.

Crank 1 then showed the ceiling is a **shape limit**, not a quantum phenomenon:

* E(θ) is **exactly independent of the dimension N** of the body λ lives on —
  the radial integral cancels identically between numerator and denominator.
  Verified symbolically (N = 3, 4, 5, 10) and by Monte Carlo (N = 2…25, error
  ≤ 4×10⁻⁴). **[D,M]** *Different perspectives give the same ratio.*
* The waveform ladder is real: triangle → **2**, cosine → **2√2**, square → **4**.
* The √2 is the ℓ¹/ℓ² ratio √d with **d = the number of terms in the
  inequality**, not the dimension of the body.
* Chained Bell: quantum/local ratio is **maximal at n = 2** (√2) and falls
  monotonically to 1 — CHSH is the extremal inscribed shape; refining the
  polygon toward the circle makes the world look *more* classical.

This left one question: the p-family ρ_p(λ|b) ∝ |λ·b|^p is a dial, and geometry
does not say which rung. **Why p = 1?**

---

## 1. Closed form for the whole ladder  **[D]**

Using |cos φ|^p sgn(cos φ) = |cos φ|^{p−1} cos φ ≡ g(φ), and the fact that g is
odd about φ = π/2 (so its full-circle integral vanishes, and the sgn factor
merely doubles the half-circle integral), differentiation under the integral
gives the endpoint terms directly:

```
    dE_p/dθ = (4/Z_p) |sin θ|^{p−1} sin θ

    E_p(θ)  = −1 + (4/Z_p) ∫₀^θ sin^p u du        (0 ≤ θ ≤ π)

    Z_p     = 2√π Γ((p+1)/2) / Γ(p/2 + 1)
```

The entire waveform ladder is a regularised incomplete beta function. Symbolic
verification at the landmarks:

| p | Z_p | E_p(θ) | |
|---|---|---|---|
| 0 | 2π | (2θ − π)/π | **matches** the triangle |
| 1 | 4 | **−cos θ** | **matches** the cosine |
| 2 | π | (2θ − sin 2θ − π)/π | — |

---

## 2. The mechanism: |x|·sgn(x) = x  **[D]**

The weighted B-response is

```
    ρ_p(λ|b) · B(b,λ)  ∝  |λ·b|^p · (−sgn(λ·b))  =  −|λ·b|^{p−1} (λ·b)
```

At **p = 1, and only there**, this collapses to exactly −(λ·b). The dichotomic
±1 outcome, weighted by how often it is sampled, **reconstitutes the linear
projection of the hidden variable onto the setting.** For p ≠ 1 a residual
|λ·b|^{p−1} survives and distorts the waveform.

This is a rectifier relation, not a coincidence of exponents: the detector
reports a *sign*, the rate at which it reports is proportional to the
*magnitude*, and sign × magnitude is the field itself.

---

## 3. Uniqueness: p = 1 is the only monochromatic rung  **[D,M]**

E_p is even and π-antisymmetric, hence a cosine series in odd harmonics.
Measured coefficients:

| p | c₁ | c₃ | c₅ | c₇ | ‖higher‖/c₁ |
|---|---|---|---|---|---|
| 0 | 0.81057 | 0.09006 | 0.03242 | 0.01654 | 1.2×10⁻¹ |
| 0.5 | 0.92880 | 0.04423 | 0.01206 | 0.00517 | 5.0×10⁻² |
| **1** | **1.00000** | −0.0000 | −0.0000 | −0.0000 | **2.4×10⁻¹⁰** |
| 1.5 | 1.04724 | −0.03879 | −0.00537 | −0.00158 | 3.7×10⁻² |
| 2 | 1.08076 | −0.07205 | −0.00618 | −0.00147 | 6.7×10⁻² |
| 16 | 1.23636 | −0.32536 | 0.12085 | −0.04128 | 2.8×10⁻¹ |

So the sharp characterisation is:

> **p = 1 ⟺ the correlation is a single harmonic ⟺ E = −cos θ ⟺ S_max = 2√2.**

Every other rung carries higher harmonics. The quantum point is the
**monochromatic** member of the family.

**And this is where the fabric enters.** The objects in this model *are*
monochromatic: Φ_a = f(r)e^{iωt}, one clock. HC-1 measured that the ball carries
at most **one** bound internal mode over the whole branch and **none** in the
working region ω ≥ 1.36. A single rotating phase, projected on a setting, yields
a single harmonic in the correlation. **p = 1 is what a monochromatic object
gives.** [D for the geometry; the fabric link is [C] until §5.3 is closed.]

---

## 4. The decisive fork — and one branch is already dead  **[D,M]**

ρ_p(λ|b) ∝ |λ·b| admits two physically different readings that produce
**identical correlations** but have **opposite experimental status**. This
distinction is load-bearing and was not drawn in crank 1.

**Reading 1 — source-side measurement dependence.** The fabric configuration
actually *produced* is correlated with the settings through the common past
(the chooser is fabric, O2). Every pair is detected. Violates (MI). **Not
excluded by any experiment.**

**Reading 2 — detector-side sampling.** The source emits uniformly; the detector
*fires* with probability ∝ |λ·b| and misses the rest. **This is exactly the
detection loophole, and it is closed.**

Reading 2 requires a detection efficiency equal to the mean weight:

| quantity | value |
|---|---|
| ⟨\|cos φ\|⟩ on the phase circle | **0.63662** (= 2/π) |
| ⟨\|λ·b\|⟩ on the sphere S² | **0.50010** (= ½) |
| CHSH detection-loophole threshold 2(√2 − 1) | **0.82843** |

Both sit far below threshold, so Reading 2 predicts an efficiency at which
loophole-free experiments could never have violated CHSH — and they did.
**Reading 2 is refuted.**

> **Hard design constraint on the fabric:** the p = 1 weight must be realised as
> a bias in *what the source produces*, correlated with the settings through the
> common past. It must **not** be realised as a detector that preferentially
> registers well-aligned configurations. Same equation, opposite verdicts.

---

## 5. Fitting it to the fabric

### 5.1 What matches, and it is more than expected

| what p = 1 needs | what v86 measured |
|---|---|
| a continuous internal phase on every object | Φ_a = f(r)e^{iωt}; ω on a branch; local clock resolved per voxel **[M]** |
| monochromatic (single harmonic) | HC-1: ≤1 bound internal mode on the branch, **none** for ω ≥ 1.36 **[M]** |
| coupling linear in the phase overlap | A6 contact law: force ∝ **cos(Δφ)**e^{−κD} **[M]** |
| a sign carried by that overlap | A6: co-phase **fuses** (×2.7), anti-phase **repels** ⇒ sgn(cos Δφ) **[M]** |
| weight = magnitude of the overlap | the same law's magnitude \|cos Δφ\| **[M]** |
| source–setting correlation via common past | O2 monism; capacity required 0.046 bits **[D]** |

The fabric already exhibits, **as a measured interaction law**,

```
    force  ∝  cos(Δφ)  =  |cos Δφ| · sgn(cos Δφ)
```

which is precisely **weight × response at p = 1**. That is not a fitted analogy;
it is A6 read literally.

### 5.2 What does not match — three gaps, stated exactly

1. **A force is not a correlation.** The contact law is a force between two
   nearby objects, not a correlation of dichotomic outcomes at spacelike
   separation. Same algebraic form, different observable. **[D]**
2. **Exponential suppression kills direct locking.** With κ ≈ 0.5117:
   e^{−κD} = 7.7×10⁻² at D = 5, **2.2×10⁻³ at D = 12** (where v86 measured the
   Adler bound < 0.06), 7.7×10⁻¹² at D = 50. Direct phase-locking across a
   Bell-relevant separation is numerically zero. The correlation must be
   established in the **common past** and survive transport — and GROK's G14
   obstruction is that the kernel's mixing works against exactly that. **[M/D]**
3. **Monochromaticity is a property of the object, not a derivation of the
   sampling weight.** §3 says a monochromatic source *would* give p = 1; it does
   not yet show that the fabric's sampling *is* the projection weight. This is
   the same "capacity is easy, architecture is the wall" gap the seats found,
   now localised to one exponent. **[C]**

### 5.3 What would close it

Derive the sampling weight from the fabric's own coupling. The chain that would
have to be shown, and each link is a well-posed calculation:

- the gauge sector couples to ρ_Q = Im(Φ̄Φ̇), which is **bilinear** in the field
  and is the *only* long-range coupling the model has; **[M]**
- an interaction whose rate is set by that bilinear form samples configurations
  at a rate ∝ the projection, i.e. p = 1 rather than p = 0; **[C — to derive]**
- the resulting bias must be **source-side** (§4), i.e. established when source
  and chooser were in causal contact, not applied at readout.

---

## 6. Status after crank 2

| question | status |
|---|---|
| Is the ceiling a shape limit, not a quantum phenomenon? | **Yes** [D,M] — dimension-invariant, ℓ¹/ℓ² ratio, extremal at CHSH |
| Closed form for the ladder? | **Yes** [D] — incomplete beta, verified symbolically |
| What singles out p = 1? | **Monochromaticity** [D] — unique single-harmonic rung; and \|x\|sgn(x) = x |
| Does the fabric have the right algebraic form? | **Yes, measured** [M] — A6's cos(Δφ) law is weight × response at p = 1 |
| Detector-side realisation? | **Refuted** [D,M] — η = 0.5–0.64 < 0.828 |
| Source-side realisation derived? | **No** [C] — the remaining wall, now one exponent wide |

The gap has narrowed from "Bell is a wall" (before) to "derive the |λ·b| sampling
weight, source-side, from the gauge bilinear" (now). That is a single, well-posed
question rather than a research programme.

**Falsifiable and cheap, unchanged from the seats' recommendation:** run the
in-kernel CHSH with two objects from a common origin and a phase-overlap
readout. Predictions: S ≤ 2 for an unbiased readout (the fabric is then a
Cauchy-slice local model and Bell applies); S → 2√2 *only* if the source
ensemble itself is setting-correlated; and if S > 2 appears with a
detector-side bias, that is the **closed detection loophole** and must be
rejected, not celebrated.
