# v87 B0-G, crank 2 — why p = 1, and how it fits the fabric
**Date:** 2026-07-26 · Follow-on to the three-seat rung (`FABLE_WORK.md`,
`GROK_WORK.md`, `GEOM_WORK.md`). Worksheets: `work/geom/nd_shape_limit.py`,
`work/geom/why_p1.py`, both run clean end to end.

Tags: **[D]** derived · **[M]** measured · **[P]** postulated · **[C]** conjecture

> **REVIEWED AND CORRECTED, 2026-07-26.** An independent grok-4.5 seat re-ran both
> worksheets and returned **SOUND WITH CORRECTIONS** with eight required fixes
> (`crank2_review.md`). All eight are applied below. The pre-review text is kept
> at `CRANK2_RESULTS_prereview.md`. Three of the corrections retract claims that
> were central to the original write-up: the monochromaticity inference (a pun),
> the A6 "read literally" identification (a tautology plus a category error), and
> the closing sentence that the gap was "one exponent wide" (it is not).

---

## 0. Where crank 1 left it

The three seats established: the ontology relaxes **(MI)**; explicit local
deterministic models reach exactly 2√2 at a cost of **0.046 bits**; and the
Tsirelson half was **not** explained, because measurement dependence gives a
smooth dial from 2 through 2√2 to 4 with nothing distinguishing the quantum
point.

Crank 1 then showed the ceiling is a **shape limit**, not a quantum phenomenon:

* E(θ) is **exactly independent of the dimension N** of the body λ lives on —
  the radial integral cancels identically. Verified symbolically (N = 3, 4, 5,
  10) and by Monte Carlo (N = 2…25, error ≤ 4×10⁻⁴). **[D,M]**
  **SCOPE (correction 1):** this holds only for responses and weights that
  *factor through the angle in the a–b plane* — sgn(a·λ) and |λ·b|^p do. A
  response sensitive to the radial component (e.g. a soft or amplitude-dependent
  readout) breaks it, and the reviewer confirmed the break numerically. The
  invariance is a property of *this response class*, not of hidden-variable
  bodies in general.
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

## 3. Uniqueness: p = 1 is the only single-harmonic rung *of this family*  **[D,M]**

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
single-angular-harmonic member of the family. **Read §3.1 before drawing any
inference from this to the fabric — the obvious one is false.**

### 3.1 The fabric inference was a PUN, and it is false  **(correction 2)**

The original text argued: the fabric's objects are monochromatic (Φ = f e^{iωt},
and HC-1 measured ≤1 bound internal mode, none in the working region), therefore
p = 1. **That inference is wrong, and it equivocates on two senses of the word:**

| sense | object | what is single |
|---|---|---|
| **temporal** | Φ_a = f(r)e^{iωt}; HC-1 modes | one **clock frequency** |
| **angular** | E_p(θ) = Σ c_k cos(kθ) | one **harmonic in setting angle** |

**Decisive counterexample (reviewer's).** Take a temporally monochromatic rotator
λ — a single phase — with responses A = sgn cos(λ−a), B = −sgn cos(λ−b) and
**uniform** sampling (p = 0). The object is perfectly monochromatic in time. The
correlation is the **triangle wave**, c₁ ≈ 0.81 with large odd harmonics, and
S_max = **2**. Classical.

So: *monochromatic object ⇒ p = 1* is **false**. Monochromatic object + uniform
λ gives the classical triangle. What selects the single **angular** harmonic is
the **sampling weight** |λ·b|, not the single temporal clock.

What survives is the geometry alone: **within the p-family**, p = 1 is the unique
single-harmonic rung **[D]**, and the rectifier identity |x|sgn(x) = x is a real
mechanism for why that weight reconstitutes a linear projection **[D]**. Neither
licenses any inference from HC-1. The fabric link is **[C], and as an implication
it is retracted.**

**Also (correction 4a):** p = 1 is unique *inside this power family only*. Other
measurement-dependent constructions outside the family (Hall-type, and related
constructions the reviewer flagged) also reach E = −cos θ. The family is a
convenient parametrisation, not an exhaustive one.

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

**Correction 6 — the fork is directionally right but was overstated as clean.**
The binary above omits hybrids. With P(detect | λ,b) = c + (1−c)|λ·b| the
reviewer measured:

| c | η_B | S |
|---|---|---|
| 0.00 | 0.50 | ≈2.83 |
| 0.70 | 0.85 | ≈2.15 |
| 0.85 | 0.93 | ≈2.07 |
| 1.00 | 1.00 | 2.00 |

with S ≈ 2.09 still obtainable at coincidence efficiency ≈ 0.90. So the correct
statement is narrower than the original: **pure |λ·b| detection cannot be the
data-generating process for high-η, near-Tsirelson results** — it is not true
that detection bias cannot produce *any* S > 2. The 0.828 figure is the threshold
for reaching the *quantum point at ideal visibility*, not a prohibition on all
violation. (Clauser–Horne/Eberhard thresholds are lower, ~2/3; pure Reading 2 at
0.5–0.637 sits below those too.) The threshold identity itself is standard but
**UNVERIFIED against a primary source this session.**

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

The fabric exhibits, as a measured interaction law, force ∝ cos(Δφ), which can
be written |cos Δφ|·sgn(cos Δφ).

**Corrections 3 and 4 — this is much weaker than the original claimed.** The
original said this "is precisely weight × response at p = 1… not a fitted
analogy; it is A6 read literally." Both halves are wrong:

* **The factorisation is a tautology.** For *any* real x, x = |x|·sgn(x). Any
  force law ∝ f(Δφ) admits the identical split. Writing cos = |cos|·sgn(cos)
  demonstrates nothing about how the fabric samples. It is reverse-engineered
  bookkeeping. **[D]**
* **Δφ is not λ·b — it is a category error.** Δφ is a *dynamical* relative clock
  phase evolving under the PDE; b is a *freely chosen setting*. One is the
  argument of a force law, the other the argument of a conditional measure
  ρ(λ|b). O2 does make choosers fabric and hence dynamical — but that means the
  chooser functional must be **constructed**, not that the two symbols may be
  silently renamed. **[D]**

The row "weight = magnitude of the overlap **[M]**" in the table above is
therefore **mis-tagged**: A6 itself is [M]; the map from A6 to a Bell
weight × response is **[C]**.

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

### 5.3 What would close it — the original chain was a non-sequitur  **(correction 5)**

The original proposed: gauge couples to ρ_Q = Im(Φ̄Φ̇), which is bilinear,
therefore the sampling rate ∝ the projection, therefore p = 1. **Retracted. It
fails twice:**

* **It points at the wrong A6 sector.** A6 splits into Coulomb/gauge, which is
  long-range and explicitly *internal-state-blind*, and phase-coherent contact,
  which carries the cos(Δφ). The structure crank 2 wants lives in **contact**,
  not in the gauge bilinear. Citing ρ_Q targets the force that by construction
  cannot see the phase. **[D]**
* **Bilinearity does not imply a sampling law.** For a monochromatic ball
  ρ_Q = ω|Φ|² — a scalar charge density, not a projection onto an external
  setting direction. Nothing about it determines the measure ρ(λ|b). **[D]**

**The real checklist**, which is a derivation programme and not a one-liner:

1. Define settings a, b as **fabric functionals** — actual chooser subsystems.
2. From the joint initial-value problem (source + two choosers + common past),
   compute or prove ρ(λ|b) ∝ |λ·b|, or an equivalent tilt giving E = −cos θ,
   **at production**, not at readout.
3. Show the tilt **survives** free evolution and scrambling out to spacelike
   measurement — this is GROK's G14 obstruction, unaddressed.
4. Map continuous fabric readouts to dichotomic ±1 outcomes **without**
   reintroducing a detection loophole.
5. Show the dynamics **forbid** nearby rungs with S > 2√2, or concede that
   Tsirelson remains imported.

### 5.4 What a "setting" would even be  **(correction, omission 7)**

The document did not define one, and the kernel does not contain one. A setting
requires at minimum: a fabric degree of freedom playing the analyzer
orientation; quasi-independent choice of it at each wing; a dichotomic readout
map A(a,λ) ∈ {±1}; and spacelike separation with local response. **None of the
four exists in A6, which is contact physics at small D.**

---

## 6. Status after crank 2

| question | status |
|---|---|
| Is the ceiling a shape limit, not a quantum phenomenon? | **Yes** [D,M] — dimension-invariant, ℓ¹/ℓ² ratio, extremal at CHSH |
| Closed form for the ladder? | **Yes** [D] — incomplete beta, verified symbolically |
| What singles out p = 1? | **Inside the p-family**: unique single-harmonic rung, plus \|x\|sgn(x) = x [D]. The monochromaticity inference is **retracted** — it was a pun (§3.1) |
| Does the fabric have the right algebraic form? | **No as a Bell correlator; yes as force algebra** [M for A6, C for the map] — and the factorisation is tautological |
| Detector-side realisation? | **Pure form refuted** [D,M] — η = 0.5–0.64 < 0.828. Hybrids survive at lower S (§4) |
| Source-side realisation derived? | **No** [C] |

**Correction 7 — "the gap is one exponent wide" is retracted as not honest.**
The residual gap is multi-dimensional. The reviewer's enumeration:

| # | gap | status |
|---|---|---|
| 1 | What *is* a setting in the kernel? | undefined |
| 2 | Source-side ρ(λ\|a,b) of quantum form | postulated |
| 3 | Survival under transport / mixing (G14) | open obstruction |
| 4 | Continuous force/phase → dichotomic outcomes | unbuilt |
| 5 | Why exactly p = 1 and not 0.9 or 1.1, from dynamics | open |
| 6 | Why not S → 4 (a Tsirelson principle) | still imported |
| 7 | Spacelike separation with local responses | only by construction |
| 8 | Gauge vs contact mis-identification | wrong target coupling |

What crank 2 actually contributed: it closed the pure detection reading, gave a
clean closed form for the ladder, and isolated the rectifier identity. It
**localised one algebraic preference inside an already-postulated family** — it
did not reduce the derivation wall to a single exponent.

**Falsifiable and cheap, unchanged from the seats' recommendation:** run the
in-kernel CHSH with two objects from a common origin and a phase-overlap
readout. Predictions: S ≤ 2 for an unbiased readout (the fabric is then a
Cauchy-slice local model and Bell applies); S → 2√2 *only* if the source
ensemble itself is setting-correlated; and if S > 2 appears with a
detector-side bias, that is the **closed detection loophole** and must be
rejected, not celebrated.
