# v88 — Symbolic refinement by fitted back-tracing

**Date:** 2026-07-27 · Seat: grok-4.5  
**Method:** measure instrument → fit secular/effective forms with residuals → back-trace minimal structure  
**Inputs:** `fabric_harmonic` EOM (frozen into `fabric_measure.c`)  
**Data:** `v88/measure_out/*.tsv`, `fit_report.json`, `fit_summary.txt`  
**Code:** `fabric_measure.c`, `fit_harmonic.py`

Tags: **[M]** measured · **[F]** fitted · **[D]** derived · **[P]** postulated (surviving)

---

## 0. Verdict up front (blunt)

The instrument **does not produce particles**. The pushbacks are correct.

| claim in GROK_V88_DESIGN | status after measurement |
|---|---|
| P1: CV(N_c) < 0.15 multi-cell peaks | **FAILS** — smoke CV ~1.5–2.9 and rising; M6: no preferred N_c |
| Multi-mode lock vectors = species | **Artifact of truncation + spectrum** — lowest-order channels of m=3 |
| θ densification traps energy | **Decorative at design params** — θ(I=0.25) ≈ −0.11, Δω/ω ≈ 4% |
| Complete-cycle blocks exterior leakage | **Only as a time average** — instantaneous hop is *independent of detune* **[M]** |
| Adjacent locked cells form clusters | **Static ΔE is repulsive**; dynamics either **merges** (in phase, short sep) or **does nothing** (far) |

**As written the model is dead for the ontology's particles.** Hopping spreads or merges; it does not bind multi-cell interiors with detuned exteriors. What follows documents that with fitted numbers, kills dead postulates, and specifies the minimal revision the back-trace forces.

---

## 1. What was measured

| ID | experiment | output |
|---|---|---|
| **M1** | Instantaneous RHS on two cells: dI, dφ vs (I₀,I₁,Δφ,θ) | `M1_transfer.tsv` |
| **M2** | Instantaneous \|dI\|_max vs detune | `M2_secular.tsv` |
| **M2b** | Free evolve two-site ring: net ΔI over T=40 vs detune | `M2b_secular_evolve.tsv` |
| **M3** | Single-cell θ relax vs forced I; scale scan | `M3_theta.tsv`, `M3_scale.tsv` |
| **M4** | Two-lump static E(sep) − 2E₁ | `M4_interaction.tsv` |
| **M4b/c** | d(sep)/dt in-phase / antiphase | `M4b_force.tsv`, `M4c_force_antiphase.tsv` |
| **M5** | Tower modes: lock residuals + 3-wave rates | `M5_tower.tsv` |
| **M6** | Seed width → final N_above / IPR | `M6_Nc_vs_width.tsv` |
| **M6b** | Isolated multi-mode phase lock evolution | `M6b_lock_evolve.tsv` |

Parameters frozen: EPS1=0.15, EPS2=0.05, GINT=0.08, SIG=0.4, THMAX=0.8, GAM=1, KTH=1, MU=0.5, ETA=0.3.

---

## 2. Fitted forms and residuals

### 2.1 Instantaneous inter-cell transfer (M1)  **[F]**

**Fit:**
```
dI₀/dt = b₀ · 2√(I₀ I₁) sin(−Δφ)
       + b₁ · 2√(I₀ I₁) cos(−Δφ)
       + b₂ · 2√(I₀ I₁) sin(−Δφ) · θ₀
       + b₃ · Δ · sin(−Δφ)
       + b₄
```

| coeff | value | interpretation |
|---|---|---|
| b₀ | **−0.150000** | = −EPS1 (live edge on the ring) |
| b₁ | −5.2×10⁻⁹ | quadrature — noise |
| b₂ | −3.4×10⁻¹¹ | θ-modulated hop — absent |
| b₃ | −6.2×10⁻¹¹ | detune in dI — absent |
| b₄ | +6.0×10⁻⁹ | offset — noise |

- **R² = 1.000000**, **RMSE = 1.6×10⁻⁸**
- **〈dI₀ + dI₁〉 = 0** exactly (exchange)

**Back-trace:** the only forced inter-cell term is **bilinear hopping**
```
H_hop ⊃ ε₁ Σ_〈ij〉,α (z_{iα}^* z_{jα} + c.c.)   [Laplacian form]
```
with action-angle power
```
dI_i/dt = 2 ε₁ √(I_i I_j) sin(φ_i − φ_j)     (+ neighbour multiplicity).
```
No fitted need for density-modulated hop or detune-dependent *instantaneous* power.

### 2.2 Detune and “complete cycles” (M2, M2b)  **[F,D]**

**Instantaneous amplitude vs detune (M2):**
```
A_sec(Δ) ≡ max_Δφ |dI₀| = c₀ + c₁ Δ² + c₂ Δ⁴
```
| c₀ | c₁ | c₂ | CV(A_sec) | R² |
|---|---|---|---|---|
| 0.075 | ~0 | ~0 | **0.000** | 1 (flat) |

A_sec = 2·EPS1·I = 0.075 exactly, **independent of detune**.

**Symbolic average (sympy):** with free relative phase `φ̇ = Δ`,
```
ΔI₀(T) = 2ε √(I₀ I₁) [cos φ₀ − cos(φ₀+ΔT)] / Δ     (Δ≠0)
        = 2ε √(I₀ I₁) T sin φ₀                        (Δ=0)
```
envelope `∝ |sin(ΔT/2)/(Δ/2)|`. Net transfer is suppressed only after averaging a joint period — and only if amplitudes stay put.

**M2b (free evolve, two-site):** net |ΔI| still O(0.01–0.16) across the whole accessible detune window (θ-driven Δ ≲ 0.33). Lorentz fit R²=0.06 (poor — oscillatory, not Lorentz). On-resonance / far ratio only **~1.65**, not a hard gate.

**Conclusion:** the ontology's "energy crosses only in complete cycles" is **not implemented** by bilinear hop. Instantaneous leakage is always on; detune only beats the average. For particle trapping you need either:

* \|Δ\| ≫ ε₁·(bandwidth factor) so the beat is fast *and* amplitudes cannot follow, or  
* a true non-bilinear / gated transfer rule.

Current densification cannot make \|Δ\| ≫ bandwidth (see §2.3).

### 2.3 Theta sector (M3)  **[F]**

**Fit of relaxed θ vs total action in mode 0:**
```
θ(I) = −a tanh(b I)
```
| a | b | R² | RMSE |
|---|---|---|---|
| 0.575 | 0.747 | **0.9993** | small |

Design was a=THMAX=0.8, b=GAM=1. Softening a,b comes from **elastic MU** pulling against a single densified site (neighbours want θ≈0).

| I | θ_relaxed **[M]** | Δω/ω̄ = −σ θ (σ=0.4) |
|---|---|---|
| 0.25 | **−0.109** | **+0.044** |
| 0.50 | ~−0.21 | ~0.084 |
| 0.60 | −0.239 | 0.096 |

**Decorative at design params.** Smoke's θ ~ −0.10 to −0.17 is exactly this regime.

**Scale to engage trapping:** want |θ| = 0.5 at I = 0.25 with THMAX=0.8:
```
0.8 tanh(GAM · 0.25) = 0.5  ⇒  GAM = atanh(0.625)/0.25 ≈ 2.93
```
Only **~3×** current GAM — not orders of magnitude. But even then
```
Δω/ω̄ = σ · 0.5 = 0.20  ⇒  Δω ~ 0.2  (mode 0)
```
while hop half-bandwidth on cubic d=3 is `2 d ε₁ = 0.90`.  
**Density-shifted frequency still sits inside the hop band.** θ cannot reflect energy at the exterior until

```
σ · |θ|_core  · ω̄  ≳  2 d ε₁
⇒  |θ|_core  ≳  2 d ε₁ / (σ ω̄) ≈ 0.90 / (0.4 · 1.03) ≈ 2.2
```
which **exceeds THMAX=0.8**. With present (σ, THMAX, ε₁, d) the density sector **cannot** produce exterior detuning outside the band. Raising GAM only saturates earlier at THMAX; it does not fix the inequality.

**Forced fix (one of):** lower ε₁, raise σ, raise THMAX past O(1), or change the hop structure — not "anneal harder."

### 2.4 Two-lump interaction (M4)  **[F]** — question (a)

**Static excess energy** ΔE(s) = E(s) − 2E₁ after short relax:

| sep | ΔE **[M]** |
|---|---|
| 2 | **+3.59** |
| 4 | +1.61 |
| 6 | +0.40 |
| 10 | +0.006 |
| ≥14 | ~0 |

**Fit** (R² = **0.996**):
```
ΔE(s) = a₀ + a₁ e^{−s/ξ} + a₂ e^{−s/ξ} cos(π s) + a₃ e^{−s/ξ} cos(2π s/3) + a₄/(s+½)
```
with ξ ≈ O(1–few) lattice spacings. **ΔE(near) > ΔE(far)** ⇒ **static interaction is REPULSIVE.**

**Dynamics:**

| protocol | 〈d(sep)/dt〉 | reading |
|---|---|---|
| in-phase (M4b) | −7.0×10⁻³ | short sep: centroids **collapse to one object** (sep 3→0.09); far: dsep/dt≈0 |
| antiphase (M4c) | >0 (repulsive) | antibonding |

In-phase "attraction" is **not binding of two particles**. It is **merger / delocalisation into one blob** when hop energy rewards phase coherence. Once separated by ≳ 12 sites there is **no force**.

**Answer (a):** fitted inter-lump potential is **short-range repulsive** (static). There is **no bound two-lump molecule**. Multi-cell *particles* cannot form by adhesion of locked cells in this H. The model is **dead as written** for composite / multi-cell assembly via cell–cell attraction.

Symbolic two-site hop energy (forced by M1):
```
E_hop = ε (I₀ + I₁ − 2√(I₀ I₁) cos Δφ)
```
Minimised by coherent equal amplitudes on neighbours — that **spreads** action onto the neighbour, opposite of isolating a cluster from exterior vacuum.

### 2.5 What sets N_c (M6)  **[F]** — question (b)

After T=20 from Gaussian seeds of width W and amplitude A:

| observation | value |
|---|---|
| mean N_above | 5.78 |
| corr(N, W) | **0.08** (weak) |
| dN/dW by amplitude | ~0, −0.17, +0.62 |
| preferred_size_exists | **False** |

At higher amplitude, objects collapse toward **N ≈ 3** (few-site breathers). At low amplitude, N tracks fluff of the seed. **Nothing selects a quantised multi-cell N_c independent of initial width/amplitude.**

**Answer (b):** N_c is set by **initial condition + dispersive/self-trapping balance of a discrete NLS-like core**, not by an integer closure condition. Fit: no stable preferred size. P1 cannot pass on this dynamics.

### 2.6 Tower choice (M5)  **[M]** — question (d)

| tower | (ω̄₀,ω̄₁,ω̄₂) | \|2ω₀−ω₁\| | \|−ω₀+2ω₁−ω₂\| | \|dI\|_rms (3-wave @ fixed phase) |
|---|---|---|---|---|
| detuned (instrument) | 1.03, 2.12, 3.27 | **0.060** | 0.060 | 2.034×10⁻³ |
| commensurate | 1, 2, 3 | **0** | 0 | 2.034×10⁻³ |
| strong detune | 1.17, 1.91, 2.53 | **0.424** | 0.130 | 2.034×10⁻³ |

Instantaneous 3-wave rates are **identical** (GINT acts on z, not on ω). Lock *residuals* track the tower. So:

* Tower choice **controls which integer channels are near-secular**, not the instantaneous vertex.  
* Labels (2,−1,0) and (−1,2,−1) are exactly the **lowest-order relations in an m=3 tower** — tension #2/#3 confirmed.  
* **m-sensitivity is expected:** change m, change the lowest channels, change the "species" histogram. That is truncation physics, not fabric physics.

**Answer (d):** commensurate vs detuned **does** control secular lock accessibility. It is **not irrelevant**. It is also **not a species mechanism** until the spectrum is derived from cell geometry rather than chosen as a tower of length m.

### 2.7 Intra-cell lock (M6b)  **[M]**

Isolated cell, hop off, multi-mode IC:
* late ψ = φ₀+φ₁−φ₂ has **std 2.44** (winds)  
* **not locked**  
* θ → −0.45 (here I is large enough for moderate densification)

GINT is too weak / detune too large for a stable relative-phase fixed point at design params. The "emergent Higgs configuration" **does not form**.

---

## 3. Back-traced minimal structure

### 3.1 Forced by the fits (keep)

| term | evidence |
|---|---|
| Complex mode amplitudes z_α with bilinear neighbour hop ε₁ | M1 R²=1, b₀=−ε₁ |
| Action-angle exchange dI ∝ √(IᵢIⱼ) sin Δφ | M1 |
| θ dynamics toward −THMAX tanh(GAM Σ w_α I_α) with elastic MU | M3 R²=0.999 |
| Free rotation ω_α(θ)=ω̄_α(1−σ θ) | used in M2/M3; consistent |

### 3.2 Redundant for current (null) phenomenology

| term | why redundant |
|---|---|
| **EPS2** (2:1 inter-cell) | M1 explained without it; smoke "species" are not 2:1 composites |
| **GINT** (3-wave) | no lock (M6b); does not set N_c (M6) |
| Multi-mode m≥2 as species engine | single-mode hop+θ already gives the same failure mode (few-site peaks) |
| Weakly detuned tower as "generic spectrum" | privileges m-dependent channels (M5) |
| Annealing schedules as discovery tool | nothing to discover until binding exists |

### 3.3 Dead weight postulates **[e]**

From DESIGN §, after back-trace:

| postulate | fate |
|---|---|
| Discrete fabric | **SURVIVES** as necessary for integer N_c *in principle* |
| Complete-cycle transfer as resonant average | **SURVIVES as math**; **FAILS as implemented trap** (instantaneous hop always on) |
| No kernel Higgs field | **SURVIVES** (prohibition still right) |
| Multi-D harmonics ⇒ integer species via lock vectors | **DEAD as currently set up** — truncation channels, no lock, no multi-cell |
| More cyclic energy → tighter cell | **SURVIVES as formula** (M3); **DEAD as trapping engine** at reachable (θ,σ,ε₁) |
| Particle = interior config + sharp exterior mismatch | **DEAD in instrument** — no interior/exterior structure measured |
| Annealing finds species | **DEAD until binding exists** |
| Weakly detuned tower | **DEAD weight** — replace with geometry-derived spectrum or drop |

**Short list that remains [P]:** discrete cells; geometry θ (preferably with Σθ=0); cyclic energy densifies; transfer should be cycle-complete in a sense strong enough to isolate a region; particles are multi-cell with exterior mismatch.  
**Everything about m-mode lock-vector species is out until rebuilt from geometry.**

### 3.4 Missing (required additions) — forced by failures

To get multi-cell particles the fits say you must add at least one of:

1. **Spectral wall:** make core–exterior detune exceed hop bandwidth  
   `σ · THMAX · ω̄ ≳ 2 d ε₁`  
   (currently 0.4·0.8·1 ≈ 0.32  ≱  0.9).

2. **True complete-cycle gate:** transfer kernel that is *structurally* O(ε²/Δ) or stroboscopic, not bilinear O(ε) instantaneous. Bilinear hop cannot implement the ontology literally.

3. **Binding energy for multi-cell cores:** static ΔE must go **negative** at small sep for in-phase densified regions. Currently ΔE>0. Candidate mechanism already known in this project: **geometric Σᵢ θᵢ = 0** (fabric_mass) so densification is globally scarce and packing of dense cells can lower elastic cost — *absent* in scalar free θ.

4. **Scale selector:** a term whose stationary condition fixes N_c (e.g. balance of surface energy ~ N_c^{(d-1)/d} against bulk cycle gain ~ N_c, or a discrete PN barrier with a preferred width in units of a). Hopping alone does not.

Without (1) or (2), exterior always leaks. Without (3)+(4), no multi-cell particle, only breathers or fluff.

---

## 4. Answers (a)–(e) compact

**(a) Adjacent locked-cell interaction?**  
Fitted static form: short-range **repulsive** ΔE(s)>0, R²=0.996. Dynamics: in-phase **merger** at contact, **no force** when separated; antiphase repulsive. **Clusters of the ontology kind cannot form. Model dead as written for multi-cell assembly.**

**(b) What sets N_c?**  
**Initial width/amplitude + NLS-like self-trapping**, not an integer closure. preferred_size_exists = **False**. Mean N_above ~ 6 with large scatter; high-amp cores → N~3 breathers.

**(c) Does θ work?**  
**Decorative** at design params (θ(0.25)=−0.11). Fit θ=−0.575 tanh(0.747 I), R²=0.999. GAM~3 gets |θ|~0.5, but **even saturated THMAX cannot push Δω outside the hop band** under current (σ,ε₁,d). θ is not a small tweak away from working — the inequality is structural.

**(d) Commensurate vs detuned tower?**  
**Controls secular lock residuals** (0 vs 0.06 vs 0.42). **Irrelevant to instantaneous vertices.** **Dangerous as species source** — labels track m and the tower choice (tension #2/#3 confirmed).

**(e) Postulate kill list?**  
Survives: discrete fabric; densify-with-cycle (formula); no imported Higgs; resonant-average reading of complete cycles (as math).  
Dead: lock-vector species from m-truncation; θ-as-trap at present scales; annealing-as-discovery; multi-mode complexity as the particle engine; bilinear hop as complete-cycle rule.

---

## 5. Revised model spec (minimal, back-trace driven)

### 5.1 What to delete

* EPS2, GINT, and m-mode "species" diagnostics as primary observables  
* Weakly detuned harmonic tower as a physics input  
* Any claim that smoke/quench histograms are particles  

### 5.2 What to keep

* Discrete lattice, spacing = physical constant  
* One primary cyclic complex field per cell **ψ** (start m=1; re-introduce internal spectrum only from a geometric operator on the cell)  
* Densification θ with **geometric** constraint from ξ: θ = div ξ, hence **Σ θ = 0** on periodic domains  
* Cycle → tightness: θ_eq driven by |ψ|², same tanh saturation, but parameters must satisfy the spectral-wall inequality below  

### 5.3 What to add (mandatory)

**A. Spectral wall inequality (design constraint, not a free guess):**
```
σ · THMAX · ω₀  ≥  κ · 2 d ε₁     with κ ≥ 1 (target 2)
```
Example feasible point: d=3, ε₁=0.05, σ=0.8, THMAX=1.2, ω₀=1 → LHS=0.96, RHS=0.30√ — wall possible.  
*Derive ε₁ from measured requirement, do not inherit 0.15.*

**B. Geometric packing (binding):** evolve ξ, not free θ. Reuse fabric_mass elastic + double-well or cycle-driven well, so multi-cell dense regions can be **energetically preferred** to N isolated densified cells (surface vs bulk). **Measure ΔE(s) again; success = ΔE(near)<0.**

**C. Transfer rule upgrade (pick one and test with M1/M2 protocol):**

| option | rule | test |
|---|---|---|
| C1 | Keep bilinear hop but enforce spectral wall (A) so exterior is evanescent | M2b: \|ΔI\|_far / \|ΔI\|_on → 0; multi-cell seed retains I |
| C2 | Replace hop with resonant normal-form only: ∂ₜIᵢ = ε² f(I) sin(ψ)/Δ_eff with hard cutoff if \|Δ\|>Δ_* | M1: instantaneous dI≈0 when detuned |
| C3 | Stroboscopic: exchange kicks only when cumulative phase re-aligns within δ | same as C2 |

Ontology prefers C2/C3; C1 is the minimal change if one insists on Hamiltonian hop.

**D. Observable that defines a particle (replace lock-vector fingerprints):**
```
particle ⇔ connected dense component with
  N_c ≥ N_min (e.g. 4),
  interior |θ| > θ_wall,
  exterior shell with hop-band mismatch,
  lifetime ≫ 2π/ω₀,
  same (N_c, I_tot/N_c) bin in ≥3 seeds.
```
Drop (2,−1,0) labels until an internal spectrum is geometry-derived and locks with measured phase variance → 0.

### 5.4 Revised first instrument (spec only; implement next)

`fabric_core.c` (proposed):

* DOF: ψ_i ∈ ℂ, ξ_i ∈ ℝ^d (or θ_i with Lagrange projection Σθ=0)  
* ε₁, σ, THMAX satisfying §5.3.A  
* No EPS2/GINT  
* Diagnostics: M1–M4 suite must pass gates before any anneal  

**Acceptance gates (quantitative):**

| gate | pass criterion |
|---|---|
| G-θ | θ(I=0.25) ≤ −0.4  **and** σ\|θ\|ω₀ ≥ 2 d ε₁ |
| G-bind | ΔE(sep=2) < 0 for two densified seeds (static) |
| G-Nc | after quench, histogram of N_c has a peak with CV<0.15 in ≥3 seeds, peak N_c≥4 |
| G-leak | seed with interior wall loses <5% action to exterior over T=100 |
| G-m | if m>1 reintroduced, species bins stable under m→m+1 |

Until G-bind and G-Nc pass, **do not run annealing programmes**.

### 5.5 What the revised theory does *not* claim

* It does not claim multi-dimensional harmonics are wrong as ontology — only that the **implemented tower truncation does not realise them**.  
* It does not claim discrete breathers are the answer — G-Nc with N_c≥4 and interior/exterior structure is required specifically to reject breather mimicry.  
* It does not lengthen the postulate list; it shortens it and adds **one inequality** and **one geometric constraint** forced by measurement.

---

## 6. Relation to pushbacks (explicit)

1. **P1 fails / CV rising / N_c~1:** confirmed; mechanism is hop-dominated single-site/few-site peaking with decorative θ.  
2. **Lock vectors = truncation:** confirmed by M5; m-sensitivity expected.  
3. **Single cells ≠ particles:** confirmed; no interior/exterior structure in smoke or M6.

---

## 7. Files

| file | role |
|---|---|
| `v88/fabric_measure.c` | controlled measurements M1–M6 of frozen EOM |
| `v88/fit_harmonic.py` | fits + sympy back-trace → `measure_out/fit_report.json` |
| `v88/measure_out/*` | raw TSV + fit_summary |
| `v88/GROK_V88_SYMBOLIC.md` | this document |
| `v88/GROK_V88_DESIGN.md` | prior design; species/annealing sections superseded until gates pass |
| `v88/fabric_harmonic.c` | cautionary forward model (keep for regression, not as particle engine) |

---

## 8. Immediate next step (only this)

Implement `fabric_core.c` with §5.3.A–B (spectral wall + Σθ=0 geometry), rerun M3/M4/M6 gates.  
**If G-bind fails under Σθ=0, say so and stop — the density packing idea is then also dead and the ontology needs a different binding mechanism.**  
Do not touch `scp_sim` / `sfa.h`. Do not anneal.
