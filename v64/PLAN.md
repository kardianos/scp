# v64 PLAN — Density, Proper Time, and the Magic Sector

**Date**: 2026-06-05
**Parents**: [`../v62/THESIS.md`](../v62/THESIS.md), [`../v63/CLOSEOUT.md`](../v63/CLOSEOUT.md),
[`../v63/PREDICTIVITY_AUDIT.md`](../v63/PREDICTIVITY_AUDIT.md), [`../v60/README.md`](../v60/README.md)
**External inputs**: Quanta, *"Entanglement Builds Space-Time. Now 'Magic' Gives It Gravity"*
(2026-06-03); HN thread 48409675 (the differential-aging drift-and-spin intuition).

---

## 0. Why this version exists (the redirection)

The v59–v63 octonionic program is **closed as a constant-derivation effort**:

- **v63 predictivity audit**: honest input/output ratio **≈ 0.36×** — on the *fit* side,
  not the predictive side. Only the lepton Koide `Q=2/3` survives as striking.
- **v62 number-type no-go**: the transcendentals (`α`, `cos 2/3`) are *provably*
  unreachable by any algebraic construction.

Both documents name the only live tracks: **dynamical / RG**, or **back to the
soliton / simulation program**. v64 takes the second, with a sharp physical target
supplied by the v64 CONCEPT and corroborated by two external sources:

> **The local density `P = φ₀φ₁φ₂` sets the local rate of proper time. A gradient in
> that rate produces drift and spin. This is the equivalence-principle reading of the
> 6-field Cosserat field — and the question is whether it is gravity, or (as every
> prior effective-metric result in this project found) a too-strong, short-range
> nuclear effect.**

## 1. The two external inputs and how they map

### 1a. The HN comment — the mechanism (right idea, wrong sign)
> *"Time moves a bit faster on the left than the right; the left side traverses space
> faster, causing a drift to the left (and a spin). Add scale and distance."*

This is the differential-aging picture of gravity: **a gradient in the rate of time
is a force, plus a torque.** Two corrections we carry forward:

- **Sign**: real gravity pulls toward *slower* time (toward the mass). The v64 lapse
  `α = 1/√(1+κP²)` gives slower time at higher `P`, so `ẍ ∝ -∂ₓ ln α` points **up**
  the density gradient — *attractive*, the correct sign. The HN comment's stated
  direction (toward faster time) is backwards; its mechanism is right.
- **Spin**: the "(and a spin)" maps onto this theory's distinctive feature — the
  `η·curl θ` torsion coupling. Frame-dragging / spin-orbit is the natural home of the
  θ-sector, not an add-on.

### 1b. The Quanta article — the organizing hypothesis
Stripped of holography: **entanglement builds the geometric scaffold, but the scaffold
is inert; it takes "magic" (non-Clifford / nonstabilizer complexity) to give spacetime
springiness — the ability to bend. "Gravity comes from the mixing of the encoded
information."** This maps onto the project's own number-type map:

| Quantum gravity (Quanta) | SCP number-type map (v62/v63) | SCP dynamics |
|---|---|---|
| Clifford / stabilizer — rigid, classically simulable, no curvature | **algebraic** column (characters, Casimirs, dims) | the **free / linear** field theory (flat, inert) |
| **magic** / non-Clifford — gives springiness | **transcendental** column (`α`, loop ratios, RG / P3) | the **nonlinearity** `κ` and the `η·curl θ` torsion coupling |
| "gravity lives in the magic, not the geometry" | v63: gravity *cannot* be algebraic; only dynamical/RG is live | gravitational response should vanish in the linear limit |

**Operational identification for v64**: the *magic sector of this model is the
nonlinearity* — `V(P)` saturation (`κ`) and the `η·curl θ` coupling. A free/linear
field theory is flat and inert; `κ` and `η` are what let `P` back-react and curve
things. This is falsifiable (Thread C).

This also re-explains **why v60's scalar `g₀₀(P)` metric was LIGO-fatal**: a static
algebraic function of density is the *inert geometry*. The springiness (the ±2 tensor
modes) must come from the dynamical/torsional ("magic") sector — exactly what v60's
Urbantke induced-metric work was reaching for.

## 1c. The two-carrier hypothesis (H1) — what force is what

The magnitude pre-check (§3) showed the lapse `α(P)=1/√(1+κP²)` is ~10³⁸–10³⁹× too
strong and Yukawa-short (~0.4–0.8 fm) to be gravity. The reframe: **that is not a
failed gravity, it is the nuclear force.** Φ/c² ~ O(1) over ~1 fm *is* the nuclear
scale (well depth ~50–260 MeV, cf. Λ_QCD ~200 MeV). The model splits the two forces by
**carrier mass**:

> **H1 (two-carrier hypothesis).** The 6-field theory carries two forces in two
> sectors distinguished by carrier mass:
> - **Strong / nuclear force = time dilation of the *massive* φ-density sector.** φ is
>   massive (`m²=2.25`) ⇒ Yukawa range `~1/m_φ ≈ 0.38` fm ⇒ short-range; the lapse
>   `α(P)` is the local proper-time well; depth ~O(100 MeV) = nuclear scale. ✓ range,
>   ✓ magnitude, ✓ attractive sign.
> - **Gravity = the *massless* sector** (θ-torsion, `m_θ²=0`, and/or the v6 conserved
>   density wave). Massless ⇒ long-range; the v6 result `□δρ = -½|ω|² ⇒ δρ ~ 1/r` is
>   the candidate gravitational tail. ✓ range (1/r).

One density field `P`, two regimes: the **near field** is the local lapse well
(~1/r⁶, Yukawa → strong/nuclear); the **far field** is the massless radiative monopole
(1/r → gravity). Same source, two observables.

**What H1 buys**: it turns the project's most-repeated "failure" (every effective-metric
result is too strong / too short for gravity) into a *correct identification in the
nuclear sector*, and it divides the two hard problems by carrier mass instead of forcing
one channel to do both.

**What H1 does NOT solve (carried as open):**
- **Hierarchy**: the long-range (gravity) candidate is *also* ~10⁴⁰× too strong (v6:
  `G_eff/G_N ≈ 3.6×10⁴⁰`). H1 relocates the strong/gravity hierarchy (~10³⁸) into "why
  is the massless-sector coupling 10³⁸× weaker than the massive-sector coupling" — a
  real question the model does not yet answer.
- **Tensor**: the 1/r density tail is a scalar monopole; θ-BLV modes gave only h=0,±1.
  Gravity needs h=±2 from the v60 induced-metric construction, not the scalar tail.
- **Confinement**: a time-dilation well reproduces the *residual nucleon–nucleon Yukawa*
  attraction, not confining QCD (no color, no linear potential, no hadron spectrum). H1
  is honest at the *nuclear-force* level, not fundamental strong interaction.

**H1 predictions / tests** (cheap, in-model):
- **T-H1a (structural range) — DONE, PASS** (`mass_range_check.py`): the force range
  scales as exactly `1/m` (log-log slope vs `m²` = **−0.500** definitional, **−0.485**
  convolved-response tail e-fold), **source-independently** (same slope for source
  widths w=0.1–1.0). The range is the carrier Compton length — short range is *forced*
  by `m_φ`, not fitted. At standard `m²=2.25`: ~0.31–0.37 fm (nuclear, not
  gravitational). Corollary: the massless sector (`m_θ²=0`) has infinite range ⇒ the
  natural long-range gravity carrier. **The mass-split of H1 is internally consistent**;
  magnitude (10⁴⁰) and tensor (h=±2) remain open. *(Note: a naive "well half-width"
  measure is mass-independent because the screened response is 1/r-singular at the
  origin — the meaningful range is the tail decay length, which is cleanly carrier-set.)*
- **T-H1b (sector segregation)**: the long-range (1/r) tail must live *only* in the
  massless sector and vanish when that sector is given a mass; the massive-φ lapse must
  stay short-range. **Stage-1 analytic precheck DONE, PASS** (`th1b_sector_range.py`):
  massless→`p=−1.000` (1/r), massive `m²=2.25`→`p=−9.6` (screened); and the tail
  requires a nonzero-monopole source (the v6 Gauss bypass). **Stage 2 GPU run DONE
  (V100), SUGGESTIVE not clean**: massless θ/φ rises to the edge (long-range *shape*),
  massive θ/φ peaks & falls (localized) — H1 direction on shape. Surprise: the sectors
  differ dynamically by **radiation (massless) vs binding (massive)**, not the static
  Coulomb-vs-Yukawa range the precheck assumed; φ is not a clean Yukawa baseline
  (extended decaying neutron). **Stage 2b clean rerun DONE (compact impulse source),
  POSITIVE**: energy-weighted radius ⟨r⟩_θ massless > massive at all times (ratio
  1.25–1.50, gap grows); φ control identical between runs (validates method); massless
  ⟨r⟩_θ out-ranges massive φ while massive θ ≈ φ — H1's split confirmed at the
  range/propagation level. Static 1/r tail, 10⁴⁰ magnitude, h=±2 tensor still open. See
  `FINDINGS_H1.md`.
- **T-H1c (tensor home)**: whether the massless sector's tensor content can reach h=±2
  via the v60 Urbantke/induced-metric route — the only path to LIGO-compatible gravity.

## 2. The three threads

Each thread is cheap, non-tautological, and confronts a documented wall.

### Thread A — Corroborate the lapse (density ↔ time, *measured* not posited)
The CONCEPT posits `α = 1/√(1+κP)`. The kernel's dynamics use `den = 1+κP²`
(`scp_sim.c:409`), so the dynamically-grounded form is `α = 1/√(1+κP²)`. **Do not
posit either** — derive `α(P)` from the **linearized 6-field dispersion relation** on
a controlled density background, using the existing BLV / null-rotor machinery
(`v3/src/nullrotor_metric.c`, the `P/m` analysis).
- **Deliverable**: measured `α(P)` from the linearized dynamics vs. both posited forms;
  resolve the `κP` / `κP²` question; state the actual functional law.

### Thread B — The HN mechanism (gradient → drift + spin), from *real* dynamics
Evolve an actual wave packet / oscillon in a 1-D `tanh` density ramp — a **real 6-field
run, no imposed metric** (imposing the metric and integrating geodesics is the
recurring `P/m=2 is tautological` trap).
- Measure **drift velocity**: is it up-gradient (correct, attractive sign)?
- Measure **angular velocity**: does the torque couple to the `η·curl θ` (torsion)
  sector — i.e. does θ-rms growth track the spin?
- **Deliverable**: drift `ẍ(∂ₓP)` and spin `ω(∂ₓP)` from dynamics, vs. the Thread-A law.

### Thread C — Is the response "magical"? + per-sector segregation (H1)
Control experiment isolating the magic sector *and* testing the two-carrier split:
- **Linear truncation** (`κ→0`, and separately `η→0`) **must give zero drift.**
- Drift magnitude must **scale with `κ`** (and separately with `η`).
- **Per-sector (T-H1b)**: the short-range lapse well must live in the **massive φ**
  sector; the long-range 1/r tail must live only in the **massless** sector and vanish
  when that sector is given a mass.
- **Deliverable**: drift vs. `κ`, `η`; and range vs. carrier mass per sector. Confirms
  or kills both "gravity = the magic sector" and H1's mass-split of the two forces.

## 3. Guardrails (the walls every prior attempt hit)

- **Magnitude pre-check FIRST, before any GPU time** (`magnitude_precheck.py`, this
  step). Prior BLV results all land on the same wall: Yukawa range ~0.5–0.8 fm,
  strength **10³⁰–10⁴⁰× too strong** — nuclear scattering, not gravity (`Φ_min/c² ≈
  -0.275`). At equilibrium core density the lapse gives `Φ/c² = ln α ~ O(1)`, which is
  ~10³⁹× the gravitational `Φ/c² ~ 10⁻³⁹` at the proton scale. A one-hour analytic
  check tells us which side of the wall v64 is on before committing the version.
- **Use the torsional/tensor sector, not scalar-only** — otherwise it is LIGO-dead
  before it starts (the v60 lesson).
- **No imposed metrics on test particles** — the drift must fall out of the field
  evolution (Thread B), or it is tautological.

## 4. Honest caveats

- The **magic ↔ transcendental** mapping is a heuristic that aligns two independent
  no-go-style conclusions; it is *organizing*, not a theorem. The Quanta result is
  holographic and "step 0.5 of 5."
- Corroborating density ↔ time does **not** by itself yield gravity — prior work
  already shows it yields a too-strong short-range force. v64's value is to (a) pin the
  density↔time law quantitatively (Thread A), (b) test whether the nonlinear/torsional
  ("magic") sector can push it toward long-range (Thread C), and (c) if not, document
  *why the magic sector is short-range* — itself a v60-style principled exit.

## 5. Decision gate

After the magnitude pre-check (§3):
- **If `Φ/c²` and the force ratio are at the nuclear wall** (expected): v64 becomes a
  *characterization-and-exit* — pin the density↔time law (Thread A), confirm the
  drift+spin mechanism is real but short-range (Thread B), and write the principled
  reason the lapse mechanism is nuclear, not gravitational.
- **If there is any regime (parameter, sector, or long-range tail) where `Φ/c²`
  approaches the gravitational scale**: escalate to the full three-thread simulation
  program with GPU runs.

## 6. Artifacts

```
v64/  CONCEPT.md             — theory document (density ↔ proper time), corrected lapse
      PLAN.md                — this file
      magnitude_precheck.py  — Guardrail §3 (run first; analytic, no GPU)
      (Thread A/B/C scripts + SFA runs follow the decision gate)
```
