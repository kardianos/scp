# v65 FINDINGS — Back-calculating a stable structure (oscillon route)

**Date**: 2026-06-05
**Question (user)**: can we *estimate* what structure / harmonics / phased frequency would
be stable, and back-calculate a seed from there — instead of searching or controlling?

## The theory (sound) — `stability_estimates.py`

- **Derrick**: no *static* soliton is stable (`E(λ)=λE₂+λ³E_V` has only a maximum). This is
  exactly the collapse/disperse instability every self-tuning run hit. Rigorous.
- **Escape = time-periodicity**: a phased frequency ω gives an effective-mass barrier
  `m²−ω²` that holds the lump open ⇒ an **oscillon**. Band: ω < m = 1.5. Size
  `R(ω)=1/√(m²−ω²)`; tuned amplitude from the core nonlinear balance.
- **Diagnosis of the old failure**: the previous seed (A=1.0, σ=1.5) is **over-driven** —
  σ=1.5 matches ω≈1.34, whose tuned amplitude is **A≈0.62, not 1.0**. The excess radiates.
- **θ loss channel**: massless θ is driven at ω with k_θ=ω ⇒ it radiates; prefer lower ω.

## The test (negative) — tuned Gaussian seeds still decay

Four seeds at κ=50, absorbing far boundary, no self-tune, fine diagnostics:

| seed | A | ω | P_max>0.01 until | outcome |
|---|---|---|---|---|
| base | 1.00 | ~1.34 | t≈4 | radiates away |
| w110 | 0.76 | 1.10 | t≈1.5 | radiates away |
| w125 | 0.69 | 1.25 | t≈3 | radiates away |
| w135 | 0.61 | 1.35 | t≈3.5 | radiates away |

**All decay within a few time units; none is a coherent oscillon.** P_max oscillates (period
~2 t.u., the 3ω product harmonic) under a decaying envelope. The amplitude-vs-ω trend is weak
and noisy (no clean "lower ω lasts longer"). **Amplitude tuning alone is not enough.**

## Interpretation

The estimate identified the right *concept and ballpark* but a **Gaussian-from-rest seed is
too crude** to land on the oscillon limit cycle (wrong radial profile, no launch velocity,
and the product structure P=φ₀φ₁φ₂ couples the three components' phases). Two possibilities,
and the next test distinguishes them:
1. The exact oscillon **exists** but needs its true profile `f(r)` + velocity launch (a
   Gaussian is far from it).
2. The **product-potential theory has no oscillon** (the P=φ₀φ₁φ₂ structure may not support a
   coherent breather), and stability needs an *explicit* ingredient.

## Where this leaves the whole self-tuning arc

Every route this session converges on one fact: **this 6-field product-potential theory has
no stable localized soliton** at standard parameters — static is a Derrick saddle, oscillons
(as seeded) radiate, and feedback (SOC / density homeostasis) cannot *manufacture* an
attractor that isn't there. The self-tuning program (X1 mechanism proven, kernel `self_tune`
implemented) is **blocked on this substrate**, not on the control logic.

## The fork (decision)

1. **Definitive oscillon test** — derive & solve the **radial oscillon BVP** for *this*
   theory (single-frequency reduction of the 6-field equations, fixed ω), get the exact
   `f(r)` + launch velocities, seed and test. Local computation (no GPU) for the profile,
   one short GPU run to test. *If even the exact profile radiates, that is strong evidence
   the theory lacks oscillons* — a real, near-conclusive result.
2. **Explicit stabilizer** (kernel physics): add a Skyrme-type `L₄` term (the project has
   L₄/L₆ analysis), or a genuine **Noether-charge Q-ball sector** (a complex/internal phase
   with a conserved current in the EOM, which *provably* evades Derrick). Bigger change;
   needs authorization and a physics choice.
3. **Reconsider the target**: whether stable *static-like* particles are the right object in
   this theory, vs. accepting the solitons are dynamical/metastable and reframing what
   "particle" means here (ties back to the v62/v63 "is a stable particle emerging" criterion).

**Recommendation**: (1) — it is the cheapest decisive test and tells us whether to invest in
(2). Honest caveat: deriving the product-potential oscillon BVP is real work, and a null
result there would point to (2)/(3) as the only ways to get a stable particle.

**Artifacts**: `stability_estimates.py`; four oscillon run logs (lifetimes above).

---

## UPDATE — the exact BVP route solved it: there is NO φ-oscillon, and we know WHY

Built a 1+1D spherically-symmetric radial solver for the φ-sector (`radial_oscillon.c`,
exact EOM `φ̈ = φ'' + (2/r)φ' − m²φ − μ φ_a∏_{b≠a}φ_b²/(1+κP²)²`) and scanned exhaustively
(local, no GPU): amplitude A∈[0.3, 3.5], width σ, frequency ω∈[1.0, 1.45], equal vs
offset δ, rest vs time-phase launch, with/without bulk cooling. **Every configuration
radiates away within a few periods** (P_max → <1% of initial). The exact-profile route did
not find an oscillon — and revealed the mechanism that forbids one:

**The cubic source P=φ₀φ₁φ₂ pumps third-harmonic radiation above the mass gap.**
- Measured core oscillation: |P| period ≈ 1.87 t.u. ⇒ field frequency `ω_φ = π/1.87 ≈ 1.68`,
  matching `ω_φ = √(m²+k²)` with `k∼1/σ` (the profile's gradient content).
- `ω_φ ≥ m = 1.5` *always* (gradient energy only raises it). A field at ω_φ feeds the cubic
  product a **3ω_φ ≈ 5 ≫ m** component — far above the mass gap, where it is a **propagating
  (radiating) mode**, not an evanescent bound one. Suppressing it would need `3ω_φ < m`, i.e.
  `ω_φ < 0.5` — unreachable (the field can't be softened a factor 3 below its mass against
  saturation).
- This is **intrinsic to the product structure** — independent of A, σ, phase, cooling — and
  it is why oscillons die in a few periods here (vs 10⁴–10⁶ periods for a φ⁴ oscillon).

**Combined with Derrick (no static soliton), the φ-sector has NO stable localized structure,
for a definite physical reason** — not a numerical artifact, not bad seeds.

## What this proves about the resolution (it is now well-motivated, not a guess)

The radiation channel is `d|P|/dt ≠ 0` driving above-gap harmonics. The structure that
**eliminates it** is a **Q-ball**: a complex/internal phase rotation `φ ~ f(r)e^{iωt}` whose
modulus `|φ|²=f(r)²` is **time-independent**, so the potential it sources is *static* — **no
harmonic radiation at all**. A conserved Noether charge then evades Derrick (Coleman). This
is exactly **option (2)** from the prior fork, now selected *for a reason*: it is the unique
fix that cures the specific loss mechanism we identified.

**Concrete next step**: introduce a U(1)/internal phase-rotation sector (a genuine conserved
current in the EOM) — the field's *modulus* localizes (static, radiation-free) while its
*phase* carries the stabilizing frequency. This is a kernel-physics change (authorization +
a definite ansatz choice), but it is now the *derived* answer to "what structure is stable,"
not a shot in the dark. Estimate of the Q-ball frequency band and profile is the next
back-calculation (same `stability_estimates.py` machinery, Coleman thin/thick-wall).

**Artifacts (update)**: `radial_oscillon.c` (the solver + exhaustive scan); the |P|-period
measurement above.

---

## UPDATE 2 — Maxima quantified the saddle; the cascade test pins what's actually needed

**Maxima (`saddle.mac`) — the 3D saddle, rigorously:**
- Derrick scaling `E(λ)=aλ+bλ³`; at any bound equilibrium **`E″=−2a<0` ⇒ SADDLE** (symbolic).
- At standard constants `a=16.1, b=37.4`; `b>0` ⇒ saturation *prevents* static binding (the
  field disperses — there isn't even a bound saddle, the minimum is λ=0).
- **Add a conserved density charge `K/λ³`** (the cascade product): `E=aλ+bλ³+K/λ³` has a
  **guaranteed stable minimum for any K>0** (all-positive). Tuning curve `K(λ)=bλ⁶+(a/3)λ⁴`,
  mass `(4/3)aλ+2bλ³`; a λ≈0.5–0.7 particle needs `K≈1–6`. **So a conserved charge doesn't
  just lift the saddle — it creates the bound particle.** This is rigorous and is the
  positive core of the cascade idea.

**Simulation (`radial_cascade.c`) — the cascade as a *slaved density* does NOT realize it.**
Built the φ-sector + a density ρ and tested every faithful coupling: dynamic cascade
(ρ̇ ∝ φ-activity), seeded conserved charge, Thomas-Fermi binding (`ρ∝φ²+μ` at fixed `∫ρ=Q_d`
with pressure `c_ρ`), attractive vs repulsive back-reaction, with and without cooling.
**All decay** (one apparent stabilization was a revival artifact at short T). The reason is
sharp: a density *slaved to φ* (`ρ∝φ²`) makes the back-reaction `∝φ³` — **self-focusing, which
adds to the collapse**, and it gives φ no independent outward pressure. Maxima's `K/λ³` is the
*charge's own* compression energy acting on the whole object; a hand-coupled slaved field
doesn't reproduce that structure.

**What this proves about the realization (decisive):** the stabilizing charge must be a
**genuine Noether current**, not a slaved density. The clean object that *does* have energy
`= a λ + bλ³ + Q²/(c λ³)` is a **Q-ball**: `φ ~ f(r)e^{iωt}`, where the U(1) charge
`Q=ω∫f²` supplies exactly the `K/λ³` pressure **and** the modulus `|φ|²=f²` is time-independent
⇒ **no 3ω radiation** (it simultaneously cures the oscillon's loss channel). Every thread now
points to the same structure.

## Final synthesis of the stability investigation

The 6-field φ-sector has **no stable localized structure**, for two independent, proven reasons:
1. **Static** → Derrick saddle (`E″=−2a`, Maxima), and saturation prevents binding (`b>0`).
2. **Oscillon** → the cubic `P=φ₀φ₁φ₂` pumps `3ω_φ≈5 ≫ m` radiation above the gap (`radial_oscillon.c`).

Both are cured by the **same** ingredient — a genuine conserved charge realized as an internal
phase rotation (**Q-ball**): it supplies the Maxima `Q²/λ³` stabilization *and* removes the
harmonic radiation. The user's cascade intuition (φ→θ→density stabilizes) is **correct in
energy bookkeeping** (Maxima confirms a charge stabilizes) but requires the density to be a
*true Noether charge with its own conserved current*, not a field slaved to φ.

**The concrete, now fully-motivated next step**: add a U(1)/internal phase-rotation sector to
the kernel so the modulus localizes (radiation-free, Derrick-evading) while the phase carries
the charge `Q`. Front it with the Coleman thin/thick-wall Q-ball estimate (frequency band +
profile) using `saddle.mac`'s `a,b`, so the seed is tuned. This is a kernel-physics change
(authorization + the U(1) ansatz), and it is the unique resolution all five investigation
threads converge on.

**Artifacts (update 2)**: `saddle.mac` (Maxima saddle + charge-stabilization + tuning curve);
`radial_cascade.c` (cascade/density-coupling tests — the slaved-density null + the diagnosis).
