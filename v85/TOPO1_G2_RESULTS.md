# TOPO-1 + G2-0 results — flux quantization gradient & dormant fiber baseline
**Date:** 2026-07-21 · **Runs:** local CPU (N=64, L=15), 6× ring T=120 + 2× ball T=60
**Protocol:** exact fixed-Q ring states via `eta_qflow` (wind=1, wtor=3.5,
profile ω=1.4500, 12k iters, maxF ≤ 2.5e-6); gauged evolution (complex_gauge=1,
bc_type=0); flux via `v73/analysis/magneton` (frame t=120, hole ρ<2).
Data: `/space/scp/v85/topo1/` (seeds, out, logs). Q600 seed reproduces the v73
ring to all printed digits (ω=1.417101, ρ_pk=3.704).

## TOPO-1 — trapped flux vs stiffness (Q) and screening (g)

**Ring branch bottom found:** Q=300 unwinds during relaxation (L_z/Q → 0.002);
the wind=1 branch bottom lies in (300, 450). Q=450–900 all stable gauged
(Q drift −0.37…−0.42% /120 t.u. at g=0.05; −1.7% at g=0.10; gauss 1e-13 floor;
flux steady t=60→120, slightly growing, matching v73 phenomenology).

**Measured (t=120, flux through hole ρ<2, units of g):**

| axis | point | flux | fraction of 2π/g |
|---|---|---|---|
| Q (g=0.05) | 450 | −0.3134 | 0.249% |
| | 600 | −0.3621 | 0.288% |
| | 760 | −0.3885 | 0.309% |
| | 900 | −0.3976 | 0.316% |
| g (Q=600) | 0.02 | −0.1485 | **0.047%** |
| | 0.05 | −0.3621 | **0.288%** |
| | 0.10 | −0.4767 | **0.759%** |

**Finding 1 — stiffness axis is NULL: geometry, not screening.** B_z(0)
ratios across Q (0.871 / 1 / 1.058 / 1.069) match the pure Biot–Savart
current-loop prediction I/R (0.866 / 1 / 1.062 / 1.072) to ≤0.5%. All flux
growth along the ring branch is classical loop geometry; the density lever is
too weak (Σ|Φ|² varies only ×1.10 across the whole branch — the κ=50 knee
caps interior density).

**Finding 2 — g axis CONFIRMS the screening law.** Fraction-of-quantum
scales as **g^2.00** on the 0.02→0.05 leg (measured ×6.10 vs g² prediction
×6.25 — 2% agreement) and g^1.4 on the 0.05→0.10 leg (entering the
nonlinear/degraded regime; the g=0.10 ring sheds 1.7%). The fabric is a
**weak superconductor**: London screening (m_A = g√(2Σ|Φ_a|²); λ_L = 34 →
14 → 7 across the g scan) is real and follows prediction. The
boundary-backwards mechanism (SLIP_NOTE §3.3) is CONFIRMED in sign and law.

**Finding 3 — the ceiling.** Extrapolating the measured g² law, full flux
quantization (fraction → 1) needs g ≈ 0.93 — **9× beyond the
no-stable-matter bound (g ≈ 0.10–0.11)**. Within the V(s) = (μ/2)s/(1+κs)
family, topological flux quantization is parametrically capped at **< 1%**
of the fluxoid quantum. The v73 "not a flux-quantizing vortex" result is now
explained quantitatively, not just observed.

**Finding 4 — the option-C design target (quantitative).** Quantization
needs m_A·a ≳ 1 ⇔ **Σ|Φ_a|²·a² ≳ 1/(2g²) ≈ 50** at g=0.10. Measured today:
Σ|Φ|²·a² ≈ 7–10 (Σ ≈ 1.1, tube a ≈ 2.5–3). **Deficit is only ×5–7** — a
stiff-modulus (vev-type) potential with |Φ|²_vev ≈ 1 and tube radius a ≈ 7
crosses the threshold. Option C is *feasible with a modest target*, but it
is a kernel-level potential change → explicit user authorization required;
design-only until granted.

## G2-0 — dormant gauge2 fiber baseline

Symmetric ball (Q=118, ω=1.45), g_gauge=0.05, with `g_gauge2=0.05`,
charges (1, 1, −2), vs baseline:

| | Q drift /60 t.u. | E | gauss | s_max |
|---|---|---|---|---|
| baseline | −0.026% | 178.6→178.5 | 7.2e-14 | 0.0485 |
| gauge2 on | −0.018% | 178.6→178.6 | 7.2e-14 | 0.0489 |

**Finding 5 — gauge2 is benign AND structurally elegant for symmetric
states:** with charges (1,1,−2) a symmetric ball is **gauge2-neutral**
(1+1−2 = 0) — the ball is a *color singlet* under the second fiber, exactly
the baryon-neutrality structure. Only the +0.8% s_max rise (mild
self-compression from fluctuation coupling) distinguishes it. Consequence:
commensurate-harmonic physics activates only for **flavored/asymmetric**
states — the next rung (G2-1) needs `gen_qball_flavored` seeds with unequal
Q_a, not symmetric balls.

## Verdicts

| Claim under test | Verdict |
|---|---|
| Flux quantization improves with interior stiffness (Q axis) | **NULL** — geometry explains everything; density lever too weak in this potential |
| Screening fraction ∝ g² (weak-superconductor law) | **CONFIRMED** (2% agreement on the clean leg) |
| Full quantization reachable in config space | **REFUTED** — capped <1%; needs g≈0.93 vs bound ≈0.11 |
| Option C (stiff-modulus potential) viability | **Sharply motivated**: ×5–7 deficit, concrete target Σ|Φ|²a² ≥ 50 |
| gauge2 fiber safe to activate | **CONFIRMED**; symmetric states are fiber-neutral (color-singlet analog) |

## G2-1 — equal-profile "flavored" seed (confounded, but two findings)

Seed: same profile all components, ω = (1.43, 1.45, 1.47). **Confound:** this
is off the flavored branch — the components **phase-locked** almost
immediately: expected relative phase 137.5° after 120 t.u., measured 7.4°
(Δω collapsed 0.02 → 0.001, ×20). Intra-object 1:1:1 Adler locking at full
overlap is strong — and explains *why* the v71 flavored branch requires
distinct amplitude profiles to hold distinct clocks (frequency partition
must be paid for in amplitude structure).
**Clean signal despite confound:** the gauge2 run differs from baseline in
exactly one number — the fiber-charge −2 component's relative phase, shifted
−0.75°/120 t.u. (δω ≈ 1.1e-4); the charge-neutral pair (Δq=0) is identical
to baseline. First dynamical signature of the second fiber.

## G2-2 — TRUE flavored-branch colored ball (the real rung)

Seed: v71 flavored profile (4-col), ω = (1.38, 1.42, 1.42), distinct
amplitudes (f₀=0.591, f₁=f₂=0.664); net fiber charge ρ₀−ρ₁ ≠ 0 (colored).
Both runs stable (Q drift −0.08%, gauss 7.9e-14).

**Finding 6 — the flavored branch holds its clocks.** Baseline: a10
accumulates −88.3° at t=120 vs −85.0° predicted from Δω = 0.04 —
**98.8% of the seeded detuning retained** (Δω_eff = 0.0395). Degenerate
pair a21 = 0.00 exactly. Dynamical confirmation of the v71 flavored branch,
and the counterpoint to G2-1's sync collapse.

**Finding 7 — the fiber splits degenerate clocks by integer charge
(color-Zeeman).** Components 1 and 2 are *identical fields* carrying fiber
charges +1 and −2. Under gauge2 they split **linearly**: a21 = −5.71°
(t=60) → −11.39° (t=120) ⇒ δω₂₁ = −1.66e-3. Meanwhile a10 (Δq = 0) is
unchanged vs baseline to 0.4°. Pattern = δω_a = q_a·g₂A₂⁰(core), a
self-sourced color-Coulomb clock shift with g₂A₂⁰ ≈ 5.5e-4:
**level splitting proportional to the integer fiber charge.**

**Why Finding 7 matters for the Q-degeneracy (§D of TOE_v85.3):** since
E = ωQ, a clock shift is a mass/energy shift — the fiber gives components
**different effective masses according to their integer harmonic number**,
with the continuous parameters (g₂, overlap) only setting the scale. This is
the first mechanism in the theory that splits mass by a *quantum number*
rather than by a continuous knob — a crack in charge≡mass degeneracy,
config-level, measured today.

## G2-3 — splitting spectroscopy (run 2026-07-21)

Same flavored-branch seed; a21 at t=120 (baseline = 0.00°):

| g₂ | charges | a21 (t=120) | ratio vs prev | g₂² predicts |
|---|---|---|---|---|
| 0.025 | (1,1,−2) | −2.94° | — | — |
| 0.050 | (1,1,−2) | −11.39° | ×3.87 | ×4.00 |
| 0.100 | (1,1,−2) | −40.38° | ×3.55 | ×4.00 |
| 0.050 | (1,−2,1) | **+11.39°** | exact sign flip | exact |

**Finding 8 — color-Zeeman law confirmed:** δω ∝ g₂² to 3% on the clean leg
(mild saturation at g₂=0.10, same softening as TOPO-1's g-axis), and the
charge-permutation symmetry is EXACT to display precision (+11.39 vs −11.39;
a10/a20 swap values perfectly). The splitting is calibrated spectroscopy of
the second fiber: δω_a = q_a·g₂A₂⁰ with the integer q_a doing exactly what
integers should.

## G2-4 / G2-5 — color force from pair trajectories (run 2026-07-21)

**G2-4 (co-phased, D=7): probe mis-design.** The coherent matter attraction
(v70 cos Δφ force) pulled both pairs into contact oscillation (D swinging
2.8–6.6); fiber differentials ±0.02 = noise. Not a null — a bad probe.

**G2-5 redesign:** quadrature phases (kills coherent force), g_gauge=0.005
(demotes ordinary Coulomb ×100 → fiber is the dominant long-range force),
same on/off × {AB, AA} matrix, T=150.

**Tool bug found & fixed (affects nothing in the repo):** my frame analyzers
mapped SFA_VELOCITY cmp 3–5 (θ velocities) as φ-imaginary velocities
(actually cmp 6–8), yielding ρ ~ ωf²sin²(phase). All trajectories and
per-component charges re-extracted with the corrected tool; corrected
D(0) = 7.10 ✓ and smooth dynamics.

**Finding 9 — color charge measured, conserved, and partition-exact.**
The flavored ball carries (Q₀, Q₁, Q₂) = (76.66, 103.70, 103.70) →
**C = Q₀+Q₁−2Q₂ = −27.04**, conserved to 0.3%/120 t.u. (drift proportional
to total radiation loss). Under gauge2 the degenerate pair sheds
*differentially by fiber charge* (Q₁−Q₂ opens to 0.02%) — the Zeeman-split
clocks radiate differently, as they should. Permuted ball: C_B = +54.1.

**Finding 10 — the long-range color field exists by construction.** The
kernel enforces the gauge2 Gauss law exactly (divE₂ = g₂ρ₂,
scp_sim.c:2404): C ≠ 0 conserved ⇒ ∮E₂·dA = g₂C at every enclosing radius
⇒ the color-Coulomb tail is present at all r, with the coupling
independently calibrated by the G2-3 Zeeman splitting. Static color
electrodynamics is established without a trajectory measurement.

**Trajectory verdict: AMBIGUOUS (first non-positive result of the series).**
dAB(on−off) = −0.14 at t=150 — right sign and order for color attraction
(Coulomb estimate ≈ −0.07) — but the same-color control dAA = −0.15 where
Coulomb predicts +0.03: a **color-independent contraction systematic**
(consistent with the G2-0 +0.8% gauge2 self-compression slowing tail-driven
expansion) is as large as the signal. The dynamical force is *not cleanly
measured*; the field's existence (Finding 10) does not substitute for it.

**Clean probes identified (not run):** (i) static-hold force balance at
fixed D (needs a pinning instrument); (ii) orbit/deflection scattering with
impact parameter (larger box, GPU-class); (iii) exporting E₂/th₂ columns to
SFA — **kernel change, requires explicit authorization** — which would make
the color field directly measurable and G2-6 a one-frame analysis.

## Series verdict (stopping point per protocol)
Positives: Findings 5–10 (fiber benign; color-singlet structure; flavored
branch holds clocks; color-Zeeman g₂²-exact; C conserved & partition-exact;
static color field by Gauss). First ambiguous: the dynamical pair force.
Per the "continue while positive" protocol, the series pauses here with the
three clean probe options above on the table.
- **Option C design doc** (paper): vev-type stiff potential hitting
  Σ|Φ|²a² ≥ 50; changes the whole branch structure — needs full re-map;
  kernel authorization gate.
- κ-scan is NOT worth running: Σ ∝ κ^(−1/3) — even κ=1 gives only ×3.7 of
  the needed ×5–7 with the entire phenomenology re-tuned; option C is the
  honest route.
