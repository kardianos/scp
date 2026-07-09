# v73 — Process-Form Stability: the Particle as Fabric Throughput

**Date**: 2026-07-08
**Question (user)**: combine the stability criteria with "fabric as particle" —
recast stability so it involves the fabric itself, its uptake and layment as the
particle moves. Particles MUST move (T > 0); they maintain themselves in the
fabric through that movement. Prove it, incorporate it, then derive a different
operating mechanic for a particle.

## 0. The claim being formalized

The classical criteria (CONCEPT §2–3: charge pressure vs Derrick, closed
radiation channels, m_θ > ω, VK, exact stationarity) are *state* criteria —
statements about an energy landscape over field configurations. This document
recasts them as *process* criteria: a particle is not a lump of fabric but a
**closed process running ON the fabric** — a pattern that continuously takes up
fabric into motion and lays it back down, and persists exactly when that ledger
closes. The fabric never travels with the particle; the process does.

## 1. The ledger variables (exact, from the field equations)

Local energy density e(x,t) and charge density ρ_Q(x,t) obey exact continuity
laws. For the complex η-coupled system (ungauged), with all 12 real fields f:

    ∂e/∂t + ∇·S = 0,     S = − Σ_f ḟ ∇f  +  η [ u̇⃗ × θ⃗_u + v̇⃗ × θ⃗_v ]
    ∂ρ_Q/∂t + ∇·J = 0,   J = Σ_a [ v_a∇u_a − u_a∇v_a ] + (Θ terms)

(the η cross-flux term follows from d/dt of the coupling energy −η θ·(∇×φ)
and the identity θ·(∇×φ̇) = ∇·(φ̇×θ) + φ̇·(∇×θ); the component index acts as
the spatial index, matching the kernel's curl convention. Implemented in
`analysis/flux_profile.c`.)

S is the **fabric ledger**: S·da through a surface is the rate at which the
process takes up (inward) or lays down (outward) fabric-motion across it.

## 2. The process-form criteria

**P1 — Motion is constitutive (T > 0 by definition).** A particle cannot be a
static arrangement of fabric: E(λ) = aλ + bλ³ has no minimum (Derrick, CONCEPT
§2) — *a motionless lump is dead*. Existence requires a nonzero internal
velocity field at every instant. For the U(1) ball the motion is the internal
clock: every fabric point runs a circle of radius f(r) in the (Re,Im) plane at
rate ω. The particle at rest is not "static with a phase label" — it is pure
motion whose *envelope* is static.

*Quantified*: the motion content is exactly the charge times the clock rate,

    2 E_kin = ω Q        (measured: 86.609 vs ωQ/2 = 86.631 → 0.03% at η=0.25;
                          89.829 vs 89.841 → 0.013% at η=0.5)

so the action the process runs through per clock period is

    ∮ 2E_kin dt = (ωQ)(2π/ω) = 2πQ = h_eff   — exactly one quantum of action
    per re-constitution cycle.

The particle *is* a periodic re-constitution of itself, and h_eff = 2πℏ_eff is
the fabric-motion cost of one cycle. This is the process reading of ℏ_eff = Q
(CONCEPT §7): quantization is metabolic — one action quantum per heartbeat.

**P2 — The ledger must close (uptake = layment).** Persistence ⟺ the
time-averaged flux through EVERY closed surface vanishes:

    ⟨S·r̂⟩_T (r) = 0  and  ⟨J·r̂⟩_T (r) = 0   for all r,  T = one clock period.

Fabric may be taken up and laid down locally (the Θ dressing is fabric ON LOAN
— driven, evanescent, returned each cycle because m_θ > ω keeps it within
reclaim range 1/√(m_θ²−ω²)), but no fabric may be laid beyond reclaim.
*Radiation is precisely unclosed layment*: the v65 3ω harmonic was fabric laid
down at a frequency the process cannot re-uptake (above the mass gap →
propagates); the June-26 "η-drain" was fabric laid into an inconsistent Θ
channel. The classical criteria map onto P2 as follows:

    time-independent |Φ|  (no harmonics)   → nothing laid above the gap
    m_θ > ω  (kinematic closure)           → the Θ loan stays within reclaim
    exact stationarity (v72 fixed-Q state) → the ledger closes to the floor

**P3 — Ledger feedback must be self-limiting (VK as supply law).** ω is the
processing RATE (re-constitutions per unit time); Q is the standing INVENTORY
of fabric held by the process. The branch relation ω(Q) is the process's
supply curve, and stability requires

    dQ/dω < 0   ⟺   dω/dQ < 0 :  adding inventory must SLOW the clock.

If extra fabric sped the clock up, a fluctuation that grabbed fabric would
process faster, grab more, and run away. VK stability is a negative-feedback
condition on the fabric ledger — measured on the η-branch (v72 §4/§5): ω falls
monotonically from 1.4679 (Q=88) to 1.4167 (Q=210), and the branch *ends*
(Q_min ∈ (80,84)) where the process can no longer hold inventory at any rate
(ω → m: the clock reaches the free-fabric band and the inventory delocalizes).

**Movement (v > 0) is the same process, conducted.** A translating particle
does not carry fabric: at the leading surface vacuum fabric is taken up into
rotation, at the trailing surface it is laid back to vacuum. For a rigidly
translating pattern every conserved density ρ must satisfy J⃗ = v⃗ρ (+ solenoidal
part), i.e. the fabric conducts the process at exactly the pattern speed:

    ⟨S_x⟩ = v · e(x),   ⟨J_x⟩ = v · ρ_Q(x)      (prediction, tested in §3)

and each fabric element's history is a CLOSED loop: from vacuum, through a
transient spin-up (uptake), a de Broglie-tilted rotation while inside, spin-down
(layment), back to vacuum. Nothing material is transported; the imperfection of
layment on a lattice (the wake) is the measured 1–5% group-velocity lag
(v70). de Broglie's phase tilt k = γωv is, in process terms, the phase LAG of
re-constitution along the direction of motion — downstream fabric is re-made
later, so the clock phase advances along −x̂·v: motion and the matter wave are
the same bookkeeping.

## 3. Measurements (all done 2026-07-08; runs fx_stat / fx_theta0 / fx_boost,
N=64, L=15, snap_dt=0.25, ungauged; tool `analysis/flux_profile.c`)

**E4 — throughput identity (P1): CONFIRMED.**
2E_kin = ωQ to 0.03% (η=0.25: 86.609 vs 86.631) and 0.013% (η=0.5: 89.829 vs
89.841). The particle runs through exactly h_eff = 2πQ of action per clock
period.

**E3 — ledger closure (P2): CONFIRMED, quantitatively.** Time-averaged outward
power P(r) = 4πr²⟨S_r⟩ over t ∈ [5,20] (`results/flux_stat.tsv`,
`flux_theta0.tsv`):

    r        stationary state      Θ=0 drain seed
    2.1      +4.5e-4               −3.6e-2   (inward: building the Θ loan)
    5.0      +4.3e-4               +5.9e-2   (outward wind begins)
    7.9      −2.2e-3               +1.27e-1  (peak shedding)
    10.7     −1.1e-4               +4.0e-2
    12.1     +6.8e-4               +1.0e-2   (entering the sponge)

The stationary state's ledger closes at every radius (|P| ≲ 2e-3, alternating
sign — zero-mean churn); the drain seed shows a sustained one-signed outward
wind, ~50–100× larger, whose far-shell value (≈0.04 at r ≈ 11) matches the
T=60 secular loss rate of the same seed (−0.047/t.u.). Radiation IS unclosed
layment, and the ledger accounts for it quantitatively. The drain profile also
resolves the mechanism: simultaneous INWARD flux at r < 4 (the ball taking up
fabric to build the Θ dressing it was seeded without) and outward shedding of
the excess — uptake and layment, both caught in the act.

**E2 — pass-through (conduction): CONFIRMED.** v=0.3 ball, fixed probes at
x = 3, 6, 9 on the flight path (`results/flux_boost.tsv.probe*.tsv`):

- The SAME waveform passes each probe: |Φ|²_peak = 1.1145 / 1.1152 / 1.1174
  (0.3% spread); ρ_Q peak 1.691 vs the boosted-ball prediction γωΣf₀² = 1.707
  (1%). The fabric at each point is taken up into exactly the ball's rotation
  state.
- Measured transit speed 0.308–0.316 vs seed 0.3 — the few-% excess is the
  known lattice-dispersion group-velocity anomaly (v70).
- **Layment completeness 99.75%**: probe-3 fabric swings to FULL ball
  amplitude during passage (max |u| = 0.61) and returns to 1.8e-2 after —
  a closed excursion. The ball ends 7.8 units downstream while the fabric it
  was "made of" at passage sits back at vacuum. The medium is not transported;
  the process is conducted through it. The 0.25% unreturned residue is the
  wake — the lattice's imperfect layment.

## 4. The different operating mechanic (design)

The U(1) ball's recirculation is entirely INTERNAL — every fabric point loops in
the (Re,Im) plane; the spatial charge current is identically zero at rest
(J = Im[Φ̄∇Φ] = 0 for real profile × e^{iωt}). The process-form view suggests a
genuinely different mechanic: **a particle whose recirculation is in REAL
space** — fabric taken up at one part of the structure, conducted through it,
laid down at another, in a closed spatial circuit. The minimal such object in
this theory is the **winding state (Q-ring)**:

    Φ_a = f(ρ,z) · e^{i(ωt + nφ)}        (all three components, same winding)

- s = Π|Φ_a|² is phase-blind → binding survives; f must vanish on the axis
  (regularity) → the object is a TORUS.
- The charge current is a real circulating flow: J_φ = n|f|²/ρ ≠ 0 — fabric
  literally circulates around the ring. Uptake and layment happen continuously
  around the circuit: this is the "operating mechanic" where the particle is a
  self-maintained fabric CIRCUIT, not a breathing point.
- It carries intrinsic angular momentum L_z = nQ exactly (for this ansatz
  L_z = ∫ω n|f|² ... = n·Q), so J/ℏ_eff = n: **spin arises as the winding
  number of the real-space circulation** — the first mechanism in this theory
  that could give particles quantized intrinsic angular momentum.
- Process criteria apply verbatim: P1 (the circuit is motion), P2 (the loop
  must close: no leakage from the ring edges), P3 (dQ/dω < 0 along the ring
  branch). Known risk: in the Q-ball literature spinning solitons exist as
  stationary solutions but often sit on higher-energy branches; a full 3D
  fixed-Q flow may unwind through the axis.

### 4.1 RESULT — the Q-ring EXISTS and is kernel-stable (2026-07-08)

**Negative first (thin ring):** a torus seed whose amplitude sits below the
binding knee (κs ≈ 0.08) does not bind — the fixed-Q flow (from two opposite
starting points, ω₀ = 3.17 and 1.24, converging to the same state — good
numerics) runs to the delocalized free-condensate branch (ω → 1.5155 > m,
s_max → 4e-8). A ring must be seeded FAT: `ring_thin_negative.log`.

**The smoke-ring seed works** (`eta_qflow` winding mode: tube cross-section =
the full ball profile around major radius R₀ = 4, axis regularizer, n = 1,
Q auto = 679.2, η = 0):

    relaxed:  ω = 1.411378, maxF → 1e-9, s_max = 0.0461 (fully bound core),
              ρ_peak = 3.77, winding@peak = 0.98, L_z = 675.5, L_z/Q = 0.9946

The winding SURVIVED the unconstrained 3D flow (no unwinding through the axis).

**Kernel verification (T=60, ungauged):** drift −0.072% — exactly the
absorbing-BC floor (identical to the ball's); s_max flat (0.0461 → 0.0461);
omega_core = 1.41132–1.41144 vs the relaxer's 1.411378 (5 decimals). Frame
analysis (`analysis/ring_check.c`) at t = 0 vs t = 60: L_z/Q = 0.9946 → 0.9946,
winding 0.979 → 0.979, J_φ(ρ_pk) = −0.287 → −0.287 — every circulation
invariant frozen to the printed digits. (Offsets from exactly 1 are azimuthal
lattice discretization; the identity is L_z = nQ.)

**The ring branch obeys the ledger law (P3):** Q = 600 → ω = 1.41710,
Q = 679 → ω = 1.41138, Q = 760 → ω = 1.40653 — dQ/dω < 0 across the branch,
with L_z/Q pinned at n (0.9943 / 0.9946 / 0.9948) and the major radius
swelling with inventory (ρ_pk = 3.70 → 3.77 → 3.94). Circulation lets the
process hold ~2.5× the ball-branch charge at the same clock rate (the ball's
branch holds Q ≈ 250–310 over ω = 1.406–1.417).

**What this is:** a stable particle whose self-maintenance is REAL-SPACE fabric
circulation — charge current J_φ ≠ 0 flowing around a torus — carrying
intrinsic angular momentum locked to its charge, L_z = nQ = n·ℏ_eff, i.e.
**J/ℏ_eff = n: spin as the winding number of the circulation.** This is the
first spin-carrying object in the theory, and the first whose operating
mechanic differs qualitatively from the Q-ball (internal-plane rotation →
spatial circuit). Open questions: the n = 2 branch (J = 2ℏ_eff, and whether
fatter/multi-wound rings fission), the gauged ring (a circulating charge is a
current loop → it should carry a MAGNETIC MOMENT through the U(1) field — the
magneton analog), the η-dressed ring, and ring–ring / ring–ball interactions.

### 4.2 The second dimension — the poloidal twist (n,m) ring (user direction)

The ring circulates around only one of the torus's two circles. Adding the
poloidal winding (phase = ωt + nφ + mχ; hollow tube, zero-set = z-axis ∪
tube-core circle = a Hopf link; circulation lines = (n,m) torus knots) was
tested at (1,1):

- **Not a fixed-Q minimum.** The flow unwinds m at both geometries tried:
  R₀ = 4 lost BOTH windings (tube core filled → ring hole closed → ball;
  `tring_11_unwound.log`); R₀ = 6 with a boosted profile kept n (a heavier
  plain ring: Q = 1233, L_z/Q = 0.996) but still released m
  (`tring_fat_n-survived.log`). The twist has no conserved charge of its own —
  only an energy barrier — and at these couplings the barrier does not hold
  against unconstrained descent.
- **Dynamically METASTABLE.** Kernel dynamics conserves E (the flow does not),
  so unwinding requires radiating the twist energy. Seeded raw (far off-shell,
  sheds 30% of E over T=40), the twist nevertheless holds 1.000 → 1.000 →
  0.967 over t = 0 → 10 → 20 before the churning geometry makes the probe
  ambiguous (`tring_k_diag.tsv`, ring_check frames). Meanwhile **n and
  L_z/Q ≈ 0.99 survive the entire violent evolution** — spin is a good quantum
  number even 30% off-shell.
- **Fair trial still open**: it needs a twist-preserving constrained
  relaxation (fixed Q plus a helicity/twist constraint, or a 2D reduced flow
  where m lives as a protected vortex of the cross-section field F(ρ,z)) —
  the concrete next build if the (n,m) taxonomy is pursued.

### 4.3 The third dimension tried — fabric swirl (η-dressed ring): FISSION

The next candidate dimension was a pattern that is NOT a phase winding: the
poloidal circulation of the Θ fabric around the tube (the literal smoke-ring
swirl), sourced by the η curl coupling. Controls first: the η-dressed ball's
Θ cloud has **identically zero swirl helicity** (H = ∫Θ⃗·(∇×Θ⃗) ≡ 0 — its
azimuthal pattern is pointwise perpendicular to its own curl), so any helicity
or coherent poloidal circulation appearing on a dressed ring would be a
genuinely new dimension. Tool: `v73/analysis/swirl_check.c` (helicity,
poloidal/toroidal Θ loop circulations at 8 azimuths, E_int, field momentum).

**Result — the fixed-Q flow does not dress the ring; it fissions it.**
At fixed Q = 679.2151 (matched to the plain ring) the flow migrated, at both
η = 0.25 and η = 0.5, to a **standing azimuthal wave**: mode decomposition
(`v73/analysis/azi_profile.c`) shows A_v/A_u = 1.04 with phase difference
180° — v = −u pointwise, arg Φ constant ≈ −46° around the ring — i.e. two
opposite-phase dressed lumps at antipodes (ρ ≈ 6.7), J_φ ≈ 2e-7, L_z/Q ≈
0.0000. (The seed reporter's "winding@peak = 1.00" on such states is an
artifact: the standing wave's two π phase-jumps at its nodes sum to 2π.)

The energy ledger says why: E(plain ring) = 1003.26; E(flow final, η=0.25) =
982.56; two FREE dressed balls at Q/2 = 339.61 each cost 2 × 491.20 = 982.41.
The final state is **unbound by +0.15 and still descending** (maxF plateau
1e-3 = the separation soft mode). Kernel verification of the two-lobe state:
drift −1.62% over T = 60 (22× the stationary floor), lobes separating
(ρ_pk 6.74 → 7.07) — fission in slow motion, not a particle. Even at η = 0
the ring is only LOCALLY stable: it sits ~18 units above its own fission
products, protected by angular momentum alone.

Reading: the fabric dressing rewards density (s_max rises 52% at fixed Q when
η = 0.25 is switched on) but at fixed Q alone the flow cashes that reward by
killing the circulation through the standing-wave channel — densification
eats velocity unless something ties them together. What ties them together
in true dynamics is angular-momentum conservation, which the fixed-Q flow
does not know about — and which the η coupling itself weakens (§5).

**No threshold.** Repeating the flow at η = 0.1 and η = 0.05 gives the same
endpoint, only slower: both reach L_z = 0.0000, winding 0.00, E at the
two-free-balls plateau for their coupling (984.6 / 985.0). Any η > 0 opens
the channel in the flow; only η = 0 holds the ring — and it holds it for an
exact reason, not a marginal one (§5.3: equivariance).

## 5. Working backwards: the fabric-density forming function
##    (user direction, 2026-07-08)

The proposed mechanism: one (or a combination) of the added dimensions makes
the fabric locally tighter — denser — which drives the particle tighter, which
raises its internal velocity, until the two effects balance at a natural
equilibrium. The requirement: that equilibrium must be a TRUE minimum in
every direction, not a saddle (the twist of §4.2 failed exactly as a saddle).
The method: do not guess shapes and relax them forward; start from the
stability equation and derive what density-forming coupling makes the
compressed state a minimum.

### 5.1 The forming function the theory already has

Integrating out the fabric at linear response (Θ ≈ η (m_θ² − ω² − ∇²)⁻¹ ∇×Φ)
turns the η coupling into an effective energy

    E_dress ≈ − (η²/2) ⟨ |∇×Φ|² ⟩_screened ,

a density-forming function whose density variable is the CURL of the matter
field — circulation density, not amplitude. This is the right shape for the
mechanism: fabric uptake rewards motion (§2's constitutive-motion criterion
P1), and it measurably compresses — at fixed Q = 679.2, switching on η = 0.25
raises the core density s_max by 52% (0.046 → 0.070). The compression→velocity
side needs a second ingredient: at fixed charge AND angular momentum J = nQ,
shrinking the ring radius forces the azimuthal circulation to speed up
(J fixed, R smaller ⇒ v_φ larger), and the n²/R² gradient pressure is the
back-reaction that halts collapse. Densification converts to velocity ONLY
if J is held.

### 5.2 The stability equation with the second charge

The fixed-Q principle (v72) generalizes. Rigidly spinning ansatz
Φ(x,t) = R(Ωt) φ(R(−Ωt)x) e^{iωt} (R rotates space AND the component index —
the Cosserat vector structure demands both), two conserved quantities, two
multipliers:

    Q = ωN − Ωℓ ,   J = ωℓ − ΩK ,
    N = ∫Σ|φ|² ,  ℓ = 2∫ uᵀĴv (per complex pair) ,  K = ∫|Ĵφ|² ,
    Ĵ = ∂_φ − ẑ×   (orbital + spin) ,

with kinetic energy E_kin = (ωQ − ΩJ)/2 and flow force gaining the Coriolis
and centrifugal terms  F_u += ω²u − 2ωΩ Ĵv − Ω² Ĵ²u  (and v ↔ −u). Fixing
(Q, J) closes the standing-wave channel BY CONSTRUCTION — a rotating→standing
deformation changes J and is no longer admissible descent. Stability becomes:
constrained Hessian ≥ 0 plus the 2×2 Vakhitov–Kolokolov matrix criterion on
∂(Q,J)/∂(ω,Ω). On an exact circulation eigenstate (Ĵφ = inφ) the pair (ω,Ω)
is gauge-degenerate with only ω_eff = ω − nΩ = Q/N physical, so the solve
needs Tikhonov regularization there; the η-dressed states of interest are
never exact eigenstates (the dressing pattern must co-rotate) and are
well-conditioned.

### 5.3 The obstruction: the fabric coupling breaks rotation symmetry

Symmetry audit of the kernel Lagrangian:

- V(s), s = Π_a|Φ_a|²  — invariant under SPATIAL rotation; NOT invariant
  under internal (component) rotation.
- η Θ̄·(∇×Φ)  — invariant under JOINT spatial+component rotation; NOT
  invariant under spatial-only rotation.

At η = 0 spatial rotation is exact ⇒ orbital L_z exactly conserved — this is
why the plain ring is kernel-stable and why even the shredding twisted ring
held L_z/Q ≈ 0.99. At η > 0 NO continuous rotation survives ⇒ no exact
angular-momentum conservation at all. One protection remains: on
equal-component states (|Φ_1| = |Φ_2| = |Φ_3| pointwise — the standard seeds)
the potential's torque vanishes to first order (δs ≡ 0 under joint rotation
when the component magnitudes are equal), so the leak is second-order small.

Symmetry also explains why η = 0 tests could never show the fission even if
energetics favor it: gradient flow and kernel dynamics are EQUIVARIANT — a
configuration invariant under (rotation by α × phase shift −nα) has an
invariant time derivative, so from the symmetric ring seed the asymmetric
standing mode simply cannot be generated at η = 0 (only lattice anisotropy
and float noise seed it, far below visibility). At η > 0 the coupling itself
manufactures the asymmetry: the Θ pattern sourced by ∇×Φ cannot co-rotate
within the broken symmetry, and the fission mode is actively driven.

Measured (ring_eta_k: plain rotating ring dropped into η = 0.25 dynamics,
T = 120): L_z/Q holds 0.995 ± 0.001 through t ≈ 80, easing to 0.984 by
t = 120 — spin-specific leak ≈ 1% per 120 t.u. But the fission channel is
dynamically OPEN: the counter-rotating (standing) admixture grows 1.4% (t=60)
→ 12% (t=120), an e-folding time of ~30 t.u. ≈ 25–50 internal clock cycles,
with the ring visibly swelling (ρ_pk 3.77 → 4.38) and the peak current down
27%. The η = 0 control run over the same geometry is a pure rotating state
to five decimals at t = 60 (A_v/A_u = 1.0000, Δphase = 90.0°). And the flow
finds no threshold: η = 0.05 and 0.1 fission completely too (§4.3).

**Conclusion: in the current kernel, fabric densification via η and spin
protection are in structural tension.** The η-dressed spinning ring is a
resonance (lifetime ~10² t.u. at η = 0.25, longer at smaller η as the
breaking is O(η²)), not a stable particle.

### 5.4 The symmetry-consistent densifier: the gauge sector

Working backwards from the requirement — a forming function that rewards
circulation density WITHOUT breaking rotation symmetry — the kernel already
contains exactly one: the U(1) gauge sector. |（∂−igA)Φ|² + F² is spatially
isotropic and touches no component index: orbital angular momentum is exactly
conserved at any g. The gauged ring is a current loop; its magnetic flux
through the hole exerts hoop stress (the vorton mechanism) — an energy
barrier against BOTH failure channels: collapsing the hole compresses the
trapped flux (E_B ~ Φ_B²/R rises), and killing the current (the standing-wave
channel) must first destroy the flux it sustains, which costs the full
field energy and is impeded by the ring's own conductivity. Gauge dressing
is the fabric-density mechanism that ties density to velocity through an
EXACT conservation law.

**Measured (gring_k: ring seed + complex_gauge = 1, g = 0.05, η = 0,
T = 120):** the gauged spinning ring is STABLE.

    L_z/Q          0.9946 → 0.9925  (−0.2% over 120 t.u., no acceleration)
    A_v/A_u (t=120)  1.0004          (zero fission-mode growth; η=0.25
                                      reached 1.28 at the same time)
    winding          0.981 constant;  J_φ steady at −0.26
    geometry         ρ_pk breathing 3.94–4.06, no secular swelling
    E_em             ≈ 4.7 sustained (~0.5% of M) — a permanent EM sector
    gauss_max        7.6e-14 throughout (exactness floor)
    Q drift          −0.26%/120 t.u. (includes the gauge-dressing transient
                     of the ungauged seed; ungauged floor is −0.14%/120)

This is the first object in the theory that simultaneously carries conserved
charge Q, integer spin J = nQ = n·ℏ_eff, and a real electromagnetic moment
(a circulating current loop with sustained field energy). The user's
densification mechanism closes here: the gauge field is fabric density that
COUPLES to circulation (it exists only because the current does), tightens
the state (hoop stress balances gradient pressure), and is protected by an
exact symmetry — a true dimensional minimum rather than a saddle, pending
the exact-state relaxation below.

**Measured (magneton, 2026-07-09; tool `analysis/magneton.c`, frames of
gring_k):** the gauged ring is a genuine magnet. It magnetizes itself from
the B = 0 seed within ~60 t.u. and holds steady: B_z(center) −0.0349 →
−0.0352 and trapped hole flux −0.374 → −0.380 (units of g) over t = 60 → 120,
with the textbook current-loop profile (sign flip at ρ ≈ 4.5, the current
radius). The gyromagnetic ratio, μ_z / (Q·L_z/2M), is

    g_factor = 1.052 / 1.059 / 1.060    (t = 0 / 60 / 120)

— a classical extended rotor (g ≈ 1) with a +6% anomaly from the mismatch
between charge density (∝ ω f²) and energy density (gradients + negative
potential concentrated at the tube), and decisively NOT a Dirac-like g = 2.
The EM energy is mostly electric (E_B ≈ 0.09 vs E_E ≈ 4.65; the tool's
E_B + E_E reproduces the kernel's E_em diagnostic to 4 decimals). The
trapped flux is far below the London/vorton quantum 2π/g ≈ 126 — this ring
is a weak current loop, not a flux-quantizing superconducting vortex.
Caveat: the relative SIGN of μ_z vs B depends on the kernel's link-phase
convention (unaudited); magnitudes and stability are convention-free.

**Open (designed, not yet built):** the two-multiplier fixed-(Q,J) flow of
§5.2 — gauged — to produce the EXACT stationary gauged spinning ring (the
kernel test above used the ungauged seed + Gauss projection, so the observed
object is the dressed ring after a small transient), and its 2×2 VK matrix
over (Q, J) to map the spin branch.

## 6. The angular-momentum sector: spin is polarization; the yrast band
##    (fixed-(Q,J) relaxer, 2026-07-09)

The two-multiplier flow of §5.2 is implemented in eta_qflow (fixJ mode:
per-iteration 2×2 solve for (ω,Ω), Coriolis + centrifugal force terms,
light-cylinder mask at ρ=11 with |Ω| ≤ 0.15, Tikhonov-guarded solve;
regression: fixJ=0 is byte-identical to the previous tool). Scanning E(J)
at fixed Q = 679.2151, η = 0, N = 64, L = 15 produced the sector's map.

### 6.1 Spin is internal polarization — and it is FREE (flat band)

The Φ triplet is a complex 3-vector: a POLARIZATION. A ball with a static
relative phase Δ between components x and y,

    Φ = f(r)·(e^{-iΔ/2}, e^{+iΔ/2}, 1)/√3 · e^{iωt},

is an exact stationary kernel solution (each component its own Q-ball,
U(1)³) carrying spin angular momentum

    J_z = ω ℓ_z ,   ℓ_z = (2/3) N sin Δ     (equal magnitudes)

with NO spatial rotation, NO moment of inertia, and NO energy cost — every
energy term is blind to static component phases. Measured (machine
precision): E(J=0) = E(J=Q/4) = E(J=Q/2) = 963.8322 with

    J target      Δ predicted = asin(J/(ω·(2/3)N))    Δ measured
    169.8038      22.02°                              22.02°
    339.6076      48.59°                              48.59°
    (ℓ_z measured 123.0212 / 246.0424 vs predicted 123.02 / 246.05)

**The SCP ball carries spin-1-like polarization structure**: the component
index is its polarization vector; elliptical polarization = intrinsic spin;
the η curl coupling acts as a spin-orbit coupling. (This is a structural
analogy at the level of representation content — no exchange/statistics
claims are made.) Unlike the ring's topological J = nQ (integer winding),
polarization spin is classically continuous — the two angular-momentum
carriers of the theory are fundamentally different objects.

### 6.2 Past saturation: magnitude trading + rigid rotation

The flat band ends at phase saturation Δ = 90°. Beyond it the flow finds,
in order of engagement: (i) magnitude trading — norm moves from component z
into the spinning (x,y) pair, raising spin capacity above (2/3)N at some
binding cost (measured ℓ_z = 361.5 > 327.6 = (2/3)N at J = 0.75Q, with
|Φ_z|²/|Φ_x|² ≈ 0.72); (ii) rigid spatial rotation Ω ≠ 0 of the deformed
envelope. Measured yrast points (Q = 679.2151):

    J/Q     E          Δ        Ω          state
    0       963.832    0°       0          ball (anchor)
    0.25    963.832    22.02°   0          ball + phase spin   (FLAT)
    0.50    963.832    48.59°   0          ball + phase spin   (FLAT)
    0.75    965.475    90.00°   -0.0569    saturated + traded + rotating
    1.00    973.938    90.00°   -0.0568    saturated + traded + rotating
    ---- the ring (isomer) band, same Q ----
    1.00    1003.258   0°       0          ring n=1 (all-orbital) — the
                                           branch point, +29.3 above yrast;
                                           constrained flow reproduces the
                                           unconstrained ring EXACTLY
    1.50    1003.258   49.30°   0          ring + phase spin   (FLAT —
                                           orbital 675.5 + spin 343.3)
    2.00    ~1016.8    88.7°    -0.0132    ring n=2 (winding 1.99 held,
                                           ρ_pk 8.34) + saturated spin;
                                           +13.5 over n=1 for +Q of J
                                           (maxF 3e-3, still descending)

Both bands are spin-broadened by the same free polarization mechanism: the
yrast (ball) band is flat over |J| ≤ ω(2/3)N ≈ 0.67Q, and the ring band is
flat over J ∈ nQ ± ω(2/3)N ≈ [0.33Q, 1.67Q] for n = 1 — overlapping bands
separated by a constant configuration gap ≈ 29.3 (the orbital circulation
cost), exactly the structure of nuclear rotational bands built on different
intrinsic states.

### 6.3 Consequence: J-protection cannot save the ring; the ring is a K-isomer

The η = 0.25 ring relaxed at fixed (Q, J = Q) does NOT stay a ring: the flow
converts orbital circulation into polarization spin (orbital L_z 675 → 217,
spin ℓ_z 0 → 317, phases forming the symmetric ladder 0°/39.5°/79°) and
fissions toward the two-lobe shape anyway (E → 983.5, ~1.1 above the free
fission plateau — the residual spin-carrying cost at η > 0). Spin and
orbital angular momentum are INTERCONVERTIBLE through the component sector,
and spin is the cheap carrier.

So the spinning ring is the theory's K-ISOMER (nuclear-physics sense): an
all-orbital, topologically wound configuration living ~3% of M above the
yrast band at the same quantum numbers, protected NOT by energetics or by J
conservation but by symmetry — exactly at η = 0 (equivariance, §5.3), and
dynamically at g > 0 by the gauge sector (§5.4). Its J = nQ is quantized;
the yrast band's polarization spin is not. If the ring is the interesting
"spin-carrying particle" of the theory, its protection is topological
+gauge, and the η fabric coupling is its decay channel (lifetime ~10² t.u.
at η = 0.25, measured §5.3).

### 6.4 Kernel exactness classification (yrast_k, the honest cut)

Dropping the constrained J = Q yrast state (spin-saturated, magnitude-traded,
Ω = −0.057) into the kernel (T = 60, η = 0) splits it into its exact and
forced parts:

    spin ℓ_z        361.6 → 349.8   (−3%, tracks the norm shed) — HOLDS:
                    trading kept |Φ_x| = |Φ_y|, so ω_x = ω_y and the z-spin
                    phase is locked; the transverse spin (l_x, l_y) precesses
                    freely (component z runs at a different clock rate)
    orbital L_z     159 → 25 → 5    — DIES in ~30 t.u. (radiated), E −3.9%

Because V = Π|Φ_a|² is not invariant under joint (space+component) rotation,
a rigidly rotating state with unequal component magnitudes is NOT an exact
kernel solution — the constrained flow's Ω ≠ 0 band segment is a transient
family, and true dynamics relaxes it back onto the flat band. The exact
stationary states at η = 0 are:

    EXACT:  balls with polarization z-spin (the flat band; equal magnitudes
            in the spinning pair, static phases, Ω = 0),
            rings n = 1, 2, ... each broadened by free polarization spin
            (yr_J150: ring + spin at E = 1003.258 exactly; yr_J200: n = 2
            winding held, E ≈ 1016.8, +13.5 over n = 1 for +Q of J)
    NOT:    the Ω ≠ 0 rigid rotors between the bands (lifetime ~30 t.u.)

The physical division of angular momentum in this theory: MECHANICAL
(orbital) angular momentum is topologically quantized — it exists in stable
form only as ring winding, J_orb = nQ = n·ℏ_eff. POLARIZATION spin is
continuous, internal, costs no energy, and is protected by the symmetry of
the spinning pair — an isospin-like internal quantum that the η coupling
(spin-orbit) converts to mechanical form. The integer-quantized carrier and
the continuous internal carrier are DIFFERENT DEGREES OF FREEDOM — the
theory's own version of the spin/orbital distinction, discovered by
constraining their sum.
