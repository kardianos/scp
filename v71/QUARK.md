# v71 QUARK — Is the Q-ball a quark or a proton? (Component-substructure verification)

**Date**: 2026-06-11. **Status**: [measured] — 5-run V100 campaign (q1–q5, all rc=0,
~25 min GPU total) + per-component analysis + flavored-baryon BVP solver.
**Question (user)**: verify whether the Q-ball is a quark rather than a proton,
and how a proton would be built from quarks.

New tools: `sfa/analysis/sfa_qcomp.c` (per-component Q/mass/centroid/rms),
`sfa/seed/gen_qball_quark.c` (component-displaced "quark" seeds),
`v71/analysis/flavored_qball.py` (3-field Newton BVP, asymmetric baryons),
`sfa_slice`/`render_slices.py` extended with per-component densities
rho2_0/1/2. Data: `/space/scp/v70/q*_*.sfa`, `v71/results/`.

## 0. Verdict

**The Q-ball is not a quark — it is the baryon.** The quark-level objects are
the three internal component fields Φ_a, and they pass every quark-style test
the theory admits:

| property | quark (QCD) | component Φ_a (this theory) | evidence |
|---|---|---|---|
| fractional charge | e/3 units | **exactly Q/3 each** (105.06×3 = 315.17) | sfa_qcomp on banked ball, stable in t |
| no isolated existence | confinement | **single-component lump cannot bind** (force ∝ Π_{b≠a}\|Φ_b\|² ≡ 0 [thm]; disperses as free KG) | q1: s ≡ 0 exactly, rms 3.4→13.2, 60% charge lost by t=100 |
| bound only in color-neutral combos | baryons + mesons | **only all-three combos bind** | q2: two components, s ≡ 0, same dispersal — **NO meson sector** (structural difference from QCD) |
| baryon = 3 bound quarks | proton = uud | symmetric ball = 3 co-located component lumps | q3/q4: assembly (below) |
| flavor structure | u vs d | **asymmetric (ω_0,ω_1,ω_2) stationary baryons exist**, Q_0/Q_1 down to 0.74 at Δω=0.04 | flavored_qball.py continuation (validated: symmetric case reproduces the v66 branch, Q=209.5 ✓) |

Caveats up front: the "color" here is a U(1)³ structure (only the diagonal is
gauged), NOT SU(3); there is no gluon sector and no linear-confinement string —
isolated components disperse rather than pull back with constant tension; and
the component index doubles as a spatial index in the η-curl term (η=0
throughout this campaign). The mapping is structural, not quantitative.

## 1. Free quark test (q1) and meson test (q2) — both die [measured + thm]

Force on Φ_a is −2Vt′(s)Φ_a Π_{b≠a}|Φ_b|²: with any component identically
zero, ALL components are free Klein–Gordon fields [thm, from the v66-verified
EOM]. Measured: a single-component lump of charge 105 (one quark's worth)
keeps s_max = 0 to machine zero for all t, spreads (rms 3.4 → 13.2 by t=100),
and loses 60% of its charge to the boundary; renders show pure dispersal
rings, the binding panel black, the Coulomb halo evaporating. Two co-located
components (q2) behave identically (phi_max 0.635 → 0.06). **Quarks and
"mesons" cannot exist; only three-component ("color-complete") objects bind.**

## 2. Proton assembly from three quarks (q3–q5) — the binding basin [measured]

Three single-component lumps at triangle radius d (each Q_a = 105, ω = 1.42
seeding, N=128, L=20, absorbing BC, full snapshot record):

| d | infall time | captured charge/quark | final object |
|---|---|---|---|
| 2 | ~15 t.u. | 75.3 of 105 (72%) | bound ball, Q ≈ 226, settled (drain → 0) |
| 4 | ~60 t.u. | 48 of 105 (46%) | bound ball, Q ≈ 144, still ringing at t=150 |
| 6 | ~120 t.u. | ~33 of 105 (31%), still falling | small remnant baryon amid dispersal fog |

Mechanism (visible frame by frame in the renders, `v71/results/q4.png`):
binding ignites ONLY in the central region where all three overlap; each
quark's core falls toward it while its outer envelope — not yet captured —
disperses as free waves. Assembly is a race between infall and dispersal:
the farther apart the quarks start, the smaller the baryon that survives.
Per-component charges remain exactly equal throughout every assembly
(e.g. 75.33/75.32/75.32 at q3 end) — capture is democratic.

**This is "how you make a proton from quarks" in this theory**: place three
component-lumps within roughly two core radii (d ≲ 4 ≈ r_half) and the baryon
self-assembles, radiating the excess as binding violence; beyond that, most
of each quark is lost to dispersal first.

## 3. Flavored baryons (proton vs neutron analog) [solved, stability open]

`flavored_qball.py`: Newton relaxation of the 3-component radial system with
unequal frequencies, Φ_a = f_a(r)e^{iω_a t} (s stays static, ansatz exact).
Continuation from the symmetric ω=1.42 solution converges smoothly to at
least Δω = 0.04 (ω = 1.38/1.42/1.42): localized stationary baryons with
**unequal charge partition** Q_a = (76.7, 103.7, 103.7) — the in-model analog
of uud vs udd flavor content. Energy rises along the branch (E 321 → 415).
Open: dynamical stability of the flavored branch (needs a 3D run with a
flavored seed — generator extension: per-component ω/profile).

## 4. Honest differences from QCD (do not oversell)

1. No mesons of any kind (the product potential binds only triples) — real
   QCD's richest sector is absent.
2. "Confinement" is non-existence (dispersal), not a linear string: pulling a
   quark out costs no growing energy; the quark simply stops being a particle.
3. Color is U(1)³/permutation, not SU(3); nothing rotates components into
   each other dynamically (at η=0), so per-component charge is separately
   conserved [verified: Q_a flat in every run].
4. The gauge charge is the diagonal U(1) — all three quarks carry the SAME
   sign of it; there is no analog of charge −1/3 vs +2/3 within one baryon.

## Files

- `v71/QUARK.md` (this), `v71/analysis/flavored_qball.py`,
  `v71/results/{q1,q4,q5}.png`, `v71/results/q*_qcomp.tsv`,
  `v71/results/flavored_profile_last.txt`
- `sfa/analysis/sfa_qcomp.c`, `sfa/seed/gen_qball_quark.c` (new tools)
- SFAs: `/space/scp/v70/q{1,2,3,4,5}_*.sfa` (per-frame visual record)

---

## UPDATE — Flavored baryon: 3D-stable, flavor conserved [measured, fb1]

The Δω=0.04 BVP solution (ω = 1.38/1.42/1.42, Q_a = 76.7/103.7/103.7) was seeded
in 3D (`gen_qball_flavored`, exact 4-column profile) and evolved at g=0, N=128,
L=20, T=300:

- **Per-component charges machine-flat**: Q_a = 76.60/103.62/103.62 at every
  frame, ratio 0.7393 constant to 4 digits. No equilibration, no shedding.
- **The three internal clocks stay distinct**: volview's inline frequency
  analysis reads W = [1.3799, 1.4199, 1.4199] at t=298 — the BVP frequencies
  to 4 decimals after 300 t.u. of full nonlinear evolution. No synchronization.
- Component geometry static: the low-ω flavor is the COMPACT one
  (rms 3.19 vs 3.37 — μ_tail = √(m²−ω_a²) is larger at lower ω).
- **Conclusion: the theory's baryons form a stable multiplet labeled by
  frequency/charge partition — the proton/neutron-analog distinction is
  dynamically real, not just a stationary point.** (Caveats: T=300, g=0,
  unperturbed; a perturbed/gauged stability test is the follow-up.)

### Viewer support (volview)

volview now renders the substructure directly (and crashed on >16-column files
before this work — fixed):
- **key 0 / -view 5 FLAVOR**: R/G/B = |Φ₀|²/|Φ₁|²/|Φ₂|² under common
  normalization (equal flavors → white; q4 assembly renders as red+green+blue
  quarks merging to a white baryon) + inline per-component clock W_a and
  charge Q_a printout per frame.
- **key C / -view 6 CLOCK**: per-voxel local frequency (R=fast, B=slow/anti-
  charge, brightness=density) — discriminates ± balls and shows weff(r).
- key 8 U(1) gauge (R=|E|, G=|Φ|², B=|A|), key 9 charge (R=+ρ_Q, B=−ρ_Q);
  field/velocity views are phase-invariant (complex moduli) on 24/30-col files.
- Honest note: at Δω=0.04 the flavored ball's color stratification is subtle
  (white core, faint cyan rim — the longer-tailed ω=1.42 pair); the sharp
  discriminators are the printed clocks and Q_a. Deeper detuning would show
  stronger color separation if the branch extends.

Artifacts: `sfa/seed/gen_qball_flavored.c`, fb1 SFA + `v71/results/fb1_qcomp.tsv`,
renders `v71/results/{vv_fb1_grid,vv_flavor_grid}.png`.

---

## UPDATE 2 — Flavored nucleus: flavor survives fusion [measured, fn1]

Deuteron-analog with flavor content: the flavored baryon (ω=1.38/1.42/1.42,
Q_a=76.6/103.6/103.6) fused with a symmetric g=0 ball (ω=1.42, Q=209.6) at
D=8, co-phase, g=0, T=300 (`gen_qball_flavored` extended to multi-ball with
mixed 2-col/4-col profiles):

- **Fusion completes by t≈100** (single object, component centroids coincide).
- **The flavor partition SURVIVES**: Q_0/Q_1 stays in [0.828, 0.849] all run
  (seed 0.832; equilibration would drive it to 1). Fusion radiation removes
  charge nearly flavor-proportionally (−9.0% / −8.7%).
- Final composite: Q_a = (147.5, 178.1, 178.1) — **a flavored nucleus**.
- **Two distinct internal clocks inside one droplet**: volview inline analysis
  reads W = [1.3655, 1.3968, 1.3968] at t=298 — the Δω split persists
  (0.040 → 0.031), both clocks redshifted from the constituents' values as the
  larger composite binds deeper (branch physics: bigger Q → lower ω).
- Visual: a single white ball with a faint cool rim (the comp-0 deficit);
  the quantitative discriminators are the printed W_a and Q_a.
  Render: `v71/results/vv_fn1_grid.png`; data `v71/results/fn1_qcomp.tsv`.

Interpretation: because the three component charges are separately conserved
and fusion radiation is roughly flavor-democratic, **flavor is an
approximately conserved nuclear attribute** in this theory — composites
remember their constituents' flavor content even after total fusion. (Caveat:
g=0, one composition, T=300; the flavor-dependent contact-force channel —
one beating + two locked components during approach — was not separately
resolved at D=8 fusion speed and would need a larger-D run.)

### Addendum — internal flavor-dipole oscillation (user observation, confirmed)

Viewing fn1's final frames in the flavor view (key 0) shows the component cores
remain DISTINCT inside the merged droplet — confirmed quantitatively: the
comp-0 vs comp-1,2 density peaks stay displaced and slosh through each other
(peak offset swinging +0.9 → −0.6 across t=275–298; centroid rms 0.088).
Dominant displacement mode ω ≈ 0.12–0.155 (period ~40–50 t.u.) — a mechanical
flavor-dipole oscillation, the two flavor fluids oscillating in each other's
binding well: the structural analog of the nuclear GIANT DIPOLE RESONANCE
(protons sloshing against neutrons), excited here by the asymmetric fusion.
Low-frequency content at the clock-beat scale Δω=0.031 is present but
resolution-limited (window = one beat period). Permanent amplitude asymmetry:
comp-0 core |Φ|=0.60 vs 0.66. A T≳1000 run would resolve mode vs beat cleanly.

### Addendum 2 — the flavor-dipole mode resolved [fn2, T=1000 restart]

Restarting from fn1's final frame for 1000 t.u. (frequency resolution 0.0063):

- **The internal sloshing is a clean collective resonance at ω ≈ 0.137–0.144**
  (period 44–46 t.u., line width ~0.01), with NO power at the clock-beat
  frequency Δω = 0.031 — a mechanical mode, not beat modulation. The giant-
  dipole-resonance identification stands on its own dynamics.
- **High-Q and damped**: displacement rms 0.088 → 0.046 over 1000 t.u.
  (damping time ~10³ t.u., Q-factor ω·τ ≈ 170) — the nucleus rings.
- **Flavor retention essentially exact long-term**: Q_0/Q_1 = 0.8284 → 0.8273
  over 1000 t.u. (−0.13%); residual drain flavor-proportional (−2.0%/−1.9%).
  Clocks at t=1000: W = [1.3721, 1.3990, 1.3990] — still split (Δω = 0.027).

Data: `v71/results/fn2_qcomp.tsv`; SFA `/space/scp/v70/fn2.sfa` (201 frames).
