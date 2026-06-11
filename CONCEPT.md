# SCP — A Gauged Complex Field Theory with Stable Composite Particles

This is a classical field theory simulated on a 3D grid. In its current form —
six complex scalar fields with a gauged diagonal U(1) symmetry — it supports
**absolutely stable particles** (Q-balls) with internal quark-like
substructure, a measured massless Coulomb force, a short-range phase-coherent
"nuclear" force, composite nuclei with mass defects, a stable flavor multiplet,
and collective nuclear modes. All of these are measured in simulation with
conserved-quantity bookkeeping at machine precision. None of them is
quantitatively mapped to real physics; the parallels are structural.

---

## 1. The Fields and Equations

Six complex scalar fields: three "position" fields Φ_a and three "angle" fields
Θ_a (a = 0, 1, 2), written Φ_a = u_a + i v_a, Θ_a = tu_a + i tv_a (12 real
fields), plus a real abelian gauge field A_i with conjugate electric field E_i
(temporal gauge A₀ = 0). The Lagrangian:

    L =  Σ_a [ ½|D_t Φ_a|² − ½|D_i Φ_a|² − (m²/2)|Φ_a|² ]
       + Σ_a [ ½|D_t Θ_a|² − ½|D_i Θ_a|² − (m_θ²/2)|Θ_a|² ]
       − Vt(s)  +  η Re[ Θ̄ · (D×Φ) ]  +  ½|E|² − ½|∇×A|²

    s     = Π_a |Φ_a|²            (triple product of moduli)
    Vt(s) = (μ/2) s / (1 + κ s)   (saturating attractive potential, μ < 0)
    D_μ   = ∂_μ + i g A_μ         (minimal coupling, charge +1 on Φ and Θ)

Standard parameters: m² = 2.25, μ = −41.345, κ = 50, m_θ = 1.6, g = 0.05.
The curl coupling η is retained in the kernel (η = 0.5 historical default);
current particle experiments run η = 0 with the Θ sector inert.

**Symmetry and charge.** The potential depends only on the moduli, so the exact
continuous symmetry is U(1)³ (one phase per component). Only the **diagonal**
U(1) is gauged. Its Noether charge

    Q = ∫ Σ_a ( u_a v̇_a − v_a u̇_a + tu_a tv̇_a − tv_a tu̇_a ) dV

is conserved exactly (discrete lattice Gauss law flat at ~10⁻¹³ without any
projection step — the lattice Hamiltonian is exactly gauge invariant). The two
**relative** phases remain global symmetries: each component's charge Q_a is
*separately* conserved. This three-fold structure is the theory's color/flavor
sector (§4).

**The real limit.** Setting v = tv = 0 and g = 0 reproduces the original real
6-field Cosserat system term-for-term. Everything the project measured before
2026-06 lives in that limit (§9).

---

## 2. Why the Fields Must Be Complex (null results that define the theory)

The real-field theory has **no stable localized structure**, for two
independent, derived reasons:

1. **Static lumps are Derrick saddles.** E(λ) = aλ + bλ³ under spatial
   rescaling has no minimum; with saturation (b > 0) even the saddle
   disappears — static blobs disperse.
2. **Time-periodic lumps (oscillons) radiate fatally.** Any localized real
   field oscillates at ω_φ ≥ m, and the cubic source P = φ₀φ₁φ₂ then pumps the
   third harmonic 3ω_φ ≥ 3m — far above the mass gap, a *propagating* wave.
   The decay is intrinsic to the product potential: independent of amplitude,
   profile, phase, or feedback control. This is why every oscillon, braid, and
   composite "particle" of the real-field era breathed and died.

Both diseases are cured by one ingredient: an internal phase. A rotating
complex field Φ_a = f(r)e^{iωt} has *time-independent* modulus — the potential
source is static (no harmonic radiation) — and the conserved charge supplies
the pressure term (∝ Q²/λ³) that converts the Derrick saddle into a true
minimum. The stable particle of this theory is necessarily a **Q-ball**: it is
not that complexification was chosen and stability followed, but that
stability *requires* the U(1).

---

## 3. The Particle: the Charged Q-Ball

Ansatz Φ_a = f(r)e^{iωt} (all components equal — but see §4), Θ = 0. The
radial profile solves f″ + (2/r)f′ = (m² − ω²)f + 2Vt′(f⁶)f⁵.

- **Existence window**: ω ∈ (ω_min, m) with ω_min² = m² + (μ/9)(2/κ)^{2/3},
  i.e. ω ∈ (1.3087, 1.5) at standard parameters. Gauging shrinks it: at
  g = 0.05 the window is (1.406, 1.5) and Coulomb self-repulsion caps the
  branch at **Q_max = 921** — no larger static ball exists (the
  end-of-periodic-table analog).
- **Stability**: the branch with dQ/dω < 0 is classically stable
  (Vakhitov–Kolokolov). The static evaporation bound E < mQ does *not* govern
  dynamics (balls with E > mQ are classically immortal). Measured: a bare ball
  retains 99.99985% of its charge over 1000 t.u. under absorbing boundaries.
- **The θ-drain and its closure**: with massless Θ and η ≠ 0, the ball's curl
  sources Θ radiation that carries charge away (drain ∝ η² at weak coupling).
  Setting m_θ = 1.6 > ω closes this channel (~400×, Yukawa-bound dressing);
  the gauge field A then supplies the long-range mediator that the mass
  removed. The standard package is therefore **gauged diagonal U(1) +
  m_θ > ω**.
- **Formation**: a uniform charged condensate fragments into Q-balls below the
  amplitude threshold κA⁶ < 1/2 (A* = (2κ)^{−1/6} = 0.464) — the Affleck–Dine
  mechanism, confirmed in both directions (above-threshold condensates are
  φ-stable; at η ≠ 0 a separate first-order-in-η parametric instability
  fragments even those).

This is the object the real-field theory could never produce: a localized
structure whose **existence density never passes through zero** — measured
flat to a few percent with no spectral line at any clock harmonic — and which
persists indefinitely. The components oscillate; the particle does not blink.

---

## 4. Substructure: Quarks, Color, and Flavor

The triple-product potential gives the ball internal structure with strong
(structural, not quantitative) QCD parallels:

- **Quarks = components.** The binding force on Φ_a carries the factor
  Π_{b≠a}|Φ_b|²: a lump in ONE component feels no potential at all and
  disperses as a free wave — isolated quarks cannot exist as particles. Two
  components together still have s ≡ 0: **this theory has no meson sector**;
  only three-component ("color-complete") objects bind.
- **Fractional charge.** The symmetric ball carries exactly Q/3 in each
  component, and each Q_a is separately conserved. The gauge field couples
  only to the diagonal sum, so the three-ness is invisible to all long-range
  physics — substructure is confined to the relative-phase sector.
- **Baryon assembly.** Three displaced single-component lumps self-assemble
  into a ball when seeded within d ≲ 4 ≈ r_half (charge captured: 72%/46%/31%
  at d = 2/4/6 — a race between infall and free dispersal). Capture is exactly
  democratic across components.
- **Flavor = frequency partition.** Asymmetric stationary baryons
  Φ_a = f_a(r)e^{iω_a t} with unequal ω_a (hence unequal Q_a) exist as a
  smooth branch and are **dynamically stable**: a Δω = 0.04 flavored ball
  holds its charge partition and its three distinct internal clocks to
  4 decimals over 300 t.u. Baryons form a stable multiplet labeled by
  frequency partition — the proton/neutron-analog distinction.

Honest differences from QCD: color is U(1)³ (not SU(3)); "confinement" is
non-existence (dispersal), not a linear string; there are no mesons; all three
quarks carry the same sign of the gauged charge.

---

## 5. Forces

Two coexisting inter-particle forces, measured separately:

1. **Long-range Coulomb (gauge sector).** The ball wears a 1/r² electric halo
   carrying its full charge (Gauss flux readout agrees with the Noether charge
   to 0.1–1%). The mediator is measured massless: halo flux flat over
   r ∈ [6,32] to ±2% after the wrapped-box jellium correction, bounding
   m_A < 10⁻³ (1500× below the matter mass). Force = gQ·E confirmed to ≤6% in
   both signs at boundary-clean separations. Because A is U(1)-neutral, the
   long-range force exists *without* a charge drain.
   **Box systematics**: the absorbing sponge acts like conductor walls;
   pair-force measurements are invalid beyond D ≈ L/2 (image forces reach
   100% and can invert signs). Quantitative force work needs L ≳ 2D.
2. **Short-range clock-interference force (the "strong" analog).** Overlapping
   tails interfere coherently: F ∝ cos(Δφ)·e^{−μ_t D}, μ_t = √(m² − ω²) ≈
   0.56. Same particles, same charges, same distance — the sign of the force
   is set by the relative clock phase alone: Δφ = 0 fuses straight through
   Coulomb repulsion; Δφ = π/2 switches the contact force off (pure Coulomb to
   6%); Δφ = π gives enhanced repulsion. Between particles of different
   frequencies the force **beats at Δω** (measured: approach, then a
   4–10×-Coulomb repulsive surge as the accumulated phase sweeps through
   anti-phase). No phase-locking down to D ≈ 11 (lock range < 0.06 in Δω).

---

## 6. The Nuclear Sector

Same-charge balls have no non-merging bound state, so nuclei are **fused
charge droplets** — the liquid-drop picture is literal:

- **Fusion and mass defect.** Two nucleons (Q = 311 each) fuse into a single
  droplet that lands ON the independently computed branch (E/Q to 0.1%) with a
  measured binding energy of **1.6% of rest mass** radiated away during the
  merger; the product rings in its quadrupole mode while settling. A
  three-light-nucleon "Li" binds at 3.5% per unit charge.
- **The charge cap acts by evaporation.** A super-critical merger product
  (Q up to 1051 > Q_max = 921) does not fission: it remains one droplet and
  radiates charge down toward the cap. A *cold* spherical super-critical ball
  sits on the symmetric saddle indefinitely; binary fission has never been
  observed (an ℓ=2-seeded cold test is open).
- **Flavor is an approximately conserved nuclear attribute.** Fusing a
  flavored with a symmetric baryon yields a composite that retains the flavor
  partition (Q₀/Q₁ drift −0.13% per 1000 t.u.; fusion losses
  flavor-proportional) and carries two distinct internal clocks.
- **Collective modes: the giant-dipole-resonance analog.** Inside a flavored
  nucleus the two flavor fluids slosh against each other — a clean mechanical
  resonance at ω ≈ 0.14 (period ~45 t.u., linewidth ~0.01, no power at the
  clock-beat frequency), high-Q (τ ~ 10³ t.u., ω·τ ≈ 170). Visible only in
  per-component views: the aggregate density is one quiet droplet while the
  interior rings.
- **The neutral sector is mortal.** Charge-neutral composites (±ω breathers)
  show 100%-depth existence oscillation — the density passes through zero
  every π/ω — and decay slowly, consistent with the real-sector theorems.
  Only charged objects are absolutely protected. Opposite-charge pairs
  annihilate slowly (charge segregation −4× over 600 t.u.), not never.

---

## 7. Kinematics

- **Exact relativistic particles**: boosted balls satisfy E² = p² + E₀², carry
  lab clock γω and de Broglie phase tilt k = γωv (the tilt is preserved
  dynamically to 4 decimals; group velocity lags nominal by 1–5% from lattice
  dispersion).
- **The action quantum is the charge**: ℏ_eff ≡ p/k = E₀/ω = Q(1 + ε) with
  ε = (E − ωQ)/ωQ ≈ 0.03 (Pohozaev identity). Measured to 3–5%. This is the
  generic U(1)-soliton identity — stability, mass, and quantization are the
  same conserved quantity — not a derivation of quantum mechanics.
- **Internal clocks are physical**: the local phase rate ω̃(r) = ω + g·a₀(r)
  carries the Coulomb potential (temporal gauge); per-component clocks are
  measurable from any single frame and label flavor.

---

## 8. What Has NOT Been Established

1. **Quantitative contact with real physics.** No dimensionless ratio from
   this theory matches a measured constant. "Proton", "quark", "He", "GDR"
   are structural analogies. g and m_θ (and historically η) are free inputs.
2. **Gravity.** Nothing in the current framework produces a long-range
   attractive tensor force. (The v33–v51 gradient-drift results belong to the
   real-field era and its dense backgrounds; a clock-gradient mechanism
   remains an open conjecture.)
3. **SU(3) / mesons.** The color structure is abelian; there is no meson
   sector and no gluon analog. The η-curl coupling treats the component index
   as a spatial index — its role at η ≠ 0 in the complex theory is largely
   unexplored.
4. **Fission.** The Q_max cap is established on the static branch and via
   evaporation dynamics; binary fission has never been observed.
5. **Long-time limits.** "Stable" means machine-flat over ~10³ t.u. and
   theorem-backed (VK + closed radiation channels), not proven for arbitrary
   times.

---

## 9. Relation to Earlier Formulations (brief)

The project's first arc used the **real** 6-field theory (v28–v53) and
produced a rich phenomenology — oscillons, braids, UUD/UDD composites, a
deuterium analog, winding-dependent forces, gradient drift — all built on
localized structures that were, without exception, **mortal**: breathing
objects whose density passed through near-zero each cycle and which radiated
away through the 3ω product channel, with lifetimes from a few to ~10³ t.u.
The v53 stability audit and the v64–v65 crisis made this precise (Derrick
saddle + above-gap harmonic radiation, §2), closing the era: *that theory has
no particles*.

The U(1) reformulation (v66) replaced "long-lived" with **absolutely stable**:
the charged Q-ball is the first object in this project whose existence density
never blinks and whose charge cannot radiate away. Gauging the U(1) (v68–v69)
restored the long-range force that stabilization had cost, and the
substructure/nuclear program (v70–v71) was built on top of it. Algebraic
side-tracks (v54–v63: octonionic mass relations, number-type maps,
induced-metric gravity) are documented in their version directories and
DISCOVERIES.md; their surviving constraints (e.g. genuine inputs = a_lepton
+ α from the v59 audit) are independent of the simulation framework.

---

## 10. Technical Infrastructure

- **Kernel (kernel-v3)**: `sfa/sim/scp_sim.c` (CPU) and `.cu` (GPU) — one
  config-driven binary covering real (12-array), complex (`complex_phi=1`,
  24-array), and gauged (`complex_gauge=1`, compact Kogut–Susskind links,
  30-column output) modes. g=0 dispatches byte-identically to the ungauged
  path; the discrete Gauss law is conserved by symplectic structure.
- **Format**: SFA (`sfa/format/sfa.h`), zstd-compressed columns; 12/24/30-
  column conventions.
- **Viewer**: `sfa/volview` (OpenGL ray marcher) — field/velocity/accel views
  are complex-aware (phase-invariant moduli); U(1)-gauge view (|E|, |A|),
  charge view (±ρ_Q), flavor view (per-component RGB + inline per-component
  clock/charge analysis), local-clock view. Headless export via `-snapshot`.
- **Seed generators** (`sfa/seed/`): `gen_qball_pair/boost/bath/multi` (whole
  balls), `gen_qball_quark` (single components), `gen_qball_flavored`
  (per-component profiles and frequencies); radial profiles from
  `radial_qball` and the gauged shooter (`v69/theory/`).
- **Analysis** (`sfa/analysis/`): `sfa_qball_track` (cluster tracking),
  `sfa_qcomp` (per-component charge/clock bookkeeping), `sfa_slice` +
  `render_slices.py` (slice/lineout visual pipeline), plus legacy real-field
  tools.
- **Execution**: the `scp-runner` MCP server (local CPU / Vast.ai GPU);
  archival to B2 via rclone (`scpsfa:scpsfa/v70/`, `/v71/`).
