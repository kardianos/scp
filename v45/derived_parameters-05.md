# Derived Parameters for SCP Field Theory — Revision 05

## Synthesis of -01 through -04

This revision brings forward what works from each prior document, corrects
errors of analysis, and proposes concrete next steps grounded in the actual
V45 data.

---

## 1. The Virial Mismatch is Misdiagnosed

### The claim (-03, carried into -04)
"E_grad ≈ 65,000, E_mass ≈ 42,000, E_pot ≈ -3,200 gives a 20:1 imbalance.
The particles are 'exploding solitons.'"

### The correction
These numbers do not appear in the V45 data. The actual V45 D=80 averages
(t=200-400) are:

| Component | Value | Background contribution |
|-----------|-------|------------------------|
| E_phi_kin | 51,839 | ~51,600 (background oscillation) |
| E_mass | 51,684 | ~51,500 (background m²φ²) |
| E_grad | 278 | ~60 (background k²φ²) |
| E_pot | -90 | 0 (P_bg≈0 due to phase offsets) |
| E_theta_kin | 92 | 0 (θ sourced only by braids) |

The dominant energy terms (E_kin ≈ E_mass ≈ 52k) are almost entirely the
**uniform background oscillation**, not the solitons. The background field
φ = A_bg cos(kz+δ) has kinetic energy ½ω²A²_bg per mode integrated over
512³ voxels ≈ 52,000. The soliton perturbation is a DEPLETION on top of
this — it reduces E_mass locally, giving negative "soliton E_mass."

Applying Derrick's virial theorem to the total energy is meaningless when
99% of the energy is background. The virial must be computed for the
perturbation δφ = φ - φ_bg, which requires per-voxel background subtraction
— not available from the global diagnostics.

**The solitons may well be in virial equilibrium.** We cannot determine
this from the current data. The "exploding soliton" diagnosis is
unsupported.

### What we CAN say from the data

The ΔE decomposition (each D relative to D=80 baseline) reveals the
actual interaction physics:

| D | ΔE_kin | ΔE_mass | ΔE_grad | ΔE_pot | ΔE_θkin | ΔE_total |
|---|--------|---------|---------|--------|---------|----------|
| 0 | +531 | +624 | +80 | -146 | +90 | +1260 |
| 5 | +408 | +477 | +70 | -113 | +86 | +1006 |
| 10 | +182 | +198 | +24 | -30 | +48 | +464 |
| 15 | -10 | -15 | +24 | -12 | +49 | +80 |
| 20 | +43 | +54 | +10 | -10 | +7 | +110 |
| 40 | +125 | +122 | +6 | -0 | +6 | +266 |

Key observations:
1. **E_pot goes MORE NEGATIVE at close range** (ΔE_pot = -146 at D=0).
   The binding potential IS attractive — V(P) favors overlap.
2. **But E_kin and E_mass increase MORE** (+531 and +624 at D=0).
   The kinetic/mass energy penalty of distorting the field outweighs
   the potential energy gain. This is classic Derrick's theorem in action.
3. **E_theta_kin increases at close range** (+90 at D=0). The θ coupling
   adds repulsive energy when braids are near each other.

**The repulsion is NOT from phase confinement (P→0).** It's from the
kinetic/mass energy cost of field distortion. This is exactly what
Derrick's theorem predicts: in 3D, gradient energy always wins over
potential energy for static solitons. Binding requires a mechanism that
the current equation doesn't have.

---

## 2. Ideas Carried Forward

### From -01: Reducing free parameters
Agree. But the derivations must be non-circular. The virial approach (-02)
is correct in principle but needs proper background subtraction to apply.

### From -02: η sweep
**CARRY FORWARD — highest priority experiment.** The V45 data shows
ΔE_θkin = +90 at D=0 — the θ coupling contributes repulsive energy.
At lower η, this repulsion decreases. At higher η, the θ radiation between
baryons is stronger, which could create a dynamic attraction through
radiation pressure (the "ponderomotive" mechanism from EM_THEORY.md).

The sweep η = {0.0, 0.1, 0.3, 0.7, 1.0, 2.0} at D=40 would reveal:
- Does the repulsive ΔE_θkin decrease at lower η? (Expected: yes)
- Does an attractive ΔE_total appear at some η? (Test: η=0 removes θ entirely)
- Does higher η create radiation-pressure binding? (Test: η=2.0)

The η=0 run is particularly important — it isolates the pure φ-sector
interaction without θ contamination.

### From -02: Massive θ test
**CARRY FORWARD.** Adding m_θ > 0 converts the 1/r θ radiation into a
Yukawa potential e^{-m_θr}/r. At m_θ ≈ 0.5 (code units), the range would
be ~2 code units — exactly the nuclear scale. This directly tests whether
a massive mediator creates binding without changing any other parameter.

### From -03: Phase incoherence diagnosis
**PARTIALLY CARRY FORWARD.** The concern about gen_deuterium seeds producing
phase-incoherent baryons is valid — we saw fragmented cores in the volview
images. However, the energy decomposition above shows the effect is
systematic across ALL separations (same seeds, same grid), so the
differential ΔE(D) comparison cancels most seed artifacts. The relative
ordering should be correct even if the absolute energies are wrong.

### From -03: Topological quench
**CARRY FORWARD with caveats.** Precipitating baryons from a hot compressed
state is physically motivated and avoids seed artifacts. However:
- V30 tried this (FRW expansion) and failed — no spontaneous braid formation.
- The specific parameter changes proposed (m²=1.0, μ=-250, κ=20) are
  underivatived guesses that would need justification.
- The Coriolis term (Ω×) modifies the equation of motion and requires
  kernel changes plus theoretical justification.

A controlled version: high A_bg=0.4 initial sphere with absorbing BC,
standard parameters, T=2000. See if dense field spontaneously fragments
into braids. No kernel modifications needed.

### From -03/-04: Breathing-mode synchronization (Adler locking)
**HIGH INTEREST — the most physically motivated binding mechanism.**

The idea: two baryons breathing at ~150t period can synchronize through
their mutual θ radiation. The synchronization energy is:

    E_sync ∝ μ × ∫ P₁(x,ψ₁) × P₂(x,ψ₂) dV × cos(ψ₁ - ψ₂)

When in-phase (ψ₁=ψ₂), this is maximally negative → binding.
When anti-phase, maximally positive → repulsion.

**The V45 runs may have missed this because T=500 ≈ 3 breathing periods
is not enough for synchronization to occur.** Adler locking in coupled
oscillators requires many cycles — typically 10-50 cycles for weak
coupling. At 150t per cycle, that's T=1500-7500.

**Test**: Run D=40 to T=5000 (33 breathing cycles). Measure the breathing
phase of each baryon over time. If the phases converge, binding will follow.

### From -04: Adiabatic effective potential
**CARRY FORWARD as analytical framework.** Averaging the fast carrier
oscillation (~4t period) to derive a slow effective potential V_eff(D)
between two baryons is the right theoretical approach. This separates
the "nuclear" timescale from the "particle" timescale. The numerical
version: compute ⟨E_total⟩ over one breathing cycle at each D, after
sufficient equilibration (T>1000).

### From -04: Coarse-graining for effective m_θ
**INTERESTING but speculative.** The idea that nonlinear φ-θ mixing
generates an effective m_θ > 0 at the collective level (even though the
bare m_θ=0) is physically motivated — this is how pions get mass from
chiral symmetry breaking. However, demonstrating this requires a proper
renormalization group calculation, not just simulation.

---

## 3. Ideas Rejected

### From -01: m² from k_bg
Conflates physical mass with box artifact. Rejected in -02.

### From -01: Exact 2π/3 phases
CMA-ES optimization found the empirical values for a reason. No basis for
replacing them with exact thirds without validation.

### From -03: Coriolis stabilization (Ω term)
Modifies the equation of motion without Lagrangian justification. The
neutron instability is real physics (free neutron decays) — "fixing" it
by adding a rotation term is fighting the theory, not deriving from it.

### From -03: Non-local integral potential
Violates relativistic causality. Rejected in -03 itself.

### From -03: Specific parameter jumps (m²=1.0, μ=-250, κ=20)
These are underivatived guesses. m²=1.0 is below the braid stability
threshold (V34: m² ≥ 1.25). μ=-250 is 6× current without justification.
κ=20 shifts the saturation point but the effect hasn't been analyzed.

### From -04: Global ∫P constraint as Lagrange multiplier
Non-local. Same objection as -03 §C.

---

## 4. The Real Question V45 Answers

V45 doesn't just say "no binding." The energy decomposition reveals WHY:

**The potential V(P) IS attractive** (ΔE_pot = -146 at D=0).
**But the field distortion cost (ΔE_kin + ΔE_mass) is 8× larger** (+1155).

This is Derrick's theorem: 3D scalar solitons cannot bind through
potential energy alone because the gradient/kinetic energy penalty always
exceeds the potential gain. In the Skyrme model, this is overcome by the
E₄ (quartic) term which scales differently under spatial rescaling.

In our Cosserat equation, what plays the role of the Skyrme E₄ term?
The **η curl coupling** between φ and θ. The θ-mediated force has
different scaling from V(P) — it's a radiation-mediated interaction that
doesn't require field deformation at the interaction point.

**This suggests**: binding must come from the θ sector, not the φ sector.
The V45 E_θkin data (+90 at D=0) shows θ is currently REPULSIVE at close
range. But this may be because η=0.5 is in the wrong regime — the
radiation pressure may become attractive at different η, or the
breathing synchronization effect may require longer timescales.

---

## 5. Proposed Experiments (v46) — Priority Order

### Priority 1: η=0 baseline at D=40 (1 run, ~2 hr)
Remove θ coupling entirely. If ΔE_total is STILL positive, the φ-sector
alone cannot bind (confirming Derrick). If ΔE_total becomes negative,
the θ coupling is what's preventing binding (surprising but informative).

### Priority 2: η sweep at D=40 (5 runs, ~10 hr)
η = {0.1, 0.3, 1.0, 2.0, 5.0}. Find the regime where ΔE_total(D=40)
relative to D=80 changes sign. This identifies whether binding is
possible at any coupling strength.

### Priority 3: Long-run D=40 at T=5000 (1 run, ~20 hr)
Test breathing-mode synchronization. 33 breathing cycles should be
enough for Adler locking if the coupling is above threshold.

### Priority 4: Massive θ at D=40 (3 runs, ~6 hr)
m_θ = {0.1, 0.5, 1.0} with standard η=0.5. Tests Yukawa binding
directly.

### Priority 5: Virial computation from SFA data
Write an analysis tool that computes the virial for the perturbation
δφ = φ - φ_bg on a per-voxel basis using the SFA snapshots. This
resolves the "exploding soliton" question definitively.

### Priority 6: Controlled quench (1 run, ~10 hr)
High A_bg=0.4 sphere, standard parameters, T=2000, absorbing BC.
See if dense field spontaneously fragments into bound structures.

---

## 6. Summary Table

| Idea | Source | Status | Rationale |
|------|--------|--------|-----------|
| Virial μ/κ constraint | -01/-02 | Needs proper data | Background subtraction required |
| η sweep | -02 | **DO FIRST** | Most direct test of binding mechanism |
| Massive θ | -02 | Do | Yukawa mediator test |
| Phase incoherence | -03 | Valid concern | But ΔE comparison cancels most artifacts |
| Topological quench | -03 | Try (controlled) | V30 failed once, needs better conditions |
| Adler breathing lock | -03/-04 | High interest | Needs T=5000 to test |
| Adiabatic V_eff | -04 | Right framework | Analytical + numerical |
| "Exploding soliton" | -03 | **WRONG** | Virial computed on total, not perturbation |
| Coriolis Ω term | -03 | Rejected | No Lagrangian basis |
| Non-local ∫P | -03/-04 | Rejected | Breaks causality |
| Parameter jumps | -03 | Rejected | Underived, m²=1.0 below stability |
| Det potential | -03 | Interesting future | Requires kernel rewrite |
