# V52 Plan: Quantitative Physics Measurements

Re-measure fundamental physical quantities under the V50/C4 equations
(Cosserat strain α=0.1, curl²-hardening β=0.5). Previous measurements
(V33/V34/V43) used the V44 base equations without these terms.

All runs use pre-converged proton templates (v43/proton_template.sfa),
NOT analytical seeds.

---

## Test 1: Force Power Law — F(D)

**Goal**: Measure the exponent n in F ∝ 1/D^n between two protons.
V33 measured n≈1.8 for braids (mixed gravity+EM). Newton requires n=2.

**Setup**:
- Two UUD protons at separations D = 15, 20, 25, 30, 40, 50
- N=384, L=60 (enough room for D=50 with boundary clearance)
- Absorbing BC, T=50 (measure initial acceleration only)
- Extract centroid acceleration a(D) = d²x/dt² from first 10 time units
- 6 runs, snap_dt=1, diag_dt=0.5

**Analysis**:
- Plot log(a) vs log(D)
- Linear fit gives exponent n
- Compare with n=2 (Newton) and n=1.8 (V33 braids)

**Expected**: C4 hardening may change the effective range by modifying
the shell structure. The proton's isotropic depletion profile should
give a cleaner power law than the anisotropic braid.

**Cost**: 6 × 10 min = 1 hour GPU

---

## Test 2: Depletion Profile — δρ(r)

**Goal**: Measure the radial depletion profile around a single proton.
V34 measured δρ ∝ 1/r^1.2 for a braid (power law, not Yukawa).

**Setup**:
- Single UUD proton, N=384, L=100
- Absorbing BC, T=300 (long equilibration)
- Extract radial profile ρ(r) from the final frame
- Compute δρ(r) = ρ(r) - ρ_bg

**Analysis**:
- Fit δρ(r) to power law 1/r^n for r = 10..50 (avoid core and boundary)
- Fit to Yukawa e^{-mr}/r for comparison
- The depletion sets the gravitational potential shape

**Expected**: The proton (isotropic) should give a cleaner profile than
the z-aligned braid. If n→2, the gravitational force automatically
becomes inverse-square.

**Cost**: 1 × 35 min GPU

---

## Test 3: Gravitational Coupling Constant — C_proton

**Goal**: Extract the coupling constant C in F = -C × ∇ρ for protons.
V33 measured C=186 for braids but that mixed gravity+EM.

**Setup**:
- Use V51 gradient data (already collected)
- Or: repeat with multiple gradient strengths (ΔA = 0.05, 0.10, 0.15)
- Measure drift velocity v_drift for each gradient
- C = m_inert × v_drift / ∇ρ

**Analysis**:
- Linear fit of v_drift vs ∇ρ gives C/m_inert
- If C/m_inert is constant across gradient strengths → F ∝ ∇ρ confirmed
- Combined with Test 1 gives G_eff

**Expected**: The V51 gradient run showed Δx≈6 over T=400 with
ΔA=0.10 (0.15-0.05). Drift rate ~0.015/t.

**Cost**: 3 × 35 min = 1.7 hours GPU (if new runs needed)

---

## Test 4: Charge-Dependent Force Ratio

**Goal**: Measure the EM force between same-chirality vs opposite-chirality
proton pairs. V34 measured +27%/-57% for same/opposite braids.

**Setup**:
- UUD+UUD pair at D=25, v=0 (same chirality, both charge +1)
- UUD+UDD pair at D=25, v=0 (opposite chirality, charge +1/-1)
- N=384, L=50, T=100, absorbing BC
- Measure infall rate (centroid approach velocity) for each

**Analysis**:
- Same chirality: gravity + EM attraction → faster infall
- Opposite chirality: gravity + EM repulsion → slower infall (or repulsion)
- Ratio quantifies EM/gravity relative strength

**Requires**: UDD (neutron) template. May need to generate one by
running gen_phase_confined with UDD chirality and equilibrating T=500.

**Cost**: 2 × 20 min + template generation = 1.5 hours GPU

---

## Test 5: Photon Speed (Already Measured)

**Status**: COMPLETE from V50 EM wave tests.

- Predicted: v = √(1 - η²/m²) = 0.9428c
- Measured: 0.9377c (plane wave, A_bg=0)
- Lean proof: zero sorrys

No additional runs needed.

---

## Test 6: Equivalence Principle — m_grav = m_inert?

**Goal**: Test whether different-mass particles experience the same
gravitational acceleration in the same density gradient.

**Setup**:
- Same gradient as Test 3 (A_high=0.15, A_low=0.05, L=100)
- Three test particles: single braid, UUD proton, deuterium (UUD+UDD)
- Each run independently in the same gradient
- Measure drift acceleration a = d²x/dt²

**Analysis**:
- If a_braid = a_proton = a_deuterium → equivalence principle holds
- If they differ → m_grav ≠ m_inert and the coupling depends on
  particle structure, not just mass

**Note**: The braid will have EM contamination (anisotropic), so the
proton-deuterium comparison is cleaner. Deuterium requires a pre-converged
template (from V42 or newly generated).

**Cost**: 3 × 35 min = 1.7 hours GPU

---

## Test 7: Emergent Gauss' Law — ∮ E·dA = Q_enclosed?

**Goal**: Test whether the θ field (EM sector) satisfies Gauss' law
around a proton. Specifically: does the surface integral of the
outward θ flux over a closed sphere equal a constant (the "charge")
independent of the sphere's radius?

**Setup**:
- Single UUD proton, equilibrated T=300, N=384, L=50
- At the final frame, compute on spherical shells at r = 5, 8, 12, 16, 20, 25, 30:
  - Φ(r) = ∮ (∂θ/∂t) · n̂ dA  (electric flux — E = -∂θ/∂t)
  - Or equivalently: Φ(r) = ∮ |curl(θ)| dA (magnetic flux analog)

**Analysis**:
- If Φ(r) = constant for all r → Gauss' law holds → 1/r² force confirmed
- If Φ(r) ∝ r^α with α≠0 → modified Gauss' law

**Theoretical expectation**: EM_THEORY.md proves ∇·E = 0 exactly (from
∇·(∇×φ) ≡ 0). This means the NET electric monopole charge is zero for
every configuration. But the RADIATION PATTERN of a composite baryon
mimics a monopole (V44 multipole analysis: l=0 carries 54-87% of power).

The test is whether the time-averaged Poynting flux ⟨S⟩ = ⟨E×B⟩ falls
as 1/r² at distances beyond the core. This is equivalent to Gauss' law
for the radiation field.

**Approach**:
1. Load equilibrated proton SFA (last frame from depletion profile run)
2. Compute E = -∂θ/∂t from two consecutive frames (finite difference)
3. Compute B = ∇×θ from spatial derivatives
4. Compute Poynting vector S = E × B at each voxel
5. Integrate |S| · r̂ over spherical shells at various radii
6. Plot Φ_S(r) — should be constant if 1/r²

**Alternative (simpler)**: Just measure θ_rms²(r) on shells. If the
radiation intensity θ_rms² ∝ 1/r² → the Poynting flux is 1/r² →
Gauss' law holds. This avoids computing time derivatives.

**Caveat**: The V34 measurement of θ_rms(r) showed r^(-0.3) decay, but
that was contaminated by periodic BC ghost images. With absorbing BC
and a larger domain (L=50 vs L=25), the measurement should be cleaner.

**Cost**: 1 run (reuse Test 2 data) + analysis tool (~1 hour coding)

---

## Priority Order

1. **Test 2** (depletion profile) — foundational, informs all others
2. **Test 7** (Gauss' law) — uses same data as Test 2, high impact
3. **Test 1** (force power law) — the headline number
4. **Test 3** (coupling constant) — uses existing V51 data
5. **Test 4** (charge ratio) — needs UDD template
6. **Test 6** (equivalence) — needs deuterium template

Tests 2+7 share one run. Tests 1 needs 6 short runs. Total GPU time
for the top 3: ~2 hours.

---

## Shared Resources

- Proton template: `v43/proton_formation/proton_template.sfa` (64³ UUD)
- Neutron template: needs generation (UDD from gen_phase_confined + T=500)
- Collision seed generator: `v51/proton_collision/gen_collision.c`
- Cluster analysis: `bin/cluster_profile`
- All runs: N=384, V50/C4 (α=0.1, β=0.5), f16 precision
