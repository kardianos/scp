# V35 Parameter Space Mapping: Hydrogen-Like Bound States

## Method

Swept 12,000 parameter combinations across:
- ℏ_sim: 30 values, logarithmic from 0.1 to 10,000
- m_eff: 25 values, logarithmic from 0.0001 to 10
- V_depth: 4 values (0.5, 1.0, 2.0, 5.0)
- n_power: 4 values (1.0, 1.1, 1.189, 1.3)

Each point: solve the 1D radial Schrödinger equation, compute bound
states for l=0 and l=1, score against hydrogen observables.

Scoring targets:
- Bohr ratio (<r>₀ / r_braid): 53,000
- Energy ratio (E₁/E₀): 0.25 (hydrogen 1/n² scaling)
- Bound states: 2-10 for l=0
- l=1 has at least 1 bound state

Runtime: 467 seconds (26 evaluations/second).

---

## Best Match to Hydrogen

**Score: 0.0004** (essentially perfect)

| Parameter | Value |
|-----------|-------|
| ℏ_sim | 126.9 |
| m_eff | 0.00287 |
| V_depth | 5.0 |
| n_power | 1.0 |

| Observable | Result | Target |
|-----------|--------|--------|
| Bohr ratio | **52,957** | 53,000 |
| E₁/E₀ | **0.2500** | 0.2500 |
| N_bound (l=0) | 8 | 2-10 |
| l=1 exists | yes | yes |

---

## Key Finding 1: One-Parameter Family of Solutions

The Bohr ratio depends on the RATIO ℏ²/m_eff, not the individual
values. Any (ℏ, m_eff) pair satisfying the constraint gives the correct
orbital size. This is identical to real physics: a₀ = ℏ²/(m_e × e²).

All viable combinations (Bohr ratio within 26,000-106,000):

| ℏ_sim | m_eff | ℏ²/m_eff | Bohr ratio |
|-------|-------|----------|------------|
| 7.9 | 0.0001 | 623,000 | 58,677 |
| 11.7 | 0.00026 | 527,000 | 49,730 |
| 17.4 | 0.00026 | 1,164,000 | 55,008 |
| 38.6 | 0.00026 | 5,698,000 | 53,840 |
| 38.6 | 0.00068 | 2,183,000 | 51,568 |
| 57.4 | 0.00287 | 1,145,000 | 54,105 |
| 85.3 | 0.00287 | 2,534,000 | 59,847 |
| **126.9** | **0.00287** | **5,605,000** | **52,957** |
| 188.7 | 0.0316 | 1,126,000 | 53,218 |
| 280.7 | 0.0316 | 2,493,000 | 58,865 |
| 417.5 | 0.0316 | 5,513,000 | 52,088 |
| 621.0 | 0.348 | 1,108,000 | 52,344 |
| 923.7 | 0.348 | 2,452,000 | 57,899 |
| 1373.8 | 0.348 | 5,422,000 | 51,234 |
| 2043.4 | 3.83 | 1,090,000 | 51,486 |
| 4520.4 | 3.83 | 5,335,000 | 50,393 |

The constraint ℏ²/m_eff clusters around **1.1M and 5.5M** — two
discrete families depending on which V_depth value is used. Within
each family, ℏ and m_eff can trade off freely.

## Key Finding 2: n_power = 1.0 Required for Exact Hydrogen Spectrum

ALL top-20 matches have n_power = 1.000 (pure 1/r Coulomb potential).
Our measured potential has n ≈ 1.189. The difference:

| n_power | Best E₁/E₀ | Hydrogen target |
|---------|-----------|----------------|
| 1.000 | 0.2500 | 0.2500 ✓ |
| 1.100 | 0.2619 | too high |
| 1.189 | 0.2719 | too high |
| 1.300 | 0.2897 | too high |

The 1/r² energy scaling (E_n ∝ 1/n²) is SPECIFIC to the 1/r potential.
Any n > 1 gives a different progression. This means:

**For an exact hydrogen spectrum, the long-range EM potential must be
pure 1/r (Coulomb), not our measured 1/r^1.2.** The 0.2 deviation
in the exponent matters for the spectrum but not for the existence of
bound states.

Possible resolution: the 1/r^1.2 measurement was from the φ-depletion
(gravity). The θ-mediated (EM) force may have a DIFFERENT exponent —
potentially closer to 1/r — since it's carried by a different mechanism
(massless wave exchange vs nonlinear depletion).

## Key Finding 3: Effective Mass is Very Small

The viable m_eff values range from 0.0001 to 3.83, but for physically
reasonable ℏ_sim (< 1000), m_eff ≈ 0.0003 to 0.03.

| m_eff / m_braid | Interpretation |
|-----------------|----------------|
| 0.0001 / 5000 = 2×10⁻⁸ | Extremely light perturbation |
| 0.003 / 5000 = 6×10⁻⁷ | Very light (like real m_e/m_p ≈ 5×10⁻⁴) |
| 0.03 / 5000 = 6×10⁻⁶ | Light |

The physical electron/proton mass ratio is m_e/m_p ≈ 1/1836 ≈ 5×10⁻⁴.
Our m_eff/m_braid ratios are 100-1000× smaller than this. However, the
braid mass in code units is E_braid ≈ 5000, so m_eff ≈ 0.003 gives
m_eff/E_braid ≈ 6×10⁻⁷ — much smaller than the physical ratio.

This suggests the "electron" in this theory is MUCH lighter relative to
the "proton" than in real physics. Or: the effective mass depends on
how the electron perturbation couples to the potential.

## Key Finding 4: 8 Bound States in Ground Configuration

All viable parameter combinations give N_bound = 8 for l=0 with
n_power=1.0. With l=0,1,2,3 angular momentum states, the total
number of orbitals is comparable to real atoms (though the exact
occupancy depends on exclusion principles not yet in the theory).

---

## Physical Interpretation

### What determines ℏ²/m_eff?

This ratio is the fundamental scale that sets the atomic size.
In real physics: ℏ²/(m_e e²) = a₀ = 53 pm.

In our theory, two things determine this:
1. **ℏ_sim**: the field quantum ε × dx (substrate granularity)
2. **m_eff**: the effective mass of a θ-perturbation in the EM potential

Neither is currently derived from the simulation. They are the
REMAINING FREE PARAMETERS of the theory. Everything else (V_depth,
force law, charge dependence) comes from the V34 simulations.

### The spectrum constraint

The energy level spectrum is EXACT hydrogen (E₁/E₀ = 0.25) only for
a pure 1/r potential. Our measured potential is 1/r^1.2. This means:
- The EM potential at long range must be closer to 1/r than the gravity
  potential (which is 1/r^1.2)
- OR: the electron orbits at such large r that the potential has
  effectively become 1/r (the exponent may approach 1.0 at very large r)
- OR: the spectrum deviates slightly from hydrogen (E₁/E₀ ≈ 0.27
  instead of 0.25) — a testable prediction

### Connection to V34 results

| V34 quantity | Used in Phase 2 as |
|-------------|-------------------|
| δρ(r) ∝ 1/r^1.2 | V_grav(r) shape |
| θ force +27% | V_em = 0.27 × V_grav |
| Braid radius r≈5 | r_braid, repulsive core |
| Winding = charge | Potential sign depends on relative W |

---

## Next Steps

1. **Determine the θ-mediated potential exponent**: Is V_em(r) closer
   to 1/r than V_grav(r)? Measure the θ_φ DC component at large r
   in the 6-field simulation.

2. **Derive m_eff from θ dynamics**: The effective mass of a test
   perturbation in the θ field comes from its dispersion relation.
   Compute: ω(k) for small-amplitude θ waves near the braid.
   m_eff = ℏ ω / c² at the bottom of the band.

3. **Field quantization**: Add ε rounding to the simulation. The
   resulting ℏ_sim = ε × dx is then DETERMINED, not free. Combined
   with the derived m_eff, this fixes the Bohr ratio with zero free
   parameters.

4. **θ self-interaction**: If θ forms solitons (θ-braids), their mass
   IS m_eff. The orbital radius is then fully predicted.

## Files

- `param_sweep.py` — parameter space sweep (12,000 evaluations)
- `data/param_sweep.tsv` — full results table
- `data/bound_states.tsv` — Phase 2 eigenvalue results
- `data/V_eff.tsv` — effective potential from V34
