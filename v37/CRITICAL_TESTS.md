# Critical Tests: Validate Soliton Feasibility Before Scaling Search

## Overview
Five experiments to run sequentially. Each is ~30-60 min. **Total: 3-4 hours.** Stop if any fails.

---

## Test A1: Convergence to Higher Resolution

**Hypothesis:** Truncated helix at N=128 is real soliton OR numerical artifact?

**Setup:**
```bash
# Run same truncated helix at 3 resolutions
./v37_compact -geom truncated -N 32 -L 15 -T 500 -o data/conv_N32
./v37_compact -geom truncated -N 64 -L 15 -T 500 -o data/conv_N64
./v37_compact -geom truncated -N 128 -L 15 -T 500 -o data/conv_N128
```

**Analysis:**
- Extract `E_pot(t)` from each `timeseries.tsv`
- Plot on same graph with labels N=32, 64, 128
- Compute: tube width in cells = R_tube / dx
  - N=32: dx=1.0, width=3.0 cells
  - N=64: dx=0.5, width=6.0 cells
  - N=128: dx=0.235, width=12.8 cells

**Pass Criteria:**
- All three curves superimpose (converged solution)
- Breathing amplitude stable or decreasing with resolution
- Final time t=500 has comparable E_pot across all N

**Fail Criteria:**
- N=32 shows strong oscillations; N=64 shows weaker oscillations; N=128 shows decay
- OR: Energy drift increases with resolution
- This means Nyquist artifact (Concern 4 confirmed)

**Expected result:** FAIL. N=32 likely shows false breathing due to 1.5-cell tube radius.

---

## Test A2: Periodic BC vs. Absorbing BC

**Hypothesis:** Is the oscillating survival real physics or periodic reflection?

**Setup:**
```bash
# Modify v37_compact.c:
# Version 1: Keep periodic BC (lines 113-115)
# Version 2: Add absorbing BC (quadratic damping, width=3 cells, rate=0.02/step)

# Run both:
./v37_compact_periodic -geom truncated -N 128 -L 15 -T 500 -o data/periodic
./v37_compact_absorbing -geom truncated -N 128 -L 15 -T 500 -o data/absorbing
```

**Analysis:**
- Plot E_pot(t) from both runs on same graph
- Compute: time when E_pot drops below 1% of initial

**Pass Criteria:**
- Absorbing BC shows similar oscillation pattern (soliton is real)
- Survival time changes by <20% (physics unchanged by boundary)

**Fail Criteria:**
- Absorbing BC shows monotonic decay to zero (periodic reflections were stabilizing)
- Survival time << 500 with absorbing BC
- This means Concern 2 confirmed

**Expected result:** FAIL. Absorbing BC will kill the "soliton" by t~150.

---

## Test B1: Radial Energy Flux

**Hypothesis:** Is the configuration actually radiating, or oscillating in place?

**Code addition to v37_compact.c (in diagnostics section):**
```c
/* Radial flux at r = L - 5 cells */
double radial_flux = 0;
for (long idx = 0; idx < N3; idx++) {
    int i = (idx / NN), j = ((idx/N) % N), k = (idx % N);
    double x = -L + i*dx, y = -L + j*dx, z = -L + k*dx;
    double r = sqrt(x*x + y*y + z*z);

    if (fabs(r - (L-5*dx)) > 2*dx) continue;  /* only r ≈ L-5 */

    for (int a = 0; a < NFIELDS; a++) {
        double dpdt = phi_vel[a][idx];
        double drx = (phi[a][...ip...] - phi[a][...im...]) / (2*dx);
        /* add x-component; similarly y, z */
        radial_flux += dpdt * drx * dV;
    }
}
fprintf(fp_flux, "%.2f\t%.6e\n", t, radial_flux);
```

**Setup:**
```bash
# Run truncated helix with flux diagnostics
./v37_compact -geom truncated -N 128 -L 15 -T 500 -o data/flux
```

**Analysis:**
- Extract `flux(t)` (radial energy flux)
- Plot flux(t) for entire run
- Integrate: total flux = ∫ flux(t) dt

**Pass Criteria:**
- Flux oscillates around zero
- Integral ≈ 0 (no net outward energy)
- Small compared to kinetic energy scale

**Fail Criteria:**
- Flux is consistently positive (energy leaving)
- Integral >> 0 (net outward flow)
- This means configuration is radiating (Concern 2 confirmed)

**Expected result:** FAIL. Flux will be positive, confirming radiation.

---

## Test C1: Energy Partition at t=0

**Hypothesis:** Is binding potential strong enough to overcome gradient pressure?

**Analysis (no new code, use existing diagnostics):**
- From RESULTS.md, truncated helix at t=0, extract:
  - E_grad (gradient energy)
  - E_mass (mass energy)
  - E_pot (binding potential)
  - E_kinetic (initial kinetic from wave initialization)

**Pass Criteria:**
- E_pot / E_grad > 0.3 (binding is comparable to gradient)
- E_pot + E_kinetic > E_grad (binding + kinetic overcome gradient)

**Fail Criteria:**
- E_pot / E_grad < 0.1 (binding is negligible)
- E_pot / E_total < 0.05 (binding is <5% of total)
- This means gradient dominates (Concern 3 confirmed)

**Expected result:** FAIL. E_pot will be ~1-2% of total; gradient will dominate.

---

## Test C2: Stronger Binding Test

**Hypothesis:** Does increasing binding strength extend survival?

**Setup:**
```bash
# Modify KAPPA and/or MU in v37_compact.c
# Current: MU = -41.345, KAPPA = 50.0
# Test 1: MU = -200.0 (4.8× stronger), KAPPA = 50.0
# Test 2: MU = -41.345, KAPPA = 10.0 (weaker feedback, V(P)_max higher)

./v37_compact -geom truncated -N 128 -L 15 -T 500 -o data/mu_strong
# (recompile with MU=-200)
```

**Analysis:**
- Plot E_pot(t) for original (MU=-41.345) and MU=-200 on same graph
- Measure: time when E_pot drops to <1% of initial (in each case)

**Pass Criteria:**
- Survival time increases significantly (e.g., t~1000 for MU=-200)
- Binding is the limiting factor (Concern 3 refuted)

**Fail Criteria:**
- Survival time unchanged or slightly shorter
- Stronger binding doesn't help
- This means gradient pressure is the limit, not binding (Concern 3 confirmed)

**Expected result:** FAIL. Stronger binding will show modest improvement (if any), suggesting gradient dominates.

---

## Test D1: Convergence of Energy with Resolution

**Hypothesis:** Do energy values converge, or do they diverge?

**Analysis (using data from Test A1):**
- From each run (N=32, 64, 128), compute final energy `E_final` at t=500
- Compute: `ΔE_percent = 100 * (E_N128 - E_N) / E_N128`

**Pass Criteria:**
- ΔE_percent < 5% across all resolutions (converged)
- Energy drift% is stable across N (numerical accuracy consistent)

**Fail Criteria:**
- ΔE_percent > 15% (non-converged)
- Energy drift increases with N (discretization error grows)

**Expected result:** MARGINAL. Energy should be nearly converged, but *structure* will not be.

---

## Test E1: Breathing Period Measurement

**Hypothesis:** Is the breathing mode genuine, and what is its period?

**Analysis (using data from Test A1, N=128):**
- Smooth E_pot(t) with moving average (window=10 steps)
- Find local maxima and minima
- Measure time between consecutive maxima: ΔT
- Breathing period = 2 × ΔT (cycle is max → min → max)

**Pass Criteria:**
- Breathing period ≈ 20-25 time units (matches literature, physical)
- Period is stable throughout evolution
- Aspect ratio oscillates in phase with E_pot (geometric breathing)

**Fail Criteria:**
- No clear periodicity (frequency content is noise)
- Period decreases over time (damping without stabilization)

**Expected result:** PASS. Breathing is real, but not proof of stability. Just oscillation.

---

## Decision Tree

```
START
  ↓
Run A1 (convergence) + A2 (periodic vs absorbing)
  ↓
  ├─→ [Both PASS]  → Proceed to C1 (energy partition)
  │
  └─→ [A1 or A2 FAIL] → STOP
        Configuration is Nyquist artifact or absorber-dependent
        → Move to Section "If Skepticism is Correct" in SKEPTICISM_REPORT.md
        → Restart with stronger binding (point 3) or topological terms (point 2)

If proceeded to C1:
  ↓
Run C1 (energy partition) + C2 (stronger binding)
  ↓
  ├─→ [C1 shows E_pot > 5% of total AND C2 shows improvement]
  │   → Binding is viable, continue field_search.py
  │
  └─→ [C1 shows E_pot < 5% AND C2 shows no improvement]
      → STOP. Binding is too weak.
      → Increase MU to -200 (or more) and restart

If all pass:
  ↓
Run field_search.py (backprop/PLAN.md Phase 1-2)
  ├─ Increased confidence in candidate validity
  ├─ Still validate candidates at higher resolution before scaling
  └─ Use Test B1 (flux) on any survivors as extra check
```

---

## Timeline
- **A1 + A2**: 1 hour
- **B1 (if needed)**: 30 min
- **C1**: 15 min (analysis only)
- **C2**: 1 hour
- **E1**: 15 min (analysis only)

**Total: 3-4 hours max.** Stop if any test fails.

---

## Success Metric for All Tests

If truncated helix passes **A1 + A2 + C1 + C2 + E1**, then:
- It's likely a genuine (or near-genuine) soliton
- Proceed to field_search.py with confidence
- Use B1 (flux) as validation on candidates

If truncated helix fails **any of A1, A2, or C2**, then:
- Current approach is fundamentally flawed
- Modify physics (stronger binding, topological terms) before scaling search
- Or accept that compact solitons don't exist in this system
