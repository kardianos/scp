# Skeptic Agent Report: Compact Soliton Search Feasibility

## Executive Summary

The approach of finding self-sustaining compact field configurations via random initialization + time evolution + evolutionary selection is confronting **four fundamental obstacles**, each of which may be a showstopper. None alone disprove the concept, but together they suggest the task may be inherently mismatched to the tools being applied.

---

## CONCERN 1: No Compact Soliton Exists in This System

### Problem
The model uses:
- Wave equation: `∂²φ/∂t² = ∇²φ - m²φ - V'(P) + coupling`
- Binding potential: `V(P) = (μ/2) P²/(1+κP²)` with μ=-41.345, κ=50
- Topology: 3 coupled scalars with product interaction P = φ₀φ₁φ₂

**Historical evidence suggests this system lacks *any* stable soliton, not just ones we haven't found yet:**

1. **Braid3 (infinite helical tube)** works ONLY under periodic BC (wrapping). The structure requires spatial periodicity to exist at all. Attempt to truncate it (Gaussian z-envelope): still depends on phase coherence across the z-range; fails at boundaries.

2. **Borromean rings, tori, loops**: All decomposed on timescale ~15-30 steps (RESULTS.md, RESULTS_crossed_v2.md). Not due to search inefficiency or initial condition tuning — they fragmented structurally.

3. **Previous v26 work (cited in MEMORY.md)**: Found that geometric solitons in 3D Skyrmion theory require EITHER:
   - Topological charge (quantized monopole/instanton)
   - Wrapping/periodic structure (infinite tube)

   The 3-field coupling here has NO topological charge. The system isn't topologically protected.

4. **Absorbing BC changes the game**: With periodic BC, radiation circles back and interferes constructively/destructively with the core — artifacts can look like stability. With absorbing BC (absorbs all outgoing radiation), there's nowhere for energy to hide. Nothing survives.

### Is This a Showstopper?
**Possibly yes.** If the system fundamentally cannot support a localized, stable configuration in 3D (only in 1D periodic tubes), then the search is chasing a phantom.

### Test It
**Experiment A1**: Take the "best" configuration found (truncated helix at t=359, RESULTS.md):
- Run the EXACT SAME initial condition with **N=256, L=30** (much finer, larger box)
- Hypothesis 1 (artifact): configuration dissolves faster at higher resolution
- Hypothesis 2 (real): it survives longer or more stably (genuine soliton)
- Hypothesis 3 (neither): it behaves identically (neutral)

Expected result: **Hypothesis 1** (artifact). The "oscillating breathing" at N=128 is likely numerical instability near the Nyquist limit (tube radius ~2 cells).

**Experiment A2**: Compare periodic vs. absorbing BC on the truncated helix:
- Run truncated helix to T=500 with **periodic BC** (current)
- Run to T=500 with **absorbing BC** (5% damping zone of width 3 cells)
- Plot E_pot(t) on same graph

Expected result: Periodic BC shows the oscillating survival we see. Absorbing BC shows monotonic decay to zero by t~150.

**Why this matters:** If the soliton is genuine, absorbing BC shouldn't matter (waves radiate but configuration stays bound). If it's an artifact of periodic reflections, absorbing BC will kill it.

---

## CONCERN 2: Absorbing BC Removes Energy Indistinguishably from Radiation

### Problem
The absorbing boundary condition applies **quadratic damping** in the outer 3-cell layer:
```c
/* Pseudo-code from v37_crossed_v2.c */
damping_rate = 0.02 / step
```

This removes energy at a rate that is **hard-wired and unconditional**. Once energy crosses the absorbing layer, it's gone.

**The issue:** We cannot distinguish between:
- **(A)** A configuration that radiates energy (inherently unstable)
- **(B)** A genuine soliton + a perfectly calibrated absorber

Consider the truncated helix (RESULTS.md):
- t=0: E_pot = -17.1
- t=359: E_pot = -0.12 (appears dead)
- But E_pot was oscillating: max +40, min -0.1, average ~-15

Did it die, or did we measure it at a breathing node during a cycle?

**Mechanism:** The damping zone absorbs ALL waves—breathing modes, radiation, everything. If the breathing mode tries to escape radially, the absorber kills it before it can re-enter. The effective "Q-factor" of the system is now dominated by the absorber, not the physics.

### Is This a Showstopper?
**Effectively yes, for validation.** We cannot use energy measurements to distinguish genuine solitons from dissipation artifacts when the dissipation is external and unconditional.

### Test It
**Experiment B1**: Measure radiation flux separately:
- Interpolate field values at r = L - 5 (just inside absorbing zone)
- Compute radial energy flux: `F_r = (∂φ/∂t)(∂φ/∂r)` summed over shell
- Plot F_r(t) on graph

For a genuine soliton: F_r should be **zero** (no net outward flux, or oscillating around zero).
For a configuration bleeding energy: F_r should be **monotonically positive** (energy leaving).

Expected result: The truncated helix will show **positive F_r** throughout, confirming it's radiating.

**Experiment B2**: Reduce absorbing strength and measure survival time:
- Run same truncated helix with damping_rate = 0.01/step (half)
- Run again with damping_rate = 0.005/step (quarter)
- Plot E_pot(t) and survival time vs. damping rate

Genuine soliton: E_pot envelope should stabilize as damping→0.
Artifact: E_pot will still decay monotonically, just slower.

Expected result: Decay time scales linearly with 1/damping_rate, confirming the dissipation is NOT the physics.

---

## CONCERN 3: Energy Balance: Binding ≠ Gradient Pressure

### Problem
On the grid (N=32–128), three energy terms are comparable in magnitude:

| Term | Per-cell scale | At φ~0.5 |
|------|---|---|
| Kinetic | ~v² | ~0.25 (v~0.5) |
| Mass | m²φ² | m²=2.25 → 0.56 |
| Gradient | φ²/dx² | At N=32, dx~1: ~0.25 |
| Binding potential | V(P)_max | |mu|/(2κ) ≈ 0.41 per cell |

**The binding potential is the *weakest* term.** V(P) provides at most ~0.41 units of holding power per cell.

For a compact object of radius R~3 cells (volume ~100 cells), binding energy is at most:
- E_bind ~ 0.41 × 100 = 41 units

But gradient energy from maintaining shape grows as:
- E_grad ~ (gradient radius)² × (N cells touched) ~ (0.5/dx)² × 100 ~ 100+ units at N=32

**Consequence:** Gradient pressure (field spreading) is ~2-3× *stronger* than the binding force. The system naturally wants to diffuse outward, not stay compact.

### Is This a Showstopper?
**Not inherently, but it explains systematic failure:** You cannot stabilize a configuration against gradient pressure using a binding force that is weaker. You'd need:
1. Much stronger V(P) (larger μ), OR
2. Lower gradient energy (softer fields, smaller |∇φ|), OR
3. Topological protection (prevents spreading)

None of these is available.

### Test It
**Experiment C1**: Energy partition analysis:
- For the truncated helix at T=0, compute:
  - E_grad / E_total
  - E_mass / E_total
  - E_pot / E_total
  - Compare: which dominates?

Expected result:
```
E_grad:  ~0.30 (dominates)
E_mass:  ~0.25
E_pot:   ~0.01 (weak)
Kinetic: ~0.40 (initial wave energy)
```

Verify from RESULTS.md, t=0 snapshot in energy column.

**Experiment C2**: Synthetic test with stronger binding:
- Create a braid3 variant with μ=-200 (instead of -41.345) — 5× stronger
- All other parameters same, same geometry, same T
- Hypothesis: stronger binding should extend survival time

Expected result: If binding is the limiting factor, survival time increases ~5× (t~1800). If gradient pressure is dominant, no improvement (still dies ~t~300).

If result is **no improvement**, then binding strength doesn't matter — the system is gradient-limited, and no localized soliton is possible.

---

## CONCERN 4: Nyquist Discretization — Features at Grid Limit

### Problem
The truncated helix has:
- Tube radius: R_tube = 3.0 (code units)
- Grid spacing at N=32: dx = 2L/(N-1) ≈ 2
- Tube radius in cells: **1.5 cells**

The **Nyquist limit** for a finite-difference grid is ~2 cells per wavelength. At 1.5 cells, we're *at the limit*. Features smaller than ~2 cells cannot be resolved.

**What happens at the Nyquist limit?**
1. Derivatives are increasingly inaccurate (stencil error grows ~Δx²)
2. Oscillations in the field can produce spurious energy
3. Numerical dispersion (different frequencies propagate at different speeds) dominates
4. The effective "diffusivity" from the discretization can be larger than the physics

**Consequence:** A configuration that "survives" at N=32 might be a **numerical artifact** that vanishes at higher resolution (N=64, N=128) where the tube is resolved by 3-6 cells.

### Is This a Showstopper?
**Yes, for reliability.** Any configuration found at N<64 is potentially fake. It may not survive convergence testing.

### Test It
**Experiment D1**: Convergence study on truncated helix:
- Run truncated helix at **N=32, L=15** (dt scales with dx)
- Run at **N=64, L=15**
- Run at **N=128, L=15**
- Plot E_pot(t), aspect ratio(t), P_int(t) on same graph

Expected result (genuine soliton): curves converge (smaller grid spacing gives same physics).
Expected result (artifact): curves diverge at N>64 (structure collapses as resolution improves).

**Likelihood:** The latter. The tube width of 2 cells is too tight.

**Experiment D2**: Smoother initialization (lower frequency content):
- Increase tube radius from 3.0 → 4.0 → 5.0
- Keep all other parameters the same
- Plot survival time vs. tube radius

Expected result: If Nyquist matters, larger tubes (fewer cell-to-cell gradients) survive longer.

---

## CONCERN 5: Breathing Mode Lifetime — T=10 Steps, Period~2, Too Short

### Problem
From RESULTS.md (truncated helix):
> E_pot oscillates between -40 and -0.1 throughout (strong breathing mode, **period ~20-25**)

At dt ≈ 0.024 (code units), one period is ~20/0.024 ≈ 830 timesteps, or ~20 time units.

The evolutionary/gradient-descent approach evaluates fitness after T=10 time units:
- This is **0.5 breathing periods**
- Less than one full cycle

**Problem:** A configuration could be oscillating healthily but measured at a breathing node (E_pot~0), appearing "dead." Conversely, a doomed configuration could be measured at a breathing crest (E_pot large), appearing "alive."

Worse: Early selection (first 10 time steps) optimizes for **ephemeral properties** (breathing amplitude), not stability. A good breather ≠ a stable soliton.

### Is This a Showstopper?
**Significant, but not absolute.** The optimizer will find locally-maximal breathes, not stable solitons. But it could stumble upon a stable configuration by chance.

### Test It
**Experiment E1**: Extend evaluation time in optimizer:
- Run field_search.py with curriculum: T_opt = [5, 20, 50, **100**] (instead of [5, 20, 50])
- Compare top candidates at T=100 with previous top candidates at T=50
- Do the same configurations remain top-ranked?

Expected result: Different configurations dominate at different timescales (temporal overfitting). This confirms T is too short.

**Experiment E2**: Multi-timescale loss:
- Modify loss function to measure stability over multiple windows:
  ```
  loss = -|E_pot(T=100)| - 0.5*|E_pot(T=50)| - 0.1*|E_pot(T=10)|
  ```
- This biases toward configurations that stay stable across all timescales

Expected result: Different candidates emerge (not just good breathers). Validation step will tell if they're better.

---

## INTEGRATION: Why All Four Concerns Matter Together

Even if one concern is false, the others compound:

1. **No compact soliton exists** (Concern 1) → search is futile
2. **Even if one exists, absorbing BC hides it** (Concern 2) → we can't verify it
3. **Even if we verify it, binding is too weak** (Concern 3) → it should decay anyway
4. **Even if it decays slowly, discretization is cheating** (Concern 4) → the config disappears at higher resolution
5. **Even if it's real at higher resolution, the optimizer measures the wrong timescale** (Concern 5) → selection finds breathers, not solitons

**The scenario most likely from available data:**
- Truncated helix at N=128 is a **breathing mode artifact** of the discrete grid
- It oscillates healthily at the resolution where tube radius ≈ 2 cells
- At N=256 or N=512, the same initial condition disperses (Concern 4)
- At absorbing BC (not periodic), it decays in ~100 steps (Concern 2)
- Binding energy is too weak to actually stabilize it against gradient pressure (Concern 3)
- The "survival" at T=359 is a measurement at a breathing node, not real stability (Concern 5)

---

## Recommendations for Next Steps

### If Skepticism is Unfounded (optimistic path)
1. **Run Experiment A1 & A2** (convergence + absorbing BC)
   - Cost: 2 × 30 min = 1 hour
   - Decision threshold: If truncated helix survives N=256 + absorbing BC, continue
   - If not, move to Recommendations 2–3 below

### If Skepticism is Correct (pessimistic path)
2. **Accept that no compact soliton exists without topological protection**
   - Shift to: lattices, periodic arrays, vortex rings (topologically nontrivial)
   - OR: add a new topological term to the Lagrangian (e.g., Hopf charge, WZW coupling)

3. **If insisting on compact configurations:**
   - Greatly strengthen binding: μ = -200 (test in Experiment C2)
   - Use lower mass: m² = 0.5 (less gradient pressure)
   - Add a confinement potential: V_conf = λ(r/R)⁴ outside radius R
   - Any of these biases the system toward compactness; test which works

4. **Re-evaluate the loss function:**
   - Current: `-E_pot + λ_compact * R_rms²`
   - Better: `-E_pot(T=100) + λ_compact * max(R_rms, 0.1*L) + λ_smooth * grad_energy`
   - Penalize spreading at multiple timescales, not just end state

### Interim Validation
Before running large-scale field_search.py (Experiment Phase 2):
1. **Run Experiment A1 + A2** to rule out numerical artifacts
2. **Run Experiment B2** to rule out absorbing BC camouflage
3. **Run Experiment C1** to quantify energy balance
4. **Only if all pass:** proceed to grid search

---

## Bottom Line

**The system's physics and numerics are conspiring against soliton existence:**
- Physics says binding is weak (Concern 3)
- Topology says no topological protection (Concern 1)
- Numerics says configurations are under-resolved (Concern 4)
- Dissipation says we can't measure stability reliably (Concern 2)
- Timescales say we're optimizing the wrong quantity (Concern 5)

The truncated helix's "oscillating survival" is likely a **breathing mode in the Nyquist-limited regime**, not a genuine soliton.

**Test this now (Experiment A1 + A2, ~1-2 hours) before committing 10+ hours to field_search.py.**

---

## References

- `/home/d/code/scp/v37/RESULTS.md` — Geometry search results
- `/home/d/code/scp/v37/RESULTS_crossed_v2.md` — Fragmentation analysis + absorbing BC effects
- `/home/d/code/scp/v37/src/v37_compact.c` — Compact braid code
- `/home/d/code/scp/v37/backprop/PLAN.md` — Field search plan (Phase 1-2)
- `/home/d/code/scp/MEMORY.md` — Prior work on topological solitons (v26, v28)
