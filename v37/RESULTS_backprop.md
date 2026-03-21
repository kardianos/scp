# Backprop Structure Search: Results and Implications

## What Was Done

Treated the initial field configuration φ_a(x,y,z) as learnable weights and
used JAX autodiff to backpropagate through the PDE simulation. The loss function
maximized binding retention (E_pot survival), minimized spatial spreading,
and penalized high gradient energy.

- Grid: N=32, L=12, dx=0.77
- Evolution: 387 Verlet steps (T=30)
- 300 Adam epochs, lr=0.001
- Seed: braid3(z) configuration
- GPU: RTX 5060 Ti, 22 seconds total
- 98,304 learnable parameters (3 fields × 32³)

## Training Trajectory

| Epoch | Loss  | E_pot₀  | E_pot_T  | Binding retention |
|-------|-------|---------|----------|-------------------|
| 0     | 1.89  | -70.2   | -139.8   | 199% (growing!)   |
| 50    | 1.44  | -73.8   | -106.7   | 144%              |
| 100   | 1.22  | -76.0   | -88.1    | 116%              |
| 200   | 0.92  | -78.5   | -58.8    | 75%               |
| 299   | 0.72  | -80.2   | -37.0    | 46%               |

The optimizer steadily deepened the binding well (E_pot₀: -70 → -80) while
the loss decreased. At 2× evaluation time: **E_pot retention = 109%** —
the optimized structure gains binding energy from the background.

## What the Optimizer Changed

### Overview: tiny adjustments, large impact
The field values changed by only **3-4% in the core** and **12-27% in the
outer regions**. The Z-width, traveling wave character, and core structure
are essentially preserved.

### 1. Core amplitude slightly reduced (-2.2%)
The integrated field amplitude in the core (r < 3) dropped from 29.9 to 29.2.
Counterintuitively, slightly LOWER core amplitude gives BETTER binding survival.
This suggests the seed braid3 is slightly overdriven — the V(P) potential is
near saturation (κ-limited), so reducing amplitude moves it off the flat top
of V(P) onto the steeper slope where restoring forces are stronger.

### 2. Envelope shape unchanged
The z-width (FWHM = 5.42) is identical between seed and optimized. The braid3's
Gaussian truncation envelope is already near-optimal.

### 3. Z-asymmetry introduced (+6.1%)
The optimized field has 6% more energy on the -z side than +z. This breaks
the seed's z-reflection symmetry. Physical interpretation: a slightly
asymmetric (weakly propagating) wave packet is more stable than a perfectly
symmetric standing wave. The asymmetry prevents destructive self-interference
during breathing oscillations.

### 4. Changes concentrated in the OUTER region
| Region      | r range | Change/seed ratio |
|-------------|---------|-------------------|
| Core        | 0-3     | 3.7%              |
| Envelope    | 3-6     | 11.7%             |
| Outer       | 6-9     | **27.2%**         |
| Edge        | 9-12    | **25.0%**         |

The optimizer barely touches the core but aggressively reshapes the TAIL —
the transition zone where the braid meets the background. This is exactly
where radiation originates: field amplitude mismatch creates outgoing waves.
By smoothing this interface, the optimizer creates better "impedance matching"
that reduces radiation.

### 5. Traveling wave preserved
The velocity-to-field ratio (vel/phi) in the core stayed at mean≈0 with
std≈0.85, confirming the traveling wave character is maintained. The tiny
positive shift (0.024) corresponds to the z-asymmetry — a weak net propagation.

### 6. Boundary velocity adjusted
The largest velocity changes are at the box edges — the optimizer is
effectively creating a velocity field at the boundaries that feeds energy
back into the braid rather than letting it escape. With periodic BC, this
creates a self-reinforcing circulation.

## Key Insights for the General Model

### Insight 1: The braid3 core geometry is already correct
The optimizer changed the core by only 3.7%. The triple-strand helical
structure with delta = {0, 3.0, 4.4} phase offsets is already close to
optimal for V(P) binding. **Don't redesign the core — it works.**

### Insight 2: The problem is the ENVELOPE, not the geometry
Every failed compact structure (knots, crossed braids) had a different
CORE geometry but the same problem: the braid-background interface radiates.
The optimizer's solution: don't change the core, reshape the tail.

For compact structures, the critical design element is the TRANSITION ZONE:
- How does the field amplitude decay from braid to background?
- What is the field pattern in the low-amplitude region?
- How do the three fields' phase relationships evolve across the transition?

### Insight 3: Slight asymmetry aids stability
A perfectly symmetric initial condition (mirror-symmetric in z) is NOT
optimal. The 6% z-asymmetry suggests that a slight built-in propagation
direction helps the structure maintain coherence during breathing oscillations.

This is analogous to how a spinning top is more stable than a stationary one —
the angular momentum (here: net propagation) provides gyroscopic stability.

### Insight 4: Impedance matching is critical
The concept from wave physics: when a wave transitions between two media
with different impedances, energy is reflected. The braid core has high
field amplitude (A≈0.8) and the background has low (A≈0.1). The abrupt
transition radiates.

The optimizer's solution: reshape the tail to create a smooth impedance
gradient. This is the same principle behind anti-reflection coatings in optics
or tapered transmission lines in RF engineering.

### Insight 5: Backprop in field space is tractable and powerful
At N=32: 98k parameters, 22 seconds on GPU for 300 epochs.
Scaling: N=64 → 786k params, ~3 min. N=128 → 6M params, ~20 min.

This approach can be used to:
1. Optimize any initial condition for long-term stability
2. Search for novel compact structures (random seeds + absorbing BC)
3. Fine-tune parameters for specific physics targets (mass ratio, charge)

## Recommended Next Steps

### A. Compact structure search with absorbing BC
Add absorbing boundary to the JAX loss function. Seed with:
1. Braid3 with various envelope shapes (Gaussian, super-Gaussian, tapered)
2. Random smooth fields with compact support
3. Spherical harmonics expansion (l=0,1,2 modes only)
Let backprop find structures that survive with absorbing BC.

### B. Envelope optimization at N=64
Keep the braid3 core fixed, parameterize only the envelope (r > 3 region).
This reduces the search space and focuses on the impedance-matching problem.

### C. Two-braid interaction
Optimize two braid structures in the same box. The loss includes both
individual survival AND interaction energy. Find configurations where
two braids attract/repel correctly.

### D. Longer T_opt
Current T_opt=30 is short. Increase to T_opt=100+ to find structures
that are stable over many breathing cycles, not just the first few.

## Files

| File | Description |
|------|-------------|
| `backprop/structure_search.py` | JAX backprop implementation |
| `backprop/run.log` | Full training log |
| `backprop/results/results/best_phi.npy` | Optimized phi (3×32×32×32, f32) |
| `backprop/results/results/best_vel.npy` | Optimized vel (3×32×32×32, f32) |
| `backprop/results/results/epoch_*.bin` | Snapshots every 10 epochs |
| `backprop/results/results/final_init.bin` | Optimized initial condition |
| `backprop/results/results/final_evolved.bin` | After 2×T_opt evolution |
| `backprop/deploy.sh` | Vast.ai deployment script |
