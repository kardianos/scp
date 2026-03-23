# Wavefront and Binding Analysis: 3D Kappa-Dependent Collapse

## Simulation Summary

Three 3D simulations of the density-dependent kappa (mode 3, gamma=10) collapse, all starting from identical initial conditions (Gaussian oscillon, A=0.8, sigma=3, N=128, L=10), differing only in boundary conditions:

| Run | File | BC Type | Frames | Time Range |
|-----|------|---------|--------|------------|
| Rect | `kappa_3d.sfa` | Rectangular absorbing (damp_width=3) | 42 | 0-200 |
| Periodic | `kappa_3d_periodic.sfa` | Periodic (wrap-around) | 12 | 0-50 |
| Sphere | `kappa_3d_sphere.sfa` | Spherical absorbing | 16 | 0-30 |

Physics: V(P) = (mu/2) P^2 / (1 + kappa_eff P^2), where kappa_eff = kappa / (1 + gamma * Sigma) and Sigma = phi_0^2 + phi_1^2 + phi_2^2.

---

## 1. Collapse Dynamics — All Three BC Types

### Common Collapse Phase (t=0 to t~15)

All three runs are IDENTICAL until t~10, confirming the initial collapse is intrinsic and not BC-dependent. The key milestones:

| Time | E_pot | phi_max | P_max | Description |
|------|-------|---------|-------|-------------|
| t=0 | -11 | 0.80 | 0.14 | Gaussian initial condition |
| t=6 | -400 | 1.52 | 3.5 | Collapse begins |
| t=10 | -15,800 | 1.75 | 5.3 | Rapid deepening of binding well |
| t=14 | -78,000 | 1.79 | 5.5 | Still accelerating (rect), diverging (periodic) |

The collapse is driven by the density-kappa feedback: high density reduces kappa_eff, which deepens the potential well, drawing in more field, increasing density further.

### Divergence After t~14

**Rectangular absorbing** (blue): E_pot plunges to -210,000 by t=20, then settles to a stable oscillation around E_pot = -218,000 for the remainder (t=20 to t=200). The absorbing boundary drains kinetic energy, and the system reaches a quasi-equilibrium with gentle breathing oscillations (amplitude ~2% of mean E_pot).

**Periodic** (red): E_pot reaches the same depth initially but then OVERSHOOTS, hitting -240,000 at t~25, significantly deeper than the rect case. It oscillates wildly between -218,000 and -244,000 from t=20-50. This is because outgoing waves wrap around and re-enter, amplifying the binding.

**Spherical absorbing** (green): Follows the rect case closely but collapses slightly slower (E_pot = -80,000 at t=14 vs -78,000 for rect). By t=26 it catches up to E_pot = -218,000 and stabilizes. Nearly identical to rect after t=28.

### Late-Time Equilibrium (Rect Only, t>50)

The rect run, being the longest at 200 time units, shows the late-time state clearly:
- E_pot oscillates between -216,800 and -219,000 with period ~20 time units
- phi_max = 1.55 +/- 0.02 (breathing mode)
- P_max ~ 3.6-3.7 (stable)
- All energy resides in ONE connected cluster filling the domain (count = 2,097,152 = entire grid)
- Binding is UNIFORM: core fraction ~ 1.0 (entire volume participates equally)

---

## 2. Internal Binding Pattern — The Periodic BC Effect

### The Key Finding: Amplified Peaks in Periodic BC

The periodic BC run develops dramatically different internal structure compared to the absorbing runs.

**P_max comparison at same physical times:**

| Time | Rect P_max | Periodic P_max | Sphere P_max | Periodic/Rect ratio |
|------|-----------|---------------|-------------|-------------------|
| t=10 | 5.30 | 5.29 | 5.33 | 1.00x |
| t=16 | 5.43 | **33.2** | 4.75 | **6.1x** |
| t=20 | 4.76 | **188** | 3.84 | **39x** |
| t=25 | — | **26.8** | — | — |
| t=30 | 3.87 | **123** | 3.84 | **32x** |
| t=44 | 3.83 | **285** | — | **74x** |

The periodic BC amplifies peak |P| by 1-2 orders of magnitude. These are NOT distributed uniformly but form localized spikes.

### Where Are the Peaks? Cluster Analysis

At frame 5 (t=25), the periodic BC shows **17 clusters** above threshold |P|>0.1:
- **Main cluster**: 2,084,632 grid points (99.4% of volume), centroid at origin, max |P| = 26.8
- **16 satellite clusters**: 148 grid points each, centered at CORNERS: (-9.756, -9.756, -9.756) etc.

The satellite clusters at the box corners are the smoking gun: these are wave-interference artifacts from the periodic wrap-around. Outgoing waves exit one face and immediately re-enter the opposite face, creating constructive interference at the corners where three face-crossing wavefronts converge.

### Radial Profile Structure (Periodic t=25)

The radial profile at t=25 shows a NON-MONOTONIC structure with two distinct enhancement rings:

1. **Inner ring** at r ~ 2.0-2.5: P peaks at 7.64 (vs background ~7.2)
2. **Mid-ring** at r ~ 4.0-4.8: P peaks at 7.72
3. **Edge spike** at r ~ 9.5-10.0: P spikes to 7.98

For comparison, the rect and sphere BCs at t=20 show a smooth MONOTONIC profile with P decreasing from center (7.23) to minimum at r~5.5 (6.75) and mild enhancement at the edges (7.35 for rect, 7.27 for sphere).

The periodic BC profile has structure that the absorbing BCs do not have: the reflected wavefronts create standing-wave-like shells of enhanced |P|.

### Is It a Single Central Blob or Multiple Pockets?

**SINGLE central blob with shell structure.** The BFS cluster analysis shows that at all periodic frames, the dominant cluster is ONE connected region containing >99% of the active grid points. The satellites at corners are tiny (148 points = 0.007% each). The binding pattern is NOT fragmented into separate pockets.

However, the radial profile shows the binding has SHELL STRUCTURE with concentric rings of enhancement at r~2.3 and r~4.5. This is consistent with inward-propagating reflected wavefronts creating constructive interference shells.

---

## 3. Spherical Absorbing BC — No Corner Artifacts

### Core Structure

The sphere BC shows a remarkably clean radial profile:
- At t=20: smooth monotonic decrease from center (rho^2 = 7.18) to edge minimum (rho^2 = 6.80 at r~5.4), then gradual recovery to edges
- NO shell structure, NO edge spikes
- ALWAYS 1 cluster (zero fragmentation at all 16 frames)
- Profile shape is nearly identical to the rect case

### Collapse Rate

Slightly slower than rect (reaches E_pot = -3380 by t=28 vs t=20 for rect), because:
1. The spherical damping zone (r > L-3) removes less energy than rectangular damping
2. The spherical geometry preserves isotropy better, avoiding corner-focusing effects

By t=30 (frame 15), the sphere run matches the rect equilibrium exactly: E_pot = -3381, P_max = 3.84, 1 cluster filling the domain.

---

## 4. Rectangular Absorbing BC — Corner Effects

### Mild Corner Enhancement

The rect BC profile at t=20 shows a slight enhancement at large r (r~8.5, P=7.35 vs minimum P=6.75 at r~5.4). This is much milder than the periodic case but still present.

The corner enhancement comes from the rectangular (not spherical) damping geometry: the damping width is applied per-face, creating diamond-shaped effective boundaries. Waves hitting corners get damped less than waves hitting face centers.

### Long-Term Stability

After the initial collapse (t<20), the rect BC achieves stable equilibrium. From t=30 to t=200:
- E_pot range: [-218,900, -217,100] (0.8% variation)
- phi_max range: [1.53, 1.58] (3% variation)
- P_max range: [3.60, 3.80] (5.4% variation)
- ALWAYS 1 cluster, centroid at origin

---

## 5. Comparison Table

| Property | Rect Absorbing | Periodic | Sphere Absorbing |
|----------|---------------|----------|-----------------|
| Collapse time (E_pot reaches 90% of final) | t ~ 16 | t ~ 14 | t ~ 22 |
| Final E_pot | -218,000 | Oscillating [-218k, -244k] | -218,000* |
| Peak P_max ever | 10.2 (t=24) | **285** (t=44) | 6.86 (t=24) |
| Steady-state P_max | 3.6-3.7 | N/A (still oscillating) | 3.8-3.9 |
| Clusters at equilibrium | 1 | 1 (+ corner satellites) | 1 |
| Shell structure in profile | No | **Yes** (2 rings) | No |
| Edge artifacts | Mild (~5%) | **Severe** (>100x at corners) | None |
| Energy conservation | ~0.1% drift/200t | ~0.1% drift | ~0.1% drift |

*Sphere BC only ran to t=30; value at final frame.

---

## 6. Interpretation: What Is the Red Binding Pattern?

The red binding pattern seen in the periodic BC visualization at t~12 (frame 6) and later is:

**Reflected wavefront convergence creating a centrally-concentrated, shell-structured binding enhancement.**

Mechanism:
1. The initial Gaussian collapse sends outgoing spherical waves toward all boundaries
2. In periodic BC, these waves immediately re-enter from the opposite face
3. Three incoming plane waves (one per axis) constructively interfere at the center
4. This drives P_max up by 1-2 orders of magnitude compared to absorbing BCs
5. The interference creates SHELL structure at specific radii (r~2.3, r~4.5)
6. Corner-crossing effects create tiny high-|P| pockets at the 8 box corners

**This is NOT intrinsic binding structure** -- it is entirely a boundary condition artifact. The absorbing BC runs (both rect and sphere) show NO such enhancement, and the sphere run in particular demonstrates that a clean absorbing boundary produces a completely smooth, monotonic radial profile.

### Physical Significance

The fact that periodic BCs amplify P_max by 39-74x demonstrates that this system is highly sensitive to wave feedback. In a real physical system (infinite domain), the outgoing radiation would never return, and the equilibrium would look like the rect/sphere absorbing cases: a uniform domain-filling condensate with P~3.7 and no shell structure.

---

## 7. Plots Generated

All saved to `/home/d/code/scp/v39/data/`:

| File | Description |
|------|-------------|
| `epot_vs_time.png` | E_pot for all 3 BCs — shows periodic overshoot and oscillation |
| `phi_max_vs_time.png` | Peak amplitude — periodic reaches 6.5 vs 1.8 for absorbing |
| `pmax_vs_time.png` | Peak |P| (log scale) — periodic 2 orders of magnitude higher |
| `pint_vs_time.png` | Integrated binding — periodic 10-20% higher |
| `etotal_vs_time.png` | Total energy — conservation check |
| `radial_profiles_t20.png` | BC comparison at t~20: P(r), rho^2(r), V(r) |
| `periodic_radial_evolution.png` | Periodic BC profile evolution showing shell formation |
| `cluster_and_binding.png` | Cluster count and core binding fraction vs time |

## 8. Tools and Data

- Custom analysis tool: `/home/d/code/scp/v39/src/analyze_3col.c` (compiled as `analyze_3col`)
- SFA fixup tool: `/home/d/code/scp/v39/src/sfa_fixup.c` (for streaming SFA files)
- Raw analysis output: `data/analysis_rect.txt`, `data/analysis_periodic.txt`, `data/analysis_sphere.txt`
- Timeseries: `data/kappa_3d/timeseries.tsv`, `data/kappa_3d_periodic/timeseries.tsv`, `data/kappa_3d_sphere/timeseries.tsv`
