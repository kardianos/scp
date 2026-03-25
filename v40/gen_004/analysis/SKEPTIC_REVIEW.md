# SKEPTIC REVIEW -- Gen 4 Structure Characterization

**Reviewer**: Skeptic Agent
**Status**: FINAL (CHECK 5 of 5) -- No research or implementation output has been produced in the analysis/ directory. Global analysis JSONs exist at the gen_004 level but do NOT answer the TASK.md questions.

---

## CRITICAL: No Output From Other Agents

After four monitoring checks, the `/home/d/code/scp/v40/gen_004/analysis/` directory contains ONLY `TASK.md`. No `RESEARCH.md`, no `.c` files, no analysis results. Either the other agents:

1. Have not been launched
2. Are stuck on an overly complex approach
3. Are running but writing output elsewhere

**Impact**: The entire Gen 4 analysis goal is unmet. No stable-vs-unstable characterization exists.

**Recommendation**: If this persists, the analysis must be built from scratch. The existing tools (`analyze_sfa`, `spatial_analysis`) compute only GLOBAL metrics and are insufficient for the per-cluster, per-region analysis that TASK.md demands.

---

## Pre-Emptive Review: Issues With Existing Tools

Even if the other agents eventually produce output, they will likely build on the existing `analyze_sfa.c` and `spatial_analysis.c`. Both have significant gaps.

### Issue 1: BFS Cluster Detection Uses Relative Threshold

**What's wrong**: `spatial_analysis.c` line 164 sets `threshold = 0.01 * Pmax`. This is a relative threshold that changes frame-to-frame as P_max oscillates. At T=10, S20's P_max drops (breathing minimum from Gen3 analysis shows P_max going from 2.34 to 0.84). A 1% threshold on P_max=0.84 is 0.0084 -- much lower than 1% of the initial P_max=2.34 (=0.023). This means the BFS will find FEWER clusters during breathing maxima (when everything is connected) and MORE during breathing minima (when connections break). The cluster count is dominated by breathing dynamics, not by actual structural fragmentation.

**Why it matters**: Comparing cluster counts across frames is meaningless if the threshold is not consistent. You cannot distinguish "structure fragmented" from "structure is at breathing minimum" without a fixed threshold.

**What to do instead**: Use a FIXED absolute threshold based on the background |P| level. The background has A_bg=0.1, so P_bg ~ 0.1^3 = 0.001. A threshold of 5x to 10x the background (0.005 to 0.01) would be physically motivated and frame-independent.

**Severity**: MAJOR -- weakens all cluster evolution conclusions in GEN2_SUMMARY.

### Issue 2: No Per-Cluster Metrics

**What's wrong**: `spatial_analysis.c::count_clusters()` returns only the NUMBER of clusters. It does not compute per-cluster mass, centroid, R_rms, aspect ratio, P_int, theta_rms, or energy. Without per-cluster metrics, you cannot track individual clusters across time, determine which clusters are stable vs unstable, or measure inter-cluster interactions.

**Why it matters**: The ENTIRE POINT of the Gen 4 analysis (per TASK.md questions 1-5) is to compare properties of STABLE regions vs UNSTABLE regions. You need to track individual clusters through time and measure their properties separately.

**What to do instead**: The BFS should label each voxel with its cluster ID. Then compute all field metrics (rho, |P|, theta_rms, velocity divergence, energy) per cluster. Output per-cluster JSON.

**Severity**: CRITICAL -- blocks the primary goal of Gen 4 analysis.

### Issue 3: Velocity and Acceleration Fields Are Not Analyzed

**What's wrong**: `spatial_analysis.c` does not use velocity data at all. `analyze_sfa.c` computes E_kin from velocities but does not compute velocity divergence (div(v)), velocity structure, or any kinematic analysis.

**Why it matters**: TASK.md Question 3 asks: "Are surviving clusters expanding or contracting (div(v) sign)?" and "Is KE localized in cores or shells?" Neither existing tool can answer this.

The velocity field contains crucial physics:
- div(v) > 0 means the region is expanding (dispersing)
- div(v) < 0 means the region is contracting (binding)
- |curl(v)| reveals rotational structure (helical braid dynamics)
- v dot r_hat (radial velocity from cluster center) shows systematic expansion/contraction

**What to do instead**: Compute div(phi_v) = d(phi_vx)/dx + d(phi_vy)/dy + d(phi_vz)/dz at each grid point. Then bin by cluster membership. Report per-cluster mean div(v), std div(v), and the fraction of volume with div(v) < 0.

**Severity**: MAJOR -- TASK.md explicitly asks for this.

### Issue 4: No Acceleration Analysis (TASK.md Q4)

**What's wrong**: The SFA files contain 12 columns: phi (3), theta (3), phi_vel (3), theta_vel (3). There are NO acceleration columns. TASK.md Question 4 asks about acceleration fields and force balance.

**Why it matters**: To answer Q4, the analysis tool must COMPUTE the acceleration from the equation of motion:

    a_phi = laplacian(phi) - m^2 * phi - dV/dphi + eta * curl(theta)

This requires computing the Laplacian and curl from the grid data, plus the V'(P) force. This is the MOST expensive per-voxel computation and the MOST physically informative. Without it, you cannot determine:
- Where forces are strongest (core vs surface)
- Whether V(P) force or curl coupling dominates
- Whether stable regions have balanced forces

**What to do instead**: Implement force decomposition: f_laplacian = laplacian(phi), f_mass = -m^2*phi, f_V = -dV/dphi, f_curl = eta*curl(theta). Report each component's magnitude per cluster. The net force |f_total| should be near zero in truly stable regions (dynamic equilibrium).

**Severity**: MAJOR -- TASK.md explicitly asks for this. Requires significant computation.

### Issue 5: No Phase Coherence Analysis

**What's wrong**: Neither existing tool measures phase relationships between the three phi fields. The triple product P = phi_0 * phi_1 * phi_2 is computed, but the PHASE STRUCTURE that determines P is not analyzed.

**Why it matters**: Per CONCEPT.md and POSTULATES.md, the braid's stability comes from specific phase offsets delta = {0, 3.0, 4.43} that maximize time-averaged |P|. If phase coherence breaks down in a region, |P| drops and binding is lost. Measuring phase coherence could PREDICT which regions will fragment before they actually do.

**What to do instead**: For each grid point or cluster, compute the local phase via:
- theta_a(x) = atan2(phi_a(x), phi_a(x+dz)) -- phase of the z-oscillation
- Phase coherence: C = |<exp(i*(theta_0 - theta_1 - theta_2))>| over the cluster
- C = 1 means perfect phase locking (strong P), C ~ 0 means random phases (no binding)

Alternatively, compute the Fourier transform along z for each phi component and measure the coherence of the dominant k_z mode.

**Severity**: MAJOR -- this is the most physically insightful metric missing. Phase coherence is the MECHANISM of binding.

### Issue 6: No Radial Profile Computation

**What's wrong**: Neither tool computes radial profiles of ANY field quantity around cluster centroids. TASK.md Question 1 specifically asks for "radial profile of rho around cluster centroids vs inter-cluster voids."

**Why it matters**: Radial profiles reveal the internal structure of each cluster -- whether it has a well-defined core with a depletion halo, or is diffuse. The comparison between stable and unstable clusters' radial profiles is the most direct way to answer "what makes a region stable?"

**What to do instead**: For each cluster centroid, compute azimuthally-averaged radial profiles of:
- rho(r) = sum(phi_a^2) binned by distance from centroid
- |P|(r) binned by distance
- theta_rms(r) binned by distance
- v_radial(r) = v dot r_hat binned by distance
- Use at least 30 radial bins out to R_max = 2 * R_rms of the cluster.

**Severity**: CRITICAL -- this is the most basic characterization requested and is entirely missing.

### Issue 7: Inter-Cluster Region Is Ignored

**What's wrong**: All analysis focuses on cluster interiors (where |P| > threshold). The space BETWEEN clusters is never analyzed.

**Why it matters**: Per POSTULATES.md, the interaction surface at r = 4-6 code units is "WHERE the braid talks to the fabric." Inter-cluster regions determine whether two clusters will merge or repel. Understanding the field properties between clusters is essential for understanding binding.

**What to do instead**: Define inter-cluster regions as voxels where |P| < threshold but the nearest cluster centroid is within 2*R_rms. Compute field statistics in these regions separately. Look for bridging structures (|P| > background but below cluster threshold) that connect clusters.

**Severity**: MAJOR -- needed for understanding composite formation.

### Issue 8: Theta Variance Computation Is Incorrect

**What's wrong**: In `spatial_analysis.c` lines 133-142, theta "variance" is computed as `trms2 - trms*trms` where `trms = sqrt(mean(theta^2))` and `trms2 = sqrt(mean(theta^4))`. This is NOT the variance of theta. The variance of theta would be `mean(theta^2) - mean(theta)^2`, which for a zero-mean field equals `mean(theta^2)` = `trms^2`. The quantity `trms2 - trms^2` is the difference between the L4 norm and the squared L2 norm, which measures kurtosis (peakedness), not spread.

**Why it matters**: The "tvar" column in the spatial analysis output is labeled as variance but measures something different. Any conclusions drawn from it about theta stability are unreliable.

**What to do instead**: If you want the spatial variation of theta field strength, compute the per-voxel theta energy density `e_theta = 0.5*(theta_x^2 + theta_y^2 + theta_z^2)` and then report `std(e_theta) / mean(e_theta)` as a coefficient of variation.

**Severity**: MINOR -- the quantity is still informative (kurtosis proxy), but should be labeled correctly.

### Issue 9: F16 Data Precision

**What's wrong**: The `analyze_sfa.c` f16 decoder (lines 171-176) handles normal f16 values but treats ALL subnormal values as zero (`if(e==0){arr[i]=0;}`). F16 subnormals cover the range ~5.96e-8 to 6.10e-5. The background field has P_bg ~ 0.001, which is well within f16 range, but fine details of the depletion halo (delta_rho ~ 1e-3 to 1e-5) may be in the subnormal range and would be zeroed.

**Why it matters**: The f16 files exist for quick analysis (`S20_f16.sfa` at 841 MB vs 11 GB). If the implementation agent uses f16 files for speed, the depletion halo and inter-cluster field properties will be inaccurate.

**What to do instead**: Either use the full f64 files (slow but accurate) or implement proper f16 subnormal decoding. The subnormal formula is: `value = (-1)^s * (m/1024) * 2^(-14)`. This is a one-line fix in the decoder.

**Severity**: MINOR if f64 files are used. MAJOR if f16 files are used for detailed analysis.

### Issue 10: Duplicate Final Frames in Gen 3 Data

**What's wrong**: The Gen 3 analysis.json files show duplicate final frames (two entries at t=30.00 with identical data). This is likely a bug in the SFA writer or the frame counting. The task says Gen 4 has 42 frames -- if the same duplication happens, there may be 41 unique frames.

**Why it matters**: Averaging over frames with duplicates skews time series statistics. If the last frame is duplicated N times, it gets N/(N+real_frames) weight in means.

**What to do instead**: Check for and deduplicate frames by time value before averaging.

**Severity**: MINOR -- affects statistical summaries slightly.

---

## What the Implementation MUST Include

For the Gen 4 analysis to actually answer the questions in TASK.md, the implementation needs AT MINIMUM:

### Tier 1 (Blocks conclusions without these):
1. **Per-cluster labeling BFS** with fixed absolute threshold
2. **Per-cluster metrics**: mass, centroid, R_rms, P_int, theta_rms, E_pot
3. **Radial profiles** around each cluster centroid (rho, |P|, theta, v_radial)
4. **Cluster tracking** across frames (nearest-centroid matching)

### Tier 2 (Weakens analysis without these):
5. **Velocity divergence** div(v) per cluster
6. **Force decomposition** (Laplacian, mass, V(P), curl) per cluster
7. **Phase coherence** metric per cluster
8. **Inter-cluster region characterization**

### Tier 3 (Nice to have):
9. **Energy flow** tracking (dE_kin/dt, dE_pot/dt per cluster)
10. **Fourier analysis** of phase structure along dominant axis
11. **Stability prediction** from t=0 data

If the implementation agent delivers only global averages (what the existing tools already compute), the Gen 4 analysis has added ZERO new information beyond Gen 3.

---

## Physics Insights the Research Agent Should Address

1. **Derrick's theorem locally**: For each cluster, compute E_2 (gradient) and E_4 (would be the quartic Skyrme term, but here it is E_V from V(P)). In the Skyrme model, E_2 = E_4 at equilibrium. If a cluster has E_2 >> E_4, it is over-extended and will contract. If E_4 >> E_2, it is too tight and will expand. This virial ratio is a direct stability predictor.

2. **The curl(theta) coupling structure**: theta is sourced by curl(phi). In a helical braid, curl(phi) has a specific structure (helical flow along the braid axis). If the theta field develops this coherent structure, the curl(theta) back-reaction on phi is STABILIZING. If theta is random, the back-reaction is noise. Measuring the alignment between curl(phi) and theta is a direct probe of electromagnetic coupling strength.

3. **P = phi_0 * phi_1 * phi_2 has NODES**: Where any one phi_a passes through zero, P = 0 regardless of the other fields. The NODE SURFACES of each phi_a partition space into binding regions and non-binding regions. The topology of these node surfaces (are they closed? do they separate clusters?) may determine stability.

4. **Energy conservation as a diagnostic**: If E_total is not conserved within 0.1% between frames, the absorbing boundary is actively draining energy. The RATE of energy loss per cluster (not just global) reveals which clusters are radiating more -- a direct stability measure.

5. **Background amplitude matters**: A_bg = 0.1 creates P_bg ~ 0.001. The cluster |P| is 10-100x this. But the RATIO P_cluster / P_bg may be more meaningful than P_cluster alone for stability analysis.

---

## Summary Assessment

| Category | Status | Notes |
|----------|--------|-------|
| Research output | NOT STARTED | No RESEARCH.md file exists |
| Implementation | NOT STARTED | No .c files in analysis/ |
| Results | NONE | No JSON, no analysis output |
| Existing tools | INSUFFICIENT | Global metrics only, no per-cluster |
| Data availability | OK | Three 11 GB SFA files present |

**Overall**: The Gen 4 analysis has not begun. When it does begin, the implementation must go SIGNIFICANTLY beyond the existing `analyze_sfa` and `spatial_analysis` tools to answer the questions posed in TASK.md. Global averages are not enough. Per-cluster, per-region analysis with radial profiles and velocity/force decomposition is the minimum useful output.

---

---

## Data-Driven Findings from Existing Global Analysis

The gen_004/ directory contains global analysis JSONs (S20_analysis.json, UDD_R4_analysis.json, CB15_analysis.json) produced by the existing `analyze_sfa` tool. Reviewing these reveals additional issues that the detailed analysis MUST address.

### Issue 11: UDD_R4 Has Stale/Duplicated Frames After t=190

**What's wrong**: In UDD_R4_analysis.json, frames at t=190, t=195, t=200, and t=200 ALL have identical data (E=3458.8, Ep=-391.87, Pint=353.48). That is 4 consecutive frames with the same values. This means either:
- The SFA writer produced duplicate frames after t=190
- The simulation stalled and produced no new data after t=190
- The output was truncated and backfilled

**Why it matters**: UDD_R4's "survival" at T=200 may be an artifact. If the simulation actually stalled at t=190, the last real data point is t=190 (not t=200). The P_int retention of 75% is computed from the duplicate t=200 value. The TRUE retention at the last real frame (t=185) is P_int=182.5/471.3 = 38.7%, dramatically worse than the reported 75%.

More fundamentally: if the simulation was producing physically correct output at t=190 with Ep=-391.9 and Pint=353.5, why would the next 3 frames be identical? This smells like a data pipeline bug, not real physics.

**What to do instead**: Check the UDD_R4_output.sfa file directly. Read frames 39-42 and compare raw field data (not just global metrics). If the raw voxel data is identical across these frames, the SFA file is corrupted/truncated. Report the LAST UNIQUE frame as the true endpoint.

**Severity**: CRITICAL -- invalidates the UDD_R4 survival claim and P_int retention metric.

### Issue 12: Strong Oscillatory Behavior Masks Trends

**What's wrong**: All three structures show large-amplitude oscillations in E_pot and P_int. For example, S20 oscillates between E_pot=-97 (t=135) and E_pot=-961 (t=45) -- a 10x range. CB15 oscillates between E_pot=-134 (t=200) and E_pot=-962 (t=10). These are breathing modes with period ~25-30 time units.

**Why it matters**: The S_final metric (score at the last frame) is dominated by WHEN the simulation stopped relative to the breathing cycle. CB15's S_final=0.82 is low because T=200 caught a trough. If the sim had stopped at T=195 or T=205, CB15 would score very differently. This makes S_final nearly useless as a stability metric.

**What to do instead**:
1. Report TIME-AVERAGED metrics over the last ~50 time units (2 breathing periods), not instantaneous values.
2. Report the ENVELOPE of the oscillation (E_pot_max, E_pot_min over the last 50 time units).
3. Report the breathing period and amplitude explicitly.
4. Use the minimum E_pot (deepest trough) in the last 50 time units as the "worst case" stability bound.

S20 last 50 time units (t=155-200): E_pot ranges from -94 to -507, mean ~ -300.
CB15 last 50 time units: E_pot ranges from -134 to -665, mean ~ -420.
UDD_R4 last 50 time units (using only unique frames up to t=190): E_pot ranges from -22 to -567, mean ~ -290.

With time-averaging, CB15 actually has the STRONGEST binding, not the weakest. The S_final ranking is misleading.

**Severity**: MAJOR -- the current ranking (S20 > UDD_R4 > CB15) is an artifact of breathing phase at T=200.

### Issue 13: Energy Loss Rate Differs Significantly Between Structures

**What's wrong**: The total energy E_total decreases steadily for all three structures (absorbing BC draining radiation). But the rates are very different:

- S20: E drops from 17420 to 3449, losing 80% of energy over T=200
- CB15: E drops from 16061 to 5273, losing 67% of energy over T=200
- UDD_R4: E drops from 12414 to 3459, losing 72% of energy over T=200

CB15 loses the LEAST energy (67%) despite starting with the most. This means CB15 radiates less energy per unit time -- it is a MORE EFFICIENT structure (less leaky).

**Why it matters**: If the analysis only looks at the final score, it misses this crucial difference. The rate of energy loss is arguably the BEST single predictor of long-term stability. A structure that radiates 80% of its energy in T=200 will be dead by T=400. A structure that radiates 67% may stabilize.

**What to do instead**: Compute dE/dt as a function of time. Fit to an exponential decay E(t) = E_inf + (E_0 - E_inf) * exp(-t/tau). The decay constant tau tells you the structure's radiation timescale. A larger tau means longer-lived.

**Severity**: MAJOR -- energy loss rate is a primary stability discriminant not currently reported.

### Issue 14: The "alive" Check Is Too Lenient

**What's wrong**: In analyze_sfa.c line 113: `m->alive = (m->P_int > 0.1 && m->E_pot < 0 && m->phi_max > 0.2)`. With P_int values ranging 40-960, a threshold of 0.1 is trivially satisfied. With phi_max ranging 0.53-1.70 and background ~0.1, a threshold of 0.2 is only 2x background. All 42 frames for all 3 structures pass this test. The "alive" metric is useless -- it never triggers.

**Why it matters**: The POSTULATES.md defines death criteria that are much more stringent: "P_int drops below 10% of initial" and "phi_max drops below 2x background." With P_int(0)=960 for S20, the threshold should be 96, not 0.1. With background phi ~0.2 (A_bg=0.1 * cos), phi_max threshold should be 0.4.

**What to do instead**: Use relative thresholds tied to initial conditions: P_int_threshold = 0.1 * P_int(0), phi_max_threshold = 2 * A_bg * sqrt(3).

**Severity**: MINOR -- does not affect the current conclusion (all survive) but would matter for weaker structures.

---

## Final Summary

After 5 monitoring checks, the Gen 4 DETAILED analysis (per-cluster, per-region characterization as described in TASK.md) has NOT been started. The only analysis available is the global metrics from analyze_sfa, which:

- Cannot answer any of the 5 questions in TASK.md
- Has a misleading S_final metric dominated by breathing phase
- Has possible data corruption in UDD_R4 (duplicate frames)
- Does not compute velocity, acceleration, phase, or radial profile information
- Does not distinguish stable vs unstable regions

The global metrics do establish that all three structures survived T=200 (which is useful), but provide NO insight into WHY they survived or WHAT distinguishes their internal structure. This is the gap that the Gen 4 analysis was supposed to fill.
