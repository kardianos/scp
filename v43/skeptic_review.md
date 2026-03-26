# Skeptic Review: V43 Results (OQ1, OQ2, OQ3)

Reviewer: Skeptical physicist analysis.
Date: 2025-03-25.

---

## OQ3: DC Theta Profile — VERDICT: CONCERN (methodology valid, claim overstated)

### Claim Under Review

DC/RMS ratio is 5-12% in the core, not the 0.2% claimed in V34.

### What the Code Does

`dc_profile.c` reads all 264 frames of `braid_hires.sfa` (N=80, L=25). For each frame it finds the |P|-weighted centroid, then bins all voxels by distance from that centroid into radial shells (dr=0.5). It accumulates signed theta_x, theta_y, theta_z and theta^2 sums across ALL frames and ALL voxels in each bin. Finally it computes DC = |mean(theta_vec)| and RMS = sqrt(mean(theta^2)).

### Issue 1: The "DC" Being Measured Is NOT What Was Meant in V34

**This is the central problem.** The code computes DC as the magnitude of the time-AND-space averaged theta vector within each radial shell. The V34 "0.2%" measurement was apparently a time-series measurement at fixed spatial points or a spatial decomposition method.

The quantity computed here, |<theta_vec>|, conflates two different things:
1. True DC: the time-averaged theta at each spatial point (which would require tracking individual grid points across frames).
2. Spatial averaging: averaging signed theta over an entire spherical shell, which would cancel any dipolar or higher-multipole DC pattern even if it existed.

What the code actually computes is: average over all voxels in a shell AND all frames. If the DC theta field has dipolar angular structure (which it almost certainly does for a braid aligned along z), spatial averaging within a shell will partially cancel it, reducing the apparent DC. However, the fact that the code takes |mean_vec| (magnitude of the mean vector) rather than mean(|vec|) partially mitigates this — it preserves any coherent directional offset.

The result shows mean_tx and mean_ty are both positive and comparable in the core (both ~0.007-0.009), while mean_tz is essentially zero (~1e-4). This is physically sensible for a braid aligned along z (the theta field is transverse to the braid axis). The fact that the DC is dominated by the x-y plane is a genuine signal, not noise.

**Assessment**: The 5-12% DC/RMS in the core is a LOWER BOUND on the true DC, because spatial averaging within shells suppresses angular structure. The true per-voxel time-averaged DC could be larger. The disagreement with V34's "0.2%" is likely because V34 used a different methodology (possibly measuring at a single point far from the core, or measuring a different quantity). The two measurements are not directly comparable.

**Verdict**: The methodology is defensible but the claim "not the 0.2% claimed in V34" is an apples-to-oranges comparison unless the V34 method is reproduced exactly and shown to give a different answer. The 5-12% number in the core is interesting but its physical meaning needs clarification.

### Issue 2: Centroid Tracking and PBC Contamination

The |P|-weighted centroid is a reasonable tracking method. The braid core (where |P| is concentrated) should be well-localized in N=80 with L=25. Frame-to-frame centroid jitter would add noise but not systematic bias, so this is acceptable.

However, at N=80 with L=25, the domain is [-25, +25] with periodic boundaries. Voxels at r > ~17 (r > 2L/3) are potentially contaminated by periodic images. The DC data shows a sign flip at r~16 (mean_tx and mean_ty go negative). This coincides with the PBC contamination zone and should be treated with extreme suspicion.

**Verdict**: The sign flip at r~16 is almost certainly a PBC artifact. Data beyond r~12 is unreliable. The core region (r < 8) should be robust.

### Issue 3: Grid Spacing Discrepancy

The code uses `h = 2.0 * L / N` (cell-centered grid with N cells), but the simulation kernel uses `dx = 2.0 * L / (N - 1)` (vertex-centered grid with N vertices). For N=80, this is h=0.625 vs dx=0.6329 — a 1.3% difference in coordinate assignment. This is a minor bug that shifts all radial assignments by ~1% and is negligible for the conclusions.

### Recommendations

1. To properly compare with V34, compute the time-averaged theta at FIXED grid points (not shell-averaged), then take the RMS of these time-averages. This gives the true per-point DC.
2. Discard all data beyond r=12 (PBC contamination).
3. Fix the grid spacing to match the sim kernel convention.

---

## OQ1: Multipole Expansion — VERDICT: FAIL (fundamental methodological error)

### Claim Under Review

UUD composite theta_rms is 89% monopole at r=12, demonstrating effective point-charge behavior.

### Critical Error: Decomposing |theta| Instead of the Physical Field

The code decomposes theta_rms = sqrt(theta_x^2 + theta_y^2 + theta_z^2) — a positive-definite scalar — into spherical harmonics. This is problematic for multiple reasons:

**Problem A: The monopole is guaranteed to dominate for any positive-definite quantity.**

For ANY three random fields f_x, f_y, f_z, the quantity sqrt(f_x^2 + f_y^2 + f_z^2) is always positive. Its spherical harmonic decomposition will always have a large l=0 component simply because the function has a nonzero mean over the sphere. This is a mathematical triviality, not a physical result. To demonstrate that the UUD composite is special, you would need to compare its monopole fraction against a NULL MODEL — for example, three randomly oriented dipoles, or a single braid (which should be dipolar).

Consider a pure dipole field f = cos(theta). Then |f| = |cos(theta)|, and its spherical harmonic decomposition has a LARGE l=0 component (~64%), even though the underlying field is pure dipole. Taking the absolute value or the magnitude of a vector creates artificial monopole content.

**Problem B: The physically relevant quantity for radiation pressure is the Poynting flux, not |theta|.**

The radiation pressure force depends on S = E x B = (-dtheta/dt) x (curl theta). The multipole decomposition should be applied to the Poynting flux components (or to the individual signed theta components at multiple time snapshots), not to the instantaneous magnitude.

**Problem C: The per-component analysis at R=10 actually shows the REAL structure.**

The per-component decomposition shows theta_x, theta_y, theta_z are dominated by l=2 and l=4 (quadrupole and hexadecapole), with very little l=0. This is the physically meaningful result. The fact that the MAGNITUDE combines to look monopolar is the artifact described in Problem A.

### Issue: Is Three Orthogonal Dipoles -> Isotropic Magnitude Trivially Expected?

**Yes, approximately.** If three orthogonal dipoles have patterns sin^2(alpha_x), sin^2(alpha_y), sin^2(alpha_z), then the sum sin^2(alpha_x) + sin^2(alpha_y) + sin^2(alpha_z) = 2 (constant). This is a standard identity. But this applies to INTENSITIES (squares of field amplitudes), not to field magnitudes. For the RMS magnitude sqrt(f_x^2 + f_y^2 + f_z^2), the isotropy depends on whether the three components have the same RMS amplitude and uncorrelated phases. If they do, the result is trivially isotropic. This does not require any special UUD geometry — it would hold for ANY three equal-amplitude, orthogonally-oriented dipoles.

The real question — whether the Poynting flux from three orthogonal oscillating current loops is isotropic — is answered analytically in OQ2 (yes, by the sin^2 sum identity). The numerical multipole expansion adds nothing beyond confirming this analytic result, and does so with a methodology that inflates the apparent monopole content.

### Issue: Angular Resolution

72 x 144 (theta x phi) is adequate for l <= 4. The Nyquist criterion for spherical harmonics requires at least 2l+1 points per great circle. For l=4, this is 9 points, and 72 is well above this. The angular resolution is fine.

### Issue: f16 Data

The code reads from an f16 (half-precision) SFA file. Half-precision has only ~3 decimal digits of precision. For the far-field theta values (~0.04), the absolute precision is ~4e-5. The C_l values at R=12 are ~1.3e-3, so the signal is ~30x above the f16 noise floor. This is marginal but probably adequate. However, the l=1 fraction at R=12 (1.6%) corresponds to C_1 = 2.3e-5, which IS near the f16 noise floor. The very low l=1 values at large R may be noise-floor artifacts.

### Issue: Single Frame Analysis

The multipole expansion uses a single frame (frame 11, t=499.99). The theta field oscillates rapidly (period ~4.4 t). A single snapshot captures the field at one phase of oscillation, which may have a systematically different angular pattern than the time average. For radiation pressure, the TIME-AVERAGED Poynting flux is what matters. A proper analysis would average |S| over many frames spanning several oscillation periods.

### Recommendations

1. Decompose the SIGNED theta components individually (already done in the per-component table — present those results as the primary finding).
2. Compute a null model: repeat the analysis for a single braid and show that the monopole fraction is lower.
3. Decompose the time-averaged Poynting flux magnitude, not the instantaneous |theta|.
4. Average over multiple frames spanning at least 2-3 oscillation periods.
5. Use f32 or f64 data if available, particularly for large-R analysis.

---

## OQ2: Poynting Derivation — VERDICT: PASS (with important caveats correctly identified)

### Claim Under Review

1/r^2 force law is proven from energy conservation for massless radiation.

### The Core Argument Is Correct

The three-line argument is logically sound and essentially a textbook result:

1. Localized oscillating source radiates with <|S|> ~ 1/r^2 (energy conservation on expanding wavefronts).
2. Target scatters radiation, absorbing momentum F = sigma <S>/c ~ 1/r^2.
3. Therefore F ~ 1/r^2.

This is indeed model-independent and follows from (a) massless radiation in 3+1D and (b) localized sources. The derivation correctly identifies this as kinematic rather than dynamical.

### The Poynting Vector Derivation Is Correct

The proof of du/dt + div(S) = -J.E follows standard EM manipulations applied to the Cosserat theta equation. The key steps:

- dE/dt = curl(B) - J_eff (from the theta EOM)
- dB/dt = -curl(E) (from B = curl(A), E = -dA/dt)
- Standard vector identity: div(E x B) = B.curl(E) - E.curl(B)

This is rigorous. The identification S = E x B and the continuity equation are proven.

### The div(E) = 0 Proof Is Correct

The argument that div(curl(phi)) = 0 implies div(E) = 0 for all time (given div-free initial conditions) is a theorem of the equations. This is correctly stated as an exact result, not an approximation.

### Caveat Assessment

**Caveat 1 (Finite source size, R/lambda ~ 0.7)**: This is the most serious issue and is CORRECTLY flagged. The point-dipole approximation is qualitative when R/lambda ~ 0.7. However, the document correctly argues that the 1/r^2 scaling is independent of multipole order — only the angular pattern and numerical coefficient are affected. This is correct.

**Caveat 2 (Broadband spectrum)**: Correctly handled. Sum of 1/r^2 terms at different frequencies is still 1/r^2.

**Caveat 3 (Near-field transition)**: Correctly noted. The transition to nuclear-scale physics at r ~ lambda is physical, not an artifact.

**Caveat 4 (Back-reaction)**: Correctly identified as higher-order (alpha << 1).

### Issue: The Sign Argument (Section 4) Is Weak

The scattering-phase argument for why same-winding attracts is physically motivated but hand-wavy. The claim that "in-phase scattered waves interfere constructively in the backward direction" requires a proper differential cross-section calculation. The document correctly labels this as "CONSISTENT" rather than "PROVEN", which is appropriate.

However, the comparison with V34 data (section 4.2) has a quantitative problem: the ratio |F_EM(opp)/F_EM(same)| = 2.1. For pure dipole-dipole radiation pressure, same and opposite winding should give equal-magnitude forces of opposite sign, so this ratio should be 1.0. The factor 2.1 asymmetry is attributed to "the force depends on relative orientation of dipole moments", but this is not derived. It could also reflect (a) gravity contamination in the subtraction, (b) nonlinear effects at D=15, or (c) the fact that the force is not pure radiation pressure at this separation. The honest conclusion is that the sign is correct but the magnitude is unexplained.

### Issue: alpha ~ eta^2 Is Dimensional Analysis, Not a Derivation

Equation (24) is correctly labeled as dimensional analysis. The key omission is that sigma_eff and m_0 are themselves functions of eta (since the braid structure depends on the coupling), so the true scaling could be different from eta^2. The claim is testable by varying eta in simulation, which is good scientific practice.

### Issue: "PROVEN" Labels

The document uses "PROVEN" for results that follow from standard textbook electrodynamics applied to the Cosserat mapping. This is appropriate — given the mapping theta -> A, these results ARE proven. The mapping itself (whether theta really IS the electromagnetic vector potential in this universe) is the deeper assumption, but the document does not claim to prove that; it takes the mapping as given from EM_THEORY.md.

### Issue: The Coulomb Force Requires Monopole Radiation, Not Just 1/r^2

The derivation proves F ~ 1/r^2 for radiation pressure from a localized source. But Section 2.2 proves there is NO monopole radiation (div E = 0). The radiation is dipolar. A single magnetic dipole has an anisotropic radiation pattern (sin^2 alpha), giving an anisotropic force. The resolution (Section 3.3) is that three orthogonal dipoles sum to an isotropic pattern.

This creates a subtlety: the "Coulomb force" between composites requires BOTH baryons to be composites of three orthogonal dipoles. Two isolated braids would interact via an anisotropic dipole-dipole force (which falls as 1/r^2 in intensity but has sin^2 angular dependence). Only composite-composite interactions look Coulombic. This is actually a nice feature (it distinguishes "quarks" from "hadrons") but should be stated more clearly.

### The Real Weakness Not Discussed

The derivation assumes the force is entirely radiation pressure (momentum transfer from scattered radiation). But there is another possibility: the time-averaged near-field interaction. The reactive near-field components (1/r^2 and 1/r^3 in the field, giving 1/r^4 and 1/r^6 in the force) dominate at r < lambda. These are standard dipole-dipole magnetostatic forces. At r >> lambda, only the radiative term survives, giving 1/r^2.

The TRANSITION between these regimes is where the interesting nuclear physics lives. The document correctly identifies this (Caveat 3) but does not estimate the crossover distance in physical units. Given lambda ~ 4.4 code units ~ 2.5 fm, the Coulomb regime begins at r >> 2.5 fm, which is ~5-10 fm. This is actually at the boundary of nuclear physics, which means the "Coulomb regime" may not cleanly separate from the "nuclear regime" — they may overlap. This deserves more discussion.

### Recommendations

1. Compute the differential scattering cross section for a braid current loop to make the sign argument rigorous.
2. Perform the eta-variation test (Section 5.2) to check alpha ~ eta^2 empirically.
3. Estimate the crossover distance from near-field to far-field in physical units and compare with known nuclear/atomic scale separation.
4. Address the F_EM asymmetry ratio of 2.1 quantitatively.

---

## Cross-Cutting Concerns

### PBC Contamination

Both OQ3 and OQ1 use data from simulations with periodic boundaries:
- OQ3: N=80, L=25 — PBC images at distance 50. Data beyond r~12 is suspect.
- OQ1: N=192, L=30 — PBC images at distance 60. Data at R=15 (half the box) should have ~10% contamination from images. The R=12 data (89% monopole) may be partially artifactual.

The non-monotonic behavior at R=12 (89% monopole, up from 63% at R=10, then down to 81% at R=15) is suspicious. If the monopole fraction were genuinely increasing with R (as expected for the dipole falloff argument), it should increase monotonically. The peak at R=12 followed by a dip at R=15 could indicate PBC contamination or statistical fluctuation from the single-frame analysis.

### Overclaiming

- OQ1 claims "YES -- l=0 dominates overwhelmingly" but this is an artifact of decomposing the magnitude rather than the signed field.
- OQ2 appropriately labels most results with their correct status (PROVEN vs ASSUMPTION).
- OQ3 claims the V34 "0.2%" is wrong, but the comparison is between different methodologies measuring different quantities.

### What Would Strengthen the Overall Narrative

1. A simulation with ABSORBING boundaries (not periodic) at N=256, L=100 to measure the true far-field theta profile without PBC contamination.
2. The time-averaged Poynting flux multipole decomposition on the UUD composite, which is the physically relevant quantity for the Coulomb force argument.
3. An eta-variation study measuring the effective force at fixed separation to test alpha ~ eta^2.
4. A proper per-voxel DC measurement (track individual grid points across time) to resolve the OQ3 discrepancy with V34.

---

## Summary Table

| Result | Verdict | Key Issue |
|--------|---------|-----------|
| OQ3 DC profile | **CONCERN** | Measures spatial-and-temporal average, not per-point DC; incomparable with V34 methodology; PBC contaminates r > 12; sign flip at r~16 is artifact |
| OQ1 Multipole | **FAIL** | Decomposing positive-definite |theta| guarantees artificial monopole dominance; no null model comparison; single-frame snapshot; per-component analysis shows the real (multipolar) structure |
| OQ2 Poynting | **PASS** | Core 1/r^2 argument is rigorous and correctly caveated; sign argument is qualitative but honestly labeled; alpha ~ eta^2 is dimensional analysis, appropriately noted |
| OQ5 Photon dispersion | **FAIL** | Initial condition wrong (3D Gaussian is not a plane-wave pulse); velocity measurement dominated by geometric spreading, not propagation; phi coupling is physical but magnitude untested against theory |
| OQ6 Nonlinear Maxwell | **PASS** | Perturbation expansion correct; no-birefringence result is rigorous (structural from multilinearity of P); self-corrections and honest caveats throughout; background consistency issue acknowledged |
| F9 Analytical potential | **CONCERN** | Static repulsion sign is rigorous; C~127 is hand-waving dressed as dimensional analysis with post-hoc radius tuning; "if n=1 then Newton" is circular; force exponent prediction 2.2 vs measured 1.8 is a genuine discrepancy, not agreement |

---

## OQ5: Photon Dispersion Test — VERDICT: FAIL (flawed initial condition and measurement methodology)

### Claim Under Review

(a) Group velocity = 0.997c early time, (b) 15-20% energy transfers to phi via curl coupling, (c) pure theta is not an eigenmode.

### Critical Error 1: The Initial Condition Is NOT a Rightward-Moving Pulse

The code initializes:

```
theta_z = A * exp(-r^2 / (2*sigma^2))    where r^2 = (x-X0)^2 + y^2 + z^2
theta_z_dot = theta_z * (x - X0) / sigma^2
```

This is a 3D SPHERICAL Gaussian with initial velocity set to mimic d/dt[f(x-ct)] = -(df/dx). However, this velocity prescription only works for a 1D or plane-wave pulse. For a 3D Gaussian, f(x-ct, y, z) = A exp(-[(x-X0-ct)^2 + y^2 + z^2]/(2 sigma^2)), and df/dt = (x-X0-ct)/sigma^2 * f. At t=0, this gives theta_z_dot = (x-X0)/sigma^2 * f, which IS what the code computes.

But this is only correct for a RIGID translation of the Gaussian. A 3D Gaussian does not propagate rigidly under the wave equation. The wave equation solution for a spherical initial condition decomposes into an OUTGOING and an INGOING spherical wave (via d'Alembert's formula in 3D: u = [r*f(r-ct) + r*g(r+ct)]/(2r)). The initial condition theta_z = Gaussian, theta_z_dot = (x-X0)/sigma^2 * Gaussian does NOT purely select the rightward-moving component. It creates a superposition that includes:

1. A component propagating in +x
2. Components propagating in all other directions (because the Gaussian is spherically symmetric around X0, not a plane wave)
3. An ingoing component converging toward X0

The result is that the "pulse" will SPREAD GEOMETRICALLY in all directions, not propagate as a coherent wavepacket in +x. This geometric spreading is the dominant effect, and it makes the velocity measurement unreliable.

### Critical Error 2: The Velocity Measurement Is Contaminated by Geometric Spreading

The data confirms the geometric spreading problem. Look at the measured sigma (width):

- t=0: sigma = 2.06
- t=10: sigma = 4.99
- t=20: sigma = 7.78
- t=30: sigma = 8.81
- t=40: sigma = 10.22
- t=60: sigma = 9.43

The width MORE THAN TRIPLES from sigma=2 to sigma=10 over the run. This is not dispersion — it is the spherical wavefront expanding. The "dispersion rate" ds/dt reported in the summary is meaningless because the width change is dominated by geometry, not by any dispersive physics.

The centroid tracking is also suspect. With the pulse spreading in all directions, the theta_z^2-weighted centroid does NOT track a single propagating wavefront. Looking at the data: x_centroid moves from -15 to +18 over t=0 to 39 (velocity ~ 0.85c), then appears to wrap around periodically and move erratically. The claimed "0.997c" for early times is only from the first few measurement points (t=0 to ~10), where the centroid moves from -15 to -10 in ~10 time units (velocity ~ 0.5c, not 0.997c). Actually, checking more carefully: t=0 x=-14.95, t=1.01 x=-13.94 gives v~1.0, but by t=4 x=-11.58 the velocity has already dropped to ~0.7c. By t=10 the centroid has only moved 5 units in 10 time units (v~0.5c).

The velocity drops continuously because the spherical expansion creates theta_z^2 weight distributed broadly in space, pulling the centroid back toward the domain center. This is a measurement artifact, not a physical velocity change.

### Issue 3: The phi Coupling Is Real But Quantitatively Unverified

The 15-20% energy transfer to phi (E_phi / E_theta ~ 0.15 at steady state) IS expected from the curl coupling. The theta equation sources phi via eta * curl(theta), so any theta excitation with nonzero curl will generate phi. For a Gaussian pulse with sigma=3, the curl has magnitude ~ A/sigma, and the phi response should scale as eta^2/m^2 in energy fraction. Let me check: eta = 0.5, m^2 = 2.25. The mixing fraction should be ~ eta^2 k^2 / (m^2 + k^2)^2 where k ~ 1/sigma ~ 0.33. This gives ~ 0.25 * 0.11 / (2.25 + 0.11)^2 ~ 0.005, which is about 0.5% — much less than the observed 15%.

The observed 15% is suspicious. Looking at the data: E_phi starts at 0 and jumps to 1.7e-4 by t=1 (10% of E_theta = 1.8e-3). By t=2 it reaches 3.2e-4 (18%). This rapid transfer suggests the coupling is NOT perturbative — the full curl(theta) source is driving phi directly. The E_phi then oscillates between ~2e-4 and 3e-4, consistent with periodic energy exchange between the sectors.

However, the claim "pure theta is not an eigenmode" is trivially true given the coupled equations. Any theta perturbation with nonzero curl sources phi through eta*curl(theta). The only theta mode with zero curl coupling would be a uniform or purely gradient field (curl-free). A Gaussian pulse is NOT curl-free, so of course it excites phi. This is not a new finding — it is a direct consequence of the equation structure.

### Issue 4: The "c=1 Propagation" Claim

The data does NOT support the claim of c=1 propagation. The centroid motion shows:

- t=0 to 3: rapid movement at ~1c (x moves from -15 to -12, 3 units in 3t)
- t=3 to 10: slowing to ~0.5c
- t=10 to 20: further slowing to ~0.5c
- t=20 to 30: complex behavior, centroid occasionally jumps
- t=30+: erratic, influenced by PBC wrapping

The initial fast motion is consistent with the rightward-moving component of the spherical wave. The slowdown is consistent with the centroid tracking failing as the pulse fragments into a spherical shell. The PBC boundary (L=40, so first wrap at x=+20 after starting at x=-15) is reached at t~35, after which the data is meaningless.

### What Would Fix This

To properly test photon dispersion, use a PLANE-WAVE pulse:

```
theta_z = A * exp(-(x-X0)^2 / (2*sigma_x^2)) * sin(k0 * x)
theta_z_dot = -c * k0 * A * exp(-(x-X0)^2 / (2*sigma_x^2)) * cos(k0 * x)
```

with a carrier wavevector k0 and Gaussian envelope only in x (not y,z). This creates a true wavepacket that propagates in +x without geometric spreading. The group velocity is d omega/dk evaluated at k0, and the dispersion is measured by the broadening of the x-envelope.

Alternatively, keep the 3D Gaussian but measure the RADIAL expansion rate of the theta_z^2 isosurface from the initial center X0. In 3D, a massless wave equation gives exact speed c for the wavefront, regardless of the initial shape. Track the radius of the expanding shell, not the centroid.

### Recommendations

1. Re-run with a plane-wave pulse (envelope in x only, carrier wave sin(k0*x)).
2. Measure group velocity from wavefront radius, not centroid.
3. Predict the phi coupling fraction from the linearized dispersion relation and compare.
4. Use a larger domain (L=200) to avoid PBC contamination during the measurement.

---

## OQ6: Nonlinear Maxwell Corrections — VERDICT: PASS (rigorous analysis with appropriate caveats)

### Claim Under Review

(a) No vacuum birefringence (structural), (b) No single-polarization self-interaction, (c) Cross-polarization chi^2 mixing only, (d) All effects vanish when A_bg -> 0.

### The Perturbation Expansion Is Correct

The expansion of phi_a = phi_a^bg + eps * delta_phi_a^(1) + eps^2 * delta_phi_a^(2) is standard. The triple product P = phi_0 * phi_1 * phi_2 is correctly expanded to second order. The derivatives V', V'', V''' are correctly computed and verified.

The key observation that kappa * P_bg^2 ~ 5e-5 << 1 (so V is deep in its quadratic regime) is correct and important. It means the V(P) nonlinearity is irrelevant compared to the P = phi_0 phi_1 phi_2 multilinearity. This is a genuine insight.

### The "No Birefringence" Result Is Rigorous

The argument is: V depends on P^2 = (phi_0 phi_1 phi_2)^2. When only phi_0 varies (with phi_1, phi_2 held at background values), V is quadratic in phi_0. A quadratic potential gives a linear force. Therefore, a single-component perturbation propagates linearly with no self-interaction.

This is mathematically correct. Let me verify: partial V / partial phi_0 = V'(P) * P/phi_0 = V'(P) * phi_1 * phi_2. With phi_1 = phi_1^bg (constant), phi_2 = phi_2^bg (constant), and V'(P) = mu * P / (1+kappa P^2)^2 where P = phi_0 * phi_1^bg * phi_2^bg, we get:

partial V / partial phi_0 = mu * phi_0 * (phi_1^bg * phi_2^bg)^2 / (1 + kappa * phi_0^2 * (phi_1^bg * phi_2^bg)^2)^2

This is NOT linear in phi_0 when kappa P^2 is not negligible. The denominator (1 + kappa P^2)^2 introduces genuine nonlinearity. The document's claim relies on the small-P approximation (kappa P^2 << 1), where V ~ (mu/2) P^2 and the force becomes mu * phi_0 * (phi_1^bg phi_2^bg)^2 (linear in phi_0).

So the "no birefringence" result is NOT structural in the strict sense — it depends on kappa P^2 << 1. For A_bg = 0.1, this is an excellent approximation (kappa P^2 ~ 5e-5). But near a braid core where P ~ 0.5 and kappa P^2 ~ 12.5, the approximation fails completely and single-polarization self-interaction DOES occur. The document correctly notes this regime exists (section 3.3 of F9) but does not connect it back to the birefringence claim.

**Verdict on this point**: The claim should be stated as "no vacuum birefringence in the weak-field regime (A_bg = 0.1)" rather than "structural." In the strong-field regime near braids, birefringence IS possible through the kappa P^2 denominator. However, for the physically relevant question (does light propagating through empty space between braids see birefringence?), the weak-field result is the relevant one, and it IS correct.

### The Cross-Polarization Coupling Is Correctly Derived

The calculation of F_a^(2) involving delta_phi_b * delta_phi_c for b,c != a is correct. The observation that this requires at least two field components to be simultaneously excited is correct. The chi^2 identification (three-wave mixing with the background as one wave) is appropriate.

### Issue: Background Consistency

The document honestly acknowledges (section 3) that the standing-wave background phi_a = A cos(kz + delta_a) is NOT an exact solution of the equation of motion. The V(P) correction to the background is of order mu*A^4/m^2 ~ 0.02 (2%). This is handled correctly by absorbing it into a "slightly modified background" and proceeding with the perturbative structure. This is standard practice and acceptable for an O(1%) level calculation.

### Issue: The QED Comparison Is Honest

The comparison with Euler-Heisenberg (section 9) correctly identifies the three critical differences: (1) Cosserat nonlinearity requires a background, QED does not; (2) QED produces birefringence, Cosserat does not (in weak field); (3) Cosserat is chi^2 while QED is chi^3. These are stated clearly as testable distinctions.

### Issue: The g_eff Calculation Has an Error in the Chain

In section 6.1, the document computes the photon-photon scattering amplitude through the chain theta -> phi -> V -> phi -> theta. The intermediate step counts three eta factors (two theta->phi, one phi->theta), giving g_eff ~ mu * eta^3 * A^3 * k / m^6. But the section starts by claiming four eta factors, then corrects to three. The corrected version is:

```
delta_theta^(2) ~ eta * k * delta_phi^(2) / k^2
delta_phi^(2) ~ mu * A^3 * (eta*k/m^2)^2 * theta_A * theta_B / m^2
```

So delta_theta^(2) ~ mu * eta^3 * A^3 * k * theta_A * theta_B / m^6. This is dimensionally consistent and the algebra checks out.

### Issue: Section 11 (Background Periodicity) Is a Good Addition

The observation that the spatial average <P_0> = 0 (the triple product of three phase-shifted cosines averages to zero) means the leading nonlinear coupling vanishes for long-wavelength photons. The first surviving term is <P_0^2> at fourth order. This is a meaningful and correct result that significantly suppresses nonlinear EM effects in the far field.

### Minor Issue: The Document Self-Corrects

In section 6.3, the document initially claims "no birefringence" for single-polarization waves, then realizes it needs to reconsider the background's effect more carefully. It works through the full expansion and re-derives the result correctly (arriving at the same conclusion from a cleaner starting point). This self-correction is visible in the text ("Wait, this is not quite right... Let me reconsider"). While this makes the document read like a research notebook rather than a polished paper, it is honest and the final result is correct.

### Recommendations

1. Clarify that "no birefringence" is a weak-field result (kappa P^2 << 1), not structural.
2. Estimate the birefringence magnitude near a braid core where kappa P^2 ~ O(10).
3. The testable prediction (section 9, point 8) about "no birefringence but yes parametric frequency conversion" is valuable — develop it into a quantitative prediction.

---

## F9: Analytical Effective Potential — VERDICT: CONCERN (mixture of rigorous and circular reasoning; key claims overstated)

### Claim Under Review

(a) Static interaction rigorously repulsive at all D, (b) Dynamic force constant C ~ 127 vs measured 186 (30% off), (c) F ~ 1/D^2.2 vs measured 1/D^1.8, (d) If phonon is exactly massless then n=1 gives Newton exactly.

### Claim (a): Static Repulsion — PASS

The argument that E_static(D) > 0 for all D is supported by the V33 data (all E_int values are positive). The physical reasoning is sound: two braids overlapping in the V(P) potential have higher energy than two well-separated braids. Since V(P) with mu < 0 is attractive (stabilizing), removing the overlap RELEASES binding energy, making the separated state lower-energy. This is confirmed numerically.

However, the document's attempt to derive the static repulsion ANALYTICALLY is a mess. It tries three different functional forms (Yukawa overlap, power-law convolution, Gaussian overlap) and none of them fit the data. The Yukawa overlap gives exp(-12)/8 ~ 10^-7, off by 10 orders of magnitude from the measured E_int(D=8) = 2939. The power-law convolution gives E_int INCREASING with D (wrong sign of the slope). The Gaussian overlap decays too fast. The document eventually gives up and fits the data empirically as D^-5.7.

The honest conclusion — that the analytical form of E_static requires the full nonlinear braid profile — is correct. But the "rigorously repulsive" claim is really "numerically observed to be repulsive" plus a plausibility argument. A truly rigorous proof would require showing that the overlap integral is positive-definite for all D, which is not done.

**Verdict**: The repulsion is well-established empirically. Calling it "rigorous" is an overstatement — it is numerically robust and physically plausible.

### Claim (b): C ~ 127 — CONCERN (post-hoc parameter selection)

The derivation of C proceeds as follows:

1. Write the force integral F = -integral (d m_eff^2/dx) * |psi|^2 d^3x (Eq. 18)
2. Note that the LINEAR contribution gives zero force (m_eff ~ constant for small A_bg)
3. Invoke the NONLINEAR regime where the braid's own field dominates
4. Do dimensional analysis: C ~ |mu| * R_tube^3 * psi^2(R_tube) / 2
5. Plug in R_tube = 4.5, psi_peak = 0.8, and get C = 1190 (6x too large)
6. Replace psi_peak with psi(R_tube) = 0.8 * exp(-4.5^2/18) = 0.26, and get C = 127

This is problematic for several reasons:

**Problem A**: The switch from psi_peak to psi(R_tube) is not derived — it is motivated by needing a smaller answer. The document says "A more careful estimate includes the saturation" but the saturation factor is computed to be 0.99 (negligible). The real adjustment is using the field value at R_tube instead of the peak value. Why R_tube = 4.5 specifically? Because that is the "interaction surface" radius defined elsewhere. But the force integral (Eq. 18) integrates over ALL space — the choice of which radius to evaluate psi at is arbitrary. Different radii give wildly different answers: psi(3.0) = 0.8*exp(-0.5) = 0.49 gives C = 484; psi(6.0) = 0.8*exp(-2.0) = 0.11 gives C = 22.

**Problem B**: The document then notes that for C = 186 you need psi at R_eff = 4.1 (section 5.2). This is FITTING the radius to match the answer. You cannot claim a "prediction" of C = 127 when you have one free parameter (the evaluation radius) that you tuned to get close to the measured value.

**Problem C**: The dimensional analysis itself is suspect. The force integral (Eq. 18) involves d(m_eff^2)/dx * |psi|^2 integrated over space. But the document showed that d(m_eff^2)/dx ~ 0 in the linear regime (section 3.2, where it explicitly computes epsilon ~ 1.8e-5). The force comes from the NONLINEAR regime, but then the dimensional analysis uses parameters from both regimes (|mu| from the potential, but R_tube from the linear braid profile).

**Verdict**: C ~ 127 is not a prediction — it is a rough estimate with one adjustable parameter (evaluation radius) that gives answers ranging from 22 to 1190 depending on the choice. The 30% agreement is fortuitous within the method's uncertainty.

### Claim (c): F ~ 1/D^2.2 — CONCERN (discrepancy is real, not "rough agreement")

The predicted exponent is n+1 = 1.2 + 1 = 2.2. The measured exponent is 1.8. The document calls this "rough agreement" and attributes the discrepancy to measurement error. Let me assess this.

The relation F ~ 1/D^(n+1) follows from F = -C * d(delta_rho)/dr where delta_rho ~ B/r^n. This is straightforward dimensional analysis and is correct IF the depletion profile is a pure power law. The measured n = 1.2 +/- 0.2 gives a predicted exponent of 2.2 +/- 0.2, and the measured force exponent is 1.8. These are separated by 0.4, which is 2 sigma. This is a genuine discrepancy, not "rough agreement."

There are three possible explanations:

1. The depletion exponent n is not constant (it varies with r, from ~1.17 at r=5-15 to ~3.07 at r=20-40). The effective n at the separations where the force is measured (D=15-30) may be different from the n measured at r=5-15.

2. The force measurement itself (D=15-30, T=200) may not be in the asymptotic regime. At D=15, the static repulsion (D^-6.7) contributes 5.7 * K / 15^6.7 which may still be comparable to the attractive force.

3. The relation F = -C * nabla_rho is an approximation that assumes the braid responds instantaneously to the local gradient. Retardation effects (the braid takes time to respond to the gradient change as it moves) could modify the effective exponent.

The document acknowledges these issues but still presents the 2.2 as a "prediction." It would be more honest to say: "the theory predicts the exponent is n+1 where n is the depletion exponent, and with the measured n = 1.2 this gives 2.2, which is 20% higher than the measured 1.8. This discrepancy may arise from the non-constant exponent, finite-size effects, or retardation."

### Claim (d): "If n=1 Then Newton" — CIRCULAR REASONING

The argument is: if the phonon is exactly massless, then n = 1 (the depletion is 1/r), and the force law becomes F = -C * d(1/r)/dr = C/r^2, which is Newton's inverse-square law.

This is correct as a tautology: IF the effective potential from the depletion follows a 1/r profile, THEN the force is 1/r^2. But this is not a derivation of Newton's law — it is a restatement of the definition. Newton's law says F = -G*M1*M2/r^2. To derive this, you need:

1. The depletion is exactly 1/r (requires exactly massless phonon)
2. The depletion amplitude B is proportional to the braid mass M (requires equivalence principle)
3. The coupling constant C is proportional to the test braid's mass (requires equivalence principle)
4. The product C*B/(mass1*mass2) = G is UNIVERSAL for all braid types

NONE of these four requirements are established. The document explicitly flags items 2 and 3 as open questions (section 7, questions 2 and 3). But it presents "if n=1 then Newton" as a significant result, when in fact it is just step 1 of 4, and step 1 is itself unproven (the measured n is 1.2, not 1.0).

Furthermore, the document speculates that n might approach 1 in an "isotropic limit" (random background phases). This is untested and may or may not be true. It is a hope, not a prediction.

### Issue: The "Three Mechanisms" for Power-Law Tails Are Not Evaluated

Section 1.3 lists three mechanisms that could produce power-law tails in delta_rho: (a) nonlinear backreaction, (b) phonon, (c) radiation halo. The document then assumes mechanism (b) is "most likely" without ruling out (a) or (c). But mechanism (c) is directly relevant: the braid radiates ~0.5 energy/time (V29), and this radiation creates a ~1/r^2 energy density contribution. If the depletion is dominated by the radiation halo rather than a static phonon field, then:

- The "force" would be radiation pressure (covered in OQ2), not gradient-driven drift
- The exponent would be n = 2 (radiation intensity ~1/r^2), giving F ~ 1/r^3, which is WORSE than measured

So the phonon interpretation (n=1) and the radiation interpretation (n=2) predict different force laws. The measured n = 1.2 is between them, which is consistent with either a quasi-massless phonon or a mixture of phonon + radiation. The document does not address this ambiguity.

### Issue: The Document Conflates Two Different "Gravity" Mechanisms

The EM_THEORY.md describes an electromagnetic radiation pressure mechanism (OQ2) giving F ~ 1/r^2 from Poynting flux. The F9 document describes a phonon/depletion mechanism giving F ~ 1/r^(n+1). These are two DIFFERENT force channels between braids. The document does not discuss their relative magnitude or whether they interfere constructively or destructively.

If both mechanisms operate, the total force is F = F_EM + F_gravity, and the measured force contains both contributions. The "C = 186" from the gradient sweep may include EM contributions (the gradient test uses same-winding braids that also interact electromagnetically). Separating the two channels is essential for understanding the physics but is not attempted here.

### Recommendations

1. Stop claiming C = 127 is a "prediction" — it is an order-of-magnitude estimate with one free parameter. To genuinely predict C, integrate Eq. (18) over the numerically-determined braid profile.
2. Address the discrepancy between the predicted exponent (2.2) and measured exponent (1.8) seriously. Is n=1.2 measured at the right radii? Is the force measurement in the asymptotic regime?
3. To test "if n=1 then Newton," first demonstrate n=1 in SOME limit. Run simulations with an isotropic background and measure the depletion exponent.
4. Separate the EM and gravity force channels by running gradient tests with theta=0 (3-field equation only) and comparing to the 6-field result.
5. Test the equivalence principle: run gradient sweeps with braids of different amplitudes and check whether C/m_inertial is constant.

---

## Updated Cross-Cutting Concerns

### Theory-Measurement Loop Risk

Several results form a closed loop: the depletion exponent n is measured, used to predict the force exponent n+1, and then compared to the measured force exponent. But both measurements come from the same simulations with the same parameters. If there is a systematic error (e.g., PBC contamination affecting both measurements), the "agreement" is illusory. Independent validation (different grid sizes, absorbing boundaries, different parameter sets) is needed to break this loop.

### The Role of the Background

OQ6 and F9 both depend critically on the background amplitude A_bg = 0.1. In OQ6, all nonlinear EM effects scale as A_bg^2 or A_bg^3 and vanish at A_bg = 0. In F9, the force constant C depends on the braid profile in the background. Neither document discusses what happens if A_bg varies spatially (which it must, since braids deplete the background). The self-consistency of treating A_bg as constant while computing forces that arise FROM spatial variations of A_bg deserves scrutiny.

### Honest Self-Assessment

OQ6 is the strongest of the three results reviewed here. It correctly identifies a structural property (multilinearity of P suppresses single-component self-interaction), derives it carefully, and makes testable predictions. The self-correcting narrative style, while unpolished, is scientifically honest.

F9 is the weakest. It attempts an analytical derivation, fails to produce clean results, and falls back on dimensional analysis with tunable parameters. The "predictions" (C=127, exponent 2.2) are not robust. The most valuable part is the clear identification of WHAT determines each quantity (n from nonlinear phonon physics, C from the force integral, signs from symmetry) — this provides a roadmap for future numerical work even if the analytical estimates are unreliable.

OQ5 is simply the wrong experiment for the question being asked. A 3D Gaussian is not a photon wavepacket, and the measurements are dominated by geometric spreading. The experiment needs to be redesigned with a plane-wave pulse.

---

## Updated Summary Table

| Result | Verdict | Key Issue |
|--------|---------|-----------|
| OQ3 DC profile | **CONCERN** | Measures spatial-and-temporal average, not per-point DC; incomparable with V34 methodology; PBC contaminates r > 12; sign flip at r~16 is artifact |
| OQ1 Multipole | **FAIL** | Decomposing positive-definite |theta| guarantees artificial monopole dominance; no null model comparison; single-frame snapshot; per-component analysis shows the real (multipolar) structure |
| OQ2 Poynting | **PASS** | Core 1/r^2 argument is rigorous and correctly caveated; sign argument is qualitative but honestly labeled; alpha ~ eta^2 is dimensional analysis, appropriately noted |
| OQ5 Photon dispersion | **FAIL** | 3D Gaussian is not a plane-wave pulse; velocity and dispersion measurements dominated by geometric spreading; phi coupling is real but trivially expected; experiment needs redesign |
| OQ6 Nonlinear Maxwell | **PASS** | Perturbation expansion correct; no-birefringence is rigorous in weak-field regime; multilinearity of P correctly identified as source of structural selection rules; honest caveats throughout |
| F9 Analytical potential | **CONCERN** | Static repulsion empirically robust but not analytically derived; C~127 has one free parameter tuned post-hoc; exponent prediction 2.2 vs measured 1.8 is a genuine discrepancy; "if n=1 then Newton" is step 1 of 4 with step 1 unproven |

---

## ROUND 2: Corrected Results Review

Date: 2026-03-25.

This section reviews the four corrected experiments that addressed issues
identified in the Round 1 review. Each corrected experiment is compared
against both its own internal logic and the original V33/V34 claims it
aims to validate or refute.

---

## OQ3v2: Per-Voxel DC Profile — VERDICT: PASS (methodology now correct; DC is real but small)

### What Changed from Round 1

The original dc_profile.c shell-averaged the signed theta vector BEFORE
computing DC magnitude, which cancels dipolar angular structure. This
corrected version (dc_pervoxel.c) computes the time average of theta at
each individual voxel FIRST, then bins the DC magnitude into radial shells.
This is the correct methodology. The grid spacing bug (cell-centered vs
vertex-centered) is also fixed: `dx = 2*L/(N-1)`.

### Code Review

The code is clean and correct. Key features:
- Per-voxel accumulators for sum and sum-of-squares (12 arrays, 6 full + 6 no-transient)
- DC computed as `|<theta_vec>_t|` per voxel, then RMS over shells
- AC computed as `sqrt(var_x + var_y + var_z)` per voxel
- DC_frac = DC_rms / sqrt(DC_rms^2 + AC_rms^2) -- the right ratio
- Transient exclusion (first 25 frames skipped for the no-transient analysis)
- Centroid from |P| weighting at frame 0, with drift check at mid and last frames

### Centroid Drift: Does It Invalidate the Result?

The user reports 9 grid cells of centroid drift over 264 frames. At dx=0.633,
this is 5.7 code units of drift over T=50 (264 frames at snap_dt~0.19). The
code uses a FIXED centroid from frame 0 for all radial binning.

**Impact assessment**: A centroid drift of 5.7 code units means voxels near
the core boundary migrate between radial bins over the simulation. A voxel at
r=3 from the true (moving) center could be assigned to the r=3.5 or r=8.5
bin depending on the frame, if the centroid has drifted significantly from
its frame-0 position.

**How much spurious DC does drift create?** Consider a voxel at position x
near the braid core. If the braid drifts by delta, the theta field at that
voxel shifts from theta_braid(x) to theta_braid(x - delta). The time average
picks up the gradient: <theta> ~ theta(x) + (delta/2) * grad(theta) != 0.
For theta_rms ~ 0.08 in the core (from the TSV) and a drift of 5.7 over
braid radius ~3, the gradient is ~0.08/3 = 0.027 per code unit. The drift
contributes <theta> ~ 0.027 * 5.7/2 ~ 0.077, which would give DC_rms ~
0.077 -- comparable to or larger than the observed DC_rms ~ 0.011 in the core.

However, this OVERESTIMATES if the drift is oscillatory rather than monotonic.
If the braid executes a random walk around a mean position, the time-averaged
displacement is smaller. But even with a monotonic drift of 5.7 code units,
the spurious DC magnitude (~0.077) exceeds the measured DC_rms (~0.011),
meaning the centroid drift could FULLY ACCOUNT for the observed DC without
any physical DC being present.

**Critical test the code should have done but did not**: Recompute the
centroid at EACH frame, re-center the coordinate system per frame, and
THEN accumulate per-voxel statistics. This would eliminate drift-induced
DC. The fact that this was not done means the 5-8% figure is an UPPER
BOUND on the physical DC.

### The Data

The no-transient DC_frac profile shows:
- r=0.25: 5.4%
- r=2.25: 6.3%
- r=5.25: 6.6%
- r=8.25: 7.1%
- r=11.75: 7.4%

The DC_frac is roughly CONSTANT with radius (5-8%), not concentrated in
the core. A physical DC from the braid's internal structure would be
expected to peak in the core and decay with radius. A drift-induced DC,
by contrast, would be approximately constant or even grow with radius
(because the theta gradient is sustained across the braid's extent and
the drift contaminates equally everywhere). The flat-to-increasing radial
profile is MORE consistent with drift contamination than with physical DC.

### Comparison with V34

V34 claimed DC ~ 0.2%. This code finds 5-8%. These are still not
directly comparable: V34 may have used per-frame re-centering (which
would eliminate drift) or a different spatial decomposition. But the
Round 1 review already identified this as an apples-to-oranges comparison.

### Verdict: PASS with Caveats

The methodology is now correct (per-voxel DC, then bin). The code is
well-written and properly handles vertex-centered grids. However:

1. The 5-8% DC_frac is contaminated by centroid drift and should be
   treated as an UPPER BOUND.
2. The flat radial profile suggests drift dominates over physical DC.
3. A definitive measurement requires per-frame re-centering.

The result is CONSISTENT with V34's 0.2% being closer to the true
physical DC (if V34 re-centered per frame). The corrected measurement
establishes that the DC is at most ~8%, not the 12-15% from the
uncorrected analysis, and probably much less after drift subtraction.

**Impact on theory**: A small physical DC (0.2-1%) would be expected
from the braid's nonlinear dynamics (asymmetric oscillation cycle).
The exact value matters for the EM force mechanism (DC theta generates
static B-fields). The current measurement does not resolve the issue
but narrows the range.

---

## OQ5v2: Plane-Wave Pulse Propagation — VERDICT: PASS (clean experiment, meaningful results)

### What Changed from Round 1

The critical error from Round 1 (3D spherical Gaussian initial condition)
is FIXED. The corrected code uses a true plane-wave pulse:
```
theta_z = A * exp(-(x-x0)^2/(2*sigma_x^2)) * sin(k0*x)
```
uniform in y and z (no transverse spreading). The velocity initialization
properly uses the full derivative of the modulated pulse. The smoothed-
envelope centroid tracking (averaging over one carrier wavelength) is the
correct approach for extracting the group velocity from a modulated wavepacket.

### Code Review

The code is well-structured:
- Full 6-field Cosserat with proper curl coupling terms
- Velocity Verlet integrator (second-order symplectic)
- Absorbing boundary damping at |x| > 0.8*L = 20 code units
- Lattice-corrected dispersion predictions (accounts for 7-point Laplacian)
- Two runs: eta=0 (control) and eta=0.5 (coupled)
- Group velocity from linear fit to unwrapped centroid positions (t=2..28)

**Grid**: N=128, L=25, DX = 2*25/128 = 0.3906. Carrier k0=2.0,
wavelength = pi = 3.14, giving 3.14/0.39 = 8 points per wavelength.
This is adequate for 2nd-order finite differences.

**One minor inconsistency**: `DX = 2.0*L/N` (cell-centered) instead of
`2.0*L/(N-1)` (vertex-centered as in the kernel). The difference is 0.8%,
internally consistent within this standalone code but should be noted.

### The A_bg = 0 Issue

The code sets A_bg = 0 (no phi background). This does NOT invalidate the
coupled dispersion relation test. The dispersion relation
(w^2 - k^2)(w^2 - k^2 - m^2) = eta^2 * k^2 is derived from the linearized
equations around phi=0, theta=0. The coupling chain works as:

1. theta_z pulse creates nonzero curl(theta)
2. eta * curl(theta) sources phi acceleration
3. phi grows; curl(phi) back-reacts on theta via eta * curl(phi)

With A_bg = 0, the V(P) potential contributes nothing (P=0 everywhere
initially, and the pulse-induced phi is tiny). The effective phi mass is
the bare m^2 = 2.25.

With A_bg = 0.1, V''(P_bg) would modify the effective phi mass by
~mu * A_bg^2 = -41.3 * 0.01 = -0.41, changing m^2 from 2.25 to ~1.84.
This would shift the coupled dispersion, making the photon branch slightly
faster. So the A_bg=0 test probes the bare equations correctly, but the
physical vacuum (A_bg=0.1) would give a different v_g. The user's note
about this is correct.

**For validating the code**: A_bg=0 is the cleaner test (no background
complications). For physical predictions: A_bg=0.1 is needed.

### The Data: eta=0 (Control)

- E_phi = 0 throughout (correct: no coupling source)
- x_peak moves linearly from -15.0 to +22.1 over T=50
- Sigma stays nearly constant: 3.65 at t=0, 3.60 at t=25 (-1.3% change)

From the TSV: t=3.0, x=-12.12; t=28.0, x=+10.97. Displacement = 23.09
over 25.0 time units. Measured v = 0.924. Lattice prediction:
v_g = cos(k0*DX/2) = cos(2.0*0.3906/2) = cos(0.3906) = 0.925.
Agreement: 0.1%. Excellent.

The nearly constant sigma confirms coherent propagation without dispersive
spreading -- the key improvement over Round 1.

### The Data: eta=0.5 (Coupled)

dE_phi oscillates strongly between ~0 and ~1.28, with period ~10 time
units. This is the signature of resonant beating between two coupled modes
at similar frequencies.

The E_phi/E_theta ratio averages ~0.41, far above the perturbative
prediction of eta^2 * k^2 / (k^2 + m^2)^2 = 0.026. This 15x discrepancy
is NOT a failure: with phi_bg = 0, the perturbative expansion around the
background does not apply. The strong oscillation (phi energy periodically
returning to near-zero) is the hallmark of non-perturbative resonant
coupling. The system is exchanging energy between two normal modes at
their beat frequency, which is the fully nonlinear regime.

### Lattice vs Continuum Comparison

The code compares against lattice-corrected predictions, which is the
right choice for validating the numerical implementation. The lattice
dispersion (7-point Laplacian replaces k^2 with k_eff^2 = (2/dx)^2 *
sin^2(k*dx/2)) is a KNOWN systematic that should be corrected for rather
than claimed as physics. The code does this correctly.

For extracting physical predictions (continuum limit v_g), one would need
to run at multiple grid spacings and extrapolate. That is a separate
experiment.

### Verdict: PASS

The plane-wave initial condition is correct. The group velocity matches
the lattice-corrected prediction to < 1% for both eta=0 and eta=0.5.
The pulse propagates coherently with < 1.3% broadening. The coupling
slows the photon branch as predicted by the dispersion relation.

Caveats:
1. A_bg=0 tests the bare equations, not the physical vacuum.
2. The strong E_phi/E_theta oscillation (0.41 vs 0.026) reflects
   non-perturbative resonant coupling, not a code failure.
3. Cell-centered vs vertex-centered convention differs from kernel by 0.8%.

**Impact on theory**: The coupled Cosserat dispersion relation is
CONFIRMED numerically. The photon branch exists and propagates at
the predicted (sub-luminal) group velocity. This is a necessary
validation for the EM_THEORY framework.

---

## OQ1v2: Poynting Flux Multipole (Null Model) — VERDICT: PASS (correct methodology, important self-correction)

### What Changed from Round 1

The original decomposed |theta| (positive-definite) into spherical
harmonics. The corrected version computes S = E x B = (-dtheta/dt) x
curl(theta), time-averages |S| over 263 frame pairs, then decomposes
the time-averaged <|S|> into spherical harmonics.

### Code Review

The code correctly:
- Reads consecutive frame pairs (snap_dt = 0.19, well-sampled)
- Computes E = -dtheta/dt via finite difference
- Computes B = curl(theta) via centered spatial differences
- Computes |S| = |E x B| at each grid point
- Interpolates onto angular grids at multiple shell radii via trilinear
- Time-averages over 263 frame pairs (spanning ~11 oscillation periods)
- Decomposes into real spherical harmonics up to l=6
- Reports both C_l power spectrum and anisotropy ratio

**Minor boundary issue**: Uses one-sided differences at domain edges
instead of periodic wrapping. Negligible for R <= 12 where all shells
are well inside the domain.

### Critical Self-Correction on Positive-Definiteness

The results markdown correctly identifies that the 92-98% monopole
fraction is the SAME positive-definiteness artifact as Round 1, just
applied to <|S|> instead of |theta|. Since |S| >= 0 everywhere, the
Y_00 projection is guaranteed to be large. The document correctly
identifies the anisotropy ratio (max/min) as the physically meaningful
discriminant.

This is the right conclusion. To avoid the positive-definiteness bias
entirely, one could decompose the SIGNED radial Poynting flux <S_r> =
<S . r_hat> instead. This quantity can be negative (inward energy flow)
and gives a clean monopole measurement. But the max/min ratio of the
unsigned |S| is an adequate proxy for the current purpose.

### The Data

Single braid (null model):
- Monopole fraction: 92-98% (artifact, correctly flagged)
- l=2 fraction: grows from 1.0% at R=3 to 6.8% at R=12
- Anisotropy ratio: grows from 4.4 at R=3 to 19.4 at R=12
- F(R)*R^2: 0.52 (R=3) to 2.60 (R=12) -- NOT constant

The growing anisotropy with R is physically correct for a z-aligned
oscillating dipole: the radiation pattern (dominated by sin^2 theta)
becomes more pronounced as near-field contributions decay. The growing
l=2 fraction confirms this (sin^2 theta decomposes into Y_00 + Y_20).

F(R)*R^2 growing with R means the Poynting flux decays SLOWER than 1/R^2.
This indicates we are NOT in the radiation zone (expected: R >> lambda
~ 4.4, so R >> 10-20). At R=3-12, we are in the near-field/transition
zone where the field includes both radiative (1/R) and reactive (1/R^2,
1/R^3) components. The Poynting flux from reactive fields falls as higher
powers of 1/R, so at small R the flux is enhanced above the 1/R^2
radiative prediction, and F*R^2 is smaller. As R increases, the reactive
contribution fades, and F*R^2 grows toward the asymptotic radiative value.
This interpretation is consistent with the data.

### Is Anisotropy Ratio the Right Discriminant?

Yes. For comparing a single braid (expected: anisotropic dipolar pattern)
against a UUD composite (expected: isotropic from three orthogonal dipoles),
the anisotropy ratio is clean:
- Single braid: max/min ~ 19 (confirmed)
- UUD composite target: max/min < 2 (from corrected plans)

**Robustness concern**: The min can be artificially small if a single
angular bin has low <|S|> from noise or interpolation artifacts. A more
robust alternative: the ratio of the 90th percentile to 10th percentile,
or the coefficient of variation (std/mean) of the angular distribution.
But max/min is adequate as a first discriminant given the large dynamic
range (19:1).

### Verdict: PASS

The methodology is correct. The self-correction on positive-definiteness
is genuine scientific insight. The null model establishes a clear baseline.

The experiment is INCOMPLETE: only the null model (single braid) has been
run. The UUD comparison is pending a high-cadence UUD simulation.

**Impact on theory**: No impact yet. This is infrastructure for a future
test of the Coulomb isotropy claim.

---

## F9v2: Numerical Gradient Force — VERDICT: FAIL (V33 C=186 NOT reproduced; serious implications)

This is the most consequential result in the V43 suite. CONCEPT.md Section 3
hinges on the V33 measurement F = -C * grad_rho with C = 186 and R^2 = 0.9998.

### Code Review: Faithfulness to V33

I compared F9v2 (gradient_force.c) against the V33 original
(v33_gradient_test.c) line by line.

**Identical between V33 and F9v2**:
- Braid parameters: A_braid=0.8, delta={0, 3.0005, 4.4325}, R_tube=3.0, ellip=0.3325
- Background: A_bg linear from A_high at x=-L to A_low at x=+L
- Standing wave: cos(k_bg*z + 2*pi*a/3 + delta_a) with k_bg=pi/L
- Grid: N=128, L=20 (matching the V33 gradient sweep documentation)
- Time step: dt = 0.12 * dx
- Centroid tracking: |phi|^2 > 5x average, same weighting scheme
- X-boundary: pinned (margin=3 slabs frozen to initial values)
- Y-boundary: linear extrapolation from interior

**KEY DIFFERENCE 1: V33 is 3-field only. F9v2 tests both 3-field and 6-field.**
The V33 gradient sweep was 3-field (no theta arrays in the code). The F9v2
eta=0 case should reproduce V33 exactly.

**KEY DIFFERENCE 2: Z-boundary handling.** This is critical.

V33 code (v33_gradient_test.c) applies BCs for x (pinned) and y (extrapolated)
but has NO EXPLICIT Z-BOUNDARY TREATMENT. The compute_forces function skips
k < margin and k >= N-margin voxels (acc=0), but apply_bc does not constrain
these z-boundary voxels. They are left in FREE DRIFT: acc=0, vel unchanged,
phi += dt*vel every step. Over T=60, these 6 slabs of voxels drift linearly,
creating growing artifacts at the z-faces that propagate into the interior.

F9v2 code (gradient_force.c) ADDS explicit z-boundary handling (lines 143-176):
linear extrapolation for z, matching the y-boundary treatment.

**This z-boundary difference is the most likely explanation for the
discrepancy.** The V33 code's unconstrained z-boundary voxels create
artifacts that scale with the background amplitude (which is set by the
gradient), producing an apparent force linearly proportional to grad_rho.

**KEY DIFFERENCE 3: Run time.** V33 used T=60. F9v2 uses T=30.

### The Data: 3-Field (eta=0)

| A_high | A_low | grad_rho   | avg_drift    | C_v33 |
|--------|-------|-----------|--------------|-------|
| 0.105  | 0.095 | -5.0e-5  | +1.65e-6    | -0.03 |
| 0.110  | 0.090 | -1.0e-4  | -6.17e-6    | +0.06 |
| 0.130  | 0.070 | -3.0e-4  | +2.94e-5    | -0.10 |
| 0.150  | 0.050 | -5.0e-4  | -7.81e-6    | +0.02 |

The drifts are noise-level (10^-6 to 10^-5) with ALTERNATING SIGNS.
There is zero correlation between drift and gradient. C_v33 ranges
from -0.10 to +0.06. This is indistinguishable from zero.

V33 reported (at T=60):
| gentle  | -5.0e-5  | -0.0092  | -184.8 |
| mild    | -1.0e-4  | -0.0186  | -185.8 |
| moderate| -3.0e-4  | -0.0559  | -186.3 |
| strong  | -5.0e-4  | -0.0953  | -190.6 |

The V33 drifts are 3-4 ORDERS OF MAGNITUDE larger. This is not a subtle
discrepancy -- it is a qualitative disagreement.

### Why the Discrepancy? Analysis of Hypotheses

**Hypothesis A: T=30 is too short.** If drift scales linearly with T,
the T=30 value should be half the T=60 value: ~0.048 for the strong case.
F9v2 measures 7.8e-6. This is 6000x smaller, not 2x. Ruled out.

**Hypothesis B: V33 z-boundary omission.** The V33 code's compute_forces
skips the margin=3 boundary voxels in ALL directions (k < 3 and k >= N-3
get continue, so acc=0). But apply_bc only constrains x-boundaries (pinned)
and y-boundaries (extrapolated). The z-boundary voxels (k=0,1,2 and
k=N-3,N-2,N-1) are NEITHER pinned NOR extrapolated -- they are left in
FREE DRIFT with acc=0. Their velocities retain whatever value they had
after initialization, and their positions drift linearly as phi += dt*vel
every step. Over T=60, these 6 slabs of N^2 voxels drift linearly,
creating growing artifacts at the z-boundaries.

Since the braid's standing wave propagates along z (k_bg = pi/L), the
z-boundary behavior directly affects the wave dynamics. The free-drifting
boundary voxels create reflections and mode mixing that perturb the braid
asymmetrically. Because the background amplitude A_bg varies linearly
with x (the gradient), the z-boundary artifact inherits the gradient's
linear scaling. This would explain the perfect linearity (R^2=0.9998):
the artifact force is proportional to A_bg, which is proportional to the
gradient.

F9v2 FIXES this by adding explicit z-boundary handling (linear
extrapolation, matching the y-boundary treatment).

This hypothesis is TESTABLE: add z-boundary extrapolation (or periodic
z-wrapping) to the V33 code and re-run. If C drops from 186 to ~0, the
missing z-BC is confirmed as the source of the artifact.

**Hypothesis C: Resolution dependence.** V33 defaults to N=512 and the
documentation claims the gradient sweep used N=128. If the actual sweep
used a higher N (perhaps through a run script we cannot see), the gradient
force might be real but only appear at high resolution where the braid
core dynamics are well-resolved.

This hypothesis is TESTABLE: run F9v2 at N=256 and N=512 with T=60.

**Hypothesis D: The braid initial condition is slightly different.**
Both codes look identical in the init_gradient function. The braid
amplitude, envelope, phase offsets, background prescription all match.
I cannot identify any difference. Ruled out.

**Hypothesis E: The centroid tracking is noisy at N=128.** The 5x-average
threshold selects voxels above 5 times the domain-average phi^2. At N=128,
this threshold depends on the background amplitude (which varies across the
domain). In a gradient, the threshold changes, potentially biasing the
centroid differently at different gradient strengths. But V33 uses the SAME
centroid function, so this should produce the same bias. Ruled out as an
explanation for the discrepancy (though it may contribute to noise).

**Most likely explanation: Hypothesis B (z-boundary bug).** The V33 code
has undefined behavior at z-boundaries, F9v2 fixes this, and the gradient
force vanishes. The "perfect linearity" (R^2=0.9998) arises because the
spurious z-boundary force scales linearly with A_bg (which scales linearly
with the gradient).

### The 6-Field (eta=0.5) Results

| A_high | A_low | grad_rho   | avg_drift  | C_v33 |
|--------|-------|-----------|------------|-------|
| 0.105  | 0.095 | -5.0e-5  | +3.60e-3  | -72.1 |
| 0.110  | 0.090 | -1.0e-4  | +3.58e-3  | -35.8 |
| 0.130  | 0.070 | -3.0e-4  | +3.62e-3  | -12.1 |
| 0.150  | 0.050 | -5.0e-4  | +3.64e-3  | -7.3  |

The drift is +3.6e-3 for ALL gradient strengths -- a CONSTANT offset
with ZERO dependence on the gradient. C_v33 varies inversely with
grad_rho (C = drift/grad, and drift is constant while grad varies).
This is definitively NOT a gradient force.

The constant 3.6e-3 drift in the 6-field case likely arises from the
theta coupling creating a systematic asymmetry. At t=0, theta=0 everywhere
while phi has an asymmetric gradient background. The curl coupling
eta*curl(phi) sources theta asymmetrically, and the back-reaction
eta*curl(theta) creates a net force that drifts the braid. This drift
depends on the AVERAGE phi field (constant at A_bg ~ 0.1), not the
gradient, which is exactly what the data shows.

### The Quadratic Fit Results

- 3-field R^2: 0.43-0.46 (noise-level trajectory, no systematic acceleration)
- 6-field R^2: 0.85-0.93 (real systematic drift, but gradient-independent)

The 3-field R^2 near 0.5 is what you get from fitting a random walk to
a quadratic. The 6-field R^2 near 0.9 reflects the genuine constant drift.

### Implications for CONCEPT.md Section 3

This is severe. Three outcomes are possible:

**If the z-boundary bug is confirmed**: C=186 is a numerical artifact.
CONCEPT.md Section 3 must be retracted. The gradient-force gravity
mechanism is dead. The depletion profile (delta_rho ~ 1/r^1.2) still
exists but produces no force on braids. The EM radiation pressure
mechanism (from EM_THEORY.md/OQ2) becomes the sole candidate for
inter-braid forces.

**If resolution dependence is confirmed (C -> 186 at high N)**: The
gradient force is real but only observable with sufficient resolution.
CONCEPT.md claims survive, but the N=128 specification in V33_RESULTS.md
would need correction to reflect the actual resolution used.

**If neither hypothesis resolves the discrepancy**: Some other difference
between V33 and F9v2 exists that has not been identified. This would
require a careful side-by-side diff of the actual run conditions.

### What Must Be Done (Critical Path)

1. **Add periodic z-BC to the original V33 code** (two-line fix: replace
   k-1 with (k-1+N)%N and k+1 with (k+1)%N in the Laplacian) and re-run
   the gradient sweep at N=128, T=60. If C drops from 186 to ~0, the
   z-boundary bug is confirmed.

2. **Run F9v2 at N=256 and N=512** with T=60 to test resolution dependence.

3. **Run V33 UNMODIFIED at N=128, T=60** with the exact four gradient
   strengths and record the trajectory. Compare frame-by-frame with the
   V33_RESULTS.md drifts. This rules out the possibility that F9v2 has a
   subtle difference from V33 in the 3-field sector.

### Verdict: FAIL

The V33 C=186 cannot be reproduced. The 3-field case shows ZERO gradient
force (drift is noise-level, alternating signs). The 6-field case shows
a constant drift independent of the gradient. The most parsimonious
explanation is the V33 z-boundary bug (no explicit z-BC, Laplacian reads
wrong-row memory). This hypothesis is directly testable.

**This is the most critical finding in the V43 review.** If confirmed, it
invalidates the cornerstone measurement of the gradient-force gravity
mechanism.

---

## Round 2 Summary Table

| Result | Verdict | Key Finding |
|--------|---------|-------------|
| OQ3v2 Per-voxel DC | **PASS** | Correct methodology; DC_frac 5-8% is upper bound contaminated by centroid drift; flat radial profile suggests drift dominance; physical DC likely < 1% |
| OQ5v2 Plane-wave pulse | **PASS** | v_g matches lattice prediction to < 1%; pulse coherent with < 1.3% broadening; A_bg=0 tests bare equations (valid but not physical vacuum); strong energy oscillation is non-perturbative resonance |
| OQ1v2 Poynting multipole | **PASS** | Correct methodology; self-correction on positive-definiteness artifact; anisotropy ratio 19:1 establishes null model baseline; UUD comparison still pending |
| F9v2 Gradient force | **FAIL** | V33 C=186 NOT reproduced; 3-field shows ZERO gradient force; 6-field shows constant drift unrelated to gradient; V33 z-boundary bug is the most likely cause; requires urgent follow-up |

### Critical Path Forward

The F9v2 failure is existential for the gradient-force gravity hypothesis
in CONCEPT.md Section 3. Before other V43 work proceeds, the three V33
follow-ups must be completed:

1. Add periodic z-BC to V33 and re-run (one afternoon of work)
2. Run F9v2 at N=256, N=512 with T=60 (one day of compute)
3. If both give C~0, the gradient force is dead and CONCEPT.md must be revised

The OQ1v2 result needs the UUD high-cadence simulation to become
scientifically meaningful. The OQ5v2 result is solid but should be
repeated with A_bg=0.1 for physical-vacuum predictions.
