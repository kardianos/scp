# Corrected Experimental Plans — V43 Redo

These plans address the four results that FAILED or drew CONCERN in
the skeptic review (`skeptic_review.md`). Each plan specifies exactly
what to build, what data to use, what to measure, and how to know if
the result is real or an artifact.

---

## OQ1 Corrected: Poynting Flux Multipole Decomposition

### Why the Original Failed

The original decomposed |theta| — a positive-definite scalar — into
spherical harmonics. Any positive-definite function has a large l=0
component by construction. The 89% monopole fraction is a mathematical
triviality, not a physical result. The physically relevant quantity
for the radiation-pressure Coulomb argument is the TIME-AVERAGED
Poynting flux magnitude on spherical shells.

### The Right Quantity

The Cosserat Poynting vector is S = E x B = (-dtheta/dt) x (curl theta).
We need to decompose **<|S|>(Omega)** — the time-averaged Poynting flux
magnitude as a function of angle on a sphere of radius R — into
spherical harmonics Y_l^m.

### Input Data

**UUD composite** (test case):
- File: `/home/d/code/scp/v41/results/stable/UUD_stable_f16.sfa`
- Grid: N=192, L=30, 12 frames, T=0..500, snap_dt=45.45
- 12 columns: phi_xyz, theta_xyz, phi_v_xyz, theta_v_xyz (all f16)
- Problem: only 12 frames at 45.45 interval. The theta oscillation
  period is ~4.4t, so consecutive frames are ~10 oscillation periods
  apart. This is FAR too sparse for computing dtheta/dt via finite
  differences — the Nyquist criterion requires at least 2 samples per
  period (snap_dt < 2.2t).

**Single braid** (null model):
- File: `/home/d/code/scp/v34/torsion_coupling/data/braid_hires.sfa`
- Grid: N=80, L=25, 264 frames, T=0..50, snap_dt=0.19
- 6 columns: phi_xyz, theta_xyz (f64)
- This HAS adequate temporal resolution (snap_dt=0.19 << 2.2t period)

**Resolution**: We need a new UUD simulation with high temporal cadence.
The single-braid hires file works as the null model.

### Step 0: Generate a High-Cadence UUD SFA

Run the UUD composite for T=50 with snap_dt=0.19 (matching braid_hires).
This gives ~263 frames over ~11 oscillation periods.

```
# Use the existing phase-confined seed generator
cd /home/d/code/scp/v41
./gen_phase_confined -chirality UUD -N 192 -L 30 -o ../v43/oq1/uud_seed.sfa -precision f32

# Then run the sim kernel (CPU is fine for N=192, T=50)
cd /home/d/code/scp/sfa/sim
./scp_sim braid_default.cfg \
    -N 192 -L 30 -T 50 -snap_dt 0.19 \
    -init sfa -seed_file ../../v43/oq1/uud_seed.sfa \
    -output ../../v43/oq1/uud_hires.sfa \
    -precision f32 -damp_width 3.0 -damp_rate 0.01
```

Expected size: N=192, 6 cols f32, 263 frames. Per frame: 192^3 * 6 * 4 =
~170 MB uncompressed, ~90 MB compressed. Total: ~24 GB. This is large
but feasible. If too large, use N=128 with L=20 (still captures the
composite at the center, frames drop to ~50 MB each, total ~13 GB).

Alternative: if the existing UUD_stable_f16.sfa was generated from a
seed that is still available, restart from the LAST frame (t=500, the
equilibrated state) and run for T=50 more with snap_dt=0.19. This
avoids the transient relaxation phase. Check:
`/home/d/code/scp/v41/seeds/UUD_phase.sfa`

Expected runtime: N=192 on CPU with 8 threads: ~3 ms/step, dt~0.008,
T=50 → 6250 steps → ~19 seconds. Very fast.

### Step 1: Compute the Poynting Flux at Each Frame Pair

For each consecutive frame pair (n, n+1) at times t_n, t_{n+1}:

1. Read theta_x, theta_y, theta_z at both frames
2. Compute E = -dtheta/dt via finite difference:
   ```
   E_a(x) = -(theta_a(x, t_{n+1}) - theta_a(x, t_n)) / (t_{n+1} - t_n)
   ```
3. Compute B = curl(theta) using centered spatial differences at frame n
   (or the average of frames n and n+1):
   ```
   B_x = dtheta_z/dy - dtheta_y/dz
   B_y = dtheta_x/dz - dtheta_z/dx
   B_z = dtheta_y/dx - dtheta_x/dy
   ```
   Use standard second-order central differences: df/dx = (f[i+1]-f[i-1])/(2*dx)
4. Compute S = E x B at each grid point:
   ```
   S_x = E_y*B_z - E_z*B_y
   S_y = E_z*B_x - E_x*B_z
   S_z = E_x*B_y - E_y*B_x
   ```
5. Compute |S| = sqrt(S_x^2 + S_y^2 + S_z^2) at each grid point

### Step 2: Interpolate onto Spherical Shells

For each shell radius R in {4, 6, 8, 10, 12}:

1. Define an angular grid: N_theta = 90, N_phi = 180 (2-degree resolution,
   adequate for l up to ~20)
2. For each (theta_i, phi_j), compute the Cartesian point:
   ```
   x = R*sin(theta_i)*cos(phi_j) + cx
   y = R*sin(theta_i)*sin(phi_j) + cy
   z = R*cos(theta_i) + cz
   ```
   where (cx, cy, cz) is the centroid of the composite (use |P|-weighted
   centroid from phi fields)
3. Trilinear interpolation from the grid to get |S|(theta_i, phi_j)
4. Store the angular map for this frame pair

### Step 3: Time-Average

Average the angular maps over ALL frame pairs. With 263 frames we get
262 Poynting maps. The average spans T=50, which is ~11 oscillation
periods — adequate for convergence.

Result: <|S|>(theta, phi; R) for each shell radius R.

### Step 4: Spherical Harmonic Decomposition

Expand <|S|>(theta, phi) in real spherical harmonics:

```
<|S|>(theta, phi) = sum_{l=0}^{l_max} sum_{m=-l}^{l} c_{lm} Y_l^m(theta, phi)
```

Compute the power spectrum:
```
C_l = (1/(2l+1)) sum_{m=-l}^{l} |c_{lm}|^2
```

Monopole fraction: f_0 = C_0 / sum_l C_l

Use l_max = 8 (sufficient — we do not expect structure beyond octupole
for three orthogonal dipoles).

The integral c_{lm} = integral <|S|> Y_l^m sin(theta) d(theta) d(phi)
is computed numerically via the angular grid (trapezoidal rule in phi,
Gauss-Legendre in theta for better accuracy, or just trapezoidal in
both since N_theta=90 is overkill for l<=8).

### Step 5: Null Model

Repeat Steps 1-4 on the single-braid data:
- File: `braid_hires.sfa` (N=80, L=25, 264 frames, f64)
- Same radii: R = {4, 6, 8, 10, 12}
- The single braid is aligned along z, so it should show a dipolar
  radiation pattern (l=1 dominant) with perhaps some l=2 from the
  elliptical cross-section

### Measurements to Report

For each shell radius R, for both UUD and single-braid:
1. C_l for l = 0, 1, 2, 3, 4
2. Monopole fraction f_0 = C_0 / sum C_l
3. Anisotropy ratio: max(<|S|>) / min(<|S|>) on each shell
4. Total Poynting flux: integral <|S|> dOmega (should fall as ~1/R^2
   if radiative; this is a separate check)

### Success Criteria

- **Pass**: UUD monopole fraction f_0 > 0.5 at R >= 8, AND single-braid
  monopole fraction f_0 < 0.2 at the same R. The UUD must be
  SIGNIFICANTLY more isotropic than a single braid.
- **Marginal**: UUD f_0 between 0.2 and 0.5. The radiation pattern is
  partially isotropized but still has significant multipole structure.
- **Fail**: UUD f_0 comparable to single-braid f_0 (within 50%). The
  three-dipole summation does not produce effective isotropy.

### Pitfalls to Avoid

1. **Do NOT decompose |S|** — decompose <|S|>_time (the time average).
   The time average removes oscillatory angular structure and isolates
   the secular radiation pattern.
2. **Do NOT use |theta| or theta_rms** — the Poynting flux S = E x B is
   the physical energy flow. |theta| is the potential, not the flux.
3. **f16 precision**: The UUD_stable file is f16. For computing spatial
   derivatives (curl), f16 gives only ~3 significant digits. At dx=0.31,
   the curl is ~(theta[i+1]-theta[i-1])/(2*0.31). If theta~0.04, the
   numerator is ~0.001, and f16 precision is ~4e-5, giving ~4% noise.
   This is marginal. Use f32 data if possible (run a new simulation).
   The braid_hires is f64, so the null model has no precision issue.
4. **PBC at large R**: At N=192, L=30, periodic images are at distance
   60. Data at R=12 is 40% of the way to the image. For R <= 12, PBC
   contamination should be < 5%. Do not trust R > 15.
5. **Centroid drift**: The composite may drift slightly over T=50. Track
   the centroid at each frame and re-center. The drift should be small
   (< 1 grid cell) for an equilibrated composite.

### Expected Runtime

- New UUD hires simulation: ~20 seconds on CPU (N=192, T=50)
- Poynting computation per frame pair: N=192 grid, 7M voxels, ~5 ops
  per voxel = 35M ops ~ 0.01s. Times 262 pairs = ~3 seconds
- Angular interpolation: 5 radii x 90x180 = 81K points per frame pair,
  262 pairs = ~21M interpolations. ~10 seconds
- SH decomposition: trivial (5 radii x 81 coefficients)
- Null model on braid_hires: same but N=80, even faster
- **Total: under 5 minutes including the simulation**

---

## OQ5 Corrected: Plane-Wave Pulse Propagation

### Why the Original Failed

The original used a 3D spherical Gaussian as the initial condition.
Under the wave equation, a spherical pulse does NOT propagate as a
coherent packet — it decomposes into outgoing and ingoing spherical
waves that spread geometrically. The measured "velocity" was the
centroid of a dispersing shell, which slowed from ~1c to ~0.5c purely
due to geometry, not physics.

### The Right Initial Condition

A plane-wave pulse: a Gaussian envelope in x modulating a carrier wave.
The pulse is uniform in y and z (plane wave), so there is no transverse
spreading. The only broadening comes from dispersion.

```
theta_z(x,y,z,0) = A * exp(-(x - x0)^2 / (2*sigma_x^2)) * sin(k0 * x)

theta_z_dot(x,y,z,0) = -v_g * A * d/dx[exp(-(x-x0)^2/(2*sigma_x^2)) * sin(k0*x)]
    = -v_g * A * exp(-(x-x0)^2/(2*sigma_x^2))
        * [-(x-x0)/sigma_x^2 * sin(k0*x) + k0 * cos(k0*x)]
```

where v_g is the predicted group velocity (see below).

All other field components are zero: theta_x = theta_y = 0, all phi = 0.

### Parameter Choices

**Carrier wavenumber**: k0 = pi (wavelength lambda = 2 in code units).
This puts the pulse well above the phi mass gap (m=1.5, so
k_mass = m = 1.5 < k0 = 3.14). The pulse is in the "photon-like"
regime.

**Envelope width**: sigma_x = 5*lambda/(2*pi) = 5/k0 = 5/pi ≈ 1.59.
This gives k0*sigma_x = 5, so ~5 wavelengths under the envelope
(narrowband enough for a well-defined group velocity).

Actually, for better separation: use sigma_x = 5*lambda = 10 code units.
Then k0*sigma_x = pi*10 ≈ 31, giving a very narrowband pulse.
The pulse FWHM in x is 2.35*sigma_x ≈ 23.5 code units.

**Amplitude**: A = 0.01 (small, to stay in the linear regime; the
background has A_bg = 0.1, so this is a 10% perturbation of the
background theta — but theta starts at zero, so there is no background
theta. Use A = 0.01 to be safe in the linear regime.)

**Initial position**: x0 = -L_x/2 + 2*sigma_x (placed near the left
edge with 2-sigma clearance)

### Predicted Group Velocity

The linearized dispersion relation for the coupled phi-theta system
(from EM_THEORY.md Section 2) has two branches. For the photon-like
branch at k = k0:

```
omega^2 = k0^2 + eta^2 * k0^2 / (k0^2 + m^2)    (approximate, for small eta)
```

More precisely, the coupled dispersion relation is:
```
(omega^2 - k^2)(omega^2 - k^2 - m^2) = eta^2 * k^2
```

This is a quadratic in omega^2. The two solutions:
```
omega^2 = k^2 + m^2/2 +/- sqrt(m^4/4 + eta^2*k^2)
```

The photon-like branch (lower omega, the minus sign):
```
omega_photon^2 = k^2 + m^2/2 - sqrt(m^4/4 + eta^2*k^2)
```

At k = k0 = pi, m^2 = 2.25, eta = 0.5:
- m^4/4 = 1.265625
- eta^2 * k0^2 = 0.25 * 9.87 = 2.4674
- sqrt(1.2656 + 2.4674) = sqrt(3.7330) = 1.9322
- omega^2 = 9.8696 + 1.125 - 1.9322 = 9.0624
- omega = 3.0104
- v_phase = omega/k = 3.0104/3.1416 = 0.9583
- v_group = d(omega)/dk. Differentiating:
  2*omega * d(omega)/dk = 2k + (eta^2 * k) / sqrt(m^4/4 + eta^2*k^2)
  d(omega)/dk = [2k + eta^2*k/sqrt(...))] / (2*omega)

  Numerator: 2*3.1416 + 0.25*3.1416/1.9322 = 6.2832 + 0.4066 = 6.6898
  v_group = 6.6898 / (2*3.0104) = 1.1108

Wait — group velocity > c? That would mean the coupling INCREASES the
group velocity above 1. Let me re-check.

Actually, for the uncoupled theta equation (eta=0), omega = k (massless),
so v_g = 1 exactly. With eta > 0, the coupling pushes the photon branch
DOWN in omega (the minus sign in the dispersion), meaning slower phase
velocity but the group velocity needs careful computation. Let me redo:

At eta=0: omega^2 = k^2 (photon), omega^2 = k^2 + m^2 (phonon).
At finite eta, the two branches repel. The photon branch gets pushed DOWN.

The exact dispersion: (omega^2 - k^2)(omega^2 - k^2 - m^2) = eta^2*k^2

Let u = omega^2 - k^2. Then u(u - m^2) = eta^2*k^2, giving
u = m^2/2 +/- sqrt(m^4/4 + eta^2*k^2).

Photon branch: u_- = m^2/2 - sqrt(m^4/4 + eta^2*k^2) < 0.
So omega^2 = k^2 + u_- < k^2. The phase velocity is < c.

v_group = d(omega)/dk. With omega^2 = k^2 + m^2/2 - sqrt(m^4/4 + eta^2*k^2):
d(omega^2)/dk = 2k - eta^2*k / sqrt(m^4/4 + eta^2*k^2)
2*omega * v_g = 2k - eta^2*k / sqrt(...)
v_g = k/omega * [1 - eta^2/(2*sqrt(m^4/4 + eta^2*k^2))]

At k=pi: k/omega = 3.1416/sqrt(9.0624) = 3.1416/3.0104 = 1.0436
Correction: eta^2/(2*sqrt(...)) = 0.25/(2*1.9322) = 0.0647
v_g = 1.0436 * (1 - 0.0647) = 1.0436 * 0.9353 = 0.9762

So v_g ≈ 0.976c for k0=pi, eta=0.5. This is the prediction.

For eta=0: v_g = 1.000c exactly (massless wave equation).

### Grid Design

The pulse propagates in x. It does not need transverse resolution:

```
N_x = 512, N_y = 16, N_z = 16
L_x = 80, L_y = 2.5, L_z = 2.5
```

Wait — the sim kernel uses cubic grids (N, L are single values).
Check the kernel to confirm. If cubic-only, use:

```
N = 256, L = 50
```

This gives dx = 2*50/255 = 0.392. The carrier wavelength is 2.0,
so we get ~5 grid points per wavelength — MARGINAL. Better:

```
N = 512, L = 50
```

dx = 0.196, ~10 points per wavelength — adequate.

But N=512 with all phi fields is (512)^3 * 12 * 8 bytes = 12.9 GB per
frame for f64. This is too large for a simple test.

**Better approach**: Since this is a plane wave, the y and z dependence
is trivial (constant). But the 3D kernel will still allocate the full
grid. Use the elongated approach IF the kernel supports non-cubic grids.
If not, use N=256, L=50 and accept ~5 pts/wavelength.

Actually, the Poynting test only needs the theta sector. With k0=pi and
sigma_x=10, the pulse occupies ~40 code units in x. At L=50 (domain
[-50,50]), the pulse starts at x0=-30, propagates to x=+30 in T=60/c=60.
With absorbing boundaries at the edges, the pulse can run for ~60 time
units before hitting the boundary. This is enough.

**Revised grid**: N=256, L=50, T=80.
dx = 2*50/255 = 0.392. Points per wavelength: 2.0/0.392 = 5.1.
This is adequate for 2nd-order finite differences (error ~(k*dx)^2/12
≈ 0.13 = 13% per wavelength — not great). With 4th-order stencils in
the kernel, this improves to ~0.2%.

Check: does the kernel use 2nd or 4th order Laplacian? The reference
implementation `v33_cosserat.c` likely uses 7-point (2nd order).
scp_sim.c probably uses the same. With 5 pts/wavelength and 2nd order,
numerical dispersion is ~5%. This is a KNOWN systematic. Either:
(a) increase to N=384 (7.5 pts/wavelength, ~2% numerical dispersion), or
(b) use a higher carrier: k0 = pi/2 (wavelength=4, 10 pts/wavelength)

**Final choice**: k0 = pi/2 (wavelength = 4, well above mass gap since
k0 = 1.57 > m = 1.5... barely. Actually k0=1.57 is almost exactly at
the mass gap. Use k0 = 2 instead (wavelength = pi ≈ 3.14, 8 pts/
wavelength at N=256 L=50). k0=2 > m=1.5, safely above gap.)

**Revised parameters**:
- k0 = 2.0 (wavelength = pi ≈ 3.14)
- sigma_x = 15.0 (k0*sigma_x = 30, very narrowband)
- N = 256, L = 50, T = 80
- dx = 0.392, pts/wavelength = 3.14/0.392 = 8.0 — adequate
- A = 0.01
- x0 = -30.0

Predicted group velocity at k0=2.0:
- m^4/4 = 1.2656, eta^2*k0^2 = 0.25*4 = 1.0
- sqrt(1.2656 + 1.0) = sqrt(2.2656) = 1.5052
- omega^2 = 4.0 + 1.125 - 1.5052 = 3.6198, omega = 1.9026
- v_phase = 1.9026/2.0 = 0.9513
- v_group = (k/omega)*[1 - eta^2/(2*sqrt(...))]
  = (2.0/1.9026)*[1 - 0.25/(2*1.5052)]
  = 1.0512 * (1 - 0.0831) = 1.0512 * 0.9169 = 0.9639

So v_g ≈ 0.964c with eta=0.5.

### Seed Generator

Write a seed generator `gen_planewave.c` in `/home/d/code/scp/sfa/seed/`:

```c
// For each grid point (i,j,k):
double x = -L + i * dx;  // (use the kernel's convention)
double env = exp(-(x - x0)*(x - x0) / (2.0 * sigma_x * sigma_x));
double theta_z     = A * env * sin(k0 * x);
double theta_z_dot = -v_g * A * env * (-(x-x0)/(sigma_x*sigma_x) * sin(k0*x) + k0 * cos(k0*x));
// All other fields = A_bg * cos(k_bg * z + delta_a) for phi, 0 for theta_x, theta_y
```

The phi fields are initialized to the standard background:
phi_a = A_bg * cos(k*z + delta_a) where k = 2*pi*m/L (one wavelength
across the box? No — k should give the background oscillation. Use the
same initialization as in the sim kernel's braid init, minus the braid:
just the A_bg * cos(k*z + delta_a) background.)

### Two Runs

**Run A (uncoupled, eta=0)**: v_g should be exactly 1.0 (pure wave equation).
This is the CONTROL.

**Run B (coupled, eta=0.5)**: v_g should be ~0.964c.

Use identical initial conditions for both (same seed SFA), just change
eta in the config file.

### Measurements

1. **Envelope peak position x_peak(t)**: At each output frame, project
   theta_z onto the x-axis by summing theta_z^2 over y,z planes:
   ```
   P(x, t) = sum_{j,k} theta_z(i,j,k,t)^2
   ```
   Then fit a Gaussian to P(x) to find the peak position and width.

2. **Group velocity**: Linear fit of x_peak(t) vs t. Measure v_g with
   uncertainty from the fit residuals.

3. **Envelope width sigma(t)**: From the Gaussian fit at each frame.
   For a non-dispersive wave, sigma(t) = sigma_x = constant. Dispersion
   causes broadening: sigma(t)^2 = sigma_x^2 + (d^2omega/dk^2)^2 * t^2 / (4*sigma_x^2).

4. **Phi energy fraction**: E_phi(t) / E_total(t). Should be ~0 for
   eta=0 and should grow to a steady value for eta=0.5. Predict the
   steady-state fraction from the linearized mixing:
   ```
   E_phi/E_theta ≈ eta^2 * k0^2 / (k0^2 + m^2)^2
   ```
   At k0=2: 0.25*4/(4+2.25)^2 = 1.0/39.06 = 0.026 = 2.6%.

5. **Carrier wavelength check**: Measure the wavelength from zero-crossings
   of theta_z along the x-axis at t=0 (should give lambda=pi) and at
   the last frame (should be unchanged if dispersion is small).

### Output Configuration

- snap_dt = 1.0 (80 frames over T=80)
- precision = f32 (adequate for this test)
- N=256 cube: each frame is 256^3 * 12 * 4 = 805 MB uncompressed.
  With compression, ~400 MB. Total for 80 frames: ~32 GB.

This is too large. Solutions:
(a) Output only every 5th frame (snap_dt=5.0, 16 frames, ~6.4 GB total)
(b) Use N=128, L=25 (pts/wavelength = 4 for k0=2 — marginal, but with
    the kernel's 7-point stencil, numerical dispersion is ~10%)
(c) Use the diag file for energy tracking (no need for full snapshots)
    and output only 5-10 frames for the Gaussian fit

**Best approach**: Output snap_dt=10.0 (8 frames) for the spatial Gaussian
fits, plus diag_dt=0.5 for continuous energy tracking. Total SFA: ~3 GB.

### Success Criteria

- **Pass (eta=0)**: v_g = 1.000 +/- 0.005. Sigma grows by < 5% over
  T=80. Confirms the massless theta wave equation works.
- **Pass (eta=0.5)**: v_g within 2% of the predicted 0.964c. Sigma
  broadening consistent with the predicted d^2omega/dk^2. E_phi/E_total
  matches the mixing prediction within 30%.
- **Fail**: v_g for eta=0 deviates from 1.0 by more than 1% (indicates
  numerical dispersion or a bug). Or v_g for eta=0.5 disagrees with
  the prediction by more than 5%.

### Pitfalls to Avoid

1. **Do NOT use a 3D Gaussian** — it spreads geometrically and the
   centroid measurement is meaningless after t ~ sigma_x/c.
2. **Carrier wavevector must be well above mass gap**: k0 > m = 1.5.
   If k0 ~ m, the photon-like branch merges with the massive branch
   and the mode identity is unclear.
3. **Numerical dispersion**: At 5 pts/wavelength with 2nd-order stencils,
   numerical dispersion is ~5%. This is a systematic error on v_g.
   Correct for it by computing the NUMERICAL dispersion relation
   (omega^2 = (2/dx^2)*[3 - 4*cos(k*dx) + cos(2*k*dx)]/... for 7-point
   stencil) and predicting v_g from that. Or simply compare Run A
   (eta=0) with the known answer v_g=1 and use the DIFFERENCE between
   Run A and Run B as the physical measurement.
4. **Absorbing boundaries**: Enable damp_width and damp_rate. Without
   damping, reflections from PBC will contaminate the measurement
   when the pulse reaches the boundary (at t ~ L/c ~ 50).
5. **Background phi initialization**: The phi fields MUST have the
   standard A_bg=0.1 background. Without it, the V(P) potential is
   inactive and the phi-theta coupling may behave differently than in
   the braid simulations. However, for the pure linearized test, the
   background provides the medium through which the mixing operates.

### Expected Runtime

- N=256, T=80: ~50000 steps at dt~0.002. At ~10 ms/step (N=256 cube
  with OpenMP): ~500 seconds = 8 minutes.
- N=128: ~12500 steps at dt~0.005, ~1 ms/step: ~12 seconds.
- Analysis (Gaussian fits on 8 frames): seconds.

---

## OQ3 Corrected: Per-Voxel DC Extraction

### Why the Original Was Problematic

The original `dc_profile.c` computed the time-AND-shell averaged theta
vector, |<theta>_{t,shell}|. This conflates temporal DC content with
spatial averaging over each shell. A dipolar DC pattern (expected for a
z-aligned braid) would partially cancel in the shell average, giving a
number that is neither the true DC amplitude nor zero. The result (5-12%
DC/RMS in the core) is a lower bound but is not comparable to V34's 0.2%
because V34 measured a different quantity.

### The Right Measurement: Per-Voxel Time Average

For each voxel at grid position (i,j,k), compute the time-averaged
theta vector:

```
<theta_a>(i,j,k) = (1/N_frames) * sum_{n=0}^{N_frames-1} theta_a(i,j,k, t_n)
```

This is a vector field: the true DC component of theta at each point in
space. It is NOT shell-averaged.

Then compute derived quantities:
- DC amplitude at each voxel: |<theta>|(i,j,k) = sqrt(<theta_x>^2 + <theta_y>^2 + <theta_z>^2)
- AC amplitude at each voxel: first compute AC_rms_a^2 = <theta_a^2> - <theta_a>^2, then
  AC_rms(i,j,k) = sqrt(AC_rms_x^2 + AC_rms_y^2 + AC_rms_z^2)
- DC/AC ratio at each voxel: |<theta>| / AC_rms

### Input Data

**Primary**: `/home/d/code/scp/v34/torsion_coupling/data/braid_hires.sfa`
- N=80, L=25, 264 frames, snap_dt=0.19, f64
- Single braid along z axis
- 264 frames spanning T=50, covering ~11 oscillation periods (period ~4.4t)
- f64 precision: no numerical noise issue

### Algorithm

1. **Allocate accumulator arrays** (per-voxel):
   - sum_theta_a[3][N^3] (double, for signed sum)
   - sum_theta2_a[3][N^3] (double, for squared sum)
   - Memory: 80^3 * 6 * 8 = 24.6 MB. Trivial.

2. **First pass — accumulate**: For each of the 264 frames:
   - Read theta_x, theta_y, theta_z
   - For each voxel: sum_theta_a += theta_a, sum_theta2_a += theta_a^2

3. **Compute per-voxel DC and AC**:
   ```
   DC_a(i,j,k) = sum_theta_a[a] / 264
   DC_mag(i,j,k) = sqrt(DC_x^2 + DC_y^2 + DC_z^2)
   AC_var_a(i,j,k) = sum_theta2_a[a]/264 - DC_a^2
   AC_rms(i,j,k) = sqrt(AC_var_x + AC_var_y + AC_var_z)
   ```

4. **Centroid tracking**: Compute the |P|-weighted centroid from the phi
   fields at frame 0. For N=80 with a single braid, the centroid should
   be near (0,0,0). Verify it doesn't drift significantly by checking
   at frames 0, 132, 263. If drift > 1 grid cell (0.63 code units),
   re-center at each frame before accumulating (shift coordinates, not
   the data — use floating-point sub-voxel centroid and adjust the
   radial binning accordingly). For braid_hires, the braid is stable
   over T=50 and drift should be negligible.

5. **Radial profiles** (shell-averaged but NOW of per-voxel quantities):
   ```
   DC_rms(r) = sqrt( <DC_mag^2>_shell )     [RMS of per-voxel DC over shell]
   AC_rms_shell(r) = sqrt( <AC_rms^2>_shell )
   DC_frac(r) = DC_rms(r) / sqrt(DC_rms(r)^2 + AC_rms_shell(r)^2)
   ```
   Use radial bins of width dr = 0.5 (matching the original analysis).

6. **Angular structure of DC field**: On shells at r = {2, 4, 6, 8, 10},
   interpolate DC_x, DC_y, DC_z onto an angular grid (30x60 in
   theta,phi). Compute the spherical harmonic decomposition of each
   component. For a z-aligned braid, expect:
   - DC_x and DC_y: dipolar (l=1) or higher, perpendicular to z
   - DC_z: small (braid is symmetric about z)
   Report the dominant multipole l for each component.

### Restricted Analysis Domain

Only report results for r < 12 (PBC-free zone). At N=80, L=25, periodic
images are at distance 50. Contamination at r=12 is from the nearest
image at distance 50-12=38, contributing ~exp(-m*38)/38 ~ negligible
for the field, but the RADIATION from images accumulates. Conservative
cutoff: r < L/2 = 12.5. Use r < 12.

### Measurements to Report

1. Table: r, DC_rms(r), AC_rms(r), DC_frac(r), DC_rms/total_rms for
   r = 0.5 to 12 in steps of 0.5
2. Angular maps of DC_vec at r = 2, 4, 6, 8, 10 (visualize as vector
   field on sphere or as component heatmaps)
3. Dominant multipole of DC field at each radius
4. Global DC fraction: sum of DC energy / sum of total theta energy
   within r < 12
5. Comparison number: what is DC_frac at r ≈ 8-10? This is the
   "far-field DC fraction" in the reliable domain. Compare to V34's
   0.2% claim. If DC_frac(r=10) ≈ 0.002, V34 was correct. If
   DC_frac(r=10) ≈ 0.05-0.10, the core result (5-12%) extends outward.

### Success Criteria

There is no pass/fail here — this is a MEASUREMENT, not a hypothesis
test. The goal is to produce the definitive DC profile that resolves
the discrepancy between V34 (0.2%) and the original V43 analysis
(5-12%). The per-voxel method is the correct methodology; the number
it produces is the answer.

However, if DC_frac is CONSTANT with radius (independent of r for
r > 4), that strongly suggests it is a numerical artifact of the
initialization or the grid (systematic bias), not a physical DC
component. A physical DC should have a specific radial and angular
structure related to the source (the braid's current loop).

### Pitfalls to Avoid

1. **Do NOT shell-average the signed theta before computing the DC**.
   Shell-averaging cancels dipolar structure. Compute DC per-voxel
   first, then take statistics over shells.
2. **Centroid drift**: Even small drift (0.5 grid cells over 264 frames)
   would smear a sharp DC structure. Check and correct if needed.
3. **Transient frames**: The braid may still be relaxing in the first
   few frames (t < 5). Optionally, exclude the first 25 frames (t < 5)
   and use only frames 26-264 for the time average. Compare the
   result with/without the first 25 frames to check sensitivity.
4. **Grid convention**: Use dx = 2*L/(N-1) = 2*25/79 = 0.6329 to match
   the simulation kernel (vertex-centered grid), NOT dx = 2*L/N = 0.625.
   The original code had this wrong (1.3% coordinate error). Small but
   fix it.

### Expected Runtime

- Reading 264 frames of N=80 f64 data: 264 * 24.6 MB = 6.5 GB total
  I/O. At ~500 MB/s SSD read: ~13 seconds.
- Accumulation: 264 * 512K voxels * 6 ops = ~800M ops. Trivial.
- Angular interpolation: 5 radii * 1800 points * 3 components = 27K
  interpolations. Trivial.
- **Total: under 30 seconds**

---

## F9 Corrected: Numerical Force Integral

### Why the Original Was Problematic

The analytical estimate C ~ |mu| * R_tube^3 * psi^2 / 2 has one free
parameter: the radius at which to evaluate psi. Different radii give
C ranging from 22 to 1190. The "prediction" C=127 chose R_tube=4.5 to
get a plausible answer. The exponent prediction (2.2 vs measured 1.8)
is a genuine discrepancy, not rough agreement.

### The Right Approach: Numerical Integration

The force on a braid in a density gradient is (from the skeptic review's
Eq. 18):

```
F = -integral (d m_eff^2/dx) * |psi|^2 d^3x
```

where m_eff^2(x) = m^2 + V''(P_bg(x)) is the local effective mass
determined by the background field, and |psi|^2 is the braid's field
intensity.

For a single braid in a linear density gradient, the force is F = -C * grad(rho).
We can extract C by:
1. Computing the braid's field profile from simulation data
2. Imposing an analytical density gradient
3. Numerically integrating the force integral
4. Comparing C to the measured value from V33 (C = 186)

### Input Data

**Braid profile**: `/home/d/code/scp/v34/torsion_coupling/data/braid_hires.sfa`
- N=80, L=25, 264 frames, f64
- Use the LAST frame (frame 263, t≈49.9) — the braid is fully relaxed

**V33 measurement**: C = 186 +/- 3 (from the gradient sweep in V33,
documented in CONCEPT.md Section 3). The force exponent is 1.8 (measured
over D=15-30).

### Algorithm

**Step 1: Extract the equilibrium braid profile**

Read the last frame of braid_hires.sfa. Extract phi_0, phi_1, phi_2
on the 80^3 grid. Locate the braid centroid (|P|-weighted). The braid
is centered near (0,0,0) and aligned along z.

**Step 2: Compute |psi|^2 and P at each voxel**

```
|psi|^2(x) = phi_0^2 + phi_1^2 + phi_2^2
P(x) = phi_0 * phi_1 * phi_2
```

**Step 3: Compute m_eff^2(x) in the UNPERTURBED background**

The potential is V(P) = (mu/2) * P^2 / (1 + kappa*P^2).

The second derivative with respect to a single phi component (say phi_0)
involves the chain rule through P:

Actually, the force equation F = -integral (dm_eff^2/dx)|psi|^2 d^3x
needs more careful derivation. The effective mass for a test particle
in the background of the braid is:

```
m_eff^2(x) = m^2 + V''(P_bg(x))
```

But V''(P) means the second derivative of V with respect to the field,
not with respect to P. Since V depends on phi_a through P = phi_0*phi_1*phi_2,
the Hessian is complicated. The simpler formulation:

The gravitational force on a braid comes from the GRADIENT of the
background density rho_bg. The braid's self-energy depends on the
local background density because m_eff^2 depends on the local field
amplitudes.

**Cleaner formulation**: The gradient force arises because the braid's
perturbation profile (its Yukawa tail) extends further on the depleted
side than the dense side. The force is:

```
F = integral f_a(x) * (d phi_a^bg / dx) d^3x
```

where f_a = laplacian(phi_a) - m^2*phi_a - dV/dP * dP/dphi_a is the
residual force (which is zero for an exact stationary solution, but
nonzero when the background is perturbed).

**Simplest correct approach**: Instead of analytically deriving the
force integral, MEASURE C directly from the simulation by the following
numerical experiment:

**Step 3 (revised): Direct force measurement on a known gradient**

Place the braid profile from braid_hires into a tilted background:
```
phi_a(x) = phi_a^braid(x) + A_bg * (1 + alpha*x) * cos(k*z + delta_a)
```

where alpha is a small fractional gradient (alpha = 0.01/L gives a 2%
density variation across the box). Then:

1. Compute the equation of motion acceleration d^2phi/dt^2 at each voxel
   using the FULL Cosserat equation (Laplacian, mass, potential, curl)
2. The net force on the braid is F_x = integral rho_braid * a_x d^3x
   where rho_braid is the excess density and a_x is the x-acceleration.
   More precisely:
   ```
   F_x = integral phi_a * (d^2phi_a/dt^2) d^3x   (sum over a)
   ```
   evaluated at t=0 (before the braid has moved).

Actually this is getting complicated. Let me use the cleanest approach:

### Revised Algorithm: Force from the Equation of Motion

The most direct approach computes the instantaneous force on the braid
from the field equations, with NO free parameters.

1. Read the last frame of braid_hires.sfa: phi_a(x), theta_a(x) on N=80.

2. Impose a LINEAR density gradient by SCALING the phi fields:
   ```
   phi_a'(x) = phi_a(x) * (1 + alpha * (x - cx))
   ```
   where alpha is a small gradient parameter (e.g., alpha = 0.001 per
   code unit, giving 5% variation over 50 code units). cx is the braid
   centroid.

3. Compute the acceleration field from the full equation of motion:
   ```
   a_phi_a(x) = laplacian(phi_a') - m^2*phi_a' - dV/dphi_a(phi_a') + eta*curl(theta)_a
   ```
   Use the same Laplacian stencil as the simulation kernel (7-point for
   2nd order, 19-point for 4th order — check which the kernel uses).

4. Compute the net force on the braid:
   ```
   F_x = sum_a integral phi_a'(x) * a_phi_a(x) d^3x
   ```
   restricted to the braid region (|r - centroid| < R_cut, where
   R_cut = 10 to capture the full braid profile but exclude far-field
   artifacts). Actually, the correct quantity is:
   ```
   F_x = sum_a integral a_phi_a(x) * x-momentum-weight d^3x
   ```
   The x-component of the total field momentum is:
   ```
   P_x = -integral sum_a (dphi_a/dt) * (dphi_a/dx) d^3x
   ```
   and dP_x/dt = F_x is the force. At t=0, dphi/dt is known (from the
   velocity fields, or zero for a stationary braid). The force is:
   ```
   F_x = -integral sum_a (d^2phi_a/dt^2) * (dphi_a/dx) + (dphi_a/dt) * d/dt(dphi_a/dx) d^3x
   ```

   For a stationary braid at t=0 (dphi/dt = 0), this simplifies to:
   ```
   F_x = -integral sum_a a_phi_a * (dphi_a'/dx) d^3x
   ```

   Here a_phi_a is the acceleration computed from the perturbed field phi_a'.

5. The force should be proportional to alpha (linear response). Verify
   by running at alpha = {0.0005, 0.001, 0.002, 0.004} and checking
   linearity.

6. Extract C = F_x / (alpha * dphi/dx_bg) where dphi/dx_bg is the
   background density gradient corresponding to the applied alpha.

Actually — the measurement from V33 was simpler. V33 placed a braid
in an externally imposed density gradient and measured the braid drift.
C = 186 is the proportionality constant in F = -C * grad(rho).

The cleanest numerical test is:

### Final Algorithm: Compute C from Instantaneous Force

1. Read last frame of braid_hires.sfa (the equilibrium braid).

2. For each voxel, compute the RHS of the phi equation of motion:
   ```
   RHS_a = laplacian(phi_a) - m^2 * phi_a - dV/dphi_a + eta * curl(theta)_a
   ```
   For a PERFECT equilibrium (d^2phi/dt^2 = 0), RHS_a = 0 everywhere.
   In practice, the braid breathes, so RHS_a is not exactly zero but
   oscillates. Its TIME AVERAGE should be ~zero.

3. Now PERTURB the background. Replace phi_a(x) with phi_a(x) * (1 + alpha*x)
   for a small alpha. Recompute RHS_a. The difference:
   ```
   delta_RHS_a = RHS_a(perturbed) - RHS_a(unperturbed)
   ```
   is the extra acceleration due to the gradient.

4. The net force (in the x-direction) is:
   ```
   F_x = sum_voxels sum_a delta_RHS_a(x) * phi_a(x) * dx^3
   ```
   (This uses the momentum formula: force = sum of field * acceleration * volume.)

   Wait, this is still not quite right. The force on a localized object
   is most cleanly extracted as the rate of change of its center of mass:
   ```
   F_x = d/dt(M * X_cm) = integral sum_a (d^2phi_a/dt^2) * phi_a * x * dx^3 / integral phi_a^2 dx^3
   ```

   This is getting messy. The cleanest approach, which avoids all these
   subtleties, is:

### SIMPLEST CORRECT APPROACH: Run the Simulation

Run the ACTUAL gradient test using the simulation kernel, exactly as
V33 did, but extract C from the initial acceleration instead of waiting
for drift. This uses NO analytical formula — just the simulation code.

1. **Create a gradient seed**: Use gen_braid.c (or manually construct)
   to place a braid in a background with a linear density gradient.
   Specifically:
   ```
   phi_a(x) = A_bg * (1 + alpha*x) * cos(k*z + delta_a) + braid_profile(x)
   ```
   Use alpha = 0.01 (1% per code unit).

2. **Run for a very short time** (T=2, about 250 steps). Output at
   snap_dt=0.1 (20 frames).

3. **Track the braid centroid** at each frame. Fit x_cm(t) = x0 + v*t + (1/2)*a*t^2.
   The acceleration a is the force F/M.

4. **Compute C**: The background density gradient is drho/dx = alpha * rho_bg.
   Then C = F / (drho/dx) = (M*a) / (alpha * rho_bg).

5. **Repeat at 4 gradient strengths**: alpha = {0.005, 0.01, 0.02, 0.04}.
   Verify F is proportional to alpha (linear response).

6. **Repeat for the 3-field equation** (eta=0): This isolates the
   gravitational force channel from the electromagnetic channel. The
   V33 measurement of C=186 used the 6-field equation, so it includes
   both gravity and EM. The 3-field run gives C_gravity alone.

### Input Files

- Seed generator: `/home/d/code/scp/sfa/seed/gen_braid.c` (modify to
  accept a gradient parameter, or write a thin wrapper)
- Sim kernel: `/home/d/code/scp/sfa/sim/scp_sim.c`
- Config template: `/home/d/code/scp/sfa/sim/braid_default.cfg`
- Use N=128, L=15, matching V33 parameters

### Measurements to Report

1. C (6-field, eta=0.5) from the centroid acceleration, with error bars
   from the 4-alpha linear fit
2. C (3-field, eta=0) from the same procedure — gravity-only channel
3. C_EM = C(6-field) - C(3-field) — electromagnetic contribution
4. Ratio C_EM / C_gravity — how much of the V33 force is EM?
5. Force exponent: repeat at multiple separations (if using two braids
   instead of a gradient) to measure F(D). But the gradient method
   directly gives C, which is the coefficient, not the exponent.

For the EXPONENT test, the gradient method gives C but not the radial
dependence. The exponent comes from the depletion profile rho(r). To
test the F ~ 1/D^(n+1) prediction:

6. Measure the depletion profile rho(r) around the braid from the last
   frame of braid_hires.sfa. Fit rho(r) = rho_bg - B/r^n for r = 5-12.
   Use this n to predict the force exponent as n+1.
7. Compare the predicted n+1 to the measured 1.8 from V33.

### Success Criteria

- **C measurement**: If C(numerical) agrees with C=186 (V33) within 20%,
  the force integral framework is validated. If C(numerical) differs by
  more than 50%, either the force integral is wrong or the V33
  measurement was contaminated by other effects.
- **Linearity**: F must be linear in alpha (R^2 > 0.99 across the 4
  gradient strengths). Nonlinearity would indicate the perturbation is
  too large.
- **3-field vs 6-field**: If C(3-field) << C(6-field), the V33 force
  is dominated by electromagnetic effects, and the "gravity" interpretation
  needs revision. If C(3-field) ≈ C(6-field), gravity dominates.
- **Exponent**: If n+1 from the measured depletion profile matches the
  force exponent 1.8 within 0.2, the prediction framework works. If
  they disagree by > 0.4, the simple F = -C*grad(rho) model is missing
  physics.

### Pitfalls to Avoid

1. **Do NOT use dimensional analysis to extract C** — the whole point
   is to compute it numerically with no free parameters.
2. **Short integration time**: Run only T=2-5. Longer runs let the braid
   drift significantly, and the gradient it experiences changes. We want
   the INSTANTANEOUS force, measured from the initial acceleration.
3. **Gradient too large**: alpha > 0.05 puts the braid in a highly
   nonlinear gradient where F ~ alpha is no longer valid. Use alpha < 0.04.
4. **Background normalization**: The density gradient must be applied to
   the BACKGROUND, not the braid. If you simply multiply all phi by
   (1+alpha*x), you also distort the braid's internal structure, which
   changes its mass. Apply the gradient only to the A_bg component.
   This requires separating the braid from the background, which is
   nontrivial. Alternative: initialize the background with a gradient
   SEPARATELY from the braid (place the braid in a pre-existing gradient
   background). The seed generator should handle this.
5. **PBC effects**: At N=128, L=15, the braid's Yukawa tail extends to
   the boundary. Gradient tests should use absorbing boundaries
   (damp_width=3, damp_rate=0.01) to suppress reflections. The gradient
   itself is linear across the box, which is consistent with PBC (the
   jump at the boundary is 2*alpha*L, which is absorbed by the damping
   layer).

### Expected Runtime

- 4 gradient strengths x 2 equation types (3-field, 6-field) = 8 runs
- Each run: N=128, T=5, ~625 steps at dt~0.008. At ~1 ms/step: <1 second.
- Total: under 10 seconds of simulation time
- Analysis (centroid tracking, quadratic fit): trivial
- **Total: under 1 minute including compilation and I/O**
