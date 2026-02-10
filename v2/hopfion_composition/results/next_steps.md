# Next Steps — Hopfion Composition Simulator

## Completed (Phase 15B)

- [x] Analytical Skyrme force (L₄) — ported from field.c, 3-pass algorithm
- [x] Consistent 9-point Laplacian — exact discrete gradient of E₂
- [x] Numerical L₆ force — finite-difference of (B⁰)² energy
- [x] Pion mass + degenerate mass forces (V_π, V_D)
- [x] Gradient verification — 6 tests, 144/144 checks pass
- [x] Skyrmion equilibrium — arrested gradient flow, E/E_FB=1.219, E₂/E₄=0.999
- [x] Effective metric — P/m = 2.000000 confirmed (sigma model, no gravity)

## Next: Finite-λ Effective Metric (highest priority)

The P/m = 2 identity proves the sigma model has no emergent gravity. The most
immediate question is: **does finite-λ (field density reduction) produce a
Schwarzschild-like metric?**

### Step 1: Load finite-λ profile into 3D grid
- Use `profile_finlam_e1_8000.dat` from hopfion_search (λ=8000, ρ(0)=0.803)
- `init_skyrmion` already supports ρ(r) profiles
- Need to NOT sigma-project (ρ(r) ≠ ρ₀)

### Step 2: Compute effective metric at finite-λ
- c₄(r) = 2ρ(r)²/e² now varies with r
- P(r)/m(r) ≠ 2 — measure the deviation
- Expected: P/m < 2 near core (where ρ < ρ₀), P/m → 2 at large r

### Step 3: Compare to Schwarzschild
- Extract g₀₀(r) and g_rr(r) from P(r), m(r)
- Fit to Schwarzschild: g₀₀ = -(1 - r_s/r), g_rr = (1 - r_s/r)⁻¹
- Extract effective r_s and compare to soliton mass M via r_s = 2GM/c²
- Key question: does the effective G match a universal constant?

### Step 4: λ-dependence
- Scan λ from 10⁸ (sigma limit) down to 8000 (maximum deviation)
- Plot P/m vs r for several λ values
- Plot effective r_s vs λ

## Next: Hopfion Topology (Phase 3)

### Hopf charge computation
- Implement CP¹ projection: q → n̂ = f/|f| (unit vector from quaternion imaginary part)
- Compute area form pullback F_ij = n̂ · (∂_i n̂ × ∂_j n̂)
- FFT-based Hopf integral: solve ∇²A = ∇×(n̂×∇n̂), then H = ∫A·(∇×A)/(4π²)
- Test on standard Hopf map (expected H=1)

### Linking number for composed states
- Initialize two hopfions at different positions
- Compute individual and total Hopf charge
- Verify H_total = H₁ + H₂ + 2×Linking(1,2)

## Next: L₆ Equilibrium (BPS approach)

### Near-BPS soliton
- Initialize B=1, set λ₆ > 0 (small) with V (Mexican hat)
- Run gradient flow with L₂ + L₄ + L₆ + V
- Expected: E → E_FB as λ₆ dominates
- CFL constraint: dt ~ h³ for L₆ (need small dt or arrested flow)

### Multi-soliton binding near BPS
- At large λ₆: multi-Skyrmion binding energy changes dramatically
- B=2 near BPS: compare to standard B=2 bound state
- Key question: does near-BPS change the multi-soliton geometry?

## Next: Dynamic Evolution on Spherical Grid

### Port scattering capability
- Leapfrog integrator already implemented (hopfion_evolve.c)
- Two-soliton initialization already implemented (init_two_skyrmions)
- Need: sponge boundary absorbing layer (already in grid structure)
- Test: B+B scattering, compare to scatter.c results

### Hopfion dynamics
- Initialize H=1 hopfion, evolve
- Monitor both B (Skyrmion charge) and H (Hopf charge) during evolution
- Key question: is H conserved under Skyrme dynamics?

## Longer Term: Atom Building

### Proton analog
- B=1 Skyrmion at finite-λ with emergent metric
- Compute effective potential for test perturbations
- Compare bound state spectrum to hydrogen-like levels

### Electron analog
- Need a light, stable topological object distinct from the Skyrmion
- Candidates: hopfion (H=1, B=0), or degenerate sector excitation
- Key challenge: mass hierarchy M_e/M_p ≈ 1/1836

### Hydrogen atom
- Two-body system: B=1 Skyrmion + light topological object
- Bound by emergent gravity (if finite-λ produces it)
- Spectrum: compare to 1/n² scaling
