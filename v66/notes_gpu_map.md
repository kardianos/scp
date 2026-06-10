# GPU kernel map (from foundations workflow, 2026-06-09)

Perfect! Now I have comprehensive information. Let me compile my findings for both parts:

## PART A: GPU KERNEL STRUCTURE SURVEY

**KERNEL FUNCTION INVENTORY** (lines 553-1139):
1. **compute_intermediates_kernel** (line 553): Computes Cosserat mismatch M and hardening Q vectors; only runs if alpha_cs>0 or beta_h>0
2. **compute_forces_kernel** (line 594): Main physics—Laplacian, masses, potential forces, curl couplings, theta saturation, chiral interactions, Cosserat strain, hardening, cross-couplings
3. **verlet_halfkick_kernel** (line 838): Velocity half-kicks in Verlet integrator
4. **verlet_drift_kernel** (line 851): Position drift in Verlet integrator
5. **absorbing_boundary_kernel** (line 864): Sponge layer damping at r>R-damp_width
6. **gradient_bc_kernel** (line 892): bc_type=1 analytical gradient boundary (no pinned arrays on GPU)
7. **downcast_f64_to_f32_kernel** (line 948): GPU-side double→float conversion for precision=1 snapshots
8. **downcast_f64_to_f16_kernel** (line 954): GPU-side double→f16 conversion for precision=0 snapshots
9. **reduce_diagnostics_kernel** (line 976): GPU reduction to compute 12 energy/diagnostics scalars (epk, etk, eg, em, ep, etg, etm, ec, phi_max, P_max, P_int, theta_rms_sum)
10. **zero_diag_kernel** (line 1127): Zero diag results buffer
11. **rescale_phi_kernel** (line 1134): Scale phi+vel by factor f for self-tuning charge projection

**FIELD ARRAY ALLOCATION** (lines 1145-1171):
- **Separate per-field pointers**, NOT one big buffer:
  - d_phi[3], d_vel_phi[3], d_acc_phi[3]: 3 phi components × 3 array types
  - d_theta[3], d_vel_theta[3], d_acc_theta[3]: 3 theta components × 3 array types
  - d_mismatch[3], d_harden_Q[3]: 6 intermediates (allocated conditionally)
- Each array is N³ doubles (8 bytes)
- Total: 18 arrays (24 with intermediates)
- gpu_blocks = ceil(N³/THREADS_PER_BLOCK) where THREADS_PER_BLOCK=256
- Memory report at line 1169: "GPU: allocated X.XX GB (physics±Cosserat)"

**ASYNC I/O SNAPSHOT ARCHITECTURE** (lines 221-244, 1612-1737):
- **SnapHookCtx structure** (line 221):
  - GPU staging buffers: d_f16_buf (f16 mode) or d_f32_buf (f32 mode), or none (f64)
  - Pinned host buffer: h_pin_buf (always, size = 12×N³×sizeof_precision)
  - Async writer thread + pthread condition variables
  - Separate cudaStream for DMA (non-blocking)
- **snap_hook function** (line 1680): Called every snap_every steps
  - Converts 12 fields (phi[3], theta[3], vel_phi[3], vel_theta[3]) to precision format on GPU (lines 1700-1718)
  - DMA to pinned host (async via stream)
  - Signals writer thread
  - Returns immediately; physics continues on default stream
- **Writer thread** (line 1572): Compresses + writes to SFA file (zstd BSS encoding, parallel per-column, lock only for file I/O)
- **Memory overhead**: 
  - f16: GPU buf + host buf = 2×12×N³×2 bytes
  - f32: GPU buf + host buf = 2×12×N³×4 bytes
  - f64: host buf only = 12×N³×8 bytes

**DIAGNOSTICS REDUCTION** (lines 976-1124, 1743-1801):
- **GPU kernel reduce_diagnostics_kernel**: Block-level reduction (shared memory) + atomic global reduction
- Computes 12 values: DIAG_NVALS=12 (line 974)
  - Indices [0-7]: sums (epk, etk, eg, em, ep, etg, etm, ec)
  - Indices [8-9]: max reductions (phi_max, P_max)
  - Index [10]: P_int (integral of |P| over volume)
  - Index [11]: theta_rms_sum (for RMS calculation)
- Results stored in pinned d_results (12 doubles) on GPU, copied to h_results on host
- **DiagHookCtx** (line 246): File pointer, config parameters, GPU/host result buffers
- **diag_hook** (line 1803): Runs reduction, prints to TSV file + stdout, tracks energy drift

**NFIELDS HARDCODING** (line 15 scp_config.h):
- `#define NFIELDS 3` — appears in:
  - Grid allocation (line 275): phi[NFIELDS], theta[NFIELDS], etc.
  - GPU upload/download loops (lines 1175-1192)
  - Diagnostics energy calc (line 1856)
  - Loop bounds in many init functions
  - **NOT hardcoded in GPU kernels** (they reference specific indices [0],[1],[2] explicitly for phi; theta same)

**SELF-TUNE INTEGRATION** (lines 3102-3176):
- Hooks in at lines 3109-3114 (setup) and 3139-3176 (per-step)
- Reads P_max and E_mass from GPU diagnostics
- Dynamically adjusts kappa via cudaMemcpyToSymbol (line 3164)
- Optionally rescales phi to maintain fixed charge Q via rescale_phi_kernel (line 3148)

**12-FIELD COMPLEXIFIED PORT ANALYSIS (complex_phi=1)**:
- **Field doubling**: φ_a → (u_a, v_a), θ_a → (tu_a, tv_a) = 12 real components
- **Kernels touched**:
  - compute_intermediates: Must handle 6 curl pairs; pairwise eta coupling
  - compute_forces: Critical—all 12 forces, s=|u|²+|v|² integral, Vt'(s) derivative chaining
  - Verlet kernels: Scale trivially (just indices 0-11 instead of 0-5)
  - Diagnostics: Energy sums expand; P → u₀u₁u₂ (product of real parts only, OR real(U₀U₁*U₂*) depending on exact definition)
  - gradient_bc, absorbing_bc: Expand to 12 field components
  - downcast kernels: Still output 12 fields to SFA (no increase in I/O footprint!)
  - rescale: Scale 12 field indices

- **Memory impact at N=384**:
  - Current 6-field: 18 arrays + 6 intermediates = 24 × 384³ × 8 = **10.87 GB physics** + 0.375-5.44 GB snapshot overhead = 11.2-16.3 GB total
  - 12-field: 36 arrays + 6 intermediates = 42 × 384³ × 8 = **19.03 GB physics** + 0.375-5.44 GB snapshot = 19.4-24.5 GB total
  - **Exceeds V100-16GB by 3.4 GB with f32 snapshots**

- **Maximum N on constrained memory**:
  - 16 GB (V100) with f32 snapshots: max N ≈ **360** (16.05 GB total)
  - 32 GB with f32 snapshots: max N ≈ **455** (32.02 GB)
  - 40 GB (A100) with f16 snapshots: max N ≈ **486** (39.93 GB)

---

## PART B: RADIAL SOLVER REUSE ASSESSMENT

**ODE INTEGRATION METHOD** (radial_oscillon.c lines 40-78):
- **Leapfrog (velocity Verlet)**, identical to GPU integrator:
  - `vel[a][i] += 0.5*dt*acc[a][i]` (half-kick)
  - `phi[a][i] += dt*vel[a][i]` (drift)
  - Force recompute (COMPUTE_ACC macro, lines 53-62)
  - `vel[a][i] += 0.5*dt*acc[a][i]` (half-kick)
- **Boundary handling**: Sponge damping at r > R-damp_width (line 60)
- **r=0 regularity**: Laplacian at i=0 uses `lap = 6(φ[1]-φ[0])/dr²` (line 57)
- **Time steps**: dt = dt_factor × dr (line 40), dt_factor=0.2 default

**CLI/CONFIG INTERFACE** (lines 25-38):
```
-N 1000 -R 30.0 -A 0.7 -sigma 1.3 -omega 1.3 -T 200 
-kappa 50 -mu -41.345 -m2 2.25 
-cool 0.0 -damp 0.05 -damp_width 6.0 
-phase 0/1 (0=from rest, 1=on limit cycle with velocity) 
-delta <d0>,<d1>,<d2> -diag_dt 1.0
```

**OUTPUT FORMAT** (lines 50, 65-72):
- stderr: header + parameters
- stdout: `t  P_max  E_core(r<5)` (header line 65, diagnostics every diag_every steps)
- Final line (82-84): summary `RESULT A σ ω phase cool : Pmax0 Pmax_final ratio lifetime [COHERENT|decayed]`

**DIRECT REUSABILITY FOR Q-BALL SHOOTING SOLVER**:
✓ **YES — very high reuse**:
1. **ODE loop structure** (lines 67-78): Copy verbatim, replace COMPUTE_ACC with complex Q-ball forces
2. **Leapfrog integrator**: Unchanged
3. **Damping/sponge**: Reuse exact code (line 60)
4. **Diagnostics pattern**: Integrate energy E and charge Q during evolution
5. **Output file format**: Extend header to include Q, output Q alongside E_core
6. **Shooting setup**: Seed φ_a(0) via amplitude parameter, search over {ω, A} → converge to threshold where oscillations stabilize

**What MUST change**:
- Force calculation: Replace `pf[a] = mu*phi_a*(prod_rest)/(1+kappa*P²)²` with:
  - `s = u₀²+v₀² + u₁²+v₁² + u₂²+v₂²`
  - `Vt'(s) = mu/(1+kappa*s)²`
  - `pf_u[a] = 2*Vt'(s)*u_a`, `pf_v[a] = 2*Vt'(s)*v_a`
  - Plus curl terms `eta*curl(v)` on u, `eta*curl(u)` on v
- Energy integral: Add `∫(u·curl(v) - v·curl(u))` term = measure Noether charge Q
- Shooting boundary condition: Match φ(r→∞) ~ 0 (exponential tail e^{-r/R} with R=1/√(m²-ω²))

**radial_cascade.c** (lines 1-88):
- Tests φ + density ρ back-reaction; not directly reusable for Q-ball (no density sector in Q-ball)
- Thomas-Fermi density (lines 74-80) is standalone; Q-ball solver doesn't need it
- Shows pattern for **conserved charge constraint** (line 27: Qd0=100.0 fixed): interpolation to maintain ∫ρ = Qd0
- **Reusable concept**: Use identical scheme to enforce ∫|φ²| or ∫(u²+v²) = Q_target during radial evolution

---

## stability_estimates.py ANALYSIS

**What it computes**:
1. **Derrick saddle confirmation** (lines 30-35): Static is a maximum in the scaling direction → must use oscillating breather (oscillon) or U(1) charge (Q-ball)
2. **Oscillon R(ω) relation** (line 42): `R = 1/√(m² - ω²)` — the exponential decay length; matches linear dispersion
3. **Amplitude-frequency balance** (lines 49-56):
   - Core force: `|μ|·g5·A⁴ / sat` where g5=cos²(δ₀)cos²(δ₁)cos²(δ₂), sat=(1+κP²)²
   - Self-consistent: `m² - ω² = |μ|·g5·A⁴ / sat` → `A = (gap/(|μ|g5))^(1/4)`
   - Returns `(A, P_core, κP²)` for each ω
4. **Theta radiation constraint** (lines 70-80): Massless theta driven by phi oscillating at ω radiates as k=ω plane waves (propagating away). Prefers **lower ω** to reduce theta drive, trading width R for longevity.
5. **Seed recommendations** (lines 84-98):
   - Tests omega ∈ {0.9, 1.1, 1.25, 1.35, 1.42, 1.46} and corresponding (σ, A)
   - Checks: A decreases toward m (correct). Previous seed (A=1.0, σ=1.5) was over-driven.
   - Recommends omega ∈ [1.0, 1.4], prints explicit -sigma -A values

**What it does NOT compute**:
- ✗ Coleman thin/thick-wall window (not applicable here; potential is nonpolynomial)
- ✗ Exact stability window (only estimates order-of-magnitude via unsaturated balance)
- ✗ Theta phase-matching analysis (heuristic: lower ω preferred, but no quantitative window)

**Reusability for Q-ball existence window**:
- **Partially reusable**, but must rederive for U(1) Q-ball:
  - Q-ball potential at finite rotation: `U(|φ|², ω) = V(|φ|²) - (ω/2)|φ|²` (rotating frame)
  - Existence window: `ω ∈ [0, m)` (all ω below mass gap allow Q-balls, unlike osci where ω ∈ (m-ε, m))
  - Amplitude from core balance: `2Vt'(s)·s = (ω² - m²)·s + ...` (different scaling)
  - Real window computed by **radial shooting** (not this script)

**SCRIPT SHOULD EMIT**:
```python
# Q-ball variant estimates
QM2_min, QM2_max = 0.0, M2  # omega ∈ [0, m)
# For each omega, compute amplitude such that shooting converges to Q-ball soliton
# (not yet in this script — would require eigenvalue solver)
```

---

## FINAL SUMMARY

**PART A — GPU Port Readiness**:
- Kernel structure is **highly modular** — 11 kernels, 6 independent force components
- Field arrays use **per-component pointers** (easy to extend from 3→6 components)
- Async I/O pipeline **decouples snapshots from physics** (stream-based DMA, separate writer thread)
- Memory doubling (6→12 fields) **reduces max N from 384→360 on 16GB**, requires f16 snapshots or larger GPU
- **Critical hotspot**: compute_forces_kernel (line 594) — all 12 field indices must be explicit; no macro loop
- Self-tune integration is **parameter-only** (cudaMemcpyToSymbol), decoupled from kernels

**PART B — Radial Solver Reuse**:
- **Leapfrog integrator is plug-and-play** — identical Verlet loop
- **Force computation changes only** — replace quintic cubic source with U(1)-invariant quadratic source
- **Output format trivially extends** — add Q column, keep diag_dt/T/R/dt_factor interface
- **stability_estimates.py is NOT existence window** — only order-of-magnitude oscillon seed estimates; Q-ball window requires radial shooting (new code)
- **Conserved-charge enforcement** — can reuse radial_cascade.c's interpolation pattern (solve for Lagrange multiplier μ to hit Q_target)

File paths and exact line numbers as requested above.