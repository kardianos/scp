# V49: Fixed-Point Leapfrog Simulation

## Summary

Same PDE. Same leapfrog integrator. Same stencil. Same physics.
Replace fp64 with fixed-point integer arithmetic.

This is NOT a new physics model. It's the existing 6-field Cosserat
simulation rewritten in integer math. Every force term, every energy
contribution, every update step has a direct 1:1 correspondence with
the current fp64 code.

## Why

1. **Exact reproducibility** — bit-identical results on any platform.
   No IEEE 754 rounding mode differences, no FMA vs separate mul+add.
2. **No float drift** — integer addition is exact. Energy conservation
   limited only by the integrator (leapfrog O(dt²)), not by accumulation
   of rounding errors across millions of steps.
3. **4x memory reduction** — int16 fields vs fp64: 24 bytes/voxel vs 96.
   At N=384: 1.3 GB vs 5.2 GB. Fits larger grids in GPU memory.
4. **SIMD throughput** — 16 × int16 ops per 256-bit register vs 4 × fp64.
   GPU int16 (dp4a) is heavily optimized for ML workloads.
5. **Compatible with unified transfer** — Variant B v2 potential works
   identically, just with integer scale factors.

## Representation

### Field values: int16 (16-bit signed integer)

    A_UNIT = 1.0 / 256.0
    phi_int = round(phi_float / A_UNIT)
    phi_float = phi_int * A_UNIT

    Range: -32768 to +32767 → -128.0 to +127.996 code units
    Resolution: 1/256 = 0.00391 code units

Typical values in the simulation:
    A_bg = 0.1      → phi_int ≈ 26
    Braid peak ≈ 0.3 → phi_int ≈ 77
    phi_max ≈ 0.9    → phi_int ≈ 230

All well within int16 range with ~8 bits of precision above the
background. Sufficient for the physics (the fp32 output mode already
truncates to ~7 decimal digits).

### Velocities: int16

    V_UNIT = A_UNIT / DT_UNIT
    phi_v_int = round(phi_v_float / V_UNIT)

Typical velocity ~0.5 code units/time → phi_v_int ≈ 128.

### Energies: int64 (64-bit signed integer)

Energy computations involve products of field values (int16 × int16 =
int32) accumulated across stencils. Use int64 for the accumulator to
prevent overflow. The energy scale factor E_UNIT absorbs all the dx²,
dt², A_UNIT² conversions.

### Scale factors: int32 (pre-computed constants)

Each physics parameter becomes an integer constant:

    M2_FIXED   = round(m² × A_UNIT² × E_SCALE / (2 × E_UNIT))
    MU_FIXED   = round(|mu| × A_UNIT⁶ × E_SCALE / (2 × E_UNIT))
    KAPPA_FIXED = round(kappa × A_UNIT⁶ × E_SCALE)
    ETA_FIXED  = round(eta × A_UNIT² × E_SCALE / (2 × dx × E_UNIT))
    GRAD_FIXED = round(A_UNIT² × E_SCALE / (2 × dx² × E_UNIT))

These are computed once at startup from the fp64 config parameters.
The simulation never uses floating point after initialization.

## The Integrator

Standard velocity-Verlet (leapfrog), translated to fixed-point.
Identical structure to the current scp_sim.c verlet_step().

```c
// ================================================================
// Fixed-point Verlet step — direct translation of fp64 version
// ================================================================

void verlet_step_fixed(GridFixed *g) {
    // Half-kick: v += (dt/2) * a
    for (int a = 0; a < 3; a++) {
        for (long idx = 0; idx < N3; idx++) {
            // int16 velocity += (int16 accel >> 1)
            // The >> 1 is the dt/2 factor (absorbed into scale)
            g->phi_v[a][idx] += g->phi_a[a][idx] >> HALFKICK_SHIFT;
            g->theta_v[a][idx] += g->theta_a[a][idx] >> HALFKICK_SHIFT;
        }
    }

    // Drift: x += dt * v
    for (int a = 0; a < 3; a++) {
        for (long idx = 0; idx < N3; idx++) {
            g->phi[a][idx] += g->phi_v[a][idx] >> DRIFT_SHIFT;
            g->theta[a][idx] += g->theta_v[a][idx] >> DRIFT_SHIFT;
        }
    }

    // Forces
    compute_forces_fixed(g);

    // Half-kick again
    for (int a = 0; a < 3; a++) {
        for (long idx = 0; idx < N3; idx++) {
            g->phi_v[a][idx] += g->phi_a[a][idx] >> HALFKICK_SHIFT;
            g->theta_v[a][idx] += g->theta_a[a][idx] >> HALFKICK_SHIFT;
        }
    }
}
```

The HALFKICK_SHIFT and DRIFT_SHIFT constants encode dt/(2*V_UNIT*A_UNIT)
and dt*V_UNIT/A_UNIT respectively. They're chosen so that the bit-shift
gives the correct scaled result without multiplication.

## Force Computation

Direct translation of the fp64 compute_forces(). Every term maps 1:1.

```c
void compute_forces_fixed(GridFixed *g) {
    for (long idx = 0; idx < N3; idx++) {
        int16_t p0 = g->phi[0][idx];
        int16_t p1 = g->phi[1][idx];
        int16_t p2 = g->phi[2][idx];

        // Laplacian: 7-point stencil, integer arithmetic
        // lap = (sum of 6 neighbors - 6*center) / dx²
        // In fixed-point: pre-multiply by GRAD_SCALE
        for (int a = 0; a < 3; a++) {
            int32_t lap = (int32_t)g->phi[a][ip] + g->phi[a][im]
                        + g->phi[a][jp] + g->phi[a][jm]
                        + g->phi[a][kp] + g->phi[a][km]
                        - 6 * (int32_t)g->phi[a][idx];
            // lap is in int32, scale to acceleration units
            int16_t lap_scaled = (int16_t)(lap >> LAP_SHIFT);

            // Mass term: -m² × phi
            int16_t mass_term = (int16_t)(-(M2_SCALED * (int32_t)g->phi[a][idx]) >> MASS_SHIFT);

            // Triple product potential: -dV/dP × dP/dphi_a
            int32_t P = (int32_t)p0 * p1;  // int32
            P = (P * (int32_t)p2) >> P_SHIFT;  // scaled triple product
            int32_t P2 = ((int64_t)P * P) >> P2_SHIFT;
            int32_t den = DENOM_ONE + ((int64_t)KAPPA_SCALED * P2 >> KAPPA_SHIFT);
            int32_t dVdP = (int32_t)((int64_t)MU_SCALED * P / ((int64_t)den * den >> DEN_SHIFT));

            int32_t dPda;
            if (a == 0) dPda = ((int32_t)p1 * p2) >> P_SHIFT;
            else if (a == 1) dPda = ((int32_t)p0 * p2) >> P_SHIFT;
            else dPda = ((int32_t)p0 * p1) >> P_SHIFT;

            int16_t pot_term = (int16_t)(((int64_t)dVdP * dPda) >> POT_SHIFT);

            // Curl coupling: eta × curl(theta)_a
            int32_t ct = curl_fixed(g->theta, a, ip,im,jp,jm,kp,km);
            int16_t curl_term = (int16_t)((ETA_SCALED * ct) >> CURL_SHIFT);

            // Total acceleration
            g->phi_a[a][idx] = lap_scaled + mass_term + pot_term + curl_term;
        }

        // Theta forces — same structure, with transfer potential
        int16_t t0 = g->theta[0][idx];
        int16_t t1 = g->theta[1][idx];
        int16_t t2 = g->theta[2][idx];
        int32_t Theta = (int32_t)t0*t0 + (int32_t)t1*t1 + (int32_t)t2*t2;
        Theta >>= THETA_SHIFT;

        // V_abs = |V_base(P)| (positive, for theta mass)
        int32_t V_abs = abs(((int64_t)MU_SCALED * P2) / (2 * den));

        // Transfer: f = eps + (1-eps) * Theta / (Theta + Theta_c)
        int32_t Theta_sum = Theta + THETA_C_SCALED;
        int32_t f_transfer = EPSILON_SCALED
            + ((int64_t)(ONE_MINUS_EPS) * Theta / Theta_sum);

        // Confinement: gamma / (Theta + Theta_c)
        int32_t confine = (int64_t)GAMMA_SCALED * ONE_SCALED / Theta_sum;

        // Theta effective mass: |V_abs| × (df/dTheta + confine) × 2
        int32_t df_dTheta = (int64_t)(ONE_MINUS_EPS) * THETA_C_SCALED
                          / ((int64_t)Theta_sum * Theta_sum >> TSUM_SHIFT);
        int32_t theta_mass = ((int64_t)V_abs * (df_dTheta + confine)) >> TMASS_SHIFT;

        for (int a = 0; a < 3; a++) {
            int32_t lapt = (int32_t)g->theta[a][ip] + g->theta[a][im]
                         + g->theta[a][jp] + g->theta[a][jm]
                         + g->theta[a][kp] + g->theta[a][km]
                         - 6 * (int32_t)g->theta[a][idx];
            int16_t lapt_scaled = (int16_t)(lapt >> LAP_SHIFT);

            int32_t cp = curl_fixed(g->phi, a, ip,im,jp,jm,kp,km);
            int16_t curl_p = (int16_t)((ETA_SCALED * cp) >> CURL_SHIFT);

            // Theta mass (confining, positive)
            int16_t mass_t = (int16_t)(-(theta_mass * (int32_t)g->theta[a][idx]) >> TMASS2_SHIFT);

            g->theta_a[a][idx] = lapt_scaled + mass_t + curl_p;
        }

        // Phi transfer modulation: dVdP × f_total
        // (Applied retroactively — multiply pot_term by f_total/f_base)
        // This is where the unified transfer affects phi.
        // In practice, compute dVdP*f_total directly above instead of
        // computing dVdP then multiplying. Left separate here for clarity.
    }
}
```

## Shift Constants

The _SHIFT values are determined at startup by the scale factor
analysis. They encode the combined effect of A_UNIT, dx, dt, and the
physics constants. Example derivation:

    Physical: lap_physical = (sum - 6*center) / dx²
    Fixed:    lap_int = sum_int16 - 6*center_int16  (range: ±6×32767 = ±196602, fits int32)
    Scaling:  lap_physical = lap_int × A_UNIT / dx²
    For acceleration: a = lap_physical (in phi equation)
    In velocity units: a_int = lap_int × A_UNIT / (dx² × V_UNIT)
    As bit shift: LAP_SHIFT = log2(dx² × V_UNIT / A_UNIT)

The startup code computes all _SHIFT values from the fp64 config and
verifies they're in valid range (0-31 for int32 shifts).

## Memory Layout

```c
typedef struct {
    int16_t *phi[3];      // 3 × N³ × 2 bytes
    int16_t *theta[3];    // 3 × N³ × 2 bytes
    int16_t *phi_v[3];    // 3 × N³ × 2 bytes
    int16_t *theta_v[3];  // 3 × N³ × 2 bytes
    int16_t *phi_a[3];    // 3 × N³ × 2 bytes  (acceleration)
    int16_t *theta_a[3];  // 3 × N³ × 2 bytes  (acceleration)
    int N;
    long N3;
    // Scale constants (int32)
    int32_t M2_SCALED, MU_SCALED, KAPPA_SCALED, ETA_SCALED;
    int32_t EPSILON_SCALED, THETA_C_SCALED, GAMMA_SCALED;
    int HALFKICK_SHIFT, DRIFT_SHIFT, LAP_SHIFT; // etc
} GridFixed;
```

Memory per voxel: 18 arrays × 2 bytes = 36 bytes
(vs current: 18 arrays × 8 bytes = 144 bytes)

At N=384: 36 × 384³ = 2.0 GB (vs 8.2 GB)
At N=512: 36 × 512³ = 4.8 GB (fits V100-16GB, current doesn't)

## Precision Analysis

With int16 (range ±32767):
- Field values: ~8 bits above background (256 × 0.1 = 26 background,
  256 × 0.9 = 230 peak, ratio 230/26 ≈ 9, ~3 bits dynamic range
  above background). This is TIGHT.

With int32 for intermediate products:
- P = phi[0] × phi[1] × phi[2]: up to 77³ = 456533, fits int32
- P² can overflow int32 (456533² > 2³¹). Use int64 for P².
- Energy sums over N³ voxels: up to ~10⁷ × 10⁹ = 10¹⁶ — fits int64.

**Concern: int16 may not have enough precision for the field values.**

The background amplitude is 26 int units, and the braid perturbation
is ~50 int units above that. The DIFFERENCE between adjacent voxels
in the gradient might be only 1-3 int units. This means the Laplacian
has only ~2 bits of precision.

**Mitigation: use int32 for field values, int16 for velocities.**

    Field:    int32 (range ±2B, A_UNIT = 1/65536, 16 fractional bits)
    Velocity: int16 (range ±32K)
    Memory:   12 × 4 + 6 × 2 = 60 bytes/voxel (vs 144 for fp64)

This gives 16 bits of fractional precision for field values — more
than enough for the Laplacian. The velocity can stay at int16 since
it's only used for the drift step (added to position via shift).

## Alternative: int32 fields throughout

If memory isn't the binding constraint (GPU has enough), use int32
for everything:

    Memory: 18 × 4 = 72 bytes/voxel (vs 144 for fp64)
    At N=384: 72 × 384³ = 4.1 GB (fits V100-16GB)
    At N=512: 72 × 512³ = 9.7 GB (fits V100-16GB)

Still 2x memory savings, with 24+ bits of precision for all values.
No overflow concerns for any intermediate computation.

This is the SAFE option: enough precision to match fp32 results,
enough range for all products, and still 2x better than fp64.

## Unified Transfer Integration

The Variant B v2 potential works directly in fixed-point:

```c
// V_base = (mu/2) × P² / (1 + kappa × P²)    — int64 computation
int64_t V_base_64 = ((int64_t)MU_HALF * P2_64) / (DENOM_ONE_64 + KAPPA_64 * P2_64);
int32_t V_abs = (int32_t)abs(V_base_64 >> VBASE_SHIFT);

// f_transfer = eps + (1-eps) × Theta / (Theta + Theta_c)
int32_t f = EPS_32 + (int32_t)(((int64_t)ONE_MINUS_EPS_32 * Theta_32) / Theta_sum_32);

// Phi force: dVdP × dPda × f_total (f_total includes ln confinement)
// Theta force: V_abs × (df/dTheta + gamma/Theta_sum) × 2 × theta[a]
// Signs: V_abs is positive, theta force is confining. Verified.
```

All the sign issues from the fp64 proposal are moot — V_abs = abs(...)
is explicitly positive. No sign errors possible.

## Implementation Plan

Phase 1: Write scp_sim_fixed.c — int32 version of scp_sim.c.
    - Same config file, same SFA output.
    - Startup: read fp64 config, compute int32 scale constants.
    - Output: convert int32 back to float for SFA writing.
    - Test: single proton, compare E_total conservation vs fp64.

Phase 2: Verify wave propagation.
    - Initialize a plane wave, measure propagation speed.
    - Compare dispersion relation to fp64 version.
    - Verify that the shift constants give correct physics.

Phase 3: Run the two-proton binding test.
    - Same setup as V48 GPU run (D=14, N=384).
    - With unified transfer (Variant B v2).
    - Compare cluster evolution and energy partition.

Phase 4: GPU kernel (scp_sim_fixed.cu).
    - Map int32 arithmetic to CUDA.
    - Profile: is the 2x memory savings worth the int32 overhead?
    - Compare wall-clock time vs fp64 version.

## What Changes vs Current fp64

| Aspect | fp64 (current) | int32 (proposed) |
|--------|---------------|-----------------|
| Field precision | 52-bit mantissa | 32-bit (24+ usable) |
| Energy drift | O(dt²) + float accumulation | O(dt²) only (no float drift) |
| Memory/voxel | 144 bytes | 72 bytes |
| N=384 memory | 8.2 GB | 4.1 GB |
| N=512 memory | 19.4 GB | 9.7 GB |
| Reproducibility | Platform-dependent | Bit-identical everywhere |
| Overflow risk | None (fp64 range huge) | Must verify intermediate products |
| Physics | Identical | Identical (same PDE, same integrator) |
| Sign errors | Possible (mu < 0 etc) | Possible (same math, same risks) |
| Integrator | Verlet (symplectic) | Verlet (symplectic) |
| Output format | SFA (f16/f32/f64) | SFA (convert int32 → f32 on write) |
