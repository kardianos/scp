# Fixed-Point Q16.16 Kernel: Known Issues

Status: BLOWS UP — 275% energy drift over T=460 while identical physics
in fp64 (variant_b) conserves energy to 0.5%. Root causes identified below.

## Issue 1 (CRITICAL): Energy diagnostic measures wrong Hamiltonian

The energy function sums `s_etr + s_ecf + s_etm` where s_etm includes a
theta mass term computed with `(df_dTheta + confine_deriv)`. But the FORCE
uses `(df_dTheta - confine_deriv)`. The energy diagnostic is measuring a
quantity that is NOT conserved by the integrator — it adds a term
`V_abs * gamma/Theta_sum * theta^2` that the force actively works against.

The fp64 version computes E_confine = |V_base| * gamma * ln(1+Theta/theta_c)
directly from the potential, not from the force coefficient. The fixed-point
should do the same.

Much of the "275% drift" may be a measurement artifact. The integrator might
be conserving the correct Hamiltonian while the diagnostic reports a different,
inconsistent quantity.

**Fix:** Replace s_etm computation with the actual potential energy:
  E_confine = |V_base_f64| * gamma * log(1 + Theta_f64 / theta_c)
(Already computed in fp64 for diagnostics — just use it directly.)

## Issue 2 (HIGH): P² underflows to zero at background

At phi = 0.1 (Q16.16 int = 6554):
  P = fp_mul(fp_mul(6554, 6554), 6554) = 65  (representing 0.00099)
  P² = fp_sq(65) = (65*65) >> 16 = 4225 >> 16 = 0  (EXACTLY ZERO)

In fp64: P² ≈ 9.8e-7 (nonzero), V_base ≈ -2.0e-5 per voxel.

The force sees V_base = 0 (no potential), but the fp64 energy diagnostic
sees V_base = -2e-5. Over ~2M background voxels, this is ~0.17 energy
units of "phantom potential" visible to diagnostics but invisible to forces.

When voxels cross the P² = 0 → 1 threshold, the energy jumps discontinuously
with no corresponding force work. This breaks symplecticity.

**Fix:** Compute P and P² in wider intermediate format:
  int64_t P_wide = ((int64_t)p0 * p1) >> 8;  // Q24.8
  P_wide = (P_wide * (int64_t)p2) >> 8;        // still Q24.8
  int64_t P2_wide = (P_wide * P_wide) >> 8;    // Q16.16-ish, no underflow
Then scale back to Q16.16 for the force computation.

## Issue 3 (HIGH): Arithmetic right-shift bias on negative products

fp_mul uses: `(int32_t)(((int64_t)a * b) >> FP_SHIFT)`

Right-shifting a negative int64 rounds toward -infinity (arithmetic shift).
Right-shifting a positive int64 rounds toward zero. This asymmetry means:
  +32767 >> 16 = 0   (rounds +0.4999 to 0)
  -32767 >> 16 = -1  (rounds -0.4999 to -1)

Since mu < 0, the potential force involves negative intermediate products.
These are systematically rounded MORE NEGATIVE than exact, making the
attractive force systematically ~0.5 ULP stronger per multiplication.

Over ~24M biased multiplications per step (2M voxels × 4 muls × 3 components),
this pumps energy into the system continuously.

**Fix:** Use rounding instead of truncation:
  return (int32_t)(((int64_t)a * b + (1 << (FP_SHIFT-1))) >> FP_SHIFT);
The +32768 before shifting rounds to nearest instead of toward -infinity.

## Issue 4 (MEDIUM): Small accelerations truncate to zero

HALF_DT = 129 (Q16.16 for dt/2 ≈ 0.002). An acceleration of 0.001
(Q16.16 int = 66) gives: fp_mul(129, 66) = 8514 >> 16 = 0.

Any acceleration below ~0.008 physical units produces ZERO velocity update.
Background voxels with weak curl coupling or gradient forces get no updates
at all. In fp64, v += hdt * a = 2e-6 — tiny but nonzero, accumulates
correctly over thousands of steps.

This breaks the reversibility of the Verlet integrator for weak-force regions.

**Fix:** Same rounding fix as Issue 3 helps. Also consider using Q8.24 for
HALF_DT (more fractional bits for the small multiplier) or accumulating
velocity updates in int64 before truncating.

## Issue 5 (MEDIUM): dVdP chain precision loss

dVdP = fp_div(fp_mul(MU, P), den2) chains two operations with truncation
at each step. At intermediate phi (0.05-0.08), the error is 1.6-2.3%
compared to fp64. This means the integrator solves a slightly different
ODE than the energy function measures.

**Fix:** Compute the numerator in int64: MU * P in int64, then divide by
den2 in int64, truncate to int32 only at the end.

## Issue 6 (LOW): fp_ln worst-case ~3% error

The 7-term Maclaurin series for ln(1+u) has ~3% error at u→1 (near
power-of-two boundaries). Affects confinement force strength.

**Fix:** Use a minimax polynomial or more terms. See TODO in fp_ln().

## Issue 7 (LOW): DT/HALF_DT rounding consistency

DT and HALF_DT are rounded independently. If DT is odd, 2*HALF_DT != DT,
making the Verlet integrator formally non-symplectic. Currently a non-issue
(values match) but a latent bug for other grid sizes.

**Fix:** Compute HALF_DT = DT >> 1 (derive from DT, don't round independently).

## Overall Assessment

Q16.16 is TOO NARROW for this physics. The dynamic range from background
(P² ~ 10⁻⁶) to braid core (P² ~ 10⁻³) spans 3 decades, and 16 fractional
bits can't represent both ends. The minimum viable format is probably:
  - Q24.8 for field values (more integer range)
  - Q16.16 for velocities and forces
  - int64 for intermediate products (P, P², V_base)

Or simply Q32.32 (int64 everywhere), which gives 32 bits of fraction —
equivalent to fp32 precision with exact integer arithmetic. Memory cost
is the same as the current fp64 simulation.
