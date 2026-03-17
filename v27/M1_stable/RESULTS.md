# V27 Milestone 1: Stable Propagating Braid — RESULTS

## Goal
Achieve BOTH |P| > 0.1 AND fc > 0.5 at t=500.

## Outcome: ALL TESTS FAIL

No configuration achieved the joint target. The best single-metric results:
- **Best |P|**: 0.574 (DynA control, lpw=0.0) — but fc=0.28
- **Best fc**: 0.47 (M1b L=10 at t=400) — but dropped to fc=0.10 by t=500

The pairwise coupling is **destructive**: it kills |P| while failing to
improve fc. The smaller domain (L=10) shows transient fc improvement but
does not sustain it.

## M1a: Pairwise Coupling Scan (L=20, k=pi/20)

| lambda_pw | fc(500) | |P|(500) | E(500) | Pz(500) | PASS? |
|-----------|---------|----------|--------|---------|-------|
| 0.0 (ctrl)| 0.276   | 0.574    | 281.3  | 35.8    | NO    |
| 0.3       | 0.112   | 0.002    | 179.9  | 23.7    | NO    |
| 0.5       | 0.185   | 0.000    | 135.7  | 17.0    | NO    |
| 0.7       | 0.225   | 0.001    | 97.2   | 20.6    | NO    |

**Finding**: Pairwise coupling catastrophically destroys the triple-product
braid structure |P|. Even lambda_pw=0.3 drops |P| from 0.574 to 0.002.
The mechanism: pairwise coupling symmetrizes the three fields (pulls them
toward phi_1 = phi_2 = phi_3), which zeros the product P = phi_1*phi_2*phi_3
when all three are similar (the cos phases overlap instead of maintaining
their 2pi/3 offset). This is the opposite of what the braid needs.

The pairwise coupling also drains energy faster (E drops from 281 to 98-180),
meaning more radiation loss to the absorbing boundaries.

## M1b: Small Domain L=10 (N=128, dx=0.157)

| Config   | fc(500) | |P|(500) | E(500) | Pz(500) | PASS? |
|----------|---------|----------|--------|---------|-------|
| L=10     | 0.102   | 0.034    | 107.9  | 22.0    | NO    |

**Finding**: The smaller domain with finer resolution (dx=0.157 vs 0.315)
and higher k (pi/10 vs pi/20) gives worse results than the baseline.
fc fluctuates wildly (0.05-0.47) during evolution and ends low.
|P| also oscillates (0.02-0.84) but ends at 0.034. The tighter
confinement from periodic BC does not help — the braid still disperses
transversely through the absorbing x,y boundaries.

**Transient peak**: At t=400, fc=0.47 and |P|=0.13 — tantalizingly close
to the joint target. But by t=500, both have degraded. The braid is
slowly losing energy to the absorbing boundaries.

## M1c: Propagation Speed Scan (L=20, lambda_pw=0.5)

| k (twists/2L) | v_g   | fc(500) | |P|(500) | Pz(500) | PASS? |
|----------------|-------|---------|----------|---------|-------|
| pi/20 (0.5)    | 0.155 | 0.185   | 0.000    | 17.0    | NO    |
| 2pi/20 (1.0)   | 0.300 | 0.182   | 0.002    | 69.8    | NO    |
| 4pi/20 (2.0)   | 0.532 | 0.197   | 0.007    | 316.7   | NO    |
| 6pi/20 (3.0)   | 0.686 | 0.166   | 0.012    | 487.8   | NO    |

**Finding**: Higher k (more twists, faster propagation) slightly improves
|P| survival but all values are still dead (< 0.02). This is because
lambda_pw=0.5 is the dominant killer — it wipes out the braid structure
regardless of propagation speed.

Higher k does maintain more z-momentum Pz (317-488 vs 17-70), indicating
the propagation is real. But the braid topology (triple product) is gone.

## Key Insight

The fundamental tension is:
- **|P| requires** the three fields to maintain distinct phases (2pi/3 offset)
- **fc requires** the energy to stay localized in the transverse plane
- **Absorbing boundaries** continuously drain transverse energy
- **Pairwise coupling** synchronizes the fields, killing the phase offset

The DynA control (lpw=0) is actually the **best configuration tested**.
It maintains |P|=0.574 but fc=0.28 because energy slowly leaks transversely.

## Recommendations for Next Steps

1. **Do NOT use pairwise coupling** for propagating braids — it kills |P|
2. **Try periodic BC in ALL directions** (not just z) to prevent energy loss
3. **Try stronger triple-product coupling** (mu=-50, kappa=50) to deepen the well
4. **Try a confining external potential** V_ext = alpha * r_perp^2 to prevent
   transverse spreading without affecting the phase structure
5. **M1b transient** (fc=0.47, |P|=0.13 at t=400) suggests that with less
   energy loss, the joint target IS reachable

## Parameters Used

All tests: mu=-20, kappa=20, m=1.0, A0=0.8, R_tube=3.0, N=128, t=500, cfl=0.20
BC: periodic in z, absorbing in x,y (R_abs_inner=0.7L, R_abs_outer=0.95L)
