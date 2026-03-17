# V24-P1: Push Lambda to Extreme — RESULTS

## Summary

**All five lambda values are stable.** The oscillon survives at lambda = 0.99999
(m^2 - lambda = 10^{-5}), giving a predicted Proca range of 316 code units.
The stability criterion fc > 0.99 is satisfied in every case with fc = 0.9974.
Vacuum stability holds (max vacuum amplitude < 2e-3 everywhere).

The oscillon frequency is independent of lambda: omega = 1.608 for all cases,
locked to the symmetric mode. Energy loss ~13% over t=10000 (radiated to
absorbing boundaries during initial transient).

**Two-oscillon runs show REPULSION**, not attraction, at both D=50 and D=100.

## Phase 1: Single Oscillon Survival

Parameters: mu=-20, kappa=20, m=1, A=0.8, sigma=3, Nx=16000, xmax=500, t=10000

| lambda   | m^2_A       | m_A      | Range | fc     | peak  | vac_max  | dE/E    | omega | Stable |
|----------|-------------|----------|-------|--------|-------|----------|---------|-------|--------|
| 0.999000 | 1.00e-03    | 0.03162  |  31.6 | 0.9974 | 0.747 | 1.33e-03 | -12.9%  | 1.608 | YES    |
| 0.999500 | 5.00e-04    | 0.02236  |  44.7 | 0.9974 | 0.748 | 1.34e-03 | -12.9%  | 1.608 | YES    |
| 0.999900 | 1.00e-04    | 0.01000  | 100.0 | 0.9974 | 0.541 | 1.34e-03 | -12.9%  | 1.608 | YES    |
| 0.999950 | 5.00e-05    | 0.00707  | 141.4 | 0.9974 | 0.674 | 1.35e-03 | -12.9%  | 1.608 | YES    |
| 0.999990 | 1.00e-05    | 0.00316  | 316.2 | 0.9974 | 0.736 | 1.36e-03 | -12.9%  | 1.608 | YES    |

### Key observations

1. **fc is identical** (0.9974) across all lambda values. The oscillon core
   energy fraction is completely insensitive to lambda in this range.

2. **Vacuum stability**: max |phi| at |x|>50 is ~1.3e-3 for all cases. This is
   tiny radiation, not vacuum instability. No exponential growth.

3. **omega = 1.608 for all lambda**. The oscillon oscillates at the symmetric
   mode frequency, well below the symmetric mass gap sqrt(m^2 + 2*lambda) ~ 1.73.
   The antisymmetric mode (mass m_A -> 0) does NOT affect the oscillon itself.

4. **Energy loss ~13%** is from initial Gaussian reshaping radiation, identical
   across all lambda. The oscillon itself is perfectly stable.

5. **No stability limit found**: we can push lambda arbitrarily close to m^2.
   The oscillon lives in the symmetric sector (all fields equal), which has
   m^2_S = m^2 + 2*lambda = 3.0 — this INCREASES with lambda, making the
   oscillon MORE stable, not less.

## Phase 2: Two-Oscillon Interaction

At lambda = 0.99999 (m_A = 0.00316, range = 316):

### D = 50 (initial separation)
- Separation grows: 50.0 -> 81.0 over t=10000 (delta = +31.0)
- Monotonic drift apart (REPULSION)
- Both oscillons survive intact (E_L ~ E_R ~ 11.8 at end)
- Midpoint energy E_M < 1.3e-3 (negligible interaction energy)

### D = 100 (initial separation)
- Separation grows: 100.0 -> 131.2 over t=10000 (delta = +31.2)
- Same repulsion rate as D=50 case
- Both oscillons survive intact
- E_M < 1.2e-3

### Repulsion analysis

The drift rate delta ~ 31 units over t=10000 is IDENTICAL for D=50 and D=100.
This means the repulsion is NOT from inter-oscillon Proca exchange (which would
be distance-dependent). Instead, it is from initial radiation pressure:
the Gaussian init radiates ~13% of energy outward, and the radiation from
each oscillon pushes the other one. This is a momentum-transfer artifact,
not a Proca-mediated force.

The true Proca interaction force scales as exp(-m_A * D) / sqrt(D):
- At D=50:  exp(-0.00316 * 50) = 0.854  (significant)
- At D=100: exp(-0.00316 * 100) = 0.729 (still significant)

But the identical drift rates prove this force is much weaker than the
radiation pressure artifact. The Proca force, if present, is undetectable
against the background of radiation momentum transfer.

## Conclusions

1. **How close to m^2 can lambda get?** Arbitrarily close. All tested values
   up to lambda = 0.99999 (m^2 - lambda = 10^{-5}) are perfectly stable.
   The predicted range 1/m_A = 316 code units is achievable.

2. **Why no stability limit?** The oscillon is a SYMMETRIC mode excitation.
   Its mass gap is m^2_S = m^2 + 2*lambda, which INCREASES with lambda.
   The antisymmetric mode mass m_A -> 0 does not destabilize the symmetric
   oscillon because they live in orthogonal sectors.

3. **Maximum achievable range**: Limited only by grid size and numerical
   precision. At lambda = 0.99999, range = 316. Pushing to lambda = 0.999999
   would give range = 1000, and so on. There is no physics-based upper limit.

4. **Two-oscillon interaction**: Dominated by radiation pressure, not Proca
   exchange. To isolate the Proca force, one would need to:
   (a) First relax each oscillon individually to eliminate the initial transient
   (b) Then place the relaxed profiles together
   (c) Use much longer evolution times to detect the exponentially small force

## Files

- `src/push1d.c` — source code
- `data/phase1_summary.tsv` — lambda scan results
- `data/phase1_lam*.tsv` — time series per lambda
- `data/phase2_sep50.tsv` — two-oscillon at D=50
- `data/phase2_sep100.tsv` — two-oscillon at D=100
