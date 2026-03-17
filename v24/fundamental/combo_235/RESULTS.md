# Combo 2+3+5: Confined Lattice with Inertial Deformation — Results

## Parameters

- Confining potential: V = -sigma*sqrt(P^2+eps^2) + (kc/2)*P^4, sigma=1.0, kc=0.01, eps=1e-6
- Standard potential: V = (mu/2)*P^2/(1+kappa*P^2), mu=-20.0, kappa=20.0
- Pairwise coupling: lambda=0.5
- Chain: 8 oscillons, d=16.0, periodic BC, Nx=2560

## Phase 1: Single Oscillon Equilibration

| Property         | Confining | Standard |
|-----------------|-----------|----------|
| Mass M          | 10.605    | 9.449    |
| Breathing omega | 1.158     | 1.260    |
| omega < m?      | NO (1.158 > 1.0) | NO (1.260 > 1.0) |

Both oscillons have breathing frequencies ABOVE the mass gap (m=1.0),
meaning they are radiative oscillons slowly losing energy. The confining
potential produces a heavier oscillon with a lower breathing frequency, as
expected from the deeper potential well.

NOTE: omega=1.158 for confining is much higher than the omega=0.69
cited in the proposal (which was from Test C with sigma=10, not sigma=1).
At sigma=1.0 the confinement is weak.

## Phase 2: Chain Phonon Spectrum

### Stability

| Metric              | Confining | Standard |
|--------------------|-----------|----------|
| Max displacement   | 63.98     | 9.01     |
| Max disp / d       | 400%      | 56%      |
| E conservation     | -1.9e-3   | -1.1e-2  |

**The confined lattice is LESS stable than the standard lattice.**
Max displacement reaches 400% of d (oscillons pass through each other),
compared to 56% for standard. The confining potential at sigma=1.0 is too
weak to provide stronger inter-oscillon binding than the saturating potential
with mu=-20, kappa=20.

### Phonon Dispersion (symmetric/compression modes)

| q | k       | omega_conf | omega_std |
|---|---------|-----------|-----------|
| 0 | 0.00000 | 0.210     | 0.473     |
| 1 | 0.04909 | 0.048     | 0.047     |
| 2 | 0.09817 | 0.090     | 0.094     |
| 3 | 0.14726 | 0.143     | 0.140     |
| 4 | 0.19635 | 0.189     | 0.190     |

Sound speeds (q=1): confining c_s=0.97, standard c_s=0.96.
The symmetric phonon branches are nearly identical -- the compression mode
is dominated by the pairwise coupling lambda=0.5 and mass gap, not the
nonlinear potential. Both show nearly linear dispersion omega~k with
c_s close to 1.

## Phase 3: Antisymmetric Perturbation

### Antisymmetric (Optical) Branch

| q | k       | omega_conf | omega_std |
|---|---------|-----------|-----------|
| 0 | 0.00000 | 0.001     | 0.426     |
| 1 | 0.04909 | 0.140     | 0.515     |
| 2 | 0.09817 | 0.144     | 0.509     |
| 3 | 0.14726 | 0.147     | 0.481     |
| 4 | 0.19635 | 0.154     | 0.494     |

**The antisymmetric gap is dramatically different:**
- Confining: omega_gap(q=0) ~ 0.001 (essentially zero -- GAPLESS)
- Standard:  omega_gap(q=0) = 0.426 (significant gap)

This is the most striking result. The confining potential produces an
essentially gapless antisymmetric mode, while the standard potential has a
clear optical branch gap at omega ~ 0.43.

### Propagation

**Confining**: The antisymmetric mode undergoes massive PARAMETRIC
AMPLIFICATION. Initial perturbation eps=0.01 grows to peak |A12| ~ 3.5
(amplification factor ~250x). The mode reaches all oscillons rapidly.
Source oscillon "decay ratio" = 472 (amplification, not decay).

**Standard**: The antisymmetric perturbation propagates coherently with
moderate amplification. Peak |A12| ~ 0.1 (amplification ~7x).
Source decay ratio = 2.2 (mild growth).

The confined chain shows dramatically stronger parametric amplification of
the antisymmetric mode, consistent with the near-zero gap frequency.

### Symmetric Displacement (during Phase 3)

Both potentials give similar symmetric phonon speeds (~0.95c), confirming
the compression mode is controlled by lambda, not the nonlinear potential.

## Key Findings

### Q1: Is the confined lattice MORE stable?
**NO.** At sigma=1.0, the confining potential provides weaker binding than
the standard mu=-20, kappa=20 potential. Max displacement is 7x larger.
The standard potential creates deeper effective wells between oscillons.

### Q2: Is the antisymmetric optical branch gap different?
**YES, dramatically.** The confining potential makes the gap essentially
zero (omega ~ 0.001 vs 0.43 for standard). This is because the confining
potential V ~ -|P| has a cusp at P=0 that creates a very flat energy
landscape for antisymmetric deformations (which push P toward zero).

### Q3: Is parametric amplification preserved?
**YES, and it is vastly enhanced.** The near-gapless antisymmetric mode
means energy can easily flow into antisymmetric excitations. Amplification
factor ~250x (confined) vs ~7x (standard).

### Q4: Breathing frequency
Confining omega=1.158, standard omega=1.260. Both are above the mass gap,
so both are radiative. The proposal's claim of omega=0.69 was for sigma=10
(strong confinement), not sigma=1.

## Interpretation

The confining potential at sigma=1.0 is in a regime where:
1. The inter-oscillon interaction is weaker than the standard potential
2. The antisymmetric mode is nearly gapless (soft mode)
3. This soft mode is highly susceptible to parametric excitation
4. The chain is unstable to positional drift

The near-zero antisymmetric gap is physically significant: it means the
"shear" mode costs almost no energy. In a higher-dimensional extension,
this would correspond to a nearly massless tensor excitation -- exactly the
graviton analog one would want. However, the instability (both positional
and amplification) suggests this regime is not physically realizable as a
stable lattice.

**To achieve stable confinement with a reduced gap**, one would need to
increase sigma (stronger confinement) while also increasing kc (stronger
quartic stabilization) to prevent the chain from dissolving. The sweet spot
likely lies at sigma ~ 5-10, kc ~ 0.1-1.0.

## Files

- `src/combo235.c` — source code
- `data/combo235_confine_*.tsv` — confining potential results
- `data/combo235_standard_*.tsv` — standard potential comparison
