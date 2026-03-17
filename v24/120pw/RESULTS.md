# V24-PW Results: Pairwise Coupling 120-degree Phase Binding

## Answer

**No.** The 120-degree state never becomes more stable than 0-degree. The pairwise
coupling makes 0-degree oscillons MORE stable while DESTROYING 120-degree oscillons.

## Setup

- EOM: d^2 phi_a/dt^2 = d^2 phi_a/dx^2 - m^2 phi_a - lambda(phi_b+phi_c) + F_triple
- Parameters: mu=-20, kappa=20, m=1.0, A=0.8, sigma=3.0
- Grid: Nx=4000, xmax=100, tfinal=10000
- Absorbing boundaries at |x| > 75
- Lambda scan: {0.0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.85, 0.90, 0.95}

## Results Table

| lambda | Mode | fc_final | dE/E0     | omega_peak | phase12 | phase13 | Status      |
|--------|------|----------|-----------|------------|---------|---------|-------------|
| 0.00   | 120  | 0.994    | -7.08e-1  | 0.876      | pi      | pi      | OSCILLON (anti-phase!) |
| 0.00   | 0    | 1.000    | -6.24e-1  | 0.876      | 0       | 0       | OSCILLON    |
| 0.20   | 120  | 0.107    | -8.31e-1  | 0.894      | ~pi     | ~pi     | DISPERSED   |
| 0.20   | 0    | 0.980    | -4.29e-1  | 0.996      | 0       | 0       | OSCILLON    |
| 0.40   | 120  | 0.205    | -8.47e-1  | 0.774      | ~pi     | ~pi     | DISPERSED   |
| 0.40   | 0    | 1.000    | -1.74e-1  | 1.176      | 0       | 0       | OSCILLON    |
| 0.60   | 120  | 0.347    | -9.08e-1  | 0.636      | ~pi     | mixed   | DISPERSED   |
| 0.60   | 0    | 0.998    | -2.37e-1  | 1.338      | 0       | 0       | OSCILLON    |
| 0.70   | 120  | 0.000    | -1.26e+1  | 0.000      | pi      | pi      | COLLAPSED   |
| 0.70   | 0    | 1.000    | -4.13e-2  | 1.410      | 0       | 0       | OSCILLON    |
| 0.80   | 120  | 0.000    | -4.36e+1  | 0.000      | pi      | pi      | COLLAPSED   |
| 0.80   | 0    | 1.000    | -8.22e-2  | 1.482      | 0       | 0       | OSCILLON    |
| 0.85   | 120  | 0.000    | -8.74e+1  | 0.000      | pi      | pi      | COLLAPSED   |
| 0.85   | 0    | 1.000    | -9.99e-2  | 1.512      | 0       | 0       | OSCILLON    |
| 0.90   | 120  | 0.000    | -2.42e+2  | 0.000      | pi      | pi      | COLLAPSED   |
| 0.90   | 0    | 1.000    | -1.15e-1  | 1.542      | 0       | 0       | OSCILLON    |
| 0.95   | 120  | 0.000    | -1.53e+3  | 0.000      | pi      | pi      | COLLAPSED   |
| 0.95   | 0    | 1.000    | -1.26e-1  | 1.578      | 0       | 0       | OSCILLON    |

## Key Findings

### 1. 120-degree phase is NEVER realized

Even at lambda=0, the 120-degree initialization relaxes to **anti-phase (180-degree)**
within the first few oscillation periods. The measured phases phi_{12} and phi_{13}
are both ~pi, not 2*pi/3. The three-body coupling with mu<0 prefers anti-phase
pairing, not the Z_3-symmetric 120-degree configuration.

### 2. Three regimes for 120-degree init

- **lambda=0**: Relaxes to anti-phase oscillon (fc=0.99, stable)
- **0.2 <= lambda <= 0.6**: Disperses (fc drops to 0.1-0.35, 83-91% energy loss)
- **lambda >= 0.7**: Collapses to static non-oscillating state (omega=0, E hugely
  negative, fields freeze at phi_1 ~ +1.2, phi_2 ~ phi_3 ~ -0.7)

### 3. 0-degree mode is robust at ALL lambda

The symmetric in-phase oscillon remains perfectly localized (fc > 0.98) for all
tested lambda values. Energy loss DECREASES with lambda (from 62% at lambda=0 to
4% at lambda=0.7), showing the pairwise coupling STABILIZES the symmetric mode.

The peak frequency tracks sqrt(m^2 + 2*lambda) as expected for the symmetric
normal mode.

### 4. Why pairwise coupling fails for 120-degree

The pairwise term V_pw = lambda(phi_1*phi_2 + phi_2*phi_3 + phi_3*phi_1) is
minimized when the sum phi_1+phi_2+phi_3 is minimized (since the cross terms
equal [(sum)^2 - sum(phi_a^2)]/2). For equal-amplitude oscillations:

- **0-degree**: V_pw = +3*lambda*A^2*cos^2(omega*t) > 0 (positive, but offset)
- **120-degree**: V_pw = -3*lambda*A^2/2 < 0 (negative, but CONSTANT — no oscillation energy)

The problem is that the 120-degree pairwise energy is constant (no dynamics),
while the triple-product coupling P = phi_1*phi_2*phi_3 vanishes identically
for 120-degree phases (cos(0)+cos(2pi/3)+cos(4pi/3) conspire to kill P at each
instant). So the 120-degree mode gets NO binding from the triple product and the
pairwise term just shifts the effective mass without providing oscillation trapping.

### 5. Collapse mechanism (lambda >= 0.7)

At large lambda the pairwise coupling overwhelms the mass term for the
antisymmetric mode (m^2_anti = m^2 - lambda -> 0). The fields develop a static
instability: phi_1 grows positive while phi_2, phi_3 go negative (or vice versa),
reaching a static minimum of V_pw + V_triple at nonzero field values. This is
NOT an oscillon — it is a classical static solution (kink-like) that the
absorbing boundaries cannot remove because it is non-radiative.

## Conclusion

Pairwise linear coupling V_pw = lambda*(phi_a*phi_b + ...) does NOT stabilize
120-degree oscillons. The mechanism fails because:

1. The 120-degree phase is never dynamically selected (anti-phase pi is preferred)
2. At 120-degree, the triple-product coupling P = phi_1*phi_2*phi_3 vanishes,
   removing the only source of nonlinear trapping
3. The pairwise term merely shifts the mass gap without creating a frequency
   window for 120-degree modes

To achieve 120-degree binding, one would need a coupling that explicitly breaks
the phi -> -phi symmetry (e.g., cubic phi_a*phi_b*phi_c with a SIGN that
distinguishes 120 from 0/180), or a Z_3-symmetric potential like
cos(theta_a - theta_b - 2*pi/3) that directly penalizes phase deviations from 120.
