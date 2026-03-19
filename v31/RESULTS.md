# V31 Results: Three Tests of the M7+c(rho_B) Gravity Model

## Model
    S (braid): d2S/dt2 = c_eff2(rho_B) nabla2 S - m2 S - dV/dS - g(SumB2)S
    B (background): d2B/dt2 = nabla2 B - m2 B - g(SumS2)B
    c_eff2 = 1 - alpha_c(1 - rho_B/rho0)
    m=1.50, mu=-41.3, kappa=50, g=0.01, A_bg=0.1

---

## T1: Radial c_eff Profile

**Setup**: N=128, L=30, T=400, alpha_c=0.2. Run A (braid+bg) vs Run B (bg only, control).
Time-averaged over T=200-400 (13 snapshots each).

**Power-law fits** (r=5..20):
| Quantity | Exponent alpha | Amplitude | Sign |
|----------|---------------|-----------|------|
| delta_rhoB | 1.275 | 4.74e-3 | POSITIVE |
| delta_c2   | 1.275 | 2.80e-2 | POSITIVE |

**Radial structure** (differential, braid effect on B):
- r < 5: delta_rhoB NEGATIVE (depletion), delta_c2 negative -> c_eff < 1 at core
  - Core c2_A = 1.056, c2_B = 1.072 -> delta = -0.015
  - Braid CONSUMES B near its center via g*S2*B coupling
- r = 5-10: crossover from negative to positive
- r > 10: delta_rhoB POSITIVE (accretion), delta_c2 positive -> c_eff slightly higher
  - At r=12: delta_c2 = +0.0017 (small but positive)
  - B-field displaced from braid core accumulates in the surrounding shell

**Interpretation**: The braid creates a c_eff profile with c < 1 at the core and c slightly > 1
in a shell. This is NOT a monotone 1/r depletion well. The sign of delta_c2 at intermediate
radius is WRONG for gravitational attraction -- the shell has c > 1, meaning test waves
propagate FASTER there, which is REPULSIVE (anti-gravity).

The exponent 1.275 is between 1/r and 1/r^1.5, consistent with V29's finding (1.21).
The power law is promising but the sign is fatal.

---

## T2: Two-Braid Attraction Test (THE GRAVITY TEST)

**Setup**: N=128, L=40, T=400. Two braids at (+-10, 0, 0). Three configs.

### Results by configuration:

**alpha_c = 0.0 (control, no c-coupling)**:
- D(t) = 20.0 +/- 0.3 throughout (stable oscillation)
- fc = 0.87-0.94 (braids survive to T=400)
- No attraction, no repulsion. Braids are NEUTRAL as expected.

**alpha_c = 0.2 (moderate c-coupling)**:
- D(t) increases from 20 to 48 over T=400 (STRONG REPULSION)
- fc drops from 0.99 to 0.02 (braids fully DISSOLVE)
- avg_c2 at midpoint hits cap of 1.5 (massive c-enhancement between braids)
- Mechanism: B accumulates between braids -> c > 1 -> braids pushed apart

**alpha_c = 0.5 (strong c-coupling)**:
- D(t) oscillates erratically between 16 and 33
- fc drops from 1.0 to 0.47 (partial dissolution)
- Same mechanism with more chaotic dynamics

### VERDICT: NO GRAVITATIONAL ATTRACTION

The c(rho_B) coupling produces REPULSION, not attraction. The mechanism:
1. Each braid consumes B at its core (via g*S2*B coupling)
2. This DEPLETES B near the braid -> c_eff < 1 at braid location
3. The midpoint between braids has NO S -> B is UNDEPLETED -> c_eff >= 1
4. The c-contrast (low at braids, high between) pushes braids apart

For attraction, we would need c_eff < 1 BETWEEN the braids (a shared depletion
well). Instead, the depletion is LOCAL to each braid, and the gap between them
is normal or enhanced.

---

## T3: Pinned Boundaries -- Accretion Dynamics

**Setup**: N=128, L=40, T=500, alpha_c=0.2. B pinned at domain edges (reservoir).
S uses normal absorbing damping.

**Accretion data**:
| Time | E_B_total | E_B(r<15) | E_S | fc |
|------|-----------|-----------|-----|-----|
| 0    | 1.77e4    | 1933      | 2190 | 1.00 |
| 100  | 1.71e4    | 1921      | 1365 | 0.87 |
| 200  | 1.62e4    | 1926      | 1522 | 0.90 |
| 300  | 1.51e4    | 2017      | 2289 | 0.77 |
| 400  | 1.40e4    | 2143      | 1341 | 0.94 |
| 500  | 1.28e4    | 2176      | 1649 | 0.92 |

E_B(r<15) increases from 1933 to 2176 (+12.6%) while E_B_total drops by 28%.
B energy flows INWARD toward the braid center -- genuine accretion.

**Profile evolution** (core r=0.25):
| Snapshot | t   | rhoB(core) | c2_eff(core) | Regime |
|----------|-----|------------|--------------|--------|
| 1        | 100 | 0.031      | 0.982        | Depletion |
| 2        | 200 | 0.029      | 0.973        | Depletion (deeper) |
| 3        | 300 | 0.034      | 0.998        | Near-unity (transition) |
| 4        | 400 | 0.042      | 1.047        | Accretion |
| 5        | 500 | 0.044      | 1.060        | Accretion (stronger) |

The core INVERTS from depletion (c=0.982) to accretion (c=1.060). The braid
pulls B toward itself, filling the initial depletion well with accreted B.

**Profile shape at T=500**: rhoB peaks at core (0.044), decreases monotonically
to ~0.019 at r=25 (far below rho0=0.034). The pinned edges maintain rho0 but
the intermediate zone is deeply depleted. The braid creates an accretion core
surrounded by a depletion shell -- the INVERSE of a gravitational well.

**NOT steady state**: Profile still evolving at T=500. Would need T>>1000 and
larger domain for true steady state.

---

## Root Cause Analysis

The M7+c(rho_B) gravity hypothesis fails because of a sign problem:

**Expected for gravity**: Braid depletes B -> low c around braid -> other braids
slow down near it -> geodesics curve inward -> attraction.

**What actually happens**: Braid depletes B locally (core), but the depleted B
re-accumulates in a shell. At distances relevant for two-braid interaction (D=20),
the B density is HIGHER than background (accretion zone), so c > 1 -> repulsion.

The core depletion IS real (T1 shows delta_rhoB < 0 at r < 5), but it's too
localized. The power-law tail (r^-1.275) is in the ACCRETION regime (positive),
not the depletion regime (negative).

**Why depletion is localized**: The g*S2*B coupling removes B only where S is
large (the braid core, r < 3). The removed B energy propagates outward as waves
and deposits in the surrounding shell. The net effect at distance is ENRICHMENT.

This is analogous to a star that heats its surroundings: the core is hot (B-depleted)
but the surrounding medium is warmer than ambient (B-enriched). For gravity, we
would need the opposite: a cold spot (B-depleted) with a long-range 1/r cold tail.

## What Would Need to Change

For the c(rho_B) mechanism to produce gravity, one of these must change:

1. **Invert the c formula**: c_eff2 = 1 + alpha_c(1 - rho_B/rho0). Then accretion
   (rhoB > rho0) gives c < 1 and depletion (rhoB < rho0) gives c > 1.
   The accretion core would then be a slow pocket (attractive).

2. **Larger separations**: At D >> accretion shell radius (~10), both braids
   would be in each other's far-field depletion zones. But the effect decays
   as r^-1.275, so it would be extremely weak.

3. **Asymmetric coupling**: If B->S conversion (absorption) were stronger than
   S->B conversion (re-emission), the net effect would be depletion, not accretion.

4. **Different field for c**: Use a separate scalar field phi that responds to
   braid topology (not just S2), creating a genuine 1/r depletion potential.

## Status

V31 T2 was NEGATIVE with the original c formula (c lower where B depleted).

---

## Sign-Fix Test: INVERTED c PRODUCES GRAVITATIONAL ATTRACTION

Tested 4 fixes for the wrong-sign problem. The INVERTED c formula works:

    c_eff² = 1 + α_c × (1 - ρ_B/ρ₀)
    (c LOWER where B enriched, HIGHER where B depleted)

| Fix | D(0) | D(final) | ΔD | fc | Verdict |
|-----|------|----------|-----|-----|---------|
| control (α_c=0) | 20 | 20.3 | +0.3 | 0.35 | Neutral |
| fix1 (D=60) | 60 | 50.1 | -9.9 | 0.01 | Approach (dead braids) |
| fix2 (asymmetric g) | 20 | 47.1 | +27.1 | 0.01 | Repulsion (dead) |
| **fix3 (inverted, α_c=0.2)** | **20** | **14.7** | **-5.3** | **0.63** | **ATTRACTION** |
| **fix4 (inverted, α_c=0.5)** | **20** | **9.3** | **-10.8** | **0.72** | **STRONG ATTRACTION** |

**FIRST GRAVITATIONAL ATTRACTION FROM FIELD DYNAMICS IN THE PROJECT.**

The inverted formula makes c LOWER where B is enriched at the braid core.
The accretion zone becomes a slow pocket → other braids curve toward it.
Stronger α_c → stronger attraction (ΔD scales with α_c). Both braids
survive with fc=0.63-0.72 (better than control's 0.35!).

### The Working Gravity Model

    S: ∂²S_a/∂t² = c_eff²(ρ_B)∇²S_a - m²S_a - ∂V/∂S_a - g(ΣB²)S_a
    B: ∂²B_a/∂t² = ∇²B_a - m²B_a - g(ΣS²)B_a
    c_eff² = 1 + α_c(1 - ρ_B/ρ₀)
    m=1.50, μ=-41.3, κ=50, g=0.01, α_c=0.2-0.5

## Files
- src/v31_t1.c, src/v31_t2.c, src/v31_t3.c, src/v31_signfix.c
- data/t1_profiles.tsv, data/t1_fit.tsv
- data/t2_separation.tsv, data/t2_timeseries.tsv
- data/t3_profiles.tsv, data/t3_accretion.tsv
- data/signfix_results.tsv
