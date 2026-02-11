# Avenue A: One-Loop B^0 p Effective Coupling

## Summary

**Question**: Does the one-loop quantum correction from pion fluctuations on
the Skyrmion background generate an effective B^0 p coupling (g_top)?

**Answer**: NO. The coupling is exactly zero at one loop, for two independent reasons:

1. The constraint vertex |q|^2 + epsilon_0^2 p^2 = rho_0^2 generates only a p^2
   (mass shift) coupling, not a linear B^0 p coupling.

2. Avenue C established that epsilon_0^2 = 0 EXACTLY (topological necessity:
   pi_3(S^7) = 0 means solitons cannot exist if epsilon_0^2 != 0). At epsilon_0^2 = 0,
   the Lagrangian factorizes and no pion-p vertex exists. The one-loop diagram
   vanishes identically regardless of the mode spectrum.

**Combined four-avenue result**: g_top = 0 at all orders in perturbation theory
in the self-consistent (N=0, epsilon_0^2=0) scenario. The WZW route (N=1) gives
10^{13} overcoupling and cannot cross-couple to the degenerate sector at
epsilon_0^2 = 0. Gravity, if it arises from this framework, must be nonperturbative
(instanton tunneling or condensate formation).

**Assessment**: The perturbative program is definitively closed. See Sections 9-12
for the combined four-avenue synthesis.

---

## 1. Setup

The Cl+(3,0,1) field Psi = (s, f1, f2, f3; j1, j2, j3, p) has bulk sector q
and degenerate sector (J, p). The constraint |q|^2 + epsilon_0^2 p^2 = rho_0^2
couples the pseudoscalar p to the bulk Skyrmion strain S_0(r).

### Parameters (e=1, rho_0=1, sigma-model)
- Profile: `data/profiles/profile_sigma_e1.dat`
- Shooting parameter: a = -f'(0) = 1.420071
- E_2 = 51.55, E_4 = 51.59, E_sol = 103.14 (virial: E_2/E_4 = 0.999)
- Lambda = 141.6 (moment of inertia)
- m_pi = 0.398 code^{-1} (physical pion mass)
- Constraint eigenvalue: E_0 = -1.487 (from varsig.c)

### Key Overlap Integrals
Computed numerically from the equilibrium profile:

| Integral | Value | Physical meaning |
|----------|-------|------------------|
| I_{SB} = int S_0 B^0 4pi r^2 dr | 2.508 | strain x baryon density overlap |
| I_{S2} = int S_0^2 4pi r^2 dr | 160.94 | strain self-overlap |
| I_{B2} = int (B^0)^2 4pi r^2 dr | 0.0433 | baryon self-overlap |
| I_{TB} = int T_0 B^0 4pi r^2 dr | 2.837 | Skyrme density x B^0 |
| Q = int B^0 4pi r^2 dr | 1.000 | topological charge (exact) |

---

## 2. Constraint Vertex Analysis

The constraint |q|^2 + epsilon_0^2 p^2 = rho_0^2 implies q = q_0 sqrt(1 - epsilon_0^2 p^2/rho_0^2).

Expanding the kinetic energy L_2 = (rho_0^2/2)|d(q/rho_0)|^2 to second order in p:

    L_2 = (rho_0^2/2) S_0(r) - (epsilon_0^2/2) S_0(r) p^2 + O(epsilon_0^4 p^4)

**The constraint vertex is EVEN in p**: it couples pion pairs to p pairs with
strength epsilon_0^2 S_0(r)/rho_0^2. This is a pi-pi-p-p four-point vertex.

At one loop, contracting the pion pair into a loop gives:

    L_{eff} = -epsilon_0^2 G_pi(r,r) S_0(r) p^2 / rho_0^2

This is a **p mass shift**, not a B^0 p coupling:

    delta m_p^2 / epsilon_0^2 = I_{S2} / (16 pi^2 m_pi^2 rho_0^2) = 6.43

**Critical observation**: The constraint vertex cannot generate a LINEAR p
coupling at one loop because the vertex is quadratic in p. A linear B^0 p
term requires a P-odd (pseudoscalar) source.

---

## 3. WZW Vertex: The Only P-odd Source

The Wess-Zumino-Witten term is the only part of the Lagrangian that is odd
under p -> -p. For the Cl+(3,0,1) theory:

    L_WZW = (Nc / 240 pi^2) epsilon^{mu nu rho sigma} B_mu J^D_{nu rho sigma}

where B_mu is the baryon current and J^D involves the degenerate sector.

**For a STATIC hedgehog**: The WZW coupling B^0 x (pion x dpion) . dp/dr
vanishes because (a) p = 0 on the equilibrium, and (b) the time-component
B^0 has no spatial current to couple to dp/dr. This was proved rigorously
in the Path 5 analysis (Track A, three independent obstructions).

**At one loop**: The pion loop could generate <pi x dpion> != 0, but this
requires a B^0 insertion vertex that already involves the coupling we're
trying to compute. The diagram is circular.

**CONCLUSION**: At one loop with the constraint coupling alone, g_top^{1-loop} = 0.
The constraint gives p^2, and the WZW vanishes on statics. No linear B^0 p
coupling is generated.

---

## 4. What One-Loop DOES Generate

Although the B^0 p coupling vanishes at one loop on a static background,
the pion fluctuations do produce physically meaningful corrections:

### (a) p-field mass shift

    delta m_p^2 / epsilon_0^2 = 6.43

This is an O(1) shift relative to the tree-level constraint eigenvalue |E_0| = 1.487.
The pion fluctuations significantly modify the p-field effective potential.

### (b) Casimir energy

    E_{Casimir} ~ -3 I_{S2} / (64 pi^2 m_pi) ~ -0.64 code units ~ -5.8 MeV

This is the zero-point energy of three pion flavors on the soliton background.
It contributes to the quantum correction of the soliton mass:

    M_{quantum} = M_{classical} + E_{Casimir} + hbar^2 J(J+1)/(2 Lambda) + ...

### (c) One-loop correction to moment of inertia

    delta Lambda / Lambda ~ (3/16 pi^2) x S_0(0) x R^3 / m_pi ~ 0.97

This is an ~97% correction, indicating that perturbation theory is marginal
at one loop. Higher-order effects and proper resummation are important.

---

## 5. g_top from Rotating Soliton (Avenue A + B Connection)

The physical nucleon has J = I = 1/2, achieved by quantizing the collective
coordinate A(t) in SU(2). The quantum angular velocity:

    Omega = sqrt(3/4) hbar / Lambda = 6.12 x 10^{-3} code^{-1}

For the rotating soliton, the WZW term becomes nonzero because the spatial
baryon current B^i is nonzero:

    B^i_{rot}(x) ~ (Omega/c) (x cross ...)^i B^0(r)

### WZW Coupling Estimates

    Nc/(240 pi^2) = 1.27 x 10^{-3}
    hbar/Lambda = 7.06 x 10^{-3}

    g_WZW = Nc/(240 pi^2) x hbar/Lambda = 8.95 x 10^{-6}

### Comparison with Gravity

    Required g_top for G_Newton: 2.81 x 10^{-17}
    WZW estimate: 8.95 x 10^{-6}
    Ratio: g_WZW / g_required = 3.2 x 10^{11}

**The WZW coupling is 10^{11} times too strong for gravity.**

### Angular Integral Suppression

The WZW coupling to spherically symmetric p(r) vanishes:
(Omega x x) . grad(p(r)) = 0 for s-wave p.

The coupling requires L=1 (dipolar) p fluctuations. This introduces additional
angular integrals and form factors that could provide significant suppression.
However, even a factor of 10^{11} suppression from angular integrals seems
implausible from a single geometric factor.

---

## 6. Channel C Near-Binding Enhancement

The K=1, Channel C (mixed L=1, I=1) angular mode has an eigenvalue at
omega^2 ~ 1.04 m_pi^2, which is 4% above the pion mass threshold.

### Propagator Enhancement

At the near-resonance, the pion propagator in the L=1 channel is enhanced:

    G_pi(omega) ~ 1/(omega^2 - m_pi^2) ~ 26/m_pi^2

This enhances the L=1 partial wave of the loop integral by a factor of ~26.

### Could Quantization Close the Gap?

If quantum corrections pushed omega^2 below m_pi^2, Channel C would become
a true bound state, creating a pole in the pion propagator. This would
dramatically enhance the loop integral.

Required shift: delta(omega^2/m_pi^2) = -0.04
Typical quantum shift: delta_omega^2 ~ hbar/(Lambda R^2) ~ 0.003-0.005

**The quantum shift is 4-8x too small to close the gap.** The near-binding is
tantalizing but insufficient for a qualitative change.

---

## 7. Feynman Diagram Structure (Detailed)

For completeness, here is the correct diagram analysis:

### Diagram 1: Constraint tadpole (VANISHES for linear p)
```
        pi loop
       /       \
p --- * -------- * --- p
      |constraint|
```
Result: p^2 mass shift. No linear B^0 p.

### Diagram 2: WZW + constraint (TWO-LOOP)
```
        pi loop
       /       \
p --- * -------- * --- B^0
      |WZW      |constraint|
```
Result: Requires two vertices and is suppressed by additional 1/(16 pi^2).
Estimate: g_top ~ epsilon_0^2 x g_WZW x I_{SB}/(16 pi^2 I_{S2}) ~ 8.8 x 10^{-10} epsilon_0^2

### Diagram 3: WZW on rotating soliton (ONE-LOOP, non-static)
```
p --- * --- B^i_{rot}
      |WZW|
```
Result: g_WZW = Nc/(240 pi^2) x hbar/Lambda ~ 10^{-5}. Too large.

### Diagram 4: Anomaly-induced (requires Avenue C)
If the Z_2 symmetry p -> -p has a quantum anomaly, a linear p term
could be generated non-perturbatively. This requires the anomaly
classification of Avenue C.

---

## 8. Summary and Implications (Original Analysis)

### Key Results

| Finding | Result |
|---------|--------|
| Constraint one-loop | Generates p mass shift (delta m^2/epsilon_0^2 = 6.43), NOT B^0 p |
| WZW on statics | Vanishes (proved in Path 5) |
| WZW rotating | g_WZW ~ 10^{-6}, too large by 10^{11} |
| Channel C enhancement | 26x in L=1 channel, insufficient for qualitative change |
| Required epsilon_0^2 | Not determined by one-loop (free parameter) |

---

## 9. Impact of Avenue C: epsilon_0^2 = 0 Exactly (Critical Update)

### 9.1 The topological obstruction

Avenue C established a hard topological constraint: if epsilon_0^2 != 0 at the
algebraic level, the norm constraint |q|^2 + epsilon_0^2 |d|^2 = rho_0^2 defines
a target space S^7 (or squashed S^7). Since pi_3(S^7) = 0, ALL topological solitons
cease to exist. Therefore:

    epsilon_0^2 = 0 EXACTLY at the fundamental (algebraic) level.

This is not a perturbative approximation — it is a topological necessity for the
existence of Skyrmions. The PGA interpretation (e_0^2 = 0 in Cl+(3,0,1)) is thus
the ONLY consistent choice.

### 9.2 Consequences for Avenue A

**The constraint mechanism is dead.** All constraint-mediated effects computed
in Sections 2 and 4 are proportional to epsilon_0^2 and therefore exactly zero:

| Effect | Dependence | Value at epsilon_0^2 = 0 |
|--------|-----------|--------------------------|
| p mass shift | ~ epsilon_0^2 I_S2/(16pi^2 m_pi^2) | **0** |
| Casimir correction to p potential | ~ epsilon_0^2 G_pi | **0** |
| Two-loop B^0 p via constraint+WZW | ~ epsilon_0^2 g_WZW / (16pi^2) | **0** |
| p-field kinetic term from constraint | ~ epsilon_0^2 | **0** |

The constraint eigenvalue E_0 = -1.487 (Section 1) exists as a mathematical
property of the S_0(r) potential, but it describes no physical mode because
at epsilon_0^2 = 0, the degenerate sector has no coupling to the bulk.

### 9.3 Does a DYNAMICALLY GENERATED effective epsilon_0^2 survive?

**Question**: Even at epsilon_0^2 = 0 fundamentally, could quantum fluctuations
generate an effective nonzero coupling?

**Answer: No, for the following reasons.**

**(a) No vertex exists to mediate the coupling.**
At epsilon_0^2 = 0, the Lagrangian factorizes: L = L_bulk[q] + L_degen[d]. There
is no term connecting the bulk sector q to the degenerate sector d = (J, p).
Without a tree-level vertex, no perturbative loop diagram can generate one — every
loop diagram requires at least one vertex to start from.

**(b) The field space factorizes.**
The field space is S^3 (bulk) x R^4 (degenerate), a direct product. The
degenerate fields live in a contractible fiber over each point of the bulk
configuration space. There are no nontrivial cycles connecting the two sectors,
so topological effects (WZW, instantons) also factorize.

**(c) <p^2> = 0 does not generate epsilon_0^2.**
Quantum fluctuations of the free degenerate scalar give <p^2> != 0 (UV divergent
vacuum energy). But this is a CONSTANT — it does not depend on the bulk field
configuration q. There is no soliton-profile-dependent enhancement of <p^2> because
p does not couple to the soliton at epsilon_0^2 = 0. The vacuum energy from <p^2>
is absorbed into the cosmological constant, not into a B^0 p vertex.

**(d) The grade mismatch is algebraically exact.**
The scalar extraction <(partial Psi)(partial Psi~)>_0 = |partial q|^2 contains
no d-dependent terms because (partial q)(partial d~) has odd grade (contains e_0
linearly). This is an algebraic identity of the Cl+(3,0,1) Clifford algebra,
exact at all loop orders. Quantum corrections cannot produce a nonzero grade-0
component where the algebra forbids one.

### 9.4 The WZW route at epsilon_0^2 = 0

The WZW term L_WZW = (N/240 pi^2) integral Tr(L^5) exists entirely in the bulk
sector at epsilon_0^2 = 0. Avenue B established that:

1. The WZW 5-form lives on S^3, not on the product S^3 x R^4.
2. Without a cross-sector term in the Lagrangian, the WZW term cannot couple
   to the degenerate scalar p.
3. The spatial coupling B^i d_i p vanishes at O(Omega) by current conservation
   (div B = 0), even if a coupling existed.
4. For SU(2): pi_5(SU(2)) = Z_2, so N in {0, 1}. At N=1, the coupling would be
   10^{13} too strong anyway.

**At epsilon_0^2 = 0, the WZW term has no mechanism to source p.**

### 9.5 The one-loop diagram at epsilon_0^2 = 0

To be maximally explicit, here is the one-loop B^0 p effective vertex at
epsilon_0^2 = 0:

```
DIAGRAM: pion loop with B^0 insertion and external p
```

The loop integral is:
    Gamma_1[B^0, p] = Sum_n integral d^4k V_{B^0-pi}(k,n) G_pi(k,n) V_{pi-p}(k,n)

At epsilon_0^2 = 0:
- V_{B^0-pi}: exists (topological current couples to pion fluctuations) = NONZERO
- G_pi: pion propagator on soliton background = NONZERO
- V_{pi-p}: pion-to-p vertex from constraint = ZERO (no coupling)

Since V_{pi-p} = 0, the ENTIRE loop vanishes regardless of the mode spectrum.

    g_top^{1-loop} = 0 at epsilon_0^2 = 0.

This is exact — there is no missing vertex, no subtlety, no way around it.

---

## 10. Impact of Avenue B: WZW Overcoupling

Avenue B independently showed that even if a coupling existed:

1. B^i d_i p vanishes at O(Omega) (current conservation).
2. The static B^0 p coupling reduces to Path 3 with g_WZW = N/(240 pi^2).
3. For N=1: g_WZW = 4.2 x 10^{-4}, giving G_eff ~ 10^{26} G_Newton (10^{13} too strong).
4. For N=0: no coupling at all.

Combined with Avenue C (epsilon_0^2 = 0), Avenue B shows the WZW term CANNOT
cross-couple to the degenerate sector. The consistent scenario is N=0 or N=1
with no cross-sector coupling.

---

## 11. Impact of Avenue D: Spectral Null Result

Avenue D showed that:

1. At epsilon_0^2 ~ 10^{-14}, all spectral effects are O(10^{-14}) — invisible.
2. The nuclear spectrum (masses, splittings) is independent of epsilon_0^2.
3. Channel C near-binding gap WIDENS with epsilon_0^2 > 0 (wrong direction).
4. No spectral self-consistency condition constrains epsilon_0^2.
5. At epsilon_0^2 = 0, the quantum mass formula has no degenerate sector contribution.

---

## 12. Synthesis: Combined Four-Avenue Result

### 12.1 All perturbative mechanisms are excluded

| Avenue | Mechanism | Result |
|--------|-----------|--------|
| A (one-loop) | Constraint vertex generates p^2, not B^0 p | NULL (p mass shift only) |
| A (one-loop at eps_0^2=0) | No pion-p vertex exists | **ZERO exactly** |
| B (WZW rotation) | B^i d_i p vanishes by div B = 0 | NULL (all partial waves) |
| B (WZW static) | B^0 p coupling at N=1 is 10^{13} too strong | OVERCOUPLED or ABSENT |
| B+C (WZW at eps_0^2=0) | No cross-sector term in factorized Lagrangian | **ZERO exactly** |
| C (anomaly) | No gauge anomalies, no chiral anomaly, no mixed anomaly | NULL (no constraint) |
| C (topology) | pi_3(S^7)=0 forces epsilon_0^2=0 | **HARD CONSTRAINT** |
| D (spectrum) | Nuclear observables blind to eps_0^2 ~ 10^{-14} | NULL (decoupled) |

### 12.2 The nonperturbative escape

At epsilon_0^2 = 0 and N=0 (the self-consistent scenario):
- g_top = 0 at ALL orders in perturbation theory.
- Gravity, if it arises from this framework, MUST be nonperturbative.

Candidate nonperturbative mechanisms:

1. **Instanton tunneling**: Euclidean configurations connecting B=0 and B=1
   sectors could generate g_top ~ e^{-S_inst}. For S_inst ~ 32:
   g_top ~ e^{-32} ~ 10^{-14}. The instanton action depends on (rho_0, e, lambda),
   so this WOULD predict G_Newton from nuclear parameters.

2. **Condensate formation**: A <p^2> condensate forming through nonperturbative
   dynamics (analogous to <qq> in QCD) could generate effective epsilon_0^2_eff.
   But at epsilon_0^2 = 0, there is no potential to drive condensation — the
   degenerate sector is completely free and decoupled.

3. **Emergent graviton**: If the BLV effective metric (P/m ≠ 2 from L_6 sextic)
   becomes dynamical at quantum level, a spin-2 mediator could emerge as a
   collective excitation of the soliton lattice, rather than from the scalar p.

### 12.3 The fundamental tension

The combined four-avenue investigation reveals a structural issue:

- **Classical level**: The B^0 p coupling (Path 3) gives exact 1/r gravity,
  but g_top is free.
- **Algebraic level**: epsilon_0^2 = 0 is forced by topology (pi_3(S^7) = 0).
  No perturbative mechanism can generate g_top.
- **Quantum level**: All loop corrections vanish at epsilon_0^2 = 0. The WZW
  term either overcouples (N=1) or is absent (N=0).

The hierarchy problem (g_top ~ 10^{-17}) appears to require physics beyond
perturbative Cl+(3,0,1). Either:
- Gravity is fundamentally nonperturbative in this framework
- Additional structure (gauge fields, extra dimensions) is needed
- The B^0 p mechanism is not the correct source of gravity

### 12.4 What the one-loop analysis DID establish

Despite the null result for g_top, the one-loop analysis produced valuable results:

1. **Correct diagram identification**: The constraint vertex is p^2 (even in p),
   not linear. Only P-odd (WZW) terms can generate B^0 p.

2. **Quantitative overlap integrals**: I_SB = 2.508, I_S2 = 160.94, I_B2 = 0.0433.
   These characterize the soliton background for any future calculation.

3. **Channel C near-binding**: 26x enhancement in the L=1 pion propagator, with
   a gap of 4-5.7% above threshold. Could be closed by L_6 or finite-lambda effects.

4. **Perturbation theory breakdown**: One-loop corrections to the WZW vertex are
   O(1) (~97%), indicating that perturbative methods are marginal for this system.

5. **Casimir energy**: E_Casimir ~ -0.64 code ~ -5.8 MeV, contributing to the
   quantum soliton mass.

---

## Code and Reproducibility

**Source**: `src/oneloop.c`
**Build**: `make oneloop`
**Run**: `./bin/oneloop -profile data/profiles/profile_sigma_e1.dat [-mpi 0.398] [-Nc 3]`

All numerical results verified against existing tools:
- E_sol = 103.14 matches radial.c
- Lambda = 141.6 matches normal_modes.c
- Q = 1.000 (topological charge exactly conserved)
- Virial E_2/E_4 = 0.999 (equilibrium verified)
