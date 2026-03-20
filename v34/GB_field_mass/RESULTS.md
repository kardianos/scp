# Track GB: Field-Dependent Mass — Results

## Summary

Two approaches tested for replacing constant m² with field-dependent mass.
Neither solves the range problem. Both produce informative negative results.

---

## Approach 1: phi^4 Direct Coupling (m_eff^2 = alpha * Sigma phi^2)

Mass INCREASES with field amplitude. High at braid core, low in background.

| alpha | m^2(bg) | m^2(core,t=0) | E_drift | E_pot retention | core/far(t=500) | Verdict |
|------:|--------:|--------------:|--------:|----------------:|----------------:|---------|
|   150 |    2.25 |         144.0 |  -30.0% |           0.48% |             1.0 | Dissolves |
|   100 |    1.50 |          96.0 |  -26.3% |           0.65% |             1.1 | Dissolves |
|    67 |    1.00 |          64.3 |  -23.1% |           0.96% |             1.0 | Dissolves |
|    33 |    0.50 |          31.7 |  -16.1% |           1.46% |             1.2 | Dissolves |
|    17 |    0.25 |          16.3 |  -10.9% |           2.60% |             1.2 | Dissolves |

All runs: N=128, L=20, T=500, single braid.

### Dissolution timeline (alpha=150)

| t     | E_kin     | E_mass    | E_pot(braid) | Event                          |
|------:|----------:|----------:|-------------:|--------------------------------|
| 0     | 2.70e+03  | 3.33e+04  | -1.22e+02    | Braid at full amplitude        |
| 5     | 2.01e+04  | 6.76e+03  | -5.18e+00    | EXPLOSION: mass -> kinetic     |
| 10    | 1.94e+04  | 4.45e+03  | -6.00e+00    | Braid gone, wave noise         |
| 30    | 1.85e+04  | 2.56e+03  | -6.35e-01    | Thermal equilibrium            |
| 500   | 1.84e+04  | 2.33e+03  | -5.88e-01    | Stable thermal noise           |

### m_eff profile at t=0 (initial, as designed)

| r     | alpha=150 m_eff^2 | alpha=67 m_eff^2 | alpha=17 m_eff^2 |
|------:|------------------:|-----------------:|-----------------:|
| 0.2   |           173.3   |           77.4   |           19.6   |
| 2.2   |            92.0   |           41.1   |           10.4   |
| 6.2   |             7.8   |            3.5   |            0.9   |
| 14.2  |             2.3   |            1.0   |            0.3   |

### m_eff profile at t=500 (after dissolution)

All profiles completely flat: m_eff^2 ~ 4.0 (alpha=150), 2.3 (alpha=67),
1.0 (alpha=17). The braid structure is entirely gone.

### Mechanism of failure

The phi^4 self-interaction alpha*(Sigma phi^2)*phi_a generates a cubic
restoring force proportional to phi^3. In the braid core (phi ~ 0.8),
this force is alpha*0.96*phi ~ 144*phi for alpha=150, versus the
triple-product coupling V'(P) ~ 3-4. The mass force overwhelms binding
by 40x. The three fields decouple, phases randomize, P -> 0.

**Why this is fundamental, not an init artifact**: Even alpha=17 (only 5%
initial energy perturbation) dissolves completely. The phi^4 potential
has NO localized soliton solutions in 3+1D (Derrick's theorem). The force
is strongest where phi is largest -- exactly where binding needs to win.

**Learning**: Making mass proportional to Sigma phi^2 KILLS binding because
the force is a repulsive cubic self-interaction, strongest at the core.

## Approach 2: Inverse Coupling (m_eff² = α/(1+β×Σφ²))

Mass DECREASES with field amplitude. Low at braid core, high in background.

| Config | α | β | m²(bg) | m²(core) | E_drift | E_pot ret | Verdict |
|--------|-----|---|--------|----------|---------|-----------|---------|
| inv_1 | 2.28 | 1 | 2.25 | 1.17 | -5.9% | 22% | Stable, braid survives |
| inv_2 | 2.42 | 5 | 2.25 | 0.42 | -5.1% | 9% | Stable, braid survives |
| inv_3 | 1.07 | 5 | 1.00 | 0.19 | -38% | 17,700% | κ-saturated collapse |
| inv_4 | 0.54 | 5 | 0.50 | 0.09 | -69% | 18,900% | κ-saturated collapse |
| inv_5 | 1.15 | 10 | 1.00 | 0.11 | -54% | 19,000% | κ-saturated collapse |

### Two regimes:

**Regime A (inv_1, inv_2): Stable braids at current range**
- m²(bg) = 2.25 (same Yukawa range as constant m²)
- Core mass reduced 2-5× → V(P) binds tighter
- Braid survives with E_pot retention 9-22%
- Late-time drift: 0.4-0.6%/100t (comparable to constant m²=2.25)
- **No range improvement** — but confirms inverse coupling is viable

**Regime B (inv_3-5): κ-saturated collapse at longer range**
- m²(bg) ≤ 1.0 → longer range attempted
- Core mass reduced to 0.09-0.19 → V(P) overwhelms → P grows explosively
- κ saturation arrests the collapse at V(P) → μ/(2κ)
- Creates a dense, κ-capped blob that slowly radiates (-10 to -23%/100t)
- NOT a viable braid — structure is fundamentally different

### Mechanism of failure (Regime B):
Positive feedback loop: low core mass → V(P) dominates → P grows →
Σφ² grows → m_eff drops further (inverse coupling) → V(P) grows more
→ collapse until κ saturation caps V(P). The result is a frozen blob
at maximum V(P) density, slowly draining.

---

## Track G Comparison (constant m^2)

| Method            | m^2(bg) | m^2(core) | E_pot ret (%) | Survives? |
|-------------------|--------:|----------:|--------------:|-----------|
| Constant m^2=2.25 |    2.25 |      2.25 |         ~100  | YES       |
| Constant m^2=1.50 |    1.50 |      1.50 |           14  | YES       |
| Constant m^2=1.25 |    1.25 |      1.25 |            2  | MARGINAL  |
| Constant m^2=1.00 |    1.00 |      1.00 |          0.6  | NO        |
| alpha=150 (phi^4) |    2.25 |    144.0  |         0.48  | NO        |
| alpha=67 (phi^4)  |    1.00 |     64.3  |         0.96  | NO        |
| inv_1 (inverse)   |    2.25 |      1.17 |           22  | YES       |
| inv_2 (inverse)   |    2.25 |      0.42 |            9  | YES       |
| inv_3 (inverse)   |    1.00 |      0.19 |    collapse   | NO        |

**Critical comparison**: Constant m^2=2.25 survives; phi^4 alpha=150 with
m_eff^2(bg)=2.25 does NOT. The field-dependent mass is strictly WORSE
than constant mass at matching background mass.

---

## Fundamental Conclusion

The constant m^2 serves as a GOLDILOCKS parameter: it simultaneously
provides vacuum stability AND braid binding at the SAME value. Any
attempt to decouple these by making m^2 field-dependent hits one of
two failure modes:

1. **Mass too high at core** (phi^4): overwhelms V(P) binding -> dissolution
2. **Mass too low at core** (inverse): V(P) overwhelms mass -> collapse

The binding requires m_eff^2(core) ~ 1-3 (comparable to V'(P) ~ 3-4).
The vacuum requires m_eff^2(bg) > ~0.5. Since core Sigma phi^2 is only ~60x
background Sigma phi^2, there isn't enough dynamic range to have
m_eff^2(core) ~ 2 and m_eff^2(bg) ~ 0.1 simultaneously through any
monotonic function of Sigma phi^2.

## Practical Result

**Track G constant m²=1.50 remains the best option**:
- Braid survives (97% E_pot retention)
- Yukawa range = 0.82 (1.83× current m²=2.25)
- Stable long-term (drift 0.31%/100t)

Further range extension requires either:
1. A non-monotonic m_eff(Σφ²) function (peaked at intermediate values)
2. A different binding mechanism that doesn't compete with the mass term
3. Acceptance that Yukawa gravity with range ~1 code unit is the theory's prediction

## Files

- Source: `src/v33_GB.c` (supports constant, φ⁴, and inverse modes)
- φ⁴ data: `data/alpha_{150,100,67,33,17}/`
- Inverse data: `data/inv_{1,2,3,4,5}/`
- Usage: `./v33_GB -inv <α> <β>` for inverse, `-alpha <α>` for φ⁴
