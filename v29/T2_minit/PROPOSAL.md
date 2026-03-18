# T2: Initialization Mass Independence

## Question
Does the bimodal synergy survive with m_init = 0 (ω = k)?

## Motivation
Gemini's strongest critique: m_init=1.5 injects excess energy into
the initial velocity field, contradicting the "emergent mass" claim.
If the braid works with purely massless dispersion, the critique falls.

## Method
1. Run bimodal t=0.85 with m_init ∈ {0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0}
2. The Lagrangian mass is ALWAYS zero (m²=0 in EOM)
3. m_init ONLY affects initial velocity: ω = √(k² + m²_init)
4. N=128, T=500, absorbing xy BC, periodic z
5. Compare all metrics: trans_l2, torsion, fc, |P|, winding, energy

## Key Observable
- If m_init=0 gives trans_l2 > 0.1 AND torsion > 0.5: claim is solid
- If metrics degrade monotonically with m_init: initialization matters
- If there's a sweet spot at intermediate m_init: physics is robust but
  the specific numbers depend on how much energy is injected

## Also test
- m_init=0 with Lagrangian mass m²=0 (the "pure" massless case)
- m_init=0 with Lagrangian mass m²=1.5² (mass in dynamics but not init)
  This separates init-mass from dynamics-mass effects.

## Grid & Runtime
- N=128, T=500: ~8 min per eval with OMP
- 7 m_init values + 2 extras = 9 runs, ~72 min total
- Can run sequentially (each uses 16 cores)

## Build
```
cd v29/T2_minit && gcc -O3 -fopenmp -o t2 src/t2.c -lm
./t2
```
