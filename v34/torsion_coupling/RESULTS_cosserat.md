# 6-Field Cosserat Simulation Results

## The Equation

3 position fields (φ) + 3 angle fields (θ), coupled through curl:

    ∂²φ_a/∂t² = ∇²φ_a - m²φ_a - V'(P) + η×curl(θ)_a
    ∂²θ_a/∂t² = ∇²θ_a - m_θ²θ_a       + η×curl(φ)_a

Parameters: m²=2.25, μ=-41.345, κ=50, A_bg=0.1
All runs: N=128, L=20, T=300, single braid, snapshots every T=50

---

## Phase 1: η Coupling Strength Scan (m_θ² = 2.25 = m²)

| η | drift | θ_rms | <E_pot> | Braid? |
|---|-------|-------|---------|--------|
| 0 | -0.28% | 0.0000 | -387 | YES (control, θ never excited) |
| 0.1 | -0.23% | 0.066 | -231 | YES |
| 0.5 | -0.25% | 0.066 | -143 | YES |
| 1.0 | -0.32% | 0.064 | -63 | YES (weaker binding) |
| 2.0 | -0.02% | 0.065 | ~0 | NO (braid dissolved) |
| 5.0 | blowup | 3.5e28 | — | EXPLODED |

**Key findings**:
1. θ fields ARE excited by the braid's helical curl (θ_rms ~ 0.065)
2. Braid survives at η=0.1 to 1.0, dissolved at η=2.0, explodes at η=5.0
3. Binding weakens monotonically with η (E_pot: -387 → -231 → -143 → -63 → 0)
4. θ_rms is remarkably constant across viable η values (~0.065)
5. Sweet spot: η = 0.1 to 0.5 (strong θ excitation, minimal braid weakening)

## Phase 2: Angle Field Mass Scan (η = 1.0)

| m_θ² | drift | θ_rms | E_pot(T=300) | Braid? |
|------|-------|-------|-------------|--------|
| 2.25 (= m²) | -0.32% | 0.064 | -63 | YES |
| 0.50 | -0.58% | 0.057 | -66 | YES |
| **0.00 (massless)** | **-0.60%** | **0.062** | **-44** | **YES** |

### CRITICAL RESULT: m_θ = 0 IS STABLE

The massless angle field (m_θ²=0) is stable for the full T=300 run:
- Braid survives with E_pot = -44 (weaker than massive case but alive)
- θ_rms = 0.062 (excited to similar level as massive case)
- Energy drift = -0.60% (slightly more than massive, but NOT growing)
- No blowup, no collapse, no runaway

**Why m_θ=0 works where m_φ=0 fails**: The φ fields have the attractive
V(P) coupling that drives tachyonic collapse at m=0. The θ fields have
NO V(P) — they're free Klein-Gordon fields sourced by curl(φ). Without
the attractive nonlinearity, there's no collapse mechanism.

## Implications

### We now have two sectors:

**Position sector (φ)**: Massive (m²=2.25), nonlinear (V(P)),
self-interacting. Creates braids. Mediates gravity through the
nonlinear depletion mechanism (power-law 1/r^1.2).

**Angle sector (θ)**: Can be MASSLESS (m_θ²=0), linear (no V(θ)),
sourced by curl(φ). Propagates at c without Yukawa decay.
Potential carrier of long-range EM-like forces.

### What makes this different from the Phase 2 torsion test:

In Phase 2, we tried to extract EM from perturbations of the existing
3 φ fields. All modes were massive (m_eff ≈ 1.5) → torsion couldn't
propagate. The failure was fundamental: you can't get a massless
mediator from a massive field.

In the Cosserat approach, the θ fields are GENUINELY NEW degrees of
freedom with their OWN mass parameter. Setting m_θ=0 makes them
massless by construction. The curl coupling to φ ensures they're
sourced by the braid's topology. This is not extracting a hidden
mode — it's adding the rotational degrees of freedom that an elastic
solid should have.

### Next steps:

1. **θ radial profile**: Does the massless θ field create a 1/r²
   pattern around the braid? (Massless field from a localized source
   should give exactly 1/r² at long range.)

2. **Two-braid θ force**: Do two braids interact through the θ field?
   Is the θ-mediated force different from the φ-mediated force?

3. **Opposite winding**: Do W=+1 and W=-1 braids create opposite
   θ patterns? Does the θ-mediated force reverse sign for opposite
   windings? (This would make winding = charge.)

4. **Moving braid**: Does a moving braid create a circular θ pattern
   (magnetic field analog)?

## Files

- `src/v33_cosserat.c` — 6-field simulation (18 arrays, curl coupling)
- `data/cosserat_eta_{0,0.1,0.5,1.0,2.0,5.0}/` — η scan results
- `data/cosserat_mt{0.5,0.0}/` — m_θ scan results
- Snapshot format: N(int), L(dbl), t(dbl), nf=6(int), φ[3][N³], θ[3][N³]
