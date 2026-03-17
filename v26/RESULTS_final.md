# V26 Final Results: Braided Solitons

## Summary

V26 introduced braided solitons — three displacement fields twisted into
knot-theory-inspired configurations (trefoil braid, Borromean rings) —
as an alternative to the breathing oscillon. The braid naturally breaks
spherical symmetry, producing 18-41% quadrupole (l=2) content.

## Phase 1-3 Results (Static Braids)

| Config | fc | |P| | Breathing? | l=2 | Status |
|--------|-----|------|-----------|-----|--------|
| Twisted tube + full | 0.07 | 0.001 | No (dispersed) | 4.7% | Disperses |
| Borromean + full | 0.92 | 0.39 | Yes (ω=1.37) | 1.6% | Survives but breathes |
| Twisted tube + triple only | 0.37 | 0.65 | **No (static)** | **6.8%** | Non-breathing |
| Oscillon control | ~1.0 | ~0.5 | Yes (ω=1.37) | 0.07% | Standard |

## Phase 4 Results (Torsion and Strain Couplings)

**41.5% l=2** from the baseline (m=1, no extra couplings) — highest ever.
Adding torsion κ_T or strain κ_S DESTROYS the braid by penalizing the
very gradients that define the twist. The simplest Lagrangian is best.

Self-consistent metric g^{ij} = δ^{ij} - 2ε_{ij} is STABLE through t=500.

## Dynamic Braid Results (V>0, T>0, δV>0)

### Completed

| Variant | m | Propagate | Rotate | fc(500) | l=2 | Verdict |
|---------|---|-----------|--------|---------|-----|---------|
| **DynAB** | 1.0 | YES | YES | 0.28 | **18%** | Pulsates (fc→0.89), damped |
| DynC | 0.0 | YES | No | 0.000 | — | Collapse in 6 t.u. |
| DynAC | 0.0 | YES | No | 0.000 | — | Collapse in 17 t.u. |
| DynABC | 0.0 | YES | YES | ~0 | 13% | Collapse, l=2 in radiation |

### Pending (still computing)

| Variant | m | Propagate | Rotate | Early signal |
|---------|---|-----------|--------|-------------|
| DynA | 1.0 | YES | No | TBD |
| DynB | 1.0 | No | YES | fc recovery to 0.65 |
| DynBC | 0.0 | No | YES | Collapse confirmed |
| DynABC | 0.0 | YES | YES | Running at N=64 |

## Key Discoveries

### 1. Highest Spin-2 Content: 41.5%

The static braided soliton (m=1, triple product only, no extra couplings)
produces 41.5% quadrupole on a spherical shell — 600× more than the
oscillon's 0.07%. This is STRUCTURAL (from the braid geometry), not
dynamic. The three-strand 2π/3 phase offset naturally creates l=2.

### 2. Pulsating Dynamic Braid

DynAB (propagating + rotating, m=1) shows damped pulsation:
- fc peaks: 0.70 → 0.89 → 0.78 → 0.35
- Period: ~100 time units
- Energy: 97.5% radiated, 2.5% persistent core
- Momenta: Lz/Pz = 6.4 (constant ratio, consistent with helical wave)

### 3. Mass Term is Necessary

m=0 collapses in ALL variants (4/4 tested). The triple product μ<0
drives fields to infinity without m²φ² restraint. No dynamical mechanism
(propagation, rotation, angular momentum) can substitute.

### 4. Torsion/Strain Couplings Hurt

Adding κ_T or κ_S UNWINDS the braid (l=2 drops from 41.5% to <5%).
The torsion term penalizes the antisymmetric gradients that define the
twist — it fights the very structure it should support.

### 5. Self-Consistent Metric Works

The strain-derived metric g^{ij} is stable for the braided soliton.
ds = c·dt is automatic from field propagation. The framework for
emergent GR is sound — what's missing is a STABLE braid to inhabit it.

## The Remaining Gap

The braided soliton produces the physics we want:
- 41% spin-2 content (gravity)
- Aspherical, non-breathing (for some configurations)
- Self-consistent metric
- Three forces from gradient decomposition

But it DISPERSES. The braid unwinds because:
1. The triple product |P| decays (three-field overlap decreases)
2. The mass term preserves energy but not topology
3. No mechanism enforces the braid crossing structure dynamically
4. Energy escapes through absorbing boundaries

The missing ingredient: a Lagrangian term that ENFORCES the braid
topology — making unwinding cost INFINITE energy (confinement) or
making the braid state a topological MINIMUM (like π₃(S³)=Z for
the Skyrmion, but using B₃ algebra instead of SU(2) homotopy).

## Equations

The V26 Lagrangian (same as v21, braided initial condition):

    L = Σ_a [½(∂_t φ_a)² - ½(∂_i φ_a)²] - ½m²Σφ_a² - (μ/2)P²/(1+κP²)

The strain (gravity): ε_{ij} = ½(∂_iφ_j + ∂_jφ_i)
The torsion (EM): ω_{ij} = ½(∂_iφ_j - ∂_jφ_i)
The metric: g_{ij} = δ_{ij} + 2ε_{ij}

Three forces from one gradient tensor ∂_iφ_j decomposed by symmetry.
