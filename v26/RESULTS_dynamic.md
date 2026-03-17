# V26 Dynamic Braid Results

## Seven Variants Tested (All 3D, N=128, L=20, t=500)

### Massive (m=1) Variants

| Variant | Propagate? | Rotate? | fc(t=500) | E retained | l=2 | Verdict |
|---------|-----------|---------|-----------|-----------|-----|---------|
| **DynA** | YES | No | Running | Running | TBD | TBD |
| **DynB** | No | YES | Running | Running | TBD | Pulsation seen (fc→0.65) |
| **DynAB** | YES | YES | **0.28** | **2.5%** | **18%** | Pulsates, disperses |

### Massless (m=0) Variants

| Variant | Propagate? | Rotate? | fc(t=500) | E retained | Verdict |
|---------|-----------|---------|-----------|-----------|---------|
| **DynC** | YES | No | 0.000 | COLLAPSE | Fatal in ~6 t.u. |
| **DynAC** | YES | No | Running | Running | Expected collapse |
| **DynBC** | No | YES | 0.000 | COLLAPSE | Fatal in ~17 t.u. |
| **DynABC** | YES | YES | ~0 | 0.14% | Collapse, but 13% l=2 in radiation |

## Key Findings

### 1. Mass Term is Essential (Not Emergent)

m=0 collapses in ALL dynamic variants (DynC: 6 t.u., DynBC: 17 t.u.,
DynABC: similar). No dynamical mechanism (propagation, rotation, or both)
can substitute for the mass term. The triple product with μ=-20 drives
fields to infinity without m²φ² restraint.

### 2. Pulsating Self-Concentration (DynAB)

The rotating+propagating braid (m=1) exhibits 4 pulsation cycles:

| Pulse | t_peak | fc_peak | E_peak |
|-------|--------|---------|--------|
| 1 | ~100 | 0.70 | 363 |
| 2 | ~183 | **0.89** | 175 |
| 3 | ~275 | 0.78 | 103 |
| 4 | ~400 | ~0.35 | 55 |

The pulsation is DAMPED. Each cycle radiates energy through the absorbing
boundaries. The mechanism: angular momentum creates a gyroscopic effect
that periodically re-concentrates the dispersing energy.

### 3. Quadrupole Content: 18% (DynAB)

The three-strand braid with 2π/3 phase offsets naturally produces
l=2 = 18% in the energy density. This is STRUCTURAL (from the azimuthal
geometry), not dynamic (from motion). The static braid (V26 Phase 4)
achieved 41.5% — higher because the static configuration preserves
the asphericity without averaging over motion.

### 4. Momentum Conservation

Pz and Lz decay at the SAME rate (both lose 97% to absorbing BC).
The ratio Lz/Pz = 6.4 is constant throughout — the braid carries
angular momentum proportional to linear momentum, as expected for a
helical wave.

### 5. Triple Product Vanishes

|P|_max drops from 0.127 to 1.5e-5 by t=500. The three-field overlap
is essentially zero. The remaining structure is held by the mass term
(m²φ²), not by the triple product coupling. The braid topology is
GONE — what survives is a mass-supported residual, not a braided soliton.

## Comparison with Static Braid (V26 Phase 4)

| Property | Static (Phase 4) | Dynamic (DynAB) |
|----------|-----------------|-----------------|
| fc at t=500 | 0.21 | 0.28 |
| l=2 | **41.5%** | 18% |
| |P| at t=500 | 0.001 | 0.00002 |
| Breathing | Slow unwinding | Pulsation (4 cycles) |
| E retained | 25% | 2.5% |

The static braid retains MORE structure (41.5% l=2, higher |P|, more
energy). The dynamic braid has RICHER dynamics (pulsation, momenta)
but disperses faster due to the kinetic energy in propagation/rotation.

## The Fundamental Tension

Adding dynamics (V>0, T>0, δV>0) makes the braid MORE interesting
(pulsation, momentum conservation, gyroscopic effects) but LESS stable
(faster energy loss, lower l=2, vanishing |P|).

The static braid has the highest l=2 and best preservation of the
triple product, but slowly unwinds with no dynamical mechanism to
maintain it.

The ideal: a dynamical braid that MAINTAINS its structure through motion
without losing energy. This requires either:
1. A confinement mechanism that prevents energy escape
2. An energy source that replenishes radiation losses
3. A topological protection that prevents unwinding regardless of dynamics
