# Field Knot Theory

## Concept

A **field knot** is a 3D+time scalar field that forms self-reinforcing wave patterns along Celtic knot curves. Unlike solitons (localized particle-like solutions), field knots are extended structures - think of a violin string shaped into a knot, vibrating with standing waves.

**Key insight**: Knots can EMERGE from the field itself, not be placed parametrically.

## Mathematical Framework

### 1. Field and Density

The fundamental object is a complex scalar field Φ(x,t):

```
Φ: R³ × R → C

Density:  ρ(x,t) = |Φ|²
Phase:    φ(x,t) = arg(Φ)
Gradient: ∇ρ = 2·Re(Φ*·∇Φ)
```

### 2. Density Gradient Dynamics

Field evolution driven by density gradients:

```
∂Φ/∂t = -γ·∇ρ·(Φ/|Φ|)·sign(ρ-ρ_eq) + i·V[Φ] + D·∇²Φ + S[Φ]
        └────────────────────────────┘   └───┬───┘  └──┬──┘  └──┬──┘
              density pressure              potential  diffusion  source
```

### 3. Emergent Structure Condition

Knots emerge where:
- |Φ| passes through zero (nodal surface)
- ∇Φ has non-zero winding
- Phase φ has non-trivial topology

**Winding number** around loop C:
```
W(C) = (1/2π) ∮_C ∇φ · dl
```

W(C) ≠ 0 indicates topological structure.

### 4. Field Consumption and Stability

Knots **consume** field density:

```
Γ_consume = ∫_knot |∂ρ/∂t| dV
Γ_refill = ∫_boundary (∇ρ · n̂) dS

Stability: Γ_consume ≈ Γ_refill
```

### 5. S³ → R³ Mapping (Hopf Fibration)

Natural topological structure from Hopf map:

```
(z₁, z₂) ∈ S³ ⊂ C² → Hopf(z₁, z₂) ∈ R³

X = 2(x₁x₃ + x₂x₄)
Y = 2(x₂x₃ - x₁x₄)
Z = x₁² + x₂² - x₃² - x₄²
```

This creates linked structures with non-trivial topology.

### 6. Coupled Knot Dynamics

Cross-coupling force between knots:

```
F₁←₂[i] = Σⱼ A₂[j] · K(‖γ₁(i) - γ₂(j)‖) · cos(φ₁(i) - φ₂(j))
```

**Coupling energy**:
```
E_coupling = ∫∫ A₁A₂K₁₂(d)cos(Δφ) ds dt'
```

- E < 0: Binding (knots attract)
- E > 0: Repulsive (knots destabilize)

### 7. Boson-Like Transfers

"Bosons" are NOT particles but **transfer patterns**:

```
B(x,t) = ∂/∂t [∫ G(x,y)·ρ(y)·cos(Δφ) dy]
```

Where G(x,y) is a Yukawa-like propagator.

**Transfer rate**:
```
Rate = κ · ∫ ρ₁·ρ₂·cos(Δφ) dV
```

- In-phase (Δφ≈0): Maximum transfer
- Quadrature (Δφ≈π/2): Minimal net transfer
- Out-of-phase (Δφ≈π): Reverse/repulsive transfer

## Phase Relationships and Stability

| Phase Difference | Coupling Type | Stability |
|------------------|---------------|-----------|
| Δφ ≈ 0 (In-phase) | Binding | Generally stable |
| Δφ ≈ π/2 (Quadrature) | Mixed | Often unstable |
| Δφ ≈ π (Out-of-phase) | Repulsive | Can be stable at low coupling |

## Stability Map

```
              Phase Offset
              0    π/2    π
           ┌────┬────┬────┐
      0.0  │ ■  │ ■  │ ■  │  (All stable - no coupling)
Coup  0.3  │ ■  │ □  │ ■  │  (Quadrature unstable)
Str   0.6  │ ■  │ □  │ ■  │  (Stronger coupling)
      0.9  │ ■  │ □  │ □  │  (High coupling destabilizes)
           └────┴────┴────┘
           ■ = Stable  □ = Unstable
```

## Key Parameters

| Parameter | Symbol | Effect |
|-----------|--------|--------|
| Wave number | n | Standing wave mode (must fit around knot) |
| Spatial decay | σ | How localized field is around curve |
| Nonlinearity | α | Saturation strength |
| Self-coupling | β | Reinforcement strength |
| Damping | γ | Energy dissipation rate |
| Cross-coupling | κ | Strength of knot-knot interaction |
| Coupling range | λ | Characteristic distance for boson transfer |

## Locality vs Non-Locality

The system exhibits both local and non-local behavior:

| Component | Nature | Description |
|-----------|--------|-------------|
| Wave propagation c²∂²A/∂s² | LOCAL | Information at speed c |
| Self-reinforcement F_self | NON-LOCAL | Instantaneous sum over curve |
| Cross-coupling | NON-LOCAL | Instantaneous phase coherence |

This duality creates quantum-like behavior (discrete stability states) from continuous input.

## Simulations

1. **sim1_density_dynamics.py**: Field consumption and density gradients
2. **sim2_boson_transfer.py**: Energy transfer between coupled knots

## Files

- `THEORY.md` - This document
- `investigation_map.py` - Phase 1: Mathematical foundations
- `investigation_map2_bosons.py` - Phase 2: Boson transfer theory
- `sim1_density_dynamics.py` - Consumption simulation
- `sim2_boson_transfer.py` - Transfer simulation
- `coupling_analysis.py` - Knot coupling stability analysis
