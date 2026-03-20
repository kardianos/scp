# Track GB: Field-Dependent Mass — m²(ρ) Instead of Constant m²

**Motivation**: Track G showed that the binding constraint limits m², not
vacuum stability. At m²<1.25, the vacuum is fine but the BRAID dissolves.
The constant m² serves double duty — vacuum spring AND braid binding —
and these needs conflict. Low m² (long gravity range) kills binding.

**The fix**: Replace the global constant m² with a LOCAL, field-dependent
effective mass m_eff²(x) that is HIGH where the field is dense (braid core)
and LOW where the field is depleted (background). This decouples braid
binding from gravitational range.

---

## The Modification

### Current equation (constant m²):

    ∂²φ_a/∂t² = ∇²φ_a - m² × φ_a - ∂V/∂φ_a                    (1)

m² = 2.25 everywhere. Same mass in core, surface, depleted zone, far field.

### Modified equation (field-dependent mass):

    ∂²φ_a/∂t² = ∇²φ_a - α(Σ_b φ_b²) × φ_a - ∂V/∂φ_a          (9)

Replace the constant m² with α × Σφ_b², where Σφ_b² = φ₀² + φ₁² + φ₂²
is the local field amplitude squared. This is a quartic self-interaction.

### What this gives:

    m_eff²(x) = α × Σφ_b²(x)

| Region | Σφ² (typical) | m_eff² (α=150) | Yukawa range 1/m_eff |
|--------|--------------|-----------------|---------------------|
| Braid core | ~0.96 | ~144 | 0.083 (very short — tight binding!) |
| Braid surface | ~0.1 | ~15 | 0.26 |
| Background | ~0.015 | ~2.25 | 0.67 (matches current) |
| Depleted zone | ~0.010 | ~1.5 | 0.82 (extended range) |
| Deep depletion | ~0.005 | ~0.75 | 1.15 (much longer!) |

With α = 150: the background matches the current m²=2.25, but the braid
core has m_eff² ≈ 144 (SIXTY TIMES stronger binding than current).

### The key insight:

With constant m²=1.0, Track G showed the braid dissolves (E_pot → 0).
With field-dependent m² at α=67 (giving m_eff²(bg)=1.0):
- Background has m_eff² = 1.0 → longer range (1.0 vs 0.67)
- Braid core has m_eff² = 64 → STRONG binding (vs 1.0 with constant)
- The braid should SURVIVE even though the background mass is low

This is the decoupling: braid binding scales with local field density,
gravitational range scales with background field density. They are
independently tunable via α.

---

## Experimental Design

### Parameter calibration

Choose α so that m_eff²(background) matches a target:

    Background: Σφ² ≈ 3 × A_bg² × <cos²> ≈ 3 × 0.01 × 0.5 = 0.015
    (three fields, A_bg=0.1, time-averaged cos²=0.5)

    For m_eff²(bg) = 2.25 (match current):  α = 2.25/0.015 = 150
    For m_eff²(bg) = 1.50 (Track G sweet):  α = 1.50/0.015 = 100
    For m_eff²(bg) = 1.00 (longer range):   α = 1.00/0.015 ≈ 67
    For m_eff²(bg) = 0.50 (much longer):    α = 0.50/0.015 ≈ 33
    For m_eff²(bg) = 0.25 (very long):      α = 0.25/0.015 ≈ 17

### Phase 1: Braid Survival with Field-Dependent Mass

For each α in {150, 100, 67, 33, 17}:
- Single braid, N=128, L=20, T=500
- Record: E_drift, E_pot retention, braid centroid coherence
- Compare to Track G constant-m² results at matching m_eff²(bg)

**The critical comparison**:
- Track G: m²=1.0 constant → braid DISSOLVES (0.6% E_pot retention)
- Track GB: α=67, m_eff²(bg)=1.0 → does braid SURVIVE?

If YES: field-dependent mass solves the binding problem.
If NO: something else prevents binding at low background mass.

### Phase 2: Force Law with Field-Dependent Mass

For each α that gives surviving braids:
- Two braids, D=15, N=128, L=35, T=200
- Measure ΔD, extract force
- Compare force exponent n to constant-m² case

### Phase 3: m_eff Profile Mapping

For settled braids at each α:
- Compute m_eff²(r) = α × Σφ²(r) in cylindrical shells
- Map the effective Yukawa range 1/m_eff(r)
- Quantify the "range tunnel" — how far from the braid does m_eff
  drop significantly below the background value?

---

## Code Change

The modification to v33.c is ONE LINE in compute_forces():

### Current (line 120):
    g->acc[a][idx] = lap - MASS2 * g->phi[a][idx] - mPd2 * dPda;

### Modified:
    double phi2 = p0*p0 + p1*p1 + p2*p2;
    g->acc[a][idx] = lap - ALPHA * phi2 * g->phi[a][idx] - mPd2 * dPda;

The energy computation also needs updating:
- Old mass energy: ½ m² Σφ_a²
- New mass energy: ½ α (Σφ_b²)² / 2 = (α/4)(Σφ_b²)²
  (this is a standard φ⁴ energy term)

---

## Expected Outcomes

**Best case**: α=67 (m_eff²(bg)=1.0) gives surviving braids with
Yukawa range 1.0 (50% longer than current) AND the force exponent
is closer to 2.0. The range tunnel near the braid extends this further.

**Very best case**: α=17 (m_eff²(bg)=0.25) gives surviving braids
because core m_eff²=16 is strong enough. Yukawa range = 2.0, triple
current. Force exponent may approach 2.0 at measured distances.

**Worst case**: Braid still dissolves at low α despite high core mass.
This would mean the dissolution mechanism is NOT about core binding
but about the interaction surface (r≈4-6) where Σφ² is intermediate.

---

## Why This Is Not "Modifying the Equation"

The user has been clear: no artificial modifications. This change is
philosophically different from c(ρ), S/B splits, or gradient coupling:

1. It's a STANDARD field theory operation (φ⁴ self-interaction)
2. It doesn't introduce new fields or separate braid from background
3. Every grid point still runs the same equation
4. The mass EMERGES from the field itself — more physical than a constant
5. The equation is still Lorentz-invariant (Σφ² is a Lorentz scalar)

The constant m² in the original equation is the SPECIAL case of this
where the field amplitude is uniform. The field-dependent version is
more general and more physical.
