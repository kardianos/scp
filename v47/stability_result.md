# V47 Stability Result

## Finding: η(P) coupling is unstable for η_max > ~3

Both linear η(P) = η₀ + η₁|P| and saturating η(P) = η₀ + η₁P²/(1+κP²)
are numerically unstable when the effective coupling exceeds ~3×η₀.

### Stability boundaries

| Form | Stable limit | η_max | Needed for binding |
|------|-------------|-------|-------------------|
| Linear η₀+η₁\|P\| | η₁ < 7 | 1.2 | η₁=115, η_max=12 |
| Saturating η₀+η₁P²/(1+κP²) | η₁ < 100 | 2.5 | η₁=375, η_max=8 |
| **Required for binding** | — | **~8** | — |

### Root cause

Positive feedback loop:
1. η(P) > η₀ → stronger curl force on φ
2. Stronger curl force → larger φ amplitude oscillation
3. Larger oscillation → larger |P| = |φ₀φ₁φ₂|
4. Larger |P| → larger η(P) → goto 1

The V(P) saturation and the P²/(1+κP²) saturation both limit P itself,
but they DON'T limit the curl FORCE, which grows proportionally to η_eff.
The curl force drives φ amplitude growth, which increases P, which
increases η_eff — the loop is fundamentally unstable above a threshold.

### Possible resolutions

1. **NEGATIVE feedback**: η(P) = η₀/(1 + η₁P²). Coupling DECREASES at
   braid cores instead of increasing. This is CONFINEMENT-like (coupling
   weakens inside, not strengthens). But this is the OPPOSITE of what we
   wanted for nuclear binding.

2. **Delayed feedback**: η depends on TIME-AVERAGED P, not instantaneous.
   This prevents the oscillation-amplification loop. Requires modifying
   the kernel to track ⟨P⟩ over a running window.

3. **Implicit coupling**: Instead of η(P) in the curl force, use η(P) only
   in the POTENTIAL V(P). This modifies the binding but not the curl force,
   breaking the feedback loop.

4. **Massive θ (m_θ > 0)**: This was the alternative from V46. A massive θ
   creates Yukawa attraction WITHOUT amplifying the curl force. It's stable
   because the mass term DAMPS θ growth instead of amplifying it.

5. **η₁ at the stable maximum**: Use η₁=100 (η_max=2.5) and accept the
   weaker coupling. It amplifies the curl force by 5× at the core, which
   may be enough for SOME binding even if below the η=8 threshold.

### Recommended path forward

**Option 5 first**: Run the V45-style separation sweep with η₁=100 on the
existing instance. If ANY separation shows lower E_total than D=80, binding
exists (even if weak). This is the fastest test.

If no binding at η₁=100, switch to **Option 4 (massive θ)** which is
theoretically cleaner and doesn't have the stability problem.
