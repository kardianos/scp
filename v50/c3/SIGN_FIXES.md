# V50/C3 Sign Fixes — Maxima-Verified Force Derivations

## The Two New Terms

V50/C3 adds two terms to the V44 Lagrangian:

    L_chiral = κ_h × P² × φ·curl(φ)        (chiral helicity bias)
    L_strain = -α × |curl(φ)/2 - θ|²       (Cosserat strain energy)

Both contribute forces to the phi equation. The theta equation gets
a force only from L_strain.

---

## Fix 1: Cosserat Strain — Phi Force Sign

### The energy

    E_strain = α Σ_a (curl(φ)_a/2 - θ_a)²

### BEFORE (wrong):

    cosserat_phi = +α × curl(M)_a

where M_a = curl(φ)_a/2 - θ_a.

### AFTER (correct):

    cosserat_phi = -α × curl(M)_a

### Why: Maxima derivation (derive_cosserat.mac)

The EL equation gives F_φ_a = ∂_j(∂E/∂(∂_j φ_a)).

For a=0, Maxima verifies:
    ∂E/∂(∂_y φ_0) = -α M_2     (confirmed: expand gives 0)
    ∂E/∂(∂_z φ_0) = +α M_1     (confirmed: expand gives 0)
    ∂E/∂(∂_x φ_0) = 0          (confirmed)

Therefore:
    F_φ_0 = ∂_y(-α M_2) + ∂_z(+α M_1)
           = α(∂_z M_1 - ∂_y M_2)
           = -α(∂_y M_2 - ∂_z M_1)
           = -α × curl(M)_0

The sign is MINUS because the EL divergence term is SUBTRACTED from
the algebraic term (which is zero here since E_strain doesn't depend
on φ_a directly, only on its spatial derivatives).

### Theta force (unchanged, was already correct):

    F_θ_a = -∂E/∂θ_a = +2α M_a = 2α(curl(φ)_a/2 - θ_a)

---

## Fix 2: Chiral Helicity — t3 (∇P² cross) Term Sign

### The Lagrangian

    L_h = κ_h × P² × (φ_0 curl_0 + φ_1 curl_1 + φ_2 curl_2)

where curl_a = curl(φ)_a and P = φ_0 φ_1 φ_2.

### The full EL force

    F_φ_a = ∂L/∂φ_a - ∂_j(∂L/∂(∂_j φ_a))

Maxima verifies:

    ∂L/∂φ_a = κ_h × (2P × dP/dφ_a × h + P² × curl(φ)_a)

    (confirmed: expand(dL/dp0 - expected) = 0 for all three components)

This gives terms t1 = 2P × dPda × h and part of t2 = P² × curl_a.

For the divergence part, Maxima verifies for a=0:

    ∂L/∂(∂_x φ_0) = 0                     (confirmed)
    ∂L/∂(∂_y φ_0) = -κ_h × P² × φ_2      (confirmed: diff = 0)
    ∂L/∂(∂_z φ_0) = +κ_h × P² × φ_1      (confirmed: diff = 0)

Therefore the divergence term is:

    ∂_j(∂L/∂(∂_j φ_0)) = ∂_y(-κ_h P² φ_2) + ∂_z(κ_h P² φ_1)
                        = κ_h × [-∂_y(P²)φ_2 - P²∂_y(φ_2) + ∂_z(P²)φ_1 + P²∂_z(φ_1)]
                        = κ_h × [-P² × curl(φ)_0 + ∂_z(P²)φ_1 - ∂_y(P²)φ_2]

The FULL force is F = ∂L/∂φ_0 - divergence:

    F_0 = κ_h × (2P dPda h + P² curl_0) - κ_h × (-P² curl_0 + dP²dz φ_1 - dP²dy φ_2)
        = κ_h × (2P dPda h + 2P² curl_0 - dP²dz φ_1 + dP²dy φ_2)

### BEFORE (wrong):

```c
if (a == 0)      t3 = dP2dz * p1 - dP2dy * p2;    // +dP²dz φ₁ - dP²dy φ₂
else if (a == 1) t3 = dP2dx * p2 - dP2dz * p0;    // +dP²dx φ₂ - dP²dz φ₀
else             t3 = dP2dy * p0 - dP2dx * p1;    // +dP²dy φ₀ - dP²dx φ₁
```

### AFTER (correct):

```c
if (a == 0)      t3 = -(dP2dz * p1 - dP2dy * p2);  // -dP²dz φ₁ + dP²dy φ₂
else if (a == 1) t3 = -(dP2dx * p2 - dP2dz * p0);  // -dP²dx φ₂ + dP²dz φ₀
else             t3 = -(dP2dy * p0 - dP2dx * p1);  // -dP²dy φ₀ + dP²dx φ₁
```

### Why:

The divergence term ∂_j(∂L/∂(∂_j φ_a)) is SUBTRACTED in the EL equation.
The code was adding the ∇P² cross terms with the WRONG SIGN because
the subtraction was not applied to the cross-product contributions.

For a=0: the correct force has (-dP²dz φ₁ + dP²dy φ₂), but the code
had (+dP²dz φ₁ - dP²dy φ₂). Same pattern for a=1 and a=2.

---

## Full Corrected Equations

### Phi equation:

    ∂²φ_a/∂t² = ∇²φ_a - m²φ_a - dV/dP × dP/dφ_a

                 + η × curl(θ)_a                              [V44 curl coupling]

                 + κ_h × (2P dPda h + 2P² curl(φ)_a + t3)    [chiral helicity]
                   where t3_0 = -(dP²dz φ₁ - dP²dy φ₂)       (and cyclic)

                 - α × curl(M)_a                               [Cosserat strain]
                   where M = curl(φ)/2 - θ

### Theta equation:

    ∂²θ_a/∂t² = ∇²θ_a - m_θ²θ_a

                 + η × curl(φ)_a                              [V44 curl coupling]

                 + 2α × M_a                                   [Cosserat strain]
                   where M = curl(φ)/2 - θ

### Verification status:

| Term | Maxima verified | Lean formalized |
|------|----------------|-----------------|
| V44: -dV/dP × dPda | Pre-existing | V44/Equations.lean |
| V44: +η curl(θ) | Pre-existing | V44/CurlCoupling.lean |
| V44: +η curl(φ) | Pre-existing | V44/CurlCoupling.lean |
| Cosserat: -α curl(M) | derive_cosserat.mac ✓ | V50C3/Cosserat.lean |
| Cosserat: +2α M | derive_cosserat.mac ✓ | V50C3/Cosserat.lean |
| Chiral: t1 = 2P dPda h | derive_chiral.mac ✓ | V50C3/ChiralHelicity.lean |
| Chiral: t2 = 2P² curl(φ) | derive_chiral.mac ✓ | V50C3/ChiralHelicity.lean |
| Chiral: t3 = -(∇P² cross) | derive_chiral.mac ✓ | TODO: update Lean |

### Note on instability at κ_h=1.0 after fix:

The corrected t3 sign caused NaN at κ_h=1.0 (t≈36). This may indicate
that the CORRECT force is actually destabilizing at large κ_h — the
chiral term with the right sign may be too strong. The wrong sign was
accidentally providing a stabilizing effect. Need to investigate:
- Does κ_h=0.1 survive?
- Is there a CFL-like stability condition on κ_h relative to dt?
- Is the chiral term fundamentally unstable, or just needs smaller dt?
