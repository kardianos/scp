# V24-MB: Three Complex Fields with Quark Charges

## Background

The most physical Maxwell integration: make all three fields complex with
charges matching quark charges. The triple product P = Ψ₁Ψ₂Ψ₃ carries total
charge e₁+e₂+e₃. If charges are chosen as +2/3, +2/3, -1/3 (up, up, down
quarks), the baryon has charge +1 (proton).

The potential V(|P|²) is automatically gauge-invariant.

## Setup

Three complex scalar fields Ψ_a = φ_a + iχ_a (6 real DOF total) plus a
U(1) gauge field A_μ.

### Charges

    e₁ = +2/3    (up)
    e₂ = +2/3    (up)
    e₃ = -1/3    (down)

Baryon charge: e₁+e₂+e₃ = +1
Anti-baryon: -1

### Lagrangian

    L = Σ_a [|D_μ^a Ψ_a|² - m_a²|Ψ_a|²] - V(|P|²) - ¼F_{μν}F^{μν}

where D_μ^a = ∂_μ - ie_aA_μ and P = Ψ₁Ψ₂Ψ₃.

V(|P|²) = (μ/2)|P|²/(1+κ|P|²) — gauge invariant since |P|² transforms as
|e^{i(e₁+e₂+e₃)α}P|² = |P|² when the total charge is... wait, that's only
gauge invariant if the total charge is zero. With e_total = +1, P picks up
a phase e^{iα} under gauge transform.

**Fix**: Use V(|P|²) which IS gauge invariant (|P|² = |Ψ₁|²|Ψ₂|²|Ψ₃|²
times phase factors that cancel in the modulus squared).

Actually: P = Ψ₁Ψ₂Ψ₃ → e^{i(e₁+e₂+e₃)α}P = e^{iα}P under A_μ → A_μ + ∂_μα.
So |P|² → |P|² (phases cancel). V(|P|²) IS gauge invariant. ✓

### EOM (1D, temporal gauge A₀=0)

For each complex field (writing Ψ_a = φ_a + iχ_a):

    ∂²φ_a/∂t² = ∂²φ_a/∂x² + 2e_aA∂_xχ_a + e_a(∂_xA)χ_a
                 + e_a²A²φ_a - m_a²φ_a - Re(∂V/∂Ψ_a*)
    ∂²χ_a/∂t² = ∂²χ_a/∂x² - 2e_aA∂_xφ_a - e_a(∂_xA)φ_a
                 + e_a²A²χ_a - m_a²χ_a - Im(∂V/∂Ψ_a*)

    ∂²A/∂t² = ∂²A/∂x² - Σ_a e_a(φ_a∂_xχ_a - χ_a∂_xφ_a + e_aA|Ψ_a|²)

The coupling derivative: ∂V/∂Ψ_a* = (μ/2)·∂|P|²/∂Ψ_a* / (1+κ|P|²)²
where ∂|P|²/∂Ψ₁* = P·Ψ₂*Ψ₃* (complex conjugate of the cofactor).

### Simplification for 1D test

Start with all REAL initial conditions (χ_a = 0). Then:
- P = φ₁φ₂φ₃ (real, same as v21)
- Current j = 0 (no imaginary parts → no current initially)
- The gauge field is not sourced initially → A stays zero

The EM field only activates when the imaginary parts χ_a develop.
This happens through the gauge coupling (A≠0 mixes φ and χ).

To SEED the EM field: add a small perturbation to A or to χ_a.

## What to Compute

### Phase 1: Oscillon formation with e=0 (control)

1. Initialize: φ_a = A·g(x), χ_a = 0, A=0
2. Verify the |P|² coupling produces a standard oscillon (real fields only)
3. Confirm oscillon ω, E, dE/dt match v21 values

### Phase 2: Turn on gauge coupling

4. Scan e_scale ∈ {0.0, 0.1, 0.3, 1.0} where charges are e_a = e_scale × {2/3, 2/3, -1/3}
5. Add small χ perturbation to seed the EM sector
6. Evolve t=10000. Measure: oscillon survival, EM energy, charge Q

### Phase 3: UUD vs UDD with EM

7. UUD: m₁=m₂=m_U=1.0, m₃=m_D=0.95, charges (+2/3,+2/3,-1/3) → Q=+1
8. UDD: m₁=m_U=1.0, m₂=m₃=m_D=0.95, charges (+2/3,-1/3,-1/3) → Q=0
9. Compare masses: does EM self-energy correct the mass ordering?

## Reference Code

- v21 1D solver: `/home/d/code/scp/v21/src/triad1d.c`

## Output

- `src/maxwell_b.c` — solver (6 real scalar fields + gauge field = 7 fields)
- `data/` — TSV output
- `RESULTS.md`

## Parameters

μ=-20, κ=20, m=1.0
Charges: e_scale × (2/3, 2/3, -1/3)
e_scale scan: {0.0, 0.1, 0.3, 1.0}
Nx=4000, xmax=100, tfinal=10000

Compile: `gcc -O3 -Wall -o maxwell_b src/maxwell_b.c -lm`
