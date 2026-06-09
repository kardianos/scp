# The Flat-Direction Law

**Law.** Representation theory pins the *symmetric skeleton* of each sector and
is constitutionally blind to the *flat (angular) directions* it leaves behind.
Every residual conjecture of v59–v61 sits on such a flat direction. Same shape
each time: **magnitude fixed, angle flat.**

Proven on two independent sectors in `flat_direction_demo.py` (self-verifying).

---

## Instance 1 — the Brannen phase

The `Z₃`-symmetric mass ansatz `√mₖ = a(1 + 2t·cos(2πk/3 + φ))` is what the
algebra forces. Splitting its invariants by `φ`-dependence (symbolic, `sympy`):

| invariant | value | `∂/∂φ` |
|---|---|---|
| `Σ√mₖ` | `3a` | `0` (flat) |
| `e₂(√m)` | `3a²(1−t²)` | `0` (flat) |
| `Σmₖ` | `3a²(1+2t²)` | `0` (flat) |
| **Koide `Q`** | `(1+2t²)/3` | `0` (flat) |
| `∏√mₖ` (`e₃`) | … `+ 2t³cos(3φ)` | `−6a³t³ sin(3φ)` |
| `Σ(√m)³` (`P₃`) | … `cos(3φ)` | `−18a³t³ sin(3φ)` |

Every **symmetric** invariant — the ones rep theory computes — is
`φ`-**independent**; `φ` enters **only** the `cos(3φ)`-class. At the lepton point
`t² = 1/2`, the flat invariants already fix `Q = 2/3` exactly. So the algebra
determines the entire skeleton and leaves `φ` as a coordinate it cannot see.

Its value is **transcendental** (`cos(3φ) = cos(2/3)`, see `../no_go/NOGO.md`),
so only dynamics can lift it.

---

## Instance 2 — the EW democracy (R1)

The O(784) Frobenius hat `V(Y) = (λ/4)(‖Y‖²_F − v₀)²` depends only on the radius
`‖Y‖_F`. At any vacuum the Hessian has exactly **1 positive mode** (the radial
Higgs) and **`n−1` zero modes** (Goldstones) — verified numerically for a demo
`n`, and formal at `n = 784` in `v61/lean/EwVevHome.lean` (`S^783`, `goldstones =
783 > 0`). "Democracy" is one point on this flat sphere, not selected by the
symmetric potential.

Its flat coordinate ranges over **rational/algebraic** points, so selecting one
is a *degeneracy* problem — closable by state-selection
(`../residual_audit/democracy_selection.py`).

---

## Why this matters

The two residuals look unrelated (a lepton mass phase; an electroweak vacuum
alignment) but are the **same structural object**: a symmetry-invariant
functional fixes a magnitude and leaves an angle flat. This:

- explains *why* every residual is residual (rep theory is blind to flat angles);
- separates the closable residuals (algebraic flat coordinate → state-selection)
  from the barred ones (transcendental flat coordinate → dynamics only);
- gives a second, independent confirmation of the number-type thesis beyond the
  phase alone.
