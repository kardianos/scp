# v61 Loop — Closeout

**Date**: 2026-05-26
**Scope**: extend the v60 dynamical Lagrangian to (1) nonlinear/curved-space gravity
backreaction and (2) a dynamical home for the EW vev `v=784a²` (R1) — the two next
steps from the v60 closeout — one new aspect per generation, each verified with ≥2
of {SymPy, Lean, Maxima, C}.
**Outcome**: **both stated goals met**, plus the LIGO motivation closed and all three
classic GR tests reproduced. Six generations, all verified and reproducible.

---

## The arc

| gen | aspect | key result | tools |
|---|---|---|---|
| 1 | curved gravity | nonlinear elimination → full GR; **Schwarzschild from `ρ_grav`**; weak-field = v60 Newtonian; deflection `4GM/b` | SymPy, Maxima ctensor, Lean |
| 2 | matter backreaction | `G^t_t=−2m'/r²` ⟹ **`m'=4πr²ρ`**; `M=∫ρ_grav`; interior↔exterior; Poisson limit | SymPy, Maxima ctensor, Lean |
| 3 | EW-vev home (R1) | **R1 home** = `End(L)` Frobenius Higgs; `v=784a²`; but O(784)-degenerate ⟹ democracy not selected (honest residual) | SymPy, Maxima, Lean |
| 4 | LIGO closure | 2 TT modes = **`h₊,h×` (helicity ±2)**; quadrupole `(32/5)Gμ²a⁴ω⁶`; GWs at `c` | SymPy, Maxima, Lean |
| 5 | perihelion | orbit GR term (coeff 3) → **`6πGM/(c²a(1−e²))`**; Mercury 42.98″/cy | SymPy, Maxima, Lean |
| 6 | synthesis | this closeout; regression **10/10**; 5 Lean modules clean | all |

---

## Goal 1 — nonlinear/curved gravity backreaction: DONE

The v60 gravity sector (linearized: TT gravitons + Newtonian OBE) is now a genuine
**curved-space general relativity**:
- the first-order action's connection elimination gives full GR nonlinearly (GEN1);
- the static spherical vacuum is **Schwarzschild**, sourced by `M = ∫ρ_grav = Σm = 9Qa²`
  (GEN1/2), with the v60 Newtonian `1/r` as its weak-field limit;
- the matter (GEN3 `T_00`) backreacts via `m'(r)=4πr²ρ` (GEN2);
- the theory passes **all three classic GR tests**: light deflection `4GM/b` (GEN1),
  gravitational waves (GEN4), perihelion precession `6πGM/(c²a(1−e²))`, Mercury
  43″/century (GEN5);
- the **LIGO motivation that began the G9 program is closed**: the 2 TT modes are the
  observed `h₊, h×` polarizations (helicity ±2), emitted via the quadrupole formula,
  traveling at `c` (GEN4).

## Goal 2 — the EW-vev home (R1): a home + a sharp residual

R1 (`v = ‖Y‖²_F = 784a²`) was the last residual without a dynamical origin. GEN3:
- **gives it a dynamical home**: an `End(L)` Frobenius Mexican-hat Higgs whose vacuum
  has `‖Y‖²_F = v₀`, realizing "v = Frobenius²"; equipartition reads `v=784a²`,
  `a=√v/dim(L)`; the `784=dim End(L)` is Burnside-forced;
- **but does not derive the value**: the O(784)-symmetric hat is degenerate (vacuum
  manifold `S^783`, 783 Goldstones), so the democratic vacuum is *not selected* —
  equipartition + the dimensionful identification remain value/symmetry conjectures
  (on the same footing as `α`). R1 goes from "homeless" to "home + isolated conjecture."

## Verification ledger

- **Python/SymPy**: 5 generation scripts, all assertions pass.
- **Maxima**: 5 cross-checks (2 via `ctensor`), all pass — *text-verified* (GEN5
  caught that `error()` exits 0 and `integrate` can hang on `assume` queries, so exit
  codes are unreliable; the harness now reads output text).
- **Lean**: 5 modules build clean against v59 Mathlib (`CurvedBackreaction`,
  `Backreaction`, `EwVevHome`, `GravitationalWaves`, `Perihelion`); genuine
  `field_simp`/`ring`/`positivity`/`Real.sqrt_sq` proofs, not just `decide`.
- **Regression**: `06_verify_all.py` → **10/10 pass**.

## Honest limitations

- The curved-gravity results use standard GR vacuum/weak-field solutions
  (Schwarzschild, mass function, geodesics) **sourced by the SCP charge `ρ_grav`**;
  they confirm the gravity sector *is* GR with the right source, not a new
  gravitational dynamics. The first-order action's elimination → GR is the v61
  content; the GR tests themselves are textbook (here tied to `ρ_grav = Σm = 9Qa²`).
- R1's value is **not derived** (the honest GEN3 result): the home exists but
  democracy/equipartition is an extra posit.
- No fully self-consistent boson-star interior (TOV with the GEN3 scalar) was solved
  numerically — GEN2 established the backreaction *law* and matching, not a numerical
  star.

## Residual conjectures (unchanged, isolated)

`α`; the EW-vev value `v=784a²` (now with a home, GEN3); the Brannen phase `φ=2/9`;
the gravity magnitude `f_g~α^{21/2}`. All are value/symmetry inputs, not dynamical gaps.

## Recommendation (v62)

- A numerical self-gravitating Koide/boson-star (TOV with the GEN3 potential),
  closing GEN2's interior solution quantitatively.
- The equipartition/democracy selection for R1 (what extra term picks the democratic
  `End(L)` vacuum) — or a proof it cannot be selected (sharper R1 negative).
- An FRW cosmology sourced by the matter sector (the one GR arena not yet touched).

The v61 loop concludes here: both goals met (curved gravity fully; R1 with a home +
sharp residual), the LIGO origin closed, the three GR tests reproduced, everything
verified and reproducible.
