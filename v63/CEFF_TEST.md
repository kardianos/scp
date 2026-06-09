# v63 — The c_eff Gating Test (NEGATIVE)

**Date**: 2026-05-28
**Status**: gating check for the "localization + density-dependent `c`" holonomy
idea. **Result: the proposed anchor is FALSIFIED.** The holonomy idea is not
killed, but it gets no parameter-free foothold from the existing density-dependent
`c`, so a holonomy v63 build is **not** yet warranted.

---

## 1. What was being tested

The inverted/integration picture (conversation following `../v62/`) would make the
Brannen phase a **holonomy**,

```
    e^{iφ·3} = e^{i·2/3} = exp(i ∮ A),     ∮A = 2/3 = c_eff · (canonical weight),
```

where a **density-dependent meta-constant `c_eff`** converts an algebraic weight
into the transcendental angle. This evades the v62 no-go (an integral, not a
finite algebraic operation), *provided* `c_eff`, the weight, and the cycle are all
pinned with **zero free parameters**. Otherwise `2/3 = c_eff·weight` is a dial and
we are back to numerology.

The tantalizing lead: the project's BLV null-rotor effective metric reports a core
velocity ratio `v_min/c ≈ 0.670 ≈ 2/3`. **Gating question:** is `c_eff = v_min/c`
(a) structurally `2/3` and (b) parameter-free (pinned by the field equation)?

---

## 2. Context pulled (the real derivation)

Source: `v3/tasks/done/gravity-null-rotor-metric.md` + `v3/src/nullrotor_metric.c`.

The BLV (Barceló–Liberati–Visser) phase velocity for a null-rotor polarization
`e_A` propagating along `j` on the hedgehog background is

```
    v_A(r)² / c² = (1 + P₄ʲʲ(r)) / (1 + m₄(r)),
    P₄ʲʲ = (c₄/2ρ₀⁴) Σ_{k≠j}|[e_A q̃₀, R_k]|²,   m₄ = (c₄/2ρ₀⁴) Σ_k|[e_A q̃₀, R_k]|².
```

`v_min/c` is computed as a **numerical minimum** of `v_A(r)/c` over radius `r` and
modes `A` (`nullrotor_metric.c` lines 479, 514–518: `if (vz[A] < vmin_all)
vmin_all = vz[A]`). It is a functional of the profile `f(r)`, **not** a closed-form
constant. The sector is explicitly a **nuclear-scale** (~250 MeV), Yukawa-range
effect with **no `h = ±2` tensor modes** — a *negative* result for gravity,
superseded by the v60/v61 `Cl(3,1)` soldering carrier.

---

## 3. Computation (reproduced this session)

Re-ran `v3/bin/nullrotor_metric` on both profiles (`ceff_test.py` re-runs the
binary and asserts agreement with these recorded values):

| soliton | `v_min/c` | gap to `2/3` |
|---|---|---|
| massless (`profile_sigma_e1.dat`) | **0.67020151** | 0.00353 (**0.53 %**) |
| massive (`profile_massive_e1_mpi0.398.dat`) | **0.62519338** | 0.04147 (**6.22 %**) |

---

## 4. Verdict — falsified on three counts

1. **Not equal to `2/3`.** Massless `v_min/c = 0.6702` misses `2/3` by 0.53 %; it
   is a numerical profile minimum with no structural reason to be `2/3`.
2. **Not density-independent.** Massless `0.6702` ≠ massive `0.6252` (a 6.7 %
   shift driven by the pion mass). A structural, parameter-free constant cannot
   move with `m_π`, so `v_min/c` is **not** a universal `c_eff = 2/3`.
3. **No consistent weight.** If `2/3 = c_eff·w`, the implied weight is `0.9947`
   (massless) vs `1.0663` (massive) — a 7 % disagreement, so there is no canonical
   weight either.

Plus the **sector caveat**: even a clean `2/3` here would be from the superseded,
nuclear-scale, no-tensor BLV metric — not the current gravity carrier.

The one **exact** structural constant in this sector is `P/m = 2` (`K=0`, machine
precision), not `2/3`.

---

## 5. Double-check (before and after)

**Before** — verified the context is the right object:
- `v_min/c` is a numerical `min` over the profile, not a closed form (read the
  loop in `nullrotor_metric.c`).
- The two distinct `0.670` values in the repo are unrelated: `v_min/c = 0.670`
  (this metric) vs `Δ = 3.02/4.51 = 0.670` (a skyrmion+hopfion anisotropy in
  `v2/.../tensor_prop.c`) — did not conflate them.

**After** — verified the conclusion:
- Numerical: `ceff_test.py` re-runs the binary, matches recorded values, and all
  six embedded checks pass.
- Formal: `lean/CeffNotStructural.lean` (no `sorry`) proves
  `vmin_massless ≠ vmin_massive` (density-dependent), `|vmin_massless − 2/3| >
  3/1000`, and the modus tollens `no_structural_ceff` / `two_thirds_not_structural_ceff`.

---

## 6. Consequence for the localization/holonomy program

The holonomy idea (v62 conversation) is **not refuted** — integration genuinely
escapes the no-go, and a holonomy *can* be a transcendental angle. What is refuted
is the specific hope that the **existing** density-dependent `c` hands us a
parameter-free `2/3`. It does not.

To proceed, `c_eff` must come from a genuinely **structural / virial-pinned**
quantity in the **current** (`Cl(3,1)`) gravity sector — or from a dimensionless
invariant of the `G₂`/`S⁷` geometry — with no free scale. Until such a `c_eff`
exists, `2/3 = c_eff·weight` remains a dial. **Recommendation: do not build a
holonomy v63 yet.** First find (or prove absent) a parameter-free `c_eff` in the
`Cl(3,1)` sector; that is the next gating check.

**Artifacts**: `ceff_test.py` (re-runs the binary, 6/6 checks),
`lean/CeffNotStructural.lean` (machine-checked, no `sorry`).
