# v63 — c_eff in the Cl(3,1) Sector (NEGATIVE, as anticipated)

**Date**: 2026-05-28
**Status**: second gating check for the holonomy idea. Follows `CEFF_TEST.md`
(which falsified the BLV `v_min/c = 2/3` lead in the *old* gravity sector). This
tests the **current** `Cl(3,1)` sector. **Result: no parameter-free `c_eff`.**

---

## 1. Question

The holonomy idea wants the phase to be `exp(i ∮A) = e^{2i/3}` with a
density-dependent meta-constant supplying `2/3 = c_eff·(weight)`. For that anchor
to be real, `c_eff` must be **parameter-free**. Is there such a `c_eff` (a speed /
metric magnitude) in the current `Cl(3,1)` gravity sector?

## 2. Context pulled (reproduced)

From `v60/gravity_recast/08_plebanski_action.{md,py}`, `05_findings.md`,
`v61/PROVEN_LEDGER.md`, `v61/CLOSEOUT.md`:

- The spacetime metric is **not fundamental** — it is the **induced** metric
  `g(B)` via the Urbantke formula, which is **cubic in the 2-form** `B` and
  densitized. Running `08_plebanski_action.py` reproduces `ĝ ∝ g` (proportional —
  recovered **up to an overall scale**), rel. error `0`/`7.9e-16`.
- The `ρ_grav` coupling normalization is **open** (`08` §5–6).
- The only magnitude on offer, `G_e = (21/16)α²¹`, is flagged **"curve-fitting,
  not derivation"** (`PROVEN_LEDGER.md` P-fit) and is `α`-dependent.

So every **magnitude** in the sector is open or fitted; what is pinned is the
**conformal (causal/angle) structure**, the 2 TT DOF, helicity `±2`, the signature.

(Aside: `PROVEN_LEDGER.md:94` already proves `phase_not_pi_rational`,
`koide_not_pi_rational`, and `cos(2/3)` transcendental — so the v62 no-go agrees
with existing v61 theorems.)

## 3. Computation (`conformal_check.py`, self-verifying)

- **Urbantke is homogeneous of degree 3**: `ĝ(λΣ) = λ³ ĝ(Σ)` to machine precision.
  An overall rescaling of `B` rescales `g`, so `B` (the matter/density 2-form)
  fixes only the **conformal class** of `g` — never an absolute length/time scale.
  A "speed `c_eff`" is an absolute-scale (magnitude) quantity ⇒ **not pinned**.
- **The conformally-invariant speed *ratio* `−g₀₀/g₁₁` is scale-free** for a fixed
  metric — but the metric depends on the matter configuration, so the ratio is a
  position/parameter-dependent **field**, not a universal constant. The BLV speed
  ratio is exactly this object, and `ceff_test.py` already showed it moves with
  the pion mass (`0.6702` vs `0.6252`). So **no universal `c_eff = 2/3`** even as
  a conformal invariant.

## 4. Verdict

**NEGATIVE.** Both gravity sectors now fail to provide a parameter-free `c_eff`:

| sector | why no parameter-free `c_eff` |
|---|---|
| BLV (old, v3) | a definite speed, but **parameter-dependent** (`0.6702`→`0.6252` with `m_π`) |
| `Cl(3,1)` (current, v60/61) | metric **conformal-class only** (Urbantke degree-3) + coupling magnitude **open/fitted** |

A density-dependent **speed** is a magnitude, and magnitudes are precisely what
the framework leaves unpinned.

## 5. Double-check (before / after)

**Before**: confirmed the relevant object — the Urbantke metric is the sector's
metric, it is cubic/densitized (read `08_plebanski_action.{md,py}`), and the
magnitude items are explicitly "open"/"fit" (read `05_findings`, `PROVEN_LEDGER`).
Reproduced `ĝ ∝ g` by running `08_plebanski_action.py`.

**After**: `conformal_check.py` 4/4 embedded checks pass — degree-3 homogeneity
(`λ³`), ratio conformal-invariance, BLV matter-dependence, and non-invariance of
the magnitude. Numbers consistent with `CEFF_TEST.md`.

## 6. Constructive conclusion

The "meta-const `c`" was the wrong object. Speeds are never pinned in the emergent
spacetime metric (conformal-class + matter-dependent). What *is* pinned
parameter-free and scale-free is the **internal `G₂`/`S⁷` geometry**: a canonical
round metric and (via triality) a canonical cycle, where a `curvature × (fractional
area)` **holonomy is a pure dimensionless number, independent of any scale**.

So if the holonomy/integration route is pursued, the phase is that **internal
geometric angle — with no `c_eff` factor at all**; the scale-freedom comes from
**conformal invariance / the fixed internal geometry**, not from a magic value of
`c`. That — computing the canonical `G₂`/`S⁷` holonomy on the triality cycle and
checking it against `2/3` with zero free parameters — is the real holonomy
experiment, and the only one worth building.

**Recommendation**: still do **not** build on a `c_eff`. The next (and genuinely
canonical) gating check is whether the internal `G₂`/`S⁷` triality holonomy is
parameter-free *and* lands on `2/3`. That is a pure-geometry computation, decoupled
from the dynamical gravity sector.

**Artifacts**: `conformal_check.py` (4/4 checks); evidence cross-referenced to
`CEFF_TEST.md` and `v60/gravity_recast/08_plebanski_action.py`.
