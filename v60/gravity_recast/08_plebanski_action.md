# The Plebański action for the SCP gravity sector + the induced-metric mechanism

**Date**: 2026-05-25 (the open item left by `07_findings.md`)
**Artifacts**: `08_plebanski_action.py` (induced-metric checks, all pass to ~1e-16),
`../lean/G9Plebanski.lean` (self-dual algebra, axiom-clean), this doc.
**Builds on**: `04`–`07` (soldering → 2 DOF → carrier forced to `Cl(3,1)` → unified algebra).

## Scope (honest, up front)

The action is **posited** Plebański-style — it is the natural 2-form gravity
action whose carrier is the `Cl(3,1)` self-dual bivector 2-form `B` that C3 forced
the graviton into, sourced by the internal second moment `ρ_grav`. It is **not
derived from the OBE** here. What is *verified* (numerically / machine-checked):
the simplicity constraint forces `B` to define a tetrad; the Urbantke formula
reconstructs the metric from `B` (the induced metric); the linearized content is
the 2 TT DOF and the trace law of `05`. What remains *posited / cited*: the full
nonlinear Einstein equation from varying `B` (standard Plebański), and the
specific normalization of the `ρ_grav` coupling.

## 1. Fields and action

Non-chiral (real `so(3,1)`) Plebański / BF form, matching the 6-field Cosserat
identification (`07`: the 6 `Cl(3,1)` bivectors `σ^{IJ}` = `3φ` rotations + `3θ`
boosts):

- `B = B^{IJ}_{μν}` — a 2-form valued in the `Cl(3,1)` bivectors `so(3,1)` (the
  soldering carrier; its self-dual half is the triple `B^i` used in `08_*.py`).
- `ω = ω^{IJ}_μ` — the `so(3,1)` (Lorentz/spin) connection; `F^{IJ} = dω^{IJ} +
  ω^{I}{}_K ∧ ω^{KJ}` its curvature.
- `Ψ_{IJKL}` — a Lagrange multiplier with Weyl symmetries (symmetric-traceless),
  enforcing simplicity.

```
S[B, ω, Ψ] = ∫_M  B^{IJ} ∧ F_{IJ}[ω]
                 − ½ Ψ_{IJKL} B^{IJ} ∧ B^{KL}              (gravity)
           +  S_source[ ρ_grav ; g(B) ]                    (matter trace, from Cl(7)_even)
```

The metric is **not fundamental** — it is `g(B)`, the induced metric (Urbantke, §3).
`ρ_grav = Tr(M†M) = Σ m_k` (the internal second moment, a Lorentz scalar) enters
`S_source` as the trace/Newtonian source, exactly the v59 coupling (`05` C2).

## 2. Euler–Lagrange equations (standard Plebański)

- **vary `Ψ`** → simplicity: `B^{IJ} ∧ B^{KL} ∝ ε^{IJKL}`. This forces `B^{IJ} =
  e^I ∧ e^J` for a tetrad `e^I_μ` ⇒ a metric `g_{μν} = η_{IJ} e^I_μ e^J_ν` emerges.
  **(verified in `08_*.py`: tetrad 2-forms satisfy it, generic ones don't.)**
- **vary `ω`** → `d_ω B^{IJ} = 0` (torsion-free): determines `ω` as the spin
  connection of `e`.
- **vary `B`** → `F^{IJ} = Ψ^{IJ}{}_{KL} B^{KL}` (+ source): Einstein's equations,
  with `Ψ` ↔ the Weyl tensor and the trace/source part fixed by `ρ_grav`.

So on-shell: `B` ⇒ tetrad ⇒ metric; `ω` ⇒ its connection; the last equation is GR
sourced by `ρ_grav`. The propagating content is the 2 TT graviton modes (`05` C1).

## 3. The induced metric (verified, `08_plebanski_action.py`)

The literal "metric emerges from the 2-form," to machine precision:

- **Simplicity does real work.** Self-dual 2-forms `Σ^i = η^i_{ab} e^a e^b` of a
  tetrad satisfy `Σ^i ∧ Σ^j ∝ δ^{ij}` (off-diagonal + trace-free residual ~1e-16);
  a *generic* 2-form triple violates it (residual ~1). So the constraint is what
  forces `B` into the tetrad/metric form.
- **Urbantke reconstruction.** `ĝ_{μν} = −(1/6) ε^{αβγδ} ε_{ijk} Σ^i_{μα} Σ^j_{βγ}
  Σ^k_{δν}` recovers the metric `g` we started from — relative error ~1e-16 for
  both flat and a random perturbed metric. The metric is a *function of `B`*.
- The three `Σ^i` are the self-dual half of the 6 `Cl(3,1)` bivectors (`3φ+3θ`) —
  the C3 carrier. (Euclidean signature for real self-dual forms; the Lorentzian
  chiral version is the analytic continuation, `*²=−1`.)

## 4. Linearized content (from `05`)

Around flat `B_0` (flat tetrad), fluctuations `δB` give: the simplicity + torsion
constraints + diffeo/Lorentz gauge leave **exactly 2 propagating modes**, helicity
±2 (`05` C1, machine-checked `graviton_dof`/`graviton_helicities`). The trace
sector reproduces `□Ω = f_g ρ_grav` (`05` C2). So the posited action's weak-field
limit is precisely the `05` result, now with `B` as the fundamental carrier.

## 5. Status

| piece | status |
|---|---|
| action written with SCP identifications (`B`∈`Cl(3,1)` bivectors, source `ρ_grav`) | **done** (§1) |
| EOM (simplicity / torsion-free / Einstein) | **stated** (§2; standard Plebański) |
| simplicity forces tetrad/metric | **verified** (§3; `08_*.py`, ~1e-16) |
| induced metric `g(B)` via Urbantke | **verified** (§3; recovers `g` to ~1e-16) |
| 2 TT DOF + trace law = v59 | **shown** (§4; `05`, Lean) |
| self-dual algebra (`η^i` orthogonality, dims) | **machine-checked** (`G9Plebanski.lean`) |
| action **derived from the OBE** | **open** (posited here, not derived) |
| nonlinear `ρ_grav` coupling normalization | **open** |

## 6. Net G9 picture (end of the `04`–`08` arc)

1. `04` — ±2 graviton requires soldering (internal index co-rotating); machine-checked.
2. `05` — soldered carrier gives **exactly 2 TT DOF** and contains the v59 scalar
   law as its trace sector.
3. `06` — Schur: no Lorentz fits inside `Cl(7)_even` ⇒ carrier **forced** to the
   spacetime factor.
4. `07` — built `Cl(3,1) ⊗ Cl(7)_even` (v59's own factorization); the carrier
   commutes with the entire internal sector + triality.
5. `08` — the carrier is a Plebański 2-form whose simplicity constraint induces a
   metric (`g(B)` via Urbantke, verified), giving GR sourced by `ρ_grav`.

**Remaining genuinely open:** deriving this action from the OBE multivector
dynamics (vs. positing it), and fixing the `ρ_grav` coupling magnitude / matching
the `α²¹` strength nonlinearly. The *structure* of the gravity sector — carrier,
its home factor, compatibility, induced metric, DOF, weak-field law — is now
established and largely machine-checked.
