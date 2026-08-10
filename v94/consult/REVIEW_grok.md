# v94 door-pair implementation review (grok)

**Reviewer:** grok  
**Date:** 2026-08-09  
**Scope:** face-1 implementation (`door_g`, pass P, χ-register settlement) vs
claude's `DESIGN_claude.md` §III.1–III.5 and programme constraints.  
**Kernels read:** `kernel/freecell.c`, `fab/step.go` (+ `sim.go`, `params.go`,
`meters.go`, `run.go`).  
**Evidence surface:** `IMPLEMENTATION_STATUS.md`, `battery.log` (ALL GREEN 87).

---

## Executive verdict

| Item | Judgment |
|---|---|
| Implementation vs design (structure) | Review written to `v94/consult/REVIEW_grok.md`.

## Summary

**Whole verdict: SOUND-WITH-FIXES**

The face-1 skeleton matches Claude’s design: Strang pass P around F, bare-pitch `g = door_g·w2`, Maxwell-shape fiber Givens, click grain on `|χ|²`, C≡Go, 87/87 at `door_g=0`, drift at the FP floor.

### Blocking fix (before M4/M5/M6)

**ΔEm solver has the wrong cross-term sign** (design §2.5 + both kernels).  
The applied Givens is correct and conserves `Em+|χ|²`. The formula used to *solve* θ is:

```text
ΔEm_code = sin²A + sin·cos·C     // wrong
ΔEm_true = sin²A − sin·cos·C     // algebra of the stated Givens
```

When `Em=0` (first clicks onto empty matter), `C=0` and both agree — which is why noise/functional tests look green. On **existing matter with phase structure**, `door_in` / `s_pull` / atom demand can disagree with the energy the rotation actually moves. Fix: flip the sign in the 256-scan and bisection (`freecell.c` ~2361/2372, `step.go` ~1414/1428) and add an algebraic unit test.

### Flags

| Flag | Verdict |
|---|---|
| (a) χ as two floats | **ACCEPTABLE** (user signoff) — mode-transit, not a fourth free mode |
| (b) `th2` write at settlement | **ACCEPTABLE** — arg of unitary Givens from χ, not QUENCH-3 heuristic |
| (c) vacuum O(g²/Δ²) | **ACCEPTABLE** — `door_g=0` byte-inert protects DS bars |
| (d) face 1 only | **ACCEPTABLE** for mechanics; M5/M6 need face-2 outtake |

### Also fix / add soon

- BH path does not eat χ (horizon leak-around if `bh_r>0` at strength)
- Design §5.1 death-flush of χ not coded (low risk while NC is fixed mid-run)
- FCS register plane missing (tooling for M4)

After the ΔEm sign fix, face 1 is ready for **E1 = M0+M1+M2** as Claude ordered; run **M3b early** before any lit-bath adoption talk.
ass U → … → pass P(τ/2) → pass F → pass P(τ/2) → pass 6 …
```

Second half recomputes `Ee = fa1²+fa2²` after the fiber rotation
(`freecell.c:2269–2270`). First half correctly leaves `Ee` stale: pass F
reads `fa1/fa2`, not `Ee`.

### 1.4 Settlement Givens component stencil — MATCH

Design / status complex form:

```
ψ_m' = cosθ · ψ_m + i sinθ · χ
χ'   = cosθ · χ   + i sinθ · ψ_m
```

Component form in both kernels (`freecell.c:2383–2386`, `step.go:1442–1445`):

```
m1'  = cc·m1  − ss·dc2
m2'  = cc·m2  + ss·dc1
dc1' = cc·dc1 − ss·m2
dc2' = cc·dc2 + ss·m1
```

These are algebraically identical. Norm `|ψ_m|² + |χ|²` is exact under the
applied rotation (not under the *solver* — see §1.5 / §5).

### 1.5 ΔEm(θ) formula used to *solve* θ — MISMATCH (design + impl)

Design §2.5 and both kernels use:

```
ΔEm(θ) = sin²θ · (|χ|² − Em) + sinθ·cosθ · C ,   C = 2·Im(ψ_m* · χ)
```

Direct expansion of the **stated** Givens gives the **minus** cross-term:

```
|ψ_m'|² = cos² Em + sin² |χ|² − 2 cos sin · Im(ψ_m* · χ)

ΔEm_true(θ) = sin²θ · (|χ|² − Em) − sinθ·cosθ · C
```

Numerical check (2000 random draws): wrong formula max error ~O(1);
corrected formula max error ~1e-15. **Special case Em = 0 (⇒ C = 0):** both
agree exactly — which is why face-1 functional tests (noise filling an empty
register, first clicks onto vacuum matter) look green and conservation still
telescopes.

**What breaks:** when settlement fires onto *existing* dense matter with a
nonzero cross-term (exactly the M4 phase-retention / M5 metabolic geometry),
the bisection solves the wrong equation. Consequences:

- `d1r` (atom-clamped demand) ≠ actual `|ψ_m'|² − Em` after the Givens;
- `door_in_total`, `cond_total`, and `dsp = s_pull·d1r` are keyed off the
  wrong amount;
- the Givens itself still conserves `Em+|χ|²` (drift meter stays honest);
- the **click grain intent** ("move exactly d1 of energy in whole atoms") is
  violated whenever C is material.

**Root cause:** design §2.5 wrote the wrong sign; the implementer copied it
faithfully and verified the *rotation* (correct) rather than the *solver
identity* (wrong). Fix: flip the cross-term in the scan and bisection only
(both kernels, one line each). Re-verify with a unit test:

```
for random (ψ_m, χ, θ): |ΔEm_formula − (|ψ_m'|² − Em)| < 1e-12
```

and a live case with Em>0, C≠0, assert `|Em_new − Em_old − d1r|` at the
bisection tolerance.

### 1.6 Click grain / beat gate / workfn — MATCH (machinery)

- Settlement only inside existing `beat_fire` (`cbeta` unchanged).
- Workfn text applied to `|χ|²` instead of `Ee` (`chi2 > econd_i`).
- `atoms_fire` / `atoms_clamp` with `eF = A₀ w1e/2π`, `eD = A₀ w2e/2π`,
  `qcnvD` credits — verbatim.
- Legacy `Ee`-condensation is `else if` (no double intake) at `door_g > 0`.
- Field `fa1/fa2/Ee` not touched at settlement (transfer already in pass P).
- Space pull `dsp = s_pull · d1r` from `Es`, same as legacy.

Atom machinery is intact; the §1.5 bug means the *amount* fed into that
machinery is not always the amount the rotation moves.

### 1.7 Register precession sign — MATCH kernel dense convention; design text off

Design §5.2 writes `χ ← e^{−i·w2e·dt/2}·χ`. Implementation does

```
dc1' = cc·dc1 − ss·dc2
dc2' = ss·dc1 + cc·dc2     // e^{+i ang}, ang = w2e·dt/2
```

which is the **same** convention as pass-U dense precession
(`freecell.c:1647–1648`: `dm1/dm2` and `th2 += w2e dt` for `ψ_m = √Em e^{i th2}`).
Field pass F uses the opposite Schrödinger sign (`e^{−i}`). The substrate
already has this field/dense sign split; matching χ to dense is the correct
choice for settlement composition with `ψ_m`. Continuum equations in design
§2.2/§3.1 idealized both sectors with one sign — not a kernel bug.

### 1.8 Deferred / not implemented (design §5, face-1 scoped)

| Design item | Status | Severity for face 1 |
|---|---|---|
| Stimulated outtake (ΔEm < 0 branch) | Deferred face 2 | In scope of deferral |
| Cell-death flush `field_inject(i, \|χ\|²)` (§5.1) | **Missing** | Low if NC fixed mid-run; see §5 |
| FCS dump of register plane (§5.1) | **Missing** | Tooling; needed before M4 viz |
| `amp_nat=0` magnitude-only settlement (§5.4) | Not done — full Givens always | **Better** than design; OK |
| Pass named "D" | Named "P" in code | Cosmetic |

---

## 2. Conservation proof — does the sum telescope?

### 2.1 Pass P — exact

Per half-step, per cell:

1. **Register precession** is U(1): `|χ|²` invariant; field untouched.
2. **Fiber Givens** is SO(2) on `(ψ_e, χ)`: `|ψ_e|² + |χ|²` invariant
   (cross terms cancel; verified algebraically and numerically).

Field precession and spatial hops in pass F preserve `Σ Ee` pairwise as
before. No path in pass P creates or destroys energy.

### 2.2 Settlement Givens — exact on the spinor pair

Applied rotation: `|ψ_m|² + |χ|²` invariant. Then:

```
dsp = min(s_pull · d1r, Es − es_floor⁺)
Em  ← |ψ_m'|² + dsp
Es  ← Es − dsp
```

So `Em + Es + |χ|²` after commit equals `Em + Es + |χ|²` before, up to the
**identity that d1r is only a label** for the space-pull fraction — the
spinor move is whatever the Givens did. Total inventory still telescopes.

### 2.3 Drift meter coverage

`total_energy()` / `totalEnergy()` include per-cell
`Es+Em+Ee+|χ|²` plus flight (`slem`) plus `Eh_total` sink. Every pass-P and
settlement write to `dc1/dc2` is inside that sum. **No silent χ store outside
the meter.**

### 2.4 Paths that conserve total energy but violate intended bookkeeping

| Path | Total energy | Physics/ledger intent |
|---|---|---|
| Wrong ΔEm solver (§1.5) | Conserved (Givens exact) | Click amount / s_pull / door_in wrong when C≠0 |
| Cap test `Em+Ee > cap` ignores `\|χ\|²` | Conserved | Large virtual dressing never trips evap; χ is transit, so arguably OK — document |
| Black hole (`bh_r>0`): eats Ee, Em, Es; **not χ** | Conserved (χ stays in sum) | Horizon **under-eats**; χ is a BH-evasion reservoir; pass P can still feed χ from neighbors' field | 
| Cell leave-roster without χ flush | Would drop energy if NC reindexed | NC fixed mid-run today; design §5.1 still requires the rule for correctness under any future roster shrink |

**M2 GREEN is real** for total-energy conservation. It does **not** certify
that `door_in` equals the energy that left the register.

### 2.5 Theorem statement (what actually holds)

Over one step, with the current code:

```
Σ_i (Em + Ee + Es + |χ|²) + Σ_slots flight + Eh
```

changes only by FP roundoff under: pass P, pass F, pass U, settlement
Givens, space-pull ledgers, legacy evap/radiance (unchanged). The atom
ledgers (`cond_total`, `door_in_total`) are **not** conserved quantities and,
with the §1.5 bug, are not equal to the spinor ΔEm either.

---

## 3. Anti-ignition — face 1 and the face-2 risk

### 3.1 Claude §4.3 four layers under face 1 (intake only)

| Layer | Still sound? | Notes |
|---|---|---|
| 1. Detuning | **Yes** | `g = door_g·w2` bare; vacuum `\|Δ\|=w2−w1` maximal; virtual admixture reversible |
| 2. Phase averaging | **Yes** | Random-phase bath → incoherent drive of χ; no secular pump in expectation |
| 3. Settlement gates | **Yes** | workfn / atoms / f_conv / cap+evap still gate *real* conversion |
| 4. Zero-throughput coupling | **Yes** | g-rotation is norm-preserving; cannot amplify total `Ee+\|χ\|` |

M3a (dark bath byte-identical `door_g=0` vs `0.3`) is the right structural
check for face 1: with no field, pass P is a no-op on empty χ, settlement
never fires. **Do not over-read it** — it does not exercise layers 1–2 under
a lit incoherent bath (that is a stronger M3a arm, worth a cheap add-on).

### 3.2 Face 2 (stimulated outtake) — new risks and gates

Face 2 adds the signed-C outtake branch: matter anti-phased with χ can dump
`Em → χ → field` coherently. That is the missing half of the metabolic cycle
and also the **amplifier** of coherent runaway (§4.3 named risk):

- **Positive feedback under phase-coherent lit bath:** more load → smaller
  `|Δ|` → larger F → more intake; outtake that is phase-matched to the bath
  can re-seed the field coherently rather than radiate away.
- **Gate recommendations before any default `door_g>0`:**
  1. Fix §1.5 first (wrong C sign makes outtake geometry nonsense).
  2. Run **M3b early** (design order): coherent lit bath over seeded clump;
     bar = Em bounded by cap/evap; **record** pumping threshold.
  3. Outtake must use **outtake atom credits** (`eF2/eD2` swap, `qcnvF`)
     exactly as design §5.3 — do not share intake credits.
  4. Prefer an explicit `door_out` knob or `door_g` with outtake behind a
     second face flag so face-1 campaigns stay reproducible.
  5. Cap/evap should continue to see Em (not χ); optional diagnostic:
     alarm if `Σ|χ|²` grows secularly without settlement.

Face 1 alone does **not** create a new ignition channel beyond virtual
optics; face 2 must be treated as the L1-C configuration v93 left unrun.

---

## 4. The four open flags

### (a) Two new persistent floats per cell (`dc1, dc2`) — **ACCEPTABLE**
**(user signoff still required by programme rule)**

χ is mode-transit bookkeeping, cousin of `slem`:

- no new pitch law, no spatial hops of χ, no new conversion endpoints;
- energy fully booked in the conserved sum;
- everything entering χ comes from Ee by unitary rotation; everything
  leaving goes to Em by atom click or back to Ee by the same pair.

It is **not** a fourth free energy mode in the programme sense (those have
independent transport and conversion ends). It *is* two new persistent
floats — that alone needs the user's explicit OK, as claude requested. My
recommendation: **accept** as door-flight.

### (b) Settlement writes `th2 = atan2(m2n, m1n)` — **ACCEPTABLE** (in spirit)

QUENCH-3 failed because **heuristic phase mixes** wrote field arg into
`th2` (`amp_door` / `qp_phase` trig blends) while the door itself was
magnitude-only — decohering, non-unitary, R2d≈0.02.

Settlement here:

- phase lives in χ between clicks (slot-borne carrier);
- `th2` is the **arg readout of a norm-preserving Givens** composing χ into
  `ψ_m` — same operation class as pass U's write-back after hops;
- field is **not** the phase source at the click.

Letter of "do not write th2 at the door" is broken; spirit of QUENCH-3 is
honored. **M4 is the measurement that decides**, not the review. Flag for
user judgment as design does; I would not block on the write itself.

### (c) Vacuum dispersion O(g²/Δ²) — **ACCEPTABLE**

- `door_g = 0`: whole pass P + settlement branch skipped → **byte-inert**
  87-bar surface (measured ALL GREEN). Protects DS optics bars completely
  for the invariant surface.
- `door_g > 0`: expected polariton dressing of the field branch; any optics
  claim at strength **must** re-run DS loci off-battery before grading.
  Not a bug — a campaign protocol.

No non-byte-inert path at default zero.

### (d) Face 1 intake-only vs full metabolic claim — **ACCEPTABLE** for
**face-1 scope; incomplete for the programme's central claim**

| Claim | Face 1 enough? |
|---|---|
| Primitive exists (virtual Ee↔χ, click χ→Em) | **Yes** |
| Conservation theorem | **Yes** (with §1.5 caveat on ledgers) |
| Anti-ignition dark | **Yes** |
| Metabolic *propagation* (M1) | **Partially** — intake cascade can propagate; true polariton cycle wants both halves |
| Antagonism release (M5) translation-by-conversion | **Needs face 2** for trailing-edge outtake |
| Metabolic fixed point I=O (M6) | **Needs face 2** |

Face 1 is a valid **mechanics** test and a correct gate before opening the
campaign arms. It is **not** a complete test of "the metabolic cycle
propagates." Orchestrator language in `IMPLEMENTATION_STATUS` should not
equate face-1 green with M1/M5/M6 readiness beyond "primitive installed."

---

## 5. What's missing or wrong (specific)

### BLOCKING for M4/M5/M6 grading

1. **ΔEm cross-term sign error** (design §2.5 + both kernels).  
   - C: `freecell.c:2361`, `2372` — change `s*s*Aa + s*c*Cc` → `s*s*Aa - s*c*Cc`.  
   - Go: `step.go:1414`, `1428` — same.  
   - Optionally redefine `C = -2 Im(ψ_m*χ)` and keep the `+` form; pick one
     and document.  
   - Add a parity unit check (C and Go) on the identity.  
   - Re-run M2 after the fix (should stay green); re-run functional noise
     test (Em=0 path unchanged).

### NEEDS-WORK (before or with early campaign)

2. **Black-hole path ignores χ** (`freecell.c:2292–2308`, Go analogue).  
   When `bh_r>0`, flush `|χ|²` into `eh_add` (or `field_inject` then eat)
   and zero `dc1/dc2`. Otherwise χ is a horizon leak-around. Battery is
   safe at `bh_r=0`; any horizon bar at `door_g>0` is wrong without this.

3. **Design §5.1 cell-death / repack rule** not coded. Today NC is fixed
   for the run and carve is init-time, so no measured leak — but the design
   theorem ("flush before leave roster") is incomplete in code. Implement
   on any path that removes or reindexes a cell; zero-init on birth already
   holds via `calloc`/`make`.

4. **FCS register plane** not added (`fab/fcs.go` has no dc1/dc2). Needed
   for M4/M5 diagnosis; not physics-blocking.

5. **M2 duration:** status used T=500; design asks T≥1000 on live matter.
   Cheap to re-run; not a code defect.

6. **Cap inventory** is `Em+Ee` only. Document that `|χ|²` is intentionally
   excluded (transit), or include it if campaign objects run huge virtual
   dressing under cap. Preference: document, don't change without a bar.

### Not bugs (do not "fix")

- Bare `w2` in θ (not `w2e`) — load-bearing.
- Full Givens settlement at `amp_nat=0` — strictly better than magnitude
  degradation; keeps phase in χ→ψ_m even without pass U.
- Pass named P not D — cosmetic.
- Precession `e^{+i}` on χ — matches dense sector.
- Outtake deferred — explicit face-2 scope.

### Constraint checklist (programme §IV)

| Constraint | Status |
|---|---|
| No background | OK — χ is dynamical per-cell transit, not a fixed fiber ledger |
| Energy never destroyed | OK as theorem of rotations; meter includes χ |
| No imported field / SI | OK — `g = door_g·w2` |
| Door stays special (atoms) | OK machinery; §1.5 undermines exact atom move when C≠0 |
| No quark road | OK — one register, no flavor |
| No sum-then-square | OK — pairwise / fiber Givens only |
| No feedback / shadow / bias | OK — unitary pair, not a want bias |
| No heuristic th2 from field | OK — th2 from Givens arg; see flag (b) |
| Ratchet / no bar softening | OK — 87 green at door_g=0; no bar edited |

---

## 6. Verdict table

### Flags

| Flag | Verdict | One-line reason |
|---|---|---|
| **(a)** Two floats / χ | **ACCEPTABLE** (user signoff) | Mode-transit cousin of `slem`, not a fourth free mode |
| **(b)** th2 at settlement | **ACCEPTABLE** | Arg-readout of unitary Givens from χ; not QUENCH-3 heuristic |
| **(c)** Vacuum O(g²/Δ²) | **ACCEPTABLE** | `door_g=0` byte-inert protects DS; re-run optics at strength |
| **(d)** Face 1 only | **ACCEPTABLE** for mechanics; incomplete for full cycle | M0–M3a yes; M5/M6 need outtake |

### Whole implementation

**SOUND-WITH-FIXES**

The primitive is installed correctly: Strang door pair, bare-pitch coupling,
Maxwell-shape fiber Givens, click grain, C≡Go, battery inertness, measured
conservation floor. Before opening M1/M4/M5/M6 as gradeable campaign arms:

1. **Must fix:** ΔEm cross-term sign (both kernels) + algebraic unit test.  
2. **Should fix with it:** BH χ flush; document cap vs χ.  
3. **Should add soon:** FCS χ plane; full §5.1 death flush for future-proofing.  
4. **Then:** M3b before any lit-bath adoption talk; face 2 outtake only after
   (1) and M3b.

After (1), face 1 is in good shape to start **E1 = M0+M1+M2** as design
§III.8 orders. Do not treat current `door_in` numbers on live matter with
Em>0 as atom-accurate until the sign fix lands.

---

## Appendix A — Site map (review anchors)

| Site | C | Go |
|---|---|---|
| knob | `freecell.c` Cfg `door_g`, default 0, `set_kv` | `params.go` `DoorG` |
| state | `dc1_`, `dc2_`, `alloc_all` | `sim.go` `dc1/dc2`, `allocAll` |
| drift | `total_energy` +`\|χ\|²` | `meters.go` `totalEnergy` |
| pass P ½ | ~2170–2199 | `step.go` ~1190–1215 |
| pass F | ~2202–2247 | ~1217–1270 |
| pass P ½ | ~2249–2271 | ~1272–1294 |
| settlement | ~2324–2411 | ~1373–1478 |
| RESULT | ~4684–4690 | `run.go` ~1832–1839 |

## Appendix B — Algebra snapshot (settlement)

With `ψ_m = m1 + i m2`, `χ = dc1 + i dc2`, `C = 2(m1 dc2 − m2 dc1)`:

```
Applied Givens (implemented, correct):
  ψ_m' = cosθ ψ_m + i sinθ χ
  χ'   = cosθ χ   + i sinθ ψ_m

True energy move:
  ΔEm(θ) = sin²θ (|χ|² − Em) − sinθ cosθ C

Implemented solver (wrong sign on C):
  ΔEm_code(θ) = sin²θ (|χ|² − Em) + sinθ cosθ C

Em = 0 ⇒ C = 0 ⇒ both agree (explains green face-1 noise tests).
```

---

*End of REVIEW_grok.md. Orchestrator: fix §1.5 before grading phase-sensitive
bars; flags (a)(b) need user eyes, not code blocks.*
