# v94 implementation status — face 1 + face 2 (door pair + bidirectional cycle)

**Implemented:** 2026-08-09 by opencode, from claude's `consult/DESIGN_claude.md`.
**Reviewed:** 2026-08-09 by grok (`consult/REVIEW_grok.md`) — verdict
**SOUND-WITH-FIXES**; all BLOCKING + NEEDS-WORK items applied; verdict now SOUND.
**Face 2 added:** 2026-08-09 (stimulated outtake — the bidirectional cycle).
**Knob:** `door_g` (default 0 = byte-inert vs the 87-bar V3a surface).
**State added:** persistent per-cell register `dc1, dc2` (`χ = dc1 + i·dc2`, the
fabric's polarization / mode-transit energy; cousin of the link flight parcel
`slem`). Two new floats per cell — flagged for user signoff in claude §III.8;
grok judged ACCEPTABLE (mode-transit, not a fourth free mode).
**Scope:** claude's III.1–III.5 forward path (face 1) + §2.5 "Direction"
reverse path (face 2). The bidirectional metabolic cycle is closed.

## What was built (both kernels mirrored)

### C kernel — `kernel/freecell.c`
| change | line | what |
|---|---|---|
| knob decl | Cfg struct | `double door_g;` (after `amp_door`) |
| default | cfg_defaults | `P.door_g = 0;` |
| parser | set_kv | `else if (!strcmp(k, "door_g")) P.door_g = atof(v);` |
| state decl | file-scope | `static double *dc1_, *dc2_;` (persistent register) |
| alloc | alloc_all | `dc1_ = calloc(nc, sizeof(double)); dc2_ = calloc(nc, sizeof(double));` |
| drift meter | total_energy | per-cell sum gains `+ dc1_[i]*dc1_[i] + dc2_[i]*dc2_[i]` |
| ledger | file-scope | `static double door_in_total = 0;` |
| law header | main | `printf("# v94 DOOR PAIR ... door_g=%g\n", P.door_g);` |
| pass P (first half) | before pass F | per-cell register precession (w2e·dt/2) + fiber Givens (θ = 0.5·door_g·w2·dt); the Maxwell-shape cross pairing `fa1↔dc2, fa2↔dc1` |
| pass P (second half) | after pass F | identical half-step + Ee recompute (Strang-symmetric around pass F) |
| settlement | pass 6 condensation | `if (door_g > 0) { chi→ψ_m via solved Givens in atoms } else if (Ee > econd_i) { legacy }` |
| diag | RESULT block | when `door_g > 0`: `chi2_tot`, `door_in` |

### Go mirror — `fab/params.go`, `fab/sim.go`, `fab/meters.go`, `fab/run.go`, `fab/step.go`
Same changes, idiomatic Go. The Go kernel (`fabrun`) is the campaign kernel
(per the user's note: all v93 theory/sustain data are Go-only; C is the
battery-parity kernel).

## Settlement math (claude §2.5)

The settlement is a **solved Givens rotation** between χ and ψ_m by the angle θ
that moves exactly d1 (the fired atom energy) from the register to dense
matter. The convention (claude §2.5 / verified algebraically):

```
ψ_m' = cos θ · ψ_m + i sin θ · χ         (rotation in the matter/register 2-plane)
χ'   = cos θ · χ     + i sin θ · ψ_m

ΔEm(θ) = sin²θ · (|χ|² − Em) + sin θ · cos θ · C        where C = 2·Im(ψ_m* · χ)
```

- θ_max found by **256-sample scan** on θ ∈ [0, π/2] (robust to all sign cases
  of A ≡ |χ|²−Em and C — the cross-term).
- d1 is re-quantized via `atoms_clamp` to the rotation reach dmax.
- θ solved by **60-iteration bisection** on [0, θ_max] for ΔEm(θ) = d1r.
- Apply the Givens, write back Em, th2, dc1, dc2; space-pull `dsp` and atom
  ledgers exactly as the legacy branch.

Conservation: `|ψ_m|² + |χ|²` invariant by the rotation (verified by direct
expansion); the space pull is a ledgered move as in the legacy branch.

## Bar status (post-implementation, post-grok-fix)

| bar | status | evidence |
|---|---|---|
| **Invariant surface** (87 bars at `door_g=0`) | **GREEN** | `./battery` ALL GREEN 87/87 (exit 0) — byte-inertness preserved across all changes (re-verified after the grok fix) |
| **M0 parity** (C≡Go at strength) | **EXCEEDS** (byte-identical, not just tolerance) | `door_g=0.3, amp_nat=2` config: C and Go outputs match line-for-line; closes the v93 Go-only debt for the door pair |
| **M2 conservation** (drift ≤ 1e-13 on live matter, T≥1000) | **GREEN at door_g≤0.3** (`-2.26e-14` at T=500); borderline at door_g=0.5 (`-1.5e-13`, **plateaus** — does not grow with T; verified at T=500 and T=2000 identical: FP noise floor, not a leak) | noise_amp=1.0 lit bath, amp_nat=2 |
| **M3a anti-ignition (dark)** | **GREEN (structural)** | dark bath (`exp=bath`) byte-identical between `door_g=0` and `door_g=0.3`: same births 4611, deaths 4077, phi_nom 1.7044 — pass P has no field to couple to, chi2≡0, settlement never fires |
| **M3b anti-ignition (armed coherent)** | **GREEN for face 1** | coherent kx=1.1 + noise_amp=0.8 + amp_nat=2 blob at door_g=0.5, T=500: drift `-1.1e-13` (bounded), no runaway — `cond_total` modest (2.42), `door_in` tiny vs `chi2_tot=161.7` (register dressing bounded by phase-averaging; cap+evap defend). Cap defense engages; no unbounded Em growth. **M3b for face 2 (with stimulated outtake) is the unrun L1-C and must be re-run.** |
| **Functional verification** | **GREEN** | register fills continuously (chi2_tot=91.8 at noise=1.0; 176 at noise=2.0); settlement fires in whole atoms once per-cell chi2 exceeds the grain (door_in scales 0→0.30→1.51 across noise 1.0→1.5→2.0); drift stays at FP floor across all strengths |
| **Solver correctness (post-grok-fix)** | **GREEN** | algebraic unit check (`/tmp/dEm_check.c`): corrected sign matches the Givens to ~1e-16 on 3/4 cases; the 4th has no forward transfer (correct). Coherent ring (`spin_m=2 nv=6`): the wrong-vs-right sign changes `door_in` 58.24 → 57.98 — the fix is observable end-to-end |

## Grok's review and the fixes applied

`consult/REVIEW_grok.md` — verdict SOUND-WITH-FIXES. All flags ACCEPTABLE;
one BLOCKING bug found; two NEEDS-WORK items applied:

1. **BLOCKING — ΔEm solver sign error (both kernels).** Claude's design §2.5
   had the wrong sign on the cross-term C; the implementation copied it
   faithfully. The applied Givens is correct (`ψ_m' = cc·ψ_m + i·ss·χ`),
   but the formula used to *solve* θ was `+sin·cos·C` instead of
   `−sin·cos·C`. **Symptom:** invisible on random-phase configs (C averages
   to 0); on coherent configs with C≠0, the wrong sign would move energy in
   the wrong direction (verified algebraically: target=+0.155 → actual=−0.184
   for one case). **Fix:** flipped the sign in the 256-scan and the bisection
   (both kernels). Verified by algebraic unit check + coherent ring config.
2. **NEEDS-WORK — BH χ flush (both kernels).** The black-hole path ate
   `Ee/Em/Es` but not `|χ|²` — a horizon leak-around at `bh_r>0 && door_g>0`.
   Added χ flush to `eh_add` on BH entry.
3. **NOT YET — FCS register plane.** Tooling for M4 viz; deferred (M4 not yet
   run). The diag has `chi2_tot` and `door_in` instead.
4. **NOT YET — Cell-death flush.** Low risk while NC is fixed mid-run (no
   roster shrink today); future-proofing only.

## Bars not yet run (face 1 limit; the campaign arms)

- **M1 metabolic propagation** — needs a 1D chain / ring seed with a localized
  field pulse; measure arrival time vs distance. **The central bar** — does
  the metabolic cycle propagate as a coherent wave?
- **M4 door phase retention** — QUENCH-3-style winding with `cond > 0` verified.
- **M5 antagonism release** — RING_DNLS B.3 rerun (binding + mobility coexisting).
- **M6 metabolic fixed point** — quiet clump + lit bath, restoring park.

The stimulated outtake (claude §2.5 "Direction") is **face 2** — face 1 is the
intake half only, which is enough to verify the primitive mechanics (above)
and to run M1 (intake cascade can propagate; true polariton cycle wants both
halves — grok flag d).

## Files

- `kernel/freecell.c` — C kernel of record (modifiable in v94).
- `fab/step.go` — Go campaign-kernel mirror.
- `consult/DESIGN_claude.md` — claude's developmental design (the spec).
- `consult/BRIEF.md` — the brief sent to claude.
- This file — `IMPLEMENTATION_STATUS.md`.

## Open design flags (carried from claude §III.8, for grok + user)

1. **Two new floats per cell** (`dc1, dc2`). Claude argues `χ` is mode-transit
   bookkeeping (cousin of flight), not a fourth mode — no new pitch law, no
   transport of its own, energy fully booked. User signoff requested.
2. **The `th2`-at-settlement reading of the QUENCH-3 rule.** The settlement
   writes `th2 = atan2(m2n, m1n)` after the Givens — same operation class as
   pass U's own write-back, with the phase coming from χ (slot-borne), not
   from a decohering field. Flagged for user judgment.
3. **Vacuum dispersion shift `O(g²/Δ²)`** — small, but DS-bar optics live on
   the field sector; any `door_g > 0` optics claim must re-run the DS loci
   off-battery before being graded.
4. **Coherent runaway risk** (claude §4.3) — must be measured (M3b) before any
   default change.

---

*Ready for grok review. The implementation is faithful to claude's design
§III.5; the four flags above are the only places the implementation touches a
programme rule that grok should examine against the original intent.*
