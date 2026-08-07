# L-1 PATH A — grok CLI review (2026-08-07)

Requested by the user before implementation. Grok's review of
`L1_PATH_A_DESIGN.md`, transcribed and condensed. **Verdict: Path A as
written is "feedback in disguise" (form 3), not the native complex-transport
completion it claims to be.**

## The core critique (e): Path A is feedback form 3

Under §3 as written:
- the magnitude path is still `gate·head·mob` **first**,
- the phase cargo is still sourced from `th2`,
- `R` is still the Phase-M composition,
- the drive is still a **multiplicative correction** on that magnitude.

So it makes the **gate** a function of composed arrivals — a better feedback,
but NOT S2's "transfers carry complex amplitude." The diagnosis in
L1_FINDINGS is sharper than Path A's fix.

**It stops being feedback-in-disguise only if ≥1 of:**
1. Parcel phase is free/composed (departure phase = `arg(R)` or integrated
   path phase, not only `m·th2`).
2. Want magnitude is constitutive of composition (`f` from `|composed
   amplitude|` or complex mobility) — not a real want with a phase garnish.
3. In-flight state is the ONE complex object that both energy and phase ride
   (not `slem` + a side shadow).

Until then, Path A is L-1 form 3, and "deeper kernel change" is not yet met.

## Chicken-and-egg not solved (L2)

Sustain ≠ ignite. With gg≈0.006 and no deposits, R≈0, boost off — same egg.
Path A can **maintain/amplify existing coherence** but not **resurrect a
dead fifth** without a separate nucleation path. L2 acceptance should say
**survival/maintenance of a live fifth**, not revival.

## Cross-cutting flaws flagged

1. **Per-cell R_j is a resolution regression.** Phase M's shadow is
   `[slot][dir]`, chart-folded, path-rotated. Collapsing to one R_j destroys
   the directionality grain (e3b needs it) and p:q chart separation (L2 needs
   it). **Use the existing per-slot `sre_/sim_` as the emission basis.**
2. **Acceptance jumps over intermediate evidence.** Need staged bars
   ("median |cos| improves and seed variance falls at small amp_coh"), not a
   cliff — two forms already "engaged but didn't deliver."
3. **L1-A (unison) and L2 (chord) fight on the anti-ignition gate.** One knob
   `amp_coh` for both is a smell — separate the faces/gates.
4. **"Smaller than it sounds" is true only if it's form 3.** If it's the S2
   completion (complex want, phase cargo, both kernels, C≡Go at strength),
   it's big. Pick one ambition; don't claim S2-native while sizing a patch.

## Anti-ignition (c): the weakest claim

The "bath starved, |R|~1/√N" math is **wrong** — Phase M measured the bath's
ρ_coh ≈ **0.77** (partially coherent, not incoherent). The ρ_coh prefactor in
the gate-boost is **dangerous** (rewards bath unison). Drop ρ_coh; put anti-
ignition on the measured alignment isotropy / birth–glow bars; pre-arm
alignment-only and/or chart-split.

## Conservation (a)

Ledger conservation OK, but "redistribution / no runaway by construction" is
**overclaimed** unless the zero-sum renorm returns. Name **coherence-runaway**
as the real risk.

## Grok's required changes before sign-off

1. Rewrite C-1: ledger yes; drop "redistribution by construction"; name
   coherence-runaway.
2. Rewrite C-3: delete the 1/√N math; cite ρ_coh~0.77; anti-ignition on
   measured alignment/birth-glow; pre-arm alignment-only + chart-split.
3. **Pick the emission form: prefer a signed zero-sum current on the PER-SLOT
   shadow (L0 deep + FLOW architecture); demote ρ_coh·gate-boost.**
4. Specify phase cargo: if outgoing phase stays `m·th2`, call it form 3; if
   `arg(R)`, write the deposit/emit rules explicitly.
5. Drop per-cell R_j; use existing slot-direction state.
6. Separate L1-A (unison current) from L2 (chord maintenance) gates/knobs.
7. Honest naming: "close the Phase-M loop with a directional emission
   coupling" — not "native complex-transport completion" — unless complex
   want + phase cargo are in scope.
8. L2 language: maintenance of live fifths; seed protocol explicit.

## S2 scorecard (grok)

| S2 requirement | Path A |
|---|---|
| dense transfers carry complex amplitude | **No** — transfers still real; side account complex |
| phase composes along paths | Yes (already Phase M) |
| translation IS the current | aspired via alignment boost; not derived |
| fifth survives coherence | unproven; chicken-egg; per-cell R risks chart scramble |
| exchange / flux-machine | deferred (ok) |
| invariant surface | plausible if inert default |

## Bottom line (grok)

Path A correctly learns that shadow-as-bias on a phase-blind magnitude want
cannot manufacture coherence. It does **not** yet implement the thing that
would: **make the dense parcel itself the complex object whose composition
*is* transport.** As a design sub-round it is useful if treated as "emission
coupling experiment on Phase-M state." As written it **oversells** that
experiment as the amplitude completion and **underweights** the bath ρ_coh
inversion that Phase L was supposed to be designed around.
