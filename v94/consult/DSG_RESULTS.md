# v94 session — DS-G results + experiment status

**Session ID:** `aebbc7eb-4716-4e30-8c22-484e115428f2`

## DS-G (the no-cancellation prediction test) — nuanced result

Claude's prediction (NOCANCELLATION_claude.md §5.2): **bright fringes are
conversion lines** — matter should accumulate where the field is bright,
stay at baseline where dark. DS-G1 tests this on the double-slit screen.

### What was run

The slit apparatus (`exp=slit L=64`) with:
- Source at x=12, wall at x=22 (two slits), screen at x=36.
- `door_g ∈ {0, 2.0}`, `door_avg=1`, `amp=2.0`, `T=60`, `e_cond=0`.
- FCS snapshots at t=60; screen Em(y) profile analyzed.

### Results

| config | door_in | door_out | screen Em | Em-Ee correlation |
|---|---|---|---|---|
| door_g=0 (legacy) | cond=58 | 0 | **accumulates** (peaks y=8, 31-33) | **r = −0.060** (≈0; not patterned) |
| door_g=2.0 (pair) | 110 | 141 | **zero** (all cycled back) | n/a |
| door_g=2.0, chi_load=1 | 110 | 141 | zero (chi_load inert at q_detune=0) | n/a |

### Findings

1. **The legacy door creates matter at the screen but NOT correlated with the
   field pattern** (Pearson r ≈ 0). The magnitude-only transfer loses the
   spatial phase information; matter appears stochastically, not at the
   bright fringes. This is the first measurement that distinguishes the
   legacy door from the coherent door pair.

2. **The door pair converts at the screen** (door_in=110) **but the coherent
   outtake cycles it all back** (door_out=141). At t=60 the screen Em is
   zero. The face-2 outtake fires because the slit's field is COHERENT —
   door_avg's synchronous detection doesn't suppress it (it only suppresses
   incoherent/random-phase outtake).

3. **chi_load is inert** in the optics regime (`q_detune=0`): the load never
   engages the pitch. For the register self-detuning to work, `q_detune > 0`
   is required — but that changes the optics (slows propagation in loaded
   regions, blurring fringes). This is a genuine experimental design tension.

### What this means for the no-cancellation prediction

- The prediction **bright fringes are conversion lines** is **not yet tested**.
  The mechanism that would confirm it (the door pair's coherent conversion
  preserving the spatial pattern) is implemented but **prevented from
  accumulating** by the face-2 outtake.
- The prediction would need **face 1 only** (intake, no stimulated outtake)
  to be directly testable. Currently face-2 always fires when door_g > 0.
- **Registered next step:** add a knob to disable the stimulated outtake
  independently (`door_out_enable=0`), re-run DS-G, and measure the Em-Ee
  correlation. If r >> 0, the prediction is confirmed: the coherent door
  pair writes a grating that the legacy door cannot.

## door_avg (claude's M4 fix) — implemented + tested

| bar | result |
|---|---|
| Byte-inertness at door_avg=0 | 87/87 GREEN |
| M4′ (lit-bath ring, door_avg=1) | R2 improved 0.24 → 0.32; door_out halved 8.56 → 4.49; below do-no-harm bar (baseline 0.62) but directional improvement |
| M3b at door_avg=1 (na=2.0) | drift FP floor; door pair self-regulates (net outflow at high energy) |

## chi_load (face-3 physics knob) — implemented + tested

| bar | result |
|---|---|
| Byte-inertness at chi_load=0 | 87/87 GREEN |
| Register self-detuning (chi_load=1, q_detune=1.2) | chi2 drops 17%, door_in/out drop ~30% — the register self-detunes when it loads the pitch |
| Inert at q_detune=0 (optics regime) | Confirmed: chi_load has no effect when q_detune=0 (the load never engages) |

## Paper format note (for claude's next revision)

The user noted: **two-column format should start AFTER the abstract and TOC,
not before.** In LaTeX: use `\onecolumn` for the title block + abstract + TOC,
then `\twocolumn` for the body. In pandoc, this requires either a custom
template or raw LaTeX commands embedded in the markdown before §1.

## Next session priorities

1. **`door_out_enable` knob** (disable stimulated outtake independently) →
   re-run DS-G face-1 only → measure Em-Ee correlation → test the
   no-cancellation prediction.
2. **M1** (metabolic propagation, v_M curve) — needs a clean chain experiment.
3. **M5** (antagonism release) + **M6** (metabolic park) at door_avg=1, chi_load=1.
4. **Resume claude** with DS-G results + the `door_out_enable` design + the
   paper format fix + grok's review items.
5. **Implement claude's `chi_load` full face-3 bars** (U0–U5) at q_detune > 0.
