# Falsifiable predictions for freecell / fabrun

**If** the quantized harmonic graph law of `UNIFIED_LAW.md` is right, these
seven statements must hold on the existing substrate. Each names the knob, the
bar and the kill.

**None of these requires a kernel change** except P2, which is flagged as such
and is the only one that adds physics. P1, P4 and P7 are pure offline analysis
of frames that already exist.

All knobs below are confirmed present in `v93/fab/params.go`:
`k_rad`, `p_rad`, `rad_clock`, `wf_on`, `q_detune`, `amp_nat`, `amp_door`,
`seed_mw`, `noise_amp`, `kx`, `hop_order`.

---

## The two that separate this law from "radiance only + noise fuel"

The brief asks for at least two sharp, implementable predictions that
distinguish this law from the incumbent. These are P2 and P5. They are the
ones to run first if only two can be run.

---

## P1 — the spectral dimension of the foam sets `p_rad`

**Claim.** `p_rad = d_s + 1` (Theorem 4.1). The adopted `p_rad = 4` is the
assertion that the free-cell contact graph has `d_s = 3`, and nothing else.

**Apparatus.** Offline. Take any settled `.fcs` frame; build the contact graph
from cell positions with the kernel's own contact rule (`r0`, `dmin`); run
`qgraph.ReturnDimension(g, 128, sources)`. That estimator reproduces exact tori
to `±0.02` and needs no diagonalisation, so foam size is not a limit.

**Bar.** `d_s = 3.00 ± 0.15` over ≥3 independent settled frames and ≥2 box
sizes (to separate finite-size effects from geometry).

**Kill.** `|d_s − 3| > 0.3`. A jammed or near-jammed contact network is exactly
the case where the vibrational spectrum is anomalous, so this is a live
possibility, not a formality. If `d_s ≈ 2`, then `p_rad` should be 3, and the
adopted table over-taxes loaded matter by a full power.

**Why it matters.** This is the only bar that can invalidate an *adopted*
constant. It should be run before anything else in this list.

---

## P2 — the dual intake law makes matter a thermometer *(kernel change)*

**Claim.** Theorem 3.4: with the intake channel written as the dual of radiance
(same `k_rad`, same `p_rad`), a hoard in a bath at load `xb` equilibrates at
`x* = xb` **exactly**, for every `xb`.

**Apparatus.** Add the spontaneous dual intake
`I = k_rad * cap * xb^p_rad` at the door, behind a knob defaulting to 0
(byte-inert vs the current surface, per the ratchet). Then sweep the ambient
field load `xb` at fixed everything else.

**Bar.** `x*` tracks `xb` with slope `1.00 ± 0.05` over `xb ∈ {0.3, 0.6, 1.0,
1.5}`, each held ≥500 t.u.

**Kill / discriminator.** *Under the incumbent law `x*` is bath-independent* —
the radiance fixed point sits where the tax balances whatever supply happens to
exist, and `RADIANCE.md` reports it as a fixed `x* = 0.62 ± 0.02`. Under this
law `x*` is a linear function of the bath. **Slope 0 kills the duality; slope 1
kills "radiance only".** There is no way for both to be right.

**Cost note.** This is the one item that touches step physics. It is a small
addition (one term, one knob, inert at default), but it is a kernel change and
needs the usual authorization and full battery.

---

## P3 — the adopted quartic under-taxes small hoards

**Claim.** §4.2: `O ~ x^4` presumes `x ~ T`, true only for a saturated stack.
Below saturation the true exponent is `p_eff = 4/(dm+1)` — as soft as 1 for a
3-dimensional internal stack. The adopted table therefore lets small hoards
survive far longer than the physics allows.

**Apparatus.** No new physics: `p_rad` is already a knob. Run the Face C
condensation ensemble at `p_rad = 4` (adopted) and `p_rad = 2`, everything else
identical, and compare the surviving mass spectra.

**Bar.** At `p_rad = 2` the small-mass tail is suppressed and the spectrum
narrows measurably (interquartile range of `Em_cmax` down ≥25%).

**Kill.** Spectra statistically indistinguishable ⇒ the outtake exponent is not
what shapes the hoard distribution, and §4.2's consequence is void even if the
algebra stands.

---

## P4 — magic numbers appear when a hoard is cooled

**Claim.** §7: the neutral occupancies are the partial capacity sums `Z_j` of
the object's own contact graph, and the mass distribution collapses onto one of
them when `Delta/T >> 1`. Face C's broad spectrum is the hot branch.

**Apparatus.** Offline + existing arms. For each hoard in an ensemble, extract
the contact graph, run `qgraph.Shells` → `Well.MagicNumbers()` to get the
*predicted* peaks **before** looking at the mass histogram. Then compare across
a temperature ladder (cooler = lower drive / higher `q_detune`).

**Bar.** As `T` falls, the `Em` histogram develops peaks within ±1 quantum of
the predicted `Z_j`, at ≥2σ above a **shuffled-graph control** (same degree
sequence, rewired) — the control is what makes this non-circular.

**Kill.** No peaks at any temperature, or peaks that do not move when the graph
changes ⇒ shells are not operating and §5–§7 are decorative.

---

## P5 — coherent intake beats hot intake at matched energy

**Claim.** Theorem 3.5: at matched `⟨Ee⟩`, a phase-matched bath holds a
structure at `x* = u*·xb` with `u* > 1`, while a random-phase bath holds it at
`x* = xb` exactly. Energy is not the variable; coherence is.

**Apparatus.** The energy-matched control discipline `PROGRESSIVE_CELL.md` §4
already established: `kx` (phase wind) vs `k0` vs phase-scrambled, all at
identical `sum Ee`. `amp_door` is the coherent-merge path. Report signed `J_x`,
not `|J|` or global `R2` (both mislead, per that note).

**Bar.** `x*(coherent) / x*(random) ≥ 1.35` at matched `⟨Ee⟩`, sustained ≥500
t.u.

**Kill / discriminator.** Ratio `1.00 ± 0.05` ⇒ the stimulated cross term is
absent or irrelevant, Theorem 3.5 is dead, and the box-fill result is the whole
story. This is the second of the two decisive tests: **"radiance only + noise
fuel" predicts exactly 1.00 here.**

---

## P6 — a wound door pumps winding by exactly `k`; a point door scrambles

**Claim.** L5 and §9: on `C12`, a uniform door allows only `Δm = 0`
(mean 0.000, sd 0.000); a door wound at `k` allows exactly `Δm = k`
(mean +2.000, sd 0.000 at `k=2`); a point-contact door allows every `Δm`
(mean 0.000, sd 3.488). Zero mean is **not** neutrality.

**Apparatus.** Seed a ring at `W = +2` (`seed_mw`), illuminate with (a) a wound
field at `kx = k`, (b) a random-phase field at matched `⟨Ee⟩`, (c) dark control.
Track the ring winding per door event.

**Bar.** Wound arm: `mean ΔW = k ± 0.2`. Random arm: `|mean ΔW| < 0.2` **with
sd > 1** — both halves required, since the mean alone cannot tell a neutral
door from a scrambling one.

**Kill.** Wound arm shows no preferred `ΔW` ⇒ the door is not a Fourier
operator on the contact set, and L5 (hence the Face B explanation) fails.

---

## P7 — charge needs cycles

**Claim.** Theorem 6.3: a hoard whose contact graph is a forest (`b1 = 0`) has
a non-degenerate spectrum, capacity 1 in every shell, no half-filled multiplet,
and therefore cannot retain winding — however much `Em` it holds.

**Apparatus.** Offline on existing runs. Compute `b1 = |E| − |V| + comps` of
each hoard's contact graph per frame; correlate with winding retention.

**Bar.** Winding retention occurs only for `b1 ≥ 1`, and retention time rises
with `b1` (or with the size of the shortest independent cycle).

**Kill.** Retention observed at `b1 = 0`, or no correlation ⇒ the topological
account of charge does not describe what the substrate actually retains.

**Note.** This is already half-measured: `RING_DNLS.md` reports the C6 ring
(`b1 = 1`) holding `W = +2` while packets and chains do not. P7 asks for that
to be turned into a correlation with a computed invariant rather than a
one-object observation.

---

## Ordering

| order | item | cost | why here |
|---|---|---|---|
| 1 | **P1** | offline | can invalidate an adopted constant; cheapest thing in the list |
| 2 | **P7** | offline | re-reads existing frames; sharpens a result already in hand |
| 3 | **P4** | offline + existing arms | pre-registered prediction with a shuffle control |
| 4 | **P5** | config only | first decisive discriminator |
| 5 | **P6** | config only | tests L5 directly; explains Face B or kills it |
| 6 | **P3** | config only (`p_rad` knob) | consequence check on §4.2 |
| 7 | **P2** | **kernel change** | second decisive discriminator; needs authorization + full battery |

Items 1–3 cost nothing but analysis time and can be done against frames that
already exist. Nothing here proposes adopting anything.
