# HBAR — the situation at fork F-ħ, and what could be missing

**Brainstorm document.** Subordinate to `PRINCIPLE.md`; companion to
`ROADMAP.md` (fork F-ħ) and `CONSONANCE.md`. Written 2026-07-27 at the
user's direction: *what fundamental thing are we missing?* Candidates are
tagged [P] postulate-grade, [G] guess-grade; nothing here is implemented.

---

## 1. The situation

Measured (E6 round 3, 42 locked pairs of the same species at the same
nominal frequency): per-cycle transferred action `act` spans
**0.084 – 1.584 — a 19× spread**, median 0.617. The kernel's transfers are
rate-proportional (dE ∝ E·dt): infinitely divisible, amplitude-scaled.
Nothing in the dynamics prefers any particular amount per cycle.

Reality is otherwise. A mode of frequency ω carries energy in steps of
**E = ħω**, with one universal slope, and this single fact anchors:

* the **photoelectric threshold** — below-threshold light ejects nothing
  at any intensity; above it, even feeble light ejects promptly;
* the **blackbody spectrum** — the UV catastrophe is cured exactly by
  ω-proportional quanta plus detailed balance;
* **shot noise** — detection variance follows Poisson counting of grains
  of size ħω;
* **atomic stability and level spacing** — orbits quantized by action;
* our own **tier-1 clicks** — whose grain size is currently a knob
  (`e_click`), which is the confession that ħ is not yet physics here.

The kernel has exactly one quantizing axiom — *a partial cycle converts
nothing* — and it quantizes **when** transfer happens, never **how much**.
That is the hole, stated in one line: **the cycle gate quantizes time
structure; nothing quantizes amount.**

## 2. The reframe that makes ħ natural: the atom of transfer is action

If every complete cycle hands over a fixed *action* A₀ — not a fixed
energy — then a channel cycling at frequency ω transfers

```
E_per_cycle = A₀ · ω / 2π      ⇒      E = ħ_eff · ω,   ħ_eff = A₀/2π.
```

Planck's relation is then not a law about energy at all; it is the
statement that **the tail-call frame has a fixed size**. Every handoff
passes one frame; faster clocks pass more energy only because their frames
are denser in time. This sits flush with the framework's own language:
conversion proceeds through complete cycles (PRINCIPLE §3), delivery is
atomic (CONSTRUCTION R3), and the user's tail-call metaphor was always a
fixed-frame metaphor — a function call passes a frame, not a fluid.

Notably, CONSTRUCTION §1.3 **R1 already postulates a deposit quantum
δ > 0** which the kernel quietly approximated away with continuous rates.
The refinement this document argues for: δ is not a fixed *energy* (that
would give E = n·δ, frequency-independent — wrong); **δ is fixed action**,
δ(ω) = A₀ω/2π.

## 3. Candidate origins of A₀ (the brainstorm)

The user's standing directive: *no Planck posited — the limit must be
computed.* So the real question is what computes A₀.

**(1) Integer harmonic content, restored [P — strongest candidate].**
CONSTRUCTION §1.1 already asserts `h(v) ∈ Z^k`: internal content is
*integer* because the internal space is compact ("free motion on a compact
space has a discrete spectrum — mode indices are integers by topology").
The species enumerator honors this; the dynamical kernel does not — we
**de-quantized the construction's own ontology** when we built continuous
amplitudes and energies. Restoring integer occupation per internal mode,
with transfer moving whole occupation steps, gives E = n·ε(ω) directly,
and ε(ω) ∝ ω if a step is one action atom. On this reading F-ħ is not a
missing mechanism but a **fidelity bug**: the kernel is a mean-field
approximation of an integer theory, and `act`'s 19× spread is exactly
what mean-field does to a quantized ladder.

**(2) Phase-space resolution — ħ as the fabric's resolving power [P].**
Tongues have finite width (Γ, p_gate): configurations closer than the
resonance acceptance are *dynamically indistinguishable*. A finite
distinguishability scale in (amplitude × phase) space is a minimal
phase-space cell — and a minimal phase-space cell **is** an action atom.
"The harmonic limit which is computed": A₀ computed from gate sharpness
and resonance width. Consequence to test: ħ_eff would be structural
(same everywhere the gates are), matching reality's universality —
but it must come out *independent* of which channel measures it, which
is a sharp internal consistency requirement.

**(3) The space-cell grain — ħ from the vacuum's own atom [P].**
The fabric already owns an action scale with no new constant:
one cell's space energy transacted over one link crossing,

```
A₀ ~ e_s0 · d̄ / C  =  1.0 × 1.15 / 1  =  1.15   (sim units).
```

The measured median per-cycle action is **0.617 ≈ half of it** — the
right order, suggestively ~A₀/2. What is missing is not the *scale* (the
foam already sets it) but the **locking**: the distribution sits at the
foam scale and spreads 19× around it; a quantizing mechanism would
collapse the spread, not move the center. This candidate predicts
ħ_eff ∝ e_s0·d̄/C — i.e., Planck's constant is a property of the vacuum
cell, and would covary with curvature (E_s depletion) — a testable, and
possibly falsifying, gravitational side-effect.

**(4) Zero-point self-consistency [G].** Reality's modes carry ħω/2 even
empty. If the fabric has an *irreducible structural churn* (not thermal —
built into the tumble/gate floor), the transfer quantum could be the
fixed point of a consistency loop: quanta must exceed the vacuum grain to
be distinguishable; the vacuum grain is set by spontaneous emission into
the churn; each defines the other. Rich, circular, dangerous — kept as a
guess.

**(5) Chirality/winding [G].** The plane-pair rotation is topologically
integer per cycle (winding). If amplitude could only change where winding
changes — transfer as integer winding exchange — quantization would be
topological. Close cousin of (1); the ± of u± would then also quantize
angular momentum in the same stroke.

**(rejected) Thermal noise floor as ħ.** A noise-floor quantum would make
ħ temperature-dependent. Reality's ħ is not. Rejected unless the noise in
question is the *structural* churn of (4).

## 4. Discriminating tests (one constant, five phenomena)

Whatever mechanism is implemented must produce a single ħ_eff that
simultaneously fits:

1. **Linearity:** grain size vs ω across modes — a straight line through
   zero (the current kernel fails this by construction).
2. **Photoelectric analog:** sub-threshold ω cannot condense a grain at
   any intensity; above threshold, prompt condensation at low intensity.
3. **Blackbody analog:** the churn bath's equilibrium spectrum — with
   ω-proportional quanta and detailed balance, occupation should go
   Bose; the spectrum should bend at the same ħ_eff.
4. **Shot noise:** tier-1 click statistics Poisson with grain ħ_eff·ω —
   replacing the `e_click` knob with physics.
5. **Cross-sector consistency:** the dense-sector lock action (`act`)
   and the field-sector grain must agree — the same A₀ on both sides of
   the melody/harmony divide.

Candidates (2) and (3) make additionally falsifiable side-predictions:
(2) ties ħ to gate structure (must come out channel-independent);
(3) ties ħ to the vacuum cell (must covary with E_s — gravitational
running of ħ, which reality constrains hard).

## 5. Recommended route

Implement **(1) + the action-atom reading of R1** together, dense sector
first: integer occupation per mode with transfer in whole action atoms
δ(ω) = A₀ω/2π, A₀ provisionally the foam grain (3) so no new constant
enters; then run the discriminating battery (§4) — beginning with the
collapse of `act`'s 19× spread to a line in ω, then the photoelectric
analog. This is a major dynamical change (integrate-and-fire transfer in
place of rate flow) touching every dense-sector result; it must be its
own campaign with full regressions, not a patch.

---

## 6. First implementation: candidate (1)+(3) tried, same day

Implemented behind a toggle (`quant_A0`; 0 = off, −1 = auto A₀ = e_s0·d̄/C
— no new constant): dense deposits become integrate-and-fire in whole
atoms ε(ω) = A₀ω/2π (credit accumulates, energy moves only in affordable
whole atoms, unfired credit lapses at 2 atoms); beat conversions floored
to whole atoms. Field sector untouched.

**Regression verdict (the question asked): none.**
* Toggle off: byte-exact `# RESULT` agreement on all eight current-era
  baselines (E1, E3a, E5, E6, D1, Tonomura, eraser, HOM — including every
  code path the change touches); the six mismatches were stale pre-repair
  baselines (E2/E3b/E4/E7/E8/E9 last run before the amplitude-field
  repair), confirmed by determinism re-runs and refreshed as the new
  unified `reg_*_R` baseline set.
* Toggle on: the field-sector experiments are **bit-identical** (fringes
  0.466, HOM 0.4342/0.5658, Bell 2.8263) — quantization is dense-only by
  design and leaks nowhere. Conservation ≤1.1e−15 in every quantized run.

**What the quantized universe did (reg_*_Q):**
* **The photoelectric threshold appeared unprompted:** in E1 the largest
  conversion demand (0.10) sits below one atom (ε = 0.27), so
  condensation stops entirely, at any intensity — the sub-quantum freeze
  this document predicted in §4.2. Dense mode now forms only where a
  whole atom's worth of demand exists.
* Static structures got *more* static (E3a blob speed 0.0004 → 0.00016 —
  sub-quantum halo trickles frozen).
* Flow physics throttled at few-quantum occupancy: blob drift 0.0088 →
  0.0012; standing-pair circulation starved (E6 in-tongue gg 0.72 → 0.16,
  flight mostly zero atoms) — pairs freeze rather than lock-and-circulate.
  The continuous-limit dense results live at 2–3 atoms per store, where
  lumpy classical transport cannot sustain coherent exchange.

**Reading, honestly:** the mechanism works, conserves exactly, and
produces the threshold physics ħ exists to produce — but it exposes the
same structural lesson the double slit taught the field sector: coherent
few-quantum behavior (locks at 2–3 atoms) needs *amplitudes*, not gated
lumps of classical energy. The dense sector's quantum-regime treatment —
amplitude-carried locks with integer occupation — is the successor
construction. Until then `quant_A0` stays **off by default**: a validated
opt-in mode, with the ħ-linearity built in (ε ∝ ω by construction) and
its behavioral consequences mapped.
