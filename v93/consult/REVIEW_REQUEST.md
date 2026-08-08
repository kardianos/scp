# v93 REVIEW REQUEST — theory, code, results, conclusions

You are reviewing a physics-simulation research programme ("SCP"). The active
version is v93. Please review the THEORY, the CODE, the RESULTS, and the
CONCLUSIONS, and suggest: (a) fixes in method or code, (b) missed aspects of
the results, (c) concrete ways forward. Be specific and critical.

## Files to read (all under /home/d/code/scp/)

- `v93/README.md` — the full design (PART I background, PART II the design,
  PART III plan, PART IV safeguards). READ THIS FIRST, top to bottom.
- `v93/L1_FINDINGS.md` — the experimental results (faces 1-3, L1-B/L1-C,
  QUENCH-3, arg(ψ) door). READ THIS SECOND.
- `v93/kernel/freecell.c` — the C kernel of record. Key regions:
  - pass U (the unitary dense hop): search "pass U: v93 UNITARY DENSE CHANNEL"
  - the hop-angle branch in pass 2: search "v93 UNITARY DENSE CHANNEL: fold"
  - the in-pass precession: search "Local clock precession IN-pass"
  - the arg(ψ) door: search "v93 arg(psi) door (§II.7)"
  - the field pass F (the template): search "pass F: unitary field hops"
- `v93/fab/step.go` — the Go A/B kernel mirror (same logic).
- `v92/consult/SUBQUARK_synthesis.md` — the 3-reviewer convergence this design
  came from (for context on what was already considered).

## The one-paragraph thesis (v93)

Give the DENSE (matter) sector the FIELD sector's transport algebra. Within-mode
dense transport becomes a product of UNITARY PAIRWISE PLANE ROTATIONS (a cousin
of the field's pass F) on the dense amplitude ψ_m = √Em·e^{i·th2}, replacing
the ADDITIVE MAGNITUDE "want" (Em[src]-=d1; Em[rcv]+=d1, phase used only to
gate). Each Givens hop conserves the two-cell norm EXACTLY, so ΣEm is conserved
to roundoff by construction (conservation becomes a THEOREM of the update, not a
patched ledger). The cross-term 2·Im(ψ_i*·ψ_j) that the additive Em-ledger
rejected IS the link current J_s = the dense momentum. Knob `amp_nat` (default
0 = byte-inert). The conversion DOOR (atoms ε=A₀ω/2π, cap, radiance) stays
discrete/quantized (NOT unitarized).

## Implementation as built (3 faces)

- **Face 1:** τ_s = amp_nat·base·g_ij_sym·√(head_i·head_j). FAILED — the blob
  loads cells to ~cap ⇒ head≈0 in the dense core ⇒ the coherent tilt is frozen.
- **Face 2:** drop head: τ_s = amp_nat·base·√(g_ij·g_ji). ENGAGED (coherent
  direction, seed-robust within an amp_nat) but the blob wobbled/melted.
- **Face 3 (current):** move the dense clock precession INTO pass U (rotate ψ_m
  by w2e·dt BEFORE the hops, mirroring pass F's precess-then-hop); skip the
  out-of-pass th2 advance in pass 6 when amp_nat>0. → BREAKTHROUGH.

The hop (pass U): reconstruct (m1,m2)=√Em·(cos th2, sin th2); precess by w2e·dt;
per slot (canonical i-side) Givens-rotate by τ_s:
  m1[i]'=cc m1[i]+ss m2[j]; m2[i]'=cc m2[i]-ss m1[j];
  m1[j]'=cc m1[j]+ss m2[i]; m2[j]'=cc m2[j]-ss m1[i]
write back Em=|ψ|², th2=arg(ψ). (Identical in form to pass F.)

**arg(ψ) door (latest face):** knob `amp_door`. At condensation (field→matter,
d1>0), instead of qp_phase's partial pull `th2 += mix·(aph−th2)`, do a COHERENT
amplitude merge: arg(ψ_m_new) = arg(√Em_old·e^{i th2} + √d1·e^{i aph}) where
aph=atan2(fa2,fa1); |ψ_m| fixed by the (already-conserved) Em.

## Results (all behind amp_nat/amp_door, byte-inert at 0, C≡Go, full battery
## ALL GREEN 87 at amp_nat=0)

- **L1-B conservation: GREEN** — drift_rel = 0.000e+00 at amp_nat 0/2.5/5.
- **L1-C anti-ignition: GREEN** — quiescent bath BYTE-IDENTICAL across amp_nat
  0-5 (births 4611, z_live 16.89). The channel is SELF-GATING: random bath
  phases ⇒ √(g_ij·g_ji)≈0 ⇒ τ_s≈0 ⇒ inert on incoherent matter.
- **L1-A translation: NEAR** — e3b (tilted blob). At amp_nat=2.0, ALL 5 seeds
  translate +x (cos 0.76-0.995, two seeds ≈0.99; seed 111 also clears speed
  0.0032>2.6e-3 → meets L1-A). But: NARROW RESONANCE (cos peaks at amp_nat≈2,
  falls off by 2.5; direction reverses between windows — a dispersion signature);
  seed-variant magnitude (0.76-0.995); speed mostly sub-threshold. blob melts
  ~3.6% over 42 t.u. (free-packet diffraction). bath=0 → zero motion.
- **QUENCH-3 spin retention: the door face.** Standing vortex (spin_m=2 field
  winding, m=+2 about axis), door carries phase, measure matter winding R2d.
  Steady-state (t=20-300): R2d≈0.02-0.03 for ALL arms (no door / qp_phase /
  arg(ψ) door / +unitary / +uniform-clock). BUT early-time (t=6-14): the arg(ψ)
  door SUSTAINS R2d~0.1 vs ~0.02 for no-door (~3× longer retention) before
  decaying to ~0.02 by t=20. The field winding itself decoheres in transit
  (RA2 1.0→~0.15).

## My conclusions

1. The v93 thesis is CONFIRMED where structural: coherent dense translation
   emerges (L1-A near, a first), conservation dissolves into a theorem (L1-B),
   anti-ignition is structural (L1-C self-gating).
2. The arg(ψ) door helps retention (slows decay ~3×) but does NOT achieve
   steady-state retention — the wall is the MATTER CHURN (births~deaths~9000
   over T=300, matter turns over) + field-winding transit decoherence, not the
   door mechanism.
3. L1-A's remaining gaps: narrow resonance (dispersion), seed-variant
   magnitude, free-packet diffraction (no binding).

## What I want from you

Critically examine, with specific reference to the code/results:

1. **Theory/code soundness:** Is pass U correct? Is the in-pass-precession
   ordering right? Is the arg(ψ) door's coherent-merge renormalization
   legitimate (energy-conserving)? Any bug or sign error in the Givens hop or
   the current J=2·τ·Im(ψ_i*ψ_j)?
2. **The narrow resonance (L1-A):** Is the amp_nat≈2 window + direction
   reversal a real dispersion relation, or an artifact (e.g. of the e3b
   centroid-fit measurement, or of τ_s coupling to the phase via √(g_ij·g_ji))?
   How would you tell? How would you WIDEN it / pin magnitude seed-robustly?
3. **The QUENCH-3 retention wall:** Is "matter churn" the right diagnosis, or
   is the imprint itself failing / the field winding the real limit? What
   experiment would distinguish these? Is there a code/method fix that could
   yield real steady-state retention?
4. **Packet binding:** the free unitary packet diffracts (no soliton). Is a
   binding companion necessary for a translating blob, or should L1-A be
   reframed (e.g. measure wave-packet group velocity, not a bound blob)?
5. **Missed aspects / alternative explanations** for any of the results.
6. **Concrete next experiments or code changes**, ranked.

Please be concrete and cite file:line where relevant. Disagreement welcome.
