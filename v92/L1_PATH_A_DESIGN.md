# L-1 PATH A — the native complex-transport design sub-round (v92)

User-directed ("design sub-round first, then implement"). This specifies
the native complex-transport completion BEFORE any kernel edit, derived
from the L-1 findings (feedback insufficient) and the S2 memo. It is a
design proposal for review, not code.

## 0. Why feedback failed, and what "native" must mean

Both L-1 forms (magnitude-bias, phase-current) are a FEEDBACK: the
magnitude want `swant` is computed as today (gate·head·mob), then a
shadow-driven BIAS perturbs it. The e3b and the chord probe proved this
cannot work: the magnitude transport carries no coherent phase, so its
shadow is incoherent, so the bias amplifies noise (L1) or finds nothing
to sustain (L2 — the dead fifth has no deposits, hence no shadow).

**Native** means: the dense transfer IS a complex amplitude whose phase
composes at the receiver, and that composition is the BASIS of the
outgoing transport — not a correction on a magnitude computed without it.
The phase must propagate through a deposit→reception→re-emission chain:
that chain IS the current.

## 1. The three constraints (non-negotiable, from the evidence)

**C-1 Conservation.** Energy is the magnitude. The debit/credit ledger
(Em[src] −= d1, Em[rcv] += f) is the invariant — it must not change. The
phase is carried ALONGSIDE the magnitude (informational), never creating
or destroying energy. The atoms machinery, the conversion door, the 87 V3a
bars are the invariant surface.

**C-2 Protected carriage (QUENCH-3 §7.7).** The phase must NOT be written
into the cell clock th2 — the naked clock churns (~0.2 rad/t.u.
differential rotation; qp_phase's retention was zero). It must be
slot/account-borne (the shadow's home), where delivery coherence is
measured 0.77–0.95. The matter clock th2 stays the matter clock.

**C-3 Anti-ignition (the ignition ledger).** The bath's transport is
incoherent (random arrival phases ⇒ the composed phase averages to zero ⇒
no coherent current). The completion must be starved in the bath BY
CONSTRUCTION (as the L-1 feedback was — byte-identical bath), not by a
statistical gate.

## 2. The complex amplitude (not a new field)

The dense cell already HAS a magnitude (Em) and a phase (th2). The complex
amplitude is the combination the sector already owns:

    ψ_i = √Em_i · e^{i·th2_i}        (the dense "wave function")

This is NOT an imported field — it is √Em and th2, both already present.
Path A makes the TRANSPORT compose ψ (phase-coherent), where today it only
GATES on th2 (magnitude transport, phase used as a closure condition).

## 3. The mechanism: deposit → coherent reception → re-emission

The transport chain becomes phase-coherent end to end:

**(a) The want carries phase (pass 2).** Today `swant[2s+dir]` is a real
magnitude w_ij. Path A: the slot also carries the departure phase
    φ_ij = q·th2_i − q·w2e_i·d/C           (chart-folded retarded source)
— the SAME quantity that is the gate argument today (freecell.c:1335), and
the SAME phase the Phase-M shadow records. So the departing amplitude is
    A_ij = w_ij · e^{i·φ_ij}                (magnitude w_ij, phase φ_ij)
The magnitude w_ij is unchanged (C-1: energy = w_ij, debit/credit as today).
The gate still gates (cos^p closure on φ_ij − p·th2_j); the phase is
ADDITIONALLY carried.

**(b) Coherent reception (pass 4 deposit).** At the receiver j, the
arriving amplitudes compose into a per-cell RECEIVED account (slot-borne,
C-2):
    R_j  ← Σ_{k→j} √w_kj · e^{i·(φ_kj − p·th2_j)}    (rotated into j's frame)
This is exactly the Phase-M shadow's delivery composition (freecell.c:1522)
— PROMOTED from meter to a state that drives emission. The coherent
magnitude |R_j| vs Σ√w (ρ_coh) measures whether j's arrivals reinforce.

**(c) Re-emission (next step's pass 2).** j's outgoing transport is biased
toward the coherent direction by ITS OWN received composition: where |R_j|
is large (coherent reception) and aligned along a link, that link's gate
g_jk is boosted by the coherent factor. The outgoing phase carries arg(R_j)
— the current propagates. **This is the native current:** coherent
reception → coherent emission → directional flow. It is not a bias on a
magnitude computed without R; R is the basis of the emission.

**Why this is native, not feedback.** In L-1, the shadow biased a want that
was already fully determined by gate·head·mob (the shadow was external). In
Path A, the emission gate READS R_j (the composed reception) — the
transport is coherently gated by what was coherently received. The chain is
closed (deposit→R→emit); phase composes along paths.

## 4. Why each constraint holds

- **C-1:** the magnitude w_ij is debited/credited exactly as today; R_j
  carries phase only. Drift untouched by construction (the energy ledger is
  the same). The emission boost is a REDISTRIBUTION (coherent links
  strengthened, incoherent weakened) within j's outflow budget — the
  outflow limiter (pass 3) caps the total at 0.98·Em, conserving energy.
- **C-2:** R_j is a per-cell received-amplitude account (slot/account-
  borne), NOT th2. th2 advances by w·dt as today, untouched. The phase
  composes in R_j, protected from the cell-clock churn that killed qp_phase.
- **C-3:** in the bath, arrival phases are random ⇒ R_j = Σ√w·e^{i·random}
  ⇒ |R_j| ≈ √(Σw)·(1/√N) (incoherent sum, small) ⇒ no coherent boost ⇒
  the bath's transport is ~unchanged (the L-1 anti-ignition result, now
  native). Starved by the random-phase structure, not a gate.

## 5. Knob, defaults, acceptance

**Knob.** `amp_coh` (default 0 = byte-inert; the reception→emission
coupling strength). R_j is accumulated whenever `amp_tau>0` (the existing
shadow window); the emission boost fires only when `amp_coh>0`. Reuses the
shadow apparatus (already in both kernels) — the change is that R drives
emission (pass 2), not just prints.

**Byte-inertness.** `amp_coh=0` ⇒ the emission boost is skipped ⇒ the want
is the V2g/V3a magnitude want ⇒ byte-identical (header line only).

**Acceptance (mapped to the now-working measurements).**
- **L1-A (transport):** e3b |cos_to_kdir| → 1, speed ≥ 2.6e-3, seed-robust
  (the phase-current that the feedback couldn't manufacture). The
  deposit→R→emit chain carries the tilt's phase gradient as a current.
- **L2 (chord):** nv=6 boundary fifth gg ≥ 0.9 at T=100, D-ring x̄ ≤ 0.35
  (the fifth's circulation sustained by coherent re-emission, not revived
  from a dead shadow).
- **L1-C (anti-ignition):** bath byte-identical (random phases starve R).
- **Invariant surface:** the 87 V3a bars; drift floor; C≡Go.

## 6. Open questions (for the design review, before implementation)

1. **The emission coupling form.** The cleanest first form: boost the gate
   `g_jk *= (1 + amp_coh·ρ_coh_j·alignment_jk)` where ρ_coh_j = |R_j|/Σ√w
   (receiver coherence) and alignment_jk = Re(R_j·e^{-i·φ_jk})/|R_j| (how
   well j's received phase aligns with the k-link's departure phase). This
   boosts emission along the coherent direction. Alternative: a phase-lock
   pull on the outgoing phase (closer to the cantus form). Recommend the
   gate-boost form first (it is a redistribution, naturally conserved).
2. **R_j's frame.** R_j is composed in j's clock frame (rotate by p·th2_j
   at reception, as the shadow does). This makes R_j a standing
   (clock-relative) coherence — robust to j's common-mode clock rotation.
   Confirm the rotation matches the shadow's (freecell.c:1490).
3. **The bath-reward edge case.** ρ_coh in the bath is ~0.77 (Phase M
   measured it, not zero). A ρ_coh boost could reward bath unison. Mitigation:
   the alignment_jk factor (the directional coherence) is ~0 for the bath's
   isotropic churn even if ρ_coh>0 — test L1-C; if it ignites, gate on
   alignment alone (drop the ρ_coh prefactor). This is the registered
   contingency.
4. **C≡Go at strength.** The feedback amplifies C/Go FP differences in
   meter lines (L-1 finding). A native coupling may amplify less (it gates
   rather than biases) or more. Verify; broaden the abx masking or move to
   tolerance-based comparison if needed.

## 7. Implementation order (after sign-off)

Land `amp_coh` + the reception→emission gate-boost in C, verify byte-
inertness + the e3b (now working), then mirror in Go + C≡Go, then the full
acceptance (L1-A e3b, L2 fifth, L1-C bath, invariant surface). The shadow
apparatus is reused; the change is one coupling in pass 2 (emission reads
R) + R already accumulated in pass 4 (the shadow). Smaller than it sounds —
the Phase-M apparatus did the heavy lifting; Path A closes the loop it
left open.
