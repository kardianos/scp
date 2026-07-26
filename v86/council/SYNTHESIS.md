# v86 council SYNTHESIS — three seats + foundational analysis
**Date:** 2026-07-26 · **Seats:** Claude Fable (resumed max session, 19
findings), Grok 4.5 (review + FOUNDATIONS follow-up), Kimi K3 (11 fixes,
first genuine K3 seat — v84's attempt fell back to glm). All four documents
in `v86/council/*/`.

## 1. Verdict matrix
| Document | Fable | Grok | Kimi |
|---|---|---|---|
| GROUNDING.md | SOUND-WITH-FIXES | SOUND-WITH-FIXES | SOUND-WITH-FIXES (5) |
| PROPOSAL.md | SOUND-WITH-FIXES | SOUND-WITH-FIXES | SOUND-WITH-FIXES (6) |

## 2. Converged findings (≥2 seats independently)
1. **GSS sign/counting error (Fable F1 ≡ Grok A1) — CRITICAL.** The stated
   criterion "stability when n(H_ω)=p(D), D=∂Q_a/∂ω_b" fails its own
   single-charge reduction (VK-stable branch: dQ/dω<0 ⇒ p(D)=0≠1=n(H)).
   Correct object: the **negative** index of ∂Q_a/∂ω_b must match n(H_ω)−…
   (see corrected GROUNDING §1). Had HC-3 been coded as written, every
   partition classification would have been inverted.
2. **Gauss projection hand-waved (F3 / A2):** the cited vortex literature
   lives in the Higgs phase; this system is Coulomb-phase with a massless
   channel reaching ω=0 — GSS is a theorem ungauged, a heuristic gauged.
   HC-3 verdicts must carry that caveat.
3. **The gauge channel voids naive gap arithmetic for l≥1 (F4 / Grok A3 /
   Kimi "multipole audit"):** any mode with a time-varying multipole
   radiates through massless A at ANY frequency. The census must classify
   modes by multipole first; arithmetic protection applies only to
   monopole/matter channels.
4. **HC-4 is the misleading-result rung (F6 / Grok B3 / Kimi width-floor):**
   the QRK-1 lines (0.018–0.126) sit BELOW the box IR cutoff (0.209 CPU /
   0.057 GPU): no open decay channel in any proposed box, while nonlinear
   frequency-pulling fakes the "confirming" ε² scaling. Hardening: width-
   floor calibration on a mode with a genuinely open channel + box-size
   pair to move the cutoff through a line.
5. **EX-1 threshold wrong by the mass factor (F7; Grok A5 partial; Kimi
   5,7–10):** correct condition v ≪ √(2ε₁/m) = α_f ≈ 0.053 (not 0.065);
   v=0.05 is 89% of threshold (predicted ~55% sudden-kick stripping), and
   Kimi adds that differential drift strips even below the bound — ramped
   boost required at the top end; EX-4 retargeted to the ball clock.

## 3. Kimi-specific adds (adopted)
- **HC-6 (converse rung):** deliberately UNSTABLE partition seeds must
  decay to their sector minima — the census needs the converse test, not
  only survivals.
- Protected-mode width-floor calibration (numeric noise floor measured on
  a mode with no open channel, so "zero width" is a calibrated statement).

## 4. FOUNDATIONS (Grok follow-up on the E=mc²→E=Qω arc) — adopted frame
**Verdict:** the v85 elevation of E=Qω repeated v84's mistake in mirror
image — CREED's own caution ("a closed formula is dualist residue unless
derived") was discarded to close the D5 units gap cosmetically. The
surviving coherent ontology (R1, adopted):
- **E** (Hamiltonian, vacuum-subtracted) is primary; **M = E/c²** — E=mc²
  returns honestly as the rest-energy relation, not a "derived shadow."
- **ω_a = ∂E/∂Q_a** are chemical potentials (Legendre exact, measured).
- **ε is not mysterious (R2):** Σ ≡ E − ωQ is the virial/surface excess
  (gradient + potential + gauge − kinetic), with three concrete closed-form
  targets: the integrated virial identity of the profile ODE (parameter-
  free), a thin-wall surface law Σ ≈ 4πR²σ_eff + g²Q²/8πR, and a scaling
  collapse ε vs ξ/R. "No closed form for ε" was an unfinished calculation,
  not an ontological failure.
- **Three-way split (R4):** charge Q_G ≠ process action 𝒜 ≠ inertia M even
  within one sector; "ħ_eff = Q" demoted from identity to measured ratio
  (v70's 3–5% and ε's 1–4% are the same residual family).
- **Degeneracy attacks without kernel surgery (R5):** ring spin-ladder
  M(n) at fixed Q (mass without charge — v73 rings as the lab), window
  widening via potential-parameter scans; cloud-mass explicitly does NOT
  solve HIER (restates it). Second sector remains motivated by carrier
  physics + universal ħ, not by rescuing E=Qω.
- **Numeric battery N1–N10** (mostly free): ε decomposition on existing
  gscan/branch data, virial residual, scaling collapse, flavored Legendre
  excess, throughput-vs-charge on breathing/flavored objects, ħ_eff triple
  comparison, N-INV weak-acceleration inertia test (ends D5′ empirically),
  ring mass-ladder, soft-window shooter scan, shell-mode E=ωQ exactness on
  X10 profiles.

## 5. Program consequence
The v86 order of work gains a Part 0: the **ε/foundation battery (N1–N4,
N10 first — all shooter/analysis, zero spend)** before the harmonic census;
HC rungs proceed with the corrected GSS object, the multipole-first
classification, HC-6, and the width-floor calibration; EX-1 runs ramped
with the corrected threshold. TOE Step 3 and §D get rewritten to the R1/R4
frame at the next TOE revision.
