# v85 FULLNESS — Harmonic fullness, the force map, and shells as response harmonics

**Date:** 2026-07-21
**Status:** design/analysis thread under `v85/PROPOSAL.md` — **no runs in this commit**
**Origin:** user threads 2026-07-21: (1) stable particles as harmonically
"full" units — cannot shed below themselves, cannot accept partial additions;
strong force = the full harmonic, weak force = partial dissolve/reform as
symmetric repulsion; (2) electron shells as **secondary response harmonics of
the field itself** — no electron particle; bonds = shared harmonics.

**Epistemic note:** same standard as PROPOSAL — qualitatively satisfying ≠
true; every claim below is bound to a ladder rung (X6–X12) or an explicit
open lock (D11–D13).

---

## 0. The two claims, sharpened

**FULL (fullness):** a stable/metastable particle is a closed process whose
binding capacity is saturated — it cannot transfer *part* of itself away (only
whole-closure units), and it cannot accept an incoming quantum (no unoccupied
harmonic slot). Contact between full objects produces short-range repulsion
(failed-exchange distortion); the "strong force" attraction is closure
exchange between **unfull** clusters; fullness is the hard core and the
saturation.

**SHELL (no-electron shells):** the atomic electron sector is not a light
soliton. It is the **bound harmonic response of the same fabric** to a full
charged nucleus — discrete gap-protected modes in the nucleus's Coulomb well.
Shells are geometry of the response spectrum; neutrality is self-limiting
screening; covalent bonds are two-center shared response modes.

**The interlock:** SHELL is only stable because of FULL — the response cloud
does not drain into the core because the core cannot accept it (saturated
capacity + Δω mismatch). Thread 1 is the protection mechanism for thread 2.

---

## 1. Fullness made precise

### 1.1 The quantitative handle (no new concepts)

On the Q-ball branch, μ(Q) ≡ dE/dQ = ω is the marginal cost of charge.
Fusion of two objects pays iff \(E(Q_1{+}Q_2) < E(Q_1) + E(Q_2)\) — liquid-drop
logic. Definitions:

- **Iron-peak analog** \(Q^\*\): minimizer of \(E/Q\). Fusion is favorable
  climbing toward \(Q^\*\), neutral/unfavorable past it.
- **Fullness measure (candidate, lock D11):**
  \[
  F(Q) \;=\; 1 - \frac{|d^2E/dQ^2|_Q}{|d^2E/dQ^2|_{Q_{\min}}}
  \]
  or equivalently the approach of μ(Q) to its thin-wall asymptote.
  \(F \to 1\): full (no marginal binding gain); \(F\) small: greedy fuser.

### 1.2 Already observed, previously unnamed

| Observation | Fullness reading |
|---|---|
| \(V(s) = (\mu/2)s/(1+\kappa s)\) κ-knee | Attraction saturates at high density — the microscopic fullness mechanism |
| \(Q_{\max} = 921\) (g=0.05) | Global capacity: no closure beyond it |
| v74 **c6_light parks at Q→650** | A droplet that *stopped absorbing* — candidate sighting of \(Q^\*\)/full state |
| v74 **c12_light Q→1411 sheds** | Super-capacity state refusing to exist |
| v71 no isolated components, no mesons | Sub-unit transfer forbidden — unit integrity |

### 1.3 Unit integrity: "cannot transfer away from what it is"

VK stability (dQ/dω < 0) + H5: a stable ball cannot drool charge continuously
because any shed fragment must itself close. **Shedding is quantized in units
of the minimal closure** (bottom of the stable branch, Q ≈ 88 on the v72 VK
branch; the ungauged window's own minimum to be read off the branch data).

**Retrodiction (cheapest test in this document):** the c12_light evaporation
(1411 → sub-critical) should show *lumpy*, quantized shedding steps near the
minimal closure in the existing cluster/diag time series — not continuous
drool. Data is already on disk (X6a).

---

## 2. Force map under fullness

| Real-physics label | Fullness reading | Status |
|---|---|---|
| Strong attraction (medium range) | Closure-compatible whole-unit exchange between **unfull** clusters (fusion pays while μ falls) | Measured — v74 fusion branch |
| Hard core / saturation | **Fullness**: both partners saturated; overlapping evanescent tails cannot commit transfer; failed-exchange distortion energy > 0 → repulsion with Yukawa-like range \(\sim 1/\sqrt{m^2-\omega^2}\) | Predicted — X6b |
| Weak force | Partial **dissolve/reform**: a full object hit by an incompatible quantum de-closes transiently, re-closes (possibly re-partitioned in flavor), shedding a minimal closure; symmetric repulsion during the excursion | Least developed — X9 |
| Neutrino | The **minimal closure**: smallest process that closes at all; aligned with almost nothing → H7-T3 transparency for free | Design target |
| Coulomb | Gauge sector; couples to charge **regardless of fullness** | Measured |

**Structural key for atoms:** Coulomb is fullness-*blind*; fusion is
fullness-*gated*. An atom is a system where the gauge force acts and the
fusion channel is forbidden.

---

## 3. Object-level harmonic gating on the current kernel
*(correction to PROPOSAL §2.1 inventory)*

\(V(s)\) is phase-blind in \(s\) — but \(s\) of a **superposition** is not.
In an overlap region \(\Phi_a = \Phi_a^{(1)} + \Phi_a^{(2)}\):

\[
|\Phi_a|^2 = |\Phi^{(1)}_a|^2 + |\Phi^{(2)}_a|^2
+ 2\,f_1 f_2 \cos\!\big((\omega_1-\omega_2)t + \delta\big)
\]

- **Δω = 0:** static, phase-dependent interaction — same-phase attracts,
  anti-phase repels (standard Q-ball physics; the v83/e1 sign_opp/sign_same
  runs probed exactly this sign structure).
- **Δω ≠ 0:** the cross term beats at Δω and **time-averages toward zero** —
  detuned objects interact weakly at contact despite the pointwise
  phase-blind potential.

So object-level H7-T1 gating is *partially emergent on the unmodified kernel*
via interference cross-terms. The §2.1 verdict "T1 violated at high density"
is softened accordingly (equal-clock dense overlap still couples via density).

---

## 4. Annihilation block: flavor charges as lepton number

The v70 pos1 slow-annihilation is the historical Stage-3 killer. Under
fullness language: annihilation happens because a conjugate pair **has** a
combined closure — vacuum + radiation. Blocking it requires a conserved
quantity conjugation does not cancel.

The fabric has one unused: at η=0 the potential \(s = \Pi_a|\Phi_a|^2\) is
invariant under **independent** phase rotation of each component → the three
flavor charges \(Q_a\) are separately conserved global U(1)s; only total
\(Q\) is gauged. **Three global U(1)s over one gauged U(1) = the
charge / lepton-number split.**

**Design principle (if a light-soliton path is retained at all):** a light
ball whose flavor partition mismatches the nucleus's conjugate cannot fully
annihilate — the residual flavor charge has no carrier state, so the channel
is blocked by conservation, not by fiat. This is the dynamical replacement
for v75 multi-fabric separation, on the existing kernel, using the v71
flavored branch. Tested at X8.

---

## 5. Shells as response harmonics (the no-electron atom)

### 5.1 Ontology

There is no electron object. A parked ball of charge \(Q_N\) sources
\(A_0 \sim g Q_N/4\pi r\). Small charged fluctuations **of the same fabric**
in that well obey a Klein–Gordon-in-Coulomb problem, which has a **discrete
bound spectrum** \(\omega_n < m\) below the mass gap — gap-protected, so
bound response modes cannot radiate away linearly. Shells are the geometry of
this spectrum: effective coupling

\[
\alpha_f = \frac{g^2 Q_N}{4\pi},
\qquad
\varepsilon_n \sim \frac{m\,\alpha_f^2}{2n^2}\ \text{(sub-critical regime)},
\qquad
a_0 \sim \frac{4\pi}{g^2 Q_N\, m}.
\]

The Bohr–Sommerfeld/H6 reading: shell radii are where the response harmonic
closes on itself around the core — integer-cycle closure of a *secondary*
motion, exactly "the field moving in secondary response to the harmonics."

### 5.2 Why the cloud does not fall in — the interlock

Classically an opposite-charge cloud should drain into the core. Two
protections, both from the FULL thread:

1. **Capacity:** the nucleus is full — no unoccupied slot below shedding
   threshold; absorbing the cloud's quanta gains no binding (κ-knee / F→1).
2. **Detuning:** cloud mode frequencies sit just under the gap
   (\(\omega_n \lesssim m\)); nuclear internal harmonics sit elsewhere — the
   §3 beat-averaging suppresses the transfer channel.

If both protections fail in-sim (cloud drains), the SHELL thread dies (F10)
— and that failure would also wound the FULL thread, since it predicts the
protection.

### 5.3 Neutrality is self-limiting (the electron count for free)

As the response cloud grows it screens the core: the well shallows, binding
weakens, uptake slows. The cloud **self-limits at net neutrality** — a
\(Z=6\) nucleus accumulates −6 of response charge because that is where the
binding that holds the cloud disappears. Atomic charge bookkeeping emerges
classically. (Integer *per-mode* occupancy does not — see D13.)

### 5.4 Bonds = shared harmonics; nucleons stay distinct

- **Covalent analog:** two neutral atoms at distance D present a two-center
  well; the two-center response problem has gerade/ungerade shared modes; a
  shared mode lowering total energy = bond. Shared **harmonic**, not shared
  particle — the user's "covalent bond, but not with electrons."
- **Nucleon distinction:** the shell harmonic is a whole-atom property of the
  far field; the nucleus's internal substructure (Stage 2B interlock, flavor
  partition) is untouched by it. Atoms bind through shared response
  harmonics while protons/neutrons retain identity — consistent with
  "harmonics don't fully merge, but are shared."

### 5.5 Honest gaps

| Gap | Statement | Lock |
|---|---|---|
| Scale separation | With g=0.05, \(Q_N\)=650, m=1.5: \(\alpha_f \approx 0.13\), \(a_0 \approx 5\) lattice units < ball radius (~10). Cloud would sit *inside* the nucleus. Need smaller \(g^2Q_N m\) product or explicit redesign before X10/X11 | **D12** |
| Shell filling | Classical cloud is bosonic — nothing enforces 2/8/18 occupancy. Self-consistent screening gives radial stratification (Thomas–Fermi/Hartree-like: filled inner cloud shallows the well, later charge binds in higher modes) but **integer occupancy / Pauli is not derived** — a T2-class gap | **D13** |
| Statistics | No exchange antisymmetry anywhere in a classical scalar cloud; periodic-table structure beyond gross shells is out of scope | D13 |

### 5.6 Stage 3/4 reinterpretation

If X10/X11 pass, the standing-goal ladder changes shape:

- **Stage 3** "light opposite-charge soliton" → **bound response spectrum of
  the existing fabric** (no new soliton species, no multi-fabric kernel).
- **Stage 4 atom** = full nucleus (2A product) + self-limited response cloud
  at neutrality, with X5 emergent orbit quantization applying to the cloud
  modes.
- The soliton-electron path (with §4 flavor protection) remains the fallback
  if F10/F11 kill the response path.

CLAUDE.md stage table update deferred until X10/X11 report (user's call).

---

## 6. Ladder extension X6–X12 (design only; cheap-first)

| Rung | Cost | Test | Kill |
|---|---|---|---|
| **X6a** | analysis only — data on disk | E(Q), E/Q, μ(Q), curvature from existing branch data → freeze F(Q) (D11), locate \(Q^\*\); test c6 park (650) against \(Q^\*\) and \(Q_{\max}\)=921. **Plus:** reanalyze c12_light shedding time series for quantized lumps near the minimal closure | **F6**: shedding is continuous drool → unit integrity fails |
| **X6b** | short runs | Collision matrix across the branch: merge/no-merge boundary vs F(Q); near-critical pairs (Q~800 each) predicted to shed-and-bounce, not merge; contact-repulsion range vs \(1/\sqrt{m^2-\omega^2}\) | **F7**: boundary uncorrelated with branch curvature → fullness ≠ saturation |
| **X7** | short runs | Δω gating scan: pair interaction strength vs detuning; phase sweep at Δω=0 (leverage v83/e1 sign runs as priors) | **F8**: no beat suppression → object-level T1 absent |
| **X8** | medium | Flavor-blocked annihilation: conjugate pair with mismatched \(Q_a\) partitions; predict stall with conserved residual flavor charge | **F9**: annihilation completes despite mismatch → flavor block fails; soliton-electron path weakened |
| **X9** | medium (after X6) | Weak analog: overdrive a full ball past capacity; look for de-close/re-close with flavor re-partition + shed minimal closure; measure transient symmetric repulsion | (exploratory — informs, does not kill) |
| **X10** | semi-analytic + ringdown run | Response spectroscopy: solve KG-Coulomb around a parked ball (radial solver); then in-sim seed small opposite-charge perturbation, measure ringdown frequencies below the gap vs prediction. Requires D12 parameter design first | **F11**: no discrete bound modes below gap → SHELL dies |
| **X11** | medium | Capture & self-limiting: bathe a full parked ball in opposite-charge condensate (AD-style); measure cloud formation, radial mode structure, uptake self-limiting at neutrality, and core integrity | **F10**: cloud drains into core → SHELL dies, FULL protection wounded |
| **X12** | large (after X10/X11) | Covalent analog: two neutral screened atoms at distance D; binding curve E(D), gerade/ungerade splitting of the shared mode | (payoff rung — molecule without electrons) |

Order: X6a → X7 → X6b → X10-design(D12) → X8 → X11 → X10-run → X9 → X12.
X6a costs nothing but analysis time and constrains everything downstream.

---

## 7. Falsifier registry (extends PROPOSAL §7)

| # | Observation | Kills |
|---|---|---|
| F6 | c12 shedding continuous, no quantized steps (X6a) | Unit integrity / minimal-closure shedding |
| F7 | Merge boundary uncorrelated with F(Q) (X6b) | Fullness = saturation identification |
| F8 | No Δω suppression of pair interaction (X7) | Object-level T1 on current kernel (§3) |
| F9 | Annihilation completes despite flavor mismatch (X8) | §4 flavor block; lepton-number analog |
| F10 | Response cloud drains into full core (X11) | SHELL thread; wounds FULL protection claim |
| F11 | No discrete bound ringdown spectrum below gap (X10) | SHELL thread (response ontology) |

FULL and SHELL are separable: F10/F11 kill SHELL while leaving FULL standing
on X6–X8; F6/F7 kill FULL's saturation identification while SHELL could
survive on detuning protection alone (weakened).

---

## 8. Files

| Path | Role |
|---|---|
| `v85/FULLNESS.md` | This document |
| `v85/PROPOSAL.md` | Parent proposal (H-axioms, X1–X5, F1–F5, D1–D10) |
| `v74/RESULTS.md` | c6 park / c12 shed — the existing fullness sightings |
| `v83/e1/` | sign_opp/sign_same — priors for X7 phase structure |

**No simulation runs in this commit.**
