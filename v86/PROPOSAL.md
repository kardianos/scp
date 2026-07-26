# v86 PROPOSAL — Harmonic census & exchange physics
**Date:** 2026-07-26 · **Status:** proposal (design + analysis first; runs
gated per rung) · **Prior:** v85 closed (STATE_OF_THE_THEORY.md is the map).
**Tooling prerequisites: DONE** (commit 6d8b267): chimera-download fix,
per-component Q_a diag (CPU verified, CUDA at next build), volview -gamma.

## Part I — The harmonic census: how many harmonics can a particle hold?

### I.1 The theoretical bodies that apply
1. **Conserved-sector stability (VK generalized).** The fabric has a lattice
   of conserved quantities: three flavor charges Q_a, ring winding n, (and
   fiber charge when gauge2 is on). A configuration is absolutely stable iff
   it minimizes E within its sector. Multi-charge stability = positivity of
   the Hessian of E(Q_0,Q_1,Q_2) on the flavored branch — a computable
   surface. Prediction: **the number of independently stable harmonics =
   the dimension of the conserved-charge lattice**; everything else decays.
2. **Resonance arithmetic (Soffer–Weinstein / Fermi golden rule for PDEs).**
   Internal modes above the clock survive only if every low-order integer
   combination n·ω_i ± m·ω_clock stays below the radiation continuum (gap
   m=1.5). Modes whose combinations cross the gap decay at computable rates
   (∝ coupling²·density-of-states). The census of long-lived harmonics is
   **gap arithmetic**, not dynamics.
3. **KAM / quasi-periodic persistence (Kuksin–Wayne).** Multi-frequency
   invariant tori in Hamiltonian PDEs persist only on thin parameter sets
   and are destroyed by coupling to the continuum — the mathematical reason
   "generic extra harmonics die," consistent with every fabric observation
   (v64 3ω decay, G2-1 sync collapse, QRK-2 pattern collapse).

### I.2 What the fabric already shows
Stable: one clock per flavor sector (flavored branch holds Δω to 98.8%);
ring winding on top of clocks (L_z = nQ). Resonant/decaying: everything
else — including the QRK-1 lines (0.018/0.036 family; flavored splittings
0.108/0.126), which are *widths waiting to be measured*.

### I.3 The correspondence claim (what "commands reason for real physics")
Real physics obeys the same rule: **stable particles are the lowest states
of their conserved-quantum-number sectors** (electron = lightest charged
lepton; proton = lightest baryon); all excited hadrons are resonances with
widths set by coupling × phase space. The fabric *derives* this rule from
survivorship + gap arithmetic rather than postulating it. Deliverable: a
written derivation + numerical confirmation = the program's cleanest export
of "reason for real physics." What most approximates real physics is then:
3 flavors + winding ≈ the fabric's quantum-number lattice, resonance
spectroscopy as its hadron physics.

### I.4 Numerical program (census rungs)
| Rung | Cost | Test |
|---|---|---|
| HC-1 | analysis (~1 day) | Ball-background linear-response (BdG) spectrum — the mode catalog the whole census scores against; also fixes X1's "alignment" definition |
| HC-2 | analysis (free) | Resonance audit: enumerate n·ω_i ± m·ω_clock vs gap for the HC-1 catalog → predicted stable/decaying classification |
| HC-3 | shooter (free) | Flavored-branch Hessian scan over partitions (Q_0,Q_1,Q_2): the stability region = how many clock partitions are truly stable |
| HC-4 | N=64 runs | Line-width spectroscopy: QRK-1 protocol vs kick amplitude; width ∝ amplitude² = golden-rule confirmation; long-T for narrow lines |
| HC-5 | N=64 runs | Overload test: seed 2-frequency content in ONE component (off-lattice of conserved charges) — measure the predicted decay to the sector minimum |

## Part II — Exchange physics: is movement simulated correctly?

### II.1 What is already verified
Translation IS fabric exchange in this kernel, and it is measured: the
moving ball's process ledger closes at **99.75%** (v73 E2 — uptake at the
front, layment at the back), with the 0.25% residue predicted ARTIFACT(dx)
(X2 refinement still pending). Boosted objects carry the correct SR
kinematics (E²=p²+E₀², γω clocks, k=γωv — v70). So single-object motion is
not "fundamentally incorrect": it is process-transfer, verified to the
stated precision, with a known 1–5% lattice group-velocity anomaly as the
honest incorrectness scale.

### II.2 What has never been tested (the suspicion is legitimate here)
1. **Composite transport.** No bound system has ever been moved. Does a
   nucleus+cloud atom drag its cloud (adiabatic co-motion), strip it, or
   shed it? Does the cloud's binding survive acceleration?
2. **Clock transport of bound modes.** A boosted ball's clock dilates (γω,
   measured); do its *bound-mode* frequencies (shell levels) transform
   consistently? If not, moving atoms detune their own shells — a real
   candidate for "movement is simulated wrong."
3. **The mode census of exchange.** WHAT carries the flow — above-gap
   matter waves vs gauge radiation vs near-field convection — has never
   been decomposed. Δω bookkeeping of a moving/absorbing object is
   unmeasured.
4. **Frame dependence.** Lab-frame vs boosted-frame physics differ at the
   lattice anomaly level; never quantified for bound systems.

### II.3 Numerical program (exchange rungs)
| Rung | Cost | Test |
|---|---|---|
| EX-1 ATOM-BOOST | GPU (~20h) | Boost the X10c end-state (ball+cloud) at v=0.02–0.05: co-motion vs stripping; cloud binding & radius in flight; the first moving atom |
| EX-2 sponge-flux spectroscopy | analysis + reuse | Decompose outflow at a shell just inside the sponge into (matter/gauge) × frequency bins for existing runs (x10c, x11pb, surv) — the transfer-mode census |
| EX-3 boost-invariance audit | CPU/GPU refinement set | Same two-ball encounter in lab vs boosted frame at 2–3 dx: quantifies "simulated correctly" as a number vs dx |
| EX-4 clock-transport check | with EX-1 | γ-scaling of the bound cloud's mode frequencies in the boosted atom |
| EX-5 (design only) | paper | Exchange bookkeeping options: continuous-classical (current, correct within lattice artifacts) vs carrier-count ledgers (analysis-level quanta) vs kernel commit rules (rejected: breaks locality/conservation unless redesigned) — pick the reporting standard for atomic-number > 1 work |

### II.4 Relation to atoms and information flow
Atomic stability under motion (EX-1/4), multi-electron structure (D13 —
whatever "occupancy" means, it must survive transport), and information
flow (EMF = gauge waves at c with measured lattice dispersion; matter
signals gap-limited) all reduce to the same audited exchange layer. The
program's rule going forward: **no Stage-4/5 claim without its transport
test** — a parked atom that cannot move is not an atom.

## Order of work
HC-1 → HC-2/HC-3 (all analysis, no spend) → EX-2 (mines existing data) →
HC-4/HC-5 (cheap CPU) → EX-1/EX-4 (the GPU centerpiece) → EX-3 → then the
Stage-4 carbon program inherits whatever the exchange audit says.

---

# AMENDMENT 1 (2026-07-26, post-council — integrates GROUNDING v2 and FOUNDATIONS)

**Part 0 (new, runs FIRST): the ε/foundation battery** — N1–N10 as
specified in GROUNDING v2 §5. N1–N6, N8–N10 are shooter/analysis (zero
spend); N7 is the one GPU-light rung and doubles as the empirical D5′
decision. Part 0 output feeds a TOE Step-3/§D rewrite to the frozen
foundation (GROUNDING v2 §0).

**Census amendments:** HC-1 also delivers n(H_ω) (needed by corrected GSS);
HC-2 classifies by multipole BEFORE arithmetic (massless-channel modes are
golden-rule, not arithmetic); HC-3 uses the corrected negative-index
criterion; HC-4 adds width-floor calibration + a box pair moving the IR
cutoff across a line (ε² scaling without cutoff sensitivity = NULL);
**HC-6 (new, Kimi):** seed deliberately GSS-unstable partitions and verify
decay to sector minima — the converse census.

**Exchange amendments:** EX-1 threshold corrected to v ≪ α_f ≈ 0.053;
ramped protocol above v=0.02; sudden-kick at 0.02 kept as a stripping
measurement; EX-4 retargeted to ball-clock transport.

**Order of work (revised):** Part 0 (N1–N4, N10) → HC-1/HC-2 → N5/N6/N8 →
HC-3/HC-6 → HC-4 (hardened) → EX-2 → N7 → EX-1/EX-4 (ramped) → EX-3.

---

# AMENDMENT 2 (2026-07-26, post N-battery review)
- **Canonical rung protocols:** `council/grok45/N_BATTERY_REVIEW.md` §1
  restores the full N1–N10 designs (formulas, discriminators, kill
  criteria) that GROUNDING v2 §5 compressed — that file governs execution.
- **Order fixes:** N7 runs BEFORE EX-1 (inertia lock precedes any boost
  interpretation); N9 restored to the sequence; a **D7-lite** refinement
  certification (one dt/dx/L pair on a representative observable) rides
  along with Part 0.
- **If halved:** cut EX-3, HC-5 (HC-6 covers the point), EX-1's v=0.05 arm,
  full EX-4 cloud transport; defer N8. Keep N1–N3, N7, N10, HC-1/3/6, EX-2.
- **Part I.1 superseded:** the harmonic-count framing stands, but all
  energy-ontology phrasing in the original body defers to GROUNDING v2 §0
  (E primary, M=E/c², Σ virial excess, ħ_eff/Q as measured ratio).
- Strategy review verdict: current plan reweighted, defended against
  pure-carbon-now and second-sector-first (see N_BATTERY_REVIEW §3).
