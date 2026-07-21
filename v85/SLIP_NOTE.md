# Side note — SLIP test (commensurate-harmonic exchange) + boundary-backwards
**Date:** 2026-07-21 · **Status:** tests run (SLIP-0 executed); verdict below.
Context: TOE_v85.3 §D; proposal was two interleaving harmonics on the shared
compact fiber, with phase-slip ("planetary conjunction") quantized exchange.

## 1. SLIP-0 — executed (scratchpad `slip0.py`)

Two weakly coupled Q-balls form a bosonic Josephson junction: charge
imbalance q = Q₁−Q₂ conjugate to relative phase δ; dq/dt = −2I_c sin δ,
dδ/dt = Δω₀ + (dω/dQ)·q + b·sin δ, with dω/dQ = −4.2e-4 from the measured
branch. Question: is charge transfer per 2π slip a fixed quantum?

**Result: NO — decisively.**
- Δq/slip ranges **0.26 → 3.6** across (Δω, I_c, b) — factor ~14, exactly
  tracking the analytic form I_c(2π/b)[1 − Δω/√(Δω²−b²)] (linear in I_c,
  divergent at the lock edge: 19.0 per slip at Δω = 1.017b).
- With measured back-reaction, **successive slips grow**: 1.65 → 6.77 over
  15 slips (each slip shifts Δω, changing the next transfer).
- What IS exact: per-slip **action** ∮(dδ/dt)dt = 2π (topological, trivially)
  and per-slip repeatability at fixed parameters (spread ~1e-6): slips are
  precise, countable **events** — with continuous, parameter-dependent
  **amounts**.

**Verdict: NOT VIABLE as an amount-quantization mechanism.** Slips give
jump *timing* (episodic, countable transfer — genuine quantum-jump
phenomenology in time), not jump *quanta*. The v86 full write-up and agent
review are **not** triggered.

## 2. Problems noted

1. Amount-per-slip is a continuous function of overlap, detuning, and pull —
   and drifts slip-to-slip under back-reaction. No universal unit.
2. The lock-edge divergence means near-resonant pairs transfer *large*
   variable amounts — the opposite of gentle quantized exchange.
3. What survives of the proposal: (i) fiber-harmonic **charge-ratio**
   quantization (gauge2 — structural, untouched by this failure);
   (ii) slip-event discreteness as the *timing* half of a jump mechanism.
   The *amount* half still requires carriers (closures) — CARR-1 remains
   pivotal, exactly as FEEDBACK ranked it.

## 3. Boundary-backwards: start at the quantum/classical wall, walk into the theory

**3.1 The wall is a theorem, not a difficulty.** A classical field theory is
a commutative observable algebra ⇒ it is a local-hidden-variable model ⇒
Bell-bounded. D6 (Bell out of scope) is *permanent* for this kernel — no
mechanism, bath, or architecture crosses it. Working backwards, this is the
one boundary that cannot move.

**3.2 Where the fabric sits relative to the wall.** In quantized field
language, a classical field configuration is a **coherent state**:
phase-definite, number-indefinite. Nature's particles are **number states**:
number-definite, phase-indefinite. The boundary variable is number–phase
conjugacy (ΔN·Δφ ≥ ½) — and the classical theory has *no floor*. Every
difficulty we have logged is this one fact wearing costumes: continuous Q
(number-indefiniteness), continuous transfer (no number jumps), 29% chirp
(phase-definite transitions), fractional interruptions.

**3.3 What CAN make integers on the classical side: topology, and only
topology.** Walking back from the wall, the classical integer-makers are:
- **Spatial winding** (vortices/flux): integers when the field *cannot reach
  zero* — stiff modulus. The v73 ring failed flux quantization (0.374 vs
  2π/g = 126) at LOW Q where the modulus is soft. **Prediction: flux/
  circulation quantization quality improves with interior stiffness κs**
  (×12 knee suppression near capacity). → **TOPO-1** (config): re-run the
  ring test on a near-capacity, thin-wall ring; measure trapped flux vs Q.
- **Fiber harmonics** (compact gauge): integer charge *ratios* by Fourier
  construction — the gauge2 sector, dormant in the kernel today. → **G2-0**.
- **Temporal winding** (slips): integer *events* — measured above; amounts
  continuous.
Dynamical mechanisms (Arnold tongues, erosion floors) yield plateaus and
floors — never integers.

**3.4 The landing.** All viable quantization in this theory is
**topological**; amount-quantization can only be carried by *closures*
(carriers), whose all-or-nothing integrity is itself topological (H5).
The reachable program: topological integers (winding, fiber ratios, slip
events) + selection floors (Q_thr) + episodic transfer. The unreachable:
noncommutativity, Bell violation, true number states. The theory's honest
ceiling is "structural analog up to the commutative wall" — now derived
from the boundary itself rather than assumed.

**3.5 Consequence for the shape question.** The Q-degeneracy remedy (second
sector) stands, but boundary-backwards adds a requirement: to quantize
charge *within* a sector, the matter field needs a topologically protected
winding — i.e., a **stiff-modulus / compact-target regime** where |Φ|
cannot pass through zero. That is a potential-shape (kernel) question —
option C below — and it is the single deepest shape candidate left.

## 4. Options going forward

| Option | Level | Action |
|---|---|---|
| A. TOPO-1 | config | Stiff-ring flux quantization vs Q — direct test of §3.3's prediction; uses v73 ring machinery |
| B. G2-0 | config | Baseline the dormant gauge2 fiber (charges 1,1,−2); then commensurability tests |
| C. Compact/stiff matter target | **kernel — needs explicit authorization** | Modulus-stiffened potential (|Φ| bounded away from 0 in bulk): makes charge winding topological. Design-only until A shows the stiffness→quantization gradient is real |
| D. CARR-1 / SURV-1 | config | Unchanged from TOE_v85.3 — carriers and the erosion floor carry the amount-quantization burden |

Slip physics is retained as the **timing** component (jump phenomenology,
line-shape sharpening) and folded into X7's deliverables (count slips,
measure transfer-per-slip on the lattice — SLIP-1 stays as an X7 extension,
now with the ODE prediction to compare against).

## 5. Files
`slip0.py` (scratchpad) — SLIP-0 model + results · `TOE_v85.3.md` §D —
parent diagnosis · v73 ring docs — TOPO-1 baseline · scp_sim.c:67–71,
1137–1141 — gauge2 (read-only cites).
