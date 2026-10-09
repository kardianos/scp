# v102 kernel step, K2a — pre-registration (the `lockq` module in the C kernel)

**Written:** 2026-10-08, before any K2a code ran. Plan: `PLAN.md`. K1:
`K1_DECIDE.md` (Go).
**Authority:** user, 2026-10-08: "Proceed with K1 and K2a, commiting
before, and between each step."

## 0. A correction to K1-4, found by reading the kernel before any K2 code

The kernel's `(fa1, fa2)` are not two transverse planes. They are the
real and imaginary parts of **one complex amplitude** per cell:

- the on-site step rotates (fa1, fa2) by `w1e·dt`;
- the hop is exp(−iτσ_x) between cells.

So kernel light carries no polarization and no isospin. **Until K2b
widens the field, every kernel photon is one fixed mode ν = 0**, which
is (S^z, I^z) = (+½, +½) under W and S^z = +½ under 3+2. Every
formation event therefore makes the same hole type. This replaces
K1-4's reading of "(fa1, fa2) as a Jones vector".

## 1. What K2a builds

A module `lockq` (file `v97/kernel/lockq.c`, included by
`freecell.c`), **off by default** (`lockq=0`).

**Knobs.**

| knob | meaning |
|---|---|
| `lockq` | 0 off, 1 on |
| `lockq_M` | 2 or 4: 3+2 or W |
| `lockq_n` | lock cells, ends included |
| `lockq_t` | the T-floor hop, in J0 units (default 0.1) |
| `lockq_g` | g/Δ (default 0.1) |
| `lockq_rate` | scale on the formation rate (default 1) |
| `lockq_every` | kernel steps per lock step (default 1) |
| `lockq_seed` | the module's own RNG seed |

**The lock.**

- A chain of `lockq_n` cells chosen at init: the cell nearest the box
  centre, then greedily the slot neighbour farthest along +x not yet
  used.
- The cells are pinned. Their classical `Ee`, `fa` and `Em` are zeroed,
  and their field energy is carried by the module.
- Every slot touching a lock cell has its lens area `sA` set to 0 each
  step. So no classical field, space or dense transport crosses a lock
  cell's links. That is R9 at the boundary, and R10 for everything
  except the formation event.
- Lock cells are skipped in the per-cell conversion pass (D1).

**Its state.**

- An exact vector over the leading-order lock model of
  `../perm/cas/lo.py`: R6 end cells (N = 1) and mobile bulk cells
  (N = 2 or 1), in a sector of (2S^z, [2I^z], N_hole).
- The Hamiltonian is the −(g²/Δ_type)A†A exchange plus the T-floor hop.
- Ground state by Lanczos at init (N_hole = 0, charges 0).
- Each lock step evolves by a Krylov exponential with dt_L =
  `lockq_every`·dt.

**The formation event** (R10, K1-1, §0):

- **Candidates.** Each pair (i, o): a bulk lock cell i, and a non-lock
  cell o that shares a live slot with i whose lens area, before the
  fence, is > 0 (a boundary link). Cell o holds n_o = floor(Ee_o/ε_o)
  whole atoms, with ε_o = A0·w1e_o/2π. A pair with n_o = 0 is not a
  candidate.
- **The rate** for each candidate is the golden rule onto the 1-D space
  band. The space band is the same lead K1 integrated out, at
  t_s = s_k:
  - r = Γ(δ)·‖J_i ψ‖²·n_o, with Γ(δ) = g²·√(4 s_k² − δ²)/s_k² for
    |δ| < 2 s_k and 0 outside the band;
  - δ = w1e_o + [E_f(2) − E_f(1)] − Ω is the space quantum's detuning
    from the link pitch.
  - The window is set on the E_f part. The exchange shift is at most a
    few J0 (about 0.03, against a band half-width of 0.12). Each event's
    exact shift is reported (G5).
- **The draw.** At most one event per lock step, as a competing-risk
  draw:
  - fire with probability P = Σ r·dt_L;
  - choose the pair with weight r.
  - If P > 0.1 the module warns that the lock step is too coarse.
- **The action:**
  - |ψ⟩ ← J|ψ⟩/‖·‖ with J = Σ_μ P[μ, 0] b_iμ, into the new sector;
  - Ee_o loses one atom (fa_o scaled);
  - Es_o gains ε_o − ħ·ΔE_lock, the exact change of the lock energy
    after the jump (rate → kernel energy by ħ = A0/2π);
  - the ledger closes by construction.
- **This departs from K1-4.** K1-4 read "the event's rate is the K1-1
  probability per resonant atom". An atom in a kernel cell sits still;
  it is not a scattering photon with a flux, so a probability per photon
  has no clock. The golden rule supplies the clock from (g, s_k) alone,
  with no fitted number. K1-1's fidelity (1.0000) is unaffected: it is a
  property of the jump state, not of the rate.

**Installation** (at t = 0, before the initial ledger `E0`):

- The lock cells' classical contents are removed (a seed). Any in-flight
  dense energy on a fenced slot is moved to its non-lock end's `Em`.
- The lock's energy is then ħ·(Σ_i E_f(N_i) + E_gs).
- Both numbers are printed.

**The ledger.** `total_energy()` adds ħ·(Σ_i E_f(N_i) + ⟨H⟩) when the
module is on. ⟨H⟩ is conserved by the exact evolution between events,
to the Krylov tolerance (≤ 1e-13 per step).

**Placement.** The chain starts at the cell nearest
(box centre + `lockq_dx` along x), with `lockq_dx` = 0 by default.

## 2. Gates (all must pass before the K2a commit)

| | gate | bar |
|---|---|---|
| **G1** | **Inert at default.** The full battery (`./battery`, both kernels) at defaults | every bar identical to the pre-change baseline. Kernel stdout for `exp=bath`, `exp=slit` and `exp=blob` (T = 20) byte-identical before/after |
| **G2** | **Hamiltonian A/B.** `lockq=1`, no light in reach, 3+2 n = 8 and W n = 6: the module's lock ground energy against `lo.py` for the same chain | to 1e-9 J0 |
| **G3** | **Dynamics A/B.** A forced event at a fixed bulk cell at step 0, then 2000 steps | the hole-density profile against Python (`lo.py` + exact exponential), to 1e-6 |
| **G4** | **The ledger.** `lockq=1` with the kernel's light reaching the lock (`exp=bath`) | the total-energy drift stays within the kernel's own floor (|drift_rel| ≤ 1e-12 over T = 20) while events fire. The number of events and the hole count are reported |
| **G5** | **Event bookkeeping.** For each event: one atom leaves Ee_o, N_hole rises by 1, the charges change by the fixed mode's, and the resonance window held | all events, by an in-kernel assertion count |

## 3. What K2a does not do

- It does not widen the kernel field (K2b).
- It does not run MPS (K2c).
- It does not add the capture rule (K1-3: about 10¹⁰ photons per pair).
- It does not touch `v97/fab`.
- It does not change any law constant.
