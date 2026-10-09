# v102 kernel step, K2a — decision (the `lockq` module in the C kernel)

**Written:** 2026-10-08, after the runs in `cas/out/k2a_*`. Gates:
`K2a_PREREG.md`. Code: `v97/kernel/lockq.c` (new), and hooks in
`v97/kernel/freecell.c` (each one guarded by `P.lockq`).

## Ending: **Go (K2b)**

The lock runs inside the classical kernel and matches the Python lock
model to every printed digit. Light makes holes in it, one atom per
event, and the energy ledger closes at the kernel's own floor.

- **One registered gate failed as first written.** G5 failed on W: the
  window was set on the E_f part only. The rule was corrected and G5
  now passes (§G5).
- **The demonstration needs K2b.** The kernel's light is one complex
  mode, so every hole it makes has the same charges, and a lock full of
  same-charge holes cannot form the all-singlet pair (§5).

## G1. Inert at default: **PASS**

- **Kernel stdout before and after,** with the module off and the V3a
  law: `exp=bath` (T = 20, 48 lines), `exp=blob` (L = 16, T = 20, 53
  lines) and `exp=slit` (the double-slit apparatus, T = 20, 210 lines)
  are all **byte-identical**.
- **The full battery, both kernels** (`cas/out/k2a_battery_base.out`,
  `k2a_battery_new.out`): **ALL GREEN (87 bars)** before and after, and
  the bar table is identical line by line.
  - 70 of the 71 run logs are byte-identical.
  - The 71st, `p2lc`, is a Go-kernel run (the Go kernel is untouched).
    It differs only in its wall-clock timings (`ms_sync`, `ms_async`).
- There are no new compiler warnings (9 before, 9 after, all
  pre-existing).

## G2. Hamiltonian A/B: **PASS**

The lock ground state in the zero-charge, no-hole sector, against
`../perm/cas/lo.py` (`cas/k2a_ab.py`):

| lock | sector dim | C module, E_gs/J0 | lo.py, E_gs/J0 |
|---|---:|---:|---:|
| 3+2, n = 8 | 534 | −18.278681579337 | −18.278681579337 |
| W, n = 6 | 8272 | −20.325778285156 | −20.325778285156 |

They agree to 1e-12 J0 (bar: 1e-9).

## G3. Dynamics A/B: **PASS**

The setup:

- a seeded jump (photon mode 0) at bulk site 3 (3+2) or site 2 (W);
- then 2000 kernel steps, to tL = 40, with formation off
  (`lockq_rate=0`);
- outputs in `cas/out/k2a_g2g3_*`.

The hole density, S^z and I^z profiles agree with Python's exact
exponential on every site to the 10 printed digits (bar: 1e-6). For
example, W n = 6 gives nh = 0, 0.0020765695, 0.9938304410,
0.0040899023, 0.0000030872, 0 in both.

Supplementary checks (not registered):

| run | profiles | ⟨H⟩ |
|---|---|---|
| 20 000 steps (tL = 400, 3+2); the hole has spread to 6% off its site | all agree to 10 digits | constant to 1e-13 |
| `lockq_every=10` (dt_L = 0.2) | identical to `lockq_every=1` | |

The Krylov error bound stays below 2e-17 per step. The norm drifts by
6e-13 over 20 000 steps; the ledger uses the normalized expectation.

## G4. The ledger: **PASS**

The setup:

- `exp=pulse` (a field packet centred on the lock, kx = 0, V3a law),
  T = 20;
- both fabrics, two brightnesses (amp 1.0/σ 2.2 and amp 4.0/σ 3.0),
  three module seeds each;
- outputs in `cas/out/k2a_g4_*`.

| lock | amp | events (seeds 1/2/3) | final N_h | max \|drift_rel\| over the run |
|---|---:|---|---|---:|
| 3+2, n = 8 | 1.0 | 1 / 1 / 2 | 1–2 | 4.9e-16 |
| 3+2, n = 8 | 4.0 | 6 / 5 / 6 | 5–6 (of 6 bulk cells) | 3.2e-15 |
| W, n = 6 | 1.0 | 1 / 1 / 1 | 1 | 4.9e-16 |
| W, n = 6 | 4.0 | 4 / 2 / 4 | 2–4 (of 4) | 3.5e-15 |

The bar is 1e-12. The same pulse with no lock drifts by 4.9e-16
(amp 1.0) and 3.6e-15 (amp 4.0), so the module, with its events, stays
at the classical kernel's own floor.

## G5. Event bookkeeping: **FAIL as registered, PASS after a rule correction**

**What was registered.** The band window was set on the E_f part of the
energy change. The exchange shift was claimed to be "at most a few J0
(about 0.03)".

**What happened** (`cas/out/k2a_registered_*`; registered rule,
amp 4.0):

- On W, 2 of 3 seeds had one event whose exact space quantum fell
  outside the band, at |δ| = 0.135 > 2 s_k = 0.12.
- The cause: a hole formed next to an existing hole loses that link's
  (2,1) exchange, about **10 J0 = 0.095** for W. That is not "a few J0".
- The ledger stayed exact, but such an event should not fire: the
  space band has no state at that energy.

**What was also found.** The registered draw (at most one event per
lock step) saturates under bright light: the expected number of events
per step reached 0.44, above the registered warning level of 0.1.

**The correction** (made before the K2a commit; both changes are in
`lockq.c`):

1. **The exact window and rate.**
   - For each bulk site i with a lit neighbour, the module maps J_i ψ
     into the next sector, which is the same for every site, because
     photon mode 0 fixes the charges. It is built once per event.
   - It then takes that state's exact mean energy change ΔE_i.
   - The window and the golden-rule density of states use
     δ = w1e_o − ΔE_i − Ω, the same energy the space quantum carries in
     the ledger.
2. **A Poisson draw** inside each lock step: waiting times at the total
   rate, with the rates recomputed after each event. This replaces the
   single competing-risk draw.

**Re-run** (the same 12 runs as G4):

- 34 events, **0 outside the band**. The largest |δ_exact| is 0.1186
  (W), inside 0.12.
- Every event removed exactly one atom (ε_o = A0·w1e_o/2π = 0.30200 at
  the vacuum pitch) and raised N_h by one.
- Every event moved the lock to the sector shifted by the photon mode's
  charges. A jump that left that sector would have stopped the run
  (`FATAL`); none did.

**The rule is now the rule.** The physics claim is unchanged: K1-1's
jump, at the golden-rule rate onto the space band. The window is now
evaluated where the claim puts it.

## 4. What one formation event looks like (reported)

A 3+2 event at the vacuum pitch (w1e_o = 1.65):

- the atom brings 0.30200;
- the lock gives up 0.21627, its E_f step minus the exchange it loses;
- the space quantum leaving into the neighbour's `Es` carries 0.51826,
  i.e. ħ(Ω + δ) with δ = −0.068.

**Formation is exothermic for the lock.** The hole is made by light,
and the lock's energy goes out into space through the boundary link
(R10).

## 5. What K2a shows the kernel cannot yet do

**Every hole the kernel makes has the same charges.** The kernel field
is one complex mode, so every photon is mode 0, (S^z, I^z) = (+½, +½);
under 3+2, S^z = +½.

- Under bright light the lock fills with same-charge holes: N_h reaches
  every bulk cell in 4 of the 6 bright runs.
- The bound nucleus is the all-singlet pair: W (S, I) = (0, 0), 3+2
  S = 0. It needs holes of opposite charges. One kernel mode cannot
  make it.
- **This is what K2b is for:** four complex modes per cell, the
  isospin and spin partners, empty by default.

## 6. Costs and limits

- **Size.** Exact sectors cap the lock at n ≤ 12 (3+2) and n ≤ 7 (W).
  The largest sector met here was 13 728 states.
- **Speed.** A lit 3+2 run (n = 8, T = 20) takes about 30 s, against
  about 14 s for the registered rule. The added cost is one map and one
  matrix product per lit site per step. Dark steps cost one Krylov step.
- **Fencing.** The module refuses to run with `bath_xb`, `door_g`,
  `bh_r`, `k_cant`, `amp_nat` or `pattern_only` set: those passes are
  not fenced.

## Next

K2b: widen the kernel field to four complex modes (`fa3`…`fa8`, or an
equivalent layout), empty by default. Its gate is the battery at
defaults, byte-identical.

---

## Correction, 2026-10-08 (found in K2b; `K2b_DECIDE.md` §B4)

**A latent heap overflow in the committed K2a code.**

- `lq_rawA` was sized by `NSLOT`, the slot table's high-water mark,
  instead of its capacity `NLMAX`. A new slot touching the lock, above
  the install-time high-water mark, was written past the end.
- It crashed one K2b run. In K2a it could have corrupted memory
  silently.
- **K2a's recorded results stand.** All 12 G4 configurations, rerun
  with the fixed binary, give event lines and results identical to
  `cas/out/k2a_g4_*.out`.
- G1 is unaffected: the fence runs only when `lockq=1`.
