# v102 kernel step, K2e — pre-registration (singlet core in the formation event)

**Written:** 2026-10-09, before any K2e code was written.
**Authority:** user, 2026-10-09: "Go with singlet core, keep going"
(`../README.md` R12).
**The ratchet governs:** the full battery, off by default, byte-inert.

## 1. Why K2e is needed

- K2d's `lockq_p11=2` (score) prices the singlet component of a (1,1)
  pair at P11_CORE = 1.0 (about 104 J0). The lock's evolution under
  that penalty is the standard way to impose the constraint: the
  excluded component appears only virtually.
- **The formation event is not faithful to R12.** Its window uses the
  jump state's mean energy change (K2a's rule). A jump beside an
  existing hole whose state has any singlet part on the new (1,1) link
  has its mean pushed out of the window, so it is refused whole.
- Under R12 only the singlet part is excluded. The other channels are
  allowed and carry no term, so their part of the jump should fire
  (`K2d_DECIDE.md`, limits).

## 2. What K2e builds

- **`lockq_proj`** (default 0). With `lockq_proj=1`, `lq_map` projects
  each jump state onto the allowed subspace (no singlet pair on any
  (1,1) link) before the window and the rate are evaluated.
  - The projector is the cyclic product of (1 − Π_l) over the links,
    repeated to convergence (von Neumann's alternating projections).
    Π_l is the singlet projector of link l on its (1,1) block,
    (1/|P|²)·A†A.
  - The convergence tolerance is 1e-15 of the norm, with at most 200
    sweeps. The largest number of sweeps used is reported.
  - The rate uses the projected weight ‖Π_allowed J ψ‖². The window
    and the ledger use the projected state's energy.
  - It applies to seeded jumps too (they go through `lq_fire`).
  - It is refused unless `lockq_p11` is 1 or 2 (a constraint to
    project onto).
- **A second seeded jump,** `lockq_force2` / `lockq_fmode2` (default 0:
  none). It is applied after the first, so the A/B can start from two
  adjacent holes.
- `cas/k2d_replay.py` gains the same projection (`proj` argument).
- The A/B reference projects independently, by the filter
  exp(−τ Σ_l Π_l) with τ = 200.

## 3. Gates

| | gate | bar |
|---|---|---|
| **E1** | **Inert at default.** The full battery, both kernels. Kernel stdout for `exp=bath`, `blob`, `slit`, `pulse` (T = 20) against the K2d binary (HEAD). K2d's D3 configurations `M2_b_p2_s1` and `M2_b_p0_s1`, full stdout, against the K2d binary | 87 ALL GREEN, bar table identical; every stdout byte-identical |
| **E2** | **A/B of the projected jump.** `lockq_p11=2 lockq_proj=1 lockq_rate=0`, 2000 steps (tL = 40), 3+2 n = 8 and W n = 6, in two cases: (a) one seeded jump at bulk site 1 (beside the end cell); (b) seeded jumps at sites 2 and 3 with complementary modes (two adjacent holes). Against the reference with the exp(−τH_c) projector | the projected weight, ⟨H⟩/J0, and the nh, S^z (I^z) profiles agree to 1e-6 |
| **E3** | **Singlet core in light.** K2d's D3 configurations at `lockq_p11=2` with `lockq_proj=1`: 3+2 n = 10, W n = 7, bright and dim, seeds 1–3 | per run: \|drift_rel\| ≤ 1e-12, `check_violations` 0, `window_viol` 0; the replay (with projection) matches the kernel's final lock line to 1e-8 |

**E3 readouts** are D3's (events, events beside a hole or an end, final
N_h and positions, adjacency, mean B, bound weight), plus the final
state's singlet weight on (1,1) links. That weight is reported, not
gated: under the penalty it is only virtual, of order (J0/104 J0)².

## 4. Endings

| ending | condition |
|---|---|
| **Go** | E1–E3's bars pass |
| **Stop** | any bar fails; the cause is found and stated first |

## 5. Then

K3, under R12, with its own pre-registration: a lock in light over long
times. Nucleons form, walk, and pair at formation, and the run reports
what the lock builds.
