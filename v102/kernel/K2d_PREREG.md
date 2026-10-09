# v102 kernel step, K2d — pre-registration (the (1,1) pricing in the kernel; a lock in light)

**Written:** 2026-10-09, before any K2d gate ran. The code (below) was
written and compiled, and timed on two configurations that are not
gates.
**Authority:** user, 2026-10-09 ("Try different mechanisms, see what
works … ideally … in a given configuration, it would also force a
specific size"), and `../bind/PREREG.md` §6. K2d puts the pricings
tested in `../bind/` into the kernel's lock, all of them, as options.
**The ratchet governs:** the full battery, off by default, byte-inert.

## 1. What K2d builds

- **`lockq_p11`** (default 0) prices the (1,1) link in `lq_two_build`,
  exactly as `lo.py` `P11` does:

  | value | name | the term on the (1,1) block |
  |---|---|---|
  | 0 | zero | none (today) |
  | 1 | core | +1.0 on every (1,1) configuration |
  | 2 | score | (1.0/\|P\|²)·A†A |
  | 3 | upper | (E₊/\|P\|²)·A†A, the field branch |
  | 4 | lower | (E₋/\|P\|²)·A†A, the space branch |

  E± = Δ₁₁/2 ± √(Δ₁₁²/4 + g²|P|²), computed from the law table at
  install.
- When `lockq_p11` = 0, the lock's code path is unchanged: the original
  coefficient expression is evaluated in the same order. A header line
  `# LOCKQ p11=…` is printed only when `lockq_p11` ≠ 0.
- **The analysis tool** `cas/k2d_replay.py` rebuilds a kernel lock's
  whole history from its event lines in `lo.py`. It checks the kernel's
  final profiles and ⟨H⟩ against it, then reports:
  - the final hole configuration;
  - the adjacency, Σ⟨n_i n_{i+1}⟩;
  - the mean binding B;
  - the weight on bound levels (B_k < −0.05 J0).
- `cas/k2a_ab.py` takes the pricing from the environment (`P11`); its
  default is unchanged.
- Nothing else changes: no law constant, no event rule, no `v97/fab`.

## 2. Gates

| | gate | bar |
|---|---|---|
| **D1** | **Inert at default.** The full battery, both kernels. Kernel stdout for `exp=bath`, `blob`, `slit` and `pulse` (T = 20, the V3a law) against the K2b binary (HEAD). With `lockq=1` and `lockq_p11=0`: K2b's B4 configurations (3+2 n = 8, F = 2, split, seed 1; W n = 6, F = 4, split, seed 1), full stdout, against the same binary | 87 bars ALL GREEN, bar table identical; every stdout byte-identical |
| **D2** | **A/B, every pricing.** `lockq_p11` = 1, 2, 3, 4; a seeded jump at bulk site 1 (next to the end cell, so a (1,1) end link) with photon mode 0; `lockq_rate=0`; 2000 steps (tL = 40); 3+2 n = 8 and W n = 6. Against `cas/k2a_ab.py` with `P11` set | ⟨H⟩/J0 and the nh, S^z (and I^z) profiles agree to 1e-6 |
| **D3** | **A lock in light** (registered readouts, ledger bars). `exp=pulse`, split light (`fseed_mode=-1`), V3a law, T = 20: 3+2 n = 10 (F = 2) and W n = 7 (F = 4); bright (amp 4, σ 3) and dim (amp 1, σ 2.2); module seeds 1, 2, 3; `lockq_p11` = 0 … 4 | per run: \|drift_rel\| ≤ 1e-12, `check_violations` = 0, `window_viol` = 0; the replay agrees with the kernel's last lock line to 1e-8 (⟨H⟩/J0 and every profile entry) |

**D3 readouts** (reported, not gated), per run and summarized per
pricing:

- the number of events, and how many formed next to an existing hole
  or to an end cell;
- the final N_h and hole positions;
- the adjacency;
- the mean B;
- the bound weight.

**The question D3 answers:** in this configuration (a lock in light of
mixed modes), does a pricing leave the lock holding bound pairs, with
no adjacent holes? And is the number of holes in its bound objects
forced (2 = the 3+2 cluster size `../bind/` measured), or not?

## 3. Endings

| ending | condition |
|---|---|
| **Go** | D1, D2 and D3's bars all pass. D3's physics is reported as measured |
| **Stop** | any bar fails. The cause is found and stated before anything else |

## 4. Not in K2d

- No new law constant. `core` (1) is a new rule, if it were ever
  adopted, and that needs the user.
- No new event rule: formation stays K2b's.
- No MPS lock. The exact sectors cap the lock at 3+2 n ≤ 12 and
  W n ≤ 7.
