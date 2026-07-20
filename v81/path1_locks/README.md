# v81 path1_locks — Locks on gauge medium (PIC-monist)

**Path:** P1 (primary Stage-3 redesign candidate)  
**Thesis:** `REAL = FreeSubstrate + Locks` (v80 shapes **4 + 10 + 6**)  
**Scope:** standalone CPU sandbox only until L0–L2 green.

## State

| Piece | Representation |
|-------|----------------|
| Free medium | 2D TE Yee lattice: \(E_x, E_y, B_z\) with fixed free \(c=1\) (CFL = engine law) |
| Locks | Array of structs `{id, type, q∈{±1}, m, E_star, x, p/u, pinned, alive}` — **not** fields |
| Charge | Integer \(q\) on locks; CIC \(\rho\); zigzag charge-conserving \(J\) deposit |
| Particle definition | Lock object (cannot evaporate by multiplet radiation) |

**Not present:** multiplet \(\phi\), Cosserat light L, multi-fab, pairwise \(1/r^2\), GRIN \(c(x)\).

## Step

1. Yee Faraday: update \(B_z\) from \(E\).  
2. Gather \((E,B)\) at lock → **Boris** push (relativistic).  
3. Advance lock positions; deposit \(J\) with **zigzag / VB-class** charge-conserving scheme.  
4. Yee Ampère: update \(E\) from \(B,J\).  
5. Diagnostics: Gauss residual, energy ledger, lock tracks.

Electrostatic ICs use periodic **SOR Poisson** so L0/L1 start on the Gauss constraint surface.

## Build

```bash
make -C v81/path1_locks
# → bin/locks_pic
```

Requires: `gcc`, OpenMP, `libm`.

## Run

```bash
./v81/path1_locks/bin/locks_pic all results   # L0–L3
./v81/path1_locks/bin/locks_pic l0 results
./v81/path1_locks/bin/locks_pic l1 results
./v81/path1_locks/bin/locks_pic l2 results
./v81/path1_locks/bin/locks_pic l3 results
# or: make -C v81/path1_locks run
```

Outputs under `results/`: `l0_*.tsv`, `l1_force.tsv`, `l2_series.tsv`, `l3_tracks.tsv`.

## Experiments

| ID | Setup | Pass |
|----|--------|------|
| **L0** | Single pinned lock | Gauss floor after Poisson; exterior \(\|E\|(r)\) decreases; short dyn hold |
| **L1** | Rest pair \(q=\pm1\), \(D\in\{8,12,16,20,24\}\) | Monotone \(a_{\mathrm{rel}}(D)\) from **medium only** (attraction sign) |
| **L2** | Soft free pair, \(T\gtrsim 2000\) | Locks alive; \(E_{\mathrm{em}}\) does not null; \(E_\star\) exact |
| **L3** | Tangential \(v_t\) (optional) | No shred-by-definition; report relative revs |

## Kill gates (OP)

Fail / kill path if: multiplet-only particle; pairwise Coulomb inserted; multi-fab L; GRIN gravity; spectroscopy-first.

## Kernel

Do **not** edit `scp_sim` until FINDINGS says GREEN for port and `KERNEL_DELTA.md` exists.

## Docs

- `FINDINGS.md` — scorecard  
- `STATUS.md` — DONE | BLOCKED | KILLED  
- `KERNEL_DELTA.md` — only if L0–L2 green  
