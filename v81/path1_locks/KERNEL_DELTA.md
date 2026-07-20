# KERNEL_DELTA — propose `scp_locks` (do not apply until accepted)

**Status:** **ACCEPTED and IMPLEMENTED** (2026-07-19) — see `KERNEL_PORT.md` for land summary and K-L* results.  
**Constraint:** Cosserat nuclear path (v74 Q-balls, complex gauged matter) remains untouched when `n_locks=0`.

---

## Goal

Add an optional **lock carrier** sector on the existing free U(1) gauge medium already present when `complex_gauge=1`, so Stage-3 light charge is a **typed lock**, not a multiplet hump / multi-fab L.

```text
when n_locks > 0:
  free medium  = existing gauge (A, E)  [already in kernel]
  locks        = new array of structs   [new]
  step add-on  = deposit J_lock → Ampère; gather E,B → Boris
```

---

## New module (proposed files)

| File | Role |
|------|------|
| `sfa/sim/scp_locks.h` | Lock struct, API: `locks_init`, `locks_deposit_rho_J`, `locks_push_boris`, `locks_energy`, `locks_write_tsv` |
| `sfa/sim/scp_locks.c` (CPU) / optional `.cu` device helpers later | CIC/Esirkepov deposit; Boris; diagnostics |
| Config keys in existing cfg parser | `n_locks`, `lock_q_i`, `lock_m_i`, `lock_x_i`, … or `locks_seed=path.tsv` |

**Not proposed:** new multiplet fabric; changes to Vt potential; multi-fab L arrays.

---

## Lock struct (kernel-facing)

```c
typedef struct {
    int    id;
    int    type;      /* 0=light charge; later: bag/nucleus-tag */
    double q;         /* ±1 (or small integer) */
    double m;         /* rest-mass label (M-chart) */
    double E_star;    /* sequestered rest energy (ledger) */
    double x[3];      /* continuous position */
    double u[3];      /* proper velocity γv (c=1 units of sim) */
    int    pinned;
    int    alive;
} ScpLock;
```

MVP: charge locks only. Bag/nucleus coupling = later (optional force from Cosserat energy density, not required for positronium-class).

---

## Step hooks (where to call)

Inside the existing complex-gauge substep, **after** free Maxwell/Ampère matter currents are formed, **before** or **as part of** gauge E update:

```text
existing:  matter Φ currents → J_matter
new:       locks_deposit_J(locks, dt) → J_lock
           J_total = J_matter + J_lock   (same staggered locations as gauge)

existing:  E update with J_total; Faraday on B/A as now

new:       interpolate E,B at lock x (same gather as sandbox)
           Boris push locks
           optional: locks_deposit_rho for Gauss diagnostic
```

**CPU (`scp_sim.c`):** single-threaded or OpenMP over locks (N_locks ≪ N³).  
**GPU (`scp_sim.cu`):** defer until CPU path proven in-kernel; locks on host or small device array.

**Critical:** when `n_locks==0`, code path byte-identical to today (no extra J, no branches that change Φ update).

---

## Config (sketch)

```text
n_locks=2
lock0 = q=+1 m=8 x=40,48,48 u=0,0,0 pinned=0
lock1 = q=-1 m=8 x=56,48,48 u=0,0,0 pinned=0
# or:
locks_file=locks_init.tsv
```

Units: same lattice units as gauged Cosserat (\(c=1\) convention already used in U(1) era).

---

## Gauss / charge

- Prefer **Esirkepov** (order-1 or order-2) matched to the discrete divergence used by the existing gauge Gauss diagnostic (`gauss_max`).  
- Sandbox zigzag is a stand-in; kernel should use the scheme that keeps `gauss_max` at the known floor when locks are pinned and after motion.  
- Tripwire: if `n_locks>0` and `gauss_max` drifts above agreed threshold → abort/flag (same culture as v69+).

---

## Energy ledger

Report each step (diag.tsv columns):

| Column | Meaning |
|--------|---------|
| `E_em` | existing gauge field energy |
| `E_locks_star` | \(\sum E_\star\) |
| `E_locks_kin` | \(\sum m(\gamma-1)\) |
| `E_matter` | existing Cosserat energy |
| `locks_alive` | count |

No fudge sinks. Boundary flux remains whatever absorbing BC already does.

---

## SFA / export

- **Prefer:** lock tracks as side-car TSV (`locks_track.tsv`) — no `sfa.h` change.  
- **Optional later:** particle-track chunk in SFA only if volview needs it; requires separate sfa.h authorization use.  
- Grid export of \(\rho_{\mathrm{lock}}\), \(J_{\mathrm{lock}}\) as optional diagnostic columns only if needed.

---

## What must NOT change

- `Vt`, \(\eta\)-coupling, Q-ball init paths, flavored profiles, carbon/nuclear generators.  
- Default configs with `n_locks=0`.  
- Multi-fab product campaigns (out of scope).  
- Pairwise Coulomb helpers.

---

## Validation after port (in-kernel)

Re-run sandbox gates with kernel binary:

1. **K-L0:** one pinned lock, `gauss_max` floor, exterior \(E\).  
2. **K-L1:** rest pair \(D\) scan — monotone \(a_{\mathrm{rel}}\) from medium.  
3. **K-L2:** long soft \(T\gtrsim2000\) — locks alive, \(E_{\mathrm{em}}\) not null.  
4. Regression: `n_locks=0` Q-ball / Gauss tripwire suite (existing v70/v71 smoke).

---

## Implementation order (when accepted)

1. Land `scp_locks.h` + CPU deposit/push with unit tests against sandbox numbers (2D reduced or 3D slab).  
2. Hook config + diag only (`n_locks=0` default).  
3. Enable J_lock in gauge update behind flag.  
4. K-L0…K-L2.  
5. Only then consider GPU and capacity/\(\ell\) (P2).

---

## Acceptance line

> Sandbox P1 is GREEN. Kernel port is **authorized by OP** but **gated on acceptance of this delta** so the nuclear stack cannot regress. This document is the contract for that port.
