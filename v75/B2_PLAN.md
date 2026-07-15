# B(2) Plan — Unlock Φ_Q from Φ_C

**Status:** draft after B1 G2/G3 PASS (F11)  
**Prerequisite:** B1 Shape β complete — private bags isolate L; Coulomb works with U^q.

---

## 1. Goal

In B1, \(\Phi_Q \equiv \Phi_C\) (β-copy). Nuclear **bag** and **EM charge** are
the same profile. B2 separates them:

| Fabric | Role | B1 | B2 target |
|--------|------|----|-----------|
| **C** | Nuclear bag / color medium | evolved; bag ON | evolved; bag ON; \(q_C=0\) preferred |
| **Q** | Charge / baryonic matter | **locked copy of C** | **independent** dynamics; \(q_Q=+g\) |
| **L** | Light opposite charge | independent; \(q_L=-g\) | unchanged |

**Science question:** Can a charged Q-lump exist *outside* or *partially
overlapping* the bag, with weak portal \(\varepsilon_{CQ}\), without
reintroducing short-range L merge?

---

## 2. What B1 already proved (do not re-litigate)

- Shared A + opposite fabric charges → Coulomb attract (G2 = E-lite ΔD).
- Private \(s_C\) / \(s_L\) → head-on **no merge** (G3 pass-through).
- Seed rule: same-sign Noether ω when \(q_L=-q_Q\).
- Transport must use \(U^q=\mathrm{e}^{iq\theta}\).

B2 is **not** needed for “does multi-fabric Coulomb work?” — that is closed.

---

## 3. Config / kernel surface

Already sketched in `MULTIFABRIC_SPEC.md`:

```
n_fabrics   = 3
mf_stage    = 2          # was 1 in B1
mf_lock_CQ  = 0          # unlock
eps_CQ      = 0.0…ε      # portal (start 0 for G4)
q_C = 0, q_Q = +1, q_L = -1
```

### Kernel work (when implementing)

1. **Verlet / force:** when `mf_lock_CQ=0`, evolve Q with own bag policy:
   - **Preferred:** bag force from \(s_C\) only still (Q is “colorless charge”
     in C’s potential, or Q has **no** self-bag and only kinetic+gauge+ε).
   - **Alt:** Q has own weak \(V_Q\) if needed for stability (document if used).
2. **Portal \(\varepsilon_{CQ}\):** bilinear or soft constraint
   \(\sim\varepsilon|\Phi_Q-\Phi_C|^2\) or mass-mixing term — **only if G4 needs
   it**. G4 first with \(\varepsilon=0\).
3. **Current:** \(\rho_{\mathrm{em}}=q_C\rho_C+q_Q\rho_Q+q_L\rho_L\) (already);
   with unlock, Q and C densities may separate.
4. **I/O:** 54-col already has room for L; Q may need columns if not locked
   (today Q is not written separately under lock — **check**: if β-copy only
   on device, Q may be absent from SFA). B2 may need **Q columns** or
   diagnostic-only Q dumps.
5. **G0:** `n_fabrics=1` and B1 path must still regress.

**Do not start kernel work until G4 seed plan is frozen** (user may prefer
lighter-L first — see §6).

---

## 4. Gate sequence (B2)

```
G4  mf_lock_CQ=0, eps_CQ=0
    Seed: single heavy lump on C+Q co-located at t=0 (copy C→Q once), L off
    Q: does Q stay bound as a charged ball, or radiate / dissolve?
      ├─ Q STABLE alone near C   → prefer bag on C, charge on Q, small ε later
      ├─ Q UNSTABLE alone        → B2-alt: charge on C, or ε-bind Q into C
      └─ Q drifts off C freely   → measure separation vs Coulomb self-energy

G4b Optional: Q only (C=0) pure charged ball under q_Q — existence window

G5  Ramp eps_CQ: force C+Q composite (nuclear “proton” = bag+charge)
    Metric: binding ΔE, Q stays inside bag radius, no L yet

G6  G3 redo with true C≠Q: heavy = C bag + Q charge, L opposite
    MUST still no-merge at head-on

G7  Hydrogenoid: slow L around heavy composite (orbit or adiabatic)
```

### Pass / fail criteria

| Gate | PASS | FAIL action |
|------|------|-------------|
| G4 | Q lump lifetime ≫ 100, Q conserved, gauss floor | keep charge on C (B2-alt) or add ε |
| G5 | C+Q bound, single center of mass for both | retune ε; do not proceed multi-L |
| G6 | D_min contact, two sectors survive | isolation broken under unlock — debug |
| G7 | bound or long-lived arc, no merge | lighter L / softer v first |

---

## 5. Seeds & campaigns (minimal)

| Run | N | T | Notes |
|-----|---|---|--------|
| G4 smoke | 64–96 | 50 | CPU or cheap GPU |
| G4 science | 192 | 200 | co-located C=Q init, free |
| G4b | 192 | 200 | Q-only seed |
| G5 scan | 192 | 200 | ε ∈ {0, 0.01, 0.05, 0.2} sketch |
| G6 head-on | 192 | 400 | same as G3 kinematics |
| G7 orbit | 192 | 400+ | vt from G2 a_rel |

Reuse `gen_mf_pair_boost` + new `gen_mf_CQ_unlock` (init Q from C profile,
optional offset).

---

## 6. Priority fork (choose one primary)

After F11, two valid nexts:

| Priority | Path | Why |
|----------|------|-----|
| **A. Nuclear honesty** | **B2 G4–G6 first** | Bag ≠ charge is the Option B story for nuclei |
| **B. Atom ASAP** | Lighter L + multi-L on **B1 lock** | G3 already isolates L; hierarchy may matter more than unlock for “carbon atom” |
| **C. Tooling** | Remote renders + Cloud Sync | Cheap; unblocks disk/cost without new physics |

**Active priority (2026-07-13): B — lighter L on B1 lock**, then multi-L /
hydrogenoid. B2 G4 unlock is deferred until hierarchy or carbon packing needs
true bag≠charge. Short **C** (archive/render) runs in parallel.

---

## 7. Out of scope for first B2 slice

- Second gauge / color SU(N)
- Fermionic exclusion
- Nested grids (unless N=192 insufficient for C≠Q separation scale)
- Self-tuning ε (only after G4/G5 manual)

---

## 8. Checklist before coding

- [ ] Freeze Q bag policy (no self-bag vs weak V_Q)
- [ ] Confirm SFA columns for free Q (or diag-only)
- [ ] G0 regression plan
- [ ] G4 seed generator
- [ ] User green-light on Priority A vs B
'''