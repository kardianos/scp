# E1b + E2 parallel options evaluation

**Date:** 2026-07-20  
**Skipped:** C4 (accept liquid-drop only)  
**Runtime note:** First attempt thrashing (8 concurrent sims); final data from sequential campaigns, N=48 (E1b scans T=280; B1 N=64 L=24 T=400; E2 N=48 T=40–100).

---

## E1b — free hybrid orbit (B1–B5)

**Gate:** band PASS = revs≥0.75 and late r∈[0.55,1.7] r₀.  
**Result: 0 / 31 band PASS.**

| Option | What | Outcome |
|--------|------|---------|
| **B1** larger box | N=64 L=24, 4 seeds | Best stability: `B1_r10_vf1p0` revs=0.34, r∈[10.0,10.5] — tight but **sub-rev**. No multi-rev. |
| **B2** seed grid | r∈{8,10,12} × vf∈{0.75…1.25} | Best revs: `B2_r8_vf0p75` **0.58** but **inspiral** r→3.75. Over-speed flings to wall. |
| **B3** mass ratio | m∈{1…16} at r=10 | Heavier → fewer revs, flatter r (kinematic). No band. |
| **B4** anti-lock bag core | bag_k∈{0.05,0.15,0.30} | **No effect** vs baseline (identical revs/r to B2 r10 vf1). Bag not load-bearing for hybrid light lock. |
| **B5** soft core | off / 1.5 / 2.5 / 3.5 | **No effect** on metrics (same as soft off). Soft never engaged or irrelevant at these r. |

### Reading

- Force law is correct (E1a); **dynamics still fail multi-rev**.  
- Under-speed → inspiral/core; over-speed → expand; circular seed from continuum F is only ~0.3 rev even when r is flattest.  
- **B1 helps slightly** (less boundary fling; flattest r) but does not promote to band.  
- **B4/B5 do not unlock** orbit under these settings.

**E1b verdict: FAIL** across all options tried. Next would need new physics (capacity co-field on matter density, radiation damping, or different light mass/g), not more seed polish alone.

---

## E2 — multi-center (C1–C3)

**Crude park** flag: two centers, |Δsep| small — **over-counts** short-T stalls at large D.

| Option | Outcome |
|--------|---------|
| **C1** D×phase map | **co** D≤8: seed already 1-cluster or merge; D=10 merge; D≥12 **stall** (sep flat). **anti** D≤10 **separate** hard; D≥12 slow separate. **mix** intermediate. |
| **C2** flavored interlock (ω=1.38/1.38/1.42) | **im1** (2-anti): separate. **im2** (1-anti): D=10 merge; D=12 mild shrink (rdsep~−2.5%); D=14 slight expand. **co** D≤12 merge/1-clust; D=14 stall. |
| **C3** cold large-D | Longer T confirms co D=14–16 nearly flat (Coulomb weak, not bound well); anti still drifts out. |

**Strict park** (|Δsep|/sep₀ < 8%, two centers, sep>4): 15 cells — almost all **D≥12** where force is weak over T=40–100. **Not** a demonstrated force-zero standoff.

### Reading

- Pattern matches v71: co merges, anti repels; **no repel-in/attract-out bound window** that parks.  
- Large-D “park” = **kinematic stall**, not binding.  
- Flavored interlock profile does **not** produce a static multi-center molecule in this smoke (im2 D=12 is the least bad, not a PASS).

**E2 verdict: NULL/FAIL** for static multi-center under C1–C3. Phase + gauge alone insufficient at g=0.05 for retained A.

---

## Combined evaluation

| Track | Result | Implications |
|-------|--------|--------------|
| E1b B1–B5 | **FAIL** multi-rev | Hybrid Coulomb works; free light lock does not park. Seed/box/mass/soft/bag-core insufficient. |
| E2 C1–C3 | **NULL** static multi-center | Merge (small D) or separate/stall (large D). No 2B park. |

**Do not** claim Stage-4 free-orbit pilot or Stage-2B multi-center from this campaign.

### Highest-leverage next (if continue)

1. **E1b:** matter-density capacity wall (true co-field), or 3D radiation-aware long run of **only** `B1_r10_vf1p0` with T≫T_orb after fixing box/sponge interaction.  
2. **E2:** only if 2B stays priority — need **new binding** (not more phase scans): e.g. second gauge, saturating attraction, or lock-unit E3.  
3. **Deprioritize** B4 pairwise bag on free hybrid light and further soft-core knobs.

### Data

- `v83/e1b/results/summary_e1b.tsv`, `master.log`  
- `v83/e2b/results/summary_e2.tsv`, `master.log`  
- runners: `v83/e1b/run_campaign.py`, `v83/e2b/run_campaign.py`
