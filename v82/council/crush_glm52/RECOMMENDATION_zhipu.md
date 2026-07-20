# v82 RECOMMENDATION — Stage-3 after P1 multi-rev miss

**Advisor:** claude · **Date:** 2026-07-19 · **Stance:** stay monist; the missing piece is a **short-range sign change**, not a new gauge charge.

---

## Verdict (one paragraph)

P1 is half-finished, not failed. The attractive channel is real and matches continuum — L1 `F_along(D=8)=1.981e-2` vs `1/(2πD)=1.989e-2`; N=256 anti-lock cuts `Δsep` from **+10.7** (base) to **+7.4** (anti). Every bag tried so far is **monotonic in r toward lock center**, so by construction it lands in exactly one of two regimes: too weak (`null`: revs `0.12`, `Δsep +5.05`, `min sep 11.8`) or strong enough to drive collapse-through (`grid`/`vt`: `min sep 0.6–1.2`, revs `≤0.29`, pass-through). A monotonic co-field cannot produce a stable radius; Bertrand-class orbit closure needs balance, balance needs a sign change. The fabric already supplies one in ontology: depletion of free measure by overlapping lock footprints (v73 process-form; OP shape-3 slot still empty). Build that.

---

## 1. Ranked next paths

| Rank | Path | Why | Citation |
|------|------|-----|----------|
| **1** | **P4-seq**: scalar `n_free(x)=n0−Σ_i w_i(x)` (CIC footprints); core force `F_core=+k_core (1−n_free/n_crit)_+ ∇n_free`. One new scalar on P1 medium; sign flips at `n_free=n_crit`. Gauss + `E_star` untouched. | S1 grid bag `min sep 0.74` collapse; vt `revs ≤0.29`; FINDINGS §L3-incomplete: *"Need capacity/saturation core (shape-3 slot)."* |
| **2** | **P2 ℓ with sign flip**: green hyperbolic ℓ (`ell_max_hist=6.32`, `t_arrival/(r/c)=0.84`, GRIN-kill PASS) is currently configured as thickener/attractor. Flip sign so `∇ℓ` pushes *out* of lock. Same physics insight as rank 1, different field, **already built** → fast fallback. | P2 E2: `Δsep_maxwell=−13.92` vs `+ℓ=−14.00` ⇒ ℓ-as-attract is **decorative-in-practice**. Sign flip is the untested hypothesis. |
| **3** | **Topological core**: phase-winding on lock struct → winding-number exclusion. Honest hard core, **non-diffusive** by construction. Bigger build; fallback if depletion diffuses badly at long T. | K-L2 `E_star` exact over T=200 ⇒ structural core consistent with ledger. |
| **4** | Pure P2 ℓ-attract (status quo). Demote. | E2 collapses identically with/without ℓ. |
| — | Multi-fab L, Cosserat Higgs, pairwise Coulomb: closed / wrong tool / kill-gate. | OP §1. |

**Drop grid bag and pairwise anti-lock as primary instruments.** Both are monotonic in r; neither can produce a band by construction.

---

## 2. First 2-week experiment (P4-seq, sandbox only, no kernel touch)

Extend `v81/path1_locks/src/locks_pic.c`. New state: `n_free(x)=n0−Σ_i w_i(x)` (CIC footprints). New force: `F_core=+k_core (1−n_free/n_crit)_+ ∇n_free` (sign chosen to push *out* of the depleted zone).

### Phase A (week 1) — minimal falsification: does the sign change exist?

Pinned probe pair, `D∈{4,6,8,12,16,20,24,30}`, measure `a_rel(D)` along line of centers. Scan `k_core∈{0.5,2,8} × n_crit/n0∈{0.3,0.6}` (6 combos).

- **Promote:** `a_rel(D)` is **non-monotonic** with a zero crossing at some `D*∈(1,20)`. Report `D*` per combo.
- **Kill (any):** `a_rel` monotonic over full range; OR zero crossing only at `D*<1` (sub-cell, unphysical); OR `gauss_max>1e-2`; OR `E_star` drift.

This is a 1D pinned scan — hours of CPU, decisive, falsifies before any orbit run. The prior bags never isolated this test.

### Phase B (week 2) — orbit, only if Phase A green

`N=96 L=24 T=200` (matches `s1_gridbag` exactly → direct numerical compare). `v_t∈{0.12,0.18,0.25}` × 3 `{k_core,n_crit}` bracketing `D*` (9 runs).

- **Promote:** `sep(t)` bounded in band with `min sep>0.6` (clears S1 `0.74` floor), `revs≥1`, `gauss_max≤2e-3`, `E_star` exact, `E_em` floor held.
- **Kill (any):** all 9 combos pass-through (`min sep<0.5`) or expand (`Δsep>+3`); OR ledger drift; OR band only at stiffness so high `revs<0.5`.
- **Promote → port:** band held full `T=200`, `revs≥1`, clean ledger → port `scp_locks` capacity hook, GPU N=256.

Wall: Phase A hours; Phase B ~9×~50 min ≈ one CPU-day at N=96.

**Why this differs from S1 grid bag.** S1 deposited a bag density `B(x)` and applied `−κ∇B` *toward* the lock — monotonic toward center by construction. P4-seq applies `+k∇n_free` *away* from the lock where `n_free` is depleted — **opposite sign in overlap**. Same scalar, opposite geometry. That is the only thing being tested.

---

## 3. Honest second scale

The four candidates are **not co-equal**:

- **Sequestration / capacity** — same object in this fabric. When overlapping footprints drive `n_free→n_crit`, rest-energy can no longer be claimed there → pressure wall. (i) Already in ontology (v73 process-form ledger; OP shape-3 empty slot). (ii) Turns on at short range by construction — the exact missing sign. (iii) Scalar ⇒ Gauss and `E_star` untouched. (iv) Produces a real second length `r_seq` where `n_free=n_crit`. **This is the honest answer** — the classical saturation analog that gives atoms a scale in non-quantum models (Lennard-Jones hard core = Pauli exclusion).
- **Magnetic (`v×B`)** — velocity-dependent, not position-dependent. Sets *gyrofrequency*, not *radius*. Necessary as angular-momentum holder once a wall exists; cannot create the radius alone. **Do not nominate as the second scale.**
- **Topological (winding core)** — honest length, honest hard core, Gauss-clean, non-diffusive. Different lock ontology (phase-charge on struct); bigger build. Rank-3 fallback.
- **Hard-core capacity** — form of sequestration; same mechanism, same rank.

**The real framing.** Real positronium has a stable radius because quantum uncertainty sets the scale against `1/r²` attraction. This fabric is classical. The "second scale" question is precisely: *what is the classical substitute for the uncertainty principle?* Saturation via depletion is the cleanest classical analog; magnetic alone is not, because it sets `ω` not `r`.

**Honest caveat to log.** Diffusive `n_free` will, at long `T`, either refill (wall softens → late collapse) or freeze (no gyration). Phase B's `T=60→200` probe must measure stiffness drift, not just final `sep`. If no `{k_core,n_crit}` holds a band across both, kill P4-seq honestly → trigger topological core (rank 3, non-diffusive by construction).

---

## 4. Open v82 — **yes, as `path4_capacity`**

- v81/NEXT decision rule is explicit: **S1–S2 FAIL → council → v82**. They failed (`revs ≤0.29`). Open v82.
- **New primary path name: yes.** Call it **`v82/path4_capacity/`**, not `path1_locks` again. The load-bearing ingredient (depletion field, sign change at short range) is new and independently falsifiable; it deserves its own scorecard so P1's clean L0–L2 GREEN (`gauss_max~1e-13` in-kernel after J-fix) stays unmuddied as the characterized attractive baseline.
- If rank-2 (ℓ sign-flip) is the variant that lands, name it `v82/path5_ell_core/` so ℓ-as-depletion is falsifiable separately from P2's ℓ-as-thickening.
- **Run rank-2 in parallel if bandwidth allows** — same physics insight, different field, pre-built. If P4-seq fails Phase A, ℓ-sign-flip is the immediate next experiment, not a topological rebuild from scratch.
- Sandbox-first. No `scp_sim` / `scp_locks.h` edit until Phase A green. Mirror v81 structure: `path4_capacity/{src,FINDINGS.md,STATUS.md,KERNEL_DELTA.md}` + this council dir.

---

## Inherited kill gates (still binding)

No multiplet-only light; no particle-as-Q-ball-profile; no inserted pairwise Coulomb; no multi-fab L; no GRIN variable-`c`; durability + ledger honesty before spectroscopy. P4-seq clears all six by construction: state = free Yee medium + lock structs + one scalar `n_free`; core force is gathered grid pressure, not `1/r²`; `c` fixed; `E_star` ledger untouched.
