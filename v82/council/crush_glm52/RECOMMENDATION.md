# v82 RECOMMENDATION — Stage-3 after P1 multi-rev miss

**Advisor:** crush (glm-5.2) · **Date:** 2026-07-19 · **Stance:** deepen P1 with a capacity core; do not pivot ontology.

---

## Verdict (one line)

P1 is not dead — its **attractive channel works** (L1 `a_rel` matches `1/(2πD)`; N=256 anti-lock cuts `Δsep` from **+10.7 → +7.4**). What it lacks is a **short-range repulsive wall**, so every "second scale" tried so far is either too weak (expand) or monotonic-attractive (collapse-through). Add a **capacity/sequestration core** that turns ON at short range. That is the one honest second scale in this fabric.

---

## 1. Ranked next paths

| Rank | Path | Why | Evidence cited |
|------|------|-----|----------------|
| **# v82 RECOMMENDATION — Stage-3 after P1 multi-rev miss

**Stance:** deepen P1 with a capacity core; do not pivot ontology.

## Verdict
P1's **attractive channel works** (L1 `a_rel` matches `1/(2πD)`; N=256 anti-lock cuts `Δsep` **+10.7 → +7.4**). What it lacks is a **short-range repulsive wall**. Every "second scale" tried is either too weak (expand) or monotonic-attractive (collapse-through). Add a **capacity/sequestration core** that turns ON at short range — the one honest second scale in this fabric.

## 1. Ranked paths

| Rank | Path | Evidence |
|------|------|----------|
| **1** | **P4-seq**: free U(1) + locks + depleted `n_free(x)` → repulsive pressure in overlap zone. First mechanism with **opposite sign at short range** = the only thing that yields a band. Scalar → Gauss + `E_star` ledger untouched. | Grid bag `min sep 0.74`, revs `0.25`; vt `revs ≤ 0.29` pass-through. L3 FINDINGS: *"Need capacity/saturation core (shape-3 slot)."* |
| **2** | **P1+P2 unified, ℓ as depletion (not thickening).** Reuse green ℓ (`ℓ_max≈6.3`, hyperbolic, GRIN-kill PASS); flip sign so `∇ℓ` pushes out of lock. | P2 E2: `Δsep_maxwell=−13.92` vs `+ℓ=−14.00` → ℓ-as-attract is decorative-in-practice. Sign flip is the untested hypothesis. |
| **3** | **Topological core** (phase-winding locks; non-overlap by winding number). Honest, Gauss-clean, bigger build. Fallback. | L3 `E_star` exact ⇒ structural core consistent with ledger. |
| **4** | Pure P2 ℓ-attract (status quo). Collapses (`Δsep −14`); demote. | P2 E2. |
| — | Multi-fab, Cosserat Higgs, pairwise Coulomb: closed / wrong tool / kill-gate. | OP §1. |

Drop grid bag & pairwise anti-lock as primary instruments — both monotonic, cannot produce a band by construction.

## 2. First 2-week experiment (P4-seq, sandbox, no kernel touch)

Extend `v81/path1_locks/src/locks_pic.c`. State: `n_free(x)=n0−Σ_i w_i` (CIC footprints). Core force `F_core = +k_core·(1−n_free/n_crit)_+·∇n_free`.

- **Scan** (≤12 runs): `k_core∈{0.5,2,8}` × `n_crit/n0∈{0.3,0.6}` × `v_t∈{0.12,0.18}`.
- **Phase A:** `N=96 L=24 T=60` (identical to s1_gridbag → direct compare to `null 0.12 / grid 0.25`). **Phase B (if band):** `N=128 L=32 T=200`.
- **Observables:** `sep(t)` band with `r_min>0.6` (clears `0.74` floor); revs `≥1`; `gauss_max ≤ 2e-3`; `E_star` drift `=0`; `E_em` floor held.
- **Kill (any):** no triple gives a bounded band (all collapse-through `min sep<0.5` or expand `Δsep>+3`); **or** `gauss_max>1e-2`; **or** `E_star` drifts; **or** band only at stiffness so high that `revs<0.5`.
- **Promote:** band held full `T=200`, revs `≥1`, `gauss≤2e-3`, `E_star` exact → port `scp_locks` core hook, GPU run.

## 3. Honest second scale

The four candidates are **not co-equal**:
- **Hard-core capacity = sequestration** (same thing here): free measure depleted by overlapping footprints → rest-energy can't be claimed → pressure wall. **This is the honest answer** — already in ontology (v73 process-form, OP shape-3), turns on at short range, scalar (Gauss-clean), produces real second length `r_seq`.
- **Magnetic:** velocity-dependent → sets gyrofrequency, not radius. Cannot create the stable radius alone.
- **Topological:** honest but bigger build; rank-3 fallback.

**Why the bag failed, cleanly:** both bags are **monotonic toward lock center**. A monotonic co-field can only be (a) too weak → `null` regime (`Δsep +5.05`, revs `0.12`), or (b) strong → drives into center → `grid`/`vt` regime (`min sep 0.6–1.2`, pass-through). No sign change ⇒ no band. A depletion core is the minimal sign-change addition without a new gauge charge.

**Honest caveat:** diffusive `n_free` will at long `T` either refill (late collapse) or freeze (no gyration). Scan must probe `T=60→200` stiffness drift. If no `{k_core,n_crit}` holds across both → honest kill → trigger topological core (rank 3, non-diffusive).

## 4. Open v82 — **yes, as `path4_capacity` (P4-seq)**

v81/NEXT rule is explicit: S1–S2 FAIL → council → v82. They failed (`revs ≤ 0.29`). Open it.

- **New primary path name: yes** — `v82/path4_capacity/`, not `path1_locks` again. The capacity field is new and independently falsifiable; giving it its own scorecard keeps P1's clean L0–L2 GREEN (`gauss ~1e-13` in-kernel) unmuddied. P1 = characterized attractive baseline; P4 owns the core.
- If the ℓ unification (rank 2) ships, name it `path5_ell_core` so ℓ-as-depletion is falsifiable separately from P2's ℓ-as-thickening.
- Sandbox-first; no `scp_sim`/`scp_locks.h` edit until Phase A green.

**Inherited kill gates still binding:** no multiplet-only light; no particle-as-Q-ball-profile; no inserted pairwise Coulomb; no multi-fab L; no GRIN variable-`c`; durability + ledger honesty before spectroscopy. P4-seq clears all six by construction.

*(Written to `v82/council/crush_glm52/RECOMMENDATION.md`; sibling seats `agy`/`claude_fable`/`kimi_k3` are currently empty.)*
-equal**:

- **Hard-core capacity = sequestration** (same thing in this fabric): when the local free measure is depleted by overlapping lock footprints, rest-energy can no longer be claimed there → a pressure wall. **This is the honest answer.** It (i) already lives in the ontology (v73 process-form ledger; OP shape-3 slot), (ii) turns on at short range by construction — the exact sign missing from the bag, (iii) is a scalar, so it preserves Gauss and the `E_star` ledger, (iv) produces a real second length scale `r_seq` where `n_free → n_crit`. **This is what was actually missing**, not a new force law.
- **Magnetic (`v×B`):** the medium already has it. It is velocity-dependent, not position-dependent, so it sets a *gyrofrequency*, not a *radius*. Useful as angular-momentum support once a core exists; **cannot create the stable radius alone.** Do not nominate it as the second scale.
- **Topological (winding core):** honest length scale and an honest hard core, and Gauss-clean. But it is a different lock ontology (phase-charge on the lock struct) and a bigger build. Fallback at rank 3, not primary.

**Why the bag failed, stated cleanly:** both the form-factor bag and the grid bag are **monotonic functions of `r` toward the lock center**. A monotonic co-field can only (a) be too weak → the `null` regime (`Δsep +5.05`, revs `0.12`), or (b) be strong enough to matter → it accelerates into the center → the `grid`/`vt` regime (`min sep 0.6–1.2`, pass-through). There is no sign change, so there is no band. **A depletion core is the minimal addition that introduces a sign change at short range without adding a new gauge charge or breaking the ledger.**

**Honest caveat to log:** a diffusive `n_free` will, at long `T`, either refill (wall softens → late collapse) or freeze (wall stiffens → no gyration). The 2-week scan must explicitly probe long-`T` stiffness drift. If no `{k_core, n_crit}` holds a band from `T=60` through `T=200`, that is the honest kill — and the trigger to try the topological core (rank 3), which has a non-diffusive core by construction.

---

## 4. Open v82? — **Yes, as `path4_capacity` (P4-seq)**

- v81/NEXT decision rule is explicit: **S1–S2 FAIL → open council → v82**. S1/S2 failed (`revs ≤ 0.29`, pass-through). Open v82.
- **New primary path name: yes.** Call it **`v82/path4_capacity/` (P4-seq)**, not `path1_locks` again. Reason: the load-bearing ingredient (the capacity field, the sign change at short range) is *new and independently falsifiable* — it deserves its own scorecard so we do not muddy P1's clean "medium Coulomb works" result (L0–L2 GREEN, `gauss_max ~1e-13` in-kernel). P1 stays as the characterized attractive baseline; P4 owns the core.
- **If** the P1+P2-ℓ unification (rank 2) is what ships, name that variant `path5_ell_core` so the ℓ-as-depletion hypothesis is separately falsifiable from ℓ-as-thickening (P2's current, decorative-in-practice reading).
- Keep v82 structure parallel to v81: `path4_capacity/{src,FINDINGS.md,STATUS.md,KERNEL_DELTA.md}` + this council dir. Sandbox-first; no `scp_sim`/`scp_locks.h` edit until Phase A green.

---

## Kill-gate reminder (inherited, still binding)

No multiplet-only light; no particle-as-Q-ball-profile; no inserted pairwise Coulomb; no multi-fab L; no GRIN variable-`c`; durability + ledger honesty before spectroscopy. P4-seq clears all six by construction: state = free Yee medium + lock structs + one scalar `n_free`; core force is gathered grid pressure, not `1/r²`; `c` fixed; `E_star` ledger untouched.
