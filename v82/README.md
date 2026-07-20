# v82 — capacity / sequestration core (after P1 multi-rev miss)

**Status:** council + Phase A **PASS** (force sign change); Phase B orbit **not yet** multi-rev  
**Baseline:** v81 P1 locks (commit `66fdfbe`) kept as attractive channel  

## Why

S1 grid bag was non-decorative but collapse-through (`min sep~0.7`, revs≪1).  
Council (2× glm-5.2): second scale = **depletion of free capacity** (repel when footprints overlap), not monotonic bag attract.

## Layout

| Path | Role |
|------|------|
| `DESIGN_BRIEF.md` | Brief to advisors |
| `council/` | glm-5.2 recommendations + synthesis |
| `path4_capacity/` | Sandbox Phase A/B |

## Results (so far)

| Gate | Result |
|------|--------|
| Phase A pinned `F_along(D)` non-monotonic | **PASS** — zero-cross near \(D^*\sim5\) for \(k_{\mathrm{core}}\ge2\) |
| Phase B free orbit multi-rev | **FAIL** (this session) — min sep held (~10–14), revs~0.18, slow expand |

## Next

Tune `{k_core,n_crit,foot_r,v_t}` / longer T; if band+revs≥1 → kernel port; if not → topological core fallback (council rank 3).
