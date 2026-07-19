# Product scorecard — three layers (no composite \(G\))

**Status:** standing rule for multi-fab product campaigns after v80  
**Replaces:** weighted \(G = \sum w_i S_i\) as a decision device  
**Principle:** integrity and Stage-3 mechanics are **AND-gates**. A strong ledger never rescues a failed force/orbit measurement.

---

## Status alphabet

| Status | Meaning |
|--------|---------|
| **PASS** | Measured; meets bar |
| **SOFT** | Incomplete or borderline — re-test that gate only; **do not** Stage-up |
| **FAIL** | Critical fail — **STOP_AND_RECONSIDER** that path |
| **N/A** | Out of scope this campaign (does not rescue others) |

---

## Layer 1 — Integrity (any FAIL → abort campaign)

| Gate | ID | PASS | FAIL |
|------|-----|------|------|
| Gauss / gauge | **L0** | `gauss_max` ~1e-13 class on physics runs | Floor leaves |
| Diagnostics present | **L1a** | Diags/runlogs for every claimed job | Missing outputs |
| Instrumentation (if force/orbit claimed) | **L1b** | Kept SFAs **and** pair/COM tracks for claimed \(D\)/\(v_t\) | SFAs pruned or tracks never run while claiming force/orbit |

**v80 overnight under L1b:** if scored as a force/orbit claim → **FAIL** (tracks pruned). Re-run with keep-SFA policy before any Stage-up.

---

## Layer 2 — Stage-3 mechanics (each is its own gate)

Advance Stage-3 only when all **in-scope** gates are **PASS**. Any **FAIL** → stop product two-body path and reconsider (e.g. representation line). Any **SOFT** → continue **that gate only**.

| Gate | ID | PASS bar | FAIL bar | SOFT / incomplete |
|------|-----|----------|----------|-------------------|
| Force (opposite) | **F** | Tracks: attract; \(a_{\mathrm{rel}}(D)\) monotone attractive over ≥3 \(D\); fit quality documented | Flat/noise \(a_{\mathrm{rel}}\); wrong sign; merge-only at all \(D\) | Jobs OK, polarity OK, **no tracks** |
| Same-sign control | **R** | Same-sign repel or clearly distinct from opposite | Same-sign attracts like opposite | Missing control |
| Pair dwell | **P** | Rest: \(Q\) held, \(E_{\mathrm{em}}\) not collapsed | Instant \(E_{\mathrm{em}}\) death / annihilation | Coexist but large E drift |
| Orbit / capture | **O** | ≥2 revs **or** clear capture on tracks; radius bounded; \(E_{\mathrm{em}}\) not v79-death | All soft \(v_t\) merge or flyby with \(E_{\mathrm{em}}\) drain | Runs complete, no multi-rev measured |
| Long-T hold (H-class) | **H** | \(Q_\phi\) held; \(E_{\mathrm{em}}\) held; **COM-aware** core/flux explained | \(E_{\mathrm{em}}\) collapse on **minimal** C+L | Global hold OK but fixed-center \(Q_{\mathrm{core}}\)/flux misleading |

**Morphology (M)** is informational only — never blocks alone.

---

## Layer 3 — Stage-4 shell (only after Layer 2 all PASS)

| Gate | ID | PASS | FAIL |
|------|-----|------|------|
| Staged multi-L | **N** | Load 1→2→… with \(E_{\mathrm{em}}\) survival per step; dynamical assembly preferred | Hand-placed multi-L park kills \(E_{\mathrm{em}}\) (v79 recipe) |

Hand-placed Z6+L6 is **FAIL** for Stage 4, not a SOFT.

---

## Verdict formula

```text
if any Layer-1 FAIL:
    VERDICT = STOP_AND_RECONSIDER   # fix method / kernel config
elif any Layer-2 FAIL (in scope):
    VERDICT = STOP_AND_RECONSIDER   # product two-body path in doubt
elif any Layer-2 SOFT (in scope):
    VERDICT = CONTINUE_GATE_ONLY    # re-measure that gate; no Stage-up
elif all Layer-2 PASS (in scope):
    VERDICT = ADVANCE               # may attempt Layer 3 or next Stage-3 notch
else:
    VERDICT = INCOMPLETE_SCOPE
```

**There is no single composite score.** Report a table of gates, not one \(G\).

---

## Ordered campaign policy (post-council)

1. **S4 COM / window re-analysis** (cheap) — confirm H is not fake charge loss  
   → see `S4_COM_WINDOW_ANALYSIS.md`  
2. **Force grid with kept SFAs + `mf_pair_track`** — gate **F**/**R**  
3. **Low-\(v_t\) orbit with tracks** — gate **O** (only if F PASS)  
4. **Staged multi-L** — gate **N** (only if F,P,O,H PASS)  
5. **Representation 2D toy** — parallel CPU; never substitutes for F/O FAIL  

---

## v80 overnight re-score (no tracks claimed)

| Layer | Gate | Status | Note |
|-------|------|--------|------|
| 1 | L0 | **PASS** | Gauss floor |
| 1 | L1a | **PASS** | Diags present |
| 1 | L1b | **FAIL** if force/orbit claimed | SFAs pruned |
| 2 | F | **SOFT** | No \(a_{\mathrm{rel}}\) |
| 2 | R | **SOFT** | Bookkeeping only |
| 2 | P | **SOFT** | Rest OK-ish; E drift |
| 2 | O | **SOFT** | No multi-rev |
| 2 | H | **PASS*** | See S4 analysis: global hold real; fixed-center core/flux are windows |
| 3 | N | **N/A** (v79 FAIL if claimed) | Do not re-park |

\*H PASS for “minimal pair keeps \(E_{\mathrm{em}}\) and \(Q_\phi\)”; not PASS for “fixed-center \(Q_{\mathrm{core}}\) is a good soliton tracker.”

**Verdict:** CONTINUE_GATE_ONLY → instrumentation + force tracks (after S4 window note). **Not** ADVANCE to multi-L.
