# v81 Path 3 — Token / Update-Budget CA (research spike)

**Status:** MVP runnable · T0–T2 executed (see `FINDINGS.md`, `STATUS.md`)  
**Scope:** Pure 2D sandbox under `v81/path3_token/` only. **Never** edits `scp_sim`.  
**Thesis refs:** `v81/OP.md` §2 P3 · `v80/SHAPES.md` shape 6 · `v81/council/SYNTHESIS.md` rank-3 probe.

---

## Ontology

| Item | Definition |
|------|------------|
| **State** | Directed token occupations `f[y,x,d]` (d ∈ {E,N,W,S}) + optional `rest[y,x]` sequestered bucket |
| **Currency** | Integer tokens = update budget |
| **Hop cap \(c\)** | Hard per-tick max tokens transferred across any directed bond (engine law) |
| **Charge** | Topological circulation: plaquette curl of velocity from directed occupations |
| **Mass** | Sequestered / cycling tokens inside a coherent vortex pattern (not a multiplet hump field) |
| **Not present** | Maxwell/Yee, multiplet \(\phi\), pairwise Coulomb, multi-fab L |

Kill-gate alignment (OP §1):

1. Primary state is free token substrate + emergent circulation locks — **not** multiplet \(\phi\).
2. “Particle” = self-sustaining vortex defect of the token flow.
3. No inserted \(1/r^2\) / springs.
4. No hand-placed multiplet L.
5. \(c\) is a hard hop cap (literal budget), not GRIN \(c_{\mathrm{eff}}\) gravity.
6. First metrics: ledger honesty + durability, not carbon spectroscopy.

---

## Step law

Each tick:

1. **Collide** — integer BGK toward D2Q4 equilibrium; mass exact (largest-remainder); optional HPP head-on scatter for isotropy. Supersonic sites projected to non-negative orthant + renormalized.
2. **Stream** — each direction transfers `min(f[d], c)` one site; excess stays. Periodic BC.
3. **Diagnostics** — total tokens, max bond transfer, vorticity, windowed centroids.

Invariant: `sum(f) + sum(rest)` is bit-exact for all steps when the mass repair path is not forced (and is forced only on float pathology).

---

## Build / run

Requires: Python 3.10+, NumPy.

```bash
cd v81/path3_token
python3 src/run_all.py                  # T0+T1+T2 → results/
python3 src/run_all.py --only T0
python3 src/run_all.py --only T1 --out results
python3 src/run_all.py --only T2
```

Outputs under `results/`:

| Path | Content |
|------|---------|
| `T0/t0_result.json`, `t0_timeseries.csv` | Conservation + hop-cap scorecard |
| `T1/t1_result.json`, `t1_*_timeseries.csv`, `t1_raw.json` | Opp vs same pair |
| `T2/t2_result.json`, `t2_timeseries.csv` | Long durability |
| `summary.json` | Aggregate pass flags |

---

## Experiments

| ID | Setup | Pass criterion |
|----|--------|----------------|
| **T0** | Uniform fill + overloaded east band; \(T=500\) | Exact token total; observed bond transfer \(=c\) and never \(>c\) |
| **T1** | Soft-core vortex pair ± vs ++; track windowed centroids | Circulation-sign-dependent interaction (attraction **or** clear signed force) |
| **T2** | Opposite pair, \(T=4000\) | Tokens exact; \(\|C\|\) and KE retained; centroids tracked (patterns do not evaporate like multiplet humps) |

Default hydro params (T1/T2): \(N=96^2\), \(\rho=24\), \(c=6\), \(\omega=1.6\), core\(=4\), \(\Gamma=8\).

---

## Source layout

```
v81/path3_token/
  README.md
  FINDINGS.md
  STATUS.md
  src/
    token_ca.py      # TokenCA core
    experiments.py   # T0–T2
    run_all.py       # CLI
  results/           # generated
```

---

## Explicit non-goals (this OP)

- Kernel / `scp_sim` / `sfa.h` changes  
- Maxwell medium (that is P1)  
- Path-measure \(\ell\) (that is P2)  
- Multi-fab product, carbon atom packaging  
- Week-scale positronium claim (council: research probe only)
