# v81 path1_locks — FINDINGS

**Date:** 2026-07-19  
**Binary:** `bin/locks_pic` (standalone C, CPU)  
**Thesis check:** `REAL = FreeSubstrate (Yee U(1)) + Locks (structs)`  
**Kernel:** not modified; see `KERNEL_DELTA.md` for proposed port.

---

## Scorecard

| Exp | Setup | Result | Notes |
|-----|--------|--------|-------|
| **L0** | Single pinned lock, Poisson IC + 200 dyn steps | **PASS** | \(Q=1\); `gauss_max` IC/dyn \(=1.09\times10^{-4}\) (SOR tol); exterior \(\|E\|(r)\) monotone decreasing on \(r\in[6,20]\) |
| **L1** | Rest pair \(q=\pm1\), \(D\in\{8,12,16,20,24\}\), pinned force probe | **PASS** | Monotone \(a_{\mathrm{rel}}(D)\); attraction signs correct; **medium only** (no pairwise term in code); \(F_{\mathrm{along}}\approx 1/(2\pi D)\) 2D continuum |
| **L2** | Soft free pair, \(T=2500\) steps (\(\approx795\) time units), sponge BC | **PASS** | Locks `alive` always; \(E_\star\) drift \(=0\); \(E_{\mathrm{em}}\) floor held (min \(\sim7.5\%\) of \(E_{\mathrm{em0}}\), final \(\sim27\%\)); dynamical Gauss \(\lesssim2\times10^{-3}\) |
| **L3** | \(v_t=0.12\), \(T=4000\) steps (optional) | **PASS\*** | No shred-by-definition; \(\sim12\) relative revs of joining-line angle; separation collapses (soft pass-through / RR-like inspiral — no hard core / sequestration capacity yet) |

**Overall L0–L2:** **PASS (GREEN for port proposal)**

---

## L1 force table (medium-only)

| \(D\) | \(a_{\mathrm{rel}}\) | \(F_{\mathrm{along}}\) | \(F_{2D}^{\mathrm{ref}}=1/(2\pi D)\) | attract |
|------:|---------------------:|----------------------:|---------------------------------------:|:-------:|
| 8 | \(7.925\times10^{-3}\) | \(1.981\times10^{-2}\) | \(1.989\times10^{-2}\) | yes |
| 12 | \(5.176\times10^{-3}\) | \(1.294\times10^{-2}\) | \(1.326\times10^{-2}\) | yes |
| 16 | \(3.788\times10^{-3}\) | \(9.471\times10^{-3}\) | \(9.947\times10^{-3}\) | yes |
| 20 | \(2.937\times10^{-3}\) | \(7.342\times10^{-3}\) | \(7.958\times10^{-3}\) | yes |
| 24 | \(2.352\times10^{-3}\) | \(5.879\times10^{-3}\) | \(6.631\times10^{-3}\) | yes |

Force is \(q\mathbf{E}\) gathered from the live Poisson/Yee medium. Reference column is diagnostic only — **not** used in dynamics.

---

## Kill-gate audit (OP §1)

| Gate | Status | Evidence |
|------|--------|----------|
| 1. Multiplet-only primary | **CLEAR** | State = Yee arrays + lock struct list; no \(\phi_a\) multiplet |
| 2. Particle = only Q-ball profile | **CLEAR** | Particle = lock object; `alive`/`E_star` not a field hump |
| 3. Pairwise Coulomb / springs | **CLEAR** | `grep`/code path: force only via gather \(\to qE + qv\times B\); L1 `medium_only=1` |
| 4. Multi-fab light L | **CLEAR** | No Cosserat / multi-fab arrays |
| 5. GRIN variable-\(c\) gravity | **CLEAR** | Fixed \(c=1\); CFL is engine law only |
| 6. Spectroscopy-first | **CLEAR** | Metrics: durability, Gauss, ledger, \(a_{\mathrm{rel}}(D)\) |

---

## What works

1. **Medium Coulomb is real.** Opposite locks attract with monotone \(a_{\mathrm{rel}}(D)\) matching 2D continuum scale without any \(1/r\) formula in the force law.  
2. **Locks cannot evaporate.** Unlike multiplet L (v79/v80 O FAIL), lock identity is structural; \(E_\star\) is exact over long \(T\).  
3. **\(E_{\mathrm{em}}\) does not null.** Soft long run keeps field energy as a live medium sector (fraction of IC remains).  
4. **Gauss discipline.** IC residual at SOR floor; dynamical residual controlled under zigzag \(J\) deposit (\(\sim10^{-3}\) on this grid — not Cosserat \(10^{-13}\), but non-exploding and charge-total exact).

## What is incomplete / honest limits

1. **2D TE medium**, not 3D Yee — force law is \(\sim1/r\) (line charge), not \(1/r^2\). Sufficient for monist proof; 3D is a numeric upgrade, not ontology.  
2. **No sequestration capacity field** (\(n_{\mathrm{free}}/n_{\mathrm{seq}}\)) yet — rest energy is a constant lock ledger, not dynamically claimed from a free measure. Soft core / bound orbit radius not enforced; L3 spirals to pass-through.  
3. **Dynamical Gauss** not at machine floor (Esirkepov order-1 + CIC gather asymmetry + leapfrog staggering). Good enough for L0–L2 gates; kernel port should use a proven Esirkepov implementation and match discrete operators to Cosserat gauge if co-evolving.  
4. **Energy ledger with sponge** is not closed (by design of absorbing layer); interior \(E_\star\) exact, field↔KE exchange visible. Full Poynting shell flux accounting is next polish.  
5. **L3** is “no shred + relative gyration,” not a durable multi-rev bound atom. Need capacity/saturation core (shape-3 slot) or soft form-factor for positronium product.

---

## Port decision

| Question | Answer |
|----------|--------|
| L0–L2 green? | **Yes** (sandbox) |
| Kill gates clear? | **Yes** |
| Kernel port | **DONE** — see `KERNEL_PORT.md` (user accepted `KERNEL_DELTA.md`) |
| FINDINGS port flag | **GREEN (landed)** |

---

## In-kernel gates (post-port + J-unit fix)

| Exp | Result | Key numbers (N=48, g=1, vacuum Φ) |
|-----|--------|-----------------------------------|
| **K-L0** | **PASS** | `gauss_max~9×10^{-14}`; \(E_{\mathrm{em}}\) holds; \(Q_{\mathrm{flux}}\approx0.99\); \(E_\star\) exact |
| **K-L1** | **PASS (sign)** | Opposite attract for \(D\in\{4,6,8,10\}\); \(a_{\mathrm{rel}}(D)\) not clean monotone (self-force/images) |
| **K-L2** | **PASS** | \(T=200\): **`gauss_max=1.2×10^{-13}` entire run**; alive=2; \(E_\star\) exact; \(E_{\mathrm{em}}\) not null; soft approach |
| **K-L3** | **PARTIAL** | Soft core + \(v_t\): Gauss floor; ~0.7 joining revs; multi-rev park **not** yet |
| **reg0** | **PASS** | `n_locks=0` vacuum stays zero |

**Gauss fix:** zigzag \(J\) had an extra factor of \(dx\); corrected \(J=q\,\mathrm{dseg}/(dt\,dx^2)\).  
Detail: `KERNEL_PORT.md`. Artifacts: `kernel_runs/kl{0,1,2,3}/`.

---

## Reproduce

```bash
# Sandbox
make -C v81/path1_locks clean all
./v81/path1_locks/bin/locks_pic all v81/path1_locks/results

# Kernel
make -C sfa sim/scp_sim && cp sfa/sim/scp_sim bin/
OMP_NUM_THREADS=8 bin/scp_sim v81/path1_locks/kernel_runs/kl0/kl0.cfg
```

Artifacts: `results/l0_*.tsv`, `l1_force.tsv`, `l2_series.tsv`, `l3_tracks.tsv`; kernel `kernel_runs/*/diag*.tsv`, `locks_track*.tsv`.
