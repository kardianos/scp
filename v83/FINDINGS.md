# v83 FINDINGS — execution of Fable-informed plan

**Date:** 2026-07-20  
**Plan:** `v83/PLAN.md`  
**Stop reason:** **E1a FAIL** (hybrid nucleus↔lock force has wrong Coulomb sign)

---

## Decision (own plan, Fable-informed)

- Run **E1 first** (droplet/ball attach), not 2B unit first.  
- Stage 3 positronium not gated on 2B.  
- Stop at first hard FAIL on critical path.

---

## E1a — Hybrid force scan (executed)

**Setup:** N=48, L=12, gauged Q-ball seed (`f_w142`, ω=1.42, Q≈315), `g=0.05`, full Cosserat+gauge (`locks_medium_only=0`), free heavy lock m=500 (pinned skips EM `f[]` in kernel — cannot use pins for force readout).

| D | q_lock | ⟨F_x⟩ | Interpretation (lock at +x) |
|---|--------|-------|-----------------------------|
| 5 | −1 | +2.12e−3 | **outward** |
| 6 | −1 | +1.49e−3 | outward |
| 8 | −1 | +7.48e−4 | outward |
| 10 | −1 | +2.92e−4 | outward |
| 8 | **+1** | **−8.92e−4** | **inward** |
| 8 | vac self | −8.7e−5 | small baseline |

- **Gauss:** `gauss_max ~ 7e−14` entire runs (floor).  
- **Nucleus:** stable; Q≈315; E_em holds.  
- **Scaling:** |F| falls with D; ratio to continuum \(g^2 Q/(4\pi D^2)\) ≈ 0.76–0.85 at D=5–8 (order-correct magnitude).  
- **Sign:** same-sign lock **attracted** to positive ball; opposite lock **repelled**. That is **inverted Coulomb**.

**Control (pure locks, prior kl2):** opposite locks attract under the same `F=-(g q)E` formula. So lock↔lock sign is the intended attract/repel; **matter ρ → E → lock force is flipped relative to lock ρ → E → lock force**.

### Gate

| Criterion | Result |
|-----------|--------|
| Hybrid coupling runs (matter + lock + Gauss) | **PASS** (implementation smoke) |
| Force magnitude Coulomb-class | **PASS** (rough) |
| Force **sign** correct for Stage-4 attach | **FAIL** |
| Overall E1a | **FAIL** |

**Plan stop:** “no clean force / coupling broken → STOP, interface design.”  
Wrong sign is coupling broken for atom attach.

---

## E1b, E2, E3

**Not run** — blocked by E1a FAIL.

---

## Implications

1. **H1 strong form not tested.** Cannot claim droplet is a good Stage-4 target until opposite lock is **attracted**.  
2. **Not a 2B issue.** This is monist hybrid bookkeeping / force convention between Cosserat charge and locks.  
3. **Likely fix locus (not applied):** sign consistency of matter \(\rho_Q\) in Gauss vs lock CIC \(\rho\), or lock force relative to total E when both sectors present. **Do not** blindly flip `F=-(gq)E` without re-checking lock–lock attract. Kernel edit needs an explicit dual regression (kl2 + E1a).  
4. **User’s 2B / anti-lock unit** remains untested; still deferred, not falsified.  
5. **Pinned locks** still omit EM force in `f[]` — instrumentation bug for future pinned scans.

---

## Files

| Path | Role |
|------|------|
| `v83/PLAN.md` | execution plan |
| `v83/e1/seeds/ball_N48.sfa` | single-ball seed |
| `v83/e1/run_e1a.sh` | force scan driver |
| `v83/e1/results/` | tracks, diags, logs, summaries |
| `v83/council/claude_fable/` | prior + Fable memo |

---

## Next (after human authorization)

1. **Fix hybrid force sign** with dual test: opposite locks attract **and** opposite lock↔Q-ball attract.  
2. Re-run E1a → E1b.  
3. Only then E2 (gauged interlock) / E3 (anti-lock unit) as parallel research.
