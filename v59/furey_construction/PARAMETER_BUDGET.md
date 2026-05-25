# v59 parameter budget — what's free, what's now dependent

*Honest accounting as of 2026-05-24.  Status tags:
**[thm]** = machine-checked Lean identity; **[emp≈X]** = empirical match of a structural number
to data at precision X; **[conj]** = conjecture, not derived; **[free]** = genuine input.*

## The cleanest win: the charged-lepton sector

The 3 charged-lepton masses are parametrised (Brannen) by `(a, t, φ)` = (scale, amplitude, phase):

| d.o.f. | v59 status | turned dependent? |
|---|---|---|
| amplitude `t` | `t²=1/2 ⟺ Koide Q=2/3 = dimG₂/dimSpin7` **[thm for the algebra; emp≈10⁻⁵ for the match]** | **yes** → `Q` |
| phase `φ` | `φ = Q/3 = 2/9` **[emp≈10⁻⁵]**; the `/3` = sedenion `S₃` **[thm-verified automorphism]**; the magnitude `Q` is residual | **yes** → `Q` |
| scale `a` | sets the overall lepton mass scale | **no — free** |

So **3 lepton masses → 1 free scale `a_lepton`** (both shape parameters collapse onto the *same*
structural `Q`).  This session's contribution: it showed `φ` is **not an independent free
phase** — its divisor `/3` is now genuinely structural (sedenion `S₃`), and its magnitude is the
*same* `Q` as the amplitude, so `φ` does not add a parameter beyond `Q`.

## What else turned dependent

| SM quantity | v59 relation | status |
|---|---|---|
| Higgs VEV `v_Higgs` | `= 28²·a_l²` (tied to the lepton scale) | **[emp≈0.07%]** → depends on `a_lepton` |
| `sin²θ_W` | `= 2/9` from Pati-Salam matching `(cW,cBL)=(5,2)` | **[thm given (5,2); cBL=2 pinned, interp. open]** |
| `cos²θ_W` | `= 7/9 = dimImO/9` | **[thm]** |
| `m_W, m_Z` | from `v_Higgs`, `α(M_Z)`, `sin²θ_W` | **[emp≈0.02–0.04%]** |
| gauge prefactor `5` in `g_W²=5√α` | `= h∨(Spin(7))` (dual Coxeter) | **[thm]** |
| quark Koide `Q_d, Q_u` | `=11/15, 23/27` (fix quark amplitudes) | **[emp≈0.3%]** → 2 dependent |
| `α(M_Z)` | `= 25/(324π²)` | **[conj, emp≈0.03%]** |
| α(0), `G_e=(21/16)α²¹` | structural forms | **[conj]** |

## What is still free

| input | status | notes |
|---|---|---|
| `a_lepton` (1 mass scale) | **[free]** | the one genuine dimensionful input of the lepton/EW/Higgs sector |
| `α` (fine structure) | **[conj→borderline]** | `α(M_Z)=25/(324π²)` conjectured; if rejected, `α` is free |
| `Q = 2/3` magnitude | **[emp, residual]** | the residual of `φ=Q/3` (the `/3` is structural; this `2/3` is not symmetry-fixed — proven not a holonomy) — *but it's the same `Q` already used for the amplitude, so not a new input* |
| quark scales `a_u, a_d` (2) | **[free]** | |
| quark phases `φ_u, φ_d` (2) | **[free]** | shown `≠ Q_q/3` — genuinely free, not structural |
| CKM (3 angles + 1 phase) | **[free / mostly unaddressed]** | only `sin²θ_C ≈ 7α` conjectured |
| neutrino masses + PMNS | **[free / unaddressed]** | |
| `θ_QCD`, `α_s` | **[unaddressed]** | |

## The count

- **Charged-lepton + EW + Higgs sector** (the part v59 actually engages):
  - *Optimistic picture* (accept the structural relations + the `α` conjecture):
    **1 free parameter — `a_lepton`** — everything else (3 lepton masses, `v_Higgs`, `sin²θ_W`,
    `m_W`, `m_Z`) is dependent.  This is the long-standing "only `a_lepton`" headline.
  - *Rigorous-only* (theorems + strong empirical matches, **drop conjectures**):
    **2 inputs — `a_lepton` and `α`** — with the same set dependent (now `α` is an input, not
    `25/(324π²)`).
- **Full Standard Model:** still **~10–12 free** — the lepton/EW reduction does **not** touch the
  quark scales (2) + quark phases (2) + CKM (4) + neutrinos (≥7) + `θ_QCD` + `α_s`.  v59 has
  genuinely reduced the **lepton + EW + Higgs** block to ~1 input; it has **not** reduced the
  quark-flavour / CKM / neutrino / strong sectors.

## What this session turned dependent (the delta)

We did **not** lower the raw count, but we:
1. **Pinned `φ` to `Q`** (so it is not an independent free phase): the `/3` is now a verified
   sedenion `S₃` automorphism, the magnitude is the same `Q` — removing a potential
   double-count and making "the lepton sector has one scale + one ratio `Q`" precise.
2. **Converted several coincidences into theorems** — `Q=dimG₂/dimSpin7`, the grade structure
   (lepton=L, `J_c∈Λ²`), `sin²θ_W=2/9` from `(5,2)`, `5=h∨(Spin7)`, `cBL=2` pinned by `2/9` —
   raising confidence that the "dependent" entries are structural, not accidental.
3. **Bounded the residual sharply**: the only genuinely un-derived piece of the lepton sector is
   the *magnitude* `Q=2/3` of the phase (= the same `Q`), and we proved it is **not** a
   symmetry/holonomy phase (not π-rational) — so it is a single mass-sector input, not a family
   of unknowns.

## Honest bottom line

- **Free now (lepton+EW+Higgs):** `a_lepton` (and `α` unless you accept its conjecture).
- **Turned dependent:** the 2 non-scale lepton mass d.o.f. (→`Q`), `v_Higgs`, `sin²θ_W`,
  `cos²θ_W`, `m_W`, `m_Z`, the `5` and `cBL=2` gauge integers, and the quark Koide ratios — at
  precisions from `10⁻⁵` (Koide/phase) to `~0.3%` (quark Koide), with the algebraic identities
  machine-checked.
- **Not touched:** quark scales/phases, CKM, neutrinos, strong sector — the bulk of the SM's
  flavour parameters remain free.
