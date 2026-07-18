# J5-β Numeric Package — Renorm Necessity (ND Round 2)

**Agent:** ND  
**Demo:** `D-DYN-j5-formfactor` / `D-DYN-J5`  
**Status:** **J5-β LOCKED as operational default** (numeric)  
**Partner theory:** TD `dynamics_package_v0.md` §3.5  
**Date:** 2026-07-18  

---

## 0. Claim (killable)

**J5-α (raw):** \(m_{\mathrm{inertial}} = \int\rho_b/c^2 = M_{\mathrm{ray}}\) at default free constants  
→ **FAIL** (form factor).

**J5-β (renorm necessity):**  
1. Naïve \(m_{\mathrm{inertial}}=\int\rho_b/c^2\) is **killed**.  
2. \(m_{\mathrm{inertial}} = \xi\,U[\psi]/c^2\) with \(\xi=1\) (boost model).  
3. \(M_{\mathrm{ray}}=M_{\mathrm{ledger}}\) closed by F1.  
4. Closing all three to one number requires **documented** renorm / ledger tag — not free-fit of \(m\).

**V77-3 bar:** non-tautology **or** renorm necessity + naïve kill.  
**ND R2 verdict:** **J5-β path MET for Dynamics residual documentation**; raw triad not claimed.

---

## 1. Baseline numbers (A=0.4, σ=1, s=κ=1, c=1)

| Symbol | Value | Role |
|--------|------:|------|
| \(M_{\mathrm{ledger}}\) | 6.29984398 | \(\int\rho_b/c^2\) |
| \(U\) | 0.89095341 | free self-energy |
| \(ff=U/Mc^2\) | 0.14142136 | form factor |
| \(m_{\mathrm{field}}=U/c^2\) | 0.89095341 | free-field mass |
| \(m_{\mathrm{boost}}\) (ξ=1) | 0.89095341 | A2 boost |
| \(M_{\mathrm{ray}}\) | 6.29984398 | F1 Born |
| \(s_*\) (per-shape) | 7.07106781 | β1 calibration for this shape |

**Form-factor law** (\(s=\kappa=1\)):
\[
ff = \frac{M}{8\pi\sigma\sqrt{\pi}}
  = 0.141421356237\times\frac{A}{0.4}\times\sigma^2.
\]

---

## 2. ND-G1 — naïve kill (PASS)

| Check | Value |
|-------|------:|
| \(\lvert m_{\mathrm{field}}-M_{\mathrm{ledger}}\rvert/M_{\mathrm{ledger}}\) | **0.8586** |
| threshold for “equal” | 0.25 |
| `kill_naive_m_eq_int_rho_b` | **True** |
| `tautology_flag` | **false** |
| protocol | A2_boost_energy + free-field \(U\) |

---

## 3. ND-G2 — renorm shape invariance (FAIL universal / PASS per-shape)

TD: *single \(s_*\) closes triad for lock family; fail if shape-by-shape fudge only.*

### 3.1 Two-width test (fixed A=0.4)

| σ | \(M\) | \(ff\) | \(s_*=1/ff\) | \(U\) |
|--:|------:|-------:|-------------:|------:|
| 1.0 | 6.29984 | 0.141421 | **7.07107** | 0.89095 |
| 1.5 | 21.26197 | 0.318198 | **3.14270** | 6.76569 |

Apply \(s_*^{(1.0)}=7.07107\) to **both** widths:

| σ | \(U(s_*^{(1)})\) | \(U/Mc^2\) | \(\lvert U/Mc^2-1\rvert\) | triad closed? |
|--:|-----------------:|-----------:|--------------------------:|:-------------:|
| 1.0 | 6.29984 | 1.000 | 0.000 | **YES** |
| 1.5 | 47.842 | 2.250 | 1.250 | **NO** |

Apply \(s_*^{(1.5)}=3.14270\) to both:

| σ | \(U/Mc^2\) | triad closed? |
|--:|-----------:|:-------------:|
| 1.5 | 1.000 | **YES** |
| 1.0 | 0.444 | **NO** |

**`shape_invariance_ok` = false**  
**β1 universal free-constant renorm: FAIL**  
**β1 per-shape (fixed form factor family): PASS**

### 3.2 Closure option tags (export)

| Tag | What it does | Numeric status |
|-----|--------------|----------------|
| `renorm_tag=R1_s` | per-shape \(s_*\) so \(U=Mc^2\) | works; **not** shape-invariant |
| `renorm_tag=ledger_as_U` (β2) | \(E_{\mathrm{ref}}^{\mathrm{(in)}}=U\) for inertia | \(m_{\mathrm{inertial}}=E_{\mathrm{ref}}/c^2\) **PASS**; \(M_{\mathrm{ray}}\) still tracks \(\int\rho_b\) → triad with single \(E_{\mathrm{ref}}\) **FAIL** unless path-cost also retagged |
| `renorm_tag=none` + raw | J5-α | **FAIL** at default |

**ND recommendation (R2):**  
- Default operational tag: **`J5-beta` + `renorm_tag=R1_s` with `shape_scope=fixed_form_factor`**  
- Or **β2** for inertia sector only (`E_star_def=free_U`) with path-cost still on bound ledger (`E_star_def=integral_rho_b`) — dual bookkeeping, monist continuum, **must be accepted by TU** as non-dualist.

---

## 4. Score \(L_{\mathrm{triad}}\)

With \(E_{\mathrm{ref}}=M_{\mathrm{ledger}}c^2\) (raw):
\[
L = \left|\frac{m_{\mathrm{boost}}}{M}-1\right| + \left|\frac{M_{\mathrm{ray}}}{M}-1\right|
  = |0.14142-1| + |1-1| = \mathbf{0.8586}.
\]

With \(E_{\mathrm{ref}}=U\) (β2 inertia):
\[
L_{\mathrm{in}} = \left|\frac{m_{\mathrm{boost}}}{U/c^2}-1\right| = 0,
\quad
\left|\frac{M_{\mathrm{ray}}}{U/c^2}-1\right| = 1/ff - 1 \approx \mathbf{6.071}.
\]

With \(s_*\) per-shape and \(E_{\mathrm{ref}}=Mc^2\):
\[
L = 0 + 0 = \mathbf{0}.
\]

---

## 5. Pass/fail table (R2 locked)

| Gate | Result |
|------|--------|
| Non-tautology | **PASS** |
| ND-G1 naïve kill | **PASS** |
| tier_A ray=ledger | **PASS** |
| tier_C boost=field | **PASS** |
| tier_B / J5-α raw | **FAIL** |
| β1 per-shape renorm | **PASS** |
| β1 shape-invariant \(s_*\) | **FAIL** |
| J5-β documentation | **PASS (locked)** |
| V77-3 full close | **PARTIAL** — needs TU accept of β1-scope or β2 dual ledger |

---

## 6. Files

- R1 numerics: `outputs/j5_formfactor_result.json`, `form_factor_scan.tsv`
- R2 package: `outputs/j5_beta_numeric.json`
- Sandbox: `sandbox_j5_formfactor.py`
