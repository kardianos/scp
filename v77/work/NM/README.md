# work/NM — Numeric Matter — multi-lock dual-source

**Agent:** NM  
**Phase 2 focus:** **RC1 co-field** (CP-RC1-NUM)

## RC1 co-field (current — CP-RC1-NUM)

| Path | Role |
|------|------|
| `sandbox_rc1_cofield.py` | Fixed multi-locks; 2D F1 ψ + NE `Maxwell2D` dynamical E,B |
| `run_rc1.py` / `run_rc1.sh` | Runner |
| `outputs/rc1_result.json` | `rc1_claim` + RG0–RG6 |
| `outputs/rc1_tm_gates.json` | TM gate package |
| `outputs/rc1_forces.tsv` | Force taxonomy |
| `outputs/rc1_summary.txt` | Human summary |

```bash
cd /home/d/code/scp/v77/work/NM
python3 run_rc1.py
# or: python3 sandbox_rc1_cofield.py --nx 48 --steps 24
```

**Composition:** same 2D grid — free Laplace \(\psi\) from \(\rho_b\) + Yee TE+TM \((\mathbf{E},\mathbf{B})\) from \(\rho_Q\) (import NE `Maxwell2D`). Not Φ-only.  
**Tags:** `E_origin=free_maxwell_full`, `em_solver=free_maxwell_yee_m1`, `phi_origin=free_capacity_f1`.

## Round 3 — unified composition

| Path | Role |
|------|------|
| `composition_r3.md` | Document: L1 dual joint + L2 full Maxwell sibling free-gauge |
| `export_r3_unified.py` | Import NE wave c + NM dual forces → one JSON |
| `outputs/r3_unified_export.json` | Unified package for TU/O |
| `outputs/r3_unified_summary.txt` | Human summary |

```bash
python3 /home/d/code/scp/v77/work/NM/export_r3_unified.py
```

**Verdict:** `composition_pass=True` — supports O unified (A) path; residual no single 3D Yee+ψ grid.

## Round 2 — joint dual-channel

| Path | Role |
|------|------|
| `sandbox_r2_dual_channel.py` | Joint ψ (F1) + Maxwell-lite Φ/E; KG7/KG8; configs N/S/O/V |
| `outputs/r2_dual_result.json` | Full package + V77-2 gates |
| `outputs/r2_dual_forces.tsv` | Force taxonomy |
| `outputs/r2_dual_radial.tsv` | ψ, Φ, E_r vs r |
| `outputs/r2_dual_gates.tsv` | Gate booleans |
| `outputs/r2_dual_deficit.tsv` | Ledger / c / forces |
| `outputs/r2_dual_summary.txt` | Human summary |

```bash
python3 /home/d/code/scp/v77/work/NM/sandbox_r2_dual_channel.py
python3 /home/d/code/scp/v77/work/NM/sandbox_r2_dual_channel.py --N 16 --iters 200  # mini SOR
```

### Laws (one medium)

| Channel | Law | Source | Tags |
|---------|-----|--------|------|
| Free capacity | \(-\sigma_0\nabla^2\psi = s\rho_b\) | \(\rho_b\) | `psi_origin=free_capacity_3d_green` |
| Maxwell-lite | \(-\varepsilon_0\nabla^2\Phi=\rho_Q\), \(\mathbf{E}=-\nabla\Phi\) | \(\rho_Q\) | `E_origin=free_maxwell_lite` |
| Shared \(c\) | \(c=1/\sqrt{\varepsilon\mu}=C_{\mathrm{LOCAL}}\) | — | JC1 |
| Budget | \(\rho_f+\rho_b=\rho_0\) | — | JC2 |
| Support | \(\mathrm{supp}\|\rho_Q\|\subseteq\mathrm{supp}\,\rho_b\) | locks | TM dual-source |

**TE-IA1:** \(\psi\neq\Phi\) (siblings, not clones).

### Key numbers (A=0.35, σ=0.9, sep=4, ε=μ=1)

| Config | \(\psi_{\mathrm{mid}}\) | \(\Phi_{\mathrm{mid}}\) | \(F_\psi\) signed | \(F_C\) signed |
|--------|------------------------:|------------------------:|------------------:|---------------:|
| vacuum | 0 | 0 | 0 | 0 |
| neutral | 0.3198 | 0 | −0.0803 | 0 |
| same-sign | 0.3198 | 0.3198 | −0.0803 | +0.0803 |
| opposite | 0.3198 | 0 | −0.0803 | −0.0803 |

| Multipole | Result |
|-----------|--------|
| exterior \(\psi\) | \(\sim 1/r\) |
| same-sign \(E_r\) | \(\sim 1/r^2\) |
| opposite \(\Phi\) | dipole \(\sim 1/r^2\) |

### V77-2 verdict

**`V77_2_joint_numeric = PASS`** (Maxwell-lite dual-channel joint numeric)

- KG7 sibling PASS, KG8 shared \(c\) PASS, JC1–JC5 PASS  
- `full_maxwell_claim=false` — static lite; E,B hook for NE full Maxwell  

---

## Round 1 — multi-lock free-capacity (baseline)

| Path | Role |
|------|------|
| `offline_r1_multilock.py` / `sandbox_r1_multilock.py` | ≥2 locks; free_gauge_lite dual stub |
| `outputs/r1_*` | G0–G3 TM gates PASS |

## Sector tags

| Tag | Meaning |
|-----|---------|
| `monist_1sector` | free-capacity ψ only (R1) |
| `multi_channel` | ψ + free_maxwell_lite Φ/E (R2 joint) |
