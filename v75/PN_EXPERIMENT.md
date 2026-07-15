# P/N (Z/N) Experiment Matrix — F17

**Date:** 2026-07-15  
**Goal:** Firm proton- vs neutron-analog distinction so nuclear \(Q_\mathrm{em}\propto Z\) and \(N\) adds mass without charge.  
**Auth:** equation/const redesign allowed for P/N.

## Physics gap (pre-experiment)

| Approach | Why it fails or may work |
|----------|---------------------------|
| **Flavored multiplet alone** (v71) | \(\Delta\omega\) partitions \(Q_a\) (uud/udd-like) but **all components same-sign EM** under diagonal U(1). Net \(Q_\mathrm{em}\neq 0\) for both branches. |
| **Opposite-\(\omega\) “cancel” multiplet** | Net Noether can cancel; risk annihilation / no stable bag. |
| **B1 lock** | \(\Phi_Q\equiv\Phi_C\) ⇒ every nuclear ball sources EM. No neutral nuclear species. |
| **B2 unlock + C-only n** | \(q_C=0\), \(q_Q=+1\): **p** = C+Q co-located; **n** = C-only. \(\rho_\mathrm{em}=q_Q\rho_Q\). **Primary path.** |

## Kernel / seed changes (this campaign)

| Item | Role |
|------|------|
| `init_sfa_Q`, `mf_seed_Q` | Seed charge fabric independently (or leave zero for neutrons) |
| `mf_lock_CQ=0` | Evolve Q free of C (B2) |
| Diag `Q_C,Q_Q,Q_L,Q_em` (CPU) | Report fabric Noether + EM charge |
| GPU Gauss B2 | Primary density from Q when unlocked (`q_C=0`) |
| `gen_pn_core` | Stamp Z on C+Q, N on C only, optional L |

## Scorecard

| Metric | Pass |
|--------|------|
| **S1** Neutron \(Q_\mathrm{em}/Q_\mathrm{bag} \ll 1\) (target \(\lvert Q_\mathrm{em}\rvert/\lvert Q_C\rvert < 0.05\)) | n is EM-neutral |
| **S2** Proton \(Q_\mathrm{em}\approx Q_Q\approx Q_C\) (same-seed) | p is charged |
| **S3** \(Q_\mathrm{em}(Z,N)\approx Q_\mathrm{em}(Z,0)\) within 5–10% | N does not add EM |
| **S4** Bag survives T≳80 (s_max, phi_max stable; no full dispersal) | species exist as lumps |
| **S5** Flavored-only: both branches \(Q_\mathrm{em}\) large | documents flavor ≠ EM Z/N |

## Run IDs

| ID | Relation | Setup | Expected |
|----|----------|-------|----------|
| F17a | flavored baseline | n_fab=1, sym ω=1.42 | large Q_em |
| F17b | flavored p-branch | ω=(1.38,1.42,1.42) | large Q_em, Q_a unequal |
| F17c | cancel multiplet | ω=(1.42,1.42,−1.42) | Q_em≈0?, unstable? |
| F17d | B2 neutron | lock=0, Q=0, C ball | Q_em≈0, bag OK |
| F17e | B2 proton | lock=0, mf_seed_Q=1 | Q_em≈Q_C |
| F17f | Z1 N1 | gen_pn_core | Q_em≈1×Q_ball |
| F17g | Z2 N0 | two p | Q_em≈2× |
| F17h | Z2 N2 | two p + two n | Q_em≈F17g; mass larger |

Grid: N=96 L=18 T=100 single; N=128 L=28 multi. g=0.05, m_θ=1.6, η=0, absorbing.

## Results (2026-07-15 GPU campaign)

| ID | \(Q_\phi\) end | \(Q_\mathrm{flux}\) end | \(E_\mathrm{em}\) | Q ret | Verdict |
|----|----------------|-------------------------|-------------------|-------|---------|
| F17a sym | 209.5 | 209.3 | 1.00 | 1.00 | charged baseline |
| F17b flav | 284.0 | 283.6 | 1.68 | 1.00 | charged (flavor ≠ Z/N) |
| F17c cancel | 69.9 | 69.8 | 0.11 | 1.00 | residual 1/3; not n |
| **F17d n** | **209.6** | **0.000** | **0** | 1.00 | **PASS S1** |
| **F17e p** | **209.6** | **209.2** | **1.00** | 1.00 | **PASS S2** |
| F17f Z1N1 | 439.5 | 205.9 | 1.13 | 0.94 | \(Q_\mathrm{em}\sim 1p\) |
| F17g Z2N0 | 439.5 | 422.7 | 3.80 | 0.94 | \(Q_\mathrm{em}\sim 2p\) |
| F17h Z2N2 | 1128 | 387.5 | 3.09 | 0.97 | flux≈Z2; E≃2.4×g |

**Frozen path:** B2 unlock + selective Q seeding. Flavored multiplet alone is **not** the isotope knob.
