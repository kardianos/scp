# V40 Search Log

## Gen 0: Baseline (Complete)
- Single V34 braid, N=128, L=15, T=200
- S_final=0.81, S_mean=0.62
- E_pot oscillates -0.2 to -174 (period ~50)
- 41/42 frames alive, P_int retention 86%

## Gen 1: Structural Variations (Complete)
- Date: 2026-03-23
- T=20 probe, N=128, L=15
- 8 candidates testing different structural modifications

### Final Results (ranked by S_final)

| Rank | ID  | Description                     | S_final | S_mean | P_ret  | E_pot   |
|------|-----|---------------------------------|---------|--------|--------|---------|
| 1    | 004 | Double amplitude (scale 1.5)    | 1.4751  | 1.0671 | 0.633  | -232.5  |
| 2    | 003 | Braid + counter-braid           | 1.2404  | 0.9192 | 0.763  | -178.5  |
| 3    | 002 | Two braids perpendicular (90d)  | 1.0822  | 0.7655 | 0.770  | -141.2  |
| 4    | 006 | Braid + oscillon                | 0.7580  | 0.5782 | 0.729  | -86.0   |
| 5    | 001 | Two braids parallel (sep=8)     | 0.7321  | 1.0898 | 0.544  | -82.2   |
| 6    | 007 | High ellipticity (0.5)          | 0.5620  | 0.7082 | 0.605  | -47.5   |
| 7    | 008 | Low background (A_bg=0.03)      | 0.5268  | 0.6729 | 0.619  | -43.2   |
| 8    | 005 | Reduced amplitude (scale 0.7)   | 0.1425  | 0.2351 | 0.224  | -3.1    |

Baseline comparison: S_final=0.81, S_mean=0.62, P_ret=0.86

### Analysis

**Winners (beat baseline S_final=0.81):**
- **004 (scale 1.5):** Highest S_final=1.48 but only 63% P_int retention. Higher amplitude
  drives deeper potential wells. S_mean=1.07 indicates consistently strong binding.
- **003 (counter-braid):** S_final=1.24 with 76% retention. Counter-chirality pair creates
  stronger interaction than same-chirality. Good stability/binding balance.
- **002 (perp braids):** S_final=1.08 with 77% retention. Orthogonal braids create
  constructive interference at the overlap region. Earlier T=10 data showed dispersion,
  but the T=20 snapshot catches a breathing peak.

**Notable observations:**
- **001 (parallel braids):** High S_mean=1.09 but low S_final=0.73. Two separated braids
  start strong but separation causes energy leakage. Interesting: S peaked at ~1.5 at t=10.
- **005 (scale 0.7):** Effectively collapsed. Amplitude too low to sustain binding (P_ret=22%).
- **007, 008:** Modified single braids underperformed baseline. Lower amplitude/changed
  geometry weakens the structure.
- **006 (braid+oscillon):** Moderate. Oscillon provides some constructive interference.

**Key insight:** Higher amplitude and multi-structure configurations dominate. The potential
V(P) = mu*P^2/(1+kappa*P^2) with mu=-41.345 rewards larger amplitudes pushing P deeper
into the attractive well.

### Survivors for Gen 2
1. **004** - Double amplitude (scale 1.5)
2. **003** - Braid + counter-braid
3. **002** - Two braids perpendicular
4. **006** - Braid + oscillon

### Gen 2 Plan
Refine the top 4:
- Scale variations around 1.5 (try 1.3, 1.7, 2.0)
- Combine winners: counter-braid at 1.5x amplitude
- Perpendicular braids at higher amplitude
- Test phase-shifted initial conditions
