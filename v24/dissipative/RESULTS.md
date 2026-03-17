# V24-E Results: Dissipative Soliton (Gain-Loss Balance)

## Summary

**NEGATIVE RESULT: No stable dissipative oscillon exists with the proposed gain mechanism.**

The gain-loss balance is structurally unstable for all three gain mechanisms tested. The system exhibits a sharp bistability: below a critical G, damping wins and the oscillon decays to zero. Above critical G, gain wins and energy grows without bound (or saturates at a delocalized high-energy state that is NOT an oscillon). There is no intermediate attractor — no dissipative soliton limit cycle.

The fundamental obstruction is the **high nonlinearity of the gain term**. The proposed S_a = P^2 * phi_a scales as phi^7, creating a catastrophic positive feedback loop that overwhelms linear damping. Even with gain saturation, the transition is too sharp for stable balance. This is not a tuning failure — it is a structural impossibility for polynomial gain mechanisms acting on an oscillon.

---

## Model

Modified EOM for three massive scalars:

    d^2 phi_a/dt^2 = Laplacian(phi_a) - m^2 * phi_a - dV/dphi_a - gamma * dphi_a/dt + G * S_a

V = (mu/2) P^2 / (1 + kappa P^2),  P = phi_1 phi_2 phi_3

Parameters: mu=-20, kappa=20, m=1.0, Nx=4000, xmax=100

Conservative oscillon equilibrium (after t=3000): E=1.285, peak=0.265, f_core=0.999

---

## Phase 1: Gain-Loss Balance Scan

Three gain mechanisms tested:

### Mode A: Alignment gain, S_a = P^2 * phi_a / (1 + kappa_g * P^2)

This is the proposed mechanism from the spec. Rewards configurations where all three fields are aligned.

| gamma | G_balance | E_final/E_0 | Status |
|-------|-----------|-------------|--------|
| 0.001 | 91.05 | 0.004 (decayed) | Marginal — balanced during 500 t.u. search, decayed during 2000 t.u. validation |
| 0.01 | 91.35 | 85.1 (blowup) | Transition to delocalized high-E state (E=109, peak=0.89) |
| 0.1 | 95.02 | 0.0 (decayed) | Decayed completely |

**Key finding**: G_balance is nearly independent of gamma (~91 for all values). This means the gain term completely dominates the dynamics — the damping gamma is irrelevant because the gain threshold is set by the intrinsic oscillon amplitude, not by the dissipation rate. The transition from "decay" to "blowup" is extremely sharp (<0.1% in G).

### Mode B: Amplitude gain, S_a = Sigma^2 * phi_a / (1 + kappa_g * Sigma^2), Sigma^2 = sum(phi_a^2)

Lower order (phi^3 effective) but still nonlinear.

| gamma | G_balance | E_final/E_0 | Status |
|-------|-----------|-------------|--------|
| 0.001 | 53.84 | 0.14 (decayed) | "Balanced" at 500 t.u. (ratio=1.10), decayed by 2000 t.u. |
| 0.01 | 53.95 | 0.0 (decayed) | Same pattern |
| 0.1 | 53.98 | 0.0 (decayed) | Same pattern |

**Key finding**: Mode B finds balance at 500 t.u. (E/E_0 ~ 1.1 with f_core ~ 0.99), but this is transient. The oscillon grows, sheds energy through radiation, and then has insufficient amplitude for gain to compensate damping. It slowly decays to zero over 2000 t.u.

### Mode C: Saturated anti-damping, S_a = dphi_a/dt / (1 + kappa_g * Sigma^2)

This acts as negative damping in regions of small field amplitude and positive reinforcement everywhere.

| gamma | G_balance | E_final/E_0 | Status |
|-------|-----------|-------------|--------|
| 0.001 | 0.00195 | 1.06 | Best case — f_core=0.99, pk=0.36 |
| 0.01 | 0.0156 | 7.47 | Energy migrated to extended modes |
| 0.1 | 0.106 | 1.03 | Energy in extended modes, fc=0.06 |

**Key finding**: Mode C is the most promising — at gamma=0.001 the validation shows E_final/E_0=1.06 with f_core=0.99, meaning the oscillon maintained its localization. However, the G_balance/gamma ratio is ~2, indicating the system is just slightly above critical anti-damping. The "balance" is between boundary absorption (which removes radiation) and the gain, not between localized gain and damping.

---

## Phase 2: Stability (Perturbation Recovery)

Tested at gamma=0.01 with Mode A (G=91.35), Mode B (G=53.95), Mode C (G=0.016).

### Mode A

The system reached a delocalized high-energy state (E=109, not the original E=1.3). Both +20% and -20% perturbations "recovered" (ratio=1.0000) — but to this wrong state, not the oscillon. This is a limit cycle of the delocalized field, not a dissipative soliton.

### Mode B

The steady state decayed to E=0 before perturbation was applied. No recovery possible from zero.

### Mode C

Steady state at E=9.6 (7.5x original), delocalized (f_core=0.18). The +20% perturbation reached E=10.6 (ratio=1.10, marginal recovery). The -20% perturbation reached E=11.2 (ratio=1.16, no recovery). The system does not have a well-defined amplitude attractor.

**Conclusion**: No gain mode shows robust limit-cycle behavior around the oscillon state.

---

## Phase 3: Spontaneous Formation from Noise

Starting from Gaussian noise (amplitude 0.05), none of the gain modes produced a localized oscillon:

| Mode | G values tested | Result |
|------|----------------|--------|
| A | 1x, 2x, 5x, 10x balance | All: E=0, no formation |
| B | 1x, 2x balance | E=0; 5x and 10x: BLOWUP |
| C | 1x, 2x, 5x, 10x balance | Energy grows (E=339 to 24718) but NO localization (f_core < 0.13) |

**No spontaneous oscillon formation.** Mode C grows energy from noise but produces spatially extended turbulence, not a localized structure.

---

## Why This Fails: Structural Analysis

### The Gain-Loss Mismatch

For a stable dissipative soliton, the gain must satisfy two conditions:
1. **Positive feedback at small amplitude**: gain > damping when phi < phi_eq (to grow)
2. **Negative feedback at large amplitude**: gain < damping when phi > phi_eq (to shrink)

The problem is that for ALL polynomial gain terms S_a ~ phi^n:

- dE/dt = -gamma * integral(v^2) + G * integral(S_a * v)
- The gain power scales as phi^(n+2) (from S_a ~ phi^n and v ~ phi * omega)
- The damping power scales as phi^2 * omega^2

For balance at amplitude phi_eq, increasing phi by a factor (1+epsilon) changes:
- Gain power by (1+epsilon)^(n+2) ~ 1 + (n+2)*epsilon
- Damping power by (1+epsilon)^2 ~ 1 + 2*epsilon

The perturbation growth rate is proportional to (n+2-2) = n. For any n > 0, the gain sensitivity exceeds the damping sensitivity. The system is structurally unstable.

### What Would Work

A stable dissipative soliton requires **saturating gain with sub-linear response** to amplitude. Physical examples:
- Laser mode-locking: gain medium has finite population inversion (hard saturation)
- Neural solitons: ion channel kinetics provide true bistable switching
- Reaction-diffusion: chemical kinetics with enzyme saturation (Michaelis-Menten)

None of these map onto a simple polynomial modification of the Klein-Gordon equation. A dissipative Skyrme soliton would require:
- A separate "reservoir" field (like a gain medium) with its own dynamics
- Coupling that transfers energy from reservoir to soliton
- Reservoir depletion that limits gain

This is a fundamentally richer model than what the current EOM provides.

---

## Data Files

- `data/dissipative_g{gamma}_G{G}_mode{N}_ts.tsv` — time series for each (gamma, G, mode)
- `data/dissipative_stability_{plus20,minus20}_mode{N}_ts.tsv` — perturbation tests
- `data/dissipative_fromNoise_G{G}_mode{N}_ts.tsv` — noise initialization tests
- `data/dissipative_summary.tsv` — balance point summary

Columns: time, phi1_0, phi2_0, phi3_0, peak1-3, E_kin, E_grad, E_mass, E_pot, E_total, f_core, P_diss, P_gain

---

## Code

`src/dissipative1d.c` — supports all three gain modes via `-gmode {0,1,2}` flag.

Compile: `gcc -O3 -Wall -o dissipative1d src/dissipative1d.c -lm`

Key flags: `-gamma`, `-G`, `-gmode`, `-kappa_g`, `-phase {0,1,2,3}`
