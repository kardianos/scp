# Phase 1 — relative entropy of a wrap on `C_n`

**Status:** executed 2026-08-19. Paper + CAS + numpy.
Scripts: `cas/srel.mac` (Maxima), `cas/srel.py` (Gaussian `S_rel`).
Logs: `cas/out/srel_mac.out`, `cas/out/srel.out`, plus TSVs.

The charter forbids importing \(S = \delta A/4\). This note asks
whether one energy, the conversion atom `ε = ω/2π`, and a proper
door force a unique `n*`. They do not.

---

## 0. Verdict in four lines

1. Exact `S_rel` of a coherent `W = 1` wrap vs the Gaussian vacuum,
   restricted to a proper arc, is well-defined and finite. Globally
   it is not (two pures).
2. That quantity is **not** `E_w` under another name, and it is
   **not** the Casini-kernel flux (discrete `S/K ~ 0.29–0.35` on a
   half-cycle). The continuum Dorau identity is not a theorem on
   these `C_n`.
3. Every conservation identity we can write without an area law
   either hoards (`n*` at the cutoff), sits on the smallest cycle,
   or reduces to v98’s unforced `ρ = 1`.
4. A grain door with a mass gap does not see `n` at all. Unique
   size is not a 3-D-door question.

---

## 1. What was computed

Free real scalar on the cycle graph `C_n`. Covariance of the
Gaussian ground state (`⟨φφ⟩`, `⟨ππ⟩`); zero mode omitted when
`m = 0`. Wrap: complex `W = 1` of amplitude `A` per site (two
chiral planes `φ = A\cosθ`, `χ = A\sinθ`), static or on-shell
traveling. Relative entropy of two Gaussians with the same
covariance and different means, on a contiguous arc `A` of `ℓ`
sites:

\[
S^{\mathrm{rel}}
= \tfrac12\,\bar\phi_A^\top X_A^{-1}\bar\phi_A
+ \tfrac12\,\bar\chi_A^\top X_A^{-1}\bar\chi_A
+ \text{(π blocks if traveling)}.
\]

This *is* the expectation of the exact modular Hamiltonian of the
reduced vacuum. Comparing it to the Casini interval kernel
`β_j = 2π j(ℓ−j)/ℓ` is a test of the *continuum Rindler*
approximation, not of Araki–Uhlmann.

Doors: 1 site (3-D grain), 2 sites (one edge), half-cycle.
Amplitudes: `A = 1` (TROUGH), `A = 1/√n` (unit L2), one quantum
(`E_w = ω`). Masses: `m = 0, 0.1, 1`.

Maxima identities: `E_w = 4n\sin^2(π/n)`, `ω = 2\sin(π/n)`,
`ε = ω/2π` (`A₀ = 1` lattice units), integer argmins, discrete
`ΔE_w + ε` and `ΔE_w + 1`.

---

## 2. `S_rel` is not `E_w`, and not `2π ×` Casini flux

Half-cycle, `A = 1`, static, `m = 0.1` (`cas/out/d1_half_casini.tsv`):

| n | `S_rel` | Casini `K` | `S/K` | `E_w` |
|---:|---:|---:|---:|---:|
| 4 | 1.85 | 6.28 | 0.294 | 8 |
| 6 | 2.42 | 8.38 | 0.289 | 6 |
| 12 | 3.06 | 9.82 | 0.312 | 3.22 |
| 24 | 3.61 | 10.21 | 0.354 | 1.64 |

`E_w` falls; `S_rel` on the half-cycle *rises*. `S/K` is not 1 and
is not 1/2π. Traveling wraps sit about 2× the static value at large
`n` (`d6_travel.tsv`) — both quadratures, not a 2π.

1-site door, `A = 1`, `m = 0.1`: `S_1` goes 0.269 (`n=3`) → 0.439
(`n=6`) → 0.693 (`n=24`). Slowly up, no interior peak.

---

## 3. Conservation identities, no area law

`ε(n) = \sin(π/n)/π`. Then `nε → 1` from below, `E_w/ε = 4π n\sin(π/n)`
increases `32.6 → 4π² ≈ 39.5`. Maxima I3–I7; numpy S1.

| total | all-n `n*` | even `n*` | reading |
|---|---|---|---|
| `E_w + n` (`ρ = 1`) | 5 | **4 and 6, tied** | v98 1B, unforced |
| `E_w + nε` | 24 (cutoff) | 24 | hoard: space is too cheap |
| `E_w/ε + n` | 3 | 4 | sits on the grain cycle |
| `S_1 + n` | 3 | 4 | space leg dominates |
| `S_{\mathrm{half}} + n` | 3 | 4 | same |
| `S_1` or `S_{\mathrm{half}}` alone | 3 | 4 | no space cost, smallest wins |
| `2π E_w + n` | 15 | 16 | imported 2π well |

Discrete walk, pay `ε(n)` to add a site (`s4_discrete_dE.tsv`):
`ΔE_w + ε` is **negative** on `n = 3..23`. Spreading always pays.
Pay 1 per site: first positive `ΔE_w+1` at `n = 5` — v98 `ρ = 1`
again, not a new law.

One-quantum wraps (`E_w = ω`): 1-site `S_rel` grows 0.052 → 0.119,
no trough. Unit-L2: 1-site `S_rel` *falls* `0.090 → 0.029` (the
wrap becomes invisible on a grain as it spreads). That is Face C,
in information language.

---

## 4. The grain door does not see width

1-site `S_rel`, `A = 1` (`d5_mass.tsv`):

| n | `m = 0` | `m = 0.1` | `m = 1` |
|---:|---:|---:|---:|
| 3 | 2.60 | 0.269 | 1.50 |
| 6 | 1.64 | 0.440 | 1.55 |
| 24 | 0.951 | 0.693 | 1.56 |

At `m = 1`, `S_1` is flat past `n = 4`. The 3-D door, with a mass
gap, reports a **fixed** distinguishability of amplitude. It cannot
pinch the fiber. `m` itself is unforced — another `ρ`.

---

## 5. Named readings that are not 1C

`ω_1(n) = 2\sin(π/n)`:

| match | `n` |
|---|---|
| `ω = 1` (grain; also `λ_1 = 1` on `C_6`) | **6** |
| `ω = T_{\rm Unruh} = 1/(2π)` at `κ = 1` | ≈ 39 |
| `ω = 2π` | off the bandwidth (`ω ≤ 2`) |
| `ω = 2` (band edge) | 3 |

`n = 6` from `ω = 1` is the CHARTEXP book of the painted vacuum
ring, recovered as a unit choice. Forbidden as a derivation.

---

## 6. What this does to 1C

- Do not open a kernel face. `S_rel` is not a dynamics.
- Do not replace `E_w` with `S_rel` in `E = αE_w + βn` and call
  `ρ` forced. The space leg is still the unforced coefficient.
- Do not import 2π, `A/4`, or `κ = 1` to move the argmin.
- The representation split survives: a grain door *cannot* select
  fiber width, which is why 3-D condensation could not. That is
  a reason the split is right, not a size law.

Reproduce:

```
python3 v99/phase1/cas/srel.py
maxima --very-quiet -b v99/phase1/cas/srel.mac < /dev/null
```
