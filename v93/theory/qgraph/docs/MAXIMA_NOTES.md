# Maxima notes — what was checked symbolically, and how to re-run it

Four scripts under `qgraph/maxima/`. All are non-interactive-safe (they avoid
`solve`/`integrate` calls with symbolic exponents, which make Maxima ask a
question and hang under `-b`). Maxima 5.46.0.

**Run all four:**

```bash
cd /home/d/code/scp/v93/theory/qgraph/maxima
for f in stefan detailed_balance shells monodispersity; do
  echo "=== $f ==="
  maxima --very-quiet -b $f.mac < /dev/null
done
```

The `< /dev/null` is required: without it a stray `ask` blocks forever instead
of failing.

---

## 1. `stefan.mac` — the Stefan exponent is a spectral dimension

| step | result |
|---|---|
| substitute `w = T*y` in `g(w) w nbar(w) dw` with `g ~ w^(ds−1)` | integrand becomes `T^(ds+1) * y^ds/(e^y − 1)` — the `T` power is **exact and geometry-free** |
| Bose sum `int y^d/(e^y−1) dy = Gamma(d+1) zeta(d+1)` | verified term-by-term for `d = 1..4` against `d!/n^(d+1)` |
| `C(ds)` table | `ds=1: pi^2/6`, `2: 2 zeta(3)`, **`3: pi^4/15 = 6.493939`**, `4: 24 zeta(5)`, `5: 8 pi^6/63` |
| conclusion | `j = (c/4) u ~ T^(ds+1)` ⇒ **`p_rad = ds + 1`**, so `p_rad = 4` ⟺ `ds = 3` |
| matter side | `x ~ T` requires Dulong–Petit (fixed, fully live mode count) |
| freeze-out | `Em ~ T^(dm+1)` ⇒ `p_eff = 4/(dm+1)`: `dm=1 → 2`, `dm=2 → 4/3`, `dm=3 → 1` |

Go cross-check: `TestStefanConstantIsPi4Over15` (`1e−9`),
`TestLowLoadRadianceExponentSoftens`.

---

## 2. `detailed_balance.mac` — emission and absorption are one law

| claim | check | result |
|---|---|---|
| CLAIM 1 | `fullratsimp(O(a→b) − I(b→a))` with one coupling `A`, Fermi `n`, Planck `N` | **identically `0`** |
| ratio identity | `[n_a(1−n_b)]/[n_b(1−n_a)] − exp(−w/T)` | `0` |
| ratio identity | `N/(1+N) − exp(−w/T)` | `0` |
| CLAIM 2 | roots of `k cap (xb^4 − x^4)` | `x* = xb`; `d(I−O)/dx = −4 k cap x^3 < 0` (restoring); `Phi* = k cap xb^4 > 0` |
| CLAIM 3 | `(I−O)/(k cap xb^4)` in `u = x/xb` | `−u^4 + g u^3 + 1`, i.e. `u^4 − g u^3 − 1 = 0` |
| CLAIM 3 roots | `realroots` at `g = 0,1,2,3` | `u* = 1, 1.380278, 2.106919, 3.035744` |

The `g = 0` root being *exactly* 1 is the box-fill theorem: a random-phase bath
can never hold matter above ambient.

Go cross-check: `TestDetailedBalanceIsExact` (`1e−12` over a `(T, mu, eps)`
grid), `TestRandomBathFixedPointIsTheBathItself`,
`TestCoherentBathHoldsMatterAboveAmbient` (agrees with the Maxima roots to
`1e−6`).

**Maxima gotcha recorded:** `fullratsimp` on the *ratio* form hangs; the
difference form and `ratsimp` on `ratio − exp(−w/T)` are both instant. The
script uses the fast forms.

---

## 3. `shells.mac` — closed shells are neutral, unclosed ones are charged

Ring `C_N`, `lambda_m = 2(1−cos t)`, `circ_m = sin t`, `t = 2 pi m/N`.

| identity | result |
|---|---|
| `lambda(t) − lambda(2 pi − t)` | **`0`** — every `m ≠ 0, N/2` is paired: eigenspaces are 2-dimensional |
| `lambda(t) − 4 sin^2(t/2)` | `0` (closed form) |
| `circ(t) + circ(2 pi − t)` | **`0`** — a **full doublet is neutral, identically, for every `m` and `N`** |
| `circ(0)`, `circ(pi)` | `0`, `0` — both singlets neutral; the zone-boundary mode is real/standing |
| half-filled doublet, `N=6` | `circ = sqrt(3)/2 = 0.866` — nonzero: **charge** |
| magic numbers `Z_j`, `N = 4..14` | `1, 3, 5, 7, ...` (singlet, then a ladder of counter-rotating doublets) |

Go cross-check: `TestFullShellIsNeutralHalfShellIsCharged` (full shells to
`1e−16`, `N = 5..12`), `TestZoneBoundaryModeIsNeutral`, `TestRingMagicNumbers`.

The Go side extends this beyond rings: `TestHypercubeShellCapacitiesAreBinomial`
(`C(d,k)` at `lambda = 2k`, `d = 2..5`) and `TestIcosahedronShellLadder`
(capacities `1,3,5,3` — an `s/p/d` ladder, magic numbers `1,4,9` = `n²`).

---

## 4. `monodispersity.mac` — one exact size is `Delta/T`

| step | result |
|---|---|
| M2: one-quantum rates about a closed shell | `R+ = exp(−(eps_j1 − mu)/T)`, `R− = exp(−(mu − eps_j)/T)`; balance at **mid-gap** |
| M2: stability window | exactly the shell gap `Delta_j = eps_(j+1) − eps_j` |
| M3: residence time at mid-gap | `tau = exp(Delta/2T)/2`; `Delta/T = 1,2,4,8,16 → 0.82, 1.36, 3.69, 27.3, 1490.5` |
| M1/M4: mass distribution on a real `C10` ladder (`omega_k = 2 sin(pi k/10)`, caps `1,2,2,2,2,1`, `mu = 1.7`) | `T=0.05`: `P(N=7) = 80%` — a **species**; `T=0.3`: broad; `T=1.5`: nearly flat — a **hierarchy** |

Go cross-check: `TestColdDistributionCollapsesOntoAMagicNumber`,
`TestSpeciesSharpnessIsTheGapOverTemperature`, `TestLargerGapGivesSharperSpecies`.

---

## 5. What was *not* done symbolically

- **No force law for `Q`.** No symbolic derivation of an interaction between
  two unclosed shells was attempted; §6.4 of `UNIFIED_LAW.md` records this as
  open, and nothing in these scripts should be read as supporting it.
- **No continuum limit of the graph door.** L5 is exact on `C_N` by Fourier;
  the general-graph statement is checked numerically (orthonormality of the
  uniform door, `TestUniformDoorIsDiagonalOnAnyGraph`) rather than proved.
- **No spectral dimension of a free-cell foam.** That is a measurement, not an
  identity, and it is the top bar in `SIM_PREDICTIONS.md`.
