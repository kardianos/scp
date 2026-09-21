# v98 Phase 1 — execution (2026-08-13)

Paper only. Maxima and FriCAS used where a claim could be
proved or refused. **1C is written:** `phase1/DECIDE.md`.
Interpretation (not a substitute for the proofs):
`phase1/INTERPRET.md`.

Full notes: `phase1/CENSUS.md`, `phase1/TROUGH.md`,
`UNFOLD.md` §9–§10.

---

## 1A — census

All simple graphs on `n ≤ 7` (NetworkX atlas, 1252 graphs).
Cycle-closed (C1–C3): 577.

- Unique smallest: **`K₃`**. Then a zoo (3, 11, 60, 502).
- Measured filter F-odd (v90 odd-ring frustration): unique
  smallest **`C₄`**. Still a zoo.
- F-odd + 3-connected: unique smallest **`K_{3,3}`**; only
  `K_{3,4}` at `n=7`. Nearest thing to a line. Breaks open
  at `n=8` (cube and `K_{4,4}` both enter).

FriCAS: exact Laplacian polynomials. `C₆` is
`t(t−1)²(t−3)²(t−4)` — the CHARTEXP book, as a check.

**Result:** closure makes sizes *discrete*. It does not make
them *unique*. Graph selection alone is not a species.

## 1B — trough

3+2 fiber = the two planes already transverse to a 3-D link.
`W ∈ ℤ`; no half-wrap (no half-proton).

Maxima, continuum:

- fiber `|∇|²` is `R`-invariant (Derrick in the plane);
- `|∇|²` + area shrinks to `R = 0`;
- log-IR + area has **no real critical point** (formal roots
  imaginary — a trough I expected is not there);
- 3-D unlocked lumps shrink as `E ~ R` (Face C, as a theorem).

Maxima, discrete `C_n` plus site cost `ρ n`:

- `ρ = 0`: no finite minimizer (hoard);
- `ρ > 0`: a trough exists and **slides** (`n* = 24 → 3`
  as `ρ` goes `0 → 2`);
- at `ρ = 1`, even `n` ties **4 and 6**.

**Result:** the only trough found needs a conversion ratio
`ρ` that this arm did not force. Writing `ρ` by hand is the
imported well under another name.

Joint observation, not a derivation: F-odd + F-3 forbids
`n = 4`, so the `ρ = 1` tie would pick **6**. That is the
cheapest way a unique smallest robust even lock currently
sits somewhere. It is not 1C.

## 1U — unfold

Numbers pinned to PDG 2024, CODATA 2022, CREMA, SK, Ji/Yang,
DIS, AME2020 (`UNFOLD.md` §9). Classes unchanged. **No veto.**

The ledger still requires, of any future lock: coarse = one
size, fine = field; p and n the same size; no half-baryon;
Higgs numbers stay on the elementary gap.

## 1C — decision (closed)

Ending A/B/C/D of the charter, as slogans for a *substrate*,
do not fire. Phase 1 succeeds as an investigation: discrete
zoo, no continuum trough, a discrete trough that slides with
unforced `ρ`, no unfold veto. **No new kernel. Phase 3 is
not authorized. Phase 2 (meter only) is.**

Full text: `phase1/DECIDE.md`.

## What is not done

- Forcing `ρ` (the remaining 1B gap).
- Any kernel edit.
- Phase 3.

## Reproduce

```
python3 v98/phase1/census.py
maxima --very-quiet -b v98/phase1/cas/scale.mac < /dev/null
maxima --very-quiet -b v98/phase1/cas/wind.mac  < /dev/null
fricas -nosman -noht < v98/phase1/cas/spectra.input
```
