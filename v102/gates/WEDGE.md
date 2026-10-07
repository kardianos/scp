# G-wedge — decision

**Written:** 2026-10-07, after `cas/out/wedge.out` and
`cas/out/two.out`. Gate: `../PREREG.md` §2.
**Status:** this page closes G-wedge. It does not edit G-two.

---

## 1. W1: the conversion operator. **PASS**

The derivation, step by step:

1. A conversion turns field quanta into space. Space carries no
   transverse polarization (v89 §2; the kernel's `Es` is a scalar).
2. R1 makes every coupling of the field SU(2)-invariant. A single
   field operator is a doublet, and there is no invariant linear form
   (`wedge.out`: 0 under SU(2), and 0 even under rotations about the
   link alone). So no single quantum can convert into a scalar mode.
3. The lowest order is therefore a pair, `M_ab b_a b_b`. SU(2) leaves
   exactly one invariant pair form: `ε`, which is antisymmetric
   (`wedge.out`).
4. The quanta are bosons (R2), so on a single cell
   `ε_ab b_{i,a} b_{i,b} = 0` identically (`wedge.out`: 0.0). The two
   quanta must come from two cells. Conversion is therefore a link
   process that takes one quantum from each end:

   `A_ij = b_{i1} b_{j2} − b_{i2} b_{j1}`.

No matrix was chosen.

**R1's scope carries this result.** Under rotations about the link
alone, two pair forms survive: `δ` (reflection-even) and `ε`
(reflection-odd). If R1 covered transport only, W1 would fail, with
the `δ/ε` ratio left free. W1 passes because R1 is the accepted "a
U(2)-invariant coupling" (`../PREREG.md` §1, amendment).

**Consistent with the record, not used in the derivation:**

- v89 §3 says conversion goes through a link and through chiral
  planes. `ε` is reflection-odd: it is the handedness of the pair.
- v89 law 2: "a channel is a joint process of its two ends."
- v90 REALITY B5: exclusion is "blocked on the exchange registry
  (two-quantum amplitudes with a sign)." `A_ij` is a two-quantum
  amplitude with a sign.

The v97 kernel converts on a single cell and without chirality. The
principle governs (`../README.md`).

## 2. W2: the form. **PASS**

On two sites with `N` quanta each, for N = 1…4 (residuals ≤ 7e-15):

- `A†A = N²/2 − 2 S_i·S_j`
- `B†B = N(N+2)/2 + 2 S_i·S_j`, where
  `B = Σ_a b_{j,a}† b_{i,a}` is transport.

Both are exactly bilinear. At leading order there is no biquadratic
term, for any `N`.

## 3. W3: the sign. **PASS**

At second order, through the converted state:

`H_eff = −(g²/Δ) Σ A†A = (2g²/Δ) Σ S_i·S_j − (g²/Δ) n N²/2`.

This is antiferromagnetic if and only if `Δ > 0`. The premise that
fixes `Δ > 0`: the lock is locked *field* (v98 §1, R5), so the field
configuration has to be the lock's ground state, not its converted
partner.

**Control, `Δ < 0`.** The ground state is mostly space: ⟨space⟩ =
3.92 on 4 links, and 5.88 on 6. That is not a lock of field. So
`Δ < 0` does not offer a ferromagnetic alternative; it removes the
lock.

**Instrument bars:**

| bar | n = 4 | n = 6 | limit |
|---|---:|---:|---:|
| spectrum-ratio miss at g/Δ = 0.1 | 2.3% | 0.044% | 5% |
| ground state total spin | 0 | 0 | 0 |

`J_eff` matches `2g²/Δ` within 2% at g/Δ = 0.05 and within 7% at
0.1; the gap is the fourth-order correction.

## 4. W4: the net sign. **PASS, resting on a licensed identity**

Transport's own virtual exchange is
`−(t²/U)(B_ij†B_ij + B_ji†B_ji) = −(4t²/U) S_i·S_j + const`. It is
ferromagnetic, and also bilinear. If it acts, the net coupling is
`J_net = 2g²/Δ − 4t²/U`, whose sign is a ratio of constants. That is
the failure route PREREG named.

It does not act, for these reasons:

1. **Quanta hoard.** The energy per quantum falls by about g²/Δ for
   each quantum added (`two.out`, n = 8: −1.41, −2.42, −3.41, −4.41).
   The held load law pushes the same way. "Load flattens pitch"
   (v89 law 2) means adding to a loaded cell costs less.
2. **Without a ceiling there is no chain.** Hoarding would pull every
   quantum into one cell. v98 §1 licenses "exclusion = an occupancy
   ceiling of a lock," and v90 PAULI-0 measured admission at cap as
   exactly zero.
3. **At a ceiling `N_c`, every cell sits at `N_c`.** A hop between two
   full cells would create `N_c + 1`, which the ceiling refuses. In
   the allowed space the hop's matrix elements are zero, and so is
   its virtual exchange. v98 §1's "locked field cannot travel as
   itself" becomes exact.
4. Therefore `J_net = 2g²/Δ > 0`.

PREREG required "a held sentence of the record … sets the locked-link
transport amplitude to zero." The ceiling is that sentence, once
step 1 holds. It is a *licensed* identity (v98 §1: "identities this
licenses, if the campaign closes them") with a measurement in v90's
substrate. v102 does not close it. W4 uses only the ceiling's
**existence**, not its value. A reader who requires a closed identity
should read W4 as conditional. That reading leads to the same
consequence in `../DECIDE.md`.

## 5. Decision

**G-wedge passes.** W1–W3 pass outright. W4 passes on the existence
of the v98 §1 ceiling. No sign, form, or matrix was chosen.

**What remains continuous:** only the `(g/Δ)²` correction. Over
g/Δ = 0.05 → 0.4, the spectrum ratio drifts by 0.6% on the n = 6 ring
and by 12% on n = 4. It is not scanned, and v102 does not use it to
place anything.

## 6. What is now true, and what is not

**True:**

1. Under R1 as accepted, plus R2, the conversion into scalar space is
   uniquely the antisymmetric (wedge) pair across a link.
2. Its leading exchange is exactly `(2g²/Δ) S_i·S_j` for every `N`.
3. It is antiferromagnetic because the lock is field.
4. Transport's competing exchange is ferromagnetic, and it vanishes
   at an occupancy ceiling.

**Not true:**

- The ceiling's value was used or derived.
- `N = 2`, `S = 1`, or ξ ≈ 6 follows from this page.
- The kernel's on-cell conversion was shown wrong. It is not the
  principle, and that is all.
