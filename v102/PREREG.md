# v102 — pre-registration

**Written:** 2026-10-07. This page comes after the scoping in
`v99/reopen/OPTIONS.md` and before any derivation text or number in
v102. The expectations in §6 are stated so that a reader can judge
how blind this was.

---

## 1. Setting

- **Fiber.** A ring `C_n` of cells, each of degree 2.
- **Field (R1, R2).** Bosons `b_{i,a}` with `a ∈ {1, 2}`, one mode
  per transverse plane. On each cell, `N_i = Σ_a n_{i,a}` and
  `S_i = ½ b_i† σ b_i` (Schwinger). A cell holding `N` quanta has
  spin `S = N/2`.
- **Space.** A mode with no transverse polarization, as in v89 §2 and
  the kernel's scalar `Es`.
- **Transport (R1).** U(2)-invariant:
  `t Σ_a (b_{i,a}† b_{j,a} + h.c.)`.
- **Every coupling of the field (R1)** is SU(2)-invariant.

*Amendment, made before any derivation:* the first draft of this
section applied R1 to transport only. The accepted wording is "a
U(2)-invariant coupling," so the scope is restored here and in the
README.

## 2. Gate G-wedge (paper, with an exact-diagonalization instrument)

**Pass requires all four:**

- **W1.** The link's conversion operator follows from
  `v89/PRINCIPLE.md` §3, the standing constraints, R1, and R2, with
  no free matrix.
- **W2.** For every `N`, the leading-order exchange that operator
  generates is affine in `S_i·S_j`, with no biquadratic term.
  Instrument: an identity check on two-site Schwinger spaces for
  N = 1…4, with residual < 1e-10.
- **W3.** The sign of that exchange follows from a stated premise,
  not a choice, and is antiferromagnetic.
- **W4.** The net leading-order exchange, including the virtual
  exchange from transport itself, has a sign fixed without comparing
  constants the premises leave unspecified. W4 passes only if a held
  sentence of the record or a ruling sets the locked-link transport
  amplitude to zero. W4 fails if the sign depends on the ratio
  `(g²/Δ)/(t²/U)`.

**Instrument bars** (they test the derivation; they are not gates):

- On a microscopic ring (cells plus one space mode per link,
  locked-link transport zero), at `g/Δ = 0.1` the ratio of the two
  lowest excitation energies matches the derived Heisenberg ring's to
  within 5%.
- The ground state is a total singlet on even `n`.

**Fail:** any of W1–W4 needs a choice.

## 3. Gate G-two (paper)

**Pass:** `N = 2` per cell follows from the premises, without using
the two-atom credit as an occupation, and without a minimality
preference ("take the smallest") that the record does not already
hold as a principle.

**Fail:** `N` is a conserved label that the premises do not fix, or
fixes only by a choice. In that case, `DECIDE.md` names the single
sentence that would fix it.

**Instrument:** the energy per quantum against `N` on derived chains,
for `S = ½, 1, 3/2, 2`.

## 4. Order and advance

G-wedge is decided first, in `gates/WEDGE.md`. G-two is decided next,
in `gates/TWO.md`, and does not edit G-wedge. Both decisions feed
`DECIDE.md`.

If both pass, Phase 1 is named, not authorized. Its bars:

- one isolated line;
- ξ measured on the derived chain and not read as a cell count;
- no second line below the two-particle threshold;
- the slit recovered from the unlocked, linearized field.

If either gate fails, v102 stops. `DECIDE.md` names the missing
premise.

## 5. What does not count

- "The relative phase is its handedness, so the coupling is the
  wedge." That is a reading. W1 needs a derivation that leaves no
  alternative.
- Choosing `Δ > 0` because it gives the answer. W3 needs the premise
  that makes `Δ > 0`.
- The two-atom credit as two quanta.
- `n = 6` from ξ ≈ 6.03.

## 6. Expectations, registered (not gates)

- W1–W3 look reachable. W4 is the risk: transport's ferromagnetic
  exchange (`v99/reopen` E1) competes.
- G-two looks likely to fail, because `N` looks like a conserved
  label.
