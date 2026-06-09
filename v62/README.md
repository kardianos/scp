# v62 — The Number-Type Map of the SCP Residuals

**Status**: spearhead complete. Python regression **4/4**; two Lean modules build
against the v59 Mathlib (`PhaseNoGo.lean` = 1 cited `sorry` for
Lindemann–Weierstrass, everything downstream proven; `EwVevHome.lean` = no
`sorry`).

**Thesis in one line**: every quantity in the SCP program is either *algebraic*
(reachable by representation theory) or *transcendental* (not), and that single
distinction reorganizes every open "residual conjecture" of v59–v61.

Start with [`THESIS.md`](THESIS.md). It supersedes
`speculative/OCTO_SPECTRAL_PROJECTION_MENTAL_MODEL.md` as the working map.

---

## Headline results

1. **The phase no-go is a theorem.** The Brannen phase invariant `cos(3φ) =
   cos(2/3)` and `cos(2/9)` are transcendental (Lindemann–Weierstrass), so the
   lepton phase `φ = 2/9` is **not producible by any algebraic construction** —
   no Casimir, character, dimension, or eigenphase of an algebraic matrix. This
   *closes* the spectral-projection program for the phase and shows the
   speculative "07 cut with better J's" is barred a priori. → [`no_go/`](no_go/)

2. **The flat-direction law.** Rep theory fixes the symmetric skeleton and is
   blind to the flat angular directions it leaves; the Brannen phase and the EW
   `S^783` democracy are the *same* phenomenon (magnitude fixed, angle flat). →
   [`flat_direction/`](flat_direction/)

3. **The residual partition.** v61's four "equal-footing" residuals split into
   three fates: an unavoidable **scale** (`a`), a **closable degeneracy**
   (democracy — selected three independent ways), and **transcendental couplings**
   (`φ`, `α`, `f_g`) that are barred from the algebra. → [`residual_audit/`](residual_audit/)

---

## File map

```
v62/
  THESIS.md                          the model (read this first)
  README.md                          this file
  verify_all.py                      regression harness (python 4/4)
  no_go/
    NOGO.md                          rigorous paper proof of the phase no-go
    transcendence_nogo.py            numerical verification + 01-06 explanation
  flat_direction/
    FLAT_DIRECTION.md                the law + its two instances
    flat_direction_demo.py           symbolic (phase) + numerical (democracy) proof
  residual_audit/
    RESIDUAL_PARTITION.md            the 3-type partition + input count
    democracy_selection.py           3 selectors all pick democracy (closable)
    phase_dynamical_home.py          the P3 loop-ratio = phase's only home
  lean/
    PhaseNoGo.lean                   machine-checked no-go (rewrites v59 LeptonPhaseEmpirical)
    EwVevHome.lean                   democracy/S^783 (carried from v61)
```

## Run

```bash
python3 v62/verify_all.py            # all python self-checks

# Lean (needs the v59 Mathlib lake env; ~1-2 min each to load Mathlib):
cd v59/furey_construction/lean && lake env lean ../../../v62/lean/PhaseNoGo.lean
cd v59/furey_construction/lean && lake env lean ../../../v62/lean/EwVevHome.lean
```

(If `lake` is not on `PATH`, use `~/.elan/bin/lake`.)

---

## Relationship to v59–v61

v62 does not introduce new physics; it **re-organizes** the existing results by
number-type and turns three of the program's standing remarks into proofs/demos:

- v59 `CLOSEOUT.md` N5 and `LeptonPhaseEmpirical.lean` *noted* `cos(2/3)` is
  transcendental → v62 *proves* it and draws the no-go (`PhaseNoGo.lean`).
- v61 `EwVevHome.lean` *established* the `S^783` democracy degeneracy → v62 shows
  it is the same flat-direction phenomenon as the phase, and that it is
  **closable** by state-selection (unlike the phase).
- v59 `algebra/v_eff_loop.py` *built* the loop `V_eff(φ)` → v62 shows the
  loop-coefficient ratio is the *only* number-type that can equal `cos(2/3)`,
  placing the phase firmly in P3.

The actionable output is an **effort allocation** (THESIS §8): stop forcing the
algebra to emit transcendentals; close the democracy degeneracy and the P4
soldering action; treat the phase as a P3/loop problem and `α` as an input.
