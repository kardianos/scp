# v100 Phase 0 — execution (2026-09-21)

Paper only. Newton on the hedgehog energy, and one derivation.
**Decision:** `phase0/DECIDE.md`. Numbers: `phase0/PROFILE.md`.
Spin: `phase0/SPIN.md`. Gates were locked first:
`phase0/PREREG.md`.

---

## Result

G1 **passes**. G2 **fails**. Phase 1 is not opened.

The corner term produces a single hedgehog minimum,
\(E = 72.923\) in units \(f_\pi/e\), with \(E_2 = E_4\).
Re-integrated half-radius and double-radius ratios are
1.25003 and 1.25001. The kill floor was 1.05.

The two-atom credit does not fix an odd spin level. Spin
\(1/2\) stays an import. The advance rule is a conjunction,
so the mesh is not written.

**User, same day:** any hedgehog is incompatible with this
theory, across the experiments. The phase-0 minimum is a
radial profile \(F(r)\). That is the object `v89/PRINCIPLE.md`
§5 retires, and the object whose index read 0 on every shell
in `v89/KIMI_REVIEW_2026-08-02.md`. The wall stays on record.
It does not become a mesh.

## Reproduce

```
python3 v100/phase0/cas/profile.py
maxima --very-quiet -b v100/phase0/cas/virial.mac < /dev/null
```
