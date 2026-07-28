# EMF — the field mode against real electromagnetism (living lab notebook)

**Living document.** Need and design for the field sector: what real
EMF phenomenology demands, what the two-plane amplitude already does,
what it measurably fails, and the design candidates. Append, date,
correct — never rewrite history.

Subordinate to `PRINCIPLE.md`; ratchet rule governs. The field mode is
one of the three established modes (space / dense / field); "EMF" here
means the field mode's correspondence with light — charge and the full
electromagnetic duality are NOT yet constructed and nothing below may
smuggle them in.

---

## 1. What already corresponds (measured, gated)

* Superposition, diffraction, two-slit fringes at parameter-free loci;
  complementarity dial; eraser ±; delayed choice (d1/t1/q2).
* Quanta: click grains on the ε(ω) grid; sub-threshold refusal —
  photoelectric structure (qt_lo/qt_hi, LIN).
* Exchange statistics ordering at a coupler (t4 HOM, partial depth).
* Momentum of light: ballistic energy centroid, momentum current
  Im[ψ*∇ψ] forward (p1); occultation by matter (g3).
* Vacuum linearity by construction (pitch load = bound energy only);
  Kerr nonlinearity is a property of loaded matter.

## 2. Need — the measured gaps (each is an experiment, not a shame)

| gap | reality | kernel today | status |
|---|---|---|---|
| **dispersion** | massless: ω = c·k, group speed k-independent | quadratic (electron-like); v_g ∝ field_J and k | OPEN — the deepest one |
| **radiation pressure** | light pushes absorbers/mirrors | ~100× recoil deficit at the conversion door (P2) | recorded; S2-full criterion |
| **polarization optics** | Malus's law, rotators, full vector optics | chirality pair (fb) exists, used for tagging/eraser only | analyzer machinery exists — cheap to attempt |
| **1/r² intensity** | point source flux law | untested (radial field-current instrument now exists) | cheap NOW |
| **Doppler** | moving source shifts frequency | untested (tilted-blob emitters translate) | cheap NOW |
| **c unification** | one ceiling for everything | dense flight has imposed C; field v_g emergent, k-dependent | tied to dispersion |
| **visibility ceiling** | fringe V ≈ 1 | V ≈ 0.46 (foam scatter + bandwidth) | engineering + foam regularity |
| **HOM depth** | full dip | g_b ≈ 0.36 (mode-match limited) | needs adiabatic coupler apparatus |

## 3. Design candidates

* **Linear dispersion.** The two-plane amplitude is second-order-like in
  its hop structure (quadratic band bottom). A massless branch wants a
  first-order (Dirac-like) coupling: the chiral PAIR (u±, CONSTRUCTION
  §2.2) hopped with an antisymmetric inter-component term could split
  the band into linear ± branches. This is a LAW change (field hop
  structure) — design on paper first, then a variant law table, full
  battery + a new dispersion experiment as its gate. If it works, c
  unification (one conversion ceiling) may fall out with it.
* **Radiation pressure** belongs to the dense side (translation IS the
  current — S2-full), not to field design; tracked in MASS/S2, cross-
  referenced here because the phenomenon is electromagnetic.
* **Polarization**: treat (fa, fb) as the two field species through real
  analyzer chains; Malus curve as the bar. No new law expected — an
  apparatus + measurement campaign.
* **Foam regularity as a physical parameter** for the visibility
  ceiling: scattering mean free path vs rjit — measures whether vacuum
  disorder is a tunable of the theory or a fixed prediction.

## 4. Experiment ladder (cheap first)

| id | question | needs |
|---|---|---|
| EM1 | ω(k): group velocity vs carrier k, both polarizations | existing pulse apparatus, k sweep |
| EM2 | 1/r² from a point emitter (radial field current vs shell) | rad_diag machinery on a field source |
| EM3 | Doppler: emitted frequency from a translating source | tilted-blob emitter + spectral read at a probe cell |
| EM4 | Malus: transmitted power vs analyzer angle | analyzer chain apparatus |
| EM5 | dispersion repair variant (linear branch) | LAW design + battery |
| EM6 | scattering mean free path vs rjit (foam as parameter) | long boxes, ensemble |

EM1–EM4 are apparatus-only and can ratchet into the battery as they
pass; EM5 is a law-table event.

## 5. Log

* **2026-07-28** — document opened. Standing table laws_V2g (21/21).
  Known deviations catalogued from ROADMAP §1 + P/G campaigns. Next
  cheap wins: EM1 dispersion scan and EM2 1/r² (instruments already in
  the kernel). Sequencing: after the speedup ladder and the MASS
  campaign start, per user direction.

* **2026-07-28 (standing practice)** — every EMF experiment renders
  frames alongside analytics: `viz/render_slice.py` (Ee panel + series
  mode). Visuals have caught what numbers missed; the requirement is in
  the MASS tech tree (D1) and the task tree.
