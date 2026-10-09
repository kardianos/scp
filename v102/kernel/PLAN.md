# v102 kernel step — plan

**Authorized:** user, 2026-10-08, "(3) approved": put R8 (the 3+4
widening), R9 (T-floor) and R10 (S-boundary) into a real simulation.
`../README.md` row K records it.

**The ratchet governs.**

- Every kernel modification runs the full battery before commit.
- New behaviour is off by default and byte-inert at its default.
- Bars are never softened.

The C kernel of record is `v97/kernel/freecell.c`. The Go kernel
`v97/fab/` is not touched in this step: the new module is C-only, and
the battery's C ↔ Go cross-checks run at defaults, where the module is
off.

---

## 1. What the kernel must carry, and what it does not have

| needed (R1–R10) | the kernel today |
|---|---|
| field quanta (R2) in four modes per cell, spin ½ × isospin ½ (R8) | a classical two-plane real amplitude per cell (`fa1`, `fa2`, `Ee = fa1² + fa2²`); conversions fire in whole action atoms ε = A0·w/2π |
| a lock: a fiber of cells at their ceiling (R6, C-lock label), whose quantum is gapped and whose open slot is a (½, ½) nucleon | nothing: a classical field cannot host the Haldane-type lock (`../lock/`, `../perm/`) |
| conversion as a coherent pair amplitude g into a pitched link oscillator Ω (S1–S3), read virtually (V) | conversion as a rate law at a door (`f_conv`, atoms), into scalar space `Es` |
| T-floor transport inside the lock (R9) | classical hops on every link |
| S-boundary: space leaves only through a lock's boundary link (R10) | `s_k` space flow everywhere |

**So the kernel step is a quantum lock module inside a classical
kernel, coupled at the lock's boundary.** The rules at that boundary
are new. They have to be pinned and validated first (K1), then built
(K2).

## 2. Stages

| stage | what | the gate to pass before the next stage |
|---|---|---|
| **K1** | **reference simulator** (Python, exact plus TeNPy): the event rules at the lock boundary and the energy ledger, validated against the exact results already measured | each rule reproduces its exact counterpart within registered bars (`K1_PREREG.md`) |
| **K2a** | C module `lockq` (off by default): an exact lock state in a dynamic charge sector (small locks: 3+2 up to about 10 cells, W up to 6), Krylov time stepping, the K1 boundary events, the ledger | battery at defaults byte-identical; with `lockq` on, A/B against K1 on shared cases |
| **K2b** | the kernel field widened to four planes (`fa3`, `fa4` = the isospin partner), empty by default | battery at defaults byte-identical |
| **K2c** | an MPS lock (TEBD, LAPACK) for long fibers, with an explicit ledger term for truncation | A/B against TeNPy |
| **K3** | the demonstration in the kernel: light in the foam meets a seeded lock; nucleons form, walk, and pair | its own pre-registration |

## 3. The questions K1 must settle

1. **Formation (R10).**
   - Is an S-boundary event (a photon converting with a lock quantum,
     the space quantum leaving) captured by a quantum jump on the lock?
     The jump would apply Σ_μ P[μ, ν] b_iμ for an incoming mode ν.
   - Does its rate match the exact S-matrix (`../perm/` X4,
     `../fabric/` F2a)?
   - What lock state does the jump leave: what fraction is the dressed
     nucleon band, and what fraction nucleon + lock quanta?
2. **The energy ledger.** The space quantum that leaves must carry the
   exact energy balance: photon + lock before = lock after + space. The
   classical kernel conserves to 1e-15; the module must state its own
   error.
3. **Binding needs a sink.**
   - Two holes bind by −0.5 J0 (W), but the lock's gap is 1.32 J0. A
     closed fiber cannot shed the binding energy into a lock quantum.
   - So two holes meeting in an isolated lock should scatter, not bind,
     unless something carries the energy away.
   - What carries it: inelastic scattering of cold light at the
     boundary (X1), or nothing in the S-law? If light does, at what rate
     per photon?
4. **Mapping the classical field to quanta.**
   - The kernel's field is classical, real, and two-plane. What a
     classical neighbour amplitude means for an S-boundary event must be
     stated: an atom count, a polarization (S^z), and under R8 an
     isospin (I^z).
   - This is the rule K2b widens.

## 4. What this plan does not do

- It does not change the law table's constants.
- It does not touch `v97/fab` or `sfa.h`.
- It does not run the kernel with the module on before K1's gates
  pass.
