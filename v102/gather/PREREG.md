# v102 gather step — pre-registration (how the background gathers)

**Written:** 2026-10-10, before any run of this step.
**Authority:** user, 2026-10-10: "keep going, report when it all
finishes. Try all 7 combinations." This follows `../vacancy/DECIDE.md`
("Then"), which named the open link of field → particle: **the
background is placed.** The standing rule: every decision is made from
experiments over all combinations of its candidates. New code is C;
LAPACK and BLAS do the linear algebra.

**A measurement in models. It rules nothing.** No kernel change, no law
change, no `v97/fab` change. Adopting any candidate needs a ruling.
Constants as in every v102 step:
- the S-law: E_f(1) = 1.46217, E_f(2) = 2.67252, Ω = 2.9;
- g/Δ = 0.1, so J0 = 2g²/Δ = 0.0095862;
- the T-floor t_in = 0.1 J0 (R9);
- R12 as the singlet core ('score', as in `../vacancy/`).

Two further readings come from the table, as used before:
- light transport t_L = 0.9 = field_J/2 (`../options/` E2 and E3);
- space transport s_k = 0.06, so the space band is |δ| < 2 s_k = 0.12
  (`../kernel/K2a_PREREG.md`).

## 0. The question and the candidates

The vacancy step found the particle is a vacancy in a background of
converted space, but that background was always placed: a set of cells
labelled "lock". How does it gather from field and space?

| | candidate | as a rule here |
|---|---|---|
| **a** | the reverse boundary event | traveling space arriving at a **boundary link** (one end in the background, one not) converts there into one quantum on each end. The outer quantum is light and may leave. This is the time-reverse of formation (R10), which is always on. |
| **b** | conversion in an unlabelled fiber | conversion acts on links with **no end in the background**. Two things follow. (1) Pairs there feel the leading-order exchange, including R12's core on (1,1). (2) Traveling space may arrive there and convert into a pair, and a pair there may convert into space that leaves. The released energy leaves as light: the light quanta carry it out through the lead. |
| **c** | membership by saturation | no label. A cell belongs to the background when it holds one quantum per link end and every neighbour holds field: N_i = z_i and N_j ≥ 1 for each chain neighbour j. Here z_i is the cell's number of links. |

With c off, membership is a fixed label Λ.

The rule sets are none, a, b, c, ab, ac, bc and abc: eight in all,
including the baseline.

## 1. The model: the leading-order fiber with events

This is the K-kernel's own structure (`../kernel/K2a_PREREG.md`): an
exact leading-order (reading V) vector, with formation and its kin as
golden-rule jumps onto bands. The kernel's lock is replaced by an
explicit 1-D chain whose membership follows the rule set.

**The chain.** It has n cells, i = 0 … n−1, joined by links (i, i+1).
- Cell 0 ends at a wall: it has one link, z_0 = 1.
- Cell n−1 couples to a semi-infinite light lead, so it has two links.
- Every other cell has two links: z_i = 2.

**Cells.** Each cell holds bosonic field quanta in M modes: M = 2 (3+2)
or M = 4 (W: spin × isospin, P = ε ⊗ ε). The occupation N_i is at most
its cap. On-site energy is E_f(N_i). There are no link oscillators
(reading V: space is never resident on a link; it arrives or leaves
as a traveling quantum).

**Membership μ_i.**
- c off: μ_i = 1 for i ∈ Λ = {0, …, k−1}, the placed region; the rest
  of the chain is light.
- c on: μ_i = [N_i = z_i and N_j ≥ 1 for every chain neighbour j]. The
  lead counts as an empty neighbour of cell n−1, so cell n−1 is never a
  member.

**Caps.**
- c off: cell 0 has cap 1, cells 1 … k−2 have cap 2, and cell k−1 has
  cap 1 (C-lock: a cell's cap counts its lock links). Light cells have
  cap 2.
- c on: cap_i = z_i, so 1 for cell 0 and 2 for the rest.

**Transport.** A hop of mode μ on link (i, j) takes configuration x to
configuration y.
- c off:
  - both ends in Λ: t_in, but refused if either end has N = 0 in x or
    in y (the T-floor may not empty or fill-from-empty a lock cell);
  - neither end in Λ: t_L;
  - one end in Λ: refused.
- c on:
  - if either end is a member in x or in y: t_in, but refused if a
    member end has N = 0 in the other configuration;
  - otherwise: t_L.

Caps are always enforced. Every rule is symmetric in (x, y), so H is
Hermitian.

**Exchange (leading order).** On link l:
- −(g²/Δ_type) A_l† A_l, with A_l = Σ_μν P_μν b_iμ b_jν;
- Δ_type is Δ on (2,2) and Ω − E_f(2) = 0.2275 on (2,1);
- on (1,1) links, the R12 core +(1.0/|P|²) A†A, as in `../vacancy/cas/vx.c`.

Which links carry it:

| | b off | b on |
|---|---|---|
| c off | links with both ends in Λ | every link |
| c on | links with at least one member end | every link |

Under c off the lock's own end links have both ends in Λ. Under c on,
the lock's light-facing end cell (N = 1 < z = 2) is not a member, so
"at least one member end" is what keeps the lock's last link.
Exchange keeps every N_i, so membership is unchanged by it.

**Link classes for events.** These are symmetric in (x, y).
- **member link:** both ends are members in x, or both are in y. Under
  c off: both ends in Λ.
- **boundary link:** not a member link, and some end is a member in x
  or in y.
- **non-member link:** otherwise.

**Events.** Golden-rule jumps between eigenstate clusters of H. A
cluster is a set of eigenvalues within 1e-9 of each other; a jump is
taken coherently within its target cluster.

| event | operator | when | rate × \|⟨k′\|O\|k⟩\|² |
|---|---|---|---|
| light out (lead) | b_{n−1,μ} | always | (1 + n_L)·√(4t_L² − ε²), ε = E_k − E_k′ − E_f(1), \|ε\| < 2t_L |
| light in (lead) | b†_{n−1,μ} | n_L > 0 | n_L·√(4t_L² − ε²), ε = E_k′ − E_k − E_f(1) |
| space out | A_l | boundary links: always (formation, R10); non-member links: b | (1 + n_s)·g²·√(4s_k² − δ²)/s_k², δ = E_k − E_k′ − Ω, \|δ\| < 2s_k |
| space in | A_l† | boundary links: a; non-member links: b | n_s·g²·√(4s_k² − δ²)/s_k², δ = E_k′ − E_k − Ω |

- Member links carry no event (R10: space leaves a lock only through a
  link with one end outside).
- The light rate is 2π t_L² ρ_end(ε), and the space rate is K2a's.
- The ledger closes by construction: every jump's energy goes to the
  quantum that leaves, or comes from the one that arrives, within its
  band.

**Environments.** n_s and n_L are flat across their bands.

| name | n_s | n_L |
|---|---|---|
| S | ν | 0 |
| L | 0 | ν |
| SL | ν | ν |

There are two regimes:
- **dilute,** ν = 1e-6: arrivals are rare against every spontaneous
  process faster than about 1e-8;
- **dense,** ν = 1: one quantum per mode.

The dark environment (ν = 0) is used for holding and light.

**Kinetics.** Kinetic Monte Carlo over clusters: competing rates, an
exponential clock, and a jump into the target cluster with the state
vector carried. Rates out of non-degenerate clusters are cached.

## 2. Systems

| fabric | n | Λ (c off) | light cells | largest sector |
|---|---|---|---|---|
| 3+2 | 6 | 0 … 3: the lock (1, 2, 2, 1) | 4, 5 | 1235 |
| W | 5 | 0 … 2: the lock (1, 2, 1) | 3, 4 | 4728 |
| 3+2, supplement | 7 | 0 … 4: (1, 2, 2, 2, 1) | 5, 6 | 6405 |

The supplement runs only if the main grid has finished.

**The lock pattern** is the configuration with Λ's cells at
(1, 2, …, 2, 1) and the light cells empty. Its weight in a state is
w_lock. The **lock state** |L0⟩ is the lowest eigenstate of the c-off,
b-off H in the sector N = 2k−2, 2S^z = 0 (2I^z = 0), restricted to the
lock pattern; it is an exact eigenstate there. **Fidelity** is
F = |⟨L0|ψ⟩|².

## 3. Experiments (every rule set; both fabrics)

| | name | start | environment | run | read |
|---|---|---|---|---|---|
| **X1** | gather | the empty chain | S, L, SL; dilute and dense | dilute: up to 100 arrivals; dense: T = 5000 | w_lock, F, member count, ⟨N_i⟩, quanta and energy in and out |
| **X2** | hold | \|L0⟩, projected on the rule set's clusters | dark; S dense; SL dense | dark: until no event can fire, or T = 1e9; dense: T = 5000 | w_lock and F at the end; survival = P(w_lock ≥ 0.9 at the end) |
| **X3** | light | two light quanta on the last two cells, modes 0 and M−1 (they can pair) | dark | until no event can fire | P(both leave as light); P(any space out); P(stuck) |
| **X4** | static | — | — | spectra | in the lock's sector: the lowest state's w_lock and F, its gap, its spontaneous rate out |
| **X5** | vacancy fate | the lock less one quantum. Two starts, each the lowest c-off, b-off state of an invariant subspace of sector N = 2k−3 (2S^z = 1, 2I^z = 1) with the light cells empty: **X5e**, the light-facing end cell empty (what formation makes in 1-D); **X5b**, one bulk cell at N = 1 (the mobile vacancy) | S dilute; S dense | as X1 | P(refilled: w_lock ≥ 0.9 at the end); time to refill |

**Trajectories:** 200 per run for 3+2, 100 for W. For dense runs, the
"end" reading is the time average over [T/2, T].

**X5 side variant (c off only: none, a, b, ab; start X5b).** In a 1-D
chain a bulk vacancy can never reach the boundary link: the T-floor may not
empty the end cell. In the kernel's 3-D lock, every bulk cell also has
boundary links to fabric cells around it. So X5 is repeated with a side
coupling on each Λ cell, i.
- **Side formation:** light from a side continuum (n_L), plus a quantum
  at i, converts into space that leaves.
- **Side a:** arriving space (n_s) converts into a quantum at i plus a
  light quantum that leaves to the side.
- Both are two-continuum golden rules:
  - rate = Σ_ν \|⟨k′\|J_iν^(†)\|k⟩\|² · K(ΔE), with J_iν = Σ_μ P_μν b_iμ;
  - K(Δ) = 2π g² ∫ dδ ρ_S(δ) ρ_L(Ω + δ − Δ − E_f(1)) (with its sign
    reversed for formation), times n or 1 + n.

Under c on, side links would count in z_i and are never converted, so
no cell could be a member. The side variant is therefore c-off only.
Side formation needs light, so an SL-dense side run is added. It is a
labelled supplement, not scored.

**Why 1-D for c.** In a 3-D fabric a cell has six links, and c's
literal reading then needs all six converted. `../options/` E4
measured that ceiling: spin-3 cells, with gap 0.010 J and ξ ≈ 640. The
1-D chain (z = 2) is c's best case.

## 4. Gates (T0, before any registered run)

1. **Hermiticity:** H − Hᵀ = 0 in every sector, under every rule set.
2. **Eigenpairs:** ‖HU − UE‖ ≤ 1e-9.
3. **LO against `../vacancy/cas/vx.c`:** with c off and b off, the lock
   sector's lowest energy, less Σ E_f, equals vx.c's LO open chain of k
   cells to 1e-10 (3+2 n = 4, W n = 3).
4. **Inertness of a and b on the placed lock (c off):** |L0⟩ is an
   eigenstate of H under none, a, b and ab with the same energy. This
   holds exactly under none and a. Under b, it holds if no exchange
   acts on the boundary link at (1, 0). That is A on an empty end, so
   it should hold.
5. **Bookkeeping:**
   - every jump's ε or δ lies inside its band;
   - N_chain + (light out − light in) + 2·(space out − space in) is
     constant along every trajectory;
   - the energy ledger E_chain + E_out − E_in is constant to 1e-9.
6. **Completeness and clusters:** for every event vector v = O|ψ⟩,
   Σ_k′ |⟨k′|v⟩|² = ‖v‖² to 1e-10, and each cluster's rate equals the
   sum of its eigenstates' rates.

## 5. Bars and verdicts

Per rule set and fabric:

| bar | passes when |
|---|---|
| **G1 gathers** | in at least one X1 environment and regime, the mean end w_lock is ≥ 0.5 |
| **G2 holds** | X2: dark survival ≥ 0.9, and the S-dense mean end w_lock ≥ 0.5 |
| **G3 same object** | among G1's X1 runs and X2's dark runs, mean F / w_lock ≥ 0.9 where w_lock ≥ 0.5 |
| **G4 light is free** | X3: P(both leave as light) ≥ 0.9 |
| **G6 no label** | c on |

A rule set **gathers the background** if G1, G2, G3 and G4 pass in both
fabrics.

**Endings:**

| ending | condition |
|---|---|
| **Gc** | a c-on set gathers the background: it gathers from field and space with no label |
| **Gl** | no c-on set does, but a c-off set does: it gathers only into a placed label |
| **H** | no set passes G1, but a c-on set passes G2–G4: a label-free background holds but does not gather |
| **N** | none of these |

X5 is reported, not scored: does a refill the vacancy, so that the
particle has a finite life in a space environment?

## 6. Disclosed before running

- **These are small chains.** The 3+2 lock has 4 cells and W's has 3.
  The light region is two cells plus a Markov lead. Every verdict is
  about these sizes.
- **The environment occupation ν is not a reading.** Two regimes are
  registered (1e-6 and 1). Nothing is scanned to place an outcome.
- **The secular jump picture** treats each cluster's internal state
  coherently, and different energies as dephased. K1 validated this
  for formation (fidelity 1.0000 against exact absorption). It is
  carried here for all events.
- **Two tool fixes, made during T0 and before any registered run.**
  - A transition within 1e-9 of a band edge has zero density of
    states. The two-cell light region has exact band-edge
    coincidences, which otherwise give round-off rates of about 1e-8;
    gate 6 caught them.
  - Rates depend on the environment, so the per-eigenstate rate cache
    is cleared at every environment change. A smoke test caught this:
    the first run's empty-light rates were reused under L.
- **What I expect, written before running.**
  - Under c off, an empty Λ cannot fill from outside. Its interior is
    reachable only by member-link events (R10 forbids them) or by
    T-floor hops from empty cells (refused). So G1 should fail for
    every c-off set.
  - b's arrivals on an empty link make a (1,1) singlet, which R12's core
    lifts by +1.0, far out of the space band. So b should never start
    a background from bare fiber.
  - Under c on, the lock's light-facing end cell is not a member, and
    it hops at t_L. So the open lock should evaporate from that end
    (G2 fail) unless a refills it.

## 7. Reproduce

Everything runs from `cas/`:
- `gcc -O3 -march=native -fopenmp -o gx gx.c -l:liblapack.so.3 -l:libopenblas.so.0 -lm`;
- then `./run_t0.sh` and `./run_grid.sh`;
- tables come from `go run tab.go`.

Outputs go to `cas/out/`.
