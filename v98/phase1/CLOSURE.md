# Phase 1A — the closure rule

**Status:** written so the census can run. Not a kernel law.
A graph that fails this rule is not a lock *candidate*. Passing
it is not a proton.

Companion: `CENSUS.md`, `cas/`, `census.py`.

---

## 0. What “closed” means here

`v89/PRINCIPLE.md` §3: conversion proceeds through complete
cycles; a partial cycle converts nothing.

v98 reads that as a statement about **which contact graphs can
host a lock**, not as a metabolic I=O balance and not as
qgraph’s “sit at magic N on a graph you were given.”

The graph is the contact set *as it stands*. It is not a stage.

---

## 1. The rule

Let `G = (V, E)` be a finite undirected simple graph.

**G is cycle-closed** iff all three hold:

| id | clause | equivalent |
|---|---|---|
| **C1** | `G` is connected | one component |
| **C2** | `G` is 2-edge-connected | no bridges; every edge lies on a cycle |
| **C3** | cycle rank `b₁ = \|E\| − \|V\| + 1 ≥ 1` | at least one independent cycle |

C3 is implied by C1+C2 for simple graphs (the smallest is `C₃`).
It is written so a single loop or a lone vertex cannot sneak
through on a wording technicality.

**Reading.** An edge that is a bridge is a conversion that
cannot close. That edge is not a lock; it is a leak or a
lead. Trees (`b₁ = 0`) are the extreme: they can transport,
they cannot lock. That matches the measured vacuum ring
(needs a cycle) and the measured open packet (diffracts).

This is the **raw** rule. It is combinatorial. It does not
mention 3-D embedding, harmonics as law, or a trough.

---

## 2. Named filters (not the rule)

The census reports these as separate columns so a later 1C
can adopt at most one without pretending it was C1–C3.

| id | filter | why it is named |
|---|---|---|
| **F-cut** | 2-vertex-connected (no cut vertex) | lock survives losing one contact |
| **F-3** | 3-vertex-connected | lock survives losing two contacts |
| **F-odd** | no odd cycle (girth even, or bipartite) | v90: odd rings are π-frustrated; even rings live |
| **F-W** | some cycle admits winding `W ∈ ℤ \ {0}` | automatic from C2; listed so occupancy is not confused with topology |

F-odd is the only filter that is already a *measurement* on
this substrate. It is still a filter, not smuggled into C2.

---

## 3. Minimality

Two notions, both used:

- **Edge-minimal cycle-closed:** no proper spanning subgraph is
  cycle-closed. These are exactly the simple cycles `C_n`
  (`n ≥ 3`).
- **Smallest by vertices** (then by edges) in a class.

Unique size, if it lives here, would be: a class whose
smallest member is unique *and* whose next members are
discrete jumps, not a zoo of the same `n`.

---

## 4. What this rule refuses to do

- It does not score `Em*`.
- It does not treat Laplacian shells as wallpaper to paint
  onto hoards. Spectra are diagnostics (`cas/spectra.input`).
- It does not claim `C₆` is a proton because a ring held
  winding in vacuum.
- It does not import a U(1). Winding is a cycle integer.

---

## 5. Pre-registered 1A questions

1. Under C1–C3 only, is there a unique smallest simple graph?
2. Under C1–C3 plus F-odd (the measured filter), same question?
3. Under F-3, same question?
4. After the smallest, is the rest a **line** (one graph per
   `n`, or a short forced ladder) or a **zoo**?
5. If it is a zoo, graph closure alone does **not** force a
   monodisperse species. That is a 1A result, and it makes
   1B load-bearing.

Answers: `CENSUS.md`.
