# v76 Multi-Agent Run Protocol

**Date**: 2026-07-18  
**Status**: active run setup  
**Problem**: [`PROBLEM.md`](PROBLEM.md)  
**Approaches**: [`APPROACHES.md`](APPROACHES.md)

---

## 0. Why this ops shape (and what we rejected)

**Chosen:** four long-lived approach tracks, each with:

- an **append-only log** under `logs/`,
- a **private work folder** under `work/`,
- one agent that **owns** that track.

Agents **read** all logs; they **write only** their own log and work folder.
Cross-ideas are written into the *author’s* log as `FOR_A` / `FOR_B` / … notes;
recipients must notice them on periodic log check-ins.

**Rejected alternatives:**

| Alternative | Why not first |
|-------------|---------------|
| Four git worktrees | Harder to share live logs; merge noise; overkill for paper/Lean |
| Single shared agent | Collapses the four epistemic stances (forward/reverse × theory/numerics) |
| Central wiki only | Loses chronological “who found what when” and race-free ownership |
| Agents writing each other’s logs | Corrupts provenance; two authors per timeline |

**Orchestration:** the parent session (or a later re-spawn) re-reads all four
logs between rounds. Agents do not delete or rewrite history.

---

## 1. Directory map

```text
v76/
  PROBLEM.md
  APPROACHES.md
  RUN.md                 ← this file
  README.md
  logs/
    A_forward_theory.log    # Approach A only writes here
    B_forward_numeric.log   # Approach B only writes here
    C_reverse_theory.log    # Approach C only writes here
    D_reverse_numeric.log   # Approach D only writes here
  work/
    A/                   # A’s scratch, Lean, notes, drafts
    B/
    C/
    D/
```

---

## 2. Approach ownership

| ID | Log file | Work dir | Mandate (from APPROACHES.md) |
|----|----------|----------|------------------------------|
| **A** | `logs/A_forward_theory.log` | `work/A/` | Forward theoretical / relational: axioms of medium → mass & warp; locality-\(c\); Lean/formal notes OK |
| **B** | `logs/B_forward_numeric.log` | `work/B/` | Forward numerical: medium evolution, free-budget deficit, rays; **not** kernel-v3 Q-balls as proof |
| **C** | `logs/C_reverse_theory.log` | `work/C/` | Reverse theoretical: from local \(c\), \(E=mc^2\), lensing → necessary medium laws |
| **D** | `logs/D_reverse_numeric.log` | `work/D/` | Reverse numerical: invert/score media; dualist adversarial baselines |

**Primitives (all agents):** field = only continuum; energy = ledger; mass =
bound field / \(c^2\); \(c\) = free-field locality (const in our frame); warp =
constant local \(c\) around locks. See `APPROACHES.md` §1–2.

---

## 3. Append-only log rules (strict)

1. **Only the owner agent appends** to its log. Never edit earlier entries.
2. **Never** write to another approach’s log.
3. **Never** truncate, rewrite, or “clean up” a log.
4. Prefer shell append:

   ```bash
   cat >> /home/d/code/scp/v76/logs/<LOGFILE> << 'EOF'
   ...entry...
   EOF
   ```

5. Or open the file, copy existing content, write **full old content + new
   entry** only if you cannot append — still no deletion of prior text.
6. Other agents **may and should read** all four logs periodically.
7. Ideas for another approach go in **your** log with a `FOR_X` tag (see §4).
   You do not implement their mandate unless it is a tiny clarifying note;
   they adopt or reject it after reading.

### 3.1 Entry format

```text
===== ENTRY <ISO-8601 or local stamp> | <short-id> | <type> =====
Agent: <A|B|C|D>
Round: <integer>
Tags: <comma-separated>
---
<title or one-line summary>

<body: findings, proposals, proofs, failed attempts, open questions>

Cross-links: <optional paths under work/ or other log entry ids>
FOR_A: <optional idea for A>
FOR_B: <optional>
FOR_C: <optional>
FOR_D: <optional>
===== END =====
```

**Types:** `bootstrap` | `proposal` | `finding` | `attempt` | `lean` |
`numeric` | `checkin` | `blocker` | `handoff` | `meta`

**short-id:** e.g. `A-003`, `B-012` (owner letter + zero-padded sequence).

### 3.2 Periodic check-in (required)

At least once per work round (and before closing a session), each agent:

1. Reads the **other three** logs (tail or full if short).
2. Appends a `checkin` entry noting:
   - what changed in other logs since last check-in,
   - any `FOR_<self>` items adopted, deferred, or rejected (with reason),
   - whether another approach’s result **blocks**, **enables**, or **contradicts** self.
3. Continues own mandate; does not take over their work folder.

---

## 4. Work folder rules

- Write freely under `work/<ID>/` (notes, Lean, Python, figures, drafts).
- Do **not** write into `work/<other>/`.
- Do **not** modify `scp_sim`, `sfa.h`, or v66–v75 kernels unless the **human**
  explicitly authorizes (project policy). v76 media sandboxes, if any, live
  under `work/B` or `work/D` as **new** code only.
- Prefer small, killable experiments over large unscoped builds.
- Point log entries at concrete paths: `work/A/axioms_v0.md`, etc.

---

## 5. Agent mandate template (what each agent is told)

1. Read `PROBLEM.md`, `APPROACHES.md`, `RUN.md`, own log, then other logs.
2. Advance **only** your approach (A/B/C/D).
3. Append findings/proposals/attempts to **your** log only.
4. Put durable artifacts in **your** work folder.
5. On cross-ideas: tag `FOR_X` in your log; do not edit their log.
6. Perform a **checkin** after substantive work and before finishing the round.
7. Stay monist-eligible: no fixed-background Q-ball “proofs”; no bolt-on
   \(T\to G\) as monism; locality-\(c\) seed is shared.
8. When blocked, log `blocker` with what you need from human or other approach.
9. End-of-round: summary entry with next steps for the **next** agent round.

---

## 6. Round lifecycle

```text
Round N:
  [spawn or resume agents A,B,C,D in parallel]
      each: read logs → work → append → checkin → stop with summary
  [parent reads all four logs]
  [human or parent decides Round N+1 focus]
```

Agents are **not** infinite daemons. Each spawn is one **round** of work
(hours of thought/tooling, not unbounded). Re-spawn for later rounds so
check-ins accumulate in the append-only logs.

---

## 7. Safety / project constraints

- **No** modification of `sfa/sim/scp_sim.c`, `scp_sim.cu`, or `sfa/format/sfa.h`
  without explicit human authorization.
- Large GPU / remote scp-runner runs are **out of scope** for v76 monism
  until a medium design is eligible (`PROBLEM.md` §7).
- Lean, local Python, small C sandboxes under `work/` are fine.
- Do not claim dark-matter or black-hole “results” without a closed medium law.

---

## 8. How to re-start a round

From parent session:

1. Skim all four logs (or `tail` recent entries).
2. Spawn four agents with `v76/RUN.md` + approach-specific focus for this round.
3. On completion, read logs and optionally write a short parent note in
   `work/ORCHESTRATION.md` (append-only if used) — **not** into A–D logs.

Optional helper:

```bash
# show last entry headers
rg -n "^===== ENTRY" /home/d/code/scp/v76/logs/*.log
```

---

## 9. Success for the run system itself

The ops setup works if:

1. Provenance is clear (who claimed what, when).
2. Cross-pollination happens via `FOR_X` + checkins, not file fights.
3. Each approach stays epistemically distinct (forward vs reverse, theory vs numerics).
4. Artifacts under `work/` are reproducible from log pointers.
5. Nothing dualist-ineligible is smuggled in as “v76 solved.”

The **physics** success criteria remain `PROBLEM.md` §7 and `APPROACHES.md` §4.
