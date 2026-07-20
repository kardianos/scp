# v78 Run Protocol

Same append-only rules as v76/v77:

- Owner writes only own log under `logs/` and work under `work/<ID>/`
- `FOR_X` for cross-ideas; checkin each round
- O writes only `logs/O_orchestrator.log`

Entry format:

```
===== ENTRY <stamp> | <ID>-00N | <type> =====
Agent: <ID>
Round: <n>
Tags: ...
---
body
FOR_X: ...
===== END =====
```

**Kernel policy:** no scp_sim/sfa.h changes without human OK.
