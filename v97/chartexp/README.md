# chartexp — forked inner-chart probe (2026-08-13)

A **copy** of `v97/fab/` that can carry a 3+1 inner-cell chart book.
The live kernel is untouched. Go only. Not a principle-face.

```
chart=0          off (default; no header, no write-path change)
chart=1          project occupations; dynamics identical to off
chart=2          clean grain: door ε from weff; w2e/beat untouched
chart=3          replace w2e so hops and door both see the chart pitch
chart=4          mode door: ε=A0|ω0−w1|/2π, residual priced; beat untouched

chart_scope=star   each live cell + its neighbours (default)
chart_scope=motif  connected live-Em component; n>24 falls back to star

chart_acct=energy  unweighted Laplacian (default)
chart_acct=mass    edge weight √(Em_i Em_j)
chart_acct=area    edge weight = slot area
```

Build: `go build -o chartexp/chartrun ./chartexp/cmd/chartrun` from `v97/`.
Runs: `chartexp/run.sh`. Results: `chartexp/RESULTS.md`.
