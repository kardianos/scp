# Bell rung B0 — how the two seats interact

THREE agents work this rung IN PARALLEL and INDEPENDENTLY, then cross-check.
(A third seat, GEOM, was added after FABLE and GROK started — if you are one of
those two, note it now and read its file.)

| seat  | model    | writes to      | extra brief    | reads |
|-------|----------|----------------|----------------|-------|
| FABLE | fable    | FABLE_WORK.md  | —              | BELL_BRIEF.md + all three work files |
| GROK  | grok-4.5 | GROK_WORK.md   | —              | BELL_BRIEF.md + all three work files |
| GEOM  | grok-4.5 | GEOM_WORK.md   | GEOM_BRIEF.md  | BELL_BRIEF.md, GEOM_BRIEF.md + all three work files |

GEOM pursues one specific hypothesis: that the metric quadratic form (the c²
relation) is the source of the 2√2 ceiling, while field monism supplies the
mechanism for exceeding 2. FABLE and GROK work the doors-and-costs angle
without that assumption. The two halves are meant to join if both succeed —
GEOM owns the BOUND, FABLE/GROK own the REACHABILITY.

RULES
1. Write ONLY your own file. Never edit the other seat's file or the brief.
2. Append as you go — do not wait until the end. The other seat reads you
   asynchronously, and a file that stays empty for three hours is useless to
   them.
3. Re-read the other seat's file periodically. When you disagree with it, say
   so IN YOUR OWN FILE, quote the passage, and give the reason. Disagreement
   that is argued is the point of running two seats; silent divergence is not.
4. Do not converge for the sake of converging. If you reach different doors or
   different bounds, that is a result — record both and state the discriminator
   that would settle it.
5. Divide labour opportunistically: if the other seat has clearly taken Door A,
   consider pushing Door B or the Tsirelson half harder, but do not leave your
   own independent judgement on the main question unstated.
6. Structure your file with: findings numbered; a claim-status tag on each
   ([D]erived / [M]easured / [P]ostulated / [C]onjecture); the symbolic
   worksheets inline or referenced by path; and a running VERDICT section you
   update as you go.
7. Maxima 5.46.0 and SymPy 1.12 are installed. Commit worksheets to
   v87/work/<seat>/ and reference them.
8. Budget: 4 hours wall clock. Land a defensible partial result rather than an
   undefended complete one. At the 3-hour mark, write your VERDICT even if the
   work is unfinished.
