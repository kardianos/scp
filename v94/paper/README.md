# v94 theory paper

- `scp_theory.md` — the paper source (pandoc markdown; session
  `aebbc7eb-4716-4e30-8c22-484e115428f2`).
- `fig_*.py` — one matplotlib script per figure; run each with `python3`
  to regenerate the corresponding `fig_*.pdf`.
- `scp_theory.pdf` / `scp_theory.tex` — built outputs (12 pp, two-column).

## Build

```sh
for f in fig_*.py; do python3 "$f"; done
pandoc scp_theory.md -o scp_theory.pdf --pdf-engine=xelatex \
  -V geometry:margin=1in --number-sections --toc
# .tex for inspection/journal submission:
pandoc scp_theory.md -o scp_theory.tex --pdf-engine=xelatex \
  -V geometry:margin=1in --number-sections --toc
```

Layout: the title block, abstract, and TOC are one-column; the body
switches to two columns via a raw `\twocolumn` at the start of the
markdown body (so do NOT pass `-V classoption=twocolumn` — that would
put the front matter in two columns too).

Notes: tables and figures are embedded as raw LaTeX (`table*`/`figure*`
for full-width floats) because pandoc's native pipe tables emit
`longtable`, which is incompatible with `twocolumn`. Math avoids a
closing `$` immediately followed by a digit (pandoc parses that as
literal text, not math).
