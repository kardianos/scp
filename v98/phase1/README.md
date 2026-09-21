# v98 Phase 1 — paper arms

Executed 2026-08-13. No kernel. CAS: Maxima 5.46.0, FriCAS 1.3.10.
Cadabra2 is installed and was not used (`../TOOLS.md`).

| arm | file | result |
|---|---|---|
| rule | `CLOSURE.md` | C1–C3 + named filters |
| **1A** | `CENSUS.md` | unique smallest exists; then a zoo |
| **1B** | `TROUGH.md` | continuum: no trough; discrete: trough slides with ρ |
| **1U** | `../UNFOLD.md` | numbers cited; no veto |
| **1C** | `DECIDE.md` | no new substrate; Phase 2 on, Phase 3 off |
| interpret | `INTERPRET.md` | how to read 1A+1B+1U as one argument |

Reproduce CAS and the census:

```
python3 v98/phase1/census.py
maxima --very-quiet -b v98/phase1/cas/scale.mac < /dev/null
maxima --very-quiet -b v98/phase1/cas/wind.mac  < /dev/null
fricas -nosman -noht < v98/phase1/cas/spectra.input
```
