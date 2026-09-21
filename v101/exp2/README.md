# v101 experiment 2 — octave

Pre-registered, then run after experiment 1's decision. No kernel.
The unison ring at twist \(2\pi/3\) is a control, not a score.

| | file |
|---|---|
| gates | `PREREG.md` |
| decision | `DECIDE.md` — G2A pass, G2C fail, **FAIL** |
| instrument | `cas/octave.py`, `cas/monodromy.mac` |
| logs | `cas/out/octave.out`, `cas/out/monodromy.out` |

```
python3 v101/exp2/cas/octave.py
maxima --very-quiet -b v101/exp2/cas/monodromy.mac < /dev/null
```
