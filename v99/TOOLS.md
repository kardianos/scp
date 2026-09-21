# v99 — computer algebra, when to reach for it

| tool | used for | use next only if |
|---|---|---|
| **numpy** | Gaussian `S_rel` (`phase1/cas/srel.py`) | covariance / linear algebra on a written finite system |
| **Maxima 5.49** | `E_w`, `ε`, argmins (`srel.mac`) | ordinary calculus, small explicit sums |
| **FriCAS** | **not used** | exact univariate polynomials |
| **Cadabra2** | **not used** | a variational identity on a *forced* 3+2 action. Writing an entropic action so Cadabra can vary it is Bianconi’s move and is refused (`PRINCIPLES.md` §5). |

NetworkX and Go (`v98/size`) were not opened. No new meter.

Reproduce only what was run:

```
python3 v99/phase1/cas/srel.py
maxima --very-quiet -b v99/phase1/cas/srel.mac < /dev/null
```
