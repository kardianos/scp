# v98 — computer algebra, when to reach for it

Installed and already used, or available. **Do not open a
CAS because it is there.** Open it when a claim is an identity
that page-and-pencil (or Go) will get wrong.

| tool | already used for | use next only if |
|---|---|---|
| **Maxima 5.46** | 1B scale/Derrick (`phase1/cas/scale.mac`); discrete `E_w+ρn` (`wind.mac`) | ordinary calculus, real/complex solve, small explicit sums |
| **FriCAS 1.3.10** | exact Laplacian characteristic polynomials (`spectra.input`) | exact univariate polynomials, integer matrices, factorizations over `ℤ` |
| **Cadabra2** | **not used** | a *variational* or *tensor* identity on a written 3+2 action: Euler–Lagrange, Derrick scaling of a real Lagrangian density, index bookkeeping. Not for graph census, histograms, or forcing `ρ` by decorating a well. |

Cadabra2 does not close the 1B gap by itself. Forcing `ρ` is a
conversion-identity question (one energy, space vs pattern).
Writing a Lagrangian so Cadabra can vary it is the imported
`V(ρ)` under another name, unless that Lagrangian is already
forced by the principle.

NetworkX (Python) enumerated graphs. Go (`v98/size`) is the
size meter. Neither is a CAS.

Reproduce only what was run:

```
maxima --very-quiet -b v98/phase1/cas/scale.mac < /dev/null
maxima --very-quiet -b v98/phase1/cas/wind.mac  < /dev/null
fricas -nosman -noht < v98/phase1/cas/spectra.input
```
