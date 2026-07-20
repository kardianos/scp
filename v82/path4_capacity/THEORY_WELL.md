# Approach A — capacity well (orbit balance)

## Effective force (2D sandbox units)

Opposite locks \(q=\pm1\), separation \(D\):

\[
F_{\mathrm{EM}}(D) \approx +\frac{1}{2\pi D}
\quad\text{(attract along joining line; continuum 2D)}
\]

Capacity / footprint overlap (pair proxy of depletion):

\[
F_{\mathrm{core}}(D) = -k_{\mathrm{core}}\,e^{-D^2/(2s^2)}
\quad\text{(repel; \(s=\) footprint scale)}
\]

Net along joining line (sign: **+ attract, − repel**):

\[
F_{\mathrm{net}}(D) = F_{\mathrm{EM}}(D) + F_{\mathrm{core}}(D).
\]

## Equilibrium \(r_*\)

Solve \(F_{\mathrm{net}}(r_*)=0\) with \(F_{\mathrm{net}}'(r_*)<0\) in the sense that displacement from \(r_*\) restores (in 1D reduced mass: \(F>0\) pulls together when \(D>r_*\)).

Approximate: \(1/(2\pi r_*) = k_{\mathrm{core}} e^{-r_*^2/(2s^2)}\).

Circular \(v_t \approx \sqrt{|F_{\mathrm{net}}|\,r_*/\mu}\) with \(\mu=m/2\) if a small residual force band; better: seed at \(r_*\) with \(v_t\) from EM-only circular estimate then adjust.

## Success

- Pinned: \(F_{\mathrm{net}}(D<r_*)<0\) (repel), \(F_{\mathrm{net}}(D>r_*)>0\) (attract).  
- Free: multi-rev ≥1, \(\mathrm{sepmin}>0.6\), not monotone expand.

## Harmonics (secondary)

Once \(r_*\) exists, radiation selects among nearby \(v_t\); unstable seeds disperse energy into the medium.
