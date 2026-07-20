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

## Static force zero vs circular orbit (do not confuse)

**Force zero** \(F_{\mathrm{net}}(r_0)=0\) only marks repel↔attract boundary.

**Circular orbit** for reduced mass \(\mu=m/2\) with angular momentum \(L=\mu r v_\theta\):

\[
F_{\mathrm{along}}(r_c)=\frac{\mu v_\theta^2}{r_c}\;>\;0
\]

(\(F_{\mathrm{along}}>0\) = still **net attract**, providing centripetal force).  
So \(r_c\) lies on the **attractive side** of \(r_0\), not at \(F_{\mathrm{net}}=0\).

Seeding \(v_t=\sqrt{F_{\mathrm{EM}} r/\mu}\) at \(F_{\mathrm{net}}=0\) is **wrong** (over-speed) — see `parallel/GAP_ANALYSIS.md` §A1.

**Stable circular:** minimum of \(V_{\mathrm{eff}}=V(r)+L^2/(2\mu r^2)\) with \(\mathrm{d}V/\mathrm{d}r=F_{\mathrm{along}}\).

## Success

- Pinned: \(F_{\mathrm{net}}(D<r_*)<0\) (repel), \(F_{\mathrm{net}}(D>r_*)>0\) (attract).  
- Free: multi-rev ≥1, \(\mathrm{sepmin}>0.6\), not monotone expand.

## Harmonics (secondary)

Once \(r_*\) exists, radiation selects among nearby \(v_t\); unstable seeds disperse energy into the medium.
