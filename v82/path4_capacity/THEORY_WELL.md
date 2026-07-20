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

**Circular orbit** for reduced mass \(\mu=m/2\), relative speed \(v_{\mathrm{rel}}\):

\[
F_{\mathrm{along}}(r_c)=\frac{\mu v_{\mathrm{rel}}^2}{r_c}\;>\;0
\]

(\(F_{\mathrm{along}}>0\) = still **net attract**). So \(r_c\) is on the **attractive side** of \(r_0\).

Equal-mass COM seed: each lock gets \(v_{\mathrm{each}}=v_{\mathrm{rel}}/2\).  
Seeding \(v_{\mathrm{each}}=\sqrt{F_{\mathrm{EM}} r/\mu}\) at \(r_0\) is **wrong** (uses bare EM and treats each as reduced mass).

**Stable circular:** \(F'(r_c)>-3F(r_c)/r_c\) (min of \(V_{\mathrm{eff}}\)).

## Numeric confirmation (2026-07-20, `parallel/`)

- Continuum \(k=0.2\), \(s=4\): \(r_0\approx8.64\); analytic orbits **PASS** (flat sep, multi-rev).  
- Live PIC: **PASS** band at \(r_c\approx9.82\), revs\(=2\).  
- CONTROL wrong seed: expands.  
- Lattice \(F_{\mathrm{EM}}\) only trustworthy at integer \(D\) (CIC noise off-grid).

## Success

- Pinned / continuum: \(F_{\mathrm{net}}(D<r_0)<0\) (repel), \(F_{\mathrm{net}}(D>r_0)>0\) (attract).  
- Free analytic or live: multi-rev ≥1, band on sep, not monotone expand.

## Harmonics (secondary)

Once \(r_*\) exists, radiation selects among nearby \(v_t\); unstable seeds disperse energy into the medium.
