# V41 Shell Structure Observations

## User Observations (from volview)

1. **Breakaway structures**: Both UUD and UDD show detached features away from the core.
   In UUD, the breakaway has mass (phi content). In UDD, the breakaway is almost
   entirely theta. Both are likely artifacts of the initial condition radiation, but noted.

2. **Bright blue and green shell**: The stable particle has a visible phi (green) and
   theta (blue) shell structure at the core boundary.

3. **Periodic blue spherical bands**: Periodic theta-dominated rings at larger radii,
   away from the core.

## Numerical Analysis

### Radial Profile Structure

**Core (r < 5):** Both UUD and UDD show identical profile shape:

| r | UUD phi_rms | UDD phi_rms | UUD θ/φ ratio | UDD θ/φ ratio |
|---|-------------|-------------|---------------|---------------|
| 0.25 | 1.154 | 1.252 | 0.120 | 0.114 |
| 2.25 | 1.112 | 1.139 | 0.080 | 0.062 |
| 4.25 | 0.835 | 0.769 | 0.087 | 0.115 |
| 6.25 | 0.420 | 0.343 | 0.151 | 0.170 |
| 8.25 | 0.129 | 0.119 | 0.216 | 0.227 |

θ/φ ratio INCREASES outward — theta becomes relatively more prominent at the
surface. This is the "blue shell" seen in the viewer: the core is phi-dominated
(green), the shell is where theta catches up (blue appears).

**The blue/green shell (r ≈ 4-7):** This is the boundary between the bound core
and the free outer region. Key features:
- φ_rms drops from ~1.0 to ~0.3 (70% decline in 3 code units)
- θ_rms drops more slowly (from 0.09 to 0.06)
- θ/φ ratio rises from ~0.08 to ~0.17
- This is where the "breathing shell" lives (from V40 structural analysis)
- Visually: green (phi) fading while blue (theta) remains → blue shell appearance

### The Periodic Blue Bands (Theta Standing Waves)

Both structures show a SECONDARY phi/theta feature at r ≈ 10-15 (UUD) and
r ≈ 12-16 (UDD). The data clearly shows this:

**UUD radial theta/phi ratio:**
```
r=8   θ/φ = 0.216  (declining from shell)
r=9   θ/φ = 0.228  (flat — entering background)
r=10  θ/φ = 0.230
r=11  θ/φ = 0.233  ← theta/phi ratio INCREASING
r=12  θ/φ = 0.184  ← dip
r=13  θ/φ = 0.168  ← SECONDARY PEAK in phi (phi bump)
r=14  θ/φ = 0.159  ← phi peak continues
r=15  θ/φ = 0.160
r=16  θ/φ = 0.161
r=17  θ/φ = 0.140  ← declining
```

There's a clear phi bump at r ≈ 12-16 (phi_rms rises from 0.072 to 0.106 then
falls back). This is NOT a monotonic decay — there's a SECONDARY CONCENTRATION
of phi field energy at ~3× the core radius.

**UDD shows the same pattern** but shifted and with more theta dominance in the
outer region:
```
r=9   θ/φ = 0.301  (higher than UUD!)
r=10  θ/φ = 0.223
r=11  θ/φ = 0.293  ← theta-dominated band
r=12  θ/φ = 0.157  ← phi recovery
r=13  θ/φ = 0.178
r=14  θ/φ = 0.129  ← another phi peak
r=22  θ/φ = 0.117  ← far-field phi dominated
```

### Breakaway Structure Characterization

**UUD breakaway structures (detected blocks > 3× background):**
- 5 **theta-dominated** blocks (θ > 3×bg but φ < 3×bg)
  - Located at r = 7.6-12.4 from centroid
  - θ/φ ratio ≈ 0.2-0.3 (θ relatively strong)
  - These are the **periodic blue bands** — theta standing waves that
    separated from the core during the initial radiation transient

- ~30 **phi+theta** or **phi-dominated** blocks
  - Located throughout r = 5-24
  - Most are the trailing edge of the core or reflected waves from the boundary

**UDD breakaway structures:**
- **10 theta-dominated blocks** — MORE than UUD (10 vs 5)
  - Located at r = 7.5-22.0 from centroid
  - Several at very large r (22.0), near the boundary
  - θ/φ ratio up to 0.5+ in some blocks (nearly pure theta)
  - This confirms the user observation: **UDD breakaway is mostly theta**

- ~30 phi-dominated or mixed blocks, similar to UUD

### The Key Difference: UDD Radiates More Theta

| Metric | UUD | UDD |
|--------|-----|-----|
| Theta-dominated breakaway blocks | 5 | **10** |
| Max distance of theta block from core | 12.4 | **22.0** |
| Background θ_rms | 0.00387 | **0.00286** |
| Core θ/φ ratio at r=9 | 0.228 | **0.301** |

UDD has MORE theta radiation to larger distances, but LOWER background theta.
This means the theta radiation in UDD is MORE CONCENTRATED in discrete structures
(breakaway blobs) rather than uniformly distributed.

**Physical interpretation:**

In UDD (one Up, two Down), the majority chirality is Down. The lone Up braid
creates a theta pattern that partially OPPOSES the two Down braids. This
opposition creates a theta DIPOLE — the theta field doesn't cancel uniformly
but instead forms standing wave patterns that can detach as discrete blobs.

In UUD (two Up, one Down), the majority chirality is Up. The theta field is
more uniform (less dipole character), so it radiates more smoothly and doesn't
form discrete breakaway structures as strongly.

**The "electron" interpretation (speculative):**

The user noted the breakaway resembles an electron-like structure. While this
is likely an artifact of the radiation transient, the physics is suggestive:
- The breakaway is theta-dominated (charge-carrying field)
- It's detached from the core (at a distance)
- In UDD, it has very little phi (almost no mass — just charge)
- The separation distance (r ≈ 10-15) is ~3× the core radius

In the Cosserat theory, a pure-theta blob would be a massless, charged excitation —
this IS what a photon or lepton would look like (no binding potential, just EM field).
However, the current breakaway structures are NOT self-sustaining — they're radiation
that will eventually disperse through the absorbing boundary.

A genuine lepton would need a mechanism to confine theta without phi — which the
current Lagrangian doesn't provide (theta is massless and uncoupled to V(P)).
This is consistent with the known physics: leptons don't participate in the
strong force.

### Summary of What The Viewer Shows

1. **Green core** (r < 4): Phi-dominated binding region, high |P|, the "quark matter"
2. **Blue-green shell** (r ≈ 4-7): Transition zone where θ/φ rises, the breathing surface
3. **Periodic blue bands** (r ≈ 10-15): Theta standing waves from initial transient
4. **Theta breakaway blobs** (r > 15): Discrete theta-dominated structures
   - More numerous and further out in UDD than UUD
   - Almost pure theta in UDD (charge without mass)
   - Mixed phi+theta in UUD (charge with some mass)
5. **Far field** (r > 20): Background oscillation, declining to boundary damping
