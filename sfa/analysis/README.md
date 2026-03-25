# SFA Analysis Tools

Suite of tools for analyzing SFA simulation output files. All tools support
f16, f32, and f64 column dtypes.

## Build

```bash
cd sfa/analysis && make
```

## Tools

### analyze_sfa
Per-frame global stability metrics: energy decomposition, P_int, theta_rms,
aspect ratio, stability score. Outputs JSON.
```bash
./analyze_sfa input.sfa --json output.json
```

### spatial_analysis
Per-frame spatial structure: centroid, R_rms, aspect ratio, cluster count
via BFS on |P|. Quick overview of structural evolution.
```bash
./spatial_analysis input.sfa
```

### cluster_profile
Per-cluster radial profiles (30 bins) of: ρ, |P|, θ, |v|, |a|, energy
decomposition, velocity divergence, phase coherence. Full force computation
from the equation of motion. JSON output.
```bash
./cluster_profile input.sfa --frame 0 --json output.json
```

### stable_vs_unstable
Compares clusters between two frames (early vs late), matches by centroid
proximity, classifies as GREW (stable) or SHRANK (unstable). Statistical
comparison of all metrics between stable and unstable groups.
```bash
./stable_vs_unstable input.sfa --early 0 --late 40 --json output.json
```

### shell_analysis
Radial shell profiles (50 bins) of phi and theta field energy centered on the
|P|-weighted centroid. Detects breakaway structures above 3× background.
Useful for examining the shell structure and theta radiation patterns.
```bash
./shell_analysis input.sfa [frame_idx]
```

### modify_sfa
Read an SFA frame, apply modifications (add braids, oscillons, perturbations,
scale, flip chirality), write a new single-frame seed SFA. Used in the
evolutionary search loop.
```bash
./modify_sfa input.sfa --add-braid 5 0 0 90 --scale 1.5 -o modified.sfa
```

## Viewer Modes

The viewer (`sfa/viewer/volview`) now supports three visualization modes:

- **Key 4**: Field mode (default) — R=|P| binding, G=φ² fabric, B=θ² angle
- **Key 5**: Velocity mode — R=|v| magnitude, G=inward velocity, B=tangential velocity
- **Key 6**: Acceleration mode — R=|a| total force, G=V(P) binding force, B=curl coupling force

Velocity mode requires 12-column SFAs (with velocity fields).
Acceleration mode computes forces from the equation of motion on each frame.
