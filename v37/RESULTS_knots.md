# Results: Knot Geometries (Trefoil + Figure-Eight)

## Concept

A single continuous tube following a knot curve, carrying all three fields
with phase offsets along the arc length. Advantages over crossed braids:
- One continuous tube — can't fragment into pieces
- Closed path — no endpoints to radiate from
- Self-crossings provide V(P) binding at multiple points

## Geometries

### Trefoil knot (2,3 torus knot)
```
x(t) = (R + r*cos(3t)) * cos(2t)
y(t) = (R + r*cos(3t)) * sin(2t)
z(t) = ±r * sin(3t)       [+ = right-handed, - = left-handed]
```
R=3.0, r=1.5. Three self-crossings. Two chiralities (mirror images).

### Figure-eight knot (4₁)
```
x(t) = scale * (2 + cos(2t)) * cos(3t)
y(t) = scale * (2 + cos(2t)) * sin(3t)
z(t) = scale * sin(4t)
```
scale=1.5. Four self-crossings. Amphichiral (same as its mirror).

## Parameters

- N=128, L=25 (larger box than v2), T=300
- m²=2.25, μ=-41.345, κ=50, η=0.5
- tube_r=2.0, A=0.8, n_osc=3 (trefoil) / 4 (figure8)
- Phase offsets: δ = {0, 3.0005, 4.4325}
- Absorbing BC: width=4.0, rate=0.005 (gentler than v2)
- Stationary start (zero braid velocity)
- Per-field directional background
- SFA output: f32 columns (lossy but compressed)

## Results

| Geometry  | E_pot₀ | E_pot_final | P_int₀ → final | Energy drift | Clusters | Survived |
|-----------|--------|-------------|----------------|--------------|----------|----------|
| trefoil_R | -45.5  | -0.003      | 76 → 1.6       | -78%         | 0        | NO       |
| trefoil_L | -45.5  | -0.005      | 76 → 2.2       | -78%         | 0        | NO       |
| figure8   | -51.2  | -0.002      | 76 → 1.3       | -78%         | 0        | NO       |

All three dissolved. Binding energy drained to zero by ~t=150-200.
Trefoil L ≈ trefoil R (mirror symmetry confirmed).
Figure-eight had slightly stronger initial binding but same outcome.

### Energy timeline (trefoil_R)
```
t=0:   E_total=5120  E_pot=-45.5  P_int=76   clust=1
t=25:  E_total=4620  E_pot=-0.3   P_int=9    clust=1
t=50:  E_total=3450  E_pot=-0.4   P_int=11   clust=1
t=100: E_total=2440  E_pot=-0.1   P_int=6    clust=1
t=200: E_total=1430  E_pot=-0.0   P_int=4    clust=1
t=300: E_total=1125  E_pot=-0.0   P_int=2    clust=0
```

## Analysis

### Why the knots failed
1. **Initial binding too weak**: E_pot=-45 (knots) vs -125 (crossed v2) vs -74 (braid3(z))
   The tube radius=2 with A=0.8 gives only moderate field overlap at crossings.

2. **Initial transient too large**: Starting from rest with V(P) far from equilibrium
   creates strong oscillations. Most binding energy radiates away in the first ~25 time units
   (E_pot drops from -45 to -0.3 by t=25).

3. **No traveling wave**: The braid3(z) works because it's a propagating soliton —
   the traveling wave maintains phase coherence. Knots are closed and stationary,
   so there's no propagation to maintain coherence. The standing wave breathes
   and radiates.

### What topology provides
The knot topology DID prevent fragmentation — the structure stayed as 1 cluster
until it dissolved into background (0 clusters at the very end). This is a genuine
improvement over the crossed braids which fragmented into 50+ pieces. But
preventing fragmentation isn't enough if the binding itself dissolves.

### What would be needed
- **Stronger initial binding**: higher amplitude, fatter tube, or different
  phase structure optimized for the knot geometry
- **Closer to equilibrium**: the initial condition should be a quasi-stationary
  state, not an arbitrary ansatz. Could use imaginary-time relaxation (gradient
  flow) to find the true ground state before switching to real-time evolution.
- **Different n_osc**: the number of oscillations around the loop affects the
  V(P) pattern. May need optimization (like v28's CMA-ES search).

## Files

| File | Size | Contents |
|------|------|----------|
| `data/trefoil_R.sfa` | 2.3 GB | f32, 61 frames |
| `data/trefoil_L.sfa` | 2.3 GB | f32, 61 frames |
| `data/figure8.sfa` | 2.3 GB | f32, 61 frames |
| `src/v37_knot.c` | ~850 lines | Simulation code |

### Build and run
```bash
gcc -O3 -march=native -fopenmp -o v37_knot src/v37_knot.c -lzstd -lm

./v37_knot -geom trefoil_R -N 128 -L 25 -T 300 -snap 5 -o data/trefoil_R
./v37_knot -geom trefoil_L -N 128 -L 25 -T 300 -snap 5 -o data/trefoil_L
./v37_knot -geom figure8   -N 128 -L 25 -T 300 -snap 5 -o data/figure8
```
