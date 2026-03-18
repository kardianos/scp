# V30 Implementation Plan

## Physical Model: FRW Expansion + Rotation

### Modified EOM (cosmological scalar field)

    ∂²φ_a/∂t² + 3H(t)·∂φ_a/∂t = (1/a(t)²)·∇²φ_a - m²φ_a - ∂V/∂φ_a

    V(P) = (μ/2)P²/(1+κP²),  P = φ₀φ₁φ₂
    μ = -41.3, κ = 50, m = 1.50

    Phase 1 (inflation, t < t_inf):
        H(t) = H₀ = const
        a(t) = exp(H₀ × t)

    Phase 2 (formation, t ≥ t_inf):
        H(t) = 0
        a(t) = a_final = exp(H₀ × t_inf) = const

### What the terms do

    3H·∂φ/∂t  — Hubble friction. Redshifts kinetic energy.
                 Velocities decay as ~exp(-3Ht/2) during inflation.
                 Effectively "cools" the field without a thermostat.

    1/a² × ∇²  — Gradient stretching. Physical wavelengths grow as a(t).
                  Short-wavelength modes are redshifted below the mass gap.
                  Only long-wavelength modes (comparable to Hubble radius)
                  survive inflation.

### Initial condition: dense rotating universe-filling field

    φ_a(x,0) = A₀ × cos(k₀z + 2πa/3 + n_wind × atan2(y,x))
    v_a(x,0) = ω₀ × A₀ × sin(k₀z + 2πa/3 + n_wind × atan2(y,x))

    A₀ = 3.0       (dense, ~4× equilibrium braid amplitude)
    k₀ = π/L       (one axial half-twist, periodic in z)
    ω₀ = √(k₀²+m²) (dispersion)
    n_wind = 0,1,2,3,5  (azimuthal phase windings = angular momentum)

    Periodic BC in ALL directions.

### Parameters

    N = 256         (grid: 16.8M points, 1.2GB for 3 fields)
    L = 100         (domain: [-100, 100]³, dx = 0.782)
    H₀ = 0.05      (Hubble rate → a_final = e^(0.05×50) ≈ 12.2)
    t_inf = 50      (inflation duration)
    T_total = 350   (50 inflation + 300 formation)
    dt = 0.15       (CFL: 0.2×dx=0.156, use 0.15 for safety)

## Numerical Method

### Modified Velocity Verlet with Hubble friction

The friction term 3Hv makes standard Verlet slightly tricky.
Use the exponential integrator approach:

    Half-kick with friction:
        v += (dt/2) × (acc - 3H×v)
      = v × (1 - 3H×dt/2) + (dt/2) × acc     [first-order in dt]

    Better: exact friction integration for the kick:
        v_new = v × exp(-3H×dt/2) + (dt/2) × acc × exp(-3H×dt/4)

    For simplicity, use the first-order version (accurate when H×dt << 1):
        H₀ × dt = 0.05 × 0.15 = 0.0075  →  3H×dt/2 = 0.011  (small, OK)

    Algorithm per step:
        1. half_kick:  v += (dt/2) × acc;  v *= (1 - 3H×dt/2)
        2. drift:      φ += dt × v
        3. compute acc: acc = (1/a²)∇²φ - m²φ - ∂V/∂φ
        4. half_kick:  v += (dt/2) × acc;  v *= (1 - 3H×dt/2)
        5. update a:   if t < t_inf: a = exp(H₀×t)

    The 1/a² in the Laplacian: just multiply the standard Laplacian by 1/a².

## Defect Detection

After inflation (t > t_inf), scan for localized structures:

    1. Compute ρ(x) = Σ_a [½v_a² + ½(∇φ_a)²/a² + ½m²φ_a² + V(P)]
    2. Smooth ρ over 3³ cells (box filter)
    3. Threshold: ρ_thresh = 3 × median(ρ)  [peaks well above background]
    4. Find connected regions above threshold (flood fill)
    5. For each region:
       - centroid (center of mass weighted by ρ)
       - total energy (integral of ρ over region)
       - winding number (phase along z through centroid)
       - spatial extent (RMS radius)
       - angular momentum proxy: ∫ (x×∇)·(φ₀,φ₁) around centroid

    Track particles across snapshots by nearest-centroid matching.

## Snapshot Strategy

### Per-diagnostic (every T=5, ~60 snapshots):
    - timeseries.tsv: t, a(t), E_total, E_kinetic, E_gradient, E_mass,
                       E_pot, max_rho, avg_rho, n_particles
    - particles_tNNN.tsv: id, x, y, z, E, winding, size (per particle)

### 2D slices (every T=10, ~30 snapshots):
    - slice_xy_tNNN.tsv: x, y, rho (at z=z_mid)
    - slice_xz_tNNN.tsv: x, z, rho (at y=y_mid)
    Format: TSV with header, one row per (i,j) pair, subsampled 2:1

    At N=256, subsampled 2:1: 128² = 16K rows × 3 cols = ~400KB per slice.
    30 snapshots × 2 slices × 400KB = ~24MB total. Trivial.

### Full 3D dumps (5 total):
    - full3d_tNNN.bin: binary, raw float64[N³] for ρ(x,y,z)
    - At T=0, t_inf, t_inf+50, t_inf+150, T_final
    - Size: 256³ × 8 = 134MB each, 670MB total
    - Header: int32 N, float64 L, float64 t, then raw data
    - Loadable in Python: np.fromfile(f, dtype=np.float64, count=N**3).reshape(N,N,N)

## Output Structure

    v30/data/
    ├── nw0/                      (n_wind=0, control)
    │   ├── timeseries.tsv
    │   ├── particles_t050.tsv    (at t_inf)
    │   ├── particles_t100.tsv
    │   ├── ...
    │   ├── slice_xy_t000.tsv
    │   ├── slice_xy_t010.tsv
    │   ├── ...
    │   ├── full3d_t000.bin
    │   ├── full3d_t050.bin
    │   ├── full3d_t100.bin
    │   ├── full3d_t200.bin
    │   └── full3d_t350.bin
    ├── nw1/                      (n_wind=1)
    │   └── (same structure)
    ├── nw2/
    ├── nw3/
    └── nw5/

## Runtime Estimate

    Per step: 16.8M points × ~30 FLOP = 504M FLOP
    At 50 GFLOP/s (16-core OMP): 0.01s/step
    Plus diagnostics overhead: ~0.05s/step average
    2333 steps × 0.05s = ~117s per run
    With snapshots + particle detection: ~3 min per run
    5 n_wind values: ~15 min total

## Code Structure (v30_expand.c)

    main():
        parse args (n_wind, H0, etc.)
        create output directory
        initialize field (dense + rotating)
        open output files

        // Phase 1: Inflation
        for step in 0..steps_inf:
            update a(t), H(t)
            verlet_step_frw(g, a, H)
            if diagnostic: measure + write

        // Phase 2: Formation (H=0, a=const)
        for step in steps_inf..steps_total:
            verlet_step_frw(g, a_final, 0)
            if diagnostic: measure + write
            if particle_detect_time: find_particles + track

        write summary
        close files
