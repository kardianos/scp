#!/usr/bin/env python3
"""Adaptive simulation controller v2 — real cluster analysis.

Reads the SFA field data directly, computes spatial statistics,
detects actual clusters vs uniform background, and aggressively
tunes parameters to create concentration.

Run on GPU host: python3 -u adaptive_run2.py <binary> <seed>
"""
import subprocess, sys, os, struct, math, json, time, zlib

BINARY = sys.argv[1] if len(sys.argv) > 1 else "/root/scp_sim_cb959c3e1d7a173e"
SEED = sys.argv[2] if len(sys.argv) > 2 else "/root/v54_seed.sfa"

params = {
    'N': 96, 'L': 15.0, 'T': 50.0,
    'm2': 400.0, 'mu': 4.0, 'kappa': 50.0,
    'eta': -4.0, 'A_bg': 0.1,
    'sigma_cross': 20.0, 'lambda_self': 1.5,
    'bc_type': 2, 'damp_width': 0, 'damp_rate': 0,
    'bc_switch_time': 0,
    'snap_dt': 50.0, 'diag_dt': 50.0, 'vec_dt': 0,
}

def write_config(path, p, init_sfa):
    with open(path, 'w') as f:
        for k, v in p.items():
            f.write(f"{k}={v}\n")
        f.write(f"init=sfa\ninit_sfa={init_sfa}\ninit_frame=-1\n")
        f.write(f"output=/root/adapt_seg.sfa\ndiag_file=/root/adapt_diag.tsv\n")

def run_segment(cfg_path):
    result = subprocess.run([BINARY, cfg_path], capture_output=True, text=True, timeout=300)
    return result.returncode

def analyze_sfa(path):
    """Run particle tracker and parse output for real cluster stats."""
    result = subprocess.run(['/root/sfa_particle_track', path],
                          capture_output=True, text=True, timeout=60)
    lines = result.stdout.strip().split('\n')

    # Parse the last frame's particles
    particles = []
    header = None
    last_t = -1
    for line in lines:
        if line.startswith('sfa_particle_track:') or line.startswith('  frame'):
            continue
        parts = line.split('\t')
        if parts[0] == 't':
            header = parts
            continue
        if len(parts) < 10:
            continue
        try:
            t = float(parts[0])
            if t >= last_t:
                if t > last_t:
                    particles = []
                    last_t = t
                particles.append({
                    't': t,
                    'pid': int(parts[1]),
                    'mass': float(parts[2]),
                    'cx': float(parts[3]), 'cy': float(parts[4]), 'cz': float(parts[5]),
                    'rms_r': float(parts[6]),
                    'phi_max': float(parts[7]),
                    'P_peak': float(parts[8]),
                    'E_pot': float(parts[9]),
                    'H_cross': float(parts[10]),
                    'nvox': int(parts[15]) if len(parts) > 15 else 0,
                })
        except (ValueError, IndexError):
            continue

    return particles

def compute_metrics(particles, N, L):
    """Compute real clustering metrics from particle data."""
    total_vox = N * N * N
    box_rms = L * math.sqrt(3) / math.sqrt(5)  # rms_r of uniform distribution in [-L,L]^3

    metrics = {
        'n_clusters': len(particles),
        'has_particle': False,
        'top_mass': 0,
        'top_rms_r': box_rms,
        'top_phi_max': 0,
        'top_P_peak': 0,
        'top_nvox': total_vox,
        'compactness': 0,       # 1.0 = point, 0.0 = fills box
        'concentration': 0,     # fraction of volume in top cluster
        'density_contrast': 0,  # peak/mean
        'total_mass': 0,
        'top_frac': 1.0,        # top cluster mass / total mass
    }

    if not particles:
        return metrics

    # Sort by mass descending
    particles.sort(key=lambda p: p['mass'], reverse=True)
    top = particles[0]
    total_mass = sum(p['mass'] for p in particles)

    metrics['n_clusters'] = len(particles)
    metrics['top_mass'] = top['mass']
    metrics['top_rms_r'] = top['rms_r']
    metrics['top_phi_max'] = top['phi_max']
    metrics['top_P_peak'] = top['P_peak']
    metrics['top_nvox'] = top['nvox']
    metrics['total_mass'] = total_mass

    # Compactness: 1 - rms_r/box_rms. 1.0 = tight point, 0.0 = fills box
    metrics['compactness'] = max(0, 1.0 - top['rms_r'] / box_rms)

    # Concentration: what fraction of volume is the top cluster
    metrics['concentration'] = top['nvox'] / total_vox

    # If top cluster fills >80% of volume, it's not a real particle
    metrics['has_particle'] = (top['nvox'] < total_vox * 0.5 and
                                top['rms_r'] < L * 0.7 and
                                len(particles) >= 1)

    # Density contrast: top cluster mass/volume vs total mass/total volume
    if top['nvox'] > 0:
        top_density = top['mass'] / top['nvox']
        mean_density = total_mass / total_vox
        metrics['density_contrast'] = top_density / (mean_density + 1e-10)

    # Top cluster fraction
    metrics['top_frac'] = top['mass'] / (total_mass + 1e-10)

    return metrics

def adjust_params(metrics, params, segment, history):
    """Aggressively adjust parameters to create clustering."""
    adjustments = []

    has_p = metrics['has_particle']
    compact = metrics['compactness']
    contrast = metrics['density_contrast']
    n_clust = metrics['n_clusters']
    conc = metrics['concentration']

    # AGGRESSIVE: always push toward concentration

    # If no real particle detected (fills the box)
    if not has_p:
        # The field is uniform. Need to break symmetry and concentrate.

        # Increase mu magnitude (stronger triple-product binding)
        if abs(params['mu']) < 200:
            params['mu'] *= 1.15
            adjustments.append(f"mu → {params['mu']:.2f} (no particle, need binding)")

        # Increase eta magnitude (stronger coupling to concentrate)
        if abs(params['eta']) < 20:
            params['eta'] *= 1.08
            adjustments.append(f"eta → {params['eta']:.3f} (no particle, need coupling)")

        # Increase sigma_cross (co-localization)
        if params['sigma_cross'] < 100:
            params['sigma_cross'] *= 1.1
            adjustments.append(f"sigma_cross → {params['sigma_cross']:.2f} (no particle)")

    else:
        # Have a particle — optimize it

        # If compactness is low (particle is loose), increase m2
        if compact < 0.3:
            params['m2'] *= 1.05
            adjustments.append(f"m2 → {params['m2']:.1f} (compact={compact:.3f}, tighten)")

        # If density contrast < 2, field is too uniform
        if contrast < 2.0:
            params['mu'] *= 1.1
            adjustments.append(f"mu → {params['mu']:.2f} (contrast={contrast:.2f}, need more)")

        # If too many clusters, increase sigma_cross to merge them
        if n_clust > 3:
            params['sigma_cross'] *= 1.1
            adjustments.append(f"sigma_cross → {params['sigma_cross']:.2f} ({n_clust} clusters, merge)")

        # If compactness improving, keep going
        if len(history) > 1 and compact > history[-1].get('compactness', 0) + 0.01:
            adjustments.append(f"(compactness improving: {history[-1].get('compactness',0):.3f} → {compact:.3f})")

    # Safety: check energy
    if len(history) > 2:
        # Track if adjustments are making things worse
        prev_compact = history[-1].get('compactness', compact)
        if compact < prev_compact - 0.05:
            adjustments.append(f"WARNING: compactness dropped {prev_compact:.3f} → {compact:.3f}")

    return adjustments

# ---- Main loop ----
segment = 0
total_t = 0.0
current_sfa = SEED
history = []
sfa_count = 0

print(f"Adaptive v2: binary={BINARY}")
print(f"Seed: {SEED}")
print(f"Initial: mu={params['mu']} eta={params['eta']} sigma_cross={params['sigma_cross']} "
      f"lambda_self={params['lambda_self']} m2={params['m2']}")
sys.stdout.flush()

# First segment: absorbing BC
first_absorb = True

while True:
    segment += 1

    if first_absorb:
        params['bc_type'] = 0
        params['damp_width'] = 5
        params['damp_rate'] = 0.01
        first_absorb = False
    else:
        params['bc_type'] = 2
        params['damp_width'] = 0
        params['damp_rate'] = 0

    cfg_path = '/root/adapt_cfg.cfg'
    out_sfa = '/root/adapt_seg.sfa'
    write_config(cfg_path, params, current_sfa)

    t_start = time.time()
    rc = run_segment(cfg_path)
    wall = time.time() - t_start

    if rc != 0:
        print(f"SEGMENT {segment} FAILED (rc={rc})")
        break

    total_t += params['T']

    # Real cluster analysis on the SFA
    particles = analyze_sfa(out_sfa)
    metrics = compute_metrics(particles, params['N'], params['L'])

    # Print comprehensive status
    print(f"\n=== Seg {segment} | t={total_t:.0f} | {wall:.1f}s ===")
    print(f"  clusters={metrics['n_clusters']}  has_particle={metrics['has_particle']}")
    print(f"  top: mass={metrics['top_mass']:.1f} rms_r={metrics['top_rms_r']:.2f} "
          f"phi_max={metrics['top_phi_max']:.3f} P_peak={metrics['top_P_peak']:.4f} "
          f"nvox={metrics['top_nvox']}")
    print(f"  compactness={metrics['compactness']:.4f} contrast={metrics['density_contrast']:.2f} "
          f"concentration={metrics['concentration']:.4f}")

    # Adjust parameters
    adjustments = adjust_params(metrics, params, segment, history)
    if adjustments:
        for a in adjustments:
            print(f"  ADJ: {a}")
    else:
        print(f"  (no adjustments)")

    print(f"  PARAMS: mu={params['mu']:.2f} eta={params['eta']:.3f} "
          f"sigma_cross={params['sigma_cross']:.2f} lambda_self={params['lambda_self']:.3f} "
          f"m2={params['m2']:.1f} kappa={params['kappa']:.1f}")

    history.append(metrics)

    # Snapshot every 4 segments (T=200)
    if segment % 4 == 0:
        sfa_count += 1
        snap = f'/root/adapt2_snap_{sfa_count:04d}.sfa'
        os.system(f'cp {out_sfa} {snap}')
        print(f"  SNAPSHOT: {snap}")

    current_sfa = out_sfa
    sys.stdout.flush()
