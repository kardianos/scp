#!/usr/bin/env python3
"""Adaptive simulation controller — runs on the GPU host.

Every T=50: analyze diag, adjust params, print status.
Every T=200: save SFA snapshot.
Runs until killed.
"""
import subprocess, sys, os, math, json, time

BINARY = sys.argv[1] if len(sys.argv) > 1 else "/root/scp_sim_cb959c3e1d7a173e"
SEED = sys.argv[2] if len(sys.argv) > 2 else "/root/v54_seed.sfa"

# Starting parameters (G4-like)
params = {
    'N': 96, 'L': 15.0, 'T': 50.0,
    'm2': 400.0, 'mu': 4.0, 'kappa': 50.0,
    'eta': -4.0, 'A_bg': 0.1,
    'sigma_cross': 20.0, 'lambda_self': 1.5,
    'bc_type': 0, 'damp_width': 5, 'damp_rate': 0.01,
    'bc_switch_time': 0,  # no switch within segments
    'snap_dt': 50.0, 'diag_dt': 2.0, 'vec_dt': 0,
}

segment = 0
total_t = 0.0
current_sfa = SEED
history = []

def write_config(path, p, init_sfa):
    with open(path, 'w') as f:
        for k, v in p.items():
            f.write(f"{k}={v}\n")
        f.write(f"init=sfa\ninit_sfa={init_sfa}\ninit_frame=-1\n")
        f.write(f"output=/root/adaptive_seg.sfa\ndiag_file=/root/adaptive_diag.tsv\n")

def run_segment(cfg_path):
    result = subprocess.run([BINARY, cfg_path], capture_output=True, text=True, timeout=300)
    return result.returncode, result.stdout, result.stderr

def parse_diag(path):
    """Read last 5 lines of diag, return dict of averaged values."""
    lines = []
    with open(path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == len(header):
                lines.append([float(x) for x in parts])
    if not lines:
        return None
    # Average last 5 lines
    n = min(5, len(lines))
    last = lines[-n:]
    avg = {}
    for i, name in enumerate(header):
        avg[name] = sum(row[i] for row in last) / n
    return avg

def analyze_and_adjust(diag, params, segment):
    """Compute metrics and adjust parameters."""
    phi_max = diag.get('phi_max', 0)
    p_max = diag.get('P_max', 0)
    p_int = diag.get('P_int', 0)
    theta_rms = diag.get('theta_rms', 0)
    e_pot = diag.get('E_pot', 0)
    e_total = diag.get('E_total', 0)
    e_coupling = diag.get('E_coupling', 0)
    e_mass = diag.get('E_mass', 0)

    # Tightness metric: P_max / phi_max^3 (how concentrated is binding vs amplitude)
    tightness = p_max / (phi_max**3 + 0.001)

    # Concentration: P_int / (phi_max * N^3 * dx^3) — fraction of volume with binding
    dx = 2*params['L'] / (params['N'] - 1)
    vol = (2*params['L'])**3
    concentration = p_int / (phi_max * vol + 0.001)

    # Energy fraction in coupling vs total
    coupling_frac = abs(e_coupling) / (abs(e_total) + 0.001)

    # Print analysis
    print(f"\n=== Segment {segment} | total_t={total_t:.0f} ===")
    print(f"  phi_max={phi_max:.3f}  P_max={p_max:.4f}  P_int={p_int:.1f}")
    print(f"  theta_rms={theta_rms:.4f}  E_pot={e_pot:.1f}  E_coupling={e_coupling:.0f}")
    print(f"  E_total={e_total:.0f}  E_mass={e_mass:.0f}")
    print(f"  tightness={tightness:.4f}  concentration={concentration:.6f}")
    print(f"  coupling_frac={coupling_frac:.3f}")

    # Parameter adjustment rules
    adjustments = []

    # Rule 1: If phi_max dropping below 0.5, particle is dissolving
    # → increase |eta| to strengthen curl coupling
    if phi_max < 0.5 and segment > 2:
        params['eta'] *= 1.1
        adjustments.append(f"eta → {params['eta']:.3f} (phi too low)")

    # Rule 2: If theta_rms growing > 0.4, theta escaping
    # → increase lambda_self to cap it harder
    if theta_rms > 0.4:
        params['lambda_self'] *= 1.3
        adjustments.append(f"lambda_self → {params['lambda_self']:.3f} (theta too high)")

    # Rule 3: If P_max < 0.1, no binding concentration
    # → increase sigma_cross to co-localize harder
    if p_max < 0.1 and segment > 2:
        params['sigma_cross'] *= 1.2
        adjustments.append(f"sigma_cross → {params['sigma_cross']:.3f} (P_max too low)")

    # Rule 4: If tightness < 0.01, binding is too diffuse
    # → increase m2 to confine phi more
    if tightness < 0.01 and phi_max > 0.3:
        params['m2'] *= 1.1
        adjustments.append(f"m2 → {params['m2']:.1f} (too diffuse)")

    # Rule 5: If phi_max > 3.0, clamp — too wild
    if phi_max > 3.0:
        params['lambda_self'] *= 1.2
        adjustments.append(f"lambda_self → {params['lambda_self']:.3f} (phi too high)")

    # Rule 6: If theta_rms < 0.05, theta is dead — loosen lambda
    if theta_rms < 0.05 and params['lambda_self'] > 0.5:
        params['lambda_self'] *= 0.8
        adjustments.append(f"lambda_self → {params['lambda_self']:.3f} (theta too dead)")

    # Rule 7: Energy drift check
    if len(history) > 2:
        e_prev = history[-2].get('E_total', e_total)
        drift = (e_total - e_prev) / (abs(e_prev) + 1)
        if abs(drift) > 0.05:
            print(f"  WARNING: energy drift {drift*100:.1f}%")

    if adjustments:
        print(f"  ADJUSTMENTS: {'; '.join(adjustments)}")
    else:
        print(f"  No adjustments needed.")

    # Print current parameter state
    print(f"  PARAMS: m2={params['m2']:.1f} mu={params['mu']:.1f} eta={params['eta']:.3f} "
          f"sigma_cross={params['sigma_cross']:.2f} lambda_self={params['lambda_self']:.3f} "
          f"kappa={params['kappa']:.1f}")
    sys.stdout.flush()

    return diag

# After absorbing phase, switch to periodic
ABSORB_SEGMENTS = 1  # First segment only is absorbing

print(f"Adaptive run starting: binary={BINARY}")
print(f"Initial seed: {SEED}")
print(f"Initial params: {json.dumps({k:v for k,v in params.items() if k not in ['N','L','snap_dt','diag_dt','vec_dt','bc_type','damp_width','damp_rate','bc_switch_time']}, indent=2)}")
sys.stdout.flush()

sfa_count = 0

while True:
    segment += 1

    # First segment: absorbing BC. Rest: periodic.
    if segment <= ABSORB_SEGMENTS:
        params['bc_type'] = 0
        params['damp_width'] = 5
        params['damp_rate'] = 0.01
    else:
        params['bc_type'] = 2  # periodic
        params['damp_width'] = 0
        params['damp_rate'] = 0

    # Write config
    cfg_path = '/root/adaptive_cfg.cfg'
    out_sfa = '/root/adaptive_seg.sfa'
    write_config(cfg_path, params, current_sfa)

    # Run segment
    t_start = time.time()
    rc, stdout, stderr = run_segment(cfg_path)
    wall = time.time() - t_start

    if rc != 0:
        print(f"SEGMENT {segment} FAILED (rc={rc})")
        print(stderr[-500:] if stderr else "no stderr")
        break

    total_t += params['T']

    # Parse diag
    diag = parse_diag('/root/adaptive_diag.tsv')
    if diag is None:
        print(f"SEGMENT {segment}: no diag data")
        break

    # Analyze and adjust
    result = analyze_and_adjust(diag, params, segment)
    history.append(result)
    print(f"  Wall: {wall:.1f}s")

    # Save SFA snapshot every 4 segments (T=200)
    if segment % 4 == 0:
        sfa_count += 1
        snap_path = f'/root/adaptive_snap_{sfa_count:04d}.sfa'
        os.system(f'cp {out_sfa} {snap_path}')
        print(f"  SNAPSHOT: {snap_path} (total_t={total_t:.0f})")

    # Use output as next seed
    current_sfa = out_sfa
    sys.stdout.flush()
