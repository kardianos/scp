#!/usr/bin/env python3
"""
viz_proton_comparison.py — Comparison plots for all proton.c experiments.

Generates:
1. Q(t) comparison across all modes and temperatures
2. Energy comparison
3. Key snapshot panels for initial/final states
"""
import os
import sys
import struct
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

DATA_ROOT = os.path.join(os.path.dirname(__file__), '..', 'data')
VIZ_DIR = os.path.join(os.path.dirname(__file__), '..', 'results', 'viz')
os.makedirs(VIZ_DIR, exist_ok=True)


def load_timeseries(subdir):
    """Load timeseries.dat from a proton run directory."""
    path = os.path.join(DATA_ROOT, subdir, 'timeseries.dat')
    if not os.path.exists(path):
        return None
    data = np.loadtxt(path, comments='#')
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return {
        'step': data[:, 0], 'time': data[:, 1],
        'E2': data[:, 2], 'E4': data[:, 3],
        'E_kin': data[:, 4], 'Q': data[:, 5],
        'T_eff': data[:, 6]
    }


def read_snapshot(fname):
    """Read binary snapshot."""
    with open(fname, 'rb') as f:
        N = struct.unpack('i', f.read(4))[0]
        L = struct.unpack('d', f.read(8))[0]
        step = struct.unpack('i', f.read(4))[0]
        time = struct.unpack('d', f.read(8))[0]
        data = np.fromfile(f, dtype=np.float64, count=N**3 * 4)
        q = data.reshape(N, N, N, 4)
    return N, L, step, time, q


def plot_Q_comparison():
    """Main comparison: Q(t) for all experiments."""
    experiments = [
        ('proton_hedgehog_T0', 'Hedgehog T=0', 'C0', '-'),
        ('proton_hedgehog_T1e-5', 'Hedgehog T=1e-5', 'C1', '-'),
        ('proton_hedgehog_T1e-4', 'Hedgehog T=1e-4', 'C2', '-'),
        ('proton_quarks_T0', 'Quarks T=0', 'C3', '--'),
        ('proton_quarks_T1e-5', 'Quarks T=1e-5', 'C4', '--'),
        ('proton_b3_T0_dt005', 'B=3 T=0', 'C5', '-.'),
        ('proton_extract_T0_R1.5', 'Extract R=1.0', 'C6', '-.'),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Proton Field Dynamics: Langevin Simulation Results\n'
                 'Skyrme model (L$_2$+L$_4$), N=128, L=8, $\\sigma$-model',
                 fontsize=14, fontweight='bold')

    # Panel 1: Q(t) — all hedgehog experiments
    ax = axes[0, 0]
    for subdir, label, color, ls in experiments[:3]:
        d = load_timeseries(subdir)
        if d is not None:
            ax.plot(d['time'], d['Q'], label=label, color=color, ls=ls, lw=2)
    ax.axhline(1.0, color='gray', ls=':', alpha=0.5, label='Q=1')
    ax.axhline(0.0, color='gray', ls=':', alpha=0.3)
    ax.set_xlabel('Time (code units)')
    ax.set_ylabel('Topological charge Q')
    ax.set_title('Hedgehog B=1: topology vs temperature')
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.2, 1.1)

    # Panel 2: Q(t) — quarks
    ax = axes[0, 1]
    for subdir, label, color, ls in experiments[3:5]:
        d = load_timeseries(subdir)
        if d is not None:
            ax.plot(d['time'], d['Q'], label=label, color=color, ls=ls, lw=2)
    ax.axhline(0.0, color='gray', ls=':', alpha=0.5)
    ax.set_xlabel('Time (code units)')
    ax.set_ylabel('Q')
    ax.set_title('Three quarks: no topological nucleation')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.05, 0.05)

    # Panel 3: Q(t) — B=3 and extract
    ax = axes[1, 0]
    for subdir, label, color, ls in experiments[5:]:
        d = load_timeseries(subdir)
        if d is not None:
            ax.plot(d['time'], d['Q'], label=label, color=color, ls=ls, lw=2)
    ax.axhline(3.0, color='gray', ls=':', alpha=0.3, label='Q=3')
    ax.axhline(1.0, color='gray', ls=':', alpha=0.3, label='Q=1')
    ax.axhline(0.0, color='gray', ls=':', alpha=0.3)
    ax.set_xlabel('Time (code units)')
    ax.set_ylabel('Q')
    ax.set_title('B=3 dissolution & quark extraction')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Panel 4: Energy — hedgehog T=0 (reference)
    ax = axes[1, 1]
    d = load_timeseries('proton_hedgehog_T0')
    if d is not None:
        ax.plot(d['time'], d['E2'], label='E$_2$', color='blue', lw=1.5)
        ax.plot(d['time'], d['E4'], label='E$_4$', color='red', lw=1.5)
        ax.plot(d['time'], d['E_kin'], label='E$_{kin}$', color='green', lw=1.5)
        ax.plot(d['time'], d['E2']+d['E4']+d['E_kin'],
                label='E$_{total}$', color='black', ls='--', lw=1)
    ax.set_xlabel('Time (code units)')
    ax.set_ylabel('Energy (code units)')
    ax.set_title('Hedgehog T=0: energy components')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    outpath = os.path.join(VIZ_DIR, 'proton_Q_comparison.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved {outpath}')


def plot_dissolution_timeline():
    """Show the dissolution sequence: snapshots at key times."""
    # Try to plot hedgehog T=0 snapshots
    snap_dir = os.path.join(DATA_ROOT, 'proton_hedgehog_T0')
    snaps = sorted([f for f in os.listdir(snap_dir) if f.endswith('.bin')])
    if len(snaps) < 2:
        print(f'Not enough snapshots in {snap_dir}')
        return

    n_panels = min(len(snaps), 5)
    fig, axes = plt.subplots(2, n_panels, figsize=(4*n_panels, 7))
    fig.suptitle('Hedgehog B=1, T=0: Dissolution Sequence (xz-plane, y=0)',
                 fontsize=13, fontweight='bold')

    for col, snap_name in enumerate(snaps[:n_panels]):
        snap_path = os.path.join(snap_dir, snap_name)
        N, L, step, time, q = read_snapshot(snap_path)
        j_mid = N // 2
        h = 2.0 * L / N
        extent = [-L, L, -L, L]

        # Row 0: pion magnitude
        pi_mag = np.sqrt(q[:, j_mid, :, 1]**2 + q[:, j_mid, :, 2]**2 + q[:, j_mid, :, 3]**2)
        ax = axes[0, col]
        im = ax.imshow(pi_mag.T, origin='lower', extent=extent, cmap='inferno',
                       vmin=0, vmax=1.0, aspect='equal')
        ax.set_title(f't={time:.1f}', fontsize=11)
        if col == 0:
            ax.set_ylabel(r'$|\vec{\pi}|$' + '\nz')
        ax.set_xlabel('x')

        # Row 1: scalar field s
        s_slice = q[:, j_mid, :, 0]
        ax = axes[1, col]
        ax.imshow(s_slice.T, origin='lower', extent=extent, cmap='RdBu',
                  vmin=-1, vmax=1, aspect='equal')
        if col == 0:
            ax.set_ylabel('s(x,z)\nz')
        ax.set_xlabel('x')

    plt.tight_layout()
    outpath = os.path.join(VIZ_DIR, 'proton_dissolution_timeline.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved {outpath}')


def plot_mode_comparison_t0():
    """Compare initial states of all four modes."""
    modes = [
        ('proton_hedgehog_T0', 'Hedgehog B=1'),
        ('proton_quarks_T0', 'Three quarks'),
        ('proton_b3_T0_dt005', 'B=3 rational map'),
        ('proton_extract_T0_R1.5', 'Extracted quark'),
    ]

    fig, axes = plt.subplots(2, 4, figsize=(16, 7))
    fig.suptitle('Initial States: Four Initialization Modes (xz-plane, y=0)',
                 fontsize=13, fontweight='bold')

    for col, (subdir, title) in enumerate(modes):
        snap_path = os.path.join(DATA_ROOT, subdir, 'snapshot_0000.bin')
        if not os.path.exists(snap_path):
            continue
        N, L, step, time, q = read_snapshot(snap_path)
        j_mid = N // 2
        extent = [-L, L, -L, L]

        # Row 0: pion magnitude
        pi_mag = np.sqrt(q[:, j_mid, :, 1]**2 + q[:, j_mid, :, 2]**2 + q[:, j_mid, :, 3]**2)
        ax = axes[0, col]
        ax.imshow(pi_mag.T, origin='lower', extent=extent, cmap='inferno',
                  vmin=0, vmax=1.0, aspect='equal')
        ax.set_title(title, fontsize=11)
        if col == 0:
            ax.set_ylabel(r'$|\vec{\pi}|$' + '\nz')

        # Row 1: scalar field
        s_slice = q[:, j_mid, :, 0]
        ax = axes[1, col]
        ax.imshow(s_slice.T, origin='lower', extent=extent, cmap='RdBu',
                  vmin=-1, vmax=1, aspect='equal')
        if col == 0:
            ax.set_ylabel('s(x,z)\nz')

        # Read timeseries to annotate
        d = load_timeseries(subdir)
        if d is not None:
            Q0 = d['Q'][0]
            E0 = d['E2'][0] + d['E4'][0]
            axes[0, col].text(0.05, 0.95, f'Q={Q0:.2f}\nE={E0:.0f}',
                              transform=axes[0, col].transAxes, fontsize=9,
                              va='top', color='white',
                              bbox=dict(boxstyle='round', facecolor='black', alpha=0.5))

    plt.tight_layout()
    outpath = os.path.join(VIZ_DIR, 'proton_initial_modes.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved {outpath}')


def plot_summary_table():
    """Create a summary figure with key numbers from all experiments."""
    experiments = [
        ('proton_hedgehog_T0', 'Hedgehog T=0, $\\gamma$=0'),
        ('proton_hedgehog_T1e-5', 'Hedgehog T=1e-5'),
        ('proton_hedgehog_T1e-4', 'Hedgehog T=1e-4'),
        ('proton_quarks_T0', 'Quarks T=0, $\\gamma$=0'),
        ('proton_quarks_T1e-5', 'Quarks T=1e-5'),
        ('proton_b3_T0_dt005', 'B=3 T=0'),
        ('proton_extract_T0_R1.5', 'Extract R=1.0'),
    ]

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.axis('off')
    fig.suptitle('Proton Langevin Simulation: Summary of Results',
                 fontsize=14, fontweight='bold')

    headers = ['Experiment', 'Q(0)', 'Q(final)', 't_final', 'E_pot(0)',
               'E_pot(final)', 'Result']
    table_data = []
    for subdir, label in experiments:
        d = load_timeseries(subdir)
        if d is None:
            continue
        Q0 = d['Q'][0]
        Qf = d['Q'][-1]
        tf = d['time'][-1]
        Ep0 = d['E2'][0] + d['E4'][0]
        Epf = d['E2'][-1] + d['E4'][-1]

        if abs(Q0) > 0.5 and abs(Qf) < 0.1:
            result = 'DISSOLVED'
        elif abs(Q0) < 0.1 and abs(Qf) < 0.1:
            result = 'NO NUCLEATION'
        elif abs(Qf) > 0.5:
            result = 'PARTIAL'
        else:
            result = 'DISSOLVED'

        table_data.append([label, f'{Q0:.3f}', f'{Qf:.3f}', f'{tf:.1f}',
                           f'{Ep0:.0f}', f'{Epf:.0f}', result])

    table = ax.table(cellText=table_data, colLabels=headers,
                     loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.5)

    # Color results column
    for i, row in enumerate(table_data):
        if row[-1] == 'DISSOLVED':
            table[i+1, 6].set_facecolor('#ffcccc')
        elif row[-1] == 'NO NUCLEATION':
            table[i+1, 6].set_facecolor('#ffffcc')

    plt.tight_layout()
    outpath = os.path.join(VIZ_DIR, 'proton_summary_table.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved {outpath}')


if __name__ == '__main__':
    print('Generating proton comparison plots...')
    plot_Q_comparison()
    plot_dissolution_timeline()
    plot_mode_comparison_t0()
    plot_summary_table()
    print('Done.')
