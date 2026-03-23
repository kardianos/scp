#!/usr/bin/env python3
"""Wavefront and binding analysis for 3D kappa-dependent collapse simulations.

Compares three BC types: rectangular absorbing, periodic, spherical absorbing.
Reads timeseries.tsv files and analysis outputs.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

DATA = '/home/d/code/scp/v39/data'

def load_ts(subdir):
    """Load timeseries.tsv from a subdirectory."""
    path = os.path.join(DATA, subdir, 'timeseries.tsv')
    data = {}
    with open(path) as f:
        header = f.readline().strip().split('\t')
        cols = {h: [] for h in header}
        for line in f:
            vals = line.strip().split('\t')
            if len(vals) != len(header):
                continue
            for h, v in zip(header, vals):
                cols[h].append(float(v))
    return {h: np.array(v) for h, v in cols.items()}

def load_analysis(filename):
    """Load SUM lines from analysis output."""
    path = os.path.join(DATA, filename)
    sums = []
    clus = []
    with open(path) as f:
        for line in f:
            if line.startswith('SUM\t'):
                parts = line.strip().split('\t')
                sums.append([float(x) for x in parts[1:]])
            elif line.startswith('CLU\t'):
                parts = line.strip().split('\t')
                clus.append(parts[1:])
    return sums, clus

# Load timeseries
ts_rect = load_ts('kappa_3d')
ts_peri = load_ts('kappa_3d_periodic')
ts_sph  = load_ts('kappa_3d_sphere')

# Load analysis summaries
sum_rect, clu_rect = load_analysis('analysis_rect.txt')
sum_peri, clu_peri = load_analysis('analysis_periodic.txt')
sum_sph,  clu_sph  = load_analysis('analysis_sphere.txt')

# ---- Plot 1: E_pot vs time ----
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(ts_rect['t'], ts_rect['E_pot'], 'b-', label='Rect absorbing', linewidth=1.5)
ax.plot(ts_peri['t'], ts_peri['E_pot'], 'r-', label='Periodic', linewidth=1.5)
ax.plot(ts_sph['t'],  ts_sph['E_pot'],  'g-', label='Sphere absorbing', linewidth=1.5)
ax.set_xlabel('Time', fontsize=13)
ax.set_ylabel('E_pot', fontsize=13)
ax.set_title('Potential Energy vs Time (density-kappa, gamma=10)', fontsize=14)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
# Zoom inset for early collapse
axins = ax.inset_axes([0.55, 0.35, 0.4, 0.35])
for ts, c, lab in [(ts_rect, 'b', 'Rect'), (ts_peri, 'r', 'Periodic'), (ts_sph, 'g', 'Sphere')]:
    mask = ts['t'] <= 30
    axins.plot(ts['t'][mask], ts['E_pot'][mask], c+'-', linewidth=1)
axins.set_xlabel('t', fontsize=9)
axins.set_ylabel('E_pot', fontsize=9)
axins.set_title('Early collapse (t<30)', fontsize=9)
axins.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(DATA, 'epot_vs_time.png'), dpi=150)
print('Saved epot_vs_time.png')

# ---- Plot 2: phi_max vs time ----
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(ts_rect['t'], ts_rect['phi_max'], 'b-', label='Rect absorbing', linewidth=1.5)
ax.plot(ts_peri['t'], ts_peri['phi_max'], 'r-', label='Periodic', linewidth=1.5)
ax.plot(ts_sph['t'],  ts_sph['phi_max'],  'g-', label='Sphere absorbing', linewidth=1.5)
ax.set_xlabel('Time', fontsize=13)
ax.set_ylabel('phi_max', fontsize=13)
ax.set_title('Peak Field Amplitude vs Time', fontsize=14)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(DATA, 'phi_max_vs_time.png'), dpi=150)
print('Saved phi_max_vs_time.png')

# ---- Plot 3: P_max vs time ----
fig, ax = plt.subplots(figsize=(10, 6))
ax.semilogy(ts_rect['t'], ts_rect['P_max'], 'b-', label='Rect absorbing', linewidth=1.5)
ax.semilogy(ts_peri['t'], ts_peri['P_max'], 'r-', label='Periodic', linewidth=1.5)
ax.semilogy(ts_sph['t'],  ts_sph['P_max'],  'g-', label='Sphere absorbing', linewidth=1.5)
ax.set_xlabel('Time', fontsize=13)
ax.set_ylabel('P_max (log scale)', fontsize=13)
ax.set_title('Peak |P| (Triple Product) vs Time', fontsize=14)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(DATA, 'pmax_vs_time.png'), dpi=150)
print('Saved pmax_vs_time.png')

# ---- Plot 4: P_int (integrated binding) vs time ----
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(ts_rect['t'], ts_rect['P_int'], 'b-', label='Rect absorbing', linewidth=1.5)
ax.plot(ts_peri['t'], ts_peri['P_int'], 'r-', label='Periodic', linewidth=1.5)
ax.plot(ts_sph['t'],  ts_sph['P_int'],  'g-', label='Sphere absorbing', linewidth=1.5)
ax.set_xlabel('Time', fontsize=13)
ax.set_ylabel('P_int (integrated |P|)', fontsize=13)
ax.set_title('Integrated Binding Measure vs Time', fontsize=14)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(DATA, 'pint_vs_time.png'), dpi=150)
print('Saved pint_vs_time.png')

# ---- Plot 5: E_total (energy conservation check) ----
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(ts_rect['t'], ts_rect['E_total'], 'b-', label='Rect absorbing', linewidth=1.5)
ax.plot(ts_peri['t'], ts_peri['E_total'], 'r-', label='Periodic', linewidth=1.5)
ax.plot(ts_sph['t'],  ts_sph['E_total'],  'g-', label='Sphere absorbing', linewidth=1.5)
ax.set_xlabel('Time', fontsize=13)
ax.set_ylabel('E_total', fontsize=13)
ax.set_title('Total Energy vs Time (conservation check)', fontsize=14)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(DATA, 'etotal_vs_time.png'), dpi=150)
print('Saved etotal_vs_time.png')

# ---- Plot 6: Radial profile comparison at t~20 ----
# Extract radial profiles from analysis files
def load_rad_profiles(filename, target_frame):
    """Load RAD lines for a specific frame."""
    path = os.path.join(DATA, filename)
    r_arr, rho2_arr, P_arr, V_arr = [], [], [], []
    with open(path) as f:
        for line in f:
            if line.startswith(f'RAD\t{target_frame}\t'):
                parts = line.strip().split('\t')
                r_arr.append(float(parts[3]))
                rho2_arr.append(float(parts[4]))
                P_arr.append(float(parts[5]))
                V_arr.append(float(parts[6]))
    return np.array(r_arr), np.array(rho2_arr), np.array(P_arr), np.array(V_arr)

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# t~20 for all three
r_r, rho2_r, P_r, V_r = load_rad_profiles('analysis_rect.txt', 4)  # frame 4 = t=20
r_p, rho2_p, P_p, V_p = load_rad_profiles('analysis_periodic.txt', 4)  # frame 4 = t=20
r_s, rho2_s, P_s, V_s = load_rad_profiles('analysis_sphere.txt', 10)  # frame 10 = t=20

ax = axes[0]
ax.plot(r_r, rho2_r, 'b-', label='Rect', linewidth=1.5)
ax.plot(r_p, rho2_p, 'r-', label='Periodic', linewidth=1.5)
ax.plot(r_s, rho2_s, 'g-', label='Sphere', linewidth=1.5)
ax.set_xlabel('r')
ax.set_ylabel('rho^2 (field amplitude)')
ax.set_title('rho^2(r) at t~20')
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes[1]
ax.plot(r_r, P_r, 'b-', label='Rect', linewidth=1.5)
ax.plot(r_p, P_p, 'r-', label='Periodic', linewidth=1.5)
ax.plot(r_s, P_s, 'g-', label='Sphere', linewidth=1.5)
ax.set_xlabel('r')
ax.set_ylabel('P (triple product)')
ax.set_title('P(r) at t~20')
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes[2]
ax.plot(r_r, V_r, 'b-', label='Rect', linewidth=1.5)
ax.plot(r_p, V_p, 'r-', label='Periodic', linewidth=1.5)
ax.plot(r_s, V_s, 'g-', label='Sphere', linewidth=1.5)
ax.set_xlabel('r')
ax.set_ylabel('V(P) potential')
ax.set_title('V(r) at t~20')
ax.legend()
ax.grid(True, alpha=0.3)

fig.suptitle('Radial Profiles at t~20: BC Comparison', fontsize=14)
fig.tight_layout()
fig.savefig(os.path.join(DATA, 'radial_profiles_t20.png'), dpi=150)
print('Saved radial_profiles_t20.png')

# ---- Plot 7: Periodic BC radial profile evolution ----
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
frames_p = [0, 2, 4, 6, 8, 10]
labels_p = ['t=0', 't=10', 't=20', 't=30', 't=40', 't=50']
colors_p = plt.cm.viridis(np.linspace(0, 1, len(frames_p)))

for ax_idx, (data_name, ylabel) in enumerate([('rho2', 'rho^2'), ('P', 'P'), ('V', 'V(P)')]):
    ax = axes[ax_idx]
    for fi, frame in enumerate(frames_p):
        r, rho2, P, V = load_rad_profiles('analysis_periodic.txt', frame)
        if len(r) == 0:
            continue
        data = {'rho2': rho2, 'P': P, 'V': V}[data_name]
        ax.plot(r, data, color=colors_p[fi], label=labels_p[fi], linewidth=1.5)
    ax.set_xlabel('r')
    ax.set_ylabel(ylabel)
    ax.set_title(f'{ylabel} evolution (periodic BC)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

fig.suptitle('Periodic BC: Radial Profile Evolution', fontsize=14)
fig.tight_layout()
fig.savefig(os.path.join(DATA, 'periodic_radial_evolution.png'), dpi=150)
print('Saved periodic_radial_evolution.png')

# ---- Plot 8: Cluster analysis — number of clusters over time ----
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Extract n_clusters from SUM lines
def get_sum_col(sums, col_idx):
    return np.array([s[col_idx] for s in sums])

ax = axes[0]
# SUM format: t, P_max, P_int, E_pot, rho2_core, rho2_outer, bind_core, bind_outer, R_half, n_clusters
if sum_rect:
    t_r = get_sum_col(sum_rect, 0)
    nc_r = get_sum_col(sum_rect, 9)
    ax.plot(t_r, nc_r, 'bo-', label='Rect', markersize=3)
if sum_peri:
    t_p = get_sum_col(sum_peri, 0)
    nc_p = get_sum_col(sum_peri, 9)
    ax.plot(t_p, nc_p, 'rs-', label='Periodic', markersize=4)
if sum_sph:
    t_s = get_sum_col(sum_sph, 0)
    nc_s = get_sum_col(sum_sph, 9)
    ax.plot(t_s, nc_s, 'g^-', label='Sphere', markersize=3)
ax.set_xlabel('Time')
ax.set_ylabel('Number of clusters (|P|>0.1)')
ax.set_title('Fragmentation: Cluster Count vs Time')
ax.legend()
ax.grid(True, alpha=0.3)

# Core vs outer binding
ax = axes[1]
if sum_rect:
    bc_r = get_sum_col(sum_rect, 6)
    bo_r = get_sum_col(sum_rect, 7)
    ax.plot(t_r, bc_r / (bc_r + bo_r + 1e-30), 'b-', label='Rect core fraction', linewidth=1.5)
if sum_peri:
    bc_p = get_sum_col(sum_peri, 6)
    bo_p = get_sum_col(sum_peri, 7)
    ax.plot(t_p, bc_p / (bc_p + bo_p + 1e-30), 'r-', label='Periodic core fraction', linewidth=1.5)
if sum_sph:
    bc_s = get_sum_col(sum_sph, 6)
    bo_s = get_sum_col(sum_sph, 7)
    ax.plot(t_s, bc_s / (bc_s + bo_s + 1e-30), 'g-', label='Sphere core fraction', linewidth=1.5)
ax.set_xlabel('Time')
ax.set_ylabel('Core binding / Total binding')
ax.set_title('Binding Concentration (core = r < L/3)')
ax.legend()
ax.grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig(os.path.join(DATA, 'cluster_and_binding.png'), dpi=150)
print('Saved cluster_and_binding.png')

print('\nAll plots saved to', DATA)
