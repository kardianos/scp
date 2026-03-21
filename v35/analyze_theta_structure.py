#!/usr/bin/env python3
"""Analyze spatial structure of theta field from final snapshot"""
import numpy as np
import sys, os

def load_field(path):
    """Load raw binary field (int32 N, then 3 x N^3 doubles)"""
    with open(path, 'rb') as f:
        N = np.frombuffer(f.read(4), dtype=np.int32)[0]
        data = np.frombuffer(f.read(), dtype=np.float64)
        data = data.reshape(3, N, N, N)
    return N, data

def analyze(dirname, label):
    theta_path = os.path.join(dirname, 'theta_final.bin')
    phi_path = os.path.join(dirname, 'phi_final.bin')

    if not os.path.exists(theta_path):
        print(f"{label}: no theta_final.bin")
        return

    N, theta = load_field(theta_path)
    _, phi = load_field(phi_path)

    L = 20.0
    dx = 2*L/(N-1)
    dV = dx**3

    # RMS values
    theta_rms = np.sqrt(np.mean(theta**2))
    phi_rms = np.sqrt(np.mean(phi**2))

    # Triple product fields
    Q = theta[0] * theta[1] * theta[2]
    P = phi[0] * phi[1] * phi[2]

    # Radial profiles (from center)
    x = np.linspace(-L, L, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    R = np.sqrt(X**2 + Y**2)  # cylindrical radius (braid along z)

    # Bin by cylindrical radius
    r_bins = np.linspace(0, L, 30)
    r_centers = 0.5*(r_bins[:-1] + r_bins[1:])

    theta2_profile = np.zeros(len(r_centers))
    phi2_profile = np.zeros(len(r_centers))
    Q2_profile = np.zeros(len(r_centers))
    count = np.zeros(len(r_centers))

    for ib in range(len(r_centers)):
        mask = (R >= r_bins[ib]) & (R < r_bins[ib+1])
        if np.sum(mask) > 0:
            theta2_profile[ib] = np.mean(np.sum(theta**2, axis=0)[mask])
            phi2_profile[ib] = np.mean(np.sum(phi**2, axis=0)[mask])
            Q2_profile[ib] = np.mean(Q[mask]**2)
            count[ib] = np.sum(mask)

    # Find peak locations
    phi_peak_r = r_centers[np.argmax(phi2_profile)]
    theta_peak_r = r_centers[np.argmax(theta2_profile)]

    # Energy densities near core vs far
    core_mask = R < 5.0
    far_mask = R > 10.0

    theta2_core = np.mean(np.sum(theta**2, axis=0)[core_mask])
    theta2_far = np.mean(np.sum(theta**2, axis=0)[far_mask])
    phi2_core = np.mean(np.sum(phi**2, axis=0)[core_mask])
    phi2_far = np.mean(np.sum(phi**2, axis=0)[far_mask])

    print(f"\n=== {label} (N={N}) ===")
    print(f"  theta_rms = {theta_rms:.4f}")
    print(f"  phi_rms   = {phi_rms:.4f}")
    print(f"  Q_max     = {np.max(np.abs(Q)):.4e}")
    print(f"  P_max     = {np.max(np.abs(P)):.4e}")
    print(f"")
    print(f"  Radial peaks: phi at r={phi_peak_r:.1f}, theta at r={theta_peak_r:.1f}")
    print(f"  Core (r<5) vs far (r>10):")
    print(f"    |phi|²:   core={phi2_core:.4e}  far={phi2_far:.4e}  ratio={phi2_core/(phi2_far+1e-30):.1f}")
    print(f"    |theta|²: core={theta2_core:.4e}  far={theta2_far:.4e}  ratio={theta2_core/(theta2_far+1e-30):.1f}")
    print(f"")

    # Radial profile table
    print(f"  r      |phi|²_avg  |theta|²_avg  Q²_avg")
    for i in range(0, len(r_centers), 3):
        print(f"  {r_centers[i]:5.1f}  {phi2_profile[i]:.4e}  {theta2_profile[i]:.4e}  {Q2_profile[i]:.4e}")

# Analyze baseline and interesting cases
base = "/home/d/code/scp/v35/data/theta_self"
cases = [
    ("mu0_kap50", "baseline (mu_t=0)"),
    ("mu-1000_kap50000", "mu_t=-1000, kap_t=50000"),
    ("mu-10000_kap50000", "mu_t=-10000, kap_t=50000"),
    ("mu-1000_kap5000", "mu_t=-1000, kap_t=5000"),
]

for dirname, label in cases:
    path = os.path.join(base, dirname)
    if os.path.exists(path):
        analyze(path, label)
