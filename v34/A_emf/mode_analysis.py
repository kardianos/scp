"""
Phase 1a: Analytical Mode Decomposition of the Background

The background is φ_a(z,t) = A × cos(kz + 2πa/3 - ωt)
with ω = sqrt(k² + m²), k = π/L.

Linearize the equation around this background. The coupling matrix
M_ab = ∂²V/∂φ_a∂φ_b determines how perturbations of the three fields
couple. Its eigenvalues give the effective masses of the three modes.

V(P) = (μ/2)P²/(1+κP²),  P = φ₀φ₁φ₂
V'(P) = μP/(1+κP²)²
V''(P) = μ(1-3κP²)/(1+κP²)³
"""

import numpy as np
from numpy.linalg import eigh
import sys

# Parameters
A = 0.1       # background amplitude
m2 = 2.25     # mass squared
mu = -41.345
kappa = 50.0
L = 20.0
k = np.pi / L  # wavenumber

omega = np.sqrt(k**2 + m2)

print("=" * 70)
print("BACKGROUND MODE ANALYSIS")
print("=" * 70)
print(f"A_bg = {A}, m² = {m2}, μ = {mu}, κ = {kappa}")
print(f"k = π/L = {k:.4f}, ω = √(k²+m²) = {omega:.4f}")
print()

# Sample the background and coupling matrix at many z values
# to understand the time-averaged mode structure
N_z = 1000
z = np.linspace(0, 2*np.pi/k, N_z)  # one full period

# We analyze at t=0 for simplicity (the structure is the same at all t
# due to the standing-wave nature — just rotated in phase)

eigenvalues_all = []
eigenvectors_all = []

for zi in z:
    # Background field values
    phi = np.array([A * np.cos(k*zi + 2*np.pi*a/3) for a in range(3)])

    P = phi[0] * phi[1] * phi[2]

    # V'(P) and V''(P)
    den = 1.0 + kappa * P**2
    Vp = mu * P / den**2
    Vpp = mu * (1 - 3*kappa*P**2) / den**3

    # ∂P/∂φ_a
    dPdphi = np.array([phi[1]*phi[2], phi[0]*phi[2], phi[0]*phi[1]])

    # ∂²P/∂φ_a∂φ_b (off-diagonal only, diagonal = 0)
    d2P = np.zeros((3,3))
    d2P[0,1] = d2P[1,0] = phi[2]
    d2P[0,2] = d2P[2,0] = phi[1]
    d2P[1,2] = d2P[2,1] = phi[0]

    # Coupling matrix M_ab = V'' × (∂P/∂φ_a)(∂P/∂φ_b) + V' × ∂²P/∂φ_a∂φ_b
    M = Vpp * np.outer(dPdphi, dPdphi) + Vp * d2P

    # Eigenvalues (effective mass² corrections from V)
    evals, evecs = eigh(M)
    eigenvalues_all.append(evals)
    eigenvectors_all.append(evecs)

eigenvalues_all = np.array(eigenvalues_all)

print("--- Coupling Matrix M_ab Eigenvalues (V contribution to m_eff²) ---")
print(f"{'':>10s}  {'λ₁ (min)':>12s}  {'λ₂ (mid)':>12s}  {'λ₃ (max)':>12s}")
print(f"{'Mean':>10s}  {eigenvalues_all[:,0].mean():+12.6f}  {eigenvalues_all[:,1].mean():+12.6f}  {eigenvalues_all[:,2].mean():+12.6f}")
print(f"{'Std':>10s}  {eigenvalues_all[:,0].std():12.6f}  {eigenvalues_all[:,1].std():12.6f}  {eigenvalues_all[:,2].std():12.6f}")
print(f"{'Min':>10s}  {eigenvalues_all[:,0].min():+12.6f}  {eigenvalues_all[:,1].min():+12.6f}  {eigenvalues_all[:,2].min():+12.6f}")
print(f"{'Max':>10s}  {eigenvalues_all[:,0].max():+12.6f}  {eigenvalues_all[:,1].max():+12.6f}  {eigenvalues_all[:,2].max():+12.6f}")

print()
print("--- Effective Mass² for Each Mode ---")
print("m_eff² = m² + <λ_i>  (time/space averaged)")
for i, label in enumerate(["Mode 1 (lightest)", "Mode 2 (middle)", "Mode 3 (heaviest)"]):
    meff2 = m2 + eigenvalues_all[:,i].mean()
    if meff2 > 0:
        meff = np.sqrt(meff2)
        yukawa_range = 1.0/meff
        print(f"  {label}: m_eff² = {meff2:+.6f}  m_eff = {meff:.4f}  range = {yukawa_range:.3f}")
    else:
        print(f"  {label}: m_eff² = {meff2:+.6f}  TACHYONIC (imaginary mass)")

print()
print("--- Mode Character Analysis ---")
print("Examining eigenvectors at z where P is maximum and zero:")

# Find z where P is maximum
phi_all = np.array([[A * np.cos(k*zi + 2*np.pi*a/3) for a in range(3)] for zi in z])
P_all = phi_all[:,0] * phi_all[:,1] * phi_all[:,2]
i_max_P = np.argmax(np.abs(P_all))
i_zero_P = np.argmin(np.abs(P_all))

for label, idx in [("P_max", i_max_P), ("P≈0", i_zero_P)]:
    phi = phi_all[idx]
    P = P_all[idx]
    evals = eigenvalues_all[idx]
    evecs = eigenvectors_all[idx]

    print(f"\n  At z={z[idx]:.3f} ({label}): φ = ({phi[0]:.4f}, {phi[1]:.4f}, {phi[2]:.4f}), P = {P:.6f}")
    for i in range(3):
        v = evecs[:,i]
        # Classify: is this mode aligned with φ (amplitude), perpendicular (phase), or mixed?
        phi_norm = phi / (np.linalg.norm(phi) + 1e-30)
        alignment = abs(np.dot(v, phi_norm))

        # Check if it changes Σφ²: δ(Σφ²) = 2Σφ_a × δφ_a = 2 φ·δφ
        delta_phi2 = 2 * np.dot(phi, v)  # proportional to Σφ² change

        # Check if it changes P: δP = Σ (∂P/∂φ_a) δφ_a
        dPdphi = np.array([phi[1]*phi[2], phi[0]*phi[2], phi[0]*phi[1]])
        delta_P = np.dot(dPdphi, v)

        mode_type = "AMPLITUDE" if alignment > 0.8 else "PHASE" if alignment < 0.3 else "MIXED"

        print(f"    Mode {i+1}: λ={evals[i]:+.6f}  v=({v[0]:+.4f},{v[1]:+.4f},{v[2]:+.4f})")
        print(f"            |φ·v|={alignment:.3f} ({mode_type})  δ(Σφ²)∝{delta_phi2:+.4f}  δP∝{delta_P:+.6f}")

# Now check: at larger A_bg, do the modes separate more?
print()
print("=" * 70)
print("AMPLITUDE DEPENDENCE: How do modes change with A_bg?")
print("=" * 70)
print(f"{'A_bg':>8s}  {'<λ₁>':>12s}  {'<λ₂>':>12s}  {'<λ₃>':>12s}  {'m_eff²(1)':>10s}  {'m_eff²(2)':>10s}  {'m_eff²(3)':>10s}  {'Δm²':>8s}")

for A_test in [0.01, 0.05, 0.1, 0.2, 0.5, 1.0]:
    evals_test = []
    for zi in z:
        phi = np.array([A_test * np.cos(k*zi + 2*np.pi*a/3) for a in range(3)])
        P = phi[0]*phi[1]*phi[2]
        den = 1.0 + kappa*P**2
        Vp = mu*P/den**2
        Vpp = mu*(1-3*kappa*P**2)/den**3
        dPdphi = np.array([phi[1]*phi[2], phi[0]*phi[2], phi[0]*phi[1]])
        d2P = np.zeros((3,3))
        d2P[0,1]=d2P[1,0]=phi[2]; d2P[0,2]=d2P[2,0]=phi[1]; d2P[1,2]=d2P[2,1]=phi[0]
        M = Vpp*np.outer(dPdphi, dPdphi) + Vp*d2P
        ev, _ = eigh(M)
        evals_test.append(ev)
    evals_test = np.array(evals_test)
    means = evals_test.mean(axis=0)
    meff = m2 + means
    split = meff[2] - meff[0]
    print(f"{A_test:8.3f}  {means[0]:+12.6f}  {means[1]:+12.6f}  {means[2]:+12.6f}  "
          f"{meff[0]:10.4f}  {meff[1]:10.4f}  {meff[2]:10.4f}  {split:8.4f}")

print()
print("=" * 70)
print("KEY QUESTION: Is any mode massless (m_eff² ≈ 0)?")
print("=" * 70)
print()
print("If all three modes have m_eff² ≈ m² (splitting negligible at A=0.1),")
print("then the 'phonon' observed in the depletion test is NOT a linearized mode.")
print("It must be a NONLINEAR collective effect — the background amplitude")
print("perturbation as a whole, not a single Fourier mode.")
print()
print("If one mode has m_eff² << m²: that mode is the massless phonon/photon")
print("candidate, and its eigenvector tells us whether it's amplitude or phase.")
