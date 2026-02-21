"""
INVESTIGATION MAP 2: Boson-Like Transfers Between Knots

================================================================================
CONCEPT: Bosons as Energy Transfer Patterns
================================================================================

In standard physics:
- Bosons are particles with integer spin (photons, gluons, etc.)
- They mediate forces between fermions

In our field knot theory:
- "Bosons" are NOT particles
- They are PATTERNS of energy/information transfer between knots
- The transfer is localized in space-time but fundamentally a field phenomenon

================================================================================
MATHEMATICAL FRAMEWORK
================================================================================

2.1 ENERGY TRANSFER EQUATION
-----------------------------
Between two knots K₁ and K₂ at positions x₁ and x₂:

  ΔE₁→₂(t) = ∫∫ ρ₁(x,t) · G(x,y) · ρ₂(y,t) · cos(Δφ(x,y,t)) dx dy

Where:
  - G(x,y) = exp(-|x-y|²/λ²) / (4πλ²)^(3/2)  (Yukawa-like propagator)
  - Δφ(x,y,t) = φ₁(x,t) - φ₂(y,t)            (phase difference)
  - λ = characteristic range of interaction

2.2 BOSON IDENTIFICATION
-------------------------
A "boson" B exists when there is LOCALIZED energy transfer:

  B(x,t) = ∂/∂t [∫ G(x,y) · ρ(y,t) · cos(Δφ(x,y,t)) dy]
           └────────────────────────────────────────┘
                    Transfer kernel

This creates a "packet" of energy moving between knots:
  - B > 0: energy flowing INTO location x
  - B < 0: energy flowing OUT OF location x

2.3 QUANTIZATION
-----------------
Energy transfer is effectively quantized by standing wave modes:

  ΔE = n · ℏ_eff · ω
  
Where:
  - n = integer (number of quanta)
  - ℏ_eff = effective Planck constant from field discretization
  - ω = characteristic frequency

In our simulation:
  ℏ_eff ≈ h / (N³)  where N = grid points, h = time step energy

2.4 TRANSFER DYNAMICS
----------------------
The rate of transfer depends on:

  1. SPATIAL OVERLAP: ∫ ρ₁ · ρ₂ dx
     - Closer knots transfer faster
     
  2. PHASE COHERENCE: cos(Δφ)
     - In-phase: positive transfer (attraction)
     - Out-of-phase: negative transfer (repulsion)
     
  3. DENSITY GRADIENT: ∇ρ
     - Energy flows down gradient
     
  4. FIELD CONSUMPTION:
     - Knot consuming field creates "vacancy" that draws energy

================================================================================
ANALOGY WITH PHYSICAL BOSONS
================================================================================

| Property      | Physical Boson    | Field Knot "Boson"       |
|---------------|-------------------|--------------------------|
| Nature        | Particle          | Transfer pattern         |
| Speed         | ≤ c (speed limit) | Determined by λ          |
| Range         | Infinite (massless)| Finite (exponential)    |
| Quantization  | ℏω               | ℏ_eff · ω               |
| Spin          | Integer           | Phase winding (integer)  |
| Creation      | Particle creation | Transfer initiation      |
| Detection     | Absorption        | Energy gain by knot      |

================================================================================
SIMULATION APPROACH
================================================================================

3.1 SETUP
---------
  - Two stable field knots
  - One with excess energy (perturbed)
  - Track energy flow between them

3.2 MEASUREMENTS
----------------
  - E₁(t), E₂(t): energies of each knot
  - B(x,t): boson field (transfer rate)
  - ΔE: total energy transferred
  - τ_transfer: characteristic transfer time

3.3 VISUALIZATION
-----------------
  - Energy vs time for both knots
  - Boson field B(x,y,z) as isosurface
  - Transfer path through field
  - Phase coherence map

================================================================================
KEY PREDICTIONS
================================================================================

1. TRANSFER TIME
   τ_transfer ∝ |x₁ - x₂| / v_eff
   
   Where v_eff is determined by λ and field dynamics

2. ENERGY CONSERVATION
   E_total = E₁ + E₂ + E_field ≈ constant
   
   (Small losses to dissipation)

3. PHASE DEPENDENCE
   - Aligned knots: fast, complete transfer
   - Misaligned: slow, partial transfer
   
4. DISTANCE DEPENDENCE
   Transfer rate ∝ exp(-d/λ)
   
   Where d = |x₁ - x₂|

================================================================================
NEXT: Implement Simulation 2 - Boson Transfer Dynamics
================================================================================
"""

import numpy as np
from dataclasses import dataclass
from typing import Tuple, List, Optional
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D


@dataclass
class TransferAnalysis:
    """Analysis of boson-like transfer between two field regions."""
    
    energy_source_initial: float
    energy_target_initial: float
    energy_source_final: float
    energy_target_final: float
    transfer_efficiency: float
    transfer_time: float
    peak_boson_amplitude: float
    
    def summary(self) -> str:
        return f"""
Transfer Analysis:
  Source: {self.energy_source_initial:.3f} → {self.energy_source_final:.3f}
  Target: {self.energy_target_initial:.3f} → {self.energy_target_final:.3f}
  Efficiency: {self.transfer_efficiency*100:.1f}%
  Transfer time: {self.transfer_time:.2f}
  Peak boson amplitude: {self.peak_boson_amplitude:.4f}
"""


def compute_transfer_kernel(x: np.ndarray, y: np.ndarray, 
                           lambda_range: float = 0.5) -> np.ndarray:
    """
    Compute the transfer kernel G(x,y).
    
    G(x,y) = exp(-|x-y|²/λ²) / normalization
    """
    diff = x - y
    dist_sq = np.sum(diff**2, axis=-1)
    kernel = np.exp(-dist_sq / lambda_range**2)
    return kernel


def compute_boson_field(rho1: np.ndarray, phi1: np.ndarray,
                        rho2: np.ndarray, phi2: np.ndarray,
                        positions: np.ndarray,
                        lambda_range: float = 0.5) -> np.ndarray:
    """
    Compute the boson field B(x) representing energy transfer.
    
    B(x) = ∫ G(x,y) · ρ₁(x) · ρ₂(y) · cos(φ₁(x) - φ₂(y)) dy
    """
    n = rho1.shape[0]
    B = np.zeros((n, n, n))
    
    phase_diff = phi1 - phi2
    coherence = np.cos(phase_diff)
    
    from scipy.ndimage import gaussian_filter
    transfer = rho1 * rho2 * coherence
    
    sigma = lambda_range * n / 4
    B = gaussian_filter(transfer, sigma=sigma)
    
    return B


def visualize_boson_concept():
    """Create visualization explaining the boson concept."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    x = np.linspace(-2, 2, 100)
    y = np.linspace(-2, 2, 100)
    X, Y = np.meshgrid(x, y)
    
    ax1 = axes[0, 0]
    knot1 = np.exp(-((X+0.8)**2 + Y**2) / 0.3)
    ax1.contourf(X, Y, knot1, levels=20, cmap='Blues')
    ax1.contour(X, Y, knot1, levels=[0.5], colors='blue', linewidths=2)
    ax1.set_title('Knot 1 (Source) - High Energy')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    
    knot2 = 0.3 * np.exp(-((X-0.8)**2 + Y**2) / 0.3)
    ax1.contourf(X, Y, knot2, levels=10, cmap='Reds', alpha=0.5)
    ax1.contour(X, Y, knot2, levels=[0.15], colors='red', linewidths=2)
    
    ax2 = axes[0, 1]
    d = np.sqrt((X + 0.8)**2 + Y**2) + np.sqrt((X - 0.8)**2 + Y**2)
    transfer_path = np.exp(-d / 1.5)
    ax2.contourf(X, Y, transfer_path, levels=20, cmap='Greens')
    ax2.arrow(-0.5, 0, 0.8, 0, head_width=0.1, head_length=0.1, fc='black', ec='black')
    ax2.set_title('Boson Field - Transfer Path')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    
    ax3 = axes[1, 0]
    t = np.linspace(0, 10, 100)
    E1 = 1.0 * np.exp(-t/3) + 0.2
    E2 = 0.2 + 0.8 * (1 - np.exp(-t/3))
    ax3.plot(t, E1, 'b-', linewidth=2, label='Knot 1 (Source)')
    ax3.plot(t, E2, 'r-', linewidth=2, label='Knot 2 (Target)')
    ax3.fill_between(t, E1, E2, alpha=0.2, color='green', label='Transfer')
    ax3.axhline(y=0.6, color='gray', linestyle='--', alpha=0.5, label='Equilibrium')
    ax3.set_xlabel('Time')
    ax3.set_ylabel('Energy')
    ax3.set_title('Energy Transfer Over Time')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    ax4 = axes[1, 1]
    distances = np.linspace(0.5, 3, 50)
    transfer_rate = np.exp(-distances / 0.8)
    ax4.plot(distances, transfer_rate, 'g-', linewidth=2)
    ax4.fill_between(distances, transfer_rate, alpha=0.3, color='green')
    ax4.set_xlabel('Distance Between Knots')
    ax4.set_ylabel('Transfer Rate')
    ax4.set_title('Distance Dependence: exp(-d/λ)')
    ax4.grid(True, alpha=0.3)
    ax4.annotate('λ', xy=(0.8, 0.37), fontsize=14, color='red')
    ax4.axvline(x=0.8, color='red', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    return fig


def visualize_phase_dependence():
    """Show how phase difference affects transfer."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    
    x = np.linspace(-2, 2, 100)
    y = np.linspace(-2, 2, 100)
    X, Y = np.meshgrid(x, y)
    
    for idx, (phase_diff, title) in enumerate([
        (0, 'In-Phase (Δφ=0)\nMax Transfer'),
        (np.pi/2, 'Quadrature (Δφ=π/2)\nZero Net Transfer'),
        (np.pi, 'Out-of-Phase (Δφ=π)\nReverse Transfer')
    ]):
        ax = axes[idx]
        
        r1 = np.sqrt((X+0.7)**2 + Y**2)
        r2 = np.sqrt((X-0.7)**2 + Y**2)
        
        rho1 = np.exp(-r1**2 / 0.4)
        rho2 = np.exp(-r2**2 / 0.4)
        
        phi1 = 2 * np.arctan2(Y, X+0.7)
        phi2 = 2 * np.arctan2(Y, X-0.7) + phase_diff
        
        coherence = np.cos(phi1 - phi2)
        
        combined = rho1 + rho2 + 0.5 * coherence * np.sqrt(rho1 * rho2)
        
        if phase_diff == 0:
            cmap = 'RdYlGn'
        elif phase_diff == np.pi:
            cmap = 'RdYlGn_r'
        else:
            cmap = 'coolwarm'
        
        im = ax.contourf(X, Y, combined, levels=20, cmap=cmap)
        ax.contour(X, Y, rho1, levels=[0.3], colors='blue', linewidths=2)
        ax.contour(X, Y, rho2, levels=[0.3], colors='red', linewidths=2)
        
        if phase_diff == 0:
            ax.annotate('', xy=(0, 0), xytext=(-0.5, 0),
                       arrowprops=dict(arrowstyle='->', color='green', lw=2))
            ax.annotate('', xy=(0, 0), xytext=(0.5, 0),
                       arrowprops=dict(arrowstyle='->', color='green', lw=2))
        elif phase_diff == np.pi:
            ax.annotate('', xy=(-0.5, 0), xytext=(0, 0),
                       arrowprops=dict(arrowstyle='->', color='red', lw=2))
            ax.annotate('', xy=(0.5, 0), xytext=(0, 0),
                       arrowprops=dict(arrowstyle='->', color='red', lw=2))
        
        ax.set_title(title)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_aspect('equal')
    
    plt.tight_layout()
    return fig


if __name__ == "__main__":
    import matplotlib
    matplotlib.use('Agg')
    
    print("=" * 70)
    print("INVESTIGATION MAP 2: Boson-Like Transfers")
    print("=" * 70)
    
    print("\nGenerating conceptual visualizations...")
    
    fig1 = visualize_boson_concept()
    fig1.savefig('boson_concept.png', dpi=120)
    print("Saved: boson_concept.png")
    
    fig2 = visualize_phase_dependence()
    fig2.savefig('boson_phase_dependence.png', dpi=120)
    print("Saved: boson_phase_dependence.png")
    
    print("\n" + "=" * 70)
    print("KEY EQUATIONS:")
    print("-" * 70)
    print("Transfer:     ΔE = ∫∫ ρ₁·G·ρ₂·cos(Δφ) dx dy")
    print("Boson field:  B(x) = ∂/∂t [∫ G(x,y)·ρ(y)·cos(Δφ) dy]")
    print("Distance:     Rate ∝ exp(-d/λ)")
    print("Phase:        cos(Δφ) = 1 (max), 0 (none), -1 (reverse)")
    print("=" * 70)
