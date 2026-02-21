"""
INVESTIGATION MAP: Emergent Field Knot Theory

Goal: Develop a field theory where:
1. Knots EMERGE from field configuration (not placed parametrically)
2. Field density gradients drive dynamics and stability
3. "Bosons" = energy/information transfers between knots
4. Knots CONSUME field, creating density gradients
5. S³ → R³ mapping creates natural topological structure

================================================================================
PHASE 1: MATHEMATICAL FOUNDATIONS
================================================================================

1.1 THE FIELD AND ITS DENSITY
------------------------------
Let Φ(x,t) be a scalar field on R³ × R (3D space + time).

DENSITY: ρ(x,t) = |Φ(x,t)|²
GRADIENT: ∇ρ = 2·Re(Φ*·∇Φ)

Key insight: High density = "vacuum", low density = "knot regions"

1.2 EMERGENT STRUCTURE CONDITION
---------------------------------
A knot emerges where the field satisfies:
  - |Φ| passes through zero (nodal surface)
  - ∇Φ has non-zero winding
  - Phase φ = arg(Φ) has non-trivial topology

The winding number around a closed loop C:
  W(C) = (1/2π) ∮_C ∇φ · dl

If W(C) ≠ 0, there's a topological structure inside C.

1.3 DENSITY GRADIENT DYNAMICS
------------------------------
The field flows toward lower density (vacuum relaxation):
  ∂Φ/∂t = -γ·∇²|Φ|²·Φ/|Φ| + i·V[Φ]  (density-driven)
           └────────────────┘   └───┘
           "pressure" term     potential

This creates:
  - Field "consumed" at high-activity regions
  - Density gradient "pressure" that stabilizes structures

1.4 STABILITY FROM CONSUMPTION
-------------------------------
A knot is STABLE when consumption rate ≈ replenishment rate:
  
  Γ_consume = ∫_knot |∂ρ/∂t| dV
  Γ_refill = ∫_boundary (∇ρ · n̂) dS
  
  Stability: Γ_consume ≈ Γ_refill

If consumption > refill → knot decays (radiates energy)
If consumption < refill → knot grows (absorbs energy)

================================================================================
PHASE 2: S³ → R³ MAPPING AND TOPOLOGY
================================================================================

2.1 THE HOPF MAP
-----------------
Map from 3-sphere S³ ⊂ R⁴ to R³:

  (z₁, z₂) ∈ S³ ⊂ C² → Hopf(z₁, z₂) ∈ R³

Where z₁ = x₁ + i·x₂, z₂ = x₃ + i·x₄ with |z₁|² + |z₂|² = 1

Explicitly:
  X = 2·(z₁·z̄₂ + z̄₁·z₂) = 2(x₁x₃ + x₂x₄)
  Y = 2·i·(z₁·z̄₂ - z̄₁·z₂) = 2(x₂x₃ - x₁x₄)
  Z = |z₁|² - |z₂|² = x₁² + x₂² - x₃² - x₄²

2.2 HOPF FIELD
---------------
Define field on S³, project to R³:
  Φ_S³(θ, φ₁, φ₂) = cos(θ)·e^(iφ₁) + sin(θ)·e^(iφ₂)
  Φ_R³(x) = Φ_S³(Hopf⁻¹(x))

This naturally creates:
  - Linked structures (Hopf fibration)
  - Non-trivial topology
  - Phase singularities = knot cores

2.3 DENSITY FROM HOPF FIELD
----------------------------
ρ(x) = |Φ_R³(x)|² creates natural density variations:
  - ρ = 0 at phase singularities (knot cores)
  - ρ = 1 far from structure
  - ∇ρ points toward/away from cores

================================================================================
PHASE 3: BOSON-LIKE TRANSFERS
================================================================================

3.1 ENERGY TRANSFER AS BOSONS
------------------------------
Between two knots K₁ and K₂, energy transfer is:
  
  ΔE = ∫∫ ρ₁(x)·G(x,y)·ρ₂(y) dx dy
  
Where G(x,y) is the "propagator" (Green's function):
  G(x,y) = exp(-|x-y|²/λ²)·cos(φ₁(x) - φ₂(y))
           └────┬────┘       └──────┬──────┘
           spatial decay     phase coherence

3.2 BOSON IDENTIFICATION
-------------------------
A "boson" is a localized energy packet in transfer:
  
  B(x,t) = ∂/∂t [∫ G(x,y)·ρ(y) dy]
  
This is NOT a particle but a PATTERN of energy flow.
Observable effects:
  - Knot 1 loses energy → Knot 2 gains energy
  - Energy is quantized by standing wave modes
  - Phase coherence determines transfer rate

3.3 EXCHANGE DYNAMICS
----------------------
K₁ ────► K₂ : boson transfer
         B₁₂ = |⟨K₁'|H_int|K₂'⟩|²
         
The "interaction Hamiltonian":
  H_int = ∫ Φ₁*·Φ₂ · V(x₁-x₂) dx₁ dx₂

================================================================================
PHASE 4: SIMULATION ROADMAP
================================================================================

4.1 SIM 1: Density Gradient Field
----------------------------------
- Start with uniform field
- Introduce perturbation
- Watch density gradient emerge
- Measure stability condition

4.2 SIM 2: Emergent Knot
-------------------------
- Field with random initial conditions
- Topological structure emerges spontaneously
- Track winding number W(C)
- Observe density consumption

4.3 SIM 3: Boson Transfer
--------------------------
- Two stable knots
- Perturb one → transfer to other
- Visualize energy flow pattern
- Measure "boson" as localized transfer

4.4 SIM 4: Full Field Dynamics
-------------------------------
- Multiple emergent knots
- Field consumption and replenishment
- Boson exchanges between knots
- Stability from balance

================================================================================
KEY EQUATIONS SUMMARY
================================================================================

Field:              Φ(x,t) ∈ C
Density:            ρ(x,t) = |Φ|²
Gradient:           ∇ρ = 2·Re(Φ*·∇Φ)
Phase:              φ = arg(Φ)
Winding:            W = (1/2π)∮∇φ·dl

Dynamics:
  ∂Φ/∂t = -γ·∇ρ·Φ/|Φ| + i·V[Φ] + D·∇²Φ
          └────┬────┘   └───┬───┘  └──┬──┘
          density      potential   diffusion
          pressure

Stability:
  Γ_consume ≈ Γ_refill
  ∫_knot |∂ρ/∂t| dV ≈ ∫_boundary (∇ρ·n̂) dS

Transfer (boson):
  B(x,t) = ∂/∂t [∫ G(x,y)·ρ(y) dy]
  
================================================================================
NEXT: Implement SIM 1 - Density Gradient Field
================================================================================
"""

import numpy as np
from dataclasses import dataclass
from typing import Tuple, Callable, Optional
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


@dataclass
class FieldConfiguration:
    """A scalar field configuration in 3D."""
    values: np.ndarray  # Complex field values (NX, NY, NZ)
    grid_spacing: float
    extent: float
    
    @property
    def density(self) -> np.ndarray:
        """Field density ρ = |Φ|²"""
        return np.abs(self.values) ** 2
    
    @property
    def phase(self) -> np.ndarray:
        """Field phase φ = arg(Φ)"""
        return np.angle(self.values)
    
    def gradient(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Gradient of field ∇Φ"""
        dx = np.gradient(self.values, self.grid_spacing, axis=0)
        dy = np.gradient(self.values, self.grid_spacing, axis=1)
        dz = np.gradient(self.values, self.grid_spacing, axis=2)
        return dx, dy, dz
    
    def density_gradient(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Gradient of density ∇ρ"""
        rho = self.density
        dx = np.gradient(rho, self.grid_spacing, axis=0)
        dy = np.gradient(rho, self.grid_spacing, axis=1)
        dz = np.gradient(rho, self.grid_spacing, axis=2)
        return dx, dy, dz
    
    def gradient_magnitude(self) -> np.ndarray:
        """|∇ρ|"""
        gx, gy, gz = self.density_gradient()
        return np.sqrt(gx**2 + gy**2 + gz**2)


def hopf_map(z1_real: float, z1_imag: float, z2_real: float, z2_imag: float) -> Tuple[float, float, float]:
    """
    Hopf map S³ → R³
    
    Input: point (z1, z2) on S³ (unit 3-sphere in C²)
    Output: point (X, Y, Z) in R³
    """
    X = 2 * (z1_real * z2_real + z1_imag * z2_imag)
    Y = 2 * (z1_imag * z2_real - z1_real * z2_imag)
    Z = z1_real**2 + z1_imag**2 - z2_real**2 - z2_imag**2
    
    return X, Y, Z


def inverse_hopf_stereographic(x: float, y: float, z: float) -> Tuple[float, float, float, float]:
    """
    Inverse of stereographic projection R³ → S³
    
    Given point in R³, find corresponding point on S³
    """
    r_sq = x**2 + y**2 + z**2
    denom = 1 + r_sq
    
    w = (1 - r_sq) / denom
    scale = 2 / denom
    
    return scale * x, scale * y, scale * z, w


def hopf_field_value(x: float, y: float, z: float, k: float = 2.0) -> complex:
    """
    Compute field value using Hopf fibration structure.
    
    The field has natural topological structure from S³ → R³ mapping.
    """
    x1, x2, x3, x4 = inverse_hopf_stereographic(x, y, z)
    
    theta = np.arctan2(np.sqrt(x3**2 + x4**2), np.sqrt(x1**2 + x2**2))
    phi1 = np.arctan2(x2, x1)
    phi2 = np.arctan2(x4, x3)
    
    phase = k * (phi1 + phi2)
    amplitude = np.exp(-0.5 * (x**2 + y**2 + z**2))
    
    return amplitude * np.exp(1j * phase)


def create_knot_field(
    n: int = 64,
    extent: float = 2.0,
    knot_type: str = "hopf",
    wave_number: int = 2
) -> FieldConfiguration:
    """
    Create a field configuration with embedded topological structure.
    """
    x = np.linspace(-extent, extent, n)
    y = np.linspace(-extent, extent, n)
    z = np.linspace(-extent, extent, n)
    
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    if knot_type == "hopf":
        values = np.zeros((n, n, n), dtype=complex)
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    values[i, j, k] = hopf_field_value(
                        X[i, j, k], Y[i, j, k], Z[i, j, k], 
                        k=wave_number
                    )
    
    elif knot_type == "trefoil":
        s = np.arctan2(Z, np.sqrt(X**2 + Y**2) - 0.5)
        t = np.arctan2(Y, X)
        
        r = np.sqrt(X**2 + Y**2 + Z**2)
        envelope = np.exp(-0.5 * (r - 0.5)**2 / 0.1**2)
        
        phase = wave_number * (2 * s + 3 * t)
        values = envelope * np.exp(1j * phase)
    
    elif knot_type == "uniform":
        values = np.ones((n, n, n), dtype=complex)
    
    else:
        values = np.ones((n, n, n), dtype=complex)
    
    return FieldConfiguration(
        values=values,
        grid_spacing=2 * extent / n,
        extent=extent
    )


def compute_winding_number(field: FieldConfiguration, center: Tuple[float, float, float], 
                           radius: float, normal: str = 'z') -> float:
    """
    Compute winding number of phase around a loop.
    
    W = (1/2π) ∮ ∇φ · dl
    
    This detects topological structure inside the loop.
    """
    n = 50
    angles = np.linspace(0, 2*np.pi, n, endpoint=False)
    
    cx, cy, cz = center
    
    phase_values = []
    for angle in angles:
        if normal == 'z':
            x = cx + radius * np.cos(angle)
            y = cy + radius * np.sin(angle)
            z = cz
        elif normal == 'y':
            x = cx + radius * np.cos(angle)
            y = cy
            z = cz + radius * np.sin(angle)
        else:
            x = cx
            y = cy + radius * np.cos(angle)
            z = cz + radius * np.sin(angle)
        
        ix = int((x + field.extent) / field.grid_spacing)
        iy = int((y + field.extent) / field.grid_spacing)
        iz = int((z + field.extent) / field.grid_spacing)
        
        n_grid = field.values.shape[0]
        if 0 <= ix < n_grid and 0 <= iy < n_grid and 0 <= iz < n_grid:
            phase_values.append(field.phase[ix, iy, iz])
        else:
            phase_values.append(0)
    
    phase_values = np.array(phase_values)
    phase_diff = np.diff(phase_values)
    phase_diff = np.where(phase_diff > np.pi, phase_diff - 2*np.pi, phase_diff)
    phase_diff = np.where(phase_diff < -np.pi, phase_diff + 2*np.pi, phase_diff)
    
    winding = np.sum(phase_diff) / (2 * np.pi)
    
    return winding


def compute_consumption_rate(field: FieldConfiguration, field_prev: FieldConfiguration, 
                            dt: float, mask: np.ndarray) -> float:
    """
    Compute rate of field consumption in a region.
    
    Γ_consume = ∫ |∂ρ/∂t| dV
    """
    rho_current = field.density
    rho_prev = field_prev.density
    
    drho_dt = (rho_current - rho_prev) / dt
    
    consumption = np.sum(np.abs(drho_dt) * mask) * field.grid_spacing**3
    
    return consumption


def compute_refill_rate(field: FieldConfiguration, mask: np.ndarray) -> float:
    """
    Compute rate of field replenishment across boundary.
    
    Γ_refill = ∫_boundary (∇ρ · n̂) dS
    """
    gx, gy, gz = field.density_gradient()
    
    from scipy.ndimage import binary_erosion
    boundary = mask & ~binary_erosion(mask)
    
    grad_mag = np.sqrt(gx**2 + gy**2 + gz**2)
    
    refill = np.sum(grad_mag * boundary) * field.grid_spacing**2
    
    return refill


def visualize_density_gradient(field: FieldConfiguration, z_slice: int = None):
    """Visualize field density and gradient."""
    n = field.values.shape[0]
    if z_slice is None:
        z_slice = n // 2
    
    x = np.linspace(-field.extent, field.extent, n)
    y = np.linspace(-field.extent, field.extent, n)
    X, Y = np.meshgrid(x, y)
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    ax1 = axes[0, 0]
    density_slice = field.density[:, :, z_slice]
    im1 = ax1.imshow(density_slice.T, extent=[-field.extent, field.extent, 
                                               -field.extent, field.extent],
                     origin='lower', cmap='viridis')
    ax1.set_title('Field Density ρ = |Φ|²')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    plt.colorbar(im1, ax=ax1)
    
    ax2 = axes[0, 1]
    phase_slice = field.phase[:, :, z_slice]
    im2 = ax2.imshow(phase_slice.T, extent=[-field.extent, field.extent,
                                            -field.extent, field.extent],
                     origin='lower', cmap='hsv', vmin=-np.pi, vmax=np.pi)
    ax2.set_title('Field Phase φ = arg(Φ)')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    plt.colorbar(im2, ax=ax2)
    
    ax3 = axes[1, 0]
    gx, gy, gz = field.density_gradient()
    gx_slice = gx[:, :, z_slice]
    gy_slice = gy[:, :, z_slice]
    
    skip = n // 16
    ax3.quiver(X[::skip, ::skip], Y[::skip, ::skip],
               gx_slice[::skip, ::skip].T, gy_slice[::skip, ::skip].T)
    ax3.set_title('Density Gradient ∇ρ')
    ax3.set_xlabel('x')
    ax3.set_ylabel('y')
    ax3.set_aspect('equal')
    
    ax4 = axes[1, 1]
    grad_mag = field.gradient_magnitude()
    grad_slice = grad_mag[:, :, z_slice]
    im4 = ax4.imshow(grad_slice.T, extent=[-field.extent, field.extent,
                                           -field.extent, field.extent],
                     origin='lower', cmap='hot')
    ax4.set_title('Gradient Magnitude |∇ρ|')
    ax4.set_xlabel('x')
    ax4.set_ylabel('y')
    plt.colorbar(im4, ax=ax4)
    
    plt.tight_layout()
    return fig


if __name__ == "__main__":
    print("=" * 70)
    print("EMERGENT FIELD KNOT THEORY - Investigation Phase 1")
    print("=" * 70)
    
    print("\n1. Creating field with Hopf structure...")
    field = create_knot_field(n=64, extent=2.0, knot_type="hopf", wave_number=2)
    
    print(f"   Grid size: {field.values.shape}")
    print(f"   Density range: [{field.density.min():.4f}, {field.density.max():.4f}]")
    
    print("\n2. Computing density gradient...")
    gx, gy, gz = field.density_gradient()
    grad_mag = field.gradient_magnitude()
    print(f"   Max |∇ρ|: {grad_mag.max():.4f}")
    print(f"   Mean |∇ρ|: {grad_mag.mean():.4f}")
    
    print("\n3. Computing winding numbers at various locations...")
    test_points = [
        ((0.0, 0.0, 0.0), "Center"),
        ((0.5, 0.0, 0.0), "Off-center X"),
        ((0.0, 0.5, 0.0), "Off-center Y"),
        ((0.0, 0.0, 0.5), "Off-center Z"),
    ]
    
    for (cx, cy, cz), name in test_points:
        W = compute_winding_number(field, (cx, cy, cz), radius=0.3, normal='z')
        print(f"   W at {name}: {W:.3f}")
    
    print("\n4. Visualizing density and gradient...")
    fig = visualize_density_gradient(field, z_slice=32)
    fig.savefig('density_gradient_field.png', dpi=120)
    print("   Saved: density_gradient_field.png")
    
    print("\n" + "=" * 70)
    print("KEY OBSERVATIONS:")
    print("- Density ρ = |Φ|² shows structure of field")
    print("- Phase φ = arg(Φ) has non-trivial topology")
    print("- Winding number W ≠ 0 indicates topological structure")
    print("- Gradient ∇ρ shows 'pressure' direction for dynamics")
    print("=" * 70)
