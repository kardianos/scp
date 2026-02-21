"""
Field Knot Mathematics Exploration

Concept: A 3D+time scalar field that forms self-reinforcing patterns along Celtic knot curves.
Unlike solitons (localized lumps), this is a "string-like" excitation - think of a violin string
shaped into a knot, vibrating with standing waves.

Mathematical Framework:
-----------------------
1. KNOT CURVE: γ(s) ∈ R³, s ∈ [0, 2π)
   Parametric representation of the knot centerline.

2. FIELD DENSITY: ρ(p, t) where p = (x, y, z)
   The field value at any point in space/time.
   
3. SEMI-LOCAL STRUCTURE:
   The field at point p is influenced by nearby portions of the knot curve,
   with phase-dependent interference (constructive/destructive).

Key Equations to Explore:
-------------------------
ρ(p, t) = ∫₀^{2π} A(s) · K(‖p - γ(s)‖) · cos(ωt + φ(s)) ds

Where:
- A(s): amplitude envelope along the knot
- K(r): spatial kernel (Gaussian or similar decay from curve)
- ω: temporal frequency
- φ(s): phase as function of arclength (creates standing wave)

For self-reinforcement, φ(s) must satisfy:
φ(s + Δs) = φ(s) + n · Δs (for some mode number n)

This creates standing wave patterns that "fit" around the knot.
"""

import numpy as np
from dataclasses import dataclass
from typing import Callable, Tuple
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


@dataclass
class KnotParams:
    p: int
    q: int
    r_major: float
    r_minor: float


def torus_knot_curve(s: np.ndarray, params: KnotParams) -> np.ndarray:
    """
    Parametric torus knot: γ(s) = (x(s), y(s), z(s))
    s: parameter in [0, 2π)
    p, q: knot type (trefoil: p=2, q=3)
    r_major, r_minor: torus radii
    """
    phi = params.p * s
    theta = params.q * s
    
    x = (params.r_major + params.r_minor * np.cos(phi)) * np.cos(theta)
    y = (params.r_major + params.r_minor * np.cos(phi)) * np.sin(theta)
    z = params.r_minor * np.sin(phi)
    
    return np.stack([x, y, z], axis=-1)


def figure_eight_knot(s: np.ndarray, scale: float = 0.5) -> np.ndarray:
    """
    Figure-8 (Lissajous) knot: a simpler Celtic-style knot
    """
    x = np.sin(s) * scale
    y = np.sin(2 * s) * scale * 0.5
    z = np.sin(3 * s) * scale * 0.3
    return np.stack([x, y, z], axis=-1)


def celtic_triple(s: np.ndarray, scale: float = 0.4) -> np.ndarray:
    """
    Celtic triple spiral / triskelion-inspired 3D knot
    Three lobes that weave over/under each other
    """
    t = s
    r = 0.3 + 0.2 * np.cos(3 * t)
    x = r * np.cos(t)
    y = r * np.sin(t)
    z = 0.15 * np.sin(3 * t)
    return np.stack([x, y, z], axis=-1) * scale / 0.4


class FieldKnot:
    """
    Represents a field knot: a scalar field with self-reinforcing
    wave-like structure along a parametric curve.
    """
    
    def __init__(
        self,
        curve_func: Callable[[np.ndarray], np.ndarray],
        n_curve_samples: int = 200,
        spatial_decay: float = 0.05,
        wave_number: int = 3,
        amplitude_envelope: Callable[[np.ndarray], np.ndarray] = None,
    ):
        self.curve_func = curve_func
        self.n_samples = n_curve_samples
        self.spatial_decay = spatial_decay
        self.wave_number = wave_number
        
        # Pre-sample the curve
        self.s_values = np.linspace(0, 2 * np.pi, n_curve_samples, endpoint=False)
        self.curve_points = curve_func(self.s_values)
        
        # Amplitude envelope (uniform by default)
        if amplitude_envelope is None:
            self.amplitude = np.ones(n_curve_samples)
        else:
            self.amplitude = amplitude_envelope(self.s_values)
    
    def spatial_kernel(self, r: np.ndarray) -> np.ndarray:
        """
        Gaussian spatial decay from the curve.
        r: distance(s) from curve point(s)
        """
        return np.exp(-r * r / (2 * self.spatial_decay ** 2))
    
    def phase(self, s: np.ndarray, t: float) -> np.ndarray:
        """
        Phase function for standing wave.
        φ(s, t) = ωt + n·s (n is wave_number)
        This creates n wavelengths around the knot.
        """
        omega = 1.0
        return omega * t + self.wave_number * s
    
    def field_density_at_point(self, p: np.ndarray, t: float = 0.0) -> float:
        """
        Compute field density at a single point p = (x, y, z) at time t.
        
        The integral is approximated as a sum over pre-sampled curve points.
        ρ(p, t) = Σᵢ A(sᵢ) · K(‖p - γ(sᵢ)‖) · cos(φ(sᵢ, t))
        """
        # Vector from point to all curve samples
        deltas = p - self.curve_points  # (n_samples, 3)
        distances = np.linalg.norm(deltas, axis=-1)  # (n_samples,)
        
        # Spatial contribution
        spatial = self.spatial_kernel(distances)  # (n_samples,)
        
        # Phase/wave contribution
        phases = self.phase(self.s_values, t)  # (n_samples,)
        wave = np.cos(phases)
        
        # Combined: sum over all curve samples
        density = np.sum(self.amplitude * spatial * wave) / self.n_samples
        
        return density
    
    def field_density_grid(
        self,
        x_range: Tuple[float, float],
        y_range: Tuple[float, float],
        z_range: Tuple[float, float],
        resolution: int = 30,
        t: float = 0.0,
    ) -> np.ndarray:
        """
        Compute field density on a 3D grid.
        """
        x = np.linspace(x_range[0], x_range[1], resolution)
        y = np.linspace(y_range[0], y_range[1], resolution)
        z = np.linspace(z_range[0], z_range[1], resolution)
        
        field = np.zeros((resolution, resolution, resolution))
        
        for i, xi in enumerate(x):
            for j, yj in enumerate(y):
                for k, zk in enumerate(z):
                    p = np.array([xi, yj, zk])
                    field[i, j, k] = self.field_density_at_point(p, t)
        
        return field, (x, y, z)


def visualize_knot_curve(curve_func, n_points=100, title="Knot Curve"):
    """Plot the 3D knot curve"""
    s = np.linspace(0, 2 * np.pi, n_points)
    points = curve_func(s)
    
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    ax.plot(points[:, 0], points[:, 1], points[:, 2], 'b-', linewidth=2)
    ax.scatter(points[0, 0], points[0, 1], points[0, 2], color='red', s=50, label='Start')
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(title)
    ax.legend()
    
    plt.tight_layout()
    return fig, ax


def visualize_field_slice(field_knot: FieldKnot, z_slice: float = 0.0, resolution: int = 50, t: float = 0.0):
    """Visualize a 2D slice of the field density"""
    x = np.linspace(-1, 1, resolution)
    y = np.linspace(-1, 1, resolution)
    X, Y = np.meshgrid(x, y)
    
    field_2d = np.zeros((resolution, resolution))
    
    for i in range(resolution):
        for j in range(resolution):
            p = np.array([X[i, j], Y[i, j], z_slice])
            field_2d[i, j] = field_knot.field_density_at_point(p, t)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(field_2d, extent=[-1, 1, -1, 1], origin='lower', cmap='RdBu')
    plt.colorbar(im, ax=ax, label='Field Density')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title(f'Field Density Slice at Z={z_slice}, t={t:.2f}')
    
    return fig, ax


if __name__ == "__main__":
    # Create different knot types
    trefoil_params = KnotParams(p=2, q=3, r_major=0.4, r_minor=0.15)
    
    print("=" * 60)
    print("FIELD KNOT MATHEMATICS EXPLORATION")
    print("=" * 60)
    
    print("\n1. Creating field knots with different geometries...")
    
    # Trefoil field knot
    trefoil_knot = FieldKnot(
        curve_func=lambda s: torus_knot_curve(s, trefoil_params),
        spatial_decay=0.08,
        wave_number=3,
    )
    
    # Figure-8 field knot
    fig8_knot = FieldKnot(
        curve_func=figure_eight_knot,
        spatial_decay=0.08,
        wave_number=4,
    )
    
    # Celtic triple
    celtic_knot = FieldKnot(
        curve_func=celtic_triple,
        spatial_decay=0.06,
        wave_number=3,
    )
    
    print("   - Trefoil (3-wave)")
    print("   - Figure-8 (4-wave)")
    print("   - Celtic triple (3-wave)")
    
    print("\n2. Testing field density at specific points...")
    
    # Test points
    test_points = [
        np.array([0.0, 0.0, 0.0]),  # Center
        np.array([0.5, 0.0, 0.0]),  # Along X
        np.array([0.0, 0.5, 0.0]),  # Along Y
        np.array([0.3, 0.3, 0.0]),  # Diagonal
    ]
    
    for i, p in enumerate(test_points):
        rho_trefoil = trefoil_knot.field_density_at_point(p, t=0.0)
        rho_fig8 = fig8_knot.field_density_at_point(p, t=0.0)
        print(f"   Point {i+1} {p}: trefoil={rho_trefoil:.4f}, fig8={rho_fig8:.4f}")
    
    print("\n3. Visualizing knot curves...")
    fig1, ax1 = visualize_knot_curve(
        lambda s: torus_knot_curve(s, trefoil_params),
        title="Trefoil Knot Curve"
    )
    fig1.savefig('trefoil_curve.png', dpi=100)
    print("   Saved: trefoil_curve.png")
    
    fig2, ax2 = visualize_knot_curve(
        figure_eight_knot,
        title="Figure-8 Knot Curve"
    )
    fig2.savefig('fig8_curve.png', dpi=100)
    print("   Saved: fig8_curve.png")
    
    print("\n4. Computing field slices...")
    fig3, ax3 = visualize_field_slice(trefoil_knot, z_slice=0.0, resolution=60, t=0.0)
    fig3.savefig('trefoil_field_z0.png', dpi=100)
    print("   Saved: trefoil_field_z0.png")
    
    print("\n" + "=" * 60)
    print("NEXT STEPS:")
    print("- Adjust spatial_decay, wave_number parameters")
    print("- Implement time evolution (dynamical equations)")
    print("- Explore self-reinforcement conditions")
    print("- Add interaction between multiple field knots")
    print("=" * 60)
    
    plt.show()
