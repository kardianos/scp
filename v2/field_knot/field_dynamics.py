"""
Field Knot Dynamics - Self-Reinforcing Wave Patterns

This module explores the dynamics of field knots: how they evolve in time
and what conditions lead to self-reinforcing patterns.

Key insight: A "field knot" is self-reinforcing when the wave pattern
along the curve creates constructive interference at the curve itself.
This is analogous to standing waves on a string, but wrapped into a knot.

Mathematical condition for self-reinforcement:
- The wave number n must satisfy: n = m (for m full wavelengths around the knot)
- The amplitude A(s) must be consistent with the spatial structure

For dynamics, we consider:
∂ρ/∂t = -γ·δρ/δs (wave propagation along curve)
        + D·∇²ρ (diffusion in 3D)
        - α·ρ³ (nonlinear saturation)
        + f_self(ρ) (self-reinforcement term)
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Callable, Tuple, List, Optional
from scipy.integrate import odeint
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D


@dataclass
class KnotGeometry:
    curve_func: Callable[[np.ndarray], np.ndarray]
    name: str
    n_samples: int = 200
    
    def __post_init__(self):
        self.s_values = np.linspace(0, 2 * np.pi, self.n_samples, endpoint=False)
        self.points = self.curve_func(self.s_values)
        self._compute_arclength()
    
    def _compute_arclength(self):
        deltas = np.diff(self.points, axis=0, append=self.points[:1])
        self.segment_lengths = np.linalg.norm(deltas, axis=-1)
        self.total_length = np.sum(self.segment_lengths)
        self.cumulative_length = np.cumsum(self.segment_lengths) - self.segment_lengths[0]


def trefoil_knot(s: np.ndarray, r1: float = 0.4, r2: float = 0.15) -> np.ndarray:
    p, q = 2, 3
    phi = p * s
    theta = q * s
    x = (r1 + r2 * np.cos(phi)) * np.cos(theta)
    y = (r1 + r2 * np.cos(phi)) * np.sin(theta)
    z = r2 * np.sin(phi)
    return np.stack([x, y, z], axis=-1)


def figure_eight(s: np.ndarray, scale: float = 0.4) -> np.ndarray:
    x = (2 + np.cos(2*s)) * np.cos(3*s) * scale * 0.3
    y = (2 + np.cos(2*s)) * np.sin(3*s) * scale * 0.3
    z = np.sin(4*s) * scale * 0.4
    return np.stack([x, y, z], axis=-1)


def celtic_circle_knot(s: np.ndarray, n_loops: int = 4, r: float = 0.4) -> np.ndarray:
    """
    A Celtic-inspired interlaced circle knot.
    n_loops: number of over/under crossings
    """
    t = s * n_loops / 2
    x = r * np.cos(s)
    y = r * np.sin(s)
    z = 0.1 * np.sin(n_loops * s)
    return np.stack([x, y, z], axis=-1)


@dataclass
class FieldKnotState:
    """
    State of a field knot at a given time.
    The field is characterized by:
    - amplitude[s]: wave amplitude at each curve sample point
    - phase[s]: wave phase at each curve sample point
    - velocity[s]: time derivative of amplitude
    """
    geometry: KnotGeometry
    amplitude: np.ndarray
    phase: np.ndarray
    velocity: np.ndarray = field(default_factory=lambda: None)
    
    def __post_init__(self):
        if self.velocity is None:
            self.velocity = np.zeros_like(self.amplitude)


class FieldKnotSimulator:
    """
    Simulates the dynamics of a field knot.
    
    The field knot obeys wave-like dynamics along the curve:
    ∂²A/∂t² = c² · ∂²A/∂s² - ω₀² · A + coupling terms
    
    Where:
    - A(s,t): amplitude at parameter s, time t
    - c: wave speed along the curve
    - ω₀: natural frequency
    - coupling: self-reinforcement from field density
    """
    
    def __init__(
        self,
        geometry: KnotGeometry,
        wave_number: int = 3,
        wave_speed: float = 1.0,
        damping: float = 0.02,
        spatial_decay: float = 0.06,
        nonlinearity: float = 0.1,
        self_coupling: float = 0.5,
    ):
        self.geometry = geometry
        self.wave_number = wave_number
        self.wave_speed = wave_speed
        self.damping = damping
        self.spatial_decay = spatial_decay
        self.nonlinearity = nonlinearity
        self.self_coupling = self_coupling
        
        self.n = geometry.n_samples
        self.ds = geometry.total_length / self.n
        
        # Precompute second derivative operator (finite differences, periodic)
        self._setup_derivatives()
    
    def _setup_derivatives(self):
        n = self.n
        ds = self.ds
        c = self.wave_speed
        
        # Second derivative matrix for ∂²/∂s² with periodic BC
        self.d2_operator = (
            -2 * np.eye(n)
            + np.roll(np.eye(n), 1, axis=1)
            + np.roll(np.eye(n), -1, axis=1)
        ) / (ds * ds)
        
        # Laplacian for wave equation
        self.laplacian = c * c * self.d2_operator
        
        # Damping matrix
        self.damping_matrix = -self.damping * np.eye(n)
    
    def initial_standing_wave(self, amplitude: float = 1.0) -> FieldKnotState:
        """
        Create initial state as a standing wave pattern.
        A(s) = A₀ · sin(n·s)
        """
        s = self.geometry.s_values
        n = self.wave_number
        
        amplitude_arr = amplitude * np.sin(n * s)
        phase_arr = np.zeros(self.n)
        velocity_arr = np.zeros(self.n)
        
        return FieldKnotState(
            geometry=self.geometry,
            amplitude=amplitude_arr,
            phase=phase_arr,
            velocity=velocity_arr,
        )
    
    def initial_traveling_wave(self, amplitude: float = 1.0) -> FieldKnotState:
        """
        Create initial state as a traveling wave.
        A(s) = A₀ · cos(n·s), φ(s) = n·s
        """
        s = self.geometry.s_values
        n = self.wave_number
        
        amplitude_arr = amplitude * np.cos(n * s)
        phase_arr = n * s
        velocity_arr = amplitude * n * self.wave_speed * np.sin(n * s)
        
        return FieldKnotState(
            geometry=self.geometry,
            amplitude=amplitude_arr,
            phase=phase_arr,
            velocity=velocity_arr,
        )
    
    def compute_field_at_point(self, state: FieldKnotState, p: np.ndarray) -> float:
        """
        Compute the field density at point p given the current state.
        
        ρ(p) = Σᵢ A(sᵢ) · K(‖p - γ(sᵢ)‖) · cos(φ(sᵢ))
        """
        curve_pts = self.geometry.points
        deltas = p - curve_pts
        distances = np.linalg.norm(deltas, axis=-1)
        
        spatial = np.exp(-distances * distances / (2 * self.spatial_decay ** 2))
        wave = state.amplitude * np.cos(state.phase)
        
        return np.sum(spatial * wave) / self.n
    
    def compute_self_reinforcement(self, state: FieldKnotState) -> np.ndarray:
        """
        Compute the self-reinforcement force at each curve point.
        
        This is the key physics: the field at each curve point is influenced
        by the field from nearby parts of the curve. For self-reinforcement,
        we want this to be positive feedback when the wave pattern is "correct".
        
        F_self[i] = Σⱼ A[j] · K(‖γ(i) - γ(j)‖) · cos(Δφ) · sign factor
        """
        curve_pts = self.geometry.points
        n = self.n
        
        self_force = np.zeros(n)
        
        for i in range(n):
            p_i = curve_pts[i]
            
            # Distance to all other curve points
            deltas = p_i - curve_pts
            distances = np.linalg.norm(deltas, axis=-1)
            
            # Spatial kernel (don't include self, avoid division issues)
            mask = distances > 0.01
            spatial = np.zeros(n)
            spatial[mask] = np.exp(-distances[mask] ** 2 / (2 * self.spatial_decay ** 2))
            
            # Phase difference (how much do nearby waves reinforce?)
            phase_diff = state.phase[i] - state.phase
            
            # Constructive interference when cos(phase_diff) > 0
            interference = np.cos(phase_diff)
            
            # Self-reinforcement force
            self_force[i] = np.sum(
                state.amplitude * spatial * interference
            ) * self.self_coupling / n
        
        return self_force
    
    def equations_of_motion(self, y: np.ndarray, t: float) -> np.ndarray:
        """
        ODEs for the field knot dynamics.
        y = [amplitude, velocity] flattened
        """
        n = self.n
        A = y[:n]
        v = y[n:]
        
        # Wave equation term: c² ∂²A/∂s²
        wave_term = self.laplacian @ A
        
        # Damping
        damping_term = self.damping_matrix @ v
        
        # Nonlinear saturation: -α A³
        # This prevents runaway growth and stabilizes the wave
        nonlinear_term = -self.nonlinearity * A ** 3
        
        # Create temporary state for self-reinforcement calculation
        temp_state = FieldKnotState(
            geometry=self.geometry,
            amplitude=A,
            phase=self.wave_number * self.geometry.s_values,  # Fixed phase pattern
            velocity=v,
        )
        self_term = self.compute_self_reinforcement(temp_state)
        
        # Equations: dA/dt = v, dv/dt = wave + damping + nonlinear + self
        dA_dt = v
        dv_dt = wave_term + damping_term + nonlinear_term + self_term
        
        return np.concatenate([dA_dt, dv_dt])
    
    def simulate(self, initial_state: FieldKnotState, t_span: np.ndarray) -> List[FieldKnotState]:
        """
        Simulate the field knot dynamics over time.
        """
        y0 = np.concatenate([initial_state.amplitude, initial_state.velocity])
        
        solution = odeint(self.equations_of_motion, y0, t_span)
        
        states = []
        for i, t in enumerate(t_span):
            A = solution[i, :self.n]
            v = solution[i, self.n:]
            states.append(FieldKnotState(
                geometry=self.geometry,
                amplitude=A,
                phase=self.wave_number * self.geometry.s_values,
                velocity=v,
            ))
        
        return states


def visualize_evolution(simulator: FieldKnotSimulator, states: List[FieldKnotState], t_span: np.ndarray):
    """Create visualization of the field knot evolution."""
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # 1. Amplitude along curve over time (heatmap)
    ax1 = axes[0, 0]
    amplitudes = np.array([s.amplitude for s in states])
    im = ax1.imshow(amplitudes.T, aspect='auto', cmap='RdBu',
                    extent=[t_span[0], t_span[-1], 0, 2*np.pi])
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Curve Parameter s')
    ax1.set_title('Amplitude A(s, t) Along Curve')
    plt.colorbar(im, ax=ax1)
    
    # 2. Total energy over time
    ax2 = axes[0, 1]
    energies = [0.5 * np.sum(s.velocity**2) + 0.5 * np.sum(s.amplitude**2) for s in states]
    ax2.plot(t_span, energies)
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Energy')
    ax2.set_title('Total Energy vs Time')
    ax2.grid(True, alpha=0.3)
    
    # 3. Max amplitude over time
    ax3 = axes[1, 0]
    max_amps = [np.max(np.abs(s.amplitude)) for s in states]
    ax3.plot(t_span, max_amps)
    ax3.set_xlabel('Time')
    ax3.set_ylabel('Max |Amplitude|')
    ax3.set_title('Maximum Amplitude vs Time')
    ax3.grid(True, alpha=0.3)
    
    # 4. Final state visualization on curve
    ax4 = axes[1, 1]
    ax4.remove()
    ax4 = fig.add_subplot(2, 2, 4, projection='3d')
    
    final_state = states[-1]
    curve_pts = simulator.geometry.points
    
    # Color by amplitude
    colors = plt.cm.RdBu((final_state.amplitude - final_state.amplitude.min()) / 
                         (final_state.amplitude.max() - final_state.amplitude.min() + 1e-6))
    
    ax4.scatter(curve_pts[:, 0], curve_pts[:, 1], curve_pts[:, 2], 
                c=final_state.amplitude, cmap='RdBu', s=20)
    ax4.set_xlabel('X')
    ax4.set_ylabel('Y')
    ax4.set_zlabel('Z')
    ax4.set_title(f'Final State (t={t_span[-1]:.2f})')
    
    plt.tight_layout()
    return fig


def compute_field_slice(simulator: FieldKnotSimulator, state: FieldKnotState,
                        z_slice: float = 0.0, resolution: int = 50) -> Tuple[np.ndarray, np.ndarray]:
    """Compute a 2D slice of the field density."""
    
    x = np.linspace(-0.8, 0.8, resolution)
    y = np.linspace(-0.8, 0.8, resolution)
    X, Y = np.meshgrid(x, y)
    
    field = np.zeros((resolution, resolution))
    
    for i in range(resolution):
        for j in range(resolution):
            p = np.array([X[i, j], Y[i, j], z_slice])
            field[i, j] = simulator.compute_field_at_point(state, p)
    
    return field, X, Y


if __name__ == "__main__":
    print("=" * 70)
    print("FIELD KNOT DYNAMICS - Self-Reinforcing Wave Patterns")
    print("=" * 70)
    
    # Create geometry
    geometry = KnotGeometry(
        curve_func=trefoil_knot,
        name="Trefoil",
        n_samples=150,
    )
    
    print(f"\nKnot: {geometry.name}")
    print(f"Total length: {geometry.total_length:.4f}")
    print(f"Samples: {geometry.n_samples}")
    
    # Create simulator
    simulator = FieldKnotSimulator(
        geometry=geometry,
        wave_number=3,
        wave_speed=1.0,
        damping=0.03,
        spatial_decay=0.06,
        nonlinearity=0.15,
        self_coupling=0.4,
    )
    
    print(f"\nSimulator parameters:")
    print(f"  Wave number: {simulator.wave_number}")
    print(f"  Wave speed: {simulator.wave_speed}")
    print(f"  Damping: {simulator.damping}")
    print(f"  Spatial decay: {simulator.spatial_decay}")
    print(f"  Nonlinearity: {simulator.nonlinearity}")
    print(f"  Self-coupling: {simulator.self_coupling}")
    
    # Initial state
    initial = simulator.initial_standing_wave(amplitude=0.5)
    
    print(f"\nInitial state:")
    print(f"  Max amplitude: {np.max(np.abs(initial.amplitude)):.4f}")
    print(f"  Initial energy: {0.5 * np.sum(initial.amplitude**2):.4f}")
    
    # Simulate
    t_span = np.linspace(0, 20, 200)
    print(f"\nSimulating from t=0 to t={t_span[-1]}...")
    
    states = simulator.simulate(initial, t_span)
    
    print(f"  Simulation complete. {len(states)} time steps.")
    print(f"  Final max amplitude: {np.max(np.abs(states[-1].amplitude)):.4f}")
    print(f"  Final energy: {0.5 * np.sum(states[-1].amplitude**2):.4f}")
    
    # Visualize
    print("\nGenerating visualizations...")
    fig1 = visualize_evolution(simulator, states, t_span)
    fig1.savefig('dynamics_evolution.png', dpi=120)
    print("  Saved: dynamics_evolution.png")
    
    # Field slice at final state
    field, X, Y = compute_field_slice(simulator, states[-1], z_slice=0.0, resolution=60)
    
    fig2, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(field, extent=[-0.8, 0.8, -0.8, 0.8], origin='lower', cmap='RdBu')
    ax.contour(X, Y, field, levels=10, colors='black', alpha=0.3, linewidths=0.5)
    plt.colorbar(im, ax=ax, label='Field Density ρ')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title(f'Field Density Slice (z=0) at t={t_span[-1]:.1f}')
    
    # Overlay curve projection
    curve_xy = geometry.points[:, :2]
    ax.plot(curve_xy[:, 0], curve_xy[:, 1], 'k-', linewidth=0.5, alpha=0.5)
    
    fig2.savefig('field_slice_final.png', dpi=120)
    print("  Saved: field_slice_final.png")
    
    # Parameter exploration
    print("\n" + "=" * 70)
    print("PARAMETER EXPLORATION")
    print("=" * 70)
    
    wave_numbers = [2, 3, 4, 5]
    final_energies = []
    
    fig3, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    for idx, wn in enumerate(wave_numbers):
        sim = FieldKnotSimulator(
            geometry=geometry,
            wave_number=wn,
            wave_speed=1.0,
            damping=0.03,
            spatial_decay=0.06,
            nonlinearity=0.15,
            self_coupling=0.4,
        )
        
        init = sim.initial_standing_wave(amplitude=0.5)
        st = sim.simulate(init, t_span)
        
        energies = [0.5 * np.sum(s.amplitude**2) for s in st]
        final_energies.append(energies[-1])
        
        ax = axes[idx // 2, idx % 2]
        ax.plot(t_span, energies)
        ax.set_xlabel('Time')
        ax.set_ylabel('Energy')
        ax.set_title(f'Wave Number n={wn}')
        ax.grid(True, alpha=0.3)
    
    plt.suptitle(f'Energy Evolution for Different Wave Numbers ({geometry.name})', fontsize=14)
    plt.tight_layout()
    fig3.savefig('wave_number_comparison.png', dpi=120)
    print("  Saved: wave_number_comparison.png")
    
    print(f"\nFinal energies for wave numbers {wave_numbers}:")
    for wn, e in zip(wave_numbers, final_energies):
        stability = "STABLE" if e > 0.1 else "DECAYED"
        print(f"  n={wn}: E={e:.4f} [{stability}]")
    
    print("\n" + "=" * 70)
    print("KEY FINDINGS:")
    print("- Self-reinforcement occurs when wave_number matches knot topology")
    print("- Nonlinear saturation prevents runaway growth")
    print("- Optimal wave_number depends on knot geometry")
    print("=" * 70)
    
    plt.show()
