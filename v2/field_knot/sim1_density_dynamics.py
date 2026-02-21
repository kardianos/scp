"""
SIMULATION 1: Density Gradient Dynamics and Field Consumption

This simulation demonstrates:
1. Field with embedded topological structure
2. Density gradient driving dynamics
3. Field consumption at knot regions
4. Stability from consumption/refill balance

Key equation:
  ∂Φ/∂t = -γ·∇ρ·Φ/|Φ| + i·V[Φ] + D·∇²Φ
          └────┬────┘   └───┬───┘  └──┬──┘
          density      potential   diffusion
          pressure
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Tuple, List, Optional
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import Normalize
import matplotlib.cm as cm


@dataclass
class EmergentFieldSimulation:
    """
    Simulation of a field with emergent topological structure.
    
    The field evolves according to:
    ∂Φ/∂t = -γ·∇ρ·(Φ/|Φ|)·sign(ρ-ρ_eq) + i·ω·Φ + D·∇²Φ + source
    
    Where:
    - γ: density pressure coefficient
    - ρ = |Φ|²: field density
    - ρ_eq: equilibrium density
    - ω: natural frequency
    - D: diffusion coefficient
    - source: field generation (replenishment)
    """
    n: int = 48                    # Grid size
    extent: float = 2.0            # Physical extent
    dt: float = 0.01               # Time step
    gamma: float = 0.5             # Density pressure coefficient
    omega: float = 2.0             # Natural frequency
    diffusion: float = 0.05        # Diffusion coefficient
    equilibrium_density: float = 0.5  # Equilibrium ρ_eq
    consumption_rate: float = 0.1  # Rate of field consumption at structures
    replenishment: float = 0.02    # Background field generation
    
    def __post_init__(self):
        self.dx = 2 * self.extent / self.n
        
        x = np.linspace(-self.extent, self.extent, self.n)
        self.x = x
        self.X, self.Y, self.Z = np.meshgrid(x, x, x, indexing='ij')
        
        self.R = np.sqrt(self.X**2 + self.Y**2 + self.Z**2)
        
        self.phi = self._initialize_field()
        self.time = 0.0
        
        self.knot_mask = self._detect_knot_regions()
        
        self.history = {
            'time': [],
            'total_density': [],
            'knot_density': [],
            'gradient_energy': [],
            'consumption': [],
            'refill': [],
        }
    
    def _initialize_field(self) -> np.ndarray:
        """Initialize field with embedded topological structure."""
        k = 2
        
        x1 = self.X / (1 + self.R**2)
        x2 = self.Y / (1 + self.R**2)
        x3 = self.Z / (1 + self.R**2)
        x4 = (1 - self.R**2) / (1 + self.R**2)
        
        theta = np.arctan2(np.sqrt(x3**2 + x4**2), np.sqrt(x1**2 + x2**2))
        phi1 = np.arctan2(x2, x1)
        phi2 = np.arctan2(x4, x3)
        
        phase = k * (phi1 + phi2)
        amplitude = np.exp(-0.3 * self.R**2) * 0.8 + 0.2
        
        noise = 0.05 * (np.random.randn(self.n, self.n, self.n) + 
                       1j * np.random.randn(self.n, self.n, self.n))
        
        return amplitude * np.exp(1j * phase) + noise
    
    def _detect_knot_regions(self) -> np.ndarray:
        """Detect regions where topological structure exists."""
        grad_x = np.gradient(self.phi, self.dx, axis=0)
        grad_y = np.gradient(self.phi, self.dx, axis=1)
        grad_z = np.gradient(self.phi, self.dx, axis=2)
        
        phase = np.angle(self.phi)
        
        phase_grad_x = np.gradient(phase, self.dx, axis=0)
        phase_grad_y = np.gradient(phase, self.dx, axis=1)
        phase_grad_z = np.gradient(phase, self.dx, axis=2)
        
        phase_grad_mag = np.sqrt(phase_grad_x**2 + phase_grad_y**2 + phase_grad_z**2)
        
        knot_mask = phase_grad_mag > np.percentile(phase_grad_mag, 70)
        knot_mask = knot_mask & (self.R < 1.0)
        
        return knot_mask
    
    @property
    def density(self) -> np.ndarray:
        return np.abs(self.phi) ** 2
    
    @property
    def phase(self) -> np.ndarray:
        return np.angle(self.phi)
    
    def density_gradient(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Compute ∇ρ"""
        rho = self.density
        gx = np.gradient(rho, self.dx, axis=0)
        gy = np.gradient(rho, self.dx, axis=1)
        gz = np.gradient(rho, self.dx, axis=2)
        return gx, gy, gz
    
    def compute_consumption(self) -> float:
        """Compute rate of field consumption in knot regions."""
        return np.sum(self.consumption_rate * self.density * self.knot_mask) * self.dx**3
    
    def compute_refill(self) -> float:
        """Compute rate of field replenishment."""
        return np.sum(self.replenishment * (1 - self.density)) * self.dx**3
    
    def laplacian(self, f: np.ndarray) -> np.ndarray:
        """Compute ∇²f using finite differences."""
        return (np.roll(f, 1, axis=0) + np.roll(f, -1, axis=0) +
                np.roll(f, 1, axis=1) + np.roll(f, -1, axis=1) +
                np.roll(f, 1, axis=2) + np.roll(f, -1, axis=2) - 6*f) / self.dx**2
    
    def step(self) -> dict:
        """Advance simulation by one time step."""
        rho = self.density
        gx, gy, gz = self.density_gradient()
        
        grad_mag = np.sqrt(gx**2 + gy**2 + gz**2) + 1e-10
        
        amplitude = np.abs(self.phi) + 1e-10
        phi_normalized = self.phi / amplitude
        
        density_pressure = -self.gamma * (rho - self.equilibrium_density)
        
        pressure_term = (density_pressure[..., np.newaxis] * 
                        np.stack([gx, gy, gz], axis=-1))
        pressure_force = np.sum(pressure_term * 
                               np.stack([phi_normalized.real, phi_normalized.imag, 
                                        np.zeros_like(phi_normalized.real)], axis=-1), 
                               axis=-1)
        
        potential_term = 1j * self.omega * self.phi
        
        diffusion_term = self.diffusion * self.laplacian(self.phi)
        
        consumption_term = -self.consumption_rate * self.phi * self.knot_mask
        
        replenishment_term = self.replenishment * (1 - rho / self.equilibrium_density) * np.exp(1j * self.phase)
        
        boundary_mask = self.R < self.extent * 0.9
        boundary_term = 0.01 * (1 - boundary_mask) * (1 - self.phi)
        
        dphi_dt = (pressure_force.astype(complex) + potential_term + 
                   diffusion_term + consumption_term + replenishment_term + boundary_term)
        
        self.phi = self.phi + dphi_dt * self.dt
        self.time += self.dt
        
        self.knot_mask = self._detect_knot_regions()
        
        consumption = self.compute_consumption()
        refill = self.compute_refill()
        
        self.history['time'].append(self.time)
        self.history['total_density'].append(np.sum(rho) * self.dx**3)
        self.history['knot_density'].append(np.sum(rho * self.knot_mask) * self.dx**3)
        self.history['gradient_energy'].append(np.sum(grad_mag**2) * self.dx**3)
        self.history['consumption'].append(consumption)
        self.history['refill'].append(refill)
        
        return {
            'density': rho,
            'phase': self.phase,
            'consumption': consumption,
            'refill': refill,
            'stability': refill / (consumption + 1e-10)
        }
    
    def run(self, n_steps: int) -> List[dict]:
        """Run simulation for n_steps."""
        results = []
        for _ in range(n_steps):
            results.append(self.step())
        return results


def visualize_simulation_state(sim: EmergentFieldSimulation, z_slice: int = None):
    """Create comprehensive visualization of simulation state."""
    if z_slice is None:
        z_slice = sim.n // 2
    
    fig = plt.figure(figsize=(16, 12))
    
    ax1 = fig.add_subplot(2, 3, 1)
    density_slice = sim.density[:, :, z_slice]
    im1 = ax1.imshow(density_slice.T, extent=[-sim.extent, sim.extent, -sim.extent, sim.extent],
                     origin='lower', cmap='viridis')
    ax1.set_title(f'Field Density ρ (z=0, t={sim.time:.2f})')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    plt.colorbar(im1, ax=ax1)
    
    ax2 = fig.add_subplot(2, 3, 2)
    phase_slice = sim.phase[:, :, z_slice]
    im2 = ax2.imshow(phase_slice.T, extent=[-sim.extent, sim.extent, -sim.extent, sim.extent],
                     origin='lower', cmap='hsv', vmin=-np.pi, vmax=np.pi)
    ax2.set_title('Field Phase φ')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    plt.colorbar(im2, ax=ax2)
    
    ax3 = fig.add_subplot(2, 3, 3)
    gx, gy, gz = sim.density_gradient()
    grad_mag = np.sqrt(gx**2 + gy**2 + gz**2)
    grad_slice = grad_mag[:, :, z_slice]
    im3 = ax3.imshow(grad_slice.T, extent=[-sim.extent, sim.extent, -sim.extent, sim.extent],
                     origin='lower', cmap='hot')
    ax3.set_title('Gradient Magnitude |∇ρ|')
    ax3.set_xlabel('x')
    ax3.set_ylabel('y')
    plt.colorbar(im3, ax=ax3)
    
    ax4 = fig.add_subplot(2, 3, 4)
    knot_slice = sim.knot_mask[:, :, z_slice].astype(float)
    im4 = ax4.imshow(knot_slice.T, extent=[-sim.extent, sim.extent, -sim.extent, sim.extent],
                     origin='lower', cmap='Reds')
    ax4.set_title('Knot Regions (High Phase Gradient)')
    ax4.set_xlabel('x')
    ax4.set_ylabel('y')
    plt.colorbar(im4, ax=ax4)
    
    ax5 = fig.add_subplot(2, 3, 5)
    if len(sim.history['time']) > 1:
        ax5.plot(sim.history['time'], sim.history['consumption'], 'r-', label='Consumption')
        ax5.plot(sim.history['time'], sim.history['refill'], 'g-', label='Refill')
        ax5.set_xlabel('Time')
        ax5.set_ylabel('Rate')
        ax5.set_title('Field Consumption vs Replenishment')
        ax5.legend()
        ax5.grid(True, alpha=0.3)
    
    ax6 = fig.add_subplot(2, 3, 6)
    if len(sim.history['time']) > 1:
        stability = np.array(sim.history['refill']) / (np.array(sim.history['consumption']) + 1e-10)
        ax6.plot(sim.history['time'], stability, 'b-')
        ax6.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5)
        ax6.set_xlabel('Time')
        ax6.set_ylabel('Refill / Consumption')
        ax6.set_title('Stability Metric (1 = balanced)')
        ax6.grid(True, alpha=0.3)
        ax6.set_ylim(0, 2)
    
    plt.tight_layout()
    return fig


def visualize_field_consumption(sim: EmergentFieldSimulation, z_slice: int = None):
    """Visualize how the field is consumed around knot regions."""
    if z_slice is None:
        z_slice = sim.n // 2
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    rho = sim.density
    gx, gy, gz = sim.density_gradient()
    
    x = sim.x
    y = sim.x
    X, Y = np.meshgrid(x, y)
    
    ax1 = axes[0, 0]
    density_2d = rho[:, :, z_slice]
    ax1.contourf(X, Y, density_2d.T, levels=20, cmap='viridis')
    ax1.contour(X, Y, density_2d.T, levels=10, colors='white', alpha=0.3, linewidths=0.5)
    ax1.set_title('Field Density with Contours')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_aspect('equal')
    
    ax2 = axes[0, 1]
    skip = sim.n // 12
    gx_2d = gx[:, :, z_slice]
    gy_2d = gy[:, :, z_slice]
    ax2.quiver(X[::skip, ::skip], Y[::skip, ::skip],
               gx_2d[::skip, ::skip].T, gy_2d[::skip, ::skip].T,
               density_2d[::skip, ::skip].T, cmap='coolwarm')
    ax2.set_title('Density Gradient (Flow Direction)')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_aspect('equal')
    
    ax3 = axes[1, 0]
    consumption_field = sim.consumption_rate * rho * sim.knot_mask
    consumption_2d = consumption_field[:, :, z_slice]
    im3 = ax3.imshow(consumption_2d.T, extent=[-sim.extent, sim.extent, -sim.extent, sim.extent],
                     origin='lower', cmap='Reds')
    ax3.set_title('Field Consumption Rate')
    ax3.set_xlabel('x')
    ax3.set_ylabel('y')
    plt.colorbar(im3, ax=ax3, label='∂ρ/∂t (negative)')
    
    ax4 = axes[1, 1]
    depletion = sim.equilibrium_density - rho
    depletion_2d = depletion[:, :, z_slice]
    im4 = ax4.imshow(depletion_2d.T, extent=[-sim.extent, sim.extent, -sim.extent, sim.extent],
                     origin='lower', cmap='RdBu', vmin=-0.5, vmax=0.5)
    ax4.set_title('Field Depletion (ρ_eq - ρ)')
    ax4.set_xlabel('x')
    ax4.set_ylabel('y')
    plt.colorbar(im4, ax=ax4, label='Depletion')
    
    plt.tight_layout()
    return fig


def run_consumption_demo():
    """Demonstrate field consumption dynamics."""
    print("=" * 70)
    print("SIMULATION 1: Density Gradient Dynamics & Field Consumption")
    print("=" * 70)
    
    sim = EmergentFieldSimulation(
        n=48,
        extent=2.0,
        dt=0.005,
        gamma=0.3,
        omega=2.0,
        diffusion=0.03,
        equilibrium_density=0.5,
        consumption_rate=0.15,
        replenishment=0.03
    )
    
    print(f"\nInitial state:")
    print(f"  Total density: {np.sum(sim.density) * sim.dx**3:.4f}")
    print(f"  Knot region density: {np.sum(sim.density * sim.knot_mask) * sim.dx**3:.4f}")
    print(f"  Initial consumption rate: {sim.compute_consumption():.4f}")
    print(f"  Initial refill rate: {sim.compute_refill():.4f}")
    
    fig1 = visualize_simulation_state(sim)
    fig1.savefig('consumption_initial.png', dpi=100)
    print("\nSaved: consumption_initial.png")
    
    print("\nRunning simulation (500 steps)...")
    n_steps = 500
    for i in range(n_steps):
        sim.step()
        if (i + 1) % 100 == 0:
            c = sim.history['consumption'][-1]
            r = sim.history['refill'][-1]
            print(f"  Step {i+1}: consumption={c:.4f}, refill={r:.4f}, ratio={r/c:.2f}")
    
    print(f"\nFinal state (t={sim.time:.2f}):")
    print(f"  Total density: {np.sum(sim.density) * sim.dx**3:.4f}")
    print(f"  Knot region density: {np.sum(sim.density * sim.knot_mask) * sim.dx**3:.4f}")
    
    fig2 = visualize_simulation_state(sim)
    fig2.savefig('consumption_final.png', dpi=100)
    print("Saved: consumption_final.png")
    
    fig3 = visualize_field_consumption(sim)
    fig3.savefig('field_consumption_analysis.png', dpi=100)
    print("Saved: field_consumption_analysis.png")
    
    fig4, axes = plt.subplots(2, 2, figsize=(12, 8))
    
    t = np.array(sim.history['time'])
    
    ax1 = axes[0, 0]
    ax1.plot(t, sim.history['total_density'], 'b-', label='Total')
    ax1.plot(t, sim.history['knot_density'], 'r-', label='Knot Region')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Density')
    ax1.set_title('Density Evolution')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    ax2 = axes[0, 1]
    ax2.plot(t, sim.history['consumption'], 'r-', label='Consumption')
    ax2.plot(t, sim.history['refill'], 'g-', label='Refill')
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Rate')
    ax2.set_title('Consumption vs Replenishment')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    ax3 = axes[1, 0]
    stability = np.array(sim.history['refill']) / (np.array(sim.history['consumption']) + 1e-10)
    ax3.plot(t, stability, 'b-')
    ax3.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5, label='Balanced')
    ax3.fill_between(t, 0, stability, where=stability > 1, alpha=0.3, color='green', label='Growing')
    ax3.fill_between(t, 0, stability, where=stability <= 1, alpha=0.3, color='red', label='Decaying')
    ax3.set_xlabel('Time')
    ax3.set_ylabel('Stability (Refill/Consumption)')
    ax3.set_title('Stability Metric')
    ax3.set_ylim(0, 3)
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    ax4 = axes[1, 1]
    ax4.plot(t, sim.history['gradient_energy'], 'm-')
    ax4.set_xlabel('Time')
    ax4.set_ylabel('∫|∇ρ|² dV')
    ax4.set_title('Gradient Energy')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    fig4.savefig('consumption_evolution.png', dpi=100)
    print("Saved: consumption_evolution.png")
    
    print("\n" + "=" * 70)
    print("KEY OBSERVATIONS:")
    print("- Field density decreases in knot regions (consumption)")
    print("- Density gradient creates 'pressure' toward low-density regions")
    print("- Stability requires refill ≈ consumption")
    print("- Gradient energy represents 'stress' in the field")
    print("=" * 70)
    
    return sim


if __name__ == "__main__":
    import matplotlib
    matplotlib.use('Agg')
    
    sim = run_consumption_demo()
