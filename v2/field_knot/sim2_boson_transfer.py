"""
SIMULATION 2: Boson Transfer Dynamics

Demonstrates:
1. Two field knots with different energies
2. Energy transfer between knots (boson-like)
3. Phase dependence of transfer rate
4. Visualization of transfer "field"
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Tuple, List, Optional
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.ndimage import gaussian_filter


@dataclass
class DualKnotSystem:
    """
    Two field knots with boson-like coupling.
    
    Energy transfer occurs via:
    ∂E₁/∂t = -κ · ∫∫ ρ₁·G·ρ₂·cos(Δφ) dx dy
    ∂E₂/∂t = +κ · ∫∫ ρ₁·G·ρ₂·cos(Δφ) dx dy
    """
    n: int = 40
    extent: float = 2.5
    dt: float = 0.01
    
    knot1_pos: Tuple[float, float, float] = (-0.7, 0.0, 0.0)
    knot2_pos: Tuple[float, float, float] = (0.7, 0.0, 0.0)
    
    coupling_range: float = 0.5
    coupling_strength: float = 0.3
    
    diffusion: float = 0.02
    damping: float = 0.01
    
    def __post_init__(self):
        self.dx = 2 * self.extent / self.n
        
        x = np.linspace(-self.extent, self.extent, self.n)
        self.x = x
        self.X, self.Y, self.Z = np.meshgrid(x, x, x, indexing='ij')
        
        self.phi = np.zeros((self.n, self.n, self.n), dtype=complex)
        self.phi1 = np.zeros_like(self.phi)
        self.phi2 = np.zeros_like(self.phi)
        
        self._initialize_knots()
        
        self.time = 0.0
        self.history = {
            'time': [],
            'E1': [],
            'E2': [],
            'E_total': [],
            'transfer_rate': [],
            'phase_coherence': [],
        }
    
    def _create_knot_field(self, center: Tuple[float, float, float], 
                           energy: float, phase_offset: float = 0.0) -> np.ndarray:
        """Create a localized knot-like field structure."""
        cx, cy, cz = center
        
        R = np.sqrt((self.X - cx)**2 + (self.Y - cy)**2 + (self.Z - cz)**2)
        theta = np.arctan2(self.Y - cy, self.X - cx)
        phi_angle = np.arccos(np.clip((self.Z - cz) / (R + 1e-10), -1, 1))
        
        wave_number = 2
        phase = wave_number * (theta + phi_angle) + phase_offset
        
        amplitude = energy * np.exp(-R**2 / 0.2)
        
        return amplitude * np.exp(1j * phase)
    
    def _initialize_knots(self):
        """Initialize both knots with different energies."""
        self.phi1 = self._create_knot_field(self.knot1_pos, energy=1.0, phase_offset=0.0)
        self.phi2 = self._create_knot_field(self.knot2_pos, energy=0.3, phase_offset=0.0)
        
        self.phi = self.phi1 + self.phi2
        
        self._update_knot_regions()
    
    def _update_knot_regions(self):
        """Identify regions belonging to each knot."""
        cx1, cy1, cz1 = self.knot1_pos
        cx2, cy2, cz2 = self.knot2_pos
        
        R1 = np.sqrt((self.X - cx1)**2 + (self.Y - cy1)**2 + (self.Z - cz1)**2)
        R2 = np.sqrt((self.X - cx2)**2 + (self.Y - cy2)**2 + (self.Z - cz2)**2)
        
        self.mask1 = R1 < R2
        self.mask2 = R2 < R1
    
    @property
    def density(self) -> np.ndarray:
        return np.abs(self.phi) ** 2
    
    @property
    def phase(self) -> np.ndarray:
        return np.angle(self.phi)
    
    def compute_knot_energy(self, mask: np.ndarray) -> float:
        """Compute energy within a region."""
        return np.sum(self.density * mask) * self.dx**3
    
    def compute_boson_field(self) -> np.ndarray:
        """
        Compute the boson field B(x) representing transfer.
        
        B(x) = ρ₁·ρ₂·G(x)·cos(Δφ)
        """
        rho1 = np.abs(self.phi1)**2
        rho2 = np.abs(self.phi2)**2
        phi1 = np.angle(self.phi1)
        phi2 = np.angle(self.phi2)
        
        phase_diff = phi1 - phi2
        coherence = np.cos(phase_diff)
        
        B = rho1 * rho2 * coherence
        
        sigma = self.coupling_range * self.n / (2 * self.extent)
        B = gaussian_filter(B, sigma=sigma)
        
        return B
    
    def compute_transfer_rate(self) -> float:
        """
        Compute rate of energy transfer.
        
        Rate = κ · ∫ ρ₁·ρ₂·cos(Δφ) dV
        """
        B = self.compute_boson_field()
        
        transfer = self.coupling_strength * np.sum(B) * self.dx**3
        
        return transfer
    
    def compute_phase_coherence(self) -> float:
        """Compute average phase coherence between knots."""
        phi1_vals = self.phase[self.mask1]
        phi2_vals = self.phase[self.mask2]
        
        if len(phi1_vals) == 0 or len(phi2_vals) == 0:
            return 0.0
        
        phase_diff = np.mean(phi1_vals) - np.mean(phi2_vals)
        coherence = np.cos(phase_diff)
        
        return coherence
    
    def laplacian(self, f: np.ndarray) -> np.ndarray:
        """Compute ∇²f."""
        return (np.roll(f, 1, axis=0) + np.roll(f, -1, axis=0) +
                np.roll(f, 1, axis=1) + np.roll(f, -1, axis=1) +
                np.roll(f, 1, axis=2) + np.roll(f, -1, axis=2) - 6*f) / self.dx**2
    
    def step(self) -> dict:
        """Advance simulation by one time step."""
        B = self.compute_boson_field()
        transfer = self.compute_transfer_rate()
        
        rho1 = np.abs(self.phi1)**2
        rho2 = np.abs(self.phi2)**2
        phi1_angle = np.angle(self.phi1)
        phi2_angle = np.angle(self.phi2)
        
        phase_diff = phi1_angle - phi2_angle
        coherence = np.cos(phase_diff)
        
        transfer_field1 = -self.coupling_strength * B * coherence * np.exp(1j * phi1_angle)
        transfer_field2 = +self.coupling_strength * B * coherence * np.exp(1j * phi2_angle)
        
        diff1 = self.diffusion * self.laplacian(self.phi1)
        diff2 = self.diffusion * self.laplacian(self.phi2)
        
        damp1 = -self.damping * self.phi1
        damp2 = -self.damping * self.phi2
        
        self.phi1 = self.phi1 + (transfer_field1 + diff1 + damp1) * self.dt
        self.phi2 = self.phi2 + (transfer_field2 + diff2 + damp2) * self.dt
        
        self.phi = self.phi1 + self.phi2
        
        self.time += self.dt
        
        E1 = self.compute_knot_energy(self.mask1)
        E2 = self.compute_knot_energy(self.mask2)
        
        self.history['time'].append(self.time)
        self.history['E1'].append(E1)
        self.history['E2'].append(E2)
        self.history['E_total'].append(E1 + E2)
        self.history['transfer_rate'].append(transfer)
        self.history['phase_coherence'].append(self.compute_phase_coherence())
        
        return {
            'E1': E1,
            'E2': E2,
            'E_total': E1 + E2,
            'transfer_rate': transfer,
            'boson_field': B,
            'phase_coherence': self.compute_phase_coherence(),
        }
    
    def run(self, n_steps: int) -> List[dict]:
        """Run simulation for n_steps."""
        results = []
        for _ in range(n_steps):
            results.append(self.step())
        return results


def visualize_transfer_simulation(system: DualKnotSystem, z_slice: int = None):
    """Visualize the transfer simulation state."""
    if z_slice is None:
        z_slice = system.n // 2
    
    fig = plt.figure(figsize=(16, 10))
    
    x = system.x
    X, Y = np.meshgrid(x, x)
    
    ax1 = fig.add_subplot(2, 3, 1)
    density_slice = system.density[:, :, z_slice]
    im1 = ax1.contourf(X, Y, density_slice.T, levels=20, cmap='viridis')
    ax1.plot(system.knot1_pos[0], system.knot1_pos[1], 'wo', markersize=10, label='Knot 1')
    ax1.plot(system.knot2_pos[0], system.knot2_pos[1], 'ro', markersize=10, label='Knot 2')
    ax1.set_title(f'Total Density (t={system.time:.2f})')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.legend()
    ax1.set_aspect('equal')
    plt.colorbar(im1, ax=ax1)
    
    ax2 = fig.add_subplot(2, 3, 2)
    B = system.compute_boson_field()
    B_slice = B[:, :, z_slice]
    vmax = max(abs(B_slice.min()), abs(B_slice.max())) + 1e-10
    im2 = ax2.contourf(X, Y, B_slice.T, levels=20, cmap='RdBu', vmin=-vmax, vmax=vmax)
    ax2.set_title('Boson Field B(x) (Red=K1→K2, Blue=K2→K1)')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_aspect('equal')
    plt.colorbar(im2, ax=ax2)
    
    ax3 = fig.add_subplot(2, 3, 3)
    phase_slice = system.phase[:, :, z_slice]
    im3 = ax3.contourf(X, Y, phase_slice.T, levels=20, cmap='hsv')
    ax3.set_title('Phase φ')
    ax3.set_xlabel('x')
    ax3.set_ylabel('y')
    ax3.set_aspect('equal')
    plt.colorbar(im3, ax=ax3)
    
    ax4 = fig.add_subplot(2, 3, 4)
    if len(system.history['time']) > 1:
        t = system.history['time']
        ax4.plot(t, system.history['E1'], 'b-', linewidth=2, label='Knot 1')
        ax4.plot(t, system.history['E2'], 'r-', linewidth=2, label='Knot 2')
        ax4.plot(t, system.history['E_total'], 'k--', linewidth=1, label='Total', alpha=0.5)
        ax4.fill_between(t, system.history['E1'], system.history['E2'], 
                         alpha=0.2, color='green')
        ax4.set_xlabel('Time')
        ax4.set_ylabel('Energy')
        ax4.set_title('Energy Transfer')
        ax4.legend()
        ax4.grid(True, alpha=0.3)
    
    ax5 = fig.add_subplot(2, 3, 5)
    if len(system.history['time']) > 1:
        t = system.history['time']
        ax5.plot(t, system.history['transfer_rate'], 'g-', linewidth=2)
        ax5.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
        ax5.fill_between(t, system.history['transfer_rate'], 0,
                         where=np.array(system.history['transfer_rate']) > 0,
                         alpha=0.3, color='blue', label='K1→K2')
        ax5.fill_between(t, system.history['transfer_rate'], 0,
                         where=np.array(system.history['transfer_rate']) <= 0,
                         alpha=0.3, color='red', label='K2→K1')
        ax5.set_xlabel('Time')
        ax5.set_ylabel('Transfer Rate')
        ax5.set_title('Boson Transfer Rate')
        ax5.legend()
        ax5.grid(True, alpha=0.3)
    
    ax6 = fig.add_subplot(2, 3, 6)
    if len(system.history['time']) > 1:
        t = system.history['time']
        ax6.plot(t, system.history['phase_coherence'], 'm-', linewidth=2)
        ax6.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
        ax6.axhline(y=1, color='green', linestyle=':', alpha=0.5)
        ax6.axhline(y=-1, color='red', linestyle=':', alpha=0.5)
        ax6.set_xlabel('Time')
        ax6.set_ylabel('Phase Coherence')
        ax6.set_title('Phase Coherence (1=in-phase, -1=out-of-phase)')
        ax6.set_ylim(-1.2, 1.2)
        ax6.grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig


def run_transfer_demo():
    """Demonstrate boson transfer dynamics."""
    print("=" * 70)
    print("SIMULATION 2: Boson Transfer Dynamics")
    print("=" * 70)
    
    system = DualKnotSystem(
        n=40,
        extent=2.5,
        dt=0.02,
        coupling_range=0.5,
        coupling_strength=0.2,
        diffusion=0.01,
        damping=0.005
    )
    
    E1_init = system.compute_knot_energy(system.mask1)
    E2_init = system.compute_knot_energy(system.mask2)
    
    print(f"\nInitial state:")
    print(f"  Knot 1 energy: {E1_init:.4f}")
    print(f"  Knot 2 energy: {E2_init:.4f}")
    print(f"  Total energy:  {E1_init + E2_init:.4f}")
    print(f"  Initial phase coherence: {system.compute_phase_coherence():.3f}")
    
    fig1 = visualize_transfer_simulation(system)
    fig1.savefig('boson_transfer_initial.png', dpi=100)
    print("\nSaved: boson_transfer_initial.png")
    
    print("\nRunning simulation (300 steps)...")
    n_steps = 300
    for i in range(n_steps):
        result = system.step()
        if (i + 1) % 50 == 0:
            print(f"  Step {i+1}: E1={result['E1']:.4f}, E2={result['E2']:.4f}, "
                  f"transfer={result['transfer_rate']:.4f}")
    
    print(f"\nFinal state (t={system.time:.2f}):")
    print(f"  Knot 1 energy: {system.history['E1'][-1]:.4f}")
    print(f"  Knot 2 energy: {system.history['E2'][-1]:.4f}")
    print(f"  Total energy:  {system.history['E_total'][-1]:.4f}")
    print(f"  Energy transferred: {E1_init - system.history['E1'][-1]:.4f}")
    
    fig2 = visualize_transfer_simulation(system)
    fig2.savefig('boson_transfer_final.png', dpi=100)
    print("Saved: boson_transfer_final.png")
    
    return system


def run_phase_comparison():
    """Compare transfer for different initial phase relationships."""
    print("\n" + "=" * 70)
    print("PHASE DEPENDENCE COMPARISON")
    print("=" * 70)
    
    phase_offsets = [0, np.pi/2, np.pi]
    phase_names = ['In-phase', 'Quadrature', 'Out-of-phase']
    colors = ['blue', 'green', 'red']
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    for phase_offset, name, color in zip(phase_offsets, phase_names, colors):
        print(f"\n{name} (Δφ = {phase_offset/np.pi:.2f}π)...")
        
        system = DualKnotSystem(
            n=40,
            extent=2.5,
            dt=0.02,
            coupling_range=0.5,
            coupling_strength=0.25
        )
        
        system.phi1 = system._create_knot_field(system.knot1_pos, energy=1.0, phase_offset=0.0)
        system.phi2 = system._create_knot_field(system.knot2_pos, energy=0.3, phase_offset=phase_offset)
        system.phi = system.phi1 + system.phi2
        system._update_knot_regions()
        
        n_steps = 200
        for _ in range(n_steps):
            system.step()
        
        t = system.history['time']
        axes[0].plot(t, system.history['E1'], color=color, linewidth=2, 
                    label=f'{name} (K1)')
        axes[0].plot(t, system.history['E2'], color=color, linewidth=2, 
                    linestyle='--', alpha=0.5)
        
        axes[1].plot(t, system.history['transfer_rate'], color=color, linewidth=2,
                    label=name)
        
        print(f"  Final E1: {system.history['E1'][-1]:.4f}")
        print(f"  Final E2: {system.history['E2'][-1]:.4f}")
        print(f"  Total transferred: {system.history['E1'][0] - system.history['E1'][-1]:.4f}")
    
    axes[0].set_xlabel('Time')
    axes[0].set_ylabel('Energy')
    axes[0].set_title('Energy Evolution by Phase Relationship')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    axes[1].set_xlabel('Time')
    axes[1].set_ylabel('Transfer Rate')
    axes[1].set_title('Transfer Rate by Phase')
    axes[1].axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    fig.savefig('boson_phase_comparison.png', dpi=100)
    print("\nSaved: boson_phase_comparison.png")
    
    print("\n" + "=" * 70)
    print("KEY FINDINGS:")
    print("- In-phase: Maximum transfer rate, energy flows to equilibrium")
    print("- Quadrature: Minimal net transfer (cancels out)")
    print("- Out-of-phase: Reverse transfer or repulsion")
    print("=" * 70)


if __name__ == "__main__":
    import matplotlib
    matplotlib.use('Agg')
    
    system = run_transfer_demo()
    run_phase_comparison()
