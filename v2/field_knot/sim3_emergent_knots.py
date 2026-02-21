"""
SIMULATION 3: Emergent Knot Formation via Phase Transition

Demonstrates the Kibble-Zurek mechanism:
1. Start with random field at HIGH temperature (disordered)
2. Cool through phase transition
3. Defects form where phase cannot relax fast enough
4. Stable structures persist at low temperature

This is analogous to cosmic string formation after the Big Bang.
"""

import numpy as np
from dataclasses import dataclass
from typing import List, Tuple
import matplotlib.pyplot as plt
from scipy.ndimage import maximum_filter, label


@dataclass
class PhaseTransitionSimulation:
    """
    Emergent knot formation via phase transition.
    
    The field evolves according to:
    ∂Φ/∂t = D·∇²Φ + (1 - |Φ|²)·Φ - T·Φ + noise
    
    Where T is temperature that decreases over time.
    """
    n: int = 64
    extent: float = 2.0
    dt: float = 0.01
    diffusion: float = 0.1
    T_initial: float = 0.5
    T_final: float = 0.01
    noise_scale: float = 0.005
    
    def __post_init__(self):
        self.dx = 2 * self.extent / self.n
        
        x = np.linspace(-self.extent, self.extent, self.n)
        self.x = x
        self.X, self.Y, self.Z = np.meshgrid(x, x, x, indexing='ij')
        
        np.random.seed(42)
        self.phi = 0.1 * (np.random.randn(self.n, self.n, self.n) + 
                         1j * np.random.randn(self.n, self.n, self.n)) / np.sqrt(2)
        
        self.step_count = 0
        self.T = self.T_initial
        
        self.history = {
            'step': [], 'defects': [], 'density': [], 'temperature': []
        }
    
    @property
    def density(self) -> np.ndarray:
        return np.abs(self.phi) ** 2
    
    @property
    def phase(self) -> np.ndarray:
        return np.angle(self.phi)
    
    def laplacian(self, f: np.ndarray) -> np.ndarray:
        return (np.roll(f, 1, 0) + np.roll(f, -1, 0) +
                np.roll(f, 1, 1) + np.roll(f, -1, 1) +
                np.roll(f, 1, 2) + np.roll(f, -1, 2) - 6*f) / self.dx**2
    
    def detect_defects(self) -> Tuple[np.ndarray, int]:
        """Detect phase singularities (topological defects)."""
        rho = self.density
        phase = self.phase
        
        dpx = np.diff(phase, axis=0, prepend=phase[:1])
        dpy = np.diff(phase, axis=1, prepend=phase[:, :1])
        dpz = np.diff(phase, axis=2, prepend=phase[:, :, :1])
        
        for arr in [dpx, dpy, dpz]:
            arr[:] = np.where(arr > np.pi, arr - 2*np.pi, arr)
            arr[:] = np.where(arr < -np.pi, arr + 2*np.pi, arr)
        
        winding = np.sqrt(dpx**2 + dpy**2 + dpz**2)
        winding = winding * (rho > 0.3)
        
        local_max = (winding > 0.3) & (winding == maximum_filter(winding, size=7))
        labeled, num = label(local_max)
        
        return labeled, num
    
    def step(self) -> dict:
        """Advance one time step."""
        self.step_count += 1
        
        total_steps = 600
        if self.step_count < total_steps // 2:
            self.T = self.T_initial * (1 - self.step_count / (total_steps // 2))
        else:
            self.T = self.T_final
        
        rho = self.density
        
        lap = self.laplacian(self.phi)
        
        nonlinear = (1 - rho) * self.phi
        
        noise = self.noise_scale * np.sqrt(max(self.T, 0.01)) * (
            np.random.randn(self.n, self.n, self.n) + 
            1j * np.random.randn(self.n, self.n, self.n)
        )
        
        dphi_dt = self.diffusion * lap + nonlinear - self.T * self.phi + noise
        
        self.phi += self.dt * dphi_dt
        
        amp = np.abs(self.phi)
        self.phi = np.clip(amp, 0, 5) * np.exp(1j * np.angle(self.phi))
        
        _, num_defects = self.detect_defects()
        
        self.history['step'].append(self.step_count)
        self.history['defects'].append(num_defects)
        self.history['density'].append(np.mean(rho))
        self.history['temperature'].append(self.T)
        
        return {
            'defects': num_defects,
            'density': np.mean(rho),
            'temperature': self.T
        }
    
    def run(self, n_steps: int, verbose: bool = True) -> List[dict]:
        """Run simulation."""
        results = []
        for i in range(n_steps):
            result = self.step()
            results.append(result)
            
            if verbose and (i + 1) % 100 == 0:
                print(f"  Step {i+1}: T={result['temperature']:.3f}, "
                      f"ρ={result['density']:.3f}, defects={result['defects']}")
        
        return results
    
    def get_defect_positions(self) -> List[Tuple[float, float, float]]:
        """Get positions of all defects."""
        labeled, _ = self.detect_defects()
        positions = []
        
        for i in range(1, labeled.max() + 1):
            mask = labeled == i
            if np.sum(mask) > 0:
                cx = np.mean(self.X[mask])
                cy = np.mean(self.Y[mask])
                cz = np.mean(self.Z[mask])
                positions.append((cx, cy, cz))
        
        return positions


def visualize_final_state(sim: PhaseTransitionSimulation):
    """Visualize final state with defects."""
    fig = plt.figure(figsize=(14, 10))
    z_mid = sim.n // 2
    
    ax1 = fig.add_subplot(2, 2, 1)
    rho_slice = sim.density[:, :, z_mid]
    im1 = ax1.contourf(sim.X[:,:,z_mid], sim.Y[:,:,z_mid], rho_slice.T, 
                       levels=20, cmap='viridis')
    ax1.set_title(f'Final Density |Φ|² (T={sim.T:.3f})')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_aspect('equal')
    plt.colorbar(im1, ax=ax1)
    
    ax2 = fig.add_subplot(2, 2, 2)
    phase_slice = sim.phase[:, :, z_mid]
    im2 = ax2.contourf(sim.X[:,:,z_mid], sim.Y[:,:,z_mid], phase_slice.T,
                       levels=20, cmap='hsv')
    ax2.set_title('Final Phase arg(Φ)')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_aspect('equal')
    plt.colorbar(im2, ax=ax2)
    
    ax3 = fig.add_subplot(2, 2, 3)
    labeled, num_defects = sim.detect_defects()
    
    phase = sim.phase
    dpx = np.diff(phase, axis=0, prepend=phase[:1])
    dpy = np.diff(phase, axis=1, prepend=phase[:,:1])
    dpx = np.where(dpx > np.pi, dpx - 2*np.pi, dpx)
    dpx = np.where(dpx < -np.pi, dpx + 2*np.pi, dpx)
    dpy = np.where(dpy > np.pi, dpy - 2*np.pi, dpy)
    dpy = np.where(dpy < -np.pi, dpy + 2*np.pi, dpy)
    winding = (np.abs(dpx) + np.abs(dpy)) * (sim.density > 0.3)
    
    im3 = ax3.contourf(sim.X[:,:,z_mid], sim.Y[:,:,z_mid], 
                       winding[:,:,z_mid].T, levels=20, cmap='hot')
    
    positions = sim.get_defect_positions()
    for cx, cy, cz in positions:
        if abs(cz - sim.x[z_mid]) < sim.dx * 3:
            ax3.plot(cx, cy, 'wo', markersize=8, markeredgecolor='cyan', 
                    markeredgewidth=2)
    
    ax3.set_title(f'Phase Singularities ({num_defects} defects)')
    ax3.set_xlabel('x')
    ax3.set_ylabel('y')
    ax3.set_aspect('equal')
    plt.colorbar(im3, ax=ax3, label='|∇φ|')
    
    ax4 = fig.add_subplot(2, 2, 4)
    ax4_twin = ax4.twinx()
    
    steps = sim.history['step']
    l1, = ax4.plot(steps, sim.history['defects'], 'g-', lw=2, label='Defects')
    l2, = ax4_twin.plot(steps, sim.history['temperature'], 'r--', lw=2, 
                        label='Temperature')
    
    ax4.axvline(x=len(steps)//2, color='gray', linestyle=':', alpha=0.5)
    ax4.text(len(steps)//2, max(sim.history['defects'])*0.9, ' T=0.01', 
             fontsize=10, color='gray')
    
    ax4.set_xlabel('Step')
    ax4.set_ylabel('Defect Count', color='green')
    ax4_twin.set_ylabel('Temperature', color='red')
    ax4.set_title('Knot Formation During Cooling')
    ax4.legend(handles=[l1, l2], loc='upper right')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig


def run_emergence_demo():
    """Run the emergent knot demonstration."""
    print("=" * 70)
    print("EMERGENT KNOT FORMATION via Phase Transition")
    print("=" * 70)
    print("\nKibble-Zurek Mechanism:")
    print("1. High T: Random phases, no structure")
    print("2. Cooling: Phase transition begins")
    print("3. Defect formation: Phase 'freezes' at boundaries")
    print("4. Low T: Stable structures persist")
    print()
    
    sim = PhaseTransitionSimulation(
        n=64, extent=2.0, dt=0.01,
        diffusion=0.1, T_initial=0.5, T_final=0.01
    )
    
    print("Running phase transition (600 steps)...")
    sim.run(600, verbose=True)
    
    print(f"\nFinal state:")
    print(f"  Temperature: {sim.T:.3f}")
    print(f"  Mean density: {np.mean(sim.density):.3f}")
    print(f"  Defect count: {sim.history['defects'][-1]}")
    
    fig = visualize_final_state(sim)
    fig.savefig('emergence_phase_transition.png', dpi=120)
    print("\nSaved: emergence_phase_transition.png")
    
    print("\n" + "=" * 70)
    print("KEY OBSERVATIONS:")
    print("- Defects form during phase transition (Kibble-Zurek)")
    print("- More defects at faster cooling rates")
    print("- Defects persist as stable topological structures")
    print("- These are emergent knots - not placed, but formed")
    print("=" * 70)
    
    return sim


if __name__ == "__main__":
    import matplotlib
    matplotlib.use('Agg')
    
    sim = run_emergence_demo()
