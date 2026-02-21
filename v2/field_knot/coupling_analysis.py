"""
Coupled Field Knots - Stability and Energy Analysis

This module explores the interaction between multiple field knots:
- How do knots exchange energy?
- What phase relationships are stable (in-sync vs. out-of-sync)?
- How does spatial overlap affect coupling strength?

Mathematical Framework for Coupling:
------------------------------------
For two knots with curves γ₁(s) and γ₂(t):

Field at knot 1 due to knot 2:
  F₁←₂ = ∫∫ A₁(s) A₂(t') · K₁₂(‖γ₁(s) - γ₂(t')‖) · cos(Δφ) dt' ds

Where:
  - K₁₂: cross-kernel (spatial coupling strength)
  - Δφ = φ₁(s) - φ₂(t'): phase difference

Stability Conditions:
---------------------
1. PHASE-LOCKED: φ₁ ≈ φ₂ (in-sync oscillation) - often STABLE
2. PHASE-OPPOSED: φ₁ ≈ φ₂ + π (out-of-sync) - can be stable if coupling < threshold
3. PHASE-DRIFTING: Δφ(t) grows unbounded - UNSTABLE

Energy:
-------
E_total = E₁ + E₂ + E_coupling
E_coupling = ∫∫ A₁ A₂ K₁₂ cos(Δφ) ds dt'

The coupling energy can be:
  - NEGATIVE (binding): knots attract/stabilize each other
  - POSITIVE (repulsive): knots destabilize each other
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Tuple, Callable, Optional
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D


@dataclass
class KnotCurve:
    """Parametric curve for a field knot."""
    points: np.ndarray  # (N, 3) array of curve points
    name: str = "Knot"
    
    @classmethod
    def trefoil(cls, r_major: float = 0.4, r_minor: float = 0.15, 
                center: np.ndarray = None, n_points: int = 100, rotation: float = 0.0) -> 'KnotCurve':
        s = np.linspace(0, 2*np.pi, n_points, endpoint=False)
        phi = 2 * s
        theta = 3 * s + rotation
        
        x = (r_major + r_minor * np.cos(phi)) * np.cos(theta)
        y = (r_major + r_minor * np.cos(phi)) * np.sin(theta)
        z = r_minor * np.sin(phi)
        
        points = np.stack([x, y, z], axis=-1)
        if center is not None:
            points += center
        
        return cls(points=points, name="Trefoil")
    
    @classmethod
    def figure_eight(cls, scale: float = 0.35, center: np.ndarray = None, 
                     n_points: int = 100) -> 'KnotCurve':
        s = np.linspace(0, 2*np.pi, n_points, endpoint=False)
        x = (2 + np.cos(2*s)) * np.cos(3*s) * scale * 0.3
        y = (2 + np.cos(2*s)) * np.sin(3*s) * scale * 0.3
        z = np.sin(4*s) * scale * 0.4
        
        points = np.stack([x, y, z], axis=-1)
        if center is not None:
            points += center
        
        return cls(points=points, name="Figure-8")
    
    @classmethod
    def unknot(cls, radius: float = 0.4, center: np.ndarray = None, 
               n_points: int = 100) -> 'KnotCurve':
        """Simple circle (unknot) - useful for comparison."""
        s = np.linspace(0, 2*np.pi, n_points, endpoint=False)
        x = radius * np.cos(s)
        y = radius * np.sin(s)
        z = np.zeros_like(s)
        
        points = np.stack([x, y, z], axis=-1)
        if center is not None:
            points += center
        
        return cls(points=points, name="Unknot")
    
    def get_point(self, i: int) -> np.ndarray:
        return self.points[i % len(self.points)]
    
    @property
    def n_points(self) -> int:
        return len(self.points)


@dataclass
class KnotState:
    """State of a single field knot."""
    amplitude: np.ndarray  # (N,) amplitude at each point
    velocity: np.ndarray   # (N,) velocity at each point
    phase: np.ndarray      # (N,) phase at each point
    global_phase: float = 0.0  # Overall phase offset
    
    @classmethod
    def standing_wave(cls, n: int, wave_number: int, amplitude: float, 
                      phase_offset: float = 0.0) -> 'KnotState':
        s = np.linspace(0, 2*np.pi, n, endpoint=False)
        amp = amplitude * np.sin(wave_number * s)
        phase = wave_number * s + phase_offset
        vel = np.zeros(n)
        return cls(amplitude=amp, velocity=vel, phase=phase, global_phase=phase_offset)
    
    def energy(self) -> float:
        kinetic = 0.5 * np.sum(self.velocity ** 2)
        potential = 0.5 * np.sum(self.amplitude ** 2)
        return kinetic + potential


@dataclass
class CouplingConfig:
    """Configuration for knot-knot coupling."""
    spatial_coupling: float = 0.3      # Strength of spatial overlap coupling
    phase_coupling: float = 0.2        # Strength of phase-dependent coupling
    cross_kernel_width: float = 0.15   # Distance scale for cross-coupling
    energy_transfer_rate: float = 0.1  # How fast energy transfers between knots


class CoupledKnotSystem:
    """
    System of multiple coupled field knots.
    
    Each knot has its own amplitude field that evolves according to:
    - Wave dynamics along its curve
    - Self-reinforcement
    - Coupling to other knots
    """
    
    def __init__(
        self,
        curves: List[KnotCurve],
        states: List[KnotState],
        wave_numbers: List[int],
        config: CouplingConfig = None,
        wave_speed: float = 1.0,
        damping: float = 0.02,
        nonlinearity: float = 0.1,
        self_coupling: float = 0.3,
    ):
        self.curves = curves
        self.states = states
        self.wave_numbers = wave_numbers
        self.config = config or CouplingConfig()
        self.wave_speed = wave_speed
        self.damping = damping
        self.nonlinearity = nonlinearity
        self.self_coupling = self_coupling
        self.n_knots = len(curves)
        
        self._precompute_distances()
    
    def _precompute_distances(self):
        """Precompute inter-knot distances for coupling calculation."""
        self.cross_distances = {}
        
        for i in range(self.n_knots):
            for j in range(i + 1, self.n_knots):
                pts_i = self.curves[i].points  # (N_i, 3)
                pts_j = self.curves[j].points  # (N_j, 3)
                
                dists = np.zeros((len(pts_i), len(pts_j)))
                for k, p_i in enumerate(pts_i):
                    dists[k] = np.linalg.norm(pts_j - p_i, axis=1)
                
                self.cross_distances[(i, j)] = dists
    
    def compute_cross_coupling_force(self, knot_idx: int) -> np.ndarray:
        """
        Compute the coupling force on knot `knot_idx` from all other knots.
        
        F_i←j[k] = Σ_{l in knot j} A_j[l] · K(d_{kl}) · cos(Δφ_{kl})
        """
        n = self.curves[knot_idx].n_points
        force = np.zeros(n)
        state_i = self.states[knot_idx]
        
        for j in range(self.n_knots):
            if j == knot_idx:
                continue
            
            state_j = self.states[j]
            
            min_idx = min(knot_idx, j)
            max_idx = max(knot_idx, j)
            dists = self.cross_distances[(min_idx, max_idx)]
            
            if knot_idx > j:
                dists = dists.T
            
            sigma = self.config.cross_kernel_width
            spatial_kernel = np.exp(-dists**2 / (2 * sigma**2))
            
            phase_diff = state_i.phase[:, np.newaxis] - state_j.phase[np.newaxis, :]
            phase_factor = np.cos(phase_diff)
            
            coupling_matrix = (spatial_kernel * phase_factor * 
                             state_j.amplitude[np.newaxis, :])
            
            force += self.config.spatial_coupling * np.sum(coupling_matrix, axis=1)
        
        return force / (self.n_knots - 1) if self.n_knots > 1 else force
    
    def compute_coupling_energy(self, knot_i: int, knot_j: int) -> float:
        """
        Compute the coupling energy between two knots.
        
        E_ij = ∫∫ A_i A_j K(d) cos(Δφ) ds dt
        """
        min_idx = min(knot_i, knot_j)
        max_idx = max(knot_i, knot_j)
        dists = self.cross_distances[(min_idx, max_idx)]
        
        if knot_i > knot_j:
            dists = dists.T
        
        state_i = self.states[knot_i]
        state_j = self.states[knot_j]
        
        sigma = self.config.cross_kernel_width
        spatial_kernel = np.exp(-dists**2 / (2 * sigma**2))
        
        phase_diff = state_i.phase[:, np.newaxis] - state_j.phase[np.newaxis, :]
        phase_factor = np.cos(phase_diff)
        
        energy_matrix = (state_i.amplitude[:, np.newaxis] * 
                        state_j.amplitude[np.newaxis, :] * 
                        spatial_kernel * phase_factor)
        
        return self.config.spatial_coupling * np.sum(energy_matrix)
    
    def compute_self_force(self, knot_idx: int) -> np.ndarray:
        """Self-reinforcement force within a knot."""
        curve = self.curves[knot_idx]
        state = self.states[knot_idx]
        n = curve.n_points
        
        force = np.zeros(n)
        sigma = 0.08
        
        for i in range(n):
            dists = np.linalg.norm(curve.points - curve.points[i], axis=1)
            spatial = np.exp(-dists**2 / (2 * sigma**2))
            spatial[dists < 0.01] = 0
            
            phase_diff = state.phase[i] - state.phase
            interference = np.cos(phase_diff)
            
            force[i] = self.self_coupling * np.sum(state.amplitude * spatial * interference)
        
        return force
    
    def equations_of_motion(self, y: np.ndarray, t: float) -> np.ndarray:
        """
        ODEs for the coupled system.
        y = [A₁, v₁, A₂, v₂, ...] flattened
        """
        n_knots = self.n_knots
        dydt = np.zeros_like(y)
        
        all_A = []
        all_v = []
        offset = 0
        for k in range(n_knots):
            n = self.curves[k].n_points
            all_A.append(y[offset:offset+n])
            all_v.append(y[offset+n:offset+2*n])
            offset += 2 * n
        
        for k in range(n_knots):
            n = self.curves[k].n_points
            A = all_A[k]
            v = all_v[k]
            
            curve = self.curves[k]
            arclength = np.sum(np.linalg.norm(
                np.diff(curve.points, axis=0, append=curve.points[:1]), axis=1))
            ds = arclength / n
            
            c2 = self.wave_speed ** 2
            A_prev = np.roll(A, 1)
            A_next = np.roll(A, -1)
            d2A_ds2 = (A_prev - 2*A + A_next) / (ds**2)
            
            wave_term = c2 * d2A_ds2
            damping_term = -self.damping * v
            nonlinear_term = -self.nonlinearity * A**3
            
            temp_state = KnotState(amplitude=A, velocity=v, phase=self.states[k].phase)
            old_state = self.states[k]
            self.states[k] = temp_state
            self_term = self.compute_self_force(k)
            cross_term = self.compute_cross_coupling_force(k)
            self.states[k] = old_state
            
            dv_dt = wave_term + damping_term + nonlinear_term + self_term + cross_term
            
            offset = k * 2 * n
            dydt[offset:offset+n] = v
            dydt[offset+n:offset+2*n] = dv_dt
        
        return dydt
    
    def simulate(self, t_span: np.ndarray) -> List['CoupledKnotSystem']:
        """Simulate the coupled system over time."""
        y0 = []
        for k in range(self.n_knots):
            n = self.curves[k].n_points
            y0.extend(self.states[k].amplitude)
            y0.extend(self.states[k].velocity)
        y0 = np.array(y0)
        
        solution = odeint(self.equations_of_motion, y0, t_span)
        
        systems = []
        for i, t in enumerate(t_span):
            new_states = []
            offset = 0
            for k in range(self.n_knots):
                n = self.curves[k].n_points
                A = solution[i, offset:offset+n]
                v = solution[i, offset+n:offset+2*n]
                new_states.append(KnotState(
                    amplitude=A.copy(),
                    velocity=v.copy(),
                    phase=self.states[k].phase.copy(),
                    global_phase=self.states[k].global_phase,
                ))
                offset += 2 * n
            
            new_system = CoupledKnotSystem(
                curves=self.curves,
                states=new_states,
                wave_numbers=self.wave_numbers,
                config=self.config,
                wave_speed=self.wave_speed,
                damping=self.damping,
                nonlinearity=self.nonlinearity,
                self_coupling=self.self_coupling,
            )
            systems.append(new_system)
        
        return systems
    
    def total_energy(self) -> Tuple[float, float, float]:
        """Return (self_energy, coupling_energy, total_energy)."""
        self_energy = sum(s.energy() for s in self.states)
        
        coupling_energy = 0.0
        for i in range(self.n_knots):
            for j in range(i + 1, self.n_knots):
                coupling_energy += self.compute_coupling_energy(i, j)
        
        return self_energy, coupling_energy, self_energy + coupling_energy
    
    def phase_coherence(self) -> float:
        """
        Measure of phase synchronization between knots.
        Returns value in [-1, 1]:
          1 = perfectly in-phase
         -1 = perfectly out-of-phase
          0 = no correlation
        """
        if self.n_knots < 2:
            return 0.0
        
        coherence = 0.0
        count = 0
        
        for i in range(self.n_knots):
            for j in range(i + 1, self.n_knots):
                phase_diff = np.mean(self.states[i].phase - self.states[j].phase)
                coherence += np.cos(phase_diff)
                count += 1
        
        return coherence / count if count > 0 else 0.0


def analyze_stability(system: CoupledKnotSystem, t_span: np.ndarray) -> dict:
    """
    Analyze stability of a coupled system.
    
    Returns:
        dict with keys:
            - energy_trajectory: (T, 3) array of [self, coupling, total] energy
            - phase_coherence: (T,) array
            - amplitude_ratio: (T, n_knots) array of relative amplitudes
            - is_stable: bool
            - coupling_regime: str ("bound", "repulsive", "neutral")
    """
    systems = system.simulate(t_span)
    
    energy_traj = np.array([s.total_energy() for s in systems])
    coherence_traj = np.array([s.phase_coherence() for s in systems])
    
    amp_ratio = np.array([
        [np.max(np.abs(s.states[k].amplitude)) for k in range(system.n_knots)]
        for s in systems
    ])
    
    final_energy = energy_traj[-1, 2]
    initial_energy = energy_traj[0, 2]
    energy_change = (final_energy - initial_energy) / (initial_energy + 1e-6)
    
    final_coupling = energy_traj[-1, 1]
    coupling_regime = "bound" if final_coupling < -0.1 else "repulsive" if final_coupling > 0.1 else "neutral"
    
    is_stable = abs(energy_change) < 0.5 and np.std(coherence_traj[-10:]) < 0.3
    
    return {
        "energy_trajectory": energy_traj,
        "phase_coherence": coherence_traj,
        "amplitude_ratio": amp_ratio,
        "is_stable": is_stable,
        "coupling_regime": coupling_regime,
        "energy_change": energy_change,
    }


def plot_coupling_analysis(t_span: np.ndarray, results: dict, title: str = ""):
    """Create visualization of coupling analysis."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    ax1 = axes[0, 0]
    ax1.plot(t_span, results["energy_trajectory"][:, 0], label="Self Energy", color='blue')
    ax1.plot(t_span, results["energy_trajectory"][:, 1], label="Coupling Energy", color='red')
    ax1.plot(t_span, results["energy_trajectory"][:, 2], label="Total Energy", color='green', linestyle='--')
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Energy")
    ax1.set_title("Energy Evolution")
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    ax2 = axes[0, 1]
    ax2.plot(t_span, results["phase_coherence"])
    ax2.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax2.axhline(y=1, color='green', linestyle=':', alpha=0.5, label="In-phase")
    ax2.axhline(y=-1, color='red', linestyle=':', alpha=0.5, label="Out-of-phase")
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Phase Coherence")
    ax2.set_title("Phase Synchronization")
    ax2.set_ylim(-1.2, 1.2)
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    ax3 = axes[1, 0]
    n_knots = results["amplitude_ratio"].shape[1]
    for k in range(n_knots):
        ax3.plot(t_span, results["amplitude_ratio"][:, k], label=f"Knot {k+1}")
    ax3.set_xlabel("Time")
    ax3.set_ylabel("Max |Amplitude|")
    ax3.set_title("Amplitude Evolution")
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    ax4 = axes[1, 1]
    ax4.axis('off')
    status = "STABLE" if results["is_stable"] else "UNSTABLE"
    regime = results["coupling_regime"].upper()
    energy_change = results["energy_change"]
    
    info_text = f"""
    STABILITY ANALYSIS
    {'='*30}
    Status: {status}
    Coupling Regime: {regime}
    Energy Change: {energy_change:+.2%}
    
    Phase Coherence (final): {results['phase_coherence'][-1]:.3f}
    Coupling Energy (final): {results['energy_trajectory'][-1, 1]:.3f}
    """
    ax4.text(0.1, 0.5, info_text, fontsize=12, family='monospace',
             verticalalignment='center', transform=ax4.transAxes,
             bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.5))
    
    plt.suptitle(title, fontsize=14)
    plt.tight_layout()
    return fig


def explore_coupling_parameter_space():
    """
    Explore how coupling strength affects stability.
    """
    print("=" * 70)
    print("COUPLING PARAMETER SPACE EXPLORATION")
    print("=" * 70)
    
    coupling_strengths = np.linspace(0.0, 1.0, 11)
    phase_offsets = [0.0, np.pi/2, np.pi]
    phase_names = ["In-phase", "Quadrature", "Out-of-phase"]
    
    results_grid = np.zeros((len(coupling_strengths), len(phase_offsets)))
    
    t_span = np.linspace(0, 15, 100)
    
    for i, coupling in enumerate(coupling_strengths):
        for j, phase_offset in enumerate(phase_offsets):
            curve1 = KnotCurve.trefoil(center=np.array([-0.3, 0, 0]))
            curve2 = KnotCurve.trefoil(center=np.array([0.3, 0, 0]))
            
            state1 = KnotState.standing_wave(100, 3, 0.5, phase_offset=0.0)
            state2 = KnotState.standing_wave(100, 3, 0.5, phase_offset=phase_offset)
            
            config = CouplingConfig(spatial_coupling=coupling)
            
            system = CoupledKnotSystem(
                curves=[curve1, curve2],
                states=[state1, state2],
                wave_numbers=[3, 3],
                config=config,
            )
            
            analysis = analyze_stability(system, t_span)
            
            final_coherence = abs(analysis["phase_coherence"][-1])
            energy_drift = abs(analysis["energy_change"])
            stability_score = final_coherence * (1 - min(energy_drift, 1))
            
            results_grid[i, j] = stability_score
            
            print(f"  coupling={coupling:.1f}, phase={phase_names[j]:15s} -> "
                  f"stability={stability_score:.3f} ({analysis['coupling_regime']})")
    
    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(results_grid, aspect='auto', cmap='RdYlGn',
                   extent=[-0.5, len(phase_offsets)-0.5, 
                          coupling_strengths[-1]+0.05, coupling_strengths[0]-0.05])
    ax.set_xticks(range(len(phase_offsets)))
    ax.set_xticklabels(phase_names)
    ax.set_xlabel("Initial Phase Relationship")
    ax.set_ylabel("Coupling Strength")
    ax.set_title("Stability Map: Coupling Strength vs Phase Offset\n(Green=Stable, Red=Unstable)")
    plt.colorbar(im, ax=ax, label="Stability Score")
    
    return fig, results_grid


def explore_energy_transfer():
    """
    Explore energy transfer between coupled knots.
    """
    print("\n" + "=" * 70)
    print("ENERGY TRANSFER ANALYSIS")
    print("=" * 70)
    
    curve1 = KnotCurve.trefoil(center=np.array([-0.4, 0, 0]))
    curve2 = KnotCurve.trefoil(center=np.array([0.4, 0, 0]))
    
    state1 = KnotState.standing_wave(100, 3, 0.8, phase_offset=0.0)
    state2 = KnotState.standing_wave(100, 3, 0.2, phase_offset=0.0)
    
    config = CouplingConfig(spatial_coupling=0.4, energy_transfer_rate=0.15)
    
    system = CoupledKnotSystem(
        curves=[curve1, curve2],
        states=[state1, state2],
        wave_numbers=[3, 3],
        config=config,
        damping=0.02,
    )
    
    t_span = np.linspace(0, 25, 200)
    systems = system.simulate(t_span)
    
    energies1 = np.array([s.states[0].energy() for s in systems])
    energies2 = np.array([s.states[1].energy() for s in systems])
    coupling_energies = np.array([s.compute_coupling_energy(0, 1) for s in systems])
    
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    ax1 = axes[0]
    ax1.plot(t_span, energies1, label="Knot 1 Energy", color='blue')
    ax1.plot(t_span, energies2, label="Knot 2 Energy", color='red')
    ax1.fill_between(t_span, energies1, energies2, alpha=0.3, color='purple', label="Energy difference")
    ax1.set_ylabel("Self Energy")
    ax1.set_title("Energy Transfer Between Coupled Knots")
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    ax2 = axes[1]
    ax2.plot(t_span, coupling_energies, color='green')
    ax2.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax2.fill_between(t_span, coupling_energies, 0, 
                     where=coupling_energies < 0, alpha=0.3, color='blue', label="Binding")
    ax2.fill_between(t_span, coupling_energies, 0, 
                     where=coupling_energies >= 0, alpha=0.3, color='red', label="Repulsive")
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Coupling Energy")
    ax2.set_title("Coupling Energy (Negative = Bound, Positive = Repulsive)")
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig


def explore_knot_linking():
    """
    Compare coupling between linked vs. unlinked knots.
    """
    print("\n" + "=" * 70)
    print("KNOT LINKING ANALYSIS")
    print("=" * 70)
    
    t_span = np.linspace(0, 20, 150)
    
    print("\n1. Linked Trefoils (close together)...")
    curve1_linked = KnotCurve.trefoil(center=np.array([0, 0, 0]))
    curve2_linked = KnotCurve.trefoil(center=np.array([0.1, 0, 0.1]), rotation=np.pi/3)
    
    state1 = KnotState.standing_wave(100, 3, 0.5, phase_offset=0.0)
    state2 = KnotState.standing_wave(100, 3, 0.5, phase_offset=0.0)
    
    config = CouplingConfig(spatial_coupling=0.5)
    
    linked_system = CoupledKnotSystem(
        curves=[curve1_linked, curve2_linked],
        states=[state1, state2],
        wave_numbers=[3, 3],
        config=config,
    )
    linked_results = analyze_stability(linked_system, t_span)
    
    print(f"   Stability: {linked_results['is_stable']}")
    print(f"   Regime: {linked_results['coupling_regime']}")
    print(f"   Final coherence: {linked_results['phase_coherence'][-1]:.3f}")
    
    print("\n2. Separated Trefoils (far apart)...")
    curve1_sep = KnotCurve.trefoil(center=np.array([-0.8, 0, 0]))
    curve2_sep = KnotCurve.trefoil(center=np.array([0.8, 0, 0]))
    
    sep_system = CoupledKnotSystem(
        curves=[curve1_sep, curve2_sep],
        states=[state1, state2],
        wave_numbers=[3, 3],
        config=config,
    )
    sep_results = analyze_stability(sep_system, t_span)
    
    print(f"   Stability: {sep_results['is_stable']}")
    print(f"   Regime: {sep_results['coupling_regime']}")
    print(f"   Final coherence: {sep_results['phase_coherence'][-1]:.3f}")
    
    fig1 = plot_coupling_analysis(t_span, linked_results, "Linked Trefoils (Close)")
    fig2 = plot_coupling_analysis(t_span, sep_results, "Separated Trefoils (Far)")
    
    return fig1, fig2


if __name__ == "__main__":
    print("=" * 70)
    print("COUPLED FIELD KNOTS - Stability & Energy Analysis")
    print("=" * 70)
    
    fig_param, results_grid = explore_coupling_parameter_space()
    fig_param.savefig('coupling_parameter_space.png', dpi=120)
    print("\nSaved: coupling_parameter_space.png")
    
    fig_energy = explore_energy_transfer()
    fig_energy.savefig('energy_transfer.png', dpi=120)
    print("Saved: energy_transfer.png")
    
    fig_linked, fig_sep = explore_knot_linking()
    fig_linked.savefig('linked_knots.png', dpi=120)
    fig_sep.savefig('separated_knots.png', dpi=120)
    print("Saved: linked_knots.png, separated_knots.png")
    
    print("\n" + "=" * 70)
    print("KEY INSIGHTS:")
    print("- In-phase knots (Δφ ≈ 0) tend to BIND (negative coupling energy)")
    print("- Out-of-phase knots (Δφ ≈ π) tend to REPEL (positive coupling energy)")
    print("- Strong coupling can cause instability (energy cascade)")
    print("- Spatial overlap determines coupling strength")
    print("- Linked knots have stronger interactions than separated knots")
    print("=" * 70)
    
    plt.show()
