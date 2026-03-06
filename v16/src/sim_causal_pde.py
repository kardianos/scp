import numpy as np
import os
import json

print("HFKT v14: Avenue 1 - Strict Causal PDE Solver")
print("=============================================")
print("Testing spontaneous resonant condensation under full c-limited propagation.")

os.makedirs("v14/data", exist_ok=True)
out_tsv = "v14/data/sim_causal_pde_raw.tsv"

class CausalPDESim:
    def __init__(self, grid_size=40, dt=0.01, c=1.0, steps=800):
        self.N = grid_size
        self.dt = dt
        self.dx = 1.0
        self.c = c
        self.steps = steps
        
        # We need R(t), R(t-dt), etc., but standard Verlet or Leapfrog is enough.
        self.R = np.zeros((4, self.N, self.N, self.N))
        self.V = np.zeros((4, self.N, self.N, self.N)) # velocity R_dot
        
        self.L = self.N * self.dx
        x = np.linspace(-self.L/2, self.L/2, self.N)
        y = np.linspace(-self.L/2, self.L/2, self.N)
        z = np.linspace(-self.L/2, self.L/2, self.N)
        self.X, self.Y, self.Z = np.meshgrid(x, y, z, indexing='ij')

    def normalize(self, R_arr):
        norm = np.linalg.norm(R_arr, axis=0)
        norm[norm == 0] = 1.0
        return R_arr / norm

    def initialize_asymmetric_pulse(self):
        """ Instantiates a raw, asymmetric, high-energy generic topological pulse.
            Wait to see if the c-limited PDE natively self-organizes it into an echo. """
        print("Injecting generic high-energy asymmetric pulse...")
        
        r3d = np.sqrt(self.X**2 + (self.Y*1.2)**2 + (self.Z*0.8)**2)
        theta = np.pi * np.exp(-r3d / 4.0)
        
        nx = self.X / (r3d + 1e-10)
        ny = self.Y / (r3d + 1e-10)
        nz = self.Z / (r3d + 1e-10)
        
        self.R[0] = np.cos(theta)
        self.R[1] = nx * np.sin(theta)
        self.R[2] = ny * np.sin(theta)
        self.R[3] = nz * np.sin(theta)
        
        self.R = self.normalize(self.R)
        
        # Inject explosive asymmetric radial momentum to trigger wave ripples
        v_mag = 10.0 * np.exp(-r3d / 2.0)
        self.V[1] = v_mag * np.sin(self.Y)
        self.V[2] = v_mag * np.cos(self.X)
        self.V[3] = -v_mag * np.sin(self.Z)
        
        # ensure velocity is tangent to S3
        R_dot_V = np.sum(self.R * self.V, axis=0)
        self.V -= R_dot_V * self.R

    def step(self):
        # Calculate Laplacian
        laplacian = np.zeros_like(self.R)
        for i in range(4):
            # 6-point 3D laplacian
            laplacian[i] = (np.roll(self.R[i], 1, axis=0) + np.roll(self.R[i], -1, axis=0) +
                            np.roll(self.R[i], 1, axis=1) + np.roll(self.R[i], -1, axis=1) +
                            np.roll(self.R[i], 1, axis=2) + np.roll(self.R[i], -1, axis=2) - 
                            6 * self.R[i]) / (self.dx**2)
        
        # The true PDE is R_ddot = c^2 * Laplacian(R) + constraint * R
        # Constraint ensures |R|^2 = 1. 
        # d^2/dt^2 (|R|^2) = 0 => R \cdot R_ddot + |R_dot|^2 = 0
        # R \cdot (c^2 Lap + constraint R) + |V|^2 = 0
        # constraint = -c^2 (R \cdot Lap) - |V|^2
        
        R_dot_lap = np.sum(self.R * laplacian, axis=0)
        v_sq = np.sum(self.V**2, axis=0)
        
        constraint = - (self.c**2)*R_dot_lap - v_sq
        constraint = np.clip(constraint, -100.0, 100.0) # Numerical stability guard
        
        acceleration = (self.c**2)*laplacian + constraint * self.R
        acceleration = np.clip(acceleration, -100.0, 100.0)
        
        self.V += acceleration * self.dt
        self.R += self.V * self.dt
        self.R = self.normalize(self.R)

    def measure_resonant_state(self):
        """ Measures if the field organizes into a stable localized geometry
            or if the topological mass radiates away ($B \to 0$). """
        core_mass = np.sum(self.R[0] < -0.5)
        kinetic_energy = np.sum(self.V**2)
        return int(core_mass), float(kinetic_energy)

    def run(self):
        print(f"Streaming causal PDE data to {out_tsv}...")
        
        with open(out_tsv, "w") as f:
            f.write("Time\tTopo_Core_Mass\tKinetic_Energy\tStatus\n")
            
            for s in range(self.steps):
                self.step()
                t = s * self.dt
                
                if s % 10 == 0:
                    cores, ke = self.measure_resonant_state()
                    status = "Condensing" if cores > 10 else ("Radiating" if cores > 0 else "Dissipated")
                    f.write(f"{t:.3f}\t{cores}\t{ke:.4f}\t{status}\n")
                    
                    if s % 100 == 0:
                        print(f"Step {s:04d} | Cores: {cores} | KE: {ke:.2f} | Status: {status}")
                        
        print("Simulation complete.")

if __name__ == "__main__":
    sim = CausalPDESim()
    sim.initialize_asymmetric_pulse()
    sim.run()
