import numpy as np
import os

print("HFKT v13: Avenue 4 - Rational Map Polynomial Evolution Simulator")
print("==================================================================")
print("Testing if a monolithic 3-lobed algebraic topology can confine quarks.")

os.makedirs("v13/results", exist_ok=True)
out_tsv = "v13/results/sim_rational_map_evolution_raw.tsv"

class RationalMapSim:
    def __init__(self, grid_size=40, dt=0.01, steps=800):
        self.N = grid_size
        self.dt = dt
        self.steps = steps
        
        self.R = np.zeros((4, self.N, self.N, self.N))
        self.R_dot = np.zeros((4, self.N, self.N, self.N))
        
        self.L = 15.0
        x = np.linspace(-self.L/2, self.L/2, self.N)
        y = np.linspace(-self.L/2, self.L/2, self.N)
        z = np.linspace(-self.L/2, self.L/2, self.N)
        self.X, self.Y, self.Z = np.meshgrid(x, y, z, indexing='ij')

    def normalize(self):
        norm = np.linalg.norm(self.R, axis=0)
        norm[norm == 0] = 1.0
        self.R /= norm

    def initialize_rational_map(self):
        print("Constructing W(z) = z^3 - a ...")
        
        r2d = np.sqrt(self.X**2 + self.Y**2)
        theta = np.arctan2(self.Y, self.X)
        
        a = 2.0  # Separation parameter
        W_r = r2d**3 * np.cos(3*theta) - a
        W_i = r2d**3 * np.sin(3*theta)
        
        W_mag2 = W_r**2 + W_i**2
        den = 1.0 + W_mag2
        den[den==0] = 1e-10
        
        nx = 2 * W_r / den
        ny = 2 * W_i / den
        nz = (1.0 - W_mag2) / den
        
        r3d = np.sqrt(self.X**2 + self.Y**2 + self.Z**2)
        f = np.pi * np.exp(-r3d / 3.0)
        
        self.R[0] = np.cos(f)
        self.R[1] = nx * np.sin(f)
        self.R[2] = ny * np.sin(f)
        self.R[3] = nz * np.sin(f)
        
        # Perturb it slightly to test stability
        self.R_dot[1] += 0.05 * np.exp(-r3d/2.0) * np.sin(3*theta)
        
        self.normalize()

    def step(self):
        laplacian = np.zeros_like(self.R)
        for i in range(4):
            laplacian[i] = (np.roll(self.R[i], 1, axis=0) + np.roll(self.R[i], -1, axis=0) +
                            np.roll(self.R[i], 1, axis=1) + np.roll(self.R[i], -1, axis=1) +
                            np.roll(self.R[i], 1, axis=2) + np.roll(self.R[i], -1, axis=2) - 
                            6 * self.R[i])
        
        R_dot_sq = np.sum(self.R_dot**2, axis=0)
        R_dot_lap = np.sum(self.R * laplacian, axis=0)
        
        constraint_mag = -(R_dot_sq + R_dot_lap)
        constraint_mag = np.clip(constraint_mag, -50.0, 50.0)
        
        acceleration = laplacian + constraint_mag * self.R
        acceleration = np.clip(acceleration, -50.0, 50.0)
        
        self.R_dot += acceleration * self.dt
        self.R += self.R_dot * self.dt
        self.normalize()

    def measure_structural_integrity(self):
        """ Proxy for B-charge topological stability over time """
        # A true topological state will have a non-zero core where R[0] < -0.5
        core_mass = np.sum(self.R[0] < -0.5)
        
        # Track the expansion/dispersion of the kinetic energy
        kinetic_energy = np.sum(self.R_dot**2)
        
        return int(core_mass), float(kinetic_energy)

    def run(self):
        print(f"Streaming topological integrity data to {out_tsv}...")
        
        with open(out_tsv, "w") as f:
            f.write("Time\tCore_Volume\tTotal_Kinetic_Energy\tStatus\n")
            
            for s in range(self.steps):
                self.step()
                t = s * self.dt
                
                if s % 10 == 0:
                    cores, ke = self.measure_structural_integrity()
                    status = "Bound" if cores > 10 else ("Dissolving" if cores > 0 else "Annihilated")
                    f.write(f"{t:.3f}\t{cores}\t{ke:.4f}\t{status}\n")
                    
                    if s % 100 == 0:
                        print(f"Step {s:04d} | Topo-Core Volume: {cores} | KE: {ke:.2f} | {status}")
                        
        print("Simulation complete.")

if __name__ == "__main__":
    sim = RationalMapSim()
    sim.initialize_rational_map()
    sim.run()
