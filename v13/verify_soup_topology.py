import numpy as np
import os
import json

print("HFKT v13: Avenue 4 - Hot Soup Topology Diagnostic")
print("==================================================")
print("Re-running the Hot Soup expansion to explicitly count distinct cores in the steady state.")

out_tsv = "v13/results/soup_topology_diagnostic.tsv"

class SoupDiagnosticSim:
    def __init__(self, grid_size=20, dt=0.02, steps=600):
        self.N = grid_size
        self.dt = dt
        self.steps = steps
        
        self.R = np.zeros((4, self.N, self.N, self.N))
        self.R_dot = np.zeros((4, self.N, self.N, self.N))
        
    def normalize(self):
        norm = np.linalg.norm(self.R, axis=0)
        norm[norm == 0] = 1.0
        self.R /= norm

    def initialize_hot_soup(self, energy_scale=50.0):
        print(f"Injecting Hot Soup (Energy Scale: {energy_scale})...")
        center = self.N // 2
        radius = self.N // 6
        
        for x in range(self.N):
            for y in range(self.N):
                for z in range(self.N):
                    dist = np.sqrt((x-center)**2 + (y-center)**2 + (z-center)**2)
                    if dist <= radius:
                        self.R[:, x, y, z] = np.random.randn(4)
                        self.R_dot[:, x, y, z] = np.random.randn(4) * energy_scale
                    else:
                        self.R[0, x, y, z] = 1.0
                        
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

    def find_distinct_cores(self):
        """ Evaluates the topological profile by looking for distinct spatial minima of R_w (the scalar part).
            In a true B=1 knot, R_w approaches -1 at the core and 1 at vacuum.
            If there are multiple knots, there should be multiple distinct R_w minima. """
            
        Rw = self.R[0]
        
        # We classify a 'core' as any grid point where Rw < -0.5
        # and we group adjacent core points to count distinct entities
        core_mask = Rw < -0.5
        
        from scipy.ndimage import label
        labeled_array, num_features = label(core_mask)
        
        # Also let's check Peak Energy just to track the generic lump volume
        volume_excited = np.sum(np.sum(self.R_dot**2, axis=0) > 0.05)
        
        return num_features, int(volume_excited)

    def run(self):
        print(f"Streaming topological counts to {out_tsv}...")
        
        with open(out_tsv, "w") as f:
            f.write("Time\tDistinct_Topological_Cores\tExcited_Volume\n")
            
            for s in range(self.steps):
                self.step()
                
                t = s * self.dt
                if s % 10 == 0:
                    cores, vol = self.find_distinct_cores()
                    f.write(f"{t:.3f}\t{cores}\t{vol}\n")
                    
                    if s % 100 == 0:
                        print(f"Step {s:04d} | Cores: {cores}, Excited Grid Vol: {vol}")
                        
        print("Diagnostic complete.")

if __name__ == "__main__":
    sim = SoupDiagnosticSim()
    sim.initialize_hot_soup()
    sim.run()
