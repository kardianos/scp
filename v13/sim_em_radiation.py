import numpy as np
import os

print("HFKT v13: Avenue 1 - Topological EM Radiation Simulator")
print("=====================================================")
print("Evolving a spinning B=1 knot and extracting emergent Far-Field E & B vectors.")

os.makedirs("v13/results", exist_ok=True)
out_tsv = "v13/results/sim_em_radiation_raw.tsv"

class EMRadiationSim:
    def __init__(self, grid_size=40, dt=0.02, steps=600):
        self.N = grid_size
        self.dt = dt
        self.steps = steps
        self.center = self.N // 2
        
        self.R = np.zeros((4, self.N, self.N, self.N))
        self.R_dot = np.zeros((4, self.N, self.N, self.N))
        
    def normalize(self):
        norm = np.linalg.norm(self.R, axis=0)
        norm[norm == 0] = 1.0
        self.R /= norm

    def initialize_spinning_knot(self):
        print("Initializing a localized B=1 topological spinning knot...")
        for x in range(self.N):
            for y in range(self.N):
                for z in range(self.N):
                    # Distance from center
                    dx, dy, dz = x - self.center, y - self.center, z - self.center
                    r = np.sqrt(dx**2 + dy**2 + dz**2)
                    
                    if r == 0:
                        self.R[:, x, y, z] = [0, 1, 0, 0] # Core flipped
                        continue
                        
                    # Profile function (Hopf-like fade to vacuum 1,0,0,0)
                    theta = np.pi * np.exp(-r / 3.0)
                    
                    # Spin axis along Z, creating an equatorial twist
                    n_x = -dy/r
                    n_y = dx/r
                    n_z = dz/r
                    
                    # Normalizing director
                    n_norm = np.sqrt(n_x**2 + n_y**2 + n_z**2)
                    if n_norm > 0:
                        n_x, n_y, n_z = n_x/n_norm, n_y/n_norm, n_z/n_norm
                        
                    self.R[0, x, y, z] = np.cos(theta)
                    self.R[1, x, y, z] = n_x * np.sin(theta)
                    self.R[2, x, y, z] = n_y * np.sin(theta)
                    self.R[3, x, y, z] = n_z * np.sin(theta)
                    
                    # Give it an initial spin (R_dot) to force oscillation/radiation
                    # Spinning in the e1-e2 plane (phase rotation)
                    omega = 0.5 * np.exp(-r/2.0)
                    self.R_dot[1, x, y, z] += -omega * self.R[2, x, y, z]
                    self.R_dot[2, x, y, z] += omega * self.R[1, x, y, z]
                    
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

    def extract_EM_fields(self, obs_r):
        """ Calculates emergent E and B magnitudes at a specific radius observation shell """
        # Let's observe at (center+obs_r, center, center) on the x-axis
        cx = self.center + obs_r
        cy, cz = self.center, self.center
        
        # Ensure within grid
        if cx >= self.N - 1: return 0, 0
        
        # Get local R and R_dot
        R_local = self.R[:, cx, cy, cz]
        R_dot_local = self.R_dot[:, cx, cy, cz]
        R_tilde = np.array([R_local[0], -R_local[1], -R_local[2], -R_local[3]])
        
        # E = (R_dot) * R~
        # E magnitude is just the norm of the bivector parts
        E_mag = np.linalg.norm([
            R_dot_local[0]*R_tilde[1] + R_dot_local[1]*R_tilde[0] + R_dot_local[2]*R_tilde[3] - R_dot_local[3]*R_tilde[2],
            R_dot_local[0]*R_tilde[2] - R_dot_local[1]*R_tilde[3] + R_dot_local[2]*R_tilde[0] + R_dot_local[3]*R_tilde[1],
            R_dot_local[0]*R_tilde[3] + R_dot_local[1]*R_tilde[2] - R_dot_local[2]*R_tilde[1] + R_dot_local[3]*R_tilde[0]
        ])
        
        # B = spatial twist. We approximate magnitude via difference to neighbor
        R_dx = (self.R[:, cx+1, cy, cz] - self.R[:, cx-1, cy, cz]) / 2.0
        B_mag = np.linalg.norm([
            R_dx[0]*R_tilde[1] + R_dx[1]*R_tilde[0], 
            R_dx[0]*R_tilde[2] + R_dx[2]*R_tilde[0], 
            R_dx[0]*R_tilde[3] + R_dx[3]*R_tilde[0]
        ])
        
        return float(E_mag), float(B_mag)

    def run(self):
        print(f"Streaming time-series EM radiation data to {out_tsv}...")
        
        # We will observe the wave passing 3 different far-field radii
        r1, r2, r3 = 5, 10, 15
        
        with open(out_tsv, "w") as f:
            f.write("Time\tE_r5\tB_r5\tE_r10\tB_r10\tE_r15\tB_r15\n")
            
            for s in range(self.steps):
                self.step()
                
                t = s * self.dt
                e1, b1 = self.extract_EM_fields(r1)
                e2, b2 = self.extract_EM_fields(r2)
                e3, b3 = self.extract_EM_fields(r3)
                
                f.write(f"{t:.3f}\t{e1:.6f}\t{b1:.6f}\t{e2:.6f}\t{b2:.6f}\t{e3:.6f}\t{b3:.6f}\n")
                
                if s % 100 == 0:
                    print(f"Step {s:04d} | Radiating (r={r2}): E={e2:.4f}, B={b2:.4f}")
        
        print("Simulation complete! Extracted field magnitudes at fixed radius shells.")

if __name__ == "__main__":
    sim = EMRadiationSim()
    sim.initialize_spinning_knot()
    sim.run()
