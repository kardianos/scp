import numpy as np
import os
import scipy.signal

print("HFKT v13: Avenue 4 - 3-Knot UUD Mutual Vibration Simulator")
print("==========================================================")
print("Evolving an entangled UUD state and explicitly tracking mutual knot distance.")

os.makedirs("v13/results", exist_ok=True)
out_tsv = "v13/results/sim_3knot_vibration_raw.tsv"

class CompositeKnotSim:
    def __init__(self, grid_size=40, dt=0.02, steps=1000):
        self.N = grid_size
        self.dt = dt
        self.steps = steps
        
        self.R = np.zeros((4, self.N, self.N, self.N))
        self.R_dot = np.zeros((4, self.N, self.N, self.N))
        self.centers = []
        
    def normalize(self):
        norm = np.linalg.norm(self.R, axis=0)
        norm[norm == 0] = 1.0
        self.R /= norm

    def quat_mult_array(self, R1, R2):
        """ Multiplies two quaternion arrays point-wise """
        w1, x1, y1, z1 = R1[0], R1[1], R1[2], R1[3]
        w2, x2, y2, z2 = R2[0], R2[1], R2[2], R2[3]
        return np.array([
            w1*w2 - x1*x2 - y1*y2 - z1*z2,
            w1*x2 + x1*w2 + y1*z2 - z1*y2,
            w1*y2 - x1*z2 + y1*w2 + z1*x2,
            w1*z2 + x1*y2 - y1*x2 + z1*w2
        ])

    def inject_knot(self, cx, cy, cz, phase_shift, isospin_sign):
        """ Injects a fractional topological twist at a given center.
            Properly composes knots via quaternion multiplication to preserve S^3 topology """
        
        K_new = np.zeros_like(self.R)
        K_new[0] = 1.0 # Default vacuum
        
        for x in range(self.N):
            for y in range(self.N):
                for z in range(self.N):
                    dx, dy, dz = x - cx, y - cy, z - cz
                    r = np.sqrt(dx**2 + dy**2 + dz**2)
                    
                    if r == 0: 
                        K_new[:, x, y, z] = [0, 1, 0, 0]
                        continue
                        
                    theta = np.pi * np.exp(-r / 3.0)
                    
                    n_x = -dy/r * isospin_sign
                    n_y = dx/r * isospin_sign
                    n_z = dz/r
                    
                    phase_angle = phase_shift * (2 * np.pi / 3.0)
                    
                    K_new[0, x, y, z] = np.cos(theta)
                    K_new[1, x, y, z] = n_x * np.sin(theta) * np.cos(phase_angle)
                    K_new[2, x, y, z] = n_y * np.sin(theta) * np.sin(phase_angle)
                    K_new[3, x, y, z] = n_z * np.sin(theta)
                    
                    self.R_dot[1, x, y, z] += 0.1 * np.exp(-r/2.0)
                    
        # Compose via exact rotation to prevent ruining the S^3 geometry
        self.R = self.quat_mult_array(self.R, K_new)

    def initialize_uud_state(self):
        print("Initializing composite UUD state (Proton Analogue)...")
        self.R[0] = 1.0 # Base vacuum
        c = self.N // 2
        
        # Place 3 knots in a triangle
        self.inject_knot(c, c+4, c, phase_shift=0, isospin_sign=1)   # U1
        self.inject_knot(c-3, c-2, c, phase_shift=1, isospin_sign=1) # U2
        self.inject_knot(c+3, c-2, c, phase_shift=2, isospin_sign=-1)# D1
        
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
        constraint_mag = np.clip(constraint_mag, -20.0, 20.0)
        
        acceleration = laplacian + constraint_mag * self.R
        acceleration = np.clip(acceleration, -20.0, 20.0)
        
        self.R_dot += acceleration * self.dt
        self.R += self.R_dot * self.dt
        self.normalize()

    def find_knot_centers(self):
        """ Finds the 3 local maxima of kinetic energy to proxy knot centers """
        energy = np.sum(self.R_dot**2, axis=0)
        
        # A simple non-maximum suppression to find top 3 peaks
        peaks = []
        # Flatten and sort indices
        flat_indices = np.argsort(energy.flatten())[::-1]
        
        for idx in flat_indices:
            pt = np.unravel_index(idx, energy.shape)
            
            # Ensure it is far enough from existing peaks to be a distinct knot
            is_valid = True
            for p in peaks:
                dist = np.sqrt(sum((pt[i] - p[i])**2 for i in range(3)))
                if dist < 3.0: # Minimum exclusion radius between cores
                    is_valid = False
                    break
                    
            if is_valid:
                peaks.append(pt)
            if len(peaks) == 3:
                break
                
        # If simulation washes out and we can't find 3, return defaults
        while len(peaks) < 3:
            peaks.append((self.N//2, self.N//2, self.N//2))
            
        return peaks

    def extract_mutual_distances(self):
        peaks = self.find_knot_centers()
        # Calculate pair-wise distances (Isosceles geometry implies 1 distinct, 2 similar)
        d12 = np.sqrt(sum((peaks[0][i] - peaks[1][i])**2 for i in range(3)))
        d23 = np.sqrt(sum((peaks[1][i] - peaks[2][i])**2 for i in range(3)))
        d13 = np.sqrt(sum((peaks[0][i] - peaks[2][i])**2 for i in range(3)))
        
        # Sort distances so we can consistently track the 'long' vs 'short' sides
        distances = sorted([d12, d23, d13])
        return distances

    def run(self):
        print(f"Streaming time-series UUD vibration data to {out_tsv}...")
        
        with open(out_tsv, "w") as f:
            f.write("Time\tD_short1\tD_short2\tD_long\n")
            
            for s in range(self.steps):
                self.step()
                
                t = s * self.dt
                dist = self.extract_mutual_distances()
                
                # dist[0], dist[1] should be similar (Isosceles legs), dist[2] is the base
                f.write(f"{t:.3f}\t{dist[0]:.4f}\t{dist[1]:.4f}\t{dist[2]:.4f}\n")
                
                if s % 100 == 0:
                    print(f"Step {s:04d} | UUD Distances: {dist[0]:.2f}, {dist[1]:.2f}, {dist[2]:.2f}")
        
        print(f"Simulation complete! Vibration tracked and explicitly logged to disk.")

if __name__ == "__main__":
    sim = CompositeKnotSim()
    sim.initialize_uud_state()
    sim.run()
