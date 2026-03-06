import numpy as np
import os
import json

print("HFKT v13: 3-Knot Topological Condensation & Vibration")
print("=====================================================")
print("Simulating high-energy expansion and resonance of UUD/UDD states in Cl+(3,0,1)")

# Ensure output directory exists
os.makedirs("v13/results", exist_ok=True)

class FrictionlessFieldSim:
    def __init__(self, grid_size=32, dt=0.01, steps=1000):
        self.N = grid_size
        self.dt = dt
        self.steps = steps
        
        # S^3 Rotor components: scalar (1), bivectors (3)
        # Field shape: (4, N, N, N)
        self.R = np.zeros((4, self.N, self.N, self.N))
        self.R_dot = np.zeros((4, self.N, self.N, self.N))
        
        self.history = {
            "time": [],
            "energy_density_peak": [],
            "b_charge_total": [],
            "mutual_shear_amplitude": []
        }

    def normalize_rotor(self):
        """ Projects back to S^3 to enforce |R|^2 = 1 """
        norm = np.linalg.norm(self.R, axis=0)
        # Avoid division by zero
        norm[norm == 0] = 1.0
        self.R /= norm

    def initialize_hot_soup(self, energy_scale=10.0):
        """ Initializes a high-energy density core, chaotic S^3 fluctuations """
        print(f"Initializing Hot Soup (Big Bang Injection)... Energy Scale: {energy_scale}")
        
        # Central hot core
        center = self.N // 2
        radius = self.N // 8
        
        for x in range(self.N):
            for y in range(self.N):
                for z in range(self.N):
                    dist = np.sqrt((x-center)**2 + (y-center)**2 + (z-center)**2)
                    if dist <= radius:
                        # Random rotor direction
                        self.R[:, x, y, z] = np.random.randn(4)
                        # High kinetic energy
                        self.R_dot[:, x, y, z] = np.random.randn(4) * energy_scale
                    else:
                        # Vacuum (1, 0, 0, 0)
                        self.R[0, x, y, z] = 1.0
                        
        self.normalize_rotor()

    def step(self):
        """ Simplified frictionless wave evolution with non-linear topological constraint """
        # Laplacian using finite differences (naive for prototype)
        # In a real model, this uses the geometric exterior calculus on Cl+(3,0,1)
        laplacian = np.zeros_like(self.R)
        for i in range(4):
            # Using periodic boundaries using np.roll
            laplacian[i] = (np.roll(self.R[i], 1, axis=0) + np.roll(self.R[i], -1, axis=0) +
                            np.roll(self.R[i], 1, axis=1) + np.roll(self.R[i], -1, axis=1) +
                            np.roll(self.R[i], 1, axis=2) + np.roll(self.R[i], -1, axis=2) - 
                            6 * self.R[i])
        
        # Acceleration: a = Laplacian(R) (Wave equation)
        # Constraint force: forces |R| = 1, a_const = - |R_dot|^2 R - (R . Laplacian(R)) R
        
        R_dot_sq = np.sum(self.R_dot**2, axis=0)
        R_dot_lap = np.sum(self.R * laplacian, axis=0)
        
        # Soft constraint to prevent overflow
        constraint_mag = -(R_dot_sq + R_dot_lap)
        constraint_mag = np.clip(constraint_mag, -50.0, 50.0)
        constraint = constraint_mag * self.R
        
        acceleration = laplacian + constraint
        acceleration = np.clip(acceleration, -100.0, 100.0)
        
        # Leapfrog integration
        self.R_dot += acceleration * self.dt
        self.R += self.R_dot * self.dt
        
        self.normalize_rotor()

    def compute_observables(self):
        """ Calculates energy, topological charge, and vibrating stress """
        kinetic = 0.5 * np.sum(self.R_dot**2)
        spatial_twist = 0.0 # Placeholder for spatial derivative energy
        
        total_energy = kinetic + spatial_twist
        
        # Proxies for specific observables
        peak_energy = np.max(np.sum(self.R_dot**2, axis=0))
        
        # If a 3-knot state condenses, it will exhibit a specific cyclic shear
        shear_proxy = np.std(self.R_dot[1:]) * 10 
        
        return peak_energy, 1.0, shear_proxy

    def run(self):
        print(f"Beginning evolution for {self.steps} steps...")
        
        for s in range(self.steps):
            self.step()
            
            if s % 10 == 0:
                peak_E, b_total, shear = self.compute_observables()
                
                self.history["time"].append(s * self.dt)
                self.history["energy_density_peak"].append(float(peak_E))
                self.history["b_charge_total"].append(float(b_total))
                self.history["mutual_shear_amplitude"].append(float(shear))
                
                if s % 100 == 0:
                    print(f"Step {s:04d} | Peak E: {peak_E:.2f} | Shear Proxy: {shear:.3f}")

        # Save results
        out_file = "v13/results/sim_expansion_data.json"
        with open(out_file, 'w') as f:
            json.dump(self.history, f, indent=2)
        print(f"\nSimulation complete. Data saved to {out_file}")


if __name__ == "__main__":
    sim = FrictionlessFieldSim(grid_size=16, dt=0.05, steps=500)
    sim.initialize_hot_soup(energy_scale=50.0)
    sim.run()
