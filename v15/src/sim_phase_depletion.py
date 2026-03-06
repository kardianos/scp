import numpy as np
import os
from scipy.stats import linregress

print("HFKT v15: Avenue 2 - Volumetric Phase Depletion")
print("===============================================")
print("Testing if a static 1/r scalar deficit induces 1/r^2 kinematic attraction.")

os.makedirs("v15/data", exist_ok=True)
out_tsv = "v15/data/sim_phase_depletion_raw.tsv"

class PhaseDepletionSim:
    def __init__(self, N=60, L=40.0):
        self.N = N
        self.L = L
        
        x = np.linspace(-L/2, L/2, N)
        y = np.linspace(-L/2, L/2, N)
        z = np.linspace(-L/2, L/2, N)
        self.X, self.Y, self.Z = np.meshgrid(x, y, z, indexing='ij')
        
        self.dx = L / (N - 1)
        self.dt = 0.05
        self.c = 1.0

    def compute_force_at_distance(self, D):
        """ Evaluates the acceleration of a resonant topology moving within a 1/r phase potential. """
        
        # Source of the Gravity (Mass M at origin) creates a 1/r scalar potential deficit
        # Phi(r) = - G*M / r
        GM = 5.0
        
        r1 = np.sqrt(self.X**2 + self.Y**2 + self.Z**2) + 1e-5
        Phi = -GM / r1
        
        # We model Particle 2 (the test mass) as a localized resonant phase echo at position D
        pos2 = np.array([D, 0, 0])
        r2 = np.sqrt((self.X - pos2[0])**2 + (self.Y - pos2[1])**2 + (self.Z - pos2[2])**2) + 1e-5
        
        # The equation of motion for a wave-envelope in a scalar potential implies
        # the potential acts as an effective index of refraction OR a direct pressure gradient.
        # F_net = - integral( rho_particle * grad(Phi) ) dV
        
        # The topological "density" or mass footprint of the resonant trap
        # Let's say it's an exponentially localized energy core
        rho_particle = np.exp(-r2 / 1.5)
        
        # Gradient of the Potential (Gravity field g = -grad(Phi))
        # grad(-GM/r) = GM/r^2 * (r_vec/r) pointing inward
        g_x = -GM / r1**2 * (self.X / r1)
        
        # Integrate the interaction over the volume of the test particle
        F_net_x = np.sum(rho_particle * g_x) * (self.dx**3)
        
        return F_net_x

    def run_distance_sweep(self):
        print(f"Streaming Phase Depletion derivations to {out_tsv}...")
        
        distances = np.linspace(5.0, 15.0, 10)
        forces = []
        
        with open(out_tsv, "w") as f:
            f.write("Distance(D)\tDepletion_Force(F)\n")
            
            for D in distances:
                F = self.compute_force_at_distance(D)
                f.write(f"{D:.4f}\t{F:.8e}\n")
                forces.append(np.abs(F))
                print(f"D = {D:.2f} | Force = {F:.4e}")
                
        # Calculate force law exponent F ~ D^n
        log_D = np.log(distances)
        log_F = np.log(forces)
        
        slope, intercept, r_value, p_value, std_err = linregress(log_D, log_F)
        print("\nDistance Sweep Complete.")
        print(f"Extracted Force Law Exponent: n = {slope:.4f} (Ideal Gravity: -2.0000)")
        print(f"Phase Depletion F ~ 1/D^{abs(slope):.4f}")

if __name__ == "__main__":
    sim = PhaseDepletionSim()
    sim.run_distance_sweep()
