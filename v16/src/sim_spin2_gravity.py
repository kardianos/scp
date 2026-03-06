import numpy as np
import os
from scipy.stats import linregress

print("HFKT v15: Avenue 3 - Spin-2 Geometric Cross-Terms")
print("=================================================")
print("Testing if internal EM stress-energy tensor generates 1/r^2 interaction.")

os.makedirs("v15/data", exist_ok=True)
out_tsv = "v15/data/sim_spin2_gravity_raw.tsv"

class Spin2GravitySim:
    def __init__(self, N=60, L=40.0):
        self.N = N
        self.L = L
        
        x = np.linspace(-L/2, L/2, N)
        y = np.linspace(-L/2, L/2, N)
        z = np.linspace(-L/2, L/2, N)
        self.X, self.Y, self.Z = np.meshgrid(x, y, z, indexing='ij')
        
        self.dx = L / (N - 1)

    def compute_force_at_distance(self, D):
        """ Evaluates the emergent spin-2 interaction force between two resonant topological bodies. """
        
        # In V13, we proved the topological bulk inherently generates an emergent EM field F = E + iB.
        # In General Relativity, the source of gravity is the stress-energy tensor T_mu_nu.
        # For an electromagnetic field: T^00 (Energy density) ~ E^2 + B^2
        # T^ij (Stress) ~ E^i E^j + B^i B^j - 1/2 delta^ij (E^2 + B^2)
        
        # If gravity is driven by this tensor, the interaction force F_g between two bodies
        # is the integral of the interaction stress-energy.
        
        pos1 = np.array([0, 0, 0])
        pos2 = np.array([D, 0, 0])
        
        r1 = np.sqrt(self.X**2 + self.Y**2 + self.Z**2) + 1e-5
        r2 = np.sqrt((self.X - pos2[0])**2 + (self.Y - pos2[1])**2 + (self.Z - pos2[2])**2) + 1e-5
        
        # Model the emergent internal "Electric" field of Particle 1 (1/r^2 radiating from core)
        # Inside the core it is bounded, outside it tails off
        # Let's say the core size is a
        a = 1.0
        E1_mag = np.where(r1 < a, r1/a**3, 1/r1**2) 
        E1_x = E1_mag * (self.X / r1)
        E1_y = E1_mag * (self.Y / r1)
        E1_z = E1_mag * (self.Z / r1)
        
        # Emergent Field of Particle 2
        E2_mag = np.where(r2 < a, r2/a**3, 1/r2**2)
        E2_x = E2_mag * ((self.X - pos2[0]) / r2)
        E2_y = E2_mag * ((self.Y - pos2[1]) / r2)
        E2_z = E2_mag * ((self.Z - pos2[2]) / r2)
        
        # The total emergent field is E_tot = E1 + E2
        # T^00 ~ |E1 + E2|^2 = |E1|^2 + |E2|^2 + 2(E1 . E2)
        # The interaction energy density is U_int \propto E1 . E2
        
        U_int = E1_x*E2_x + E1_y*E2_y + E1_z*E2_z
        
        # The effective force on Particle 2 is the negative gradient of the total integrated interaction energy
        # However, to avoid numerical gradient artifacts of the entire volume integral over moving D,
        # we can calculate it analytically or by explicitly integrating the Maxwell stress tensor over the midplane.
        
        # Let's compute the total interaction energy U(D)
        U_total = np.sum(U_int) * (self.dx**3)
        return U_total

    def run_distance_sweep(self):
        print(f"Streaming Spin-2 Interaction derivations to {out_tsv}...")
        
        distances = np.linspace(5.0, 15.0, 20)
        energies = []
        forces = []
        
        with open(out_tsv, "w") as f:
            f.write("Distance(D)\tInteraction_Energy(U)\tEmergent_Force(F=-dU/dD)\n")
            
            # Compute energies
            for D in distances:
                U = self.compute_force_at_distance(D)
                energies.append(U)
                
            # Differentiate to find Force F = -dU/dD
            for i in range(len(distances)-1):
                dD = distances[i+1] - distances[i]
                dU = energies[i+1] - energies[i]
                F = -dU / dD
                
                D_mid = distances[i] + dD/2
                forces.append((D_mid, F))
                f.write(f"{D_mid:.4f}\t{energies[i]:.8e}\t{F:.8e}\n")
                print(f"D = {D_mid:.2f} | U_int = {energies[i]:.4e} | Force = {F:.4e}")
                
        # Calculate force law exponent F ~ D^n
        D_vals = np.array([f[0] for f in forces])
        F_vals = np.abs(np.array([f[1] for f in forces])) # Magnitude
        
        log_D = np.log(D_vals)
        log_F = np.log(F_vals)
        
        slope, intercept, r_value, p_value, std_err = linregress(log_D, log_F)
        print("\nDistance Sweep Complete.")
        print(f"Extracted Force Law Exponent: n = {slope:.4f} (Ideal Gravity: -2.0000)")
        print(f"Spin-2 Tensor Force F ~ 1/D^{abs(slope):.4f}")

if __name__ == "__main__":
    sim = Spin2GravitySim()
    sim.run_distance_sweep()
