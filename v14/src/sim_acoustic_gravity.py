import numpy as np
import os
from scipy.stats import linregress

print("HFKT v14: Avenue 3 - Acoustic Gravity Simulator")
print("===============================================")
print("Testing if continuous high-frequency resonances generate 1/r^2 macro-forces.")

os.makedirs("v14/data", exist_ok=True)
out_tsv = "v14/data/sim_acoustic_gravity_raw.tsv"

class AcousticGravitySim:
    def __init__(self, N=60, L=40.0, f_res=10.0):
        self.N = N
        self.L = L
        
        x = np.linspace(-L/2, L/2, N)
        y = np.linspace(-L/2, L/2, N)
        z = np.linspace(-L/2, L/2, N)
        self.X, self.Y, self.Z = np.meshgrid(x, y, z, indexing='ij')
        
        self.dx = L / (N - 1)
        self.omega = 2 * np.pi * f_res
        self.k = self.omega / 1.0 # arbitrary c=1
        
        self.rho = 1.0 # Medium density

    def compute_force_at_distance(self, D):
        """ Computes the hydrodynamic/acoustic average force on Sink 2 due to vibrating Source 1. """
        
        # Position 1 at origin
        pos1 = np.array([0, 0, 0])
        # Position 2 at D along x-axis
        pos2 = np.array([D, 0, 0])
        
        r1 = np.sqrt((self.X - pos1[0])**2 + (self.Y - pos1[1])**2 + (self.Z - pos1[2])**2) + 1e-5
        r2 = np.sqrt((self.X - pos2[0])**2 + (self.Y - pos2[1])**2 + (self.Z - pos2[2])**2) + 1e-5
        
        # We time-average over 1 full period.
        steps = 20
        dt = (2 * np.pi / self.omega) / steps
        
        F_net_x = 0.0
        
        # Radii to define the "body" of particle 2 to integrate force over
        R_core = 1.5
        mask2 = r2 <= R_core
        
        # We need the gradient of the time-averaged pressure
        # P_static = - 0.5 * rho * < |v1 + v2|^2 >
        
        v2_sq_sum = np.zeros((self.N, self.N, self.N))
        
        for i in range(steps):
            t = i * dt
            
            # v = grad( A/r * cos(kr - wt) ) 
            # dv/dr = -A/r^2 cos() - A k/r sin()
            
            # Source 1
            A1 = 5.0
            phase1 = self.k * r1 - self.omega * t
            dv1_dr = -A1/r1**2 * np.cos(phase1) - A1*self.k/r1 * np.sin(phase1)
            v1_x = dv1_dr * (self.X - pos1[0]) / r1
            v1_y = dv1_dr * (self.Y - pos1[1]) / r1
            v1_z = dv1_dr * (self.Z - pos1[2]) / r1
            
            # Source 2 (vibrating in phase, resonance locking)
            A2 = 5.0
            phase2 = self.k * r2 - self.omega * t
            dv2_dr = -A2/r2**2 * np.cos(phase2) - A2*self.k/r2 * np.sin(phase2)
            v2_x = dv2_dr * (self.X - pos2[0]) / r2
            v2_y = dv2_dr * (self.Y - pos2[1]) / r2
            v2_z = dv2_dr * (self.Z - pos2[2]) / r2
            
            v_tot_x = v1_x + v2_x
            v_tot_y = v1_y + v2_y
            v_tot_z = v1_z + v2_z
            
            v2_sq = v_tot_x**2 + v_tot_y**2 + v_tot_z**2
            v2_sq_sum += v2_sq
            
        v2_sq_avg = v2_sq_sum / steps
        P_static = -0.5 * self.rho * v2_sq_avg
        
        # Force is the integral of the pressure gradient over the volume of particle 2
        dP_dx = np.gradient(P_static, self.dx, axis=0)
        
        # Integral over the core
        F_net_x = np.sum(-dP_dx[mask2]) * (self.dx**3)
        
        return F_net_x

    def run_distance_sweep(self):
        print(f"Streaming Acoustic Gravity derivation to {out_tsv}...")
        
        distances = np.linspace(5.0, 15.0, 10)
        forces = []
        
        with open(out_tsv, "w") as f:
            f.write("Distance(D)\tAcoustic_Force(F)\n")
            
            for D in distances:
                F = self.compute_force_at_distance(D)
                f.write(f"{D:.4f}\t{F:.8e}\n")
                forces.append(F)
                print(f"D = {D:.2f} | Force = {F:.4e}")
                
        # Calculate force law exponent F ~ D^n
        log_D = np.log(distances)
        # Taking the absolute force as we just want the magnitude decay law
        # The force calculation here yields an attractive (negative) force if in phase
        log_F = np.log(np.abs(forces))
        
        slope, intercept, r_value, p_value, std_err = linregress(log_D, log_F)
        print("\nDistance Sweep Complete.")
        print(f"Extracted Force Law Exponent: n = {slope:.4f} (Ideal Gravity: -2.0000)")
        print(f"Acoustic Gravity F ~ 1/D^{abs(slope):.4f}")

if __name__ == "__main__":
    sim = AcousticGravitySim()
    sim.run_distance_sweep()
