import numpy as np
import os

print("HFKT v13: Avenue 4 - Causal Retarded-Time Phase Echo Simulator")
print("================================================================")
print("Testing if c*dt delayed wave interference natively prevents topological decay.")

os.makedirs("v13/results", exist_ok=True)
out_tsv = "v13/results/sim_retarded_time_resonance_raw.tsv"

class CausalEchoSim:
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

    def initialize_phase_echo(self):
        """ Initializes a field using the interference pattern of 3 spherical emitters """
        print("Constructing Retarded-Time Interference pattern...")
        
        k = 2.0
        
        # 3 distinct emission centers (Isosceles)
        d1 = np.sqrt((self.X-3)**2 + self.Y**2 + self.Z**2)
        d2 = np.sqrt((self.X+1.5)**2 + (self.Y-2.6)**2 + self.Z**2)
        d3 = np.sqrt((self.X+1.5)**2 + (self.Y+2.6)**2 + self.Z**2)
        
        # Simulate c*dt delay visually matching the phi = kr - wt wave constraint 
        interference = np.cos(k*d1) + np.cos(k*d2 + 2*np.pi/3) + np.cos(k*d3 + 4*np.pi/3)
        
        # We also need to give it an explicit internal R_dot (oscillation) 
        # that drives the resonance forward in time
        w = 3.0
        d_interf_dt = w*np.sin(k*d1) + w*np.sin(k*d2 + 2*np.pi/3) + w*np.sin(k*d3 + 4*np.pi/3)
        
        # Exonential localized boundary
        theta = np.pi * np.exp(-np.sqrt(self.X**2 + self.Y**2 + self.Z**2) / 4.0)
        
        # Initial Spatial gradients mapped to the vector components
        nx = np.gradient(interference, axis=0)
        ny = np.gradient(interference, axis=1)
        nz = np.gradient(interference, axis=2)
        
        norm = np.sqrt(nx**2 + ny**2 + nz**2) + 1e-10
        self.R[0] = np.cos(theta)
        self.R[1] = (nx/norm) * np.sin(theta)
        self.R[2] = (ny/norm) * np.sin(theta)
        self.R[3] = (nz/norm) * np.sin(theta)
        
        # The dynamic velocity is driven by the time-derivative of the echo
        self.R_dot[1] = d_interf_dt * 0.1 * np.exp(-np.sqrt(self.X**2 + self.Y**2)/3.0)
        
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
        core_mass = np.sum(self.R[0] < -0.5)
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
                    status = "Resonating" if cores > 10 else ("Dissolving" if cores > 0 else "Annihilated")
                    f.write(f"{t:.3f}\t{cores}\t{ke:.4f}\t{status}\n")
                    
                    if s % 100 == 0:
                        print(f"Step {s:04d} | Topo-Core Volume: {cores} | KE: {ke:.2f} | {status}")
                        
        print("Simulation complete.")

if __name__ == "__main__":
    sim = CausalEchoSim()
    sim.initialize_phase_echo()
    sim.run()
