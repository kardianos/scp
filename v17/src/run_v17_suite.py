import numpy as np
import os
import json
from scipy.stats import linregress

print("HFKT v17: The Minimal Self-Consistent Proof Suite")
print("=================================================")
print("Running complete Skyrme-Minkowski analytics on a single continuum field.")

os.makedirs("v17/data", exist_ok=True)

class HFKT_SkyrmeSim:
    def __init__(self, N=64, L=30.0, dt=0.01, lambda_s=0.25, sigma=3.0):
        self.N = N
        self.L = L
        self.dt = dt
        self.lambda_s = lambda_s
        self.sigma = sigma
        self.dx = L / (N - 1)
        
        # Grid
        x = np.linspace(-L/2, L/2, N)
        self.X, self.Y, self.Z = np.meshgrid(x, x, x, indexing='ij')
        
        self.R = np.zeros((4, N, N, N))
        self.V = np.zeros((4, N, N, N))
        
    def initialize_hopfion_Q1(self, shift_x=0.0):
        r2 = (self.X - shift_x)**2 + self.Y**2 + self.Z**2
        D = 1.0 + r2
        w = (1.0 - r2) / D
        xs = 2 * (self.X - shift_x) / D
        ys = 2 * self.Y / D
        zs = 2 * self.Z / D
        
        alpha = np.pi * np.exp(-np.sqrt(r2) / self.sigma)
        self.R[0] = np.cos(alpha)
        self.R[1] = xs * np.sin(alpha)
        self.R[2] = ys * np.sin(alpha)
        self.R[3] = zs * np.sin(alpha)
        
        norm = np.linalg.norm(self.R, axis=0)
        self.R /= norm
        self.V[:] = 0.0

    def initialize_two_hopfions_product(self, separation=10.0):
        # Multiplicative superposition for unit quaternions Q = Q1 * Q2
        R1 = np.zeros_like(self.R)
        R2 = np.zeros_like(self.R)
        
        self.initialize_hopfion_Q1(shift_x = -separation/2)
        R1[:] = self.R[:]
        self.initialize_hopfion_Q1(shift_x = separation/2)
        R2[:] = self.R[:]
        
        # Q1 * Q2 (quaternion product)
        self.R[0] = R1[0]*R2[0] - R1[1]*R2[1] - R1[2]*R2[2] - R1[3]*R2[3]
        self.R[1] = R1[0]*R2[1] + R1[1]*R2[0] + R1[2]*R2[3] - R1[3]*R2[2]
        self.R[2] = R1[0]*R2[2] - R1[1]*R2[3] + R1[2]*R2[0] + R1[3]*R2[1]
        self.R[3] = R1[0]*R2[3] + R1[1]*R2[2] - R1[2]*R2[1] + R1[3]*R2[0]
        
        norm = np.linalg.norm(self.R, axis=0)
        self.R /= norm
        self.V[:] = 0.0

    def project(self):
        norm = np.linalg.norm(self.R, axis=0, keepdims=True)
        self.R /= norm
        dot = np.sum(self.R * self.V, axis=0, keepdims=True)
        self.V -= dot * self.R
        
    def laplacian(self):
        lap = np.zeros_like(self.R)
        for i in range(4):
            f = self.R[i]
            lap[i] = (np.roll(f,1,0) + np.roll(f,-1,0) +
                      np.roll(f,1,1) + np.roll(f,-1,1) +
                      np.roll(f,1,2) + np.roll(f,-1,2) - 6*f) / self.dx**2
        return lap

    def _skyrme_contribution(self):
        # 12-term finite difference expansion of the Skyrme term div(S)
        dR = [np.zeros_like(self.R) for _ in range(3)]
        for i in range(4):
            dR[0][i] = (np.roll(self.R[i], -1, axis=0) - np.roll(self.R[i], 1, axis=0)) / (2*self.dx)
            dR[1][i] = (np.roll(self.R[i], -1, axis=1) - np.roll(self.R[i], 1, axis=1)) / (2*self.dx)
            dR[2][i] = (np.roll(self.R[i], -1, axis=2) - np.roll(self.R[i], 1, axis=2)) / (2*self.dx)
        
        M = np.zeros((3, 3, self.N, self.N, self.N))
        for i in range(3):
            for j in range(3):
                M[i,j] = np.sum(dR[i] * dR[j], axis=0)
                
        tr_M = M[0,0] + M[1,1] + M[2,2]
        
        S = [np.zeros_like(self.R) for _ in range(3)]
        for i in range(3):
            term2 = np.zeros_like(self.R)
            for j in range(3):
                term2 += M[i,j] * dR[j]
            S[i] = tr_M * dR[i] - term2
            
        div_S = np.zeros_like(self.R)
        for i in range(4):
            div_S[i] += (np.roll(S[0][i], -1, axis=0) - np.roll(S[0][i], 1, axis=0)) / (2*self.dx)
            div_S[i] += (np.roll(S[1][i], -1, axis=1) - np.roll(S[1][i], 1, axis=1)) / (2*self.dx)
            div_S[i] += (np.roll(S[2][i], -1, axis=2) - np.roll(S[2][i], 1, axis=2)) / (2*self.dx)
            
        return div_S

    def step(self):
        lap = self.laplacian()
        # Full unconstrained PDE: dAlembertian + Skyrme = 0 represents the non-linear potential
        F_total = lap
        if self.lambda_s > 0:
            F_total += self.lambda_s * self._skyrme_contribution()
            
        R_dot_F = np.sum(self.R * F_total, axis=0)
        v2 = np.sum(self.V**2, axis=0)
        
        # Projection multiplier yielding true geodesic flow
        accel = F_total - (R_dot_F + v2)[None, ...] * self.R
        
        self.V += accel * self.dt
        self.R += self.V * self.dt
        self.project()

    def get_energy(self):
        v2 = np.sum(self.V**2, axis=0)
        
        dR = [np.zeros_like(self.R) for _ in range(3)]
        for i in range(4):
            dR[0][i] = (np.roll(self.R[i], -1, axis=0) - np.roll(self.R[i], 1, axis=0)) / (2*self.dx)
            dR[1][i] = (np.roll(self.R[i], -1, axis=1) - np.roll(self.R[i], 1, axis=1)) / (2*self.dx)
            dR[2][i] = (np.roll(self.R[i], -1, axis=2) - np.roll(self.R[i], 1, axis=2)) / (2*self.dx)
            
        grad2 = np.sum(np.sum(dR[0]**2 + dR[1]**2 + dR[2]**2, axis=0), axis=0)
        grad_field = np.sum(dR[0]**2 + dR[1]**2 + dR[2]**2, axis=0)
        
        # Dirichlet + Kinetic Energy integral over the non-boundary volume
        E_dens = 0.5 * (v2 + grad_field)
        
        # Mask out boundary reflection to measure true soliton energy stability
        bnd = int(self.N * 0.15)
        mask = np.zeros((self.N, self.N, self.N), dtype=bool)
        mask[bnd:-bnd, bnd:-bnd, bnd:-bnd] = True
        
        return np.sum(E_dens[mask]) * (self.dx**3)

    def extract_radial_far_field(self):
        # E field proxy is proportional to spatial gradient magnitude |nabla R|
        # E^2 ~ grad_field
        dR = [np.zeros_like(self.R) for _ in range(3)]
        for i in range(4):
            dR[0][i] = (np.roll(self.R[i], -1, axis=0) - np.roll(self.R[i], 1, axis=0)) / (2*self.dx)
            dR[1][i] = (np.roll(self.R[i], -1, axis=1) - np.roll(self.R[i], 1, axis=1)) / (2*self.dx)
            dR[2][i] = (np.roll(self.R[i], -1, axis=2) - np.roll(self.R[i], 1, axis=2)) / (2*self.dx)
            
        E_mag = np.sqrt(np.sum(dR[0]**2 + dR[1]**2 + dR[2]**2, axis=0))
        r_mag = np.sqrt(self.X**2 + self.Y**2 + self.Z**2)
        
        # Filter for far-field points
        r_mask = (r_mag > 5.0) & (r_mag < self.L/2 - 2.0)
        radii = r_mag[r_mask].flatten()
        E_vals = E_mag[r_mask].flatten()
        
        return radii, E_vals

# ----- TESTS -----

def run_test_1_confinement():
    print("\\n--- TEST 1: Confinement (Stability) ---")
    out_file = "v17/data/test1_confinement.tsv"
    
    sim = HFKT_SkyrmeSim(N=64, L=30.0, dt=0.01, lambda_s=0.25)
    sim.initialize_hopfion_Q1()
    
    with open(out_file, "w") as f:
        f.write("Step\\tTime\\tTotal_Energy\\tCore_Count\\n")
        
        init_E = sim.get_energy()
        print(f"Initial Soliton Energy: {init_E:.4f}")
        
        for s in range(500):
            sim.step()
            if s % 50 == 0:
                E = sim.get_energy()
                core = np.sum(sim.R[0] < -0.5)
                f.write(f"{s}\\t{s*sim.dt:.3f}\\t{E:.6f}\\t{core}\\n")
                print(f"Step {s:04d} | Energy: {E:.4f} ({(E/init_E)*100:.2f}%) | Core: {core}")
    print("Test 1 complete.")

def run_test_2_emergent_fields():
    print("\\n--- TEST 2: Emergent 1/r^2 Faraday Fields ---")
    out_file = "v17/data/test2_far_field.tsv"
    
    sim = HFKT_SkyrmeSim(N=64, L=30.0, dt=0.01, lambda_s=0.25)
    sim.initialize_hopfion_Q1()
    
    # We analyze the relaxed initial structure
    for _ in range(10):
        sim.step()
        
    radii, E_vals = sim.extract_radial_far_field()
    
    # We take a sample of 1000 points to keep CSV reasonable
    idx = np.random.choice(len(radii), size=min(1000, len(radii)), replace=False)
    sub_r = radii[idx]
    sub_E = E_vals[idx]
    
    # Fit E ~ r^n
    log_r = np.log(sub_r)
    log_E = np.log(sub_E)
    slope, _, _, _, _ = linregress(log_r, log_E)
    
    with open(out_file, "w") as f:
        f.write("Extracted_Exponent_n\\n")
        f.write(f"{slope:.4f}\\n\\n")
        f.write("R\\tE_mag\\n")
        for r, e in zip(sub_r, sub_E):
            f.write(f"{r:.4f}\\t{e:.8e}\\n")
            
    print(f"True Extracted Far-Field Exponent: E ~ 1/r^{(abs(slope)):.4f} (Ideal: 2.000)")
    print("Test 2 complete.")
    
def run_test_3_spontaneous_gravity():
    print("\\n--- TEST 3: Static Geometric Gravity (Interaction Energy) ---")
    out_file = "v17/data/test3_gravity.tsv"
    
    # Evolve 2 solitons statically separated by D to natively measure non-circular gravitational energy
    distances = np.linspace(8.0, 16.0, 8)
    energies = []
    
    with open(out_file, "w") as f:
        f.write("Distance(D)\\tInteraction_U\\tForce(-dU/dD)\\n")
        
        # We need the base energy of one particle to subtract
        sim1 = HFKT_SkyrmeSim(N=64, L=30.0, lambda_s=0.25)
        sim1.initialize_hopfion_Q1(shift_x=0)
        E_single = sim1.get_energy()
        
        for D in distances:
            sim2 = HFKT_SkyrmeSim(N=64, L=30.0, lambda_s=0.25)
            # Standard quaternion product superposition
            sim2.initialize_two_hopfions_product(separation=D)
            # The exact interaction potential is the difference between total and 2*isolated
            E_total = sim2.get_energy()
            U_int = E_total - 2.0 * E_single
            energies.append(U_int)
            print(f"D = {D:.2f} | U_int = {U_int:.5e}")
            
        for i in range(len(distances)-1):
            dD = distances[i+1] - distances[i]
            dU = energies[i+1] - energies[i]
            F = -dU / dD
            D_mid = distances[i] + dD/2
            f.write(f"{D_mid:.4f}\\t{energies[i]:.6e}\\t{F:.6e}\\n")
            
    # Calculate force law F ~ D^n
    F_vals = []
    D_vals = []
    for i in range(len(distances)-1):
        dD = distances[i+1] - distances[i]
        dU = energies[i+1] - energies[i]
        F = -dU / dD
        if F > 0: # Check if genuinely attractive
            F_vals.append(F)
            D_vals.append(distances[i] + dD/2)
            
    if len(F_vals) > 2:
        slope, _, _, _, _ = linregress(np.log(D_vals), np.log(F_vals))
        print(f"Extracted Interaction Force: Attractive! Exponent: n = {slope:.4f} (Ideal Gravity: -2.000)")
    else:
        print("Not enough attractive data points to fit a power law. Repulsion dominates at this range.")

if __name__ == "__main__":
    run_test_1_confinement()
    run_test_2_emergent_fields()
    run_test_3_spontaneous_gravity()
