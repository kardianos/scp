import numpy as np
import os

print("HFKT v14: Avenue 2 - Origin of Inertial Mass")
print("============================================")
print("Deriving F=ma from internal Doppler-shifted radiation pressure.")

os.makedirs("v14/data", exist_ok=True)
out_tsv = "v14/data/sim_inertial_mass_raw.tsv"

class InertialMassSim:
    def __init__(self, L0=1.0, f0=1e15, c=3e8):
        self.L0 = L0          # Rest length of the resonant trap
        self.f0 = f0          # Base trapped frequency
        self.c = c            # Speed of light
        self.E_rest = f0      # Proportional to rest energy (h*f0)
        self.m_rest = self.E_rest / (self.c**2)

    def simulate_acceleration(self, a, t_max=1e-5):
        t = 0.0
        pos = 0.0 
        f_current = self.f0
        direction = 1 # 1: forward, -1: backward
        
        net_momentum_transfer = 0.0
        time_bins = []
        force_bins = []
        bounces = 0
        last_t = 0.0
        last_p = 0.0
        
        while t < t_max:
            # We solve 0.5 * a * t_new^2 - v_p * t_new + (offset) = 0
            # v_p = c if dir=1, -c if dir=-1
            v_p = self.c * direction
            
            # The photon equation: x(t_new) = pos + v_p * (t_new - t)
            # The boundary equation: 
            # If dir = 1 (hitting front): x_F(t_new) = self.L0 + 0.5 * a * t_new^2
            # If dir = -1 (hitting rear): x_R(t_new) = 0.5 * a * t_new^2
            
            if direction == 1:
                C = self.L0 - pos + v_p * t
            else:
                C = 0.0 - pos + v_p * t
                
            # 0.5 * a * t_new^2 - v_p * t_new + C = 0
            A = 0.5 * a
            B = -v_p
            
            # Quadratic formula: t_new = (-B +/- sqrt(B^2 - 4AC)) / 2A
            disc = B**2 - 4*A*C
            if disc < 0:
                break # Light can't reach boundary (accelerating too fast - unruh horizon)
            
            if a != 0:
                t1 = (-B + np.sqrt(disc)) / (2*A)
                t2 = (-B - np.sqrt(disc)) / (2*A)
                # We want the smallest t_new > t
                t_new = t2 if (t2 > t and t2 < t1) or t1 <= t else t1
            else:
                t_new = -C / B
            
            if t_new <= t:
                break
                
            t = t_new
            pos = 0.5 * a * t**2 + (self.L0 if direction == 1 else 0.0)
            v_trap = a * t
            beta = v_trap / self.c
            
            # Momentum p = E/c = hf/c. We track p directly.
            # Initial momentum striking the boundary
            p_in = (f_current / self.c) * direction 
            
            if direction == 1:
                # Hitting front -> Redshift 
                f_new = f_current * np.sqrt((1 - beta)/(1 + beta))
                p_out = -(f_new / self.c)
            else:
                # Hitting rear -> Blueshift
                f_new = f_current * np.sqrt((1 + beta)/(1 - beta))
                p_out = (f_new / self.c)
                
            # Force exerted BY photon ON boundary: dp = p_in - p_out
            dp = p_in - p_out
            net_momentum_transfer += dp
            f_current = f_new
            direction *= -1
            bounces += 1
            
            if bounces % 100 == 0:
                dt_bin = t - last_t
                dp_bin = net_momentum_transfer - last_p
                F_avg = dp_bin / dt_bin
                time_bins.append(t)
                force_bins.append(F_avg)
                last_t = t
                last_p = net_momentum_transfer
                
        # We need the average force opposing the acceleration
        avg_resisting_force = -np.mean(force_bins) if len(force_bins) > 0 else 0.0
        return time_bins, force_bins, avg_resisting_force

    def run_accleration_sweep(self):
        print(f"Streaming Inertial mass derivation to {out_tsv}...")
        
        # Test linear relationship F = ma
        accelerations = [1e4, 5e4, 1e5, 5e5, 1e6]
        results = []
        
        with open(out_tsv, "w") as f:
            f.write("Acceleration(a)\tNet_Resisting_Force(F)\tTheoretical_m(E/c^2)\tCalculated_m(F/a)\tError_Margin\n")
            
            for a in accelerations:
                print(f"Simulating trap under acceleration a = {a:e} m/s^2 ...")
                t_bins, F_bins, F_net = self.simulate_acceleration(a, t_max=1e-5)
                
                # The net force calculated is the force exerted BY the wave.
                # Inertial resistance is the negative of that. F_resisting = - F_net
                F_resisting = -F_net
                
                calculated_m = F_resisting / a
                error = abs((calculated_m - self.m_rest) / self.m_rest) * 100
                
                f.write(f"{a:e}\t{F_resisting:e}\t{self.m_rest:e}\t{calculated_m:e}\t{error:.4f}%\n")
                results.append((a, F_resisting, calculated_m, error))
                
        print("\nAcceleration Sweep Complete.")
        for r in results:
            print(f"a = {r[0]:e} | F = {r[1]:e} | m_calc = {r[2]:e} | Error: {r[3]:.4f}%")
        print(f"Theoretical Rest Mass (E/c^2): {self.m_rest:e}")

if __name__ == "__main__":
    sim = InertialMassSim(L0=0.1, f0=1e14, c=3e8) # Scaled to allow more bounces in short dt
    sim.run_accleration_sweep()
