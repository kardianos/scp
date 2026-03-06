import numpy as np
import os
from scipy.stats import linregress

print("HFKT v15: Avenue 1 - Refraction Gravity")
print("=======================================")
print("Testing if a gradient in c_eff natively accelerates a resonant wave-trap.")

os.makedirs("v15/data", exist_ok=True)
out_tsv = "v15/data/sim_refraction_gravity_raw.tsv"

class RefractionGravitySim:
    def __init__(self, L0=0.5, f0=1e14, c_base=3e8):
        self.L0 = L0
        self.f0 = f0
        self.c_base = c_base
        
        # We model a 1D resonant trap (bounce chamber) placed in a spatial gradient of c
        # A massive body configures a non-local strain, causing c_eff to decrease closer to it.
        # c(x) = c_base * (1 - GM / (c^2 r)) ~ c_base * (1 - alpha/x)
        self.alpha = 10.0 # Strength of the refractive gradient

    def c_eff(self, x):
        """ The effective speed of light at position x.
            It drops slightly as x -> 0 (closer to the massive gravitational source) """
        return self.c_base * (1.0 - self.alpha / x)

    def simulate_refractive_acceleration(self, initial_x=100.0, t_max=1e-5):
        """ Tracks a wave bouncing between two structural boundaries of a trap
            positioned in a c_eff gradient. """
            
        t = 0.0
        # The trap has a rear boundary x_R and front boundary x_F.
        # It begins at rest: v_trap = 0.
        x_R = initial_x
        x_F = initial_x + self.L0
        v_trap = 0.0
        
        # A photon starts at the rear heading forward
        pos = x_R
        direction = 1
        
        bounces = 0
        time_bins = []
        trap_pos_bins = []
        trap_vel_bins = []
        
        # Unlike the v14 inertial mass where the trap was *forced* to accelerate uniformly,
        # here the trap is *free*. It accelerates based entirely on the net momentum 
        # transferred by the photon bouncing off the walls.
        # We track the trap's mass to calculate its velocity change.
        m_trap = (self.f0 / self.c_base**2) * 2.0 # Arbitrary structural mass + photon energy
        f_current = self.f0
        
        while t < t_max:
            # We must numerically integrate the photon's path because c varies with x.
            # dx/dt = direction * c(x)
            dt_step = 1e-12
            
            # Move the photon and the trap
            t += dt_step
            pos += direction * self.c_eff(pos) * dt_step
            
            x_R += v_trap * dt_step
            x_F += v_trap * dt_step
            
            hit = False
            if direction == 1 and pos >= x_F:
                hit = True
            elif direction == -1 and pos <= x_R:
                hit = True
                
            if hit:
                # Photon hits a boundary. 
                # Energy E = hf. Momentum p = E/c. 
                # But c is local! p = hf / c_local
                c_local = self.c_eff(pos)
                p_in = (f_current / c_local) * direction
                
                # Relativistic Doppler shift off moving boundary
                beta = v_trap / c_local
                if direction == 1:
                    # Hitting front -> chasing it.
                    f_new = f_current * np.sqrt((1 - beta)/(1 + beta)) if beta > -1 else f_current
                    p_out = -(f_new / c_local)
                else:
                    # Hitting rear -> head on if trap moving forward, or chasing if trap moving backward.
                    # Since gravity pulls trap left (backward), rear boundary is running away from photon.
                    f_new = f_current * np.sqrt((1 + beta)/(1 - beta)) if beta < 1 else f_current
                    p_out = (f_new / c_local)
                    
                dp_photon = p_out - p_in
                # Conservation of momentum: trap absorbs -dp_photon
                # Wait: p_in applies force in its direction.
                # dp_trap = p_in - p_out
                dp_trap = p_in - p_out
                
                # Update trap velocity v = p / m
                v_trap += dp_trap / m_trap
                
                f_current = f_new
                direction *= -1
                bounces += 1
                
                if bounces % 10 == 0:
                    time_bins.append(t)
                    trap_pos_bins.append(x_R)
                    trap_vel_bins.append(v_trap)
                    
        # Calculate resulting macroscopic acceleration of the trap
        # a = dv / dt
        if len(trap_vel_bins) > 1:
            dt = time_bins[-1] - time_bins[0]
            dv = trap_vel_bins[-1] - trap_vel_bins[0]
            avg_accel = dv / dt
        else:
            avg_accel = 0.0
            
        return time_bins, trap_pos_bins, avg_accel

    def run_sweep(self):
        print(f"Streaming Refraction Gravity derivations to {out_tsv}...")
        
        # Place the trap at various distances from the gravitational source
        distances = [100.0, 200.0, 300.0, 400.0, 500.0]
        results = []
        
        with open(out_tsv, "w") as f:
            f.write("Distance(X)\tInduced_Acceleration(a)\n")
            
            for D in distances:
                print(f"Simulating trap released at X = {D} ...")
                t_bins, x_bins, a = self.simulate_refractive_acceleration(initial_x=D, t_max=0.5e-5)
                f.write(f"{D}\t{a:e}\n")
                results.append((D, a))
                
        # Analyze the acceleration curve: a ~ 1/X^n
        X_vals = np.array([r[0] for r in results])
        # Acceleration should be negative (pulling toward source at X=0)
        a_mags = np.abs(np.array([r[1] for r in results]))
        
        log_X = np.log(X_vals)
        log_a = np.log(a_mags)
        slope, intercept, _, _, _ = linregress(log_X, log_a)
        
        print("\nSweep Complete.")
        for r in results:
            print(f"Distance = {r[0]:.1f} | Induced a = {r[1]:e} m/s^2")
        print(f"Extracted Gravity Law: a ~ 1/X^{abs(slope):.4f} (Newtonian is 2.0000)")

if __name__ == "__main__":
    sim = RefractionGravitySim()
    sim.run_sweep()
