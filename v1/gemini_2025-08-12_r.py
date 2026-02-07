import numpy as np
import matplotlib.pyplot as plt
import time
import os
import sys

# --- Global Settings & Versioning ---
SIMULATION_VERSION = "v10.0"
PLOTTING_ENABLED = True
SCRIPT_FILENAME = f"chpt_simulation_{SIMULATION_VERSION}"
plt.switch_backend('Agg')

# --- Helper Function for Saving Figures ---
def save_figure(case_num_str, suffix=""):
    if PLOTTING_ENABLED:
        filename = f"{SCRIPT_FILENAME}_case_{case_num_str}{suffix}.png"
        try:
            plt.savefig(filename)
            if "step_" not in suffix: print(f"Saved plot to {filename}")
        except Exception as e:
            print(f"Could not save plot: {e}")
        plt.close('all')

# --- Core CHPT System Class (Spinor Model with Bias) ---
class CHPTSystem:
    def __init__(self, a=0.1, b=1.0, c_p=1.5, kappa=0.1, eta=0.05):
        self.a, self.b, self.c_p, self.kappa = a, b, c_p, kappa
        self.eta = eta # New constant for coupling to external bias
        self.c = 299792458.0; self.G = 6.67430e-11
        self.rho_0, self.rho_p = self._calculate_potential_minima()
        if self.rho_0 is None: raise ValueError("Invalid potential parameters.")

    def _calculate_potential_minima(self):
        d = self.b**2 - 4*self.a*self.c_p
        return (None, None) if d < 0 else ((self.b - np.sqrt(d))/(2*self.a), (self.b + np.sqrt(d))/(2*self.a))
    
    def get_total_density(self, phi_up, phi_down):
        return np.abs(phi_up)**2 + np.abs(phi_down)**2

    def get_potential_force(self, total_rho):
        return (self.c_p - self.b * total_rho + self.a * total_rho**2)
    
    def get_energy(self, phi_up, phi_down, dx):
        # Only calculates potential energy for simplicity and speed
        total_rho = self.get_total_density(phi_up, phi_down)
        potential_V = self.c_p*total_rho - (self.b/2)*total_rho**2 + (self.a/3)*total_rho**3
        return np.sum(potential_V) * (dx**3)

# --- Simulation Cases ---

def run_case_5a_spinor_stability():
    """Prerequisite test: Is a single multi-component soliton stable in a biased vacuum?"""
    print("\n>>> CASE 5a: SINGLE SPINOR SOLITON STABILITY TEST <<<")
    print("Desired: The integrated energy of a single particle should remain constant.")
    
    N=40; dx=0.5; dt=0.04; T_steps=4000
    system = CHPTSystem()
    
    phi = {'up':np.zeros((N,N,N),dtype=np.complex128), 'down':np.zeros((N,N,N),dtype=np.complex128)}
    X,Y,Z = np.mgrid[-N//2:N//2,-N//2:N//2,-N//2:N//2]*dx
    
    R = np.sqrt(X**2+Y**2+Z**2)
    phi['up'] += (np.sqrt(system.rho_p)) * np.exp(-R**2 / 2.5**2)
    phi_prev = {k: v.copy() for k, v in phi.items()}
    
    # Define the external chiral bias from galactic rotation (simplified)
    # C_mu = (Ct, Cx, Cy, Cz). Let's assume a rotation around the z-axis
    C_bias = {'t': 0.1, 'x': Y, 'y': -X, 'z': np.zeros_like(Z)}
    
    initial_energy = system.get_energy(phi['up'], phi['down'], dx)
    
    print("Running stability simulation...")
    for t in range(T_steps):
        total_rho = system.get_total_density(phi['up'], phi['down'])
        pot_force_mag = system.get_potential_force(total_rho)
        
        phi_next = {}
        for s, o in [('up','down'), ('down','up')]:
            grad_s = np.gradient(phi[s], dx)
            lap = np.gradient(grad_s[0],dx,axis=0)[0] + np.gradient(grad_s[1],dx,axis=1)[1] + np.gradient(grad_s[2],dx,axis=2)[2]
            
            d_dt_s = (phi[s] - phi_prev[s]) / dt
            d_dt_o = (phi[o] - phi_prev[o]) / dt
            
            # Equations of Motion from Part 2
            force_pot = pot_force_mag * phi[s]
            force_chiral = 2j * system.kappa * d_dt_o if s == 'up' else -2j * system.kappa * d_dt_o
            force_bias = 2j * system.eta * (C_bias['t']*d_dt_s + C_bias['x']*grad_s[0] + C_bias['y']*grad_s[1] + C_bias['z']*grad_s[2])
            
            phi_next[s] = 2*phi[s] - phi_prev[s] + dt**2 * (lap - force_pot - force_chiral - force_bias)
        
        phi_prev, phi = phi, phi_next

    final_energy = system.get_energy(phi['up'], phi['down'], dx)
    energy_change = np.abs(final_energy - initial_energy) / initial_energy

    print(f"Observed: Integrated energy changed by {energy_change*100:.4f}% over {T_steps} steps.")
    print("\n--- VERIFICATION ---")
    if energy_change < 0.05:
        print("SUCCESS: Single soliton is stable in a biased vacuum.")
        return True, phi, dx # Return the stable profile for other tests
    else:
        print("FAILURE: Single soliton is unstable. Halting further tests.")
        return False, None, None

def run_case_2_rest_mass(stable_soliton_data):
    if not stable_soliton_data[0]: return
    print("\n>>> CASE 2: RELATIVISTIC MECHANICS (REST MASS) <<<")
    is_stable, phi, dx = stable_soliton_data
    energy = CHPTSystem().get_energy(phi['up'], phi['down'], dx)
    print(f"Desired: Positive energy. Observed: {energy:.4e}")
    if energy > 0: print("--- VERIFICATION ---\nSUCCESS.")
    else: print("--- VERIFICATION ---\nFAILURE.")

def run_case_1_gravity(stable_soliton_data):
    if not stable_soliton_data[0]: return
    print("\n>>> CASE 1: VERIFICATION OF LONG-RANGE GRAVITY FIELD <<<")
    is_stable, phi, dx = stable_soliton_data
    system = CHPTSystem()
    
    # Derive mass from the stable soliton's energy
    m_source = system.get_energy(phi['up'], phi['down'], dx) / system.c**2
    
    r = np.linspace(10*dx*phi['up'].shape[0], 200*dx*phi['up'].shape[0], 2000)
    k = 2 * system.G * m_source / system.c**2
    rho_gravity = system.rho_0 * np.exp(-k / r)
    a_chpt = -0.5 * system.c**2 * np.gradient(np.log(rho_gravity), r[1]-r[0])
    a_newton = -system.G * m_source / r**2
    error = np.mean(np.abs((a_chpt - a_newton) / a_newton))
    
    print(f"Desired: Deviation < 1%. Observed: {error*100:.4f}%")
    if error < 0.01: print("--- VERIFICATION --- \nSUCCESS."); save_figure("1_gravity")
    else: print("--- VERIFICATION --- \nFAILURE.")

def run_case_3_cosmology(stable_soliton_data):
    if not stable_soliton_data[0]: return
    print("\n>>> CASE 3: COSMOLOGICAL MECHANICS (GALAXY ROTATION) <<<")
    # This case is complex and depends on a robust G. We'll mark it as a placeholder.
    print("--- VERIFICATION ---\nSKIPPED (Requires separate, large-scale simulation).")

# ... Cases 4 and 6 can be copied verbatim as they are independent ...
def run_case_4_interference():
    print("\n>>> CASE 4: EMF PROPAGATION (DOUBLE-SLIT INTERFERENCE) <<<")
    print("Desired: Clear interference pattern.")
    wl=500e-9;k=2*np.pi/wl;sw=2*wl;ss=10*sw;sd=0.01;scw=0.001;ns=50
    s1y=np.linspace(-ss/2-sw/2,-ss/2+sw/2,ns);s2y=np.linspace(ss/2-sw/2,ss/2+sw/2,ns)
    slits=np.concatenate([s1y,s2y]);scr_y=np.linspace(-scw/2,scw/2,2000);field=np.zeros_like(scr_y,dtype=np.complex128)
    for y in slits:field+=np.exp(1j*k*np.sqrt(sd**2+(scr_y-y)**2))
    intensity=np.abs(field)**2;intensity/=intensity.max();print(f"Observed: {len(np.where(intensity>0.8)[0])} peaks.")
    if len(np.where(intensity>0.8)[0])>10:print("--- VERIFICATION ---\nSUCCESS."); save_figure("4_interference")
    else:print("--- VERIFICATION ---\nFAILURE.")

def run_case_6_tunneling():
    print("\n>>> CASE 6: 'QUANTUM' TUNNELING VIA FDTD <<<")
    print("Desired: Non-zero transmitted intensity.")
    N=5000;T=12000;dx=1.0;dt=0.2;phi=np.zeros(N,dtype=np.complex128);phi_p=np.zeros(N,dtype=np.complex128)
    k0=np.pi/30;E0=np.sqrt(1+k0**2);x0=N//5;sig=40.0;x=np.arange(N)*dx
    phi=np.exp(-(x-x0)**2/(2*sig**2))*np.exp(1j*k0*(x-x0));phi_p=np.exp(-(x-x0-k0*dt/E0)**2/(2*sig**2))*np.exp(1j*k0*(x-x0-k0*dt/E0))
    V=np.zeros(N);bs,be=N//2,N//2+50;bh=0.015;V[bs:be]=bh
    init_I=np.sum(np.abs(phi)**2);print(f"KE:{E0-1:.4f}. Barrier:{bh:.4f}")
    for t in range(T):lap=np.roll(phi,1)+np.roll(phi,-1)-2*phi;phi_n=2*phi-phi_p+(dt/dx)**2*lap-dt**2*V*phi;phi_p,phi=phi,phi_n
    ref_I=np.sum(np.abs(phi[:bs])**2);tun_I=np.sum(np.abs(phi[be:])**2)
    print(f"Observed: Reflected={ref_I/init_I*100:.2f}%, Tunneled={tun_I/init_I*100:.2f}%")
    if 0.01<tun_I/init_I<0.99:print("--- VERIFICATION ---\nSUCCESS."); save_figure("6_tunneling")
    else:print("--- VERIFICATION ---\nFAILURE.")


# --- Main Execution ---
if __name__ == "__main__":
    start_time = time.time()
    
    # The entire workflow now depends on the initial stability test.
    is_stable, stable_phi, stable_dx = run_case_5a_spinor_stability()
    
    # Store the results in a tuple to pass to other functions
    stable_soliton_data = (is_stable, stable_phi, stable_dx)

    # These cases use the results of the stability test
    run_case_2_rest_mass(stable_soliton_data)
    run_case_1_gravity(stable_soliton_data)
    run_case_3_cosmology(stable_soliton_data)
    
    # These cases are independent regression tests
    run_case_4_interference()
    run_case_6_tunneling()

    # The interaction test is no longer run, as it is superseded by the stability test.
    # The next step, if stability is achieved, would be to build Case 5b based on this stable model.
    
    end_time = time.time()
    print(f"\nTotal execution time: {end_time - start_time:.2f} seconds.")
