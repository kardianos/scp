import numpy as np
import matplotlib.pyplot as plt
import time
import os
import sys

# --- Global Settings & Versioning ---
SIMULATION_VERSION = "v13.1"
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

# --- Core CHPT System Class ---
class CHPTSystem:
    def __init__(self, a=0.1, b=1.0, c_p=1.5, kappa_base=0.1, j_stable=0.1, beta=10.0):
        self.a, self.b, self.c_p = a, b, c_p
        self.kappa_base, self.j_stable, self.beta = kappa_base, j_stable, beta
        
        # --- FIX: Re-instated the missing calculation ---
        self.rho_0, self.rho_p = self._calculate_potential_minima()
        if self.rho_0 is None: raise ValueError("Invalid potential parameters.")
        # --- END FIX ---

    def _calculate_potential_minima(self):
        d = self.b**2 - 4*self.a*self.c_p
        return (None, None) if d < 0 else ((self.b - np.sqrt(d))/(2*self.a), (self.b + np.sqrt(d))/(2*self.a))

    def get_potential_force_magnitude(self, total_rho):
        return (self.c_p - self.b * total_rho + self.a * total_rho**2)

    def get_energy(self, phi, pi, K2, dx):
        total_rho = np.abs(phi['up'])**2 + np.abs(phi['down'])**2
        V_density = self.c_p*total_rho - (self.b/2)*total_rho**2 + (self.a/3)*total_rho**3
        
        phi_k_up = np.fft.fftn(phi['up']); phi_k_down = np.fft.fftn(phi['down'])
        grad_energy = np.sum(K2 * (np.abs(phi_k_up)**2 + np.abs(phi_k_down)**2)) / len(K2.ravel())
        pi_energy = np.sum(np.abs(pi['up'])**2 + np.abs(pi['down'])**2)
        
        j0 = 1j * (np.conj(phi['up'])*pi['down'] - np.conj(phi['down'])*pi['up'])
        V_chiral = (self.beta / 2.0) * (np.real(j0) - self.j_stable)**2

        total_energy = np.sum(V_density + V_chiral) + grad_energy + pi_energy
        return total_energy * (dx**3)

# --- The Definitive Stability Test ---
def run_case_5a_stability_v13():
    print("\n>>> CASE 5a: STABILITY TEST (v13 IMPLICIT SOLVER) <<<")
    print("Desired: The integrated energy must remain constant.")
    
    N=32; dx=0.5; dt=0.01; T_steps=10000
    system = CHPTSystem(j_stable=0.05, beta=50.0)
    
    phi = {'up': np.zeros((N,N,N), dtype=np.complex128), 'down': np.zeros((N,N,N), dtype=np.complex128)}
    pi = {'up': np.zeros_like(phi['up']), 'down': np.zeros_like(phi['down'])}
    
    X,Y,Z = np.mgrid[-N//2:N//2, -N//2:N//2, -N//2:N//2] * dx
    R = np.sqrt(X**2+Y**2+Z**2)
    profile = np.sqrt(system.rho_p) * np.exp(-R**2 / 2.5**2)
    phi['up'] += profile
    pi['down'] = 1j * system.j_stable * phi['up'] 
    
    kx=2*np.pi*np.fft.fftfreq(N,d=dx);ky=2*np.pi*np.fft.fftfreq(N,d=dx);kz=2*np.pi*np.fft.fftfreq(N,d=dx)
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij'); K2 = KX**2 + KY**2 + KZ**2
    
    energies = []
    print("Running stability simulation with Implicit Operator Splitting solver...")
    for t in range(T_steps):
        
        # KICK 1 (dt/2)
        total_rho = np.abs(phi['up'])**2 + np.abs(phi['down'])**2
        force_mag = system.get_potential_force_magnitude(total_rho)
        pi['up'] -= (force_mag * phi['up']) * (dt / 2.0)
        pi['down'] -= (force_mag * phi['down']) * (dt / 2.0)
        
        j0 = 1j * (np.conj(phi['up'])*pi['down'] - np.conj(phi['down'])*pi['up'])
        kappa_eff = system.kappa_base * np.maximum(0, 1 - system.beta * (np.real(j0) - system.j_stable)**2)
        
        C_up = 2j * kappa_eff * dt / 2.0; C_down = -2j * kappa_eff * dt / 2.0
        det = 1 - C_up * C_down
        pi_up_next = (pi['up'] - C_up * pi['down']) / det
        pi_down_next = (pi['down'] - C_down * pi_up_next) # Denominator is 1 after substitution
        pi['up'], pi['down'] = pi_up_next, pi_down_next
        
        # DRIFT (dt)
        for s in ['up', 'down']:
            phi_k = np.fft.fftn(phi[s]); pi_k = np.fft.fftn(pi[s])
            k_abs = np.sqrt(K2)
            phi_k_next = phi_k*np.cos(k_abs*dt) + pi_k*np.sin(k_abs*dt)/(k_abs+1e-9)
            pi_k_next = -phi_k*k_abs*np.sin(k_abs*dt) + pi_k*np.cos(k_abs*dt)
            phi[s] = np.fft.ifftn(phi_k_next); pi[s] = np.fft.ifftn(pi_k_next)

        # KICK 2 (dt/2)
        total_rho = np.abs(phi['up'])**2 + np.abs(phi['down'])**2
        force_mag = system.get_potential_force_magnitude(total_rho)
        pi['up'] -= (force_mag * phi['up']) * (dt / 2.0)
        pi['down'] -= (force_mag * phi['down']) * (dt / 2.0)
        j0 = 1j * (np.conj(phi['up'])*pi['down'] - np.conj(phi['down'])*pi['up'])
        kappa_eff = system.kappa_base * np.maximum(0, 1 - system.beta * (np.real(j0) - system.j_stable)**2)
        C_up = 2j*kappa_eff*dt/2.0; C_down = -2j*kappa_eff*dt/2.0
        det = 1 - C_up * C_down
        pi_up_next = (pi['up'] - C_up * pi['down']) / det
        pi_down_next = (pi['down'] - C_down * pi_up_next)
        pi['up'], pi['down'] = pi_up_next, pi_down_next
        
        if t % 200 == 0:
            current_energy = system.get_energy(phi, pi, K2, dx)
            energies.append(current_energy)
            print(f"Step {t}/{T_steps}, Total Energy: {current_energy:.5e}")

    energies = np.array(energies)
    initial_energy = energies[0]
    final_energy = energies[-1]
    energy_change_percent = np.abs(final_energy - initial_energy) / initial_energy * 100

    print(f"Observed: Integrated energy changed by {energy_change_percent:.6f}% over {T_steps} steps.")
    print("\n--- VERIFICATION ---")
    if energy_change_percent < 0.01 and np.isfinite(final_energy):
        print("SUCCESS: Single spinor soliton is STABLE. The topological tension model and implicit solver are valid.")
    else:
        print("FAILURE: Soliton remains UNSTABLE. The model is still flawed.")

    if PLOTTING_ENABLED:
        plt.figure(figsize=(10, 6))
        plt.plot(np.linspace(0, T_steps, len(energies)), energies)
        plt.title("Case 5a: Energy Conservation (v13.1 Implicit Solver)")
        plt.xlabel("Time Step"); plt.ylabel("Total System Energy")
        plt.grid(True);
        if energies.size > 0 and np.all(np.isfinite(energies)):
             plt.ylim(energies.min()*(1-0.001*np.sign(energies.min())), energies.max()*(1+0.001*np.sign(energies.max())))
        save_figure("5a_stability")

# --- Main Execution ---
if __name__ == "__main__":
    start_time = time.time()
    run_case_5a_stability_v13()
    end_time = time.time()
    print(f"\nTotal execution time: {end_time - start_time:.2f} seconds.")
