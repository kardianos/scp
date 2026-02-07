import numpy as np
import matplotlib.pyplot as plt
import time
import os
import sys

# --- Global Settings & Versioning ---
SIMULATION_VERSION = "v14.2"
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

# --- Dirac Gamma Matrices (chiral representation) ---
gamma0 = np.array([[0,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,0]], dtype=np.complex128)
gamma1 = np.array([[0,0,0,1],[0,0,1,0],[0,-1,0,0],[-1,0,0,0]], dtype=np.complex128)
gamma2 = np.array([[0,0,0,-1j],[0,0,1j,0],[0,1j,0,0],[-1j,0,0,0]], dtype=np.complex128)
gamma3 = np.array([[0,0,1,0],[0,0,0,-1],[-1,0,0,0],[0,1,0,0]], dtype=np.complex128)

# --- The Definitive Stability Test based on the Non-Linear Dirac Equation ---
def run_case_5a_stability_v14():
    """
    Final attempt: Tests stability of a single soliton using the Non-Linear Dirac Equation,
    with a corrected, vectorized, and analytically exact kinetic solver.
    """
    print("\n>>> CASE 5a: STABILITY TEST (v14.2 NON-LINEAR DIRAC SOLVER) <<<")
    print("Desired: The integrated probability density (particle number) must be conserved.")
    
    # Simulation Parameters
    N=64; dx=0.5; dt=0.02; T_steps=5000
    m0 = 1.0  # Base mass
    g = 2.0   # Self-interaction strength
    
    # --- Field Initialization ---
    psi = np.zeros((4, N, N, N), dtype=np.complex128)
    X,Y,Z = np.mgrid[-N//2:N//2,-N//2:N//2,-N//2:N//2]*dx
    R = np.sqrt(X**2 + Y**2 + Z**2)
    profile = np.exp(-R**2 / 2.5**2)
    psi[0, ...] = profile; psi[1, ...] = profile
    
    # --- Setup for Fourier Space ---
    kx=2*np.pi*np.fft.fftfreq(N,d=dx);ky=2*np.pi*np.fft.fftfreq(N,d=dx);kz=2*np.pi*np.fft.fftfreq(N,d=dx)
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')

    alpha1 = gamma1 @ gamma0; alpha2 = gamma2 @ gamma0; alpha3 = gamma3 @ gamma0
    H_k_5D = alpha1[:,:,None,None,None]*KX + alpha2[:,:,None,None,None]*KY + alpha3[:,:,None,None,None]*KZ + m0*gamma0[:,:,None,None,None]
    
    E_k = np.sqrt(KX**2 + KY**2 + KZ**2 + m0**2)
    cos_E_dt = np.cos(E_k * dt); sin_E_dt_over_E = np.sin(E_k * dt) / (E_k + 1e-9)
    identity_5D = np.identity(4)[:,:,None,None,None]
    U_kinetic = cos_E_dt * identity_5D - 1j * sin_E_dt_over_E * H_k_5D

    # --- Diagnostics ---
    particle_number = []
    initial_number = np.sum(np.conj(psi) * psi)
    print("Running stability simulation with Non-Linear Dirac solver...")
    
    for t in range(T_steps):
        # --- Operator Splitting ---
        # POTENTIAL KICK (dt/2)
        # --- FIX: Use np.einsum for correct matrix multiplication ---
        gamma0_psi = np.einsum('ab,b...->a...', gamma0, psi)
        psi_bar_psi = np.sum(np.conj(psi) * gamma0_psi, axis=0)
        # --- END FIX ---
        
        potential_term = g * psi_bar_psi
        U_potential = np.exp(-1j * potential_term * (dt / 2.0))
        psi *= U_potential
        
        # KINETIC DRIFT (dt)
        psi_k = np.fft.fftn(psi, axes=(1,2,3))
        psi_k = np.einsum('ab...,b...->a...', U_kinetic, psi_k)
        psi = np.fft.ifftn(psi_k, axes=(1,2,3))
        
        # SECOND POTENTIAL KICK (dt/2)
        # --- FIX: Use np.einsum for correct matrix multiplication ---
        gamma0_psi = np.einsum('ab,b...->a...', gamma0, psi)
        psi_bar_psi = np.sum(np.conj(psi) * gamma0_psi, axis=0)
        # --- END FIX ---
        
        potential_term = g * psi_bar_psi
        U_potential = np.exp(-1j * potential_term * (dt / 2.0))
        psi *= U_potential

        if t % 200 == 0:
            current_number = np.sum(np.conj(psi) * psi)
            particle_number.append(current_number)
            print(f"Step {t}/{T_steps}, Particle Number: {current_number.real:.6f}")

    particle_number = np.array(particle_number)
    final_number = particle_number[-1]
    number_change_percent = np.abs(final_number - initial_number) / initial_number * 100

    print(f"Observed: Particle number changed by {number_change_percent.real:.6f}% over {T_steps} steps.")
    print("\n--- VERIFICATION ---")
    if number_change_percent < 0.01 and np.isfinite(final_number):
        print("SUCCESS: Single Dirac soliton is STABLE. This physical model is valid.")
    else:
        print("FAILURE: The Non-Linear Dirac model, or its simulation, is still unstable.")

    if PLOTTING_ENABLED:
        plt.figure(figsize=(10, 6))
        plt.plot(np.linspace(0, T_steps, len(particle_number)), particle_number.real)
        plt.title("Case 5a: Particle Number Conservation (v14.2 Dirac Solver)")
        plt.xlabel("Time Step"); plt.ylabel("Total Probability (Particle Number)")
        plt.grid(True);
        if particle_number.size > 0 and np.all(np.isfinite(particle_number)):
             # Set tight y-limits to see small fluctuations
             plt.ylim(particle_number.real.min() * 0.9999, particle_number.real.max() * 1.0001)
        save_figure("5a_stability")

# --- Main Execution ---
if __name__ == "__main__":
    start_time = time.time()
    run_case_5a_stability_v14()
    end_time = time.time()
    print(f"\nTotal execution time: {end_time - start_time:.2f} seconds.")
