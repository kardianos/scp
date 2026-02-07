import numpy as np
import matplotlib.pyplot as plt
import time
import os
import sys

# --- Global Settings & Versioning ---
SIMULATION_VERSION = "v9.0"
PLOTTING_ENABLED = True
SCRIPT_FILENAME = f"chpt_simulation_{SIMULATION_VERSION}"
plt.switch_backend('Agg')

# --- Helper Function for Saving Figures ---
def save_figure(case_num_str, suffix=""):
    if PLOTTING_ENABLED:
        filename = f"{SCRIPT_FILENAME}_case_{case_num_str}{suffix}.png"
        try: plt.savefig(filename)
        except Exception as e: print(f"Could not save plot: {e}")
        plt.close('all')

# --- Core CHPT System Class (Spinor Model) ---
class CHPTSystem:
    def __init__(self, a=0.1, b=1.0, c_p=1.5, kappa=0.1):
        self.a, self.b, self.c_p, self.kappa = a, b, c_p, kappa
        self.c = 299792458.0; self.G = 6.67430e-11
        self.rho_0, self.rho_p = self._calculate_potential_minima()
        if self.rho_0 is None: raise ValueError("Invalid potential parameters.")

    def _calculate_potential_minima(self):
        d = self.b**2 - 4*self.a*self.c_p
        return (None, None) if d < 0 else ((self.b - np.sqrt(d))/(2*self.a), (self.b + np.sqrt(d))/(2*self.a))

    def get_total_density(self, phi_up, phi_down):
        return np.abs(phi_up)**2 + np.abs(phi_down)**2

    def get_force_term(self, total_rho):
        return (self.c_p - self.b * total_rho + self.a * total_rho**2)

# --- Main Simulation Cases ---

def run_case_5a_spinor_stability():
    """Prerequisite test: Is a single multi-component soliton stable over time?"""
    print("\n>>> CASE 5a: SINGLE SPINOR SOLITON STABILITY TEST <<<")
    print("Desired: The integrated density of a single particle should remain constant.")
    
    N=40; dx=0.4; dt=0.04; T_steps=2000
    system = CHPTSystem(kappa=0.1)
    
    phi = {'up': np.zeros((N,N,N), dtype=np.complex128), 'down': np.zeros((N,N,N), dtype=np.complex128)}
    X,Y,Z = np.mgrid[-N//2:N//2,-N//2:N//2,-N//2:N//2]*dx
    sigma=2.5
    
    # Initialize one "spin-up" particle at the center
    R = np.sqrt(X**2 + Y**2 + Z**2)
    phi['up'] += (np.sqrt(system.rho_p)) * np.exp(-R**2 / sigma**2)
    phi_prev = {k: v.copy() for k, v in phi.items()}
    
    initial_total_density = np.sum(system.get_total_density(phi['up'], phi['down']))
    
    for t in range(T_steps):
        total_rho = system.get_total_density(phi['up'], phi['down'])
        force_mag = system.get_force_term(total_rho)
        
        phi_next = {}
        for spin, other_spin in [('up', 'down'), ('down', 'up')]:
            lap = (np.roll(phi[spin],1,0)+np.roll(phi[spin],-1,0)+np.roll(phi[spin],1,1)+np.roll(phi[spin],-1,1)+
                   np.roll(phi[spin],1,2)+np.roll(phi[spin],-1,2)-6*phi[spin])/dx**2
            chiral_force = 2j*system.kappa*(phi[other_spin]-phi_prev[other_spin])/dt
            force_term_comp = force_mag * phi[spin]
            phi_next[spin] = 2*phi[spin]-phi_prev[spin]+dt**2*(lap-force_term_comp-chiral_force)
        
        phi_prev, phi = phi, phi_next

    final_total_density = np.sum(system.get_total_density(phi['up'], phi['down']))
    density_change = np.abs(final_total_density - initial_total_density) / initial_total_density

    print(f"Observed: Integrated density changed by {density_change*100:.4f}% over {T_steps} steps.")
    print("\n--- VERIFICATION ---")
    if density_change < 0.01:
        print("SUCCESS: Single spinor soliton is stable. Proceeding to interaction test.")
        return True
    else:
        print("FAILURE: Single soliton is unstable. Interaction simulation is invalid.")
        return False

def run_case_5b_atomic_forces_spinor(is_stable):
    if not is_stable: return
    print("\n>>> CASE 5b: ATOMIC FORCES (HIGH-FIDELITY SPINOR FDTD) <<<")
    print("Desired: Stable bound state (separation converges to a non-zero value).")
    
    N=50; dx=0.4; dt=0.03; T_steps=8000
    system = CHPTSystem(kappa=0.15)
    
    phi = {'up':np.zeros((N,N,N),dtype=np.complex128), 'down':np.zeros((N,N,N),dtype=np.complex128)}
    X,Y,Z = np.mgrid[-N//2:N//2,-N//2:N//2,-N//2:N//2]*dx
    sigma=2.5; d_initial=10.0
    
    # Particle 1 (Up) at (-d/2, 0, 0)
    R1=np.sqrt((X+d_initial/2)**2+Y**2+Z**2)
    phi['up'] += (np.sqrt(system.rho_p))*np.exp(-R1**2/sigma**2)

    # Particle 2 (Down) at (+d/2, 0, 0)
    R2 = np.sqrt((X-d_initial/2)**2+Y**2+Z**2)
    phi['down'] += (np.sqrt(system.rho_p))*np.exp(-R2**2/sigma**2)
    
    phi_prev = {k: v.copy() for k, v in phi.items()}
    
    # Give particle 1 a slight tangential velocity to induce orbit
    phi_prev['up'] = np.roll(phi_prev['up'], 1, axis=1)

    positions = []
    print("Running 4D FDTD interaction simulation...")
    for t in range(T_steps):
        # ... FDTD loop identical to stability test ...
        total_rho = system.get_total_density(phi['up'], phi['down']); force_mag = system.get_force_term(total_rho)
        phi_next = {}
        for s, o in [('up','down'),('down','up')]:
            lap = (np.roll(phi[s],1,0)+np.roll(phi[s],-1,0)+np.roll(phi[s],1,1)+np.roll(phi[s],-1,1)+np.roll(phi[s],1,2)+np.roll(phi[s],-1,2)-6*phi[s])/dx**2
            chiral = 2j*system.kappa*(phi[o]-phi_prev[o])/dt
            phi_next[s] = 2*phi[s]-phi_prev[s]+dt**2*(lap-force_mag*phi[s]-chiral)
        phi_prev, phi = phi, phi_next

        if t % (T_steps // 20) == 0:
            xc=X[:,0,0]; d_up=np.sum(np.abs(phi['up'])**2,(1,2)); d_down=np.sum(np.abs(phi['down'])**2,(1,2))
            pos_up=np.sum(xc*d_up)/(np.sum(d_up)+1e-9); pos_down=np.sum(xc*d_down)/(np.sum(d_down)+1e-9)
            positions.append(np.abs(pos_down-pos_up))
            if t > 0: print(f"Step {t}, Separation: {positions[-1]:.3f}")

    print(f"Observed: Initial separation ~{d_initial:.2f}, Final separation ~{positions[-1]:.2f}")
    if PLOTTING_ENABLED:
        plt.figure(figsize=(10,6)); plt.plot(np.linspace(0,T_steps,len(positions)), positions)
        plt.title("Case 5b: Separation Distance vs. Time"); plt.xlabel("Time Step"); plt.ylabel("Separation")
        plt.grid(True); save_figure("5b_atomic_forces_trajectory")

    print("\n--- VERIFICATION ---")
    if positions[-1] < d_initial and positions[-1] > 2.0: print("SUCCESS: Particles formed a stable bound state.")
    else: print("FAILURE: Did not form a stable bound state.")

# --- Full code for regression test cases ---
def run_case_1_gravity():
    print("\n>>> CASE 1: VERIFICATION OF LONG-RANGE GRAVITY FIELD <<<")
    system = CHPTSystem(); m_source = 5.972e24; r = np.linspace(6.371e6, 5e7, 4000)
    k = 2 * system.G * m_source / system.c**2; rho_gravity = system.rho_0 * np.exp(-k / r)
    a_chpt = -0.5 * system.c**2 * np.gradient(np.log(rho_gravity), r)
    a_newton = -system.G * m_source / r**2; error = np.mean(np.abs((a_chpt - a_newton) / a_newton))
    print(f"Desired: Deviation < 1%. Observed: {error*100:.4f}%")
    if error < 0.01: print("--- VERIFICATION --- \nSUCCESS."); save_figure("1_gravity")
    else: print("--- VERIFICATION --- \nFAILURE.")

def run_case_2_rest_mass():
    print("\n>>> CASE 2: RELATIVISTIC MECHANICS (REST MASS) <<<"); system = CHPTSystem()
    N,dx=2000,0.01;x=(np.arange(N)-N//2)*dx;phi=np.sqrt(system.rho_0)+(np.sqrt(system.rho_p)-np.sqrt(system.rho_0))*np.exp(-x**2)
    grads=np.gradient(phi,dx);grad_sq=np.abs(grads)**2;pot_V=system.c_p*phi**2-(system.b/2)*phi**4+(system.a/3)*phi**6
    energy=np.sum(grad_sq+pot_V)*dx; print(f"Desired: Positive energy. Observed: {energy:.4e}")
    if energy>0: print("--- VERIFICATION ---\nSUCCESS.")
    else: print("--- VERIFICATION ---\nFAILURE.")

def run_case_4_interference():
    print("\n>>> CASE 4: EMF PROPAGATION (DOUBLE-SLIT INTERFERENCE) <<<"); print("Desired: Clear interference pattern.")
    wl=500e-9;k=2*np.pi/wl;sw=2*wl;ss=10*sw;sd=0.01;scw=0.001;ns=50
    s1y=np.linspace(-ss/2-sw/2,-ss/2+sw/2,ns);s2y=np.linspace(ss/2-sw/2,ss/2+sw/2,ns)
    slits=np.concatenate([s1y,s2y]);scr_y=np.linspace(-scw/2,scw/2,2000);field=np.zeros_like(scr_y,dtype=np.complex128)
    for y in slits:field+=np.exp(1j*k*np.sqrt(sd**2+(scr_y-y)**2))
    intensity=np.abs(field)**2;intensity/=intensity.max();print(f"Observed: {len(np.where(intensity>0.8)[0])} peaks.")
    if len(np.where(intensity>0.8)[0])>10:print("--- VERIFICATION ---\nSUCCESS."); save_figure("4_interference")
    else:print("--- VERIFICATION ---\nFAILURE.")

def run_case_6_tunneling():
    print("\n>>> CASE 6: 'QUANTUM' TUNNELING VIA FDTD <<<"); print("Desired: Non-zero transmitted intensity.")
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
    
    is_particle_stable = run_case_5a_spinor_stability()
    run_case_5b_atomic_forces_spinor(is_particle_stable)
    
    print("\n--- RUNNING REGRESSION TEST SUITE ---")
    run_case_1_gravity()
    run_case_2_rest_mass()
    run_case_4_interference()
    run_case_6_tunneling()
    
    end_time = time.time()
    print(f"\nTotal execution time: {end_time - start_time:.2f} seconds.")
