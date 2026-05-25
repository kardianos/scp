import numpy as np
from scipy.special import eval_gegenbauer
import math

def volume_S7():
    return (math.pi**4) / 3.0

def degeneracy(l):
    # D(l, d) = (2l+d-1)/(l+d-1) * (l+d-1 choose d-1)
    # For S^7, d=7
    return (2 * l + 6) / 6.0 * math.comb(l + 5, 5)

def eigenvalue(l):
    # E(l) = l(l+d-1)
    return l * (l + 6)

def C3_l_1(l):
    # C_l^{(3)}(1) = (2*3)_l / l! = (l+5)! / (5! l!) = choose(l+5, 5)
    return math.comb(l + 5, 5)

def green_function(theta, l_max=1000):
    G = 0.0
    for l in range(1, l_max + 1):
        num = degeneracy(l) * eval_gegenbauer(l, 3, np.cos(theta))
        den = eigenvalue(l) * C3_l_1(l)
        G += num / den
    return G / volume_S7()

if __name__ == "__main__":
    # The family shift angle (Z3 triality)
    theta_shift = 2.0 * np.pi / 3.0
    
    # We evaluate the interaction strength lambda (which comes from the cross term G(theta_shift))
    # and the self energy mu (which comes from G(0)). 
    # G(0) is divergent in 7D, so we regularize it by evaluating at a small cutoff angle 
    # corresponding to the discrete "Planck" length of the 64-dimensional algebra.
    # The algebra has 64 dimensions, so a geometric cutoff might be theta ~ 1/sqrt(64) = 1/8
    theta_cutoff = 1.0 / 8.0
    
    l_max_cutoff = 100 # Smooth cutoff for the infinite sum to avoid Gegenbauer blowups at high l
    
    # Lambda is the overlap at the shift
    lam_raw = green_function(theta_shift, l_max=l_max_cutoff)
    
    # Mu is the regularized self-energy
    mu_raw = green_function(theta_cutoff, l_max=l_max_cutoff)
    
    print(f"G(theta=2pi/3) [lambda proxy]: {lam_raw}")
    print(f"G(theta=cutoff) [mu proxy]: {mu_raw}")
    print(f"Ratio lambda / |mu|: {abs(lam_raw / mu_raw)}")
    print(f"Empirical Target Ratio (0.012 / 41.345): {0.012 / 41.345}")
