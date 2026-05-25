from fractions import Fraction

def compute_green_sum(x_frac, lmax):
    c0 = Fraction(1)
    c1 = Fraction(6) * x_frac
    w1 = Fraction(4, 21)
    
    if lmax == 0: return Fraction(0)
    if lmax == 1: return w1 * c1
    
    sum_val = w1 * c1
    c_prev2 = c0
    c_prev1 = c1
    
    for l in range(2, lmax + 1):
        ml = Fraction(l)
        cl = (Fraction(2) * (ml + 2) * x_frac * c_prev1 - (ml + 4) * c_prev2) / ml
        w = (ml + 3) / (Fraction(3) * ml * (ml + 6))
        sum_val += w * cl
        c_prev2 = c_prev1
        c_prev1 = cl
        
    return sum_val

lam_raw = compute_green_sum(Fraction(-1, 2), 100)
mu_raw = compute_green_sum(Fraction(127, 128), 100)

print(f"lam_raw: {lam_raw}")
print(f"mu_raw: {mu_raw}")
