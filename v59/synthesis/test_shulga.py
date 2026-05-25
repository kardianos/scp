def greenSumLoop(steps, c_prev2, c_prev1, sum_val, ml, x):
    if steps == 0:
        return sum_val
    cl = (2 * (ml + 2) * x * c_prev1 - (ml + 4) * c_prev2) / ml
    w = (ml + 3) / (3 * ml * (ml + 6))
    return greenSumLoop(steps - 1, c_prev1, cl, sum_val + w * cl, ml + 1, x)

def computeGreenSum(x, lmax):
    if lmax == 0: return 0
    if lmax == 1:
        c1 = 6 * x
        w1 = 4 / 21
        return w1 * c1
    c0 = 1
    c1 = 6 * x
    w1 = 4 / 21
    sum1 = w1 * c1
    return greenSumLoop(lmax - 1, c0, c1, sum1, 2, x)

print("lambda_raw (-1/2) =", computeGreenSum(-0.5, 100))
print("mu_raw (127/128) =", computeGreenSum(127/128, 100))
print("geometric_ratio =", computeGreenSum(-0.5, 100) / computeGreenSum(127/128, 100))
