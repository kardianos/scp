import numpy as np

from scipy.optimize import minimize

data = [
    {'A': 2, 'Z': 1, 'E': 2.224},
    {'A': 3, 'Z': 1, 'E': 8.482},
    {'A': 3, 'Z': 2, 'E': 7.718},
    {'A': 4, 'Z': 2, 'E': 28.3},
]

def get_r(A):

    return 1.2 * A**(1/3)

def get_S(A, Z):

    return 0.5 * ((2*Z - A)**2 - A)

def E_model (A, Z, lambda_, k):

    r = get_r(A)

    S = get_S(A, Z)

    return k * (-S) * np.exp(- r / lambda_) / r

def loss(params):

    k, lambda_ = params

    l = 0

    for d in data:

        e_pred = E_model(d['A'], d['Z'], lambda_, k)

        l += (e_pred - d['E'])**2

    return l

initial = [10.0, 1.0] # k, lambda_

result = minimize(loss, initial, bounds = (0, None), (0.01, None))

print(result);
