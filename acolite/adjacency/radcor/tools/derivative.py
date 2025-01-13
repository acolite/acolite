## def derivative
## computes 1st order derivative with 3 points
##
## QV 2025-01-13 based on deprecated scipy code
##
## modifications:

def derivative(func, x0, dx=1.0, args=(), ):
    import numpy as np

    n = 1
    order = 3
    weights = np.array([-1, 0, 1]) / 2.0

    val = 0.0
    ho = order >> 1
    for k in range(order):
        val += weights[k] * func(x0 + (k - ho) * dx, *args)
    return(val / np.prod((dx,) * n, axis=0))
