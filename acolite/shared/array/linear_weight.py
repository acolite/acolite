## def linear_weight
## find indices and weights for linear scaling of an array and given value
## written by Quinten Vanhellemont, RBINS
## 2026-06-01
## modifications:

def linear_weight(x, x_val):
    import numpy as np
    idx = np.interp(x_val, x, np.arange(len(x)))
    idx_0 = int(np.floor(idx))
    w = 1 - (idx - idx_0)
    return(idx_0, w, idx_0+1, (1-w))
