## Predict annular cumulative distribution function
## QV 2021-07-06 from AC v0.9 R function
#### 2021-07-06 (QV) default pressure = 1 (normalised to 1013.25)

def pred_annular_cdf(r, cf, pressure=1):
    import numpy as np

    c1 = cf[0] * 1.0
    c2 = cf[1] * pressure ** -1
    c3 = cf[2] * pressure ** -1
    c4 = (c1 - c2) * cf[3]
    c5 = cf[4] * 1.0
    c6 = (c1 - c2) * (1 - cf[3])
    c7 = cf[5]
    res = c1 - (c2 * np.exp(c3 * r) + c4 * np.exp(c5 * r) + c6 * np.exp(c7 * r))

    if type(res) is np.ndarray:
        res[np.isinf(res)] = 0
    else:
        if not np.isfinite(res): res = 0
    return(res)
