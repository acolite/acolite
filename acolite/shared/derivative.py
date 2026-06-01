## def der and der2
## returns derivative and second derivative
## written by Quinten Vanhellemont, RBINS
## 2025-05-07
## modifications: 2026-06-01 (QV) added gradient based functions derg and derg2
##

def der(x, y):
    import numpy as np
    x = np.asarray(x)
    y = np.asarray(y)
    dy = np.diff(y,1)
    dx = np.diff(x,1)
    y1 = dy/dx
    x1 = 0.5*(x[:-1]+x[1:])
    return(x1, y1)

def der2(x, y):
    xd1, yd1 = der(x, y)
    xd2, yd2 = der(xd1, yd1)
    return(xd2, yd2)

def derg(x, y):
    import numpy as np
    res = np.gradient(y, x)
    return(res)

def derg2(x, y):
    res = derg(x, y)
    res2 = derg(x, res)
    return(res2)
