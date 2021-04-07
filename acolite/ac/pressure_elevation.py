## def pressure_elevation
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-07-17
## modifications: 2021-04-07 (QV) changed numpy import

def pressure_elevation(x, ratio=False, temperature=None, to_elevation=False):
    import numpy as np

    h = np.asarray(x) # elevation in m

    # from http://en.wikipedia.org/wiki/Atmospheric_pressure
    p0 = 101325 # standard pressure at sea level (Pa)
    L = 0.0065 # temperature lapse rate (K/m)
    Cp = 1007 # constant pressure specific heat (J/kg*K)
    T0 = 288.15 # standard temperature at sea level (K)
    g = 9.80665 # gravitational acceleration (m/s^2)
    M = 0.0289644 # molar mass of dry air (kg/mol)
    R = 8.31447 # universal gas constant (J/mol*K)

    T_ = T0 if temperature is None else float(temperature)

    h0 = (R*T_)/(M*g) # scale height in m

    if not to_elevation:
        p1 = p0 * np.exp(-1 * (h / h0))

        pressure = p1 / 100
        if ratio: return(pressure / 1013.25)
        else: return(pressure)
    else:
        p1 = float(h) * 100 ## pressure is given instead of elevation
        h = h0 * np.log(p0/p1)
        return(h)
