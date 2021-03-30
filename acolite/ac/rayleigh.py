# function ray_tau
# computes Rayleigh optical thickness for wl in microns
# QV 2016-12-14

def ray_tau(wl, Patm=1013.25):
    tau_ray = Patm/1013.25*(0.008569*pow(wl,-4)*(1.+0.0113*pow(wl,-2)+0.00013*pow(wl,-4)))
    return tau_ray

# function ray_phase
# computes Rayleigh phase function for given geometry
# QV 2016-12-14

def ray_phase(theta_0,theta_v,phi_0, phi_v):
    from math import cos, sin, pow
    costheta_min = -1. * cos(theta_0)*cos(theta_v) - sin(theta_0)*sin(theta_v)*cos(abs(phi_0-phi_v))
    costheta_plus = 1. * cos(theta_0)*cos(theta_v) - sin(theta_0)*sin(theta_v)*cos(abs(phi_0-phi_v))
    phase_r= (0.75*(1.+pow(costheta_min,2.))) \
              + (sky_refl(theta_0)+sky_refl(theta_v)) * (0.75*(1.+pow(costheta_plus,2.)))
               
    return phase_r


# function ray_tr
# computes Rayleigh transmittance for given geometry
# QV 2016-12-14

def ray_tr(wl, theta_0, theta_v, Patm=1013.25):
    from math import cos, exp
    tau_ray = ray_tau(wl, Patm=Patm)
    ray_tr = (1.+exp(-1.*tau_ray/cos(theta_v))) * (1.+exp(-1.*tau_ray/cos(theta_0))) / 4.
    return ray_tr


# function ray_refl
# computes Rayleigh reflectance for given geometry
# QV 2016-12-14

def ray_refl(wl, theta_0, theta_v, phi_0, phi_v, Patm=1013.25, tau_ray=None):
    from math import cos
    if tau_ray is None: tau_ray = ray_tau(wl, Patm=Patm)
    phase_ray = ray_phase(theta_0,theta_v,phi_0, phi_v)
    rho_ray = (tau_ray * phase_ray) / (4. * cos(theta_0)*cos(theta_v))
    return rho_ray

# function sky_refl
# computes diffuse sky reflectance
# QV 2016-12-14

def sky_refl(theta, n_w=1.34):
    from numpy import arcsin, sin, tan, power
    # angle of transmittance theta_t for air incident rays (Mobley, 1994 p156)
    theta_t = arcsin(1./n_w*sin(theta))
    r_int=0.5*(power(sin(theta-theta_t)/sin(theta+theta_t),2)+\
              power(tan(theta-theta_t)/tan(theta+theta_t),2))
    return r_int


# function ray_phase_nosky
# computes Rayleigh phase function for given geometry (no diffuse sky reflectance)
# QV 2016-12-14

def ray_phase_nosky(theta_0,theta_v,phi_0, phi_v):
    from math import cos, sin, pow
    costheta_min = -1. * cos(theta_0)*cos(theta_v) - sin(theta_0)*sin(theta_v)*cos(abs(phi_0-phi_v))
    phase_r= (0.75*(1.+pow(costheta_min,2.)))
    return phase_r

# function ray_refl_nosky
# computes Rayleigh reflectance for given geometry (no diffuse sky reflectance)
# QV 2016-12-14

def ray_refl_nosky(wl, theta_0, theta_v, phi_0, phi_v, Patm=1013.25, tau_ray=None):
    from math import cos
    if tau_ray is None: tau_ray = ray_tau(wl, Patm=Patm)
    phase_ray = ray_phase_nosky(theta_0,theta_v,phi_0, phi_v)
    rho_ray = (tau_ray * phase_ray) / (4. * cos(theta_0)*cos(theta_v))
    return rho_ray


# function ray_phase_onlysky
# computes Rayleigh phase function for given geometry - only diffuse sky reflectance
# QV 2017-11-14

def ray_phase_onlysky(theta_0,theta_v,phi_0, phi_v):
    from math import cos, sin, pow
    costheta_plus = 1. * cos(theta_0)*cos(theta_v) - sin(theta_0)*sin(theta_v)*cos(abs(phi_0-phi_v))
    phase_r= (sky_refl(theta_0)+sky_refl(theta_v)) * (0.75*(1.+pow(costheta_plus,2.)))
    return phase_r

# function ray_refl_onlysky
# computes Rayleigh reflectance for given geometry - only diffuse sky reflectance
# QV 2017-11-14

def ray_refl_onlysky(wl, theta_0, theta_v, phi_0, phi_v, Patm=1013.25, tau_ray=None):
    from math import cos
    if tau_ray is None: tau_ray = ray_tau(wl, Patm=Patm)
    phase_ray = ray_phase_onlysky(theta_0,theta_v,phi_0, phi_v)
    rho_ray = (tau_ray * phase_ray) / (4. * cos(theta_0)*cos(theta_v))
    return rho_ray
