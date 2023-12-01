import numpy as np



## Functions
def Stefan_lambda_upper_crust(lambda_0):
    LHS = np.exp(-1*lambda_0**2) / lambda_0 / erf(lambda_0)
    return LHS

def Stefan_lambda_lower_crust(lambda_0):
    LHS = np.exp(-1*lambda_0**2) / lambda_0 / (1+erf(lambda_0))
    return LHS

def crustal_thickness(l, kappa, t):
    return 2 * l * (kappa * t)**0.5

def temp_lower_crust(Tlava, Tsub, l, kappa, thickness, time):
    m = (Tlava - Tsub) / (1+erf(l))
    n = 1 + erf(thickness / 2 / (kappa*time)**0.5)
    return ((m*n) + Tsub)

def temp_upper_crust(Tlava, Tsurf, l, kappa, thickness, time):
    m = (Tlava - Tsurf) / erf(l)
    n = erf(thickness / 2 / (kappa*time)**0.5)
    return (m*n) + Tsurf


def density(T, rho0, beta):
    rho = rho0 * (1 + beta/(T-1450))
    return rho

def heat_capacity(T):
    if T > 1010:
        return 1100
    elis T<= 1010:
        return 1211-(1.12e5/T)
    
# def conductivity(T, phi_ves = 0.5):
#     ## ignoring k_radiation of the vesicles
#     krad = 0
    
#     kbas = 0.848 + 1.1e3 * T
#     kcond = kbas * ((2*(1-phi_ves)*kbas) + ((1 + 2*phi_ves)*kgas)) / ()
    
    
    