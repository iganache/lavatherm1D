import numpy as np
from scipy.special import erf



## Functions
def Stefan_lambda_upper_crust(lambda_0):
    LHS = np.exp(-1*lambda_0**2) / lambda_0 / erf(lambda_0)
    return LHS

def Stefan_lambda_basal_crust(lambda_0):
    LHS = np.exp(-1*lambda_0**2) / lambda_0 / (1+erf(lambda_0))
    return LHS

def crust_thickness(l, kappa, t):
    return 2 * l * (kappa * t)**0.5

def temp_basal_crust(Tlava, Tsub, l, kappa, thickness, time):
    m = (Tlava - Tsub) / (1+erf(l))
    n = 1 + erf(thickness / 2 / (kappa*time)**0.5)
    return ((m*n) + Tsub)

def temp_upper_crust(Tlava, Tsurf, l, kappa, thickness, time):
    m = (Tlava - Tsurf) / erf(l)
    n = erf(thickness / 2 / (kappa*time)**0.5)
    return (m*n) + Tsurf


def density(T, rho0 = 2600, beta=1.5e-5):
    rho = rho0 * (1 + beta/(T-1450))
    return rho

def heat_capacity(T):
    if T > 1010:
        return 1100
    elif T<= 1010:
        return 1211-(1.12e5/T)
    
def diffusivity(k, rho, c):
    return k/rho/c
    
def Radiation_flux(Ts, Ta, emis=0.95):
    SB_const = 5.67e-8
    Qrad = SB_const * emis * (Ts - Ta)**4    
    return Qrad

def Convection_flux(Ts, Ta, hf_coef=50):
    return hf_coef * (Ts - Ta)

def Snyder_flux(T):
    Qsnyder = 1.07 * 1e-13 * T**4.85 * 1000  
    return Qsnyder
    
    
def getTime(Ts, Tlava, Q, l, k, kappa):
    t_sqrt = k * (Tlava - Ts) / erf(l) / (np.pi*kappa)**0.5 / Q
    return t_sqrt**2
    
# def conductivity(T, phi_ves = 0.5):
#     ## ignoring k_radiation of the vesicles
#     krad = 0
    
#     kbas = 0.848 + 1.1e3 * T
#     kcond = kbas * ( 2*(1-phi_ves)*kbas + (1 + 2*phi_ves)*kgas ) / ( (2+phi_ves)*kbas + (1-phi_ves)*kgas )
    
    
    