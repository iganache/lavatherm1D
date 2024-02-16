import numpy as np

def viscosity(Terupt, Tcurrent, phi_crystal, phi_crystal_max, eta0):
    ## from Harris and Rowland (2001)
    a = 0.04  #in K-1
    eta = (1-(phi_crystal / phi_crystal_max))**-2.5  * eta0 * np.exp(a * (Terupt - Tcurrent))
    return eta

def yield_strength(Terupt, Tcurrent, phi_crystal):
    ## from Harris and Rowland (2001)
    a = 1e-2 ## in PA
    b = 0.08 ## in K-1
    ys = (6500*phi_crystal**2.85) + a*np.exp(b*(Terupt-Tcurrent)-1)
    return ys

def basal_shear_stress(rho, g, h, theta_deg):
    ## from Fagents, Gregg, and Lopes
    return rho * g * h * np.sin(np.rad2deg(theta_deg))

def velocity_Newtonian(rho, g, h, eta, theta_deg, C = 3):
    ## Jeffreys equation
    ## C = 3 for a filled channel
    ## C=8 for a wide flow with a moving top
    ## C = 12 for a sheet flow
    vel = rho*g*h**2*np.sin(np.deg2rad(theta_deg)) / C / eta
    return vel

def velocity_Bingham(rho, g, h, eta, theta_deg, yield_strength, c1=3, c2=3/2, c3=1/2, c4=3):
    ## from Fagents, Gregg, and Lopes
    basal_stress = basal_shear_stress(rho, g, h, theta_deg)
    tau = yield_strength/basal_stress
    vel = rho*g*h**2*np.sin(np.deg2rad(theta_deg)) / (c1*eta) * (1 - c2*tau + c3*tau**c4)
    return vel

def levee_width(rho, g, eta, theta_deg, yield_strength, effusion_rate):
    ## from Davies (1996)
    wb = yield_strength / (2 * g * rho * np.deg2rad(theta_deg)**2)
    return wb
    
def levee_thickness(rho, g, eta, theta_deg, yield_strength, effusion_rate):
    ## from Head and Wilson (1986)
    h = yield_strength / g / rho/ np.deg2rad(theta_deg)
    return h

def channel_width(rho, g, eta, theta_deg, yield_strength, effusion_rate, wb):
    ## from Davies (1996)
    E= effusion_rate
    wc1 = (24 * E * eta / yield_strength / np.deg2rad(theta_deg)**2) ** (1/3)
    wc2 = ((24 * E * eta)**4 * g * rho / yield_strength**5 / np.deg2rad(theta_deg)**6) ** (1/11)
   
    if wc1 <= 2*wb:
        return wc1
    elif wc2 > 2*wb:
        return wc2
    else:
        print("channel width computation error")
        

def channel_thickness(rho, g, eta, theta_deg, yield_strength, effusion_rate):
    wb = levee_width(rho, g, eta, theta_deg, yield_strength, effusion_rate)
    wc = channel_width(rho, g, eta, theta_deg, yield_strength, effusion_rate, wb)
    wavg = 0.5 *(2* wb + wc)
    h = (yield_strength * wavg / g / rho) ** 0.5
    return h


def cooling_timescale(rho, h, Tstart, Tend, C, flux):
    t_cool = (rho*C*h*(Tstart-Tend)) / flux
    return t_cool

def advancement_timescale(w, h, l, E):
    E= effusion_rate
    volume = l*w*h
    t_adv = volume / E
    return t_adv  

def flow_length():
    pass


def radiative_flux(Tsurf, Tatm, emis_surf = 0.95):
    sigmaB = 5.670374419e-8                    ## Stefan Boltzmann constant
    q = sigmaB * emis_surf * (Tsurf**4 - Tatm**4)
    return q

def convection_flux(Tsurf, Tatm, hconv = 50):
    q = hconv * (Tsurf - Tatm)
    return q

def conductive_flux():
    pass

def flynn_flux(Tsurface):
    ## from Flynn et al. (2023)
    q = 1.07e-13 * Tsurface**4.84 * 1e3
    return q

def coupled_flux(Tsurf, Tatm, k, emis_surf = 0.95):
    ## using formula from Siegel and Howell chapter 7
    ## 2 solutions are presented; the additive solution (eq 7-20) is used
    
    g = 8.87
    sigmaB = 5.670374419e-8                    ## Stefan Boltzmann constant
    
    k_air = 5.5e-2                             ## thermal conductivity of atmosphere Wm-1K-1 (Snyder, 2002)
    C_air = 1.15e3                             ## specific heat capacity of atmosphere Jkg-1K-1 (Snyder, 2002)
    rho_air = 65                               ## density of the atmosphere in kgm-3 (Snyder, 2002)
    eta_air = 3.4e-5                           ## dynamic viscosity of the atmosphere in kgm-3 (Snyder, 2002)
    beta_air = 1.3e-3                          ## thermal expansion coefficient (Snyder, 2002)
    kappa_air =  k_air / rho_air / C_air       ## thermal diffusivity of the atmosphere
    emis_air = 1                               ## emissivity of the atmosphere (Snyder, 2002)
    # emis_surf = 1                              ## emissivity of the surface (Snyder, 2002)
    
    H = 1                                   ## from Kesthelyi and Denlinger, 1995
    a = 0.1                                    ## gray absorption coefficient in m-1 (check if right)
    
    N = k*a / 4 / sigmaB / Tsurf**3             ## conduction parameter
    
    # Ra_crit = 2040 * 7100 / 2.3e-4           ## Snyder (2002) use N = 2.3e-4
    Ra_crit = 2040 * 7100 / N
    Ra = g * beta_air * (Tsurf - Tatm) * H**3 / eta_air / kappa_air
    D = (Ra_crit / Ra)**(1/3)
    
    qcond = a * (Tsurf-Tatm) / D
    qrad = sigmaB*(Tsurf**4-Tatm**4) / (0.74*D*a + emis_air**-1 + emis_surf**-1 - 1)
    qtotal = qrad + qcond
    return qtotal