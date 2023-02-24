"""
Testing differential equation solutions

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import pandas as pd


# read data files
df_cure_depth = pd.read_csv("data/curedepth.csv")
cure_depth_measurement = np.asarray(df_cure_depth)

df_dark_cure = pd.read_csv("data/darkcure.csv")
dark_cure_measurement = np.asarray(df_dark_cure)

depth = cure_depth_measurement[:, 0]  * 1e-6
time_to_cure = cure_depth_measurement[:, 1]


# formulation
percent_PI = 0.0333                               #|   wt.%  | weight percent of photo initiator
percent_PEA = 0.15                                #|   wt.%  | weight percent of PEA
percent_HDDA = 0.0168                             #|   wt.%  | weight percent of HDDA
percent_8025D = 0.4084                            #|   wt.%  | weight percent of 8025D
percent_8025E = 0.0408                            #|   wt.%  | weight percent of 8025E
percent_E4396 = 0.3507                            #|   wt.%  | weight percent of HDDA
percent_M = percent_PEA + percent_HDDA            #|   wt.%  | weight percent of monomer

# physical parameters
rho_PEA = 1.02e3                                  #| kg/m^3  | density of PEA (estimated)
rho_HDDA = 1.02e3                                 #| kg/m^3  | density of HDDA (estimated)
rho_E4396 = 1.10e3                                #| kg/m^3  | density of EBECRYL 4396 
rho_M = 0.899 * rho_PEA + 0.101 * rho_HDDA        #| kg/m^3  | weighted average density of monomer
rho_P = 0.03 * rho_HDDA + \
        0.29 * rho_PEA + \
        0.68 * rho_E4396                          #| kg/m^3  | weighted average density of polymer
# rho_P = rho_E4396                                 #| kg/m^3  | weighted average density of polymer
rho_UGAP = 1840                                   #| kg/m^3  | estimated density of UGAP

mw_PEA = 0.19221                                  #|  kg/mol | molecular weight of PEA
mw_HDDA = 0.226                                   #|  kg/mol | molecular weight of HDDA
mw_M = 0.899 * mw_PEA + 0.101 * mw_HDDA           #|  kg/mol | weighted average molecular weight of monomer
mw_PI = 0.4185                                    #|  kg/mol | molecular weight of photo initiator

basis_wt = 0.05                                   #|   kg    | arbtirary starting ink weight
basis_vol = basis_wt / rho_UGAP                   #|   m^3   | arbitrary starting ink volume
mol_PI = basis_wt * percent_PI / mw_PI            #|   mol   | required PI for basis weight
mol_M = basis_wt * percent_M / mw_M               #|   mol   | required monomer for basis weight

c_PI0 = mol_PI / basis_vol                        #| mol/m^3 | initial concentration of PI
c_M0 = mol_M / basis_vol                          #| mol/m^3 | initial concentration of monomer
c_Mrad0 = 0                                       #| mol/m^3 | initial concentral monomer radicals

I_0 = 6.0                                         #|  W/m^2  | initial laser intensity
L = depth[-1]                                     #|    m    | dimension of RVE
Dm0 = 2.36e-6                                     #|  m^2/s  | diffusion constant pre-exponential, monomer
Am = 0.66                                         #| unitless| diffusion constant parameter, monomer

# SHANNA PARAMS
shanna_c_PI0 = 150                                #| mol/m^3 | initial PI concentration
shanna_c_M0 = 8.25e3                              #| mol/m^3 | initial monomer concentration 
shanna_mw_M = 130.14e-3                           #| mol/kg  | mw of monomer
shanna_mw_PI = 256.301e-3                         #| mol/kg  | mw of PI

# c_PI0 = shanna_mw_PI                   
# c_M0 = shanna_c_M0                          
# mw_M = shanna_mw_M
# mw_PI = shanna_mw_PI

# compute reaction rate constants
def rxn_rate_constant(i, c_M, theta):
    """
        inputs:
            c_M - concentration of monomer (mol / m^3)
            theta - temperature (K)
        outputs:
            k_p - polymerization rate constant 
            k_t - termination rate constant
        
    """
    # system parameters
    Rg = 8.3145                                   #| J/mol K | universal gas constant
    alpha_P = 0.000075                            #|   1/K   | coefficent of thermal expansion, polymerization (lit.)
    alpha_M = 0.0005                              #|   1/K   | coefficent of thermal expansion, monomer (lit.)
    theta_gP = 236.75                             #|    K    | glass transition temperature, polymer UGAP
    theta_gM = 328                                #|    K    | glass transition temperature, monomer (lit.)

    # constants polymerization 
    k_P0 = 1.6e3                                  #|m^3/mol s| true kinetic constant, polymerization (lit.)
    E_P = 18.23e3                                 #|  J/mol  | activation energy, polymerization (lit.)
    A_Dp = 0.5                                     #| unitless| diffusion parameter, polymerization (lit.)
    f_cp = 0.042                                  #| unitless| critical free volume, polymerization (lit.)
    
    # constants for termination
    k_T0 = 3.6e3                                  #|m^3/mol s| true kinetic constant, termination (lit.)
    E_T = 2.94e3                                  #|  J/mol  | activation energy, termination (lit.)
    A_Dt = 1.2                                    #| unitless| activation energy, termination (lit.)
    f_ct = 0.060                                  #| unitless| critical free volume, termination (lit.)
    R_rd = 4                                      #|  1/mol  | reaction difussion parameter (lit.)

    # compute the total fractional free volume
    vT = c_M * mw_M / rho_M + (c_M0 - c_M) * mw_M / rho_P
    phi_M = c_M * mw_M / rho_M / vT
    phi_P = (c_M0 - c_M) * mw_M / rho_P / vT
    f = 0.025 + alpha_M * phi_M * (theta - theta_gM) + alpha_P * phi_P * (theta - theta_gP)

    # # compute rate constants BOWMAN - TEMPERATURE
    const = (1 / f - 1 / f_cp)
    denom_p = (1 + np.exp(A_Dp * const))
    k_P = k_P0 * np.exp(-E_P / Rg / theta) / denom_p
    k_tr = R_rd * k_P * c_M
    k_T = k_T0 * np.exp(-E_T / Rg / theta) / (1 + 1 / (k_tr / (k_T0 * np.exp(-E_T / Rg / theta)) + np.exp(-A_Dt * (1 / f - 1 / f_ct))))

    # compute rate constants BOWMAN - NO TEMPERATURE
    # k_P = k_P0  / (1 + np.exp(A_Dp * (1 / f - 1 / f_cp)))
    # k_T = k_T0 / (1 + (R_rd * k_P * c_M / k_T0 + np.exp(-A_Dt * (1 / f - 1 / f_ct))))

    # k_T = k_T0 * np.ones(len(f))
    # k_P = k_P0 * np.ones(len(f))
    return k_P, k_T, f, denom_p
