"""
Testing differential equation solutions


next steps:
    - add heat equation
    - use GA for further parameter tuning
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import pandas as pd
import copy

# read data files
# df_cure_depth = pd.read_csv("data/curedepth.csv")
# cure_depth_measurement = np.asarray(df_cure_depth)

# df_dark_cure = pd.read_csv("data/darkcure.csv")
# dark_cure_measurement = np.asarray(df_dark_cure)

# depth = cure_depth_measurement[:, 0]  * 1e-6
# time_to_cure = cure_depth_measurement[:, 1]


# formulation - wt. %
percent_PI = 0.0333                               #|   wt.%  | weight percent of photo initiator
percent_PEA = 0.15                                #|   wt.%  | weight percent of PEA
percent_HDDA = 0.0168                             #|   wt.%  | weight percent of HDDA
percent_8025D = 0.4084                            #|   wt.%  | weight percent of 8025D
percent_8025E = 0.0408                            #|   wt.%  | weight percent of 8025E
percent_E4396 = 0.3507                            #|   wt.%  | weight percent of HDDA

# # # NEW formulation - wt. %
# percent_PI = 0.01                               #|   wt.%  | weight percent of photo initiator
# percent_PEA = 0.15                                #|   wt.%  | weight percent of PEA
# percent_HDDA = 0.020                             #|   wt.%  | weight percent of HDDA
# percent_8025D = 0.420                            #|   wt.%  | weight percent of 8025D
# percent_8025E = 0.04                            #|   wt.%  | weight percent of 8025E
# percent_E4396 = 0.36                            #|   wt.%  | weight percent of HDDA

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
# mw_M = mw_PEA

basis_wt = 0.5                                   #|   kg    | arbtirary starting ink weight
basis_vol = basis_wt / rho_UGAP                   #|   m^3   | arbitrary starting ink volume
mol_PI = basis_wt * percent_PI / mw_PI            #|   mol   | required PI for basis weight
mol_M = basis_wt * percent_M / mw_M               #|   mol   | required monomer for basis weight

# inital concentrations
c_PI0 = mol_PI / basis_vol                        #| mol/m^3 | initial concentration of PI
c_M0 = mol_M / basis_vol                          #| mol/m^3 | initial concentration of monomer
c_Mrad0 = 0                                       #| mol/m^3 | initial concentral monomer radicals

I_0 = 6.0                                        #|  W/m^2  | initial laser intensity
L = depth[-1]                                     #|    m    | dimension of RVE
Dm0 = 1.08e-6                                     #|  m^2/s  | diffusion constant pre-exponential, monomer
Am = 0.66                                         #| unitless| diffusion constant parameter, monomer
k_i = 4.8e-5                                      #|  s^-1   | primary radical rate constant
dHp = 54.8e3                                      #|  W/mol  | heat of polymerization of acrylate monomers
Cp_nacl = 880.0                                   #| J/kg/K  | heat capacity of NaCl
Cp_pea = 180.3                                    #| J/mol/K | heat capacity of PEA @ 298K - https://polymerdatabase.com/polymer%20physics/Cp%20Table.html
Cp_pea /= mw_PEA                                  #| J/kg/K  | convert units
Cp_hdda = 202.9                                   #| J/mol/K | solid heat capacity of HDDA - https://webbook.nist.gov/cgi/cbook.cgi?ID=C629118&Units=SI&Mask=1F
Cp_hdda /= mw_HDDA                                #| J/kg/K  | convert units

# SHANNA PARAMS
shanna_c_PI0 = 150                                #| mol/m^3 | initial PI concentration
shanna_c_M0 = 8.25e3                              #| mol/m^3 | initial monomer concentration 
shanna_mw_M = 130.14e-3                           #| mol/kg  | mw of monomer
shanna_mw_PI = 256.301e-3                         #| mol/kg  | mw of PI
Cp_shanna = 1700                                  #| J/kg/K  | shanna's heat capacity
K_thermal_shanna = 0.2                            #| W/m/K   | shanna's thermal conductivity

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
    theta_gM = 313.6                              #|    K    | glass transition temperature, monomer (Taki lit.)

    # constants polymerization 
    k_P0 = 1.145e2                                #|m^3/mol s| true kinetic constant, polymerization (lit.)
    E_P = 10.23e3                                 #|  J/mol  | activation energy, polymerization (lit.)
    A_Dp = 0.05                                   #| unitless| diffusion parameter, polymerization (lit.)
    f_cp = 5.17e-2                                #| unitless| critical free volume, polymerization (lit.)
    
    # constants for termination
    k_T0 = 1.337e3                                #|m^3/mol s| true kinetic constant, termination (lit.)
    # k_T0 = 1.5e2                                #|m^3/mol s| true kinetic constant, termination (lit.)

    E_T = 2.94e3                                  #|  J/mol  | activation energy, termination (lit.)
    A_Dt = 1.2                                    #| unitless| activation energy, termination (lit.)
    f_ct = 5.81e-2                                #| unitless| critical free volume, termination (lit.)
    R_rd = 11                                     #|  1/mol  | reaction difussion parameter (lit.)

    # compute the total fractional free volume
    vT = c_M * mw_M / rho_M + (c_M0 - c_M) * mw_M / rho_P
    phi_M = c_M * mw_M / rho_M / vT
    phi_P = (c_M0 - c_M) * mw_M / rho_P / vT
    f = 0.025 + alpha_M * phi_M * (theta - theta_gM) + alpha_P * phi_P * (theta - theta_gP)

    # # compute rate constants - TEMPERATURE DEPENDENT

    # BOWMAN 1  eqn 17
    const1 = (1 / f - 1 / f_cp)
    denom_p = (1 + np.exp(A_Dp * const1))
    k_P = k_P0 * np.exp(-E_P / Rg / theta) / denom_p

    # BOWMAN 1  eqn 18
    k_tr = R_rd * k_P * c_M
    const2 = -A_Dt * (1 / f - 1 / f_ct)
    denom_t = (k_tr / (k_T0 * np.exp(-E_T / Rg / theta)) + np.exp(const2))
    k_T = k_T0 * np.exp(-E_T / Rg / theta) / (1 + 1 / denom_t)

    # compute rate constants BOWMAN - NO TEMPERATURE
    # these equations were used in Taki et al
    # k_P = k_P0  / (1 + np.exp(A_Dp * (1 / f - 1 / f_cp)))
    # k_T = k_T0 / (1 + (R_rd * k_P * c_M / k_T0 + np.exp(-A_Dt * (1 / f - 1 / f_ct))))

    # k_T = k_T0 * np.ones(len(f))
    # k_P = k_P0 * np.ones(len(f))
    return k_P, k_T, f, f

eps = 9.66e-1                                 #|m^3/mol m| initiator absorbtivity
phi = 1.0                                     #| unitless| quantum yield inititation
dt = 1e-3                                     #|    s    | time step
T = 15.0                                      #|    s    | total simulation time

# discretize time domain
t_step = np.arange(0, T + dt, dt)

# discretize spacial domain
node = 50
zspace = np.linspace(0, L, node)

# allocate memory for 
z_absorb = np.zeros(len(zspace))
dz = zspace[1]

# INITIAL CONDITIONS
c_PI_space = np.ones(len(zspace)) * c_PI0
c_PIdot_space = np.zeros(len(zspace))
c_Mdot_space = np.zeros(len(zspace))
c_M_space = np.ones(len(zspace)) * c_M0
# c_P_space = np.ones(len(zspace)) 
theta_space = np.ones(len(zspace)) * 303.15

# A matrix
A1 = np.zeros((node, node))
for i in range(node):
    if i == 0:
        A1[i, i] = 0
    elif i == node - 1:
        A1[i, i] = 0
    else:
        A1[i, i] = -2
        A1[i, i - 1] = 1
        A1[i, i + 1] = 1


# solution arrays
time = []
sol_intensity_profile = []
sol_c_PI = []

sol_c_PIdot = []
sol_c_Mdot = []
sol_c_M = []
sol_kt = []
sol_kp = []
sol_test = []
sol_consume = []
sol_diffuse = []
sol_theta = []

tol = 5e-4
thresh = 50
for i, t in enumerate(t_step):

    # compute energy intensity and reaction rate constants
    k_P, k_T, f, _ = rxn_rate_constant(i, c_M_space, theta_space)
    

    z_intensity = np.zeros(len(zspace))
    z_intensity[0] = I_0   
    for j in range(len(z_intensity) - 1):
        z_intensity[j + 1] = z_intensity[j] - dz * eps * c_PI_space[j] * z_intensity[j]
    
    # first compute reaction rate of photo initiator
    rhs_1 = lambda z_intensity, c_PI: -(phi / 2 * eps * z_intensity * c_PI)
    # c_PI_space_ = c_PI_space + dt * rhs_1(z_intensity, c_PI_space)
    err_1 = 100
    iter_1 = 0
    sol1_0 = copy.deepcopy(c_PI_space)
    while err_1 > tol:
        iter_1 +=1 
        sol1_1 = c_PI_space + dt * rhs_1(z_intensity, sol1_0)
        err_1 = np.linalg.norm(sol1_0 - sol1_1, ord=2)
        sol1_0 = copy.deepcopy(sol1_1)
        
        if iter_1 > thresh:
            raise Exception("--- PDE 1 did not converge | err: {} ---".format(err_1))

    c_PI_space_ = copy.deepcopy(sol1_0)

    # second compute reaction rate of PIdot
    rhs_2 =  lambda z_intensity, c_PI, c_PIdot, c_M: phi * eps * z_intensity * c_PI - k_i * c_PIdot * c_M
    
    err_2 = 100 
    iter_2 = 0
    sol2_0 = copy.deepcopy(c_PIdot_space)
    while err_2 > tol:
        iter_2 += 1
        sol2_1 = c_PIdot_space + dt * rhs_2(z_intensity, c_PI_space, c_PIdot_space, c_M_space)
        err_2 = np.linalg.norm(sol2_0 - sol2_1, ord=2)
        sol2_0 = copy.deepcopy(sol2_1)
        
        if iter_2 > thresh:
            raise Exception("--- PDE 2 did not converge | err: {} ---".format(err_2))

    c_PIdot_space_ = copy.deepcopy(sol2_0)

    # third compute reaction rate of monomer radical
    rhs_3 = lambda z_intensity, c_PIdot, c_M, c_Mdot: (k_i * c_PIdot * c_M - k_T * c_Mdot ** 2 + Dm0 * np.exp(-Am / f) * A1 @ c_Mdot / dz ** 2)
    # c_Mdot_space_ = c_Mdot_space + dt * rhs_3(z_intensity, c_PI_space, c_Mdot_space)

    err_3 = 100
    iter_3 = 0
    sol3_0 = copy.deepcopy(c_Mdot_space)
    while err_3 > tol:
        iter_3 += 1
        sol3_1 = c_Mdot_space + dt * rhs_3(z_intensity, c_PIdot_space, c_M_space, sol3_0)
        err_3 = np.linalg.norm(sol3_0 - sol3_1, ord=2)
        sol3_0 = copy.deepcopy(sol3_1)

        if iter_3 > thresh:
            raise Exception("--- PDE 3 did not converge | err: {} ---".format(err_3))
    c_Mdot_space_ = copy.deepcopy(sol3_0) #c_Mdot_space + dt * rhs_3(z_intensity, sol3_0, c_Mdot_space)
    test = Dm0 * np.exp(-Am / f) * A1 @ c_Mdot_space / dz ** 2
    # fourth compute rate of consumption of monomer

    err_4 = 100
    iter_4 = 0 
    sol4_0 = copy.deepcopy(c_M_space)
    while err_4 > tol:
        iter_4 +=1
        Dm = Dm0 * np.exp(-Am / f)
        diffuse = Dm / dz**2 * A1 @ sol4_0
        consume = k_P * sol4_0 * c_Mdot_space + k_i * c_PIdot_space * sol4_0
        sol4_1 = c_M_space + (diffuse - consume) * dt

        err_4 = np.linalg.norm(sol4_0 - sol4_1, ord=2)

        sol4_0 = copy.deepcopy(sol4_1)
        
        if iter_4 > thresh:
            raise Exception("--- PDE 4 did not converge | err: {} ---".format(err_4))

    c_M_space_ = copy.deepcopy(sol4_0) #c_M_space + dt * (Dm / dz**2 * A1 @ sol4_0 - k_P * sol4_0 * c_Mdot_space)


    # enforce neumann bc's
    c_M_space_[-1] = c_M_space_[-2]
    c_M_space_[0] = c_M_space_[1]

    # # fifth compute rate of temperature change
    err_5 = 100
    iter_5 = 0
    sol5_0 = copy.deepcopy(theta_space)
    
    term2 = k_P * c_M_space * c_Mdot_space * dHp
    term3 = eps * c_PI_space * z_intensity
    while err_5 > tol:
        iter_5 += 1
        
        term1 = K_thermal_shanna / dz ** 2 * A1 @ sol5_0 
        sol5_1 = theta_space + (term1 + term2 + term3) * dt / rho_UGAP / Cp_shanna
        
        err_5 = np.linalg.norm(sol5_0 - sol5_1, ord=2)
        
        sol5_0 = copy.deepcopy(sol5_1)
        
        if iter_5 > thresh: 
            raise Exception("--- PDE 5 did not converge | err: {} ---".format(err_5))
        
    theta_space_ = copy.deepcopy(sol5_0)

    # enforce neumann bc's
    theta_space_[-1] = theta_space_[-2]
    theta_space_[0] = theta_space_[1]

    if i % 1000 == 0:
        time.append(t)
        sol_intensity_profile.append(z_intensity)
        sol_c_PI.append(c_PI_space_)
        sol_c_PIdot.append(c_PIdot_space_)
        sol_c_Mdot.append(c_Mdot_space_)
        sol_c_M.append(c_M_space_)
        sol_theta.append(theta_space_)
        sol_kt.append(k_T)
        sol_kp.append(k_P)
        sol_consume.append(consume)
        sol_diffuse.append(diffuse)
        sol_test.append(test)

    c_PI_space = c_PI_space_
    c_PIdot_space = c_PIdot_space_
    c_Mdot_space = c_Mdot_space_
    c_M_space = c_M_space_
    theta_space = theta_space_


    
time = np.asarray(time)
sol_intensity_profile = np.asarray(sol_intensity_profile)
sol_c_PI = np.asarray(sol_c_PI)
sol_c_PIdot = np.asarray(sol_c_PIdot)
sol_c_Mdot = np.asarray(sol_c_Mdot)
sol_c_M = np.asarray(sol_c_M)
sol_kt = np.asarray(sol_kt)
sol_kp = np.asarray(sol_kp)
sol_diffuse = np.asarray(sol_diffuse)
sol_consume = np.asarray(sol_consume)
sol_test = np.asarray(sol_test)
print("np.amin(sol_c_M): ", np.amin(sol_c_M))

# plt.figure(figsize=(15, 10))
# for i in range(15):

#     plt.plot(zspace, sol_c_PI[i], label="time: {}".format(np.round(time[i], 2)), lw=5)

# plt.title('PI concentration profile over time', fontsize=fs)
# plt.xlabel(r'z ($m$)', fontsize=fs)
# plt.ylabel(r'concentration of PI ($mol/m^3$)', fontsize=fs)
# plt.xticks(fontsize=fs-5)
# plt.yticks(fontsize=fs-5)
# plt.xlim(0, 0.0012)
# plt.legend(loc='best', fontsize=fs-5)
# plt.show()


print('\n-----------------------------------------\n')
# plot time evolution of PI concentration at the top/bottom surfaces:

fs_ = 30
plt.figure(figsize=(10, 7.5))
plt.scatter(time, sol_c_PI[:, 1], color='r', s=150, label='top layer')
plt.scatter(time, sol_c_PI[:, -2], color='b', s=150, label='bottom layer')
plt.legend(fontsize=fs_-5)
plt.title(r'concentration of PI ($mol/m^3$)', fontsize=fs_)
plt.ylabel(r'concentration ($mol/m^3$)', fontsize=fs_-5)
plt.xlabel('time ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)
plt.show()

# plot time evolution free radical concentration
plt.figure(figsize=(10, 7.5))
plt.scatter(time, sol_c_PIdot[:, 1], color='r', s=150, label='top layer')
plt.scatter(time, sol_c_PIdot[:, -2], color='b', s=150, label='bottom layer')
plt.legend(fontsize=fs_-5)
plt.title(r'concentration of PI$\cdot$ ($mol/m^3$)', fontsize=fs_)
plt.ylabel(r'concentration ($mol/m^3$)', fontsize=fs_-5)
plt.xlabel('time ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)
plt.show()

# plot time evolution of monomer dot concentration at the top/bottom surfaces:
plt.figure(figsize=(10, 7.5))
plt.scatter(time, sol_c_Mdot[:, 0], color='r', s=150, label='top layer')
plt.scatter(time, sol_c_Mdot[:, -1], color='b', s=150, label='bottom layer')
plt.legend(fontsize=fs_-5)
plt.title(r'concentration of M$\cdot$ ($mol/m^3$)', fontsize=fs_)
plt.ylabel(r'concentration ($mol/m^3$)', fontsize=fs_-5)
plt.xlabel('time ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)
plt.show()

# plot time evolution of M concentration at the top/bottom surfaces:
plt.figure(figsize=(10, 7.5))
plt.scatter(time, sol_c_M[:, 1], color='r', s=150, label='top layer')
plt.scatter(time, sol_c_M[:, -2], color='b', s=150, label='bottom layer')
plt.legend(fontsize=fs_-5)
plt.title(r'concentration of M ($mol/m^3$) ', fontsize=fs_)
plt.ylabel(r'concentration monomer ($mol/m^3$)', fontsize=fs_-5)
plt.xlabel('time ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)
plt.show()

# diffusion of M over time
plt.figure(figsize=(10, 7.5))
plt.scatter(time, sol_diffuse[:, 1] * 1e6, color='r', s=150, label='top layer')
plt.scatter(time, sol_diffuse[:, -2] * 1e6, color='b', s=150, label='bottom layer')
plt.legend(fontsize=fs_-5)
plt.title(r'diffusion of monomer over time', fontsize=fs_)
plt.ylabel(r'diffusion ($\mu mol/s^2$)', fontsize=fs_-5)
plt.xlabel('time ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)
plt.show()

# consumption of M over time
plt.figure(figsize=(10, 7.5))
plt.scatter(time[:], sol_consume[:, 1], color='r', s=150, label='top layer')
plt.scatter(time[:], sol_consume[:, -2], color='b', s=150, label='bottom layer')
plt.legend(fontsize=fs_-5)
plt.title(r'rate of consumption of monomer', fontsize=fs_)
plt.ylabel(r'rate of consumption ($mol/s$)', fontsize=fs_-5)
plt.xlabel('time ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)
plt.show()

plt.figure()
plt.scatter(time, sol_test[:, 0], color='r', label='$k_T$')
# plt.scatter(time, sol_test, color='r', label='$k_T$')
plt.legend()
plt.title(r'sol_test over time')
plt.ylabel(r'sol_test')
plt.xlabel('time ($s$)')
plt.show()

# plot time evolution of rate constant:
plt.figure(figsize=(10, 7.5))
plt.scatter(time, sol_kp[:, 1], color='b', s=150, label='$k_P$')
plt.legend(fontsize=fs_-5)
plt.title(r'$k_P$ over time', fontsize=fs_)
plt.ylabel(r'rate constant', fontsize=fs_-5)
plt.xlabel('time ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)
plt.show()

plt.figure(figsize=(10, 7.5))
plt.scatter(time, sol_kt[:, 1], color='r', s=150, label='$k_T$')
plt.legend(fontsize=fs_-5)
plt.title(r'$k_T$ over time', fontsize=fs_)
plt.ylabel(r'rate constant', fontsize=fs_-5)
plt.xlabel('time ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)
plt.show()

print('\n-----------------------------------------\n')



# %%
