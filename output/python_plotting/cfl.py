
#%%
import numpy as np

thermal_cond = 0.2
heat_cap     = 800
rho          = 1000
kappa        = thermal_cond / heat_cap / rho
D            = 1.08e-6
nodes = [11, 15, 21, 25, 31, 41, 51]
h     = [0.00084 / (node - 1) for node in nodes]
dt    = [.5 * h_i ** 2 / D for h_i in h]

print('dt: ', dt)


# %%
