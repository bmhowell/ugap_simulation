#%%

import glob
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
import matplotlib.cm as cm
import os
from matplotlib.ticker import ScalarFormatter

skip = 10

# central comuting 
# dir_path = '/home/brian/Documents/berkeley/ugap_simulation/output/python_plotting/'

# macbook computing
dir_path = '/Users/brianhowell/Desktop/Berkeley/MSOL/ugap_simulation/output/python_plotting/convergence'

time            = []
avg_tot_cPI     = []
avg_tot_cPIdot  = []
avg_tot_cMdot   = []
avg_tot_cM      = []
avg_tot_theta   = []
avg_free_volume = []
avg_k_t         = []
avg_k_p         = []
tot_df          = []    


counter = 1; 
csv_files = [filename for filename in os.listdir(dir_path) if filename.endswith('.csv')]
csv_files_sorted = sorted(csv_files)
for filename in csv_files_sorted:
    # import data with pandas
    file_path = os.path.join(dir_path, filename)
    print(file_path)
    df = pd.read_csv(file_path, sep=',')
    tot_df.append(df)
    
    time.append(df['time'].tolist()[::skip])

    avg_tot_cPI.append(    df[' avg_tot_cPI'].tolist()[::skip])
    avg_tot_cPIdot.append( df[' avg_tot_cPIdot'].tolist()[::skip])
    avg_tot_cMdot.append(  df[' avg_tot_cMdot'].tolist()[::skip])
    avg_tot_cM.append(     df[' avg_tot_cM'].tolist()[::skip])
    avg_tot_theta.append(  df[' avg_tot_theta'].tolist()[::skip])
    avg_free_volume.append(df[' avg_free_volume'].tolist()[::skip])
    avg_k_t.append(        df[' avg_k_t'].tolist()[::skip])
    avg_k_p.append(        df[' avg_k_p'].tolist()[::skip])

    counter += 1
# %%
# loop through all data and plot conversion
fs_   = 25
bump  = 1
exact = -1

DT    = [5e-3, 2.5e-3, 1e-3, 5e-4, 1e-4, 5e-5] 
NODES = [11,   15,     21,   31,   41,   51]
DX    = [0.00084 / (node - 1) for node in NODES]


plt.figure(figsize=(10, 7.5)) 
# plt.title(r'convergence: conversion of M', fontsize=fs_)
plt.title(r'convergence: cM', fontsize=fs_)
plt.ylabel(r'$L_2$ norm error', fontsize=fs_-5)
plt.xlabel('mesh size $dx$ ($m$)', fontsize=fs_-5)


plot_err = []
plot_dx  = []

exact_cM    = np.array(avg_tot_cM[exact])
for i in range(len(time)-bump):
    avg_tot_cM_plot = avg_tot_cM[i]

    err  = np.linalg.norm(avg_tot_cM_plot[:60] - exact_cM[:60])
    print('error: ', err)
    plot_dx.append(DX[i])
    plot_err.append(err)

plt.scatter(np.array(plot_dx), plot_err, s=150, label='nodes 21')
   
plt.legend(fontsize=10)
plt.yscale('log')
plt.xscale('log')
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)
# plt.savefig("convergence/converge_conversion.png")
plt.show()

# %%

plot_err = []
plot_dx  = []

# loop through all data and plot theta
plt.figure(figsize=(10, 7.5)) 
plt.title(r'error: average temperature', fontsize=fs_)
plt.ylabel(r'$L_2$ norm error', fontsize=fs_-5)
plt.xlabel('time step $k$ ($s$)', fontsize=fs_-5)
exact_temp = np.array(avg_tot_theta[exact])
for i in range(len(time) - bump):

    color = 'g'
    label = 'nodes = 21'
    avg_tot_theta_plot = avg_tot_theta[i]

    
    err  = np.linalg.norm(avg_tot_theta_plot[:60] - exact_temp[:60])
    print('error: ', err)
    plot_dx.append(DX[i])
    plot_err.append(err)
    


plt.scatter(plot_dx, plot_err, s=150, label='nodes 21')


plt.yscale('log')
plt.xscale('log')
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)

plt.legend(fontsize=fs_-5)
plt.savefig("dt_21node_sweep/dt_sweep_figs/converge_theta.png")
plt.show()

# %%
# loop through all data and plot error for cPI
plt.figure(figsize=(10, 7.5))
plt.title(r'error: average cPI', fontsize=fs_)
plt.ylabel(r'$L_2$ norm error', fontsize=fs_-5)
plt.xlabel('time step $k$ ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)

plot_err = []
plot_dx  = []

# compute 'exact' cPI

for i in range(len(time)-bump):
    
    color = 'g'
    label = 'nodes = 21'
    avg_tot_cPI_plot = avg_tot_cPI[i]

    # compute 'exact' cPI
    exact_cPI = np.array(avg_tot_cPI[exact])

    # compute error
    err = np.linalg.norm(avg_tot_cPI_plot[:60] - exact_cPI[:60])
    print('error {}: '.format(NODES[i]), err)
    plot_dx.append(DT[i])
    plot_err.append(err)
    
plt.scatter(plot_dx, plot_err, s=150, label='nodes 21')

plt.yscale('log')
plt.xscale('log')
plt.legend(fontsize=fs_-5)
plt.savefig("dt_21node_sweep/dt_sweep_figs/converge_cPI.png")
plt.show()

# %%

# loop through all data and plot cPIdot
plt.figure(figsize=(10, 7.5))
plt.title(r'error: average cPIdot', fontsize=fs_)
plt.ylabel(r'$L_2$ norm error', fontsize=fs_-5)
plt.xlabel('time step $k$ ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)

plot_err = []
plot_dx  = []

# compute 'exact' temperature
exact_cPIdot = np.array(avg_tot_cPIdot[exact])

for i in range(len(time)-bump):
    color = 'g'
    label = 'nodes = 21'

    # compute error
    err = np.linalg.norm(avg_tot_cPIdot[i][:60] - exact_cPIdot[:60])
    print('error {}: '.format(NODES[i]), err)
    plot_dx.append(DX[i])
    plot_err.append(err)

plt.scatter(plot_dx, plot_err, s=150, label='nodes 21')

plt.legend(fontsize=fs_-5)
plt.yscale('log')
plt.xscale('log')
# plt.savefig("dt_21node_sweep/dt_sweep_figs/converge_cPIdot.png")
plt.show()

# %%
# loop through all data and plot cMdot
plt.figure(figsize=(10, 7.5))
plt.title(r'error: average cMdot', fontsize=fs_)
plt.ylabel(r'$L_2$ norm error', fontsize=fs_-5)
plt.xlabel('time step $k$ ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)

plot_err = []
plot_dx  = []

# compute 'exact' cMdot
exact_cMdot = np.array(avg_tot_cMdot[exact])

for i in range(len(time)-bump):
    color = 'g'
    label = 'nodes = 21'

    # compute error
    err = np.linalg.norm(avg_tot_cMdot[i][:60] - exact_cMdot[:60])
    print('error {}: '.format(NODES[i]), err)
    plot_dx.append(DX[i])
    plot_err.append(err)

plt.scatter(plot_dx, plot_err, s=150, label='nodes 21')

plt.legend(fontsize=fs_-5)
plt.yscale('log')
plt.xscale('log')
plt.savefig("dt_21node_sweep/dt_sweep_figs/converge_cMdot.png")
plt.show()
# %%
