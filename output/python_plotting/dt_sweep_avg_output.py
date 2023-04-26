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
dir_path = '/Users/brianhowell/Desktop/Berkeley/MSOL/ugap_simulation/output/python_plotting/dt_21node_sweep'

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

#%%
# loop through all data and plot conversion
fs_ = 25
bump = 1
exact_21 = -1
exact_51 = -1

DT = [0.001, 0.0005, 0.0001, 5e-05, 1e-05, 1e-4]
len_data = len(time[-1])

plt.figure(figsize=(10, 7.5)) 
plt.title(r'error: conversion of M', fontsize=fs_)
plt.ylabel(r'$L_2$ norm error', fontsize=fs_-5)
plt.xlabel('time step $k$ ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)

# compute 'exact' converion
plot_err_51 = []
plot_dt_51  = []
plot_err_21 = []
plot_dt_21  = []

for i in range(len(time)-bump):
    color = 'g'
    label = 'nodes = 21'
    avg_tot_cM_plot = avg_tot_cM[i]

    exact_conv = (1 - np.array(avg_tot_cM[exact_21]) / avg_tot_cM[exact_21][0])

    conv = (1 - np.array(avg_tot_cM_plot) / avg_tot_cM_plot[0])
    err  = np.linalg.norm(conv[2:len_data] - exact_conv[2:len_data])
    print('error 21: ', err)
    plot_dt_21.append(DT[i])
    plot_err_21.append(err)

plt.scatter(plot_dt_21, plot_err_21, s=150, label='nodes 21')
plt.scatter(plot_dt_51, plot_err_51, s=150, label='nodes 51')
    
plt.legend(fontsize=10)
plt.yscale('log')
plt.xscale('log')
# plt.savefig("dt_21node_sweep/dt_sweep_figs/converge_conversion.png")
plt.show()

#%%

plot_err_51 = []
plot_dt_51  = []
plot_err_21 = []
plot_dt_21  = []

# loop through all data and plot theta
plt.figure(figsize=(10, 7.5)) 
plt.title(r'error: average temperature', fontsize=fs_)
plt.ylabel(r'$L_2$ norm error', fontsize=fs_-5)
plt.xlabel('time step $k$ ($s$)', fontsize=fs_-5)

for i in range(len(time) - bump):

    color = 'g'
    label = 'nodes = 21'
    avg_tot_theta_plot = avg_tot_theta[i]

    exact_temp = np.array(avg_tot_theta[exact_21])

    # err  = np.linalg.norm(avg_tot_theta_plot[2:len_data] - exact_temp[2:len_data])
    err = np.linalg.norm(avg_tot_theta_plot - exact_temp))
    print('error 21: ', err)
    plot_dt_21.append(DT[i])
    plot_err_21.append(err)
    


plt.scatter(plot_dt_21, plot_err_21, s=150, label='nodes 21')
plt.scatter(plot_dt_51, plot_err_51, s=150, label='nodes 51')

plt.yscale('log')
plt.xscale('log')
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)

plt.legend(fontsize=fs_-5)
plt.savefig("dt_21node_sweep/dt_sweep_figs/converge_theta.png")
plt.show()

#%%
# loop through all data and plot error for cPI
plt.figure(figsize=(10, 7.5))
plt.title(r'error: average cPI', fontsize=fs_)
plt.ylabel(r'$L_2$ norm error', fontsize=fs_-5)
plt.xlabel('time step $k$ ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)

plot_err_51 = []
plot_dt_51  = []
plot_err_21 = []
plot_dt_21  = []

# compute 'exact' cPI

for i in range(len(time)-bump):
    
    color = 'g'
    label = 'nodes = 21'
    avg_tot_cPI_plot = avg_tot_cPI[i]

    # compute 'exact' cPI
    exact_cPI = np.array(avg_tot_cPI[exact_21])

    # compute error
    # err = np.linalg.norm(avg_tot_cPI_plot[:len(exact_cPI)] - exact_cPI)
    err = np.linalg.norm(avg_tot_cPI_plot - exact_cPI)
    print('error 21: ', err)
    plot_dt_21.append(DT[i])
    plot_err_21.append(err)
    
plt.scatter(plot_dt_21, plot_err_21, s=150, label='nodes 21')
# plt.scatter(plot_dt_51, plot_err_51, s=150, label='nodes 51')

plt.yscale('log')
plt.xscale('log')
plt.legend(fontsize=fs_-5)
plt.savefig("dt_21node_sweep/dt_sweep_figs/converge_cPI.png")
plt.show()

# # loop through all data and plot cPIdot
# plt.figure(figsize=(10, 7.5))
# plt.title(r'error: average cPIdot', fontsize=fs_)
# plt.ylabel(r'$L_2$ norm error', fontsize=fs_-5)
# plt.xlabel('time step $k$ ($s$)', fontsize=fs_-5)
# plt.xticks(fontsize=fs_-5)
# plt.yticks(fontsize=fs_-5)

# # compute 'exact' temperature
# exact_cPIdot = np.array(avg_tot_cPIdot[exact])

# for i in range(len(time)-bump):
#     time_plot           = time[i]
#     avg_tot_cPIdot_plot = avg_tot_cPIdot[i]
#     dt                  = DT[i]
#     err                 = np.linalg.norm(avg_tot_cPIdot_plot - exact_cPIdot)

#     print('cPIdot err: ', err)
#     plt.scatter(dt, err, 
#                 s=150, c='g', label='time step $k$: {}$s$'.format(dt))

# plt.legend(fontsize=fs_-5)
# plt.yscale('log')
# plt.xscale('log')
# # plt.savefig("dt_21node_sweep/dt_sweep_figs/converge_cPIdot.png")
# plt.show()

# # loop through all data and plot cMdot
# plt.figure(figsize=(10, 7.5))
# plt.title(r'error: average cMdot', fontsize=fs_)
# plt.ylabel(r'$L_2$ norm error', fontsize=fs_-5)
# plt.xlabel('time step $k$ ($s$)', fontsize=fs_-5)
# plt.xticks(fontsize=fs_-5)
# plt.yticks(fontsize=fs_-5)

# # compute 'exact' temperature
# exact_cMdot = np.array(avg_tot_cMdot[exact])

# for i in range(len(time)-bump):
#     time_plot       = time[i]
#     avg_tot_cMdot_plot = avg_tot_cMdot[i]
#     dt                  = DT[i]
#     err                 = np.linalg.norm(avg_tot_cMdot_plot - exact_cMdot)

#     print('cMdot err: ', err)
#     plt.scatter(dt, err, 
#                 s=150, c='g', label='time step $k$: {}$s$'.format(dt))
    
# plt.legend(fontsize=fs_-5)
# plt.yscale('log')
# plt.xscale('log')
# plt.savefig("dt_21node_sweep/dt_sweep_figs/converge_cMdot.png")
# plt.show()




# # %%

# %%
