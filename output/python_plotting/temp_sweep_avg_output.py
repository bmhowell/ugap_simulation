#%%

import glob
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
import matplotlib.cm as cm
import os
from matplotlib.ticker import ScalarFormatter

#%%
skip = 1

# central comuting 
# dir_path = '/home/brian/Documents/berkeley/ugap_simulation/output/python_plotting/temp_sweep'

# macbook computing
dir_path = '/Users/brianhowell/Desktop/Berkeley/MSOL/ugap_simulation/output/python_plotting/temp_sweep'

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

    avg_tot_cPI.append(df[' avg_tot_cPI'].tolist()[::skip])
    avg_tot_cPIdot.append(df[' avg_tot_cPIdot'].tolist()[::skip])
    avg_tot_cMdot.append(df[' avg_tot_cMdot'].tolist()[::skip])
    avg_tot_cM.append(df[' avg_tot_cM'].tolist()[::skip])
    avg_tot_theta.append(df[' avg_tot_theta'].tolist()[::skip])
    avg_free_volume.append(df[' avg_free_volume'].tolist()[::skip])
    avg_k_t.append(df[' avg_k_t'].tolist()[::skip])
    avg_k_p.append(df[' avg_k_p'].tolist()[::skip])

    counter += 1

# loop through all data and plot conversion
fs_ = 25
bump = 0
I_UV = [2, 6, 10, 50, 100]
temp_amb = [298.15, 303.15, 315.15, 325.15, 335.15]
plt.figure(figsize=(10, 7.5)) 
plt.title(r'conversion of M', fontsize=fs_)
plt.ylabel(r'conversion (%)', fontsize=fs_-5)
plt.xlabel('time ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)
plt.ylim(0, 1)
for i in range(len(time)-bump):
    time_plot       = time[i]
    avg_tot_cM_plot = avg_tot_cM[i]
    
    plt.scatter(time_plot, (1 - np.array(avg_tot_cM_plot) / avg_tot_cM_plot[0]), 
                s=150, label=r'$\theta_0$: {}$K$'.format(temp_amb[i]))
    
    
plt.legend(fontsize=fs_-5)
plt.savefig("temp_sweep/temp_sweep_fig/multi_conversion.png")
plt.show()


# loop through all data and plot theta
plt.figure(figsize=(10, 7.5)) 
plt.title(r'Average temperature', fontsize=fs_)
plt.ylabel(r'Temperature ($K$)', fontsize=fs_-5)
plt.xlabel('time ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)
# plt.ylim(303.15, 303.2)
for i in range(len(time) - bump):
    time_plot       = time[i]
    avg_tot_theta_plot = avg_tot_theta[i]
    
    plt.scatter(time_plot, avg_tot_theta_plot, 
                s=150, label=r'$\theta_0$: {}$K$'.format(temp_amb[i]))
    
    
plt.legend(fontsize=fs_-5)
plt.savefig("temp_sweep/temp_sweep_fig/multi_theta.png")
plt.show()

# loop through all data and plot cPI
plt.figure(figsize=(10, 7.5))
plt.title(r'Average cPI', fontsize=fs_)
plt.ylabel(r'cPI ($mol/m^3$)', fontsize=fs_-5)
plt.xlabel('time ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)

for i in range(len(time)-bump):
    time_plot       = time[i]
    avg_tot_cPI_plot = avg_tot_cPI[i]
    
    plt.scatter(time_plot, np.array(avg_tot_cPI_plot), 
                s=150, label=r'$\theta_0$: {}$K$'.format(temp_amb[i]))
    
plt.legend(fontsize=fs_-5)
plt.savefig("temp_sweep/temp_sweep_fig/multi_cPI.png")
plt.show()

# loop through all data and plot cPIdot
plt.figure(figsize=(10, 7.5))
plt.title(r'Average cPIdot', fontsize=fs_)
plt.ylabel(r'cPIdot ($mol/m^3$)', fontsize=fs_-5)
plt.xlabel('time ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)

for i in range(len(time)-bump):
    time_plot       = time[i]
    avg_tot_cPIdot_plot = avg_tot_cPIdot[i]
    
    plt.scatter(time_plot, np.array(avg_tot_cPIdot_plot), 
                s=150, label=r'$\theta_0$: {}$K$'.format(temp_amb[i]))

plt.legend(fontsize=fs_-5)
plt.savefig("temp_sweep/temp_sweep_fig/multi_cPIdot.png")
plt.show()

# loop through all data and plot cMdot
plt.figure(figsize=(10, 7.5))
plt.title(r'Average cMdot', fontsize=fs_)
plt.ylabel(r'cMdot ($mol/m^3$)', fontsize=fs_-5)
plt.xlabel('time ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)

for i in range(len(time)-bump):
    time_plot       = time[i]
    avg_tot_cMdot_plot = avg_tot_cMdot[i]
    
    plt.scatter(time_plot, np.array(avg_tot_cMdot_plot), 
                s=150, label=r'$\theta_0$: {}$K$'.format(temp_amb[i]))
    
plt.legend(fontsize=fs_-5)
plt.savefig("temp_sweep/temp_sweep_fig/multi_cMdot.png")
plt.show()

# loop through all data and plot cMdot
plt.figure(figsize=(10, 7.5))
plt.title(r'Average cM', fontsize=fs_)
plt.ylabel(r'cM ($mol/m^3$)', fontsize=fs_-5)
plt.xlabel('time ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)

for i in range(len(time)-bump):
    time_plot       = time[i]
    avg_tot_cM_plot = avg_tot_cM[i]
    
    plt.scatter(time_plot, np.array(avg_tot_cM_plot), 
                s=150, label=r'$\theta_0$: {}$K$'.format(temp_amb[i]))
    
plt.legend(fontsize=fs_-5)
plt.savefig("temp_sweep/temp_sweep_fig/multi_cM.png")
plt.show()

# loop through all data and plot k_t
plt.figure(figsize=(10, 7.5))
plt.title(r'Average $k_t$', fontsize=fs_)
plt.ylabel(r'$k_t$ ($m^3/mol/s$)', fontsize=fs_-5)
plt.xlabel('time ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)

for i in range(len(time)-bump):
    time_plot       = time[i]
    avg_k_t_plot = avg_k_t[i]
    
    plt.scatter(time_plot, np.array(avg_k_t_plot), 
                s=150, label=r'$\theta_0$: {}$K$'.format(temp_amb[i]))
    
plt.legend(fontsize=fs_-5)
plt.savefig("temp_sweep/temp_sweep_fig/multi_k_t.png")
plt.show()

# loop through all data and plot k_p
plt.figure(figsize=(10, 7.5))
plt.title(r'Average $k_p$', fontsize=fs_)
plt.ylabel(r'$k_p$ ($m^3/mol/s$)', fontsize=fs_-5)
plt.xlabel('time ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)

for i in range(len(time)-bump):
    time_plot       = time[i]
    avg_k_p_plot = avg_k_p[i]
    
    plt.scatter(time_plot, np.array(avg_k_p_plot), 
                s=150, label=r'$\theta_0$: {}$K$'.format(temp_amb[i]))

plt.legend(fontsize=fs_-5)
plt.savefig("temp_sweep/temp_sweep_fig/multi_k_p.png")
plt.show()


# %%
