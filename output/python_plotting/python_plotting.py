#%%

import glob
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
import matplotlib.cm as cm
import os
from matplotlib.ticker import ScalarFormatter


# %%
skip = 10

df = pd.read_csv("avg_concentration.csv", sep=',')
time = df['time'].tolist()[::skip]
avg_top_cPI = df[' avg_top_cPI'].tolist()[::skip]
avg_tot_cPI = df[' avg_tot_cPI'].tolist()[::skip]
avg_bot_cPI = df[' avg_bot_cPI'].tolist()[::skip]

avg_top_cPIdot = df[' avg_top_cPIdot'].tolist()[::skip]
avg_tot_cPIdot = df[' avg_tot_cPIdot'].tolist()[::skip]
avg_bot_cPIdot = df[' avg_bot_cPIdot'].tolist()[::skip]

avg_top_cMdot = df[' avg_top_cMdot'].tolist()[::skip]
avg_tot_cMdot = df[' avg_tot_cMdot'].tolist()[::skip]
avg_bot_cMdot = df[' avg_bot_cMdot'].tolist()[::skip]

avg_top_cM = df[' avg_top_cM'].tolist()[::skip]
avg_tot_cM = df[' avg_tot_cM'].tolist()[::skip]
avg_bot_cM = df[' avg_bot_cM'].tolist()[::skip]

avg_tot_theta = df[' avg_tot_theta'].tolist()[::skip]
avg_free_volume = df[' avg_free_volume'].tolist()[::skip]
avg_k_t = df[' avg_k_t'].tolist()[::skip]
avg_k_p = df[' avg_k_p'].tolist()[::skip]

avg_diff_pdot = df[' avg_diff_pdot'].tolist()[::skip]
avg_diff_mdot = df[' avg_diff_mdot'].tolist()[::skip]
avg_diff_m = df[' avg_diff_m'].tolist()[::skip]
avg_diff_theta = df[' avg_diff_theta'].tolist()[::skip]

# plot avg theta with time
fs_ = 25
plt.figure(figsize=(10, 7.5))
# plt.scatter(time, avg_top_theta, color='b', s=150, label='average top layer')
# plt.scatter(time, avg_bot_theta, color='g', s=150, label='average bottom layer')
plt.scatter(time, avg_tot_theta, color='r', s=150, label='average total')
plt.legend(fontsize=fs_-5)
plt.title(r'Average temperature', fontsize=fs_)
plt.ylabel(r'Temperature ($K$)', fontsize=fs_-5)
plt.xlabel('time ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)
# formatter = ScalarFormatter(useMathText=False)
# plt.gca().yaxis.set_major_formatter(formatter)
# plt.ylim(min(avg_tot_theta)-0.1, max(avg_tot_theta)+0.1)
# plt.ylim(303.15, 303.2)
# plt.savefig("figs/c_PI.png")
plt.show()

# plot avg cPI with time

plt.figure(figsize=(10, 7.5))
plt.scatter(time, avg_top_cPI, color='b', s=150, label='average top layer')
plt.scatter(time, avg_bot_cPI, color='g', s=150, label='average bottom layer')
plt.scatter(time, avg_tot_cPI, color='r', s=150, label='average total')
plt.legend(fontsize=fs_-5)
plt.title(r'concentration of PI ($mol/m^3$)', fontsize=fs_)
plt.ylabel(r'concentration ($mol/m^3$)', fontsize=fs_-5)
plt.xlabel('time ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)
# plt.savefig("figs/c_PI.png")  
plt.show()

# plot avg cPIdot with time
plt.figure(figsize=(10, 7.5))
plt.scatter(time, avg_top_cPIdot, color='b', s=150, label='average top layer')
plt.scatter(time, avg_bot_cPIdot, color='g', s=150, label='average bottom layer')
plt.scatter(time, avg_tot_cPIdot, color='r', s=150, label='average total')
plt.legend(fontsize=fs_-5)
plt.title(r'concentration of PIdot ($mol/m^3$)', fontsize=fs_)
plt.ylabel(r'concentration ($mol/m^3$)', fontsize=fs_-5)
plt.xlabel('time ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)
# plt.savefig("figs/c_PI.png")
plt.show()

# plot avg cMdot with time
plt.figure(figsize=(10, 7.5))
plt.scatter(time, avg_top_cMdot, color='b', s=150, label='average top layer')
plt.scatter(time, avg_bot_cMdot, color='g', s=150, label='average bottom layer')
plt.scatter(time, avg_tot_cMdot, color='r', s=150, label='average total')
plt.legend(fontsize=fs_-5)
plt.title(r'concentration of cMdot ($mol/m^3$)', fontsize=fs_)
plt.ylabel(r'concentration ($mol/m^3$)', fontsize=fs_-5)
plt.xlabel('time ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)
# plt.savefig("figs/c_PI.png")
plt.show()

# plot avg cM with time
plt.figure(figsize=(10, 7.5))
plt.scatter(time, avg_top_cM, color='b', s=150, label='average top layer')
plt.scatter(time, avg_bot_cM, color='g', marker='*', s=300, label='average bottom layer')
plt.scatter(time, avg_tot_cM, color='r', s=150, label='average total')
plt.legend(fontsize=fs_-5)
plt.title(r'concentration of cM ($mol/m^3$)', fontsize=fs_)
plt.ylabel(r'concentration ($mol/m^3$)', fontsize=fs_-5)
plt.xlabel('time ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)
plt.ylim(0, 1650)
# plt.savefig("figs/c_PI.png")
plt.show()

# plot avg free volume with time
plt.figure(figsize=(10, 7.5))
plt.scatter(time[1:], avg_free_volume[1:], color='r', s=150, label='average free volume')
plt.legend(fontsize=fs_-5)
plt.title(r'average free volume', fontsize=fs_)
plt.ylabel(r'free volume', fontsize=fs_-5)
plt.xlabel('time ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)
# plt.savefig("figs/c_PI.png")
plt.show()

# plot avg k_t with time
plt.figure(figsize=(10, 7.5))
plt.scatter(time, avg_k_t, color='r', s=150, label='average k_t')
plt.legend(fontsize=fs_-5)
plt.title(r'average $k_t$', fontsize=fs_)
plt.ylabel(r'$k_t$ ($m^3/mol/s$)', fontsize=fs_-5)
plt.xlabel('time ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)
# plt.savefig("figs/c_PI.png")
plt.show()

# plot avg k_p with time
plt.figure(figsize=(10, 7.5))
plt.scatter(time, avg_k_p, color='r', s=150, label='average k_p')
plt.legend(fontsize=fs_-5)
plt.title(r'average $k_p$', fontsize=fs_)
plt.ylabel(r'$k_p$ ($m^3/mol/s$)', fontsize=fs_-5)
plt.xlabel('time ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)
# plt.savefig("figs/c_PI.png")
plt.show()

# plot avg diff_pdot with time
plt.figure(figsize=(10, 7.5))
plt.scatter(time[1:], avg_diff_pdot[1:], color='r', s=150, label='average diff_pdot')
plt.legend(fontsize=fs_-5)
plt.title(r'average diff_pdot', fontsize=fs_)
plt.ylabel(r'diff_pdot', fontsize=fs_-5)
plt.xlabel('time ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)
# plt.savefig("figs/c_PI.png")
plt.show()

# plot avg diff_mdot with time
plt.figure(figsize=(10, 7.5))
plt.scatter(time[1:], avg_diff_mdot[1:], color='r', s=150, label='average diff_mdot')
plt.legend(fontsize=fs_-5)
plt.title(r'average diff_mdot', fontsize=fs_)
plt.ylabel(r'diff_mdot', fontsize=fs_-5)
plt.xlabel('time ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)
# plt.savefig("figs/c_PI.png")
plt.show()

# plot avg diff_M with time
plt.figure(figsize=(10, 7.5))
plt.scatter(time[1:], avg_diff_m[1:], color='r', s=150, label='average diff_M')
plt.legend(fontsize=fs_-5)
plt.title(r'average diff_M', fontsize=fs_)
plt.ylabel(r'diff_M', fontsize=fs_-5)
plt.xlabel('time ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)
# plt.savefig("figs/c_PI.png")
plt.show()

# plot avg diff_theta with time
plt.figure(figsize=(10, 7.5))
plt.scatter(time[1:], avg_diff_theta[1:], color='r', s=150, label='average diff_theta')
plt.legend(fontsize=fs_-5)
plt.title(r'average diff_theta', fontsize=fs_)
plt.ylabel(r'diff_theta ', fontsize=fs_-5)
plt.xlabel('time ($s$)', fontsize=fs_-5)
plt.xticks(fontsize=fs_-5)
plt.yticks(fontsize=fs_-5)
# plt.savefig("figs/c_PI.png")
plt.show()












# %%
