#-*- coding: utf-8 -*-
"""
Created on Thu Jan 23 17:13:01 2020
@author: jj5718
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
plt.style.use('seaborn-whitegrid')

from matplotlib.ticker import StrMethodFormatter
from matplotlib import ticker, cm
from matplotlib.cm import ScalarMappable

plt.rcParams.update({
    'font.size': 12,
    "text.usetex": False,
    "font.family": "Helvetica",
    'figure.autolayout': True})

tag = '_Aug2021'
folder1 = 'Empty/'
folder2 = 'Full/'

plt.figure(dpi = 400, figsize = (5,5))
print("Reading data...")
#for i in range(122,141):
#for i in range(141,160):
# inputlist1 = pd.read_csv(f'inputlist_old.csv')

inputlist1 = pd.read_csv(folder1+'inputlist.csv')
inputlist1["Avg # TC bonds empty"] = np.nan
inputlist1["Avg # total objs empty"] = np.nan
inputlist1["Avg # Polymers empty"] = np.nan
inputlist1["Avg # TC bonds full"] = np.nan
inputlist1["Avg # total objs full"] = np.nan
inputlist1["Avg # Polymers full"] = np.nan


#%%
avgslist = np.zeros((6,inputlist1.shape[0]))
for i in range(inputlist1.shape[0]):
#for i in range(122,141):
    print(i)
    tempdf0 = pd.read_csv(folder1+f'{i+1}/traj.csv',engine='python',usecols=['Time (s)','# TC bonds','# total','# Polymers'])
    tempdf2 = pd.read_csv(folder2+f'{i+1}/traj.csv',engine='python',usecols=['Time (s)','# TC bonds','# total','# Polymers'])
    avgslist[:,i] = np.array([tempdf0["# TC bonds"].mean(),
                                     tempdf0["# total"].mean(),
                                     tempdf0["# Polymers"].mean(),
                                     tempdf2["# TC bonds"].mean(),
                                     tempdf2["# total"].mean(),
                                     tempdf2["# Polymers"].mean()])

    #cutoff = -1
    #plt.plot(tempdf['Time (s)'][0:cutoff],tempdf['# TC bonds'][0:cutoff],alpha = 0.6,linewidth = 0.5)
#%%

inputlist1["Avg # TC bonds empty"] = avgslist[0,:]
inputlist1["Avg # total objs empty"] = avgslist[1,:]
inputlist1["Avg # Polymers empty"] = avgslist[2,:]
inputlist1["Avg # TC bonds full"] = avgslist[3,:]
inputlist1["Avg # total objs full"] = avgslist[4,:]
inputlist1["Avg # Polymers full"] = avgslist[5,:]

l10 = inputlist1[inputlist1['length_limit']==10]
l30 = inputlist1[inputlist1['length_limit']==30]
avgl10 = l10.groupby(["G0"]).agg({"Avg # TC bonds empty"   : ['mean', 'std', 'count', 'sem'], 
                                  "Avg # total objs empty" : ['mean', 'std', 'count', 'sem'],
                                  "Avg # Polymers empty"   : ['mean', 'std', 'count', 'sem'],
                                  "Avg # TC bonds full"    : ['mean', 'std', 'count', 'sem'],
                                  "Avg # total objs full"  : ['mean', 'std', 'count', 'sem'],
                                  "Avg # Polymers full"    : ['mean', 'std', 'count', 'sem']})

avgl30 = l30.groupby(["G0"]).agg({"Avg # TC bonds empty"   : ['mean', 'std', 'count', 'sem'], 
                                  "Avg # total objs empty" : ['mean', 'std', 'count', 'sem'],
                                  "Avg # Polymers empty"   : ['mean', 'std', 'count', 'sem'],
                                  "Avg # TC bonds full"    : ['mean', 'std', 'count', 'sem'],
                                  "Avg # total objs full"  : ['mean', 'std', 'count', 'sem'],
                                  "Avg # Polymers full"    : ['mean', 'std', 'count', 'sem']})

avgl10.reset_index(level = 0,inplace = True)
avgl30.reset_index(level = 0,inplace = True)

#%%
n=20
idx0 = 461

G3 = inputlist1[inputlist1['G0']==3]
G2_75 = inputlist1[inputlist1['G0']==2.75]
#%%
fig,ax = plt.subplots(nrows = n,ncols = 2,dpi = 400, figsize = (10,10),sharex =True)
for j,i in enumerate(range(idx0,idx0+n)):
    traj3_10_empty = pd.read_csv(folder1+f'{i}/traj.csv',engine='python',usecols=['Time (s)','# TC bonds','# total','# Polymers'])
    traj3_10_full = pd.read_csv(folder2+f'{i}/traj.csv',engine='python',usecols=['Time (s)','# TC bonds','# total','# Polymers'])

    cutoff = -1
    ax[j,0].plot(traj3_10_empty['Time (s)'][0:cutoff],traj3_10_empty['# TC bonds'][0:cutoff],alpha = 0.6,label = f'Empty IC {i}')
    ax[j,1].plot(traj3_10_full['Time (s)'][0:cutoff],traj3_10_full['# TC bonds'][0:cutoff],alpha = 0.6,label = f'Full IC {i}')
    
    ax[j,0].set_ylabel("# CT bonds")

    
ax[n-1,0].set_xlabel("Time")
ax[n-1,1].set_xlabel("Time")
#ax.legend()
plt.show()
#%%
fig,ax = plt.subplots(nrows = 2, ncols = 1,dpi = 400, figsize = (6,6))
# make a plot
ax[0].errorbar(avgl10['G0'].values[:],avgl10[('Avg # TC bonds empty','mean')].values[:]/10,avgl10[('Avg # TC bonds empty','sem')].values[:]/10, 
            color = 'r',ls='-',marker = 'x',alpha = 0.5,label = 'L = 10 IC empty')
ax[0].errorbar(avgl10['G0'].values[:],avgl10[('Avg # TC bonds full','mean')].values[:]/10,avgl10[('Avg # TC bonds full','sem')].values[:]/10, 
            color = 'r',ls='-',marker = 'o',mfc = 'none', alpha = 0.5,label = 'L = 10 IC full')
ax[0].errorbar(avgl30['G0'].values[:],avgl30[('Avg # TC bonds empty','mean')].values[:]/30,avgl30[('Avg # TC bonds empty','sem')].values[:]/30, 
            color = 'b',ls='-',marker = 'x',alpha = 0.5,label = 'L = 30 IC empty')
ax[0].errorbar(avgl30['G0'].values[:],avgl30[('Avg # TC bonds full','mean')].values[:]/30,avgl30[('Avg # TC bonds full','sem')].values[:]/30, 
            color = 'b',ls='-',marker = 'o',mfc = 'none',alpha = 0.5,label = 'L = 30 IC full')

L = 30
edg = np.exp(avgl10['G0'].values[:]-np.log(45))

ax[0].set_xlabel(r"$\Delta G_{spec}$",fontsize=14)
ax[0].set_ylabel("Avg # TC bonds/L",fontsize=14)
ax[0].set_ylim([0,1.1])
#ax[0].set_yscale("log")
ax[0].legend()

ax[1].errorbar(avgl10['G0'].values[:],avgl10[('Avg # total objs empty','mean')].values[:],avgl10[('Avg # total objs empty','sem')].values[:],
             color="r",ls = "-",marker = 'x',alpha = 0.7,label = 'L = 10 IC empty')
ax[1].errorbar(avgl10['G0'].values[:],avgl10[('Avg # total objs full','mean')].values[:],avgl10[('Avg # total objs full','sem')].values[:],
             color="r",ls = "-",marker = 'o',mfc = 'none',alpha = 0.7,label = 'L = 10 IC full')
ax[1].errorbar(avgl30['G0'].values[:],avgl30[('Avg # total objs empty','mean')].values[:],avgl30[('Avg # total objs empty','sem')].values[:], 
             color="b",ls = "-",marker = 'x',alpha = 0.7,label = 'L = 30 IC empty')
ax[1].errorbar(avgl30['G0'].values[:],avgl30[('Avg # total objs full','mean')].values[:],avgl30[('Avg # total objs full','sem')].values[:], 
             color="b",ls = "-",marker = 'o',mfc = 'none',alpha = 0.7,label = 'L = 30 IC full')
ax[1].set_ylabel("Avg # objs",fontsize=14)
ax[1].set_ylim([0,1.3])

ax[1].legend()

ax[1].set_xlabel(r"$\Delta G_{spec}$",fontsize=14)
plt.show()

fig.savefig("SimpleAvgNumTCBonds_Aug2021.pdf", bbox_inches="tight")
