# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 12:01:58 2020

@author: jordy
"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import colorsys
#from matplotlib.transforms import offset_copy

plt.rcParams.update({
    'font.size': 12,
    "text.usetex": False,
    "font.family": "Arial",
    'figure.autolayout': True})

df = pd.read_csv('CombinedData.csv')
#df = pd.read_csv('data.csv')
#df['meanlength'] = df['meanlength'].fillna(0)
#df=df[df['k']==100]

labels = ['k', 'G1', 'Conc1', 'Gbb']
ids = list(range(len(labels)))
colorid = 1  # Choose which variable to label colours
markerid = 2  # Choose which variable to assign markers

# Reject groups where the mean number of pols produced was less than 10 per run
minpolcutoff = 10

ids.remove(colorid)
ids.remove(markerid)

uniques = [df[i].unique() for i in labels]

totavgdf = df.groupby(labels).agg({'meanlength': [
    'mean', 'std', 'count', 'sem'], 'numpols': ['mean', 'std', 'count', 'sem']})

#totavgdf[totavgdf.numpols['mean'].isna()] = 0
avgdf = totavgdf[totavgdf.numpols['mean'] > minpolcutoff]

avgdf = avgdf.reset_index()

plt.figure()
nandf = totavgdf[totavgdf.numpols['mean'] > minpolcutoff]

#%%

labels = ['k', 'G1', 'Conc1', 'Gbb']
axeslabels = ['k', r'$\Delta G_{spec}$', '[M]', r'$\Delta G_{BB}$']
ids = list(range(len(labels)))
subplotrowvarid = 2
subplotcolvarid = 1
colorid = 0  # Choose which variable to label colours

collabels = [
    f'{axeslabels[subplotcolvarid]} = {i}' for i in uniques[subplotcolvarid]]
rowlabels = [
    f'{axeslabels[subplotrowvarid]} = {i}' for i in uniques[subplotrowvarid]]

Ncolours = len(uniques[colorid])
HSV_tuples = [(0.0+x*0.9/Ncolours, 0.9, 0.75) for x in range(Ncolours)]
RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
colours = list(RGB_tuples)
a = 0.7

fig6 = plt.figure(figsize=(15, 8), dpi=500)
ax6 = fig6.subplots(nrows=len(rowlabels), ncols=len(
    collabels), sharex=True, sharey=True)
fig7 = plt.figure(figsize=(15, 8), dpi=500)
ax7 = fig7.subplots(nrows=len(rowlabels), ncols=len(
    collabels), sharex=True, sharey=True)
#plt.setp(ax6.flat)#, xlabel=r'$-\frac{\Delta Gpol}{k_{B}T}$', ylabel=r'$\langle L \rangle$')

pad = 1  # in points

xtickvals = [-2, 0, 2, 4, 6, 8, 10, 12, 14, 16]
lineformat = "--"
markerstyles = ['^', 'o', 's', 'P', 'd','x']
markersize = '4'

for row, subrow in enumerate(uniques[subplotrowvarid]):
    mdf = avgdf[avgdf[labels[subplotrowvarid]] == subrow]
    for col, subcol in enumerate(uniques[subplotcolvarid]):
        mdf2 = mdf[mdf[labels[subplotcolvarid]] == subcol]

        for colourvar, c, ms in zip(uniques[colorid], colours, markerstyles):
            avgsubdf = mdf2[mdf2[labels[colorid]] == colourvar]

            if (row == 1) & (col == 4):
                # Plot errorbar
                ax6[row, col].errorbar(-avgsubdf['Gbb']+np.log(avgsubdf['Conc1'])-np.log(100),
                                        avgsubdf.meanlength['mean'],
                                        avgsubdf.meanlength['sem'],
                                        color=c, alpha=a,
                                        label=fr'{colourvar:.2e}',
                                        #label=f'{colourvar}',
                                        fmt=lineformat,
                                        marker=ms,
                                        ms=markersize,
                                        elinewidth=1.5)
                
                ax7[row, col].errorbar(-avgsubdf['Gbb']+np.log(avgsubdf['Conc1'])-np.log(100),
                                       avgsubdf.numpols['mean'],
                                       avgsubdf.numpols['sem'],
                                       color=c, alpha=a,
                                       label=fr'{colourvar:.2e}',
                                       #label=f'{colourvar}',
                                       fmt=lineformat,
                                       marker=ms,
                                       ms=markersize,
                                       elinewidth=1.5)
                #Plot without errorbar
                #ax6[row,col].plot(-avgsubdf['Gbb']+np.log(avgsubdf['Conc1'])-np.log(100),avgsubdf.meanlength['mean'],color=c,alpha=0.6,label=f'{axeslabels[colorid]} = {colourvar}')
            else:
                #Plot errorbar
                ax6[row, col].errorbar(-avgsubdf['Gbb']+np.log(avgsubdf['Conc1'])-np.log(100),
                                       avgsubdf.meanlength['mean'],
                                       avgsubdf.meanlength['sem'],
                                       color=c,
                                       alpha=a,
                                       fmt=lineformat,
                                       marker=ms,
                                       ms=markersize,
                                       elinewidth=1.5)
                
                ax7[row, col].errorbar(-avgsubdf['Gbb']+np.log(avgsubdf['Conc1'])-np.log(100),
                                       avgsubdf.numpols['mean'],
                                       avgsubdf.numpols['sem'],
                                       color=c,
                                       alpha=a,
                                       fmt=lineformat,
                                       marker=ms,
                                       ms=markersize,
                                       elinewidth=1.5)
                #Plot without errorbar
                #ax6[row,col].plot(-avgsubdf['Gbb']+np.log(avgsubdf['Conc1'])-np.log(100),avgsubdf.meanlength['mean'],color=c,alpha=0.6)

        ax6[row, col].set_ylim([1.5, 4.5])
        ax6[row, col].set_yticks([2, 2.5, 3, 3.5, 4])
        ax6[row, col].grid(1)

        ax6[row, col].set_xticks(xtickvals)
        ax6[row, col].set_xticklabels([])
        
        ax7[row, col].grid(1)
        ax7[row, col].set_xticks(xtickvals)
        ax7[row, col].set_xticklabels([])
        
        if col == 0:
            ax6[row, col].set_ylabel(
                r'$\langle L \rangle$'+'\nMean length', fontsize=16)
            ax7[row, col].set_ylabel(
                r'$\langle \# pols \rangle$'+'\nAvg polymer count', fontsize=16)

        if row == len(rowlabels)-1:
            ax6[row, col].set_xlabel(
                '$-\Delta G_{pol}$'+'\nDriving strength', fontsize=16)
            ax6[row, col].set_xticks(xtickvals)
            ax6[row, col].set_xticklabels(xtickvals)
            ax7[row, col].set_xlabel(
                '$-\Delta G_{pol}$'+'\nDriving strength', fontsize=16)
            ax7[row, col].set_xticks(xtickvals)
            ax7[row, col].set_xticklabels(xtickvals)

for ax00, col00 in zip(ax6[0], collabels):
    ax00.annotate(col00, xy=(0.5, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')

for ax00, row00 in zip(ax6[:, 0], rowlabels):
    ax00.annotate(row00, xy=(0, 0.7), xytext=(-ax00.yaxis.labelpad - pad, 0),
                xycoords=ax00.yaxis.label, textcoords='offset points',
                size='large', ha='right', va='center')
    
for ax01, col01 in zip(ax7[0], collabels):
    ax01.annotate(col01, xy=(0.5, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')

for ax01, row01 in zip(ax7[:, 0], rowlabels):
    ax01.annotate(row01, xy=(0, 0.7), xytext=(-ax01.yaxis.labelpad - pad, 0),
                xycoords=ax01.yaxis.label, textcoords='offset points',
                size='large', ha='right', va='center')


fig6.subplots_adjust(wspace=0.0001)
ax6[1, 4].legend(loc=(1.1, 0), title=axeslabels[colorid]+" :", fontsize=15)

fig7.subplots_adjust(wspace=0.0001)
ax7[1, 4].legend(loc=(1.1, 0), title=axeslabels[colorid]+" :", fontsize=15)

fig6.savefig("SimpleModelParamSweepAugGspecvsM.pdf", bbox_inches="tight")
fig6.savefig("SimpleModelParamSweepAugGspecvsM.png", bbox_inches="tight")

fig7.savefig("SimpleModelParamSweepAugGspecvsMNumpols.pdf", bbox_inches="tight")
fig7.savefig("SimpleModelParamSweepAugGspecvsMNumpols.png", bbox_inches="tight")
