# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 13:46:52 2020
@author: jordy
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import os
from operator import itemgetter
def fE(inp,a,b):
    E,L = inp
    
    return 1/(1/a/np.exp(E) + 1/b )


def funcratelim(inp,r,a,b,d,e):
    E,L = inp
    
    return 1/(1+r*L/fE(inp,a,b)+fE(inp,d,e)/fE(inp,a,b))

def funcratelimnr(inp,a,b,d,e):
    E,L = inp
    
    return 1/(1+L/fE(inp,a,b)+fE(inp,d,e)/fE(inp,a,b))

def fitfunc(f,inp,p,b):
    
    fitval, fitcov = curve_fit(f, inp, p,bounds=b)#,bounds=(0, [10000, 10000, 100]))
    
    fit = f(inp,*fitval)
    return ((fitval,fitcov),fit)

def fEbig(inp,a,b,c):
    E,L = inp
    
    return 1/(1/a/np.exp(E) + 1/b + 1/c/np.exp(-E))


def funcratelimbig(inp,r,a,b,c,d,e,f):
    E,L = inp
    
    return 1/(1 +r*L/fEbig(inp,a,b,c)+fEbig(inp,d,e,f)/fEbig(inp,a,b,c))

def ratepolymerprod(inp,m,k0,ggen,gx,k,meff):
    E,L = inp
    eX = np.exp(gx)
    eGX = np.exp(gx+ggen)
    eEGX = np.exp(E+gx+ggen)
    Rform = k*m/(m+eEGX)#1/(1+(k0*eX/(k*m/(m+eEGX))))#m/(m+eEGX) #1/(1+(k0*eGX/(k*m/(m+eEGX))))
    Rend = k0*eEGX
    Rpen = k0*m*eX/(0*m+meff)
    #Rpen = 1/((np.exp(-gx) * k * (np.exp(gx) + m + meff) + k0 * (np.exp(ggen + gx) + m + np.exp(ggen) * meff))/(k * k0 * m))
    
    
    Rfall = 1/(1/Rend+1/Rpen)
    Rcomp = 1/(1/Rform+1/Rfall)
    fall = k0/(1/eX+1/eEGX)
    uniform = (L-1)*k0*eGX
    
    dimerform = 0#/(1+eEGX)
    dimerfall = k0/(1/eX+1/eEGX)
    
    pcomp = Rcomp/(Rcomp+uniform)
    return pcomp#(form*fall/(form*fall + uniform + dimerform*dimerfall))

def drazinInv(K):
    
    u,v = np.linalg.eig(K)
    nzeros = np.sum(u==0)
    Dm = np.zeros(np.shape(v))
    
    for i in range(len(u)):
        if u[i] == 0:
            Dm[i,i] = 0
        else:
            Dm[i,i] = 1/u[i]
    kdraz = np.array(v.dot(Dm.dot(np.linalg.inv(v))))
    return kdraz


def ratepolymerprod2(inp,m,k0,ggen,gx,k,meff):
    #Doesn't quite fit
    
    E,L = inp
    eX = np.exp(gx)
    eGX = np.exp(gx+ggen)
    eEGX = np.exp(E+gx+ggen)
    eEG = np.exp(E+ggen)
    
    rateform = (k0*m*k)/(k0*m+k0*eEGX+k) #Rate of completing penultimate
    pform = 1#rateform/(rateform+k0*eEGX)
    fall = k0/(1/(m*eX/(0*m+meff))+0/(eX)+1/eEGX) #Rate of complete falling off
    comprate = 1/(1/fall + 1/rateform+0/(k0*eX*m/(0*m+meff)/L)) #total rate
    
    uniform = 1*(L**1)*k0*eGX
    #uniform = 0.7*L/(1/(1*k0*eGX)+1/(k0*eX*m/(m+meff))) #Rate of release of mid length polymers
    
    
    R1 = k0*m*eEGX/meff/(1+eEG+meff*eEG/m) # Rate of going from flat end of template to two monomers at end
    R2 = k*m/(m+eEGX)
    
    R1 = k0*m*eEGX/meff # Rate of going from flat end of template to two monomers at end
    R2 = k0*m*eX/meff
    
    dimerrateform = 1/(1/R1+1/R2)
    #dimerrateform = (k0*m*k)/(k0*m+k0*eEGX+k) #Rate of forming dimer
    dimerrateform = 1/(1/R1+1/R2+1/(k/(m/(m+eEGX))))
    
    dimerfall = k0/(0/(m*eX/(m+meff))+1/(eX)+1/eEGX) # Rate of dimer release
    
    penultfillrate = 1*k0*m*eX/meff#k0*eEGX*eX*m/(eEGX*(eX+meff+m)+meff*(meff+m)+eX*m)
    dimerrate = 1/(1/dimerfall + 1/dimerrateform+ 0/penultfillrate)
    pdimerform = dimerrate/(dimerrate+1*fall+0*m*eX/(m+meff))
    pdimerform = 0#*dimerfall/(dimerfall+1*fall+0*m*eX/(m+meff))
    
    return (pform*comprate/(pform*comprate + uniform + pdimerform*dimerrate))
              

def ratepolymerprod3(inp,m,k0,ggen,gx,k,meff):
    #Doesn't quite fit
    E,L = inp
    eX = np.exp(gx)
    eGX = np.exp(gx+ggen)
    eEGX = np.exp(E+gx+ggen)
    eEG = np.exp(E+ggen)
    
    tfall2Drazin = np.zeros(len(E))
    for i,e in enumerate(E):
        K = np.zeros((6,6))
        teX = np.exp(gx)
        teGX = np.exp(gx+ggen)
        teEGX = np.exp(e+gx+ggen)
        
        #K[to,from]
        K[1,0] = teX
        K[3,0] = teEGX
        K[0,1] = k0*meff
        K[2,1] = k0*meff#m
        K[5,1] = teEGX
        K[1,2] = teX#k0*np.exp(gX)#+gG)teGX
        K[5,2] = teEGX
        
        K[0,3] = k0*meff
        K[4,3] = k0*m
        K[5,3] = teX
        K[3,4] = teEGX
        K[5,4] = teX
    
        K = K + np.diag(-np.sum(K,axis=0))
        
        kdraztemp2 = drazinInv(K)
        tfall2Drazin[i] = kdraztemp2[5,0]
    
    rateform = (k0*m*k)/(0*k0*m+k0*eEGX+0*k) #Rate of completing penultimate
    
    fall = k0/(1/(m*eX/(m+meff))+0/(eX)+1/eEGX) #Rate of complete falling off
    #fall = 1/tfall2Drazin
    comprate = 1/(1/fall + 1/rateform+0/(k0*eX*m/(m+0*meff)/L)) #total rate
    
    uniform = (L**1)*k0*eGX
    #uniform = 0.7*L/(1/(1*k0*eGX)+1/(k0*eX*m/(m+meff))) #Rate of release of mid length polymers
    
    
    R1 = k0*m*eEGX/meff/(1+eEG+meff*eEG/m) # Rate of going from flat end of template to two monomers at end
    R2 = k*m/(m+eEGX)
    
    R1 = k0*m*eEGX/meff # Rate of going from flat end of template to two monomers at end
    #R2 = k0*m*eX/meff
    
    dimerrateform = 1/(1/R1+1/R2)
    dimerrateform = 1/(1/m+1/R2)
    #dimerrateform = (k0*m*k)/(k0*m+k0*eEGX+k) #Rate of forming dimer
    dimerrateform = 1/(1/R1+1/R2+0/(k/(m/(m+eEGX))))
    
    dimerfall = k0/(0/(m*eX/(m+meff))+1/(eX)+1/eEGX) # Rate of dimer release
    dimerfall = 1/tfall2Drazin # Rate of dimer release
    
    dimerrate = 1/(1/dimerfall + 1/dimerrateform)
    
    pdimerform = dimerrate/(dimerrate+1*fall+0*m*eX/(m+meff))
    #pdimerform = dimerfall/(dimerfall+1*fall)
    
    return (comprate/(comprate + uniform + pdimerform*dimerrate))
              
    
    
rootdir = r'C:\Users\jordy\OneDrive - Imperial College London\CPlusPlus\PPC Simulation 2020\Oct2020\Gend_T_2'

plt.rcParams.update({
    'font.size': 14,
    "text.usetex": False,
    'lines.linewidth': 1.5,
    "font.sans-serif": ["Arial"],
    'figure.autolayout': True})

alldata = pd.DataFrame([])
inputparamnames = []
inputparams = []
datasets = []

datalist = []

for subdir, dirs, files in os.walk(rootdir):
    #print(subdir)
    for dir in dirs:
        try:        
            inputlist = pd.read_csv(os.path.join(dir,'inputlist.csv'),nrows=1,usecols=(['length_limit','k0','k','G1','ConcEff','Conc1','Gbb','Ggen']))
            inputparams.append(inputlist)
            avgdata = pd.read_csv(os.path.join(dir,'Avgdata.csv'))        
            #inputparams.append(list(inputlist.values[0]))
            #datasets.append(avgdata)
            templist = list(inputlist.values[0])
            templist.append(avgdata)
            datalist.append(templist)
        except:
            continue
            
k = 1
k0 = 1
G1 = -4
Meff = 100
    
inputparamnames = inputlist.columns.tolist()
#print(f for f in files for files in os.walk(rootdir))
#df = pd.concat(map(pd.read_csv,glob.glob('*/Avgdata.csv')))

datalist.sort(key=itemgetter(5,7,0))
uniqueparams = list(list(np.unique(a)) for a in np.transpose(datalist)[0:8])

fig, axes = plt.subplots(3,3,dpi=400,figsize=(12, 8),sharex=True,sharey=True)

plt.setp(axes.flat)
fig.subplots_adjust(hspace=0.7,wspace=0.3)
pad = 5 # in points

labels = ['k','G1','Conc1','Gbb']
axeslabels = ['k',r'$\Delta G_{spec}$','[M]',r'$\Delta G_{BB}$']
pad = 4

fig2, ax2 = plt.subplots(1,1,dpi=400,figsize=(10, 8))


collabels = [r'$\Delta G_{gen} = -12 $',r'$\Delta G_{gen} = -10 $',r'$\Delta G_{gen} = -8 $']
rowlabels = [r'$[M] = 1$',r'$[M] = 3$',r'$[M] = 10$']


params_data=zip(inputparams,datasets)
mfocus = 1
ggenfocus = -12

linecols = ["#23EAFF","#0377C3","#00064D"]


for row,M in enumerate(uniqueparams[5]):
    for col,Ggen in enumerate(uniqueparams[7]):
        axes[row,col].set_xlim([0,40])
        axes[row,col].set_ylim([0,1.05])
        axes[row,col].grid(1)
        #axes[row,col].set_title("[M] = %d , Ggen = %d kbT" % (M,Ggen))
        axes[row,col].set_xticklabels([])
        axes[row,col].set_yticks([0,0.2,0.4,0.6,0.8,1])    
        axes[row,col].set_xticks([0,5,10,15,20,25,30,35,40])
        
        
        if col == 0:
            axes[row,col].set_ylabel(r'P(complete)',fontsize=18)
            
        if row == 2:
            axes[row,col].set_xlabel(r'$\Delta G_{end}$',fontsize=18)
            axes[row,col].set_xticklabels([0,5,10,15,20,25,30,35,40])
        
        
        if (M==mfocus) & (Ggen== ggenfocus):
            ax2.set_xlim([0,20])
            ax2.set_ylim([0,1])
            ax2.set_title("[M] = %d , Ggen = %d kbT" % (M,Ggen))
            
            #ax2.plot([-Ggen,-Ggen],[0,1])
            #ax2.plot([np.log(M)-Ggen-G1,np.log(M)-Ggen-G1],[0,1],'r')
            #ax2.plot([np.log(Meff)-Ggen-G1,np.log(Meff)-Ggen-G1],[0,1],'r')
            #ax2.plot([np.log(Meff)-Ggen-2*G1,np.log(Meff)-Ggen-2*G1],[0,1],'g')
            #ax2.plot([np.log(M*k/k0)-2*Ggen-2*G1,np.log(M*k/k0)-2*Ggen-2*G1],[0,1],'b')
        
        
        p_multi = np.array([])
        E_multi = np.array([])
        T_multi = np.array([])
        
        for data in datalist:
            if (data[5]==M) & (data[7]== Ggen):
                T = data[0]
                #Template length based line colour
                if (T==10):
                    linecol = linecols[0]
                elif (T==30):
                    linecol = linecols[1]
                elif (T==100):
                    linecol = linecols[2]
                
                
                if (data[5]==mfocus) & (data[7]== ggenfocus):
                    ax2.plot(data[8]['Gend'],data[8]['Complete probability'],linewidth=3,color=linecol)
                    cutoff = int(-6*Ggen+8)
                    
                    E = data[8]['Gend'].values[0:cutoff] 
                    p = data[8]['Complete probability'].values[0:cutoff] 
                    inp = (E,np.repeat(T,len(E)))
                    #ax2.plot(E,ratepolymerprod(inp,M,k0,Ggen,G1,k,Meff),':',linewidth=5)
                                        
                    ax2.set_xlabel("Gend, kT ",fontsize=18)
                    ax2.set_ylabel("P(Complete)",fontsize=18)
                    
                    ax2.grid(1)
                    
                cutoff = int(-6*Ggen+8)
                
                E = data[8]['Gend'].values[0:cutoff] 
                p = data[8]['Complete probability'].values[0:cutoff] 
                inp = (E,np.repeat(T,len(E)))
                
                p_multi = np.append(p_multi,p)
                E_multi = np.append(E_multi,E)
                T_multi = np.append(T_multi,np.repeat(T,len(E)))    
                
                if (T==10):
                    l1, = axes[row,col].plot(data[8]['Gend'],data[8]['Complete probability'],color=linecol,alpha=0.6)
                    #f1, = axes[row,col].plot(E,ratepolymerprod3(inp,M,k0,Ggen,G1,k,Meff),':',color=linecol)
                    f1, = axes[row,col].plot(E,ratepolymerprod(inp,M,k0,Ggen,G1,k,Meff),':',color=linecol)
                elif (T==30):
                    l2, = axes[row,col].plot(data[8]['Gend'],data[8]['Complete probability'],color=linecol,alpha=0.6)
                    #f2, = axes[row,col].plot(E,ratepolymerprod3(inp,M,k0,Ggen,G1,k,Meff),':',color=linecol)
                    f2, = axes[row,col].plot(E,ratepolymerprod(inp,M,k0,Ggen,G1,k,Meff),':',color=linecol)
                elif (T==100):
                    l3, = axes[row,col].plot(data[8]['Gend'],data[8]['Complete probability'],color=linecol,alpha=0.6)
                    #f3, = axes[row,col].plot(E,ratepolymerprod3(inp,M,k0,Ggen,G1,k,Meff),':',color=linecol)
                    f3, = axes[row,col].plot(E,ratepolymerprod(inp,M,k0,Ggen,G1,k,Meff),':',color=linecol)
                
                #axes[row,col].plot([0,40],np.repeat(1/(1 + T*Meff*np.exp(Ggen)/M),2),color=linecol)
                #axes[row,col].plot(np.repeat(np.log(T),2),[0,1],color=linecol)
                #axes[row,col].plot(np.repeat(np.log(M*k/(T-1)/k0)-2*Ggen-2*G1,2),[0,1],color=linecol)
                try:
                    1/0
                    #((fitval,fitcov),fit) = fitfunc(funcratelimnr,inp,p,(0,[100,1000,100,100]))
                    #((fitval,fitcov),fit) = fitfunc(funcratelim,inp,p,(0,[1000,1000,10000,1000,1000]))
                    #((fitval,fitcov),fit) = fitfunc(funcratelimbig,inp,p,(0,[0.01,0.01,10,np.exp(35),0.00001,10,np.exp(27)]))
                    #axes[row,col].plot(E,fit,'k:',label=None)
                
                    #print('M = {}, Ggen = {}, T = {} , fit = {}'.format(M,Ggen,T,fitval))
                except:
                    #print('Fit not found.')
                    continue
        
        #c1, = axes[row,col].plot(np.repeat(np.log(T),2),[0,1])
        c2, = axes[row,col].plot(np.repeat(np.log(M/Meff*np.exp(G1))-Ggen-G1,2),[0,1],'r')
        c3, = axes[row,col].plot(np.repeat(np.log(Meff)-Ggen-G1,2),[0,1],'g')
        #c4, = axes[row,col].plot([np.log(Meff)-Ggen-2*G1,np.log(Meff)-Ggen-2*G1],[0,1],'k')
        #c5, = axes[row,col].plot([np.log(M*k/k0)-2*Ggen-2*G1,np.log(M*k/k0)-2*Ggen-2*G1],[0,1],'k')
        # c1, = axes[row,col].plot([-Ggen,-Ggen],[0,1])
        # c2, = axes[row,col].plot([np.log(M)-Ggen-G1,np.log(M)-Ggen-G1],[0,1],'k')
        # c3, = axes[row,col].plot([np.log(Meff)-Ggen-G1,np.log(Meff)-Ggen-G1],[0,1],'k')
        # c4, = axes[row,col].plot([np.log(Meff)-Ggen-2*G1,np.log(Meff)-Ggen-2*G1],[0,1],'k')
        # c5, = axes[row,col].plot([np.log(M*k/k0)-2*Ggen-2*G1,np.log(M*k/k0)-2*Ggen-2*G1],[0,1],'k')
        
        try:
            inp_multi = (E_multi,T_multi)
            #((fitval_multi,fitcov_multi),fit_multi) = fitfunc(funcratelimbig,inp_multi,p_multi,[[0,0,0,100,0,0,1000],[0.01,0.01,30,10**16,0.000001,30,10**14]])
            #fit_10 = fit_multi[0:cutoff]
            #fit_30 = fit_multi[cutoff:2*cutoff]
            #fit_100 = fit_multi[2*cutoff:3*cutoff]
            
            #axes[row,col].plot(E,fit_10,'r',label=None,alpha=0.7)
            #axes[row,col].plot(E,fit_30,'r',label=None,alpha=0.7)
            #axes[row,col].plot(E,fit_100,'r',label=None,alpha=0.7)
            
            
            #print('[M] = {}, Ggen = {}, fit_multi = {}'.format(M,Ggen,fitval_multi))
        except:
            #print('Fit not found.')
            continue
        
    for ax, col in zip(axes[0], collabels):
        ax.annotate(col, xy=(0.5, 1), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')

    for ax, row in zip(axes[:,0], rowlabels):
        ax.annotate(row, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                    xycoords=ax.yaxis.label, textcoords='offset points',
                    size='large', ha='right', va='center')
                
'''    
for data in datalist:
    for row,M in enumerate(uniqueparams[5]):
        for col,Ggen in enumerate(uniqueparams[7]):
            if (data[5]==M) & (data[7]== Ggen):
                
                axes[row,col].plot(data[8]['Gend'],data[8]['Complete probability'])
                axes[row,col].set_xlim([0,40])
                axes[row,col].set_ylim([0,1])
                axes[row,col].set_title("M = %d , Ggen = %d" % (M,Ggen))
                
                cutoff = int(-2*Ggen+8)
                T = data[0]
                E = data[8]['Gend'].values[0:cutoff] 
                p = data[8]['Complete probability'].values[0:cutoff] 
                inp = (E,np.repeat(T,len(E)))
                
                try:
                    #((fitval,fitcov),fit) = fitfunc(funcratelimnr,inp,p,([0,0,0,0],[100,1000,100,100]))
                    ((fitval,fitcov),fit) = fitfunc(funcratelim,inp,p,([0,0,0,0,0],[2,1,40,1,40]))
                    #((fitval,fitcov),fit) = fitfunc(funcratelimbig,inp,p,([0,0,0,np.exp(11),0,0,np.exp(8)],[0.01,0.01,10,np.exp(35),0.00001,10,np.exp(27)]))
                    axes[row,col].plot(E,fit,'k:',label=None)
                
                    print('M = {}, Ggen = {}, T = {} , fit = {}'.format(M,Ggen,T,fitval))
                except:
                    print('Fit not found.')
                    continue
                
'''                                
fig.legend((l1, f1,l2, f2,l3, f3,c2,c3),
           labels=['L = 10','','L = 30','','L = 100','Analytic curve',r'$\tau_{end} = \tau_{disp}$',r'$\Delta G_{end} = -\Delta G_{gen}-\Delta G_{spec}+\ln [M_{eff}]$'],   # The labels for each line
           bbox_to_anchor=(1.35, 0.7, 0, 0),   # Position of legend
           borderaxespad=0.1,    # Small spacing around legend box
           title="Legend"  # Title for the legend
           )

#fig2.legend(labels=['L = 10','L = 30','L = 100'],title="Template lengths")

plt.subplots_adjust(right=0.85)
plt.show()
fig.savefig('GENDMultipleMGgen.pdf', bbox_inches = "tight")  
#fig.savefig('GENDMultipleMGgen.png', bbox_inches = "tight")  
#fig2.savefig('FNANOEND.eps')  
