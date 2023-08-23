#-*- coding: utf-8 -*-
"""
Created on Thu Jan 23 17:13:01 2020
@author: jj5718
"""
import sys
import csv
import os
from os.path import dirname, abspath
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def filter_by_length(dataframe,L):
    # Function to select rows of dataframe that contain L non-NaN entries
   if (L < 1 ):
       print("L is too short")
       return
   if (L>dataframe.shape[1]):
       print("L is greater than columns in dataframe")
       return

   return dataframe[dataframe.count(axis=1)==L]

def get_start_dist(df):
    length = df.shape[1]
    cols = df.columns
    tempn = []
    tempn.append(df[df[cols[0]].isna()==False].shape[0] )
    for i in range(1,length):

        tempn.append(df[( df[cols[i-1]].isna() ) & ( df[cols[i]].isna()==False )].shape[0])
    return tempn

def get_start_dist_matrix(df):
    normedmat = []
    mat = []
    length = df.shape[1]
    for i in range(2,length+1):
        temp = np.array(get_start_dist(filter_by_length(df,i)))
        if( temp.sum() != 0):

            normedmat.append(temp/temp.sum())
            mat.append(temp)
        else:
            normedmat.append(temp)
            mat.append(temp)
    return normedmat,mat

try:
    plt.style.use('seaborn-whitegrid')
    
    filename = sys.argv[1]
    
    dirpath = dirname(abspath(filename))
    titlestr = dirpath[dirpath.rfind('/',0,dirpath.rfind('/'))+1::]
    print(titlestr)
    df = pd.read_csv(filename,engine='python')
    
    entries = df.shape[0]
    length = df.shape[1]


    if(entries != 0):
        

        lengths = []
        nsamples = []
        errors = []
        stderrors = []
    
        
        for i in range(2,length+1):
    
            dfi = filter_by_length(df,i)
            lengths.append( i)
            nsamples.append(dfi.shape[0])
            errors.append(dfi.mean(axis=1).rsub(1).mean())
            stderrors.append(dfi.mean(axis=1).rsub(1).std()/np.sqrt(dfi.shape[0]))
        
        
        plt.figure(1)
        plt.title('Samples distribution \n'+ titlestr )
        plt.plot(lengths, nsamples, linewidth=3)
        plt.xlabel('Copy length')
        plt.ylabel('Count')
    
        plt.figure(2)
        plt.title('Length distribution \n'+ titlestr )
        plt.plot(lengths, nsamples/np.sum(nsamples), linewidth=3)
        plt.xlabel('Copy length, L')
        plt.ylabel('P(Length = L)')
        
        try:
            cutoff = 10  
            
        finally:
            plt.figure(1)
            plt.savefig("LengthCountDist.pdf".format(filename), bbox_inches='tight')
    
            plt.figure(2)
            plt.ylim(0,1)
            plt.savefig("LengthProbDist.pdf".format(filename), bbox_inches='tight')
            
        nsamples = np.array(nsamples)
        arr = np.empty(np.size(nsamples))
        arr[:] = np.nan
        arr[0] = np.sum(lengths*nsamples/np.sum(nsamples))

        arr2 = np.empty(np.size(nsamples))
        arr2[:] = np.nan
        arr2[0] = np.sum(nsamples)

        CSVdata = np.transpose(np.array([lengths,nsamples, nsamples/np.sum(nsamples),errors,stderrors,arr,arr2]))
        header = ['Length','Samples','Probability','Error','Std Dev of Error','Mean length','Total samples']

    else:
        nanarr = np.empty(length-1)
        nanarr[:] = np.nan
        CSVdata = np.transpose(np.array([np.zeros(length-1),np.zeros(length-1), nanarr,nanarr,nanarr,nanarr,nanarr]))
    
        header = ['Length','Samples','Probability','Error','Std Dev of Error','Mean length','Total samples']
    
    with open('OutputAnalysis.csv','w') as csvfile:
         writer = csv.writer(csvfile, delimiter=',')
    # Gives the header name row into csv
         writer.writerow(header)
    #Data add in csv file
         for row in CSVdata:
             writer.writerow(row)

except FileNotFoundError as fnf_error:
    print(fnf_error)
except IndexError as error:
    print(error)
else:
    print('Ran with no exceptions.')
