#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 17:40:07 2024

@author: cwuhank
"""

import numpy as np
from datetime import date as dt
import netCDF4 as nc
from math import*
import matplotlib.pyplot as plt
#!!!
test=1
should_save=0
#!!!
method=1
if method==1:
    from sklearn.cluster import KMeans as cluster
if method==2 or method==0:
    from sklearn.cluster import AgglomerativeClustering as cluster
if method==3:
    from sklearn.mixture import GaussianMixture as cluster
from sklearn.metrics import silhouette_score
import gap
from scipy import stats
from classifiability import kmeanscluster

tslat=np.load('tslat.npy')
tdlon=np.load('tdlon.npy')



tdtimestr=np.load('tdtimestr.npy', allow_pickle=True)
tdhour=np.load('tdhour.npy', allow_pickle=True)
name=np.load('name.npy', allow_pickle=True)
LMIR=np.load('LMIR.npy')

tdsws=np.load('tdsws.npy')
tdses=np.load('tdses.npy')
tdeses=np.load('tdeses.npy')
tdnn=np.load('tdnn.npy')
tdwsw=np.load('tdwsw.npy')
tdsw=np.load('tdsw.npy')
tdnw=np.load('tdnw.npy')
tdne=np.load('tdne.npy')
tdse=np.load('tdse.npy')
tdese=np.load('tdese.npy')
tddp=np.load('tddp_0.npy')
tdvo=np.load('tdvo.npy')
tdvoe=np.load('tdvoe.npy')
#tdsh=np.load('tdsh.npy')
tdwc=np.load('tdwc.npy')
tdtb=np.load('tdtb.npy')
tdqv=np.load('tdqv.npy')
tdqv_5_10=np.load('tdqv_5_10.npy')
tdrh=np.load('tdrh.npy')

tdew_w=np.load('tdew_w.npy')
tdew_m=np.load('tdew_m.npy')

"""

for tt in range(tdtimestr.size):
    if tdsw[tt]>-7 and tdsw[tt]<23 and tddp[tt]>=8 and tdwc[tt]>=10: # and tddp[tt]>0 and tddp[tt]<30 and tdwc[tt]>0 and tdwc[tt]<160 and tdvo[tt]>50 and tdvo[tt]<1200
        tdew[tt]=-2

"""



title_wm=['method_w','method_m']
title_num=[['a','b'],['c','d'],['e','f']]
file_wm=['_w','_m']
ew_num=[1,-1,-1.5,-2]
ew_n=['EW','MC','MS','MD']
var=5
var_n=['tdsw','tdse','tdese','tdnn','tdqv']

for wm in range(2):
    
    if wm==0:
        tdew=tdew_w
    else:
        tdew=tdew_m
        
    print(title_wm[wm],'------------------')
    num=tdtimestr.size

    X = np.zeros((num,var))
    NX = np.zeros((num,var))

    for tt in range(tdtimestr.size):

        X[tt,0]=tdsw[tt]#(tdese[tt]+tdeses[tt])/2#tdese[tt]#
        X[tt,1]=tdse[tt]#(tdsw[tt]+tdsws[tt])/2#tdsw[tt]#
        X[tt,2]=tdese[tt]#(tdse[tt]+tdses[tt])/2#tdse[tt]#
        
        X[tt,3]=tdnn[tt]
        X[tt,4]=tdqv[tt]



    for ew in range(4):
        print(ew_n[ew],'---------')
        for v in range(var):
            if v==4:
                r1=4
                r2=6
            else:
                r1=2
                r2=2
            print(var_n[v])
            print(round(np.mean(X[:,v][tdew==ew_num[ew]]),r1))
            print(round(np.std(X[:,v][tdew==ew_num[ew]]),r2))
            print('---------')
    

print('all---------------------')
num=0
for tt in range(tdtimestr.size):
    if tdsw[tt]>-7 and tdsw[tt]<23 and tddp[tt]>0 and tddp[tt]<3000 and tdlon[tt]>0:# and tdtb[tt]>0: # and tdwc[tt]>=0 and tdwc[tt]<200 and tdvo[tt]>50 and tdvo[tt]<1200
        num+=1
        
    #if tddp[tt]>=500:
    #    tdew[tt]=-2
#num=201 #more than 10yrs better
print(num) #376


X = np.zeros((num,var))
NX = np.zeros((num,var))
DX = np.zeros((num,2))
tdnum = np.zeros(num)
uu=0
for tt in range(tdtimestr.size):
    if tdsw[tt]>-7 and tdsw[tt]<23 and tddp[tt]>0 and tddp[tt]<3000 and tdlon[tt]>0:# and tdtb[tt]>0:

        X[uu,0]=tdsw[tt]#(tdese[tt]+tdeses[tt])/2#tdese[tt]#
        X[uu,1]=tdse[tt]#(tdsw[tt]+tdsws[tt])/2#tdsw[tt]#
        X[uu,2]=tdese[tt]#(tdse[tt]+tdses[tt])/2#tdse[tt]#
        
        X[uu,3]=tdnn[tt]
        X[uu,4]=tdqv[tt]
        DX[uu,0]=tdnn[tt]
        DX[uu,1]=tdqv[tt]
        #DX[uu,1]=tdvo[tt]-tdvoe[tt]
        #DX[uu,0]=(tdese[tt]+tdse[tt]+tdsw[tt])/3
        #X[uu,3]=tdtb[tt]
        
        #X[uu,5]=tdnw[tt]
        #X[uu,3]=abs(tdnw[tt])+abs(tdsw[tt])
        tdnum[uu]=tt
        uu+=1
        if uu==num:
            break

NX[:,0]=( X[:,0] - np.mean(X[:,0]) ) / ( np.std(X[:,0]) )
NX[:,1]=( X[:,1] - np.mean(X[:,1]) ) / ( np.std(X[:,1]) )
NX[:,2]=( X[:,2] - np.mean(X[:,2]) ) / ( np.std(X[:,2]) )
NX[:,3]=( X[:,3] - np.mean(X[:,3]) ) / ( np.std(X[:,3]) )
NX[:,4]=( X[:,4] - np.mean(X[:,4]) ) / ( np.std(X[:,4]) )
    


for v in range(var):
    if v==4:
        r1=4
        r2=6
    else:
        r1=2
        r2=2
    print(var_n[v])
    print(round(np.mean(X[:,v]),r1))
    print(round(np.std(X[:,v]),r2))
    print('---------')