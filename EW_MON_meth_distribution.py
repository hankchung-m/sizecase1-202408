#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 16:37:56 2023

@author: cwuhank
"""


import numpy as np
from datetime import date as dt
import netCDF4 as nc
from math import*
import matplotlib.pyplot as plt
#!!!
test=0
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

tdlat_ny=np.load('tdlat_ny.npy')
tdlon_nx=np.load('tdlon_nx.npy')
tdroci=np.load('tdroci.npy')
tdlat=np.load('tdlat.npy')
tdlon=np.load('tdlon.npy')
tdtimestr=np.load('tdtimestr.npy', allow_pickle=True)
tdhour=np.load('tdhour.npy', allow_pickle=True)
tdyear=np.load('tdyear.npy',allow_pickle=True)
name=np.load('name.npy', allow_pickle=True)
LMIR=np.load('LMIR.npy')

tdew=np.load('tdew.npy')

tdsw=np.load('tdsw.npy')
tdnw=np.load('tdnw.npy')
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

typeew4=[1,-1,-1.5,-2]
colors = ['red','purple','green','blue']

for var in range(3):
    for n in range(4):
        if var==0:
            ys = tddp[ tdew==typeew4[n] ]
            xs = tdese[ tdew==typeew4[n] ]
        elif var==1:
            #remove outlier
            tdqv[139]=np.nan
            ys = tdqv[ tdew==typeew4[n] ]
            xs = tdsw[ tdew==typeew4[n] ]
        else:
            #remove outlier
            tdqv_5_10[139]=np.nan
            tdse[84]=np.nan
            ys = tdqv_5_10[ tdew==typeew4[n] ]
            xs = tdse[ tdew==typeew4[n] ]
    
        plt.scatter(xs, ys, s=10, color=colors[n])
    
    #print(n,xs.shape)

    plt.subplots_adjust(left=0.25,right=0.75)    
    if var==0:
        plt.ylabel('DOCI (km)')
        plt.xlabel('u_ESE (m/s)')
        plt.savefig('fig/ews_ese_doci.png', dpi=300)
    elif var==1:
        plt.ylabel('q_inner (kg/kg)')
        plt.xlabel('u_SW (m/s)')
        plt.savefig('fig/ews_sw_q_inner.png', dpi=300)
    else:
        plt.ylabel('q_outer (kg/kg)')
        plt.xlabel('u_SE (m/s)')
        plt.savefig('fig/ews_se_q_outer.png', dpi=300)
    plt.show()
    
"""
#EW
for tt in range(tdew.size):
    if tdew[tt]==1 and tdese[tt]>5:
        print(tdese[tt],name[tt],tdyear[tt],tt)
        print(tdlat[tt],tdlon[tt],tdtimestr[tt])
"""
"""
#EW
for tt in range(tdew.size):
    if tdew[tt]==1 and tdqv_5_10[tt]<0.01:
        print(tdqv_5_10[tt],name[tt],tdyear[tt],tt)
        print(tdlat[tt],tdlon[tt],tdtimestr[tt])
"""
"""
#MD
for tt in range(tdew.size):
    if tdew[tt]==-2 and tddp[tt]<500:
        print(tddp[tt],name[tt],tdyear[tt],tt)
        print(tdlat[tt],tdlon[tt],tdtimestr[tt])
"""  
        
        