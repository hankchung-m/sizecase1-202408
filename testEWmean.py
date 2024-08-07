#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 13:55:00 2024

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

tdlat_ny=np.load('tdlat_ny.npy')
tdlon_nx=np.load('tdlon_nx.npy')
tdroci=np.load('tdroci.npy')

tdtimestr=np.load('tdtimestr.npy', allow_pickle=True)
tdhour=np.load('tdhour.npy', allow_pickle=True)
name=np.load('name.npy', allow_pickle=True)
LMIR=np.load('LMIR.npy')


tdxx=np.load('tdese.npy')

"""

tdew=np.load('tdew_w.npy')

print(np.mean(tdxx[tdew==1]))
print(np.mean(tdxx[tdew==-1]))
print(np.mean(tdxx[tdew==-1.5]))
print(np.mean(tdxx[tdew==-2]))

"""

tdew_w=np.load('tdew_w.npy')
tdew_m=np.load('tdew_m.npy')
ms_m=0
ms_mw=0
for tt in range(tdew_w.size):
    if tdew_m[tt]==-1.5:
        ms_m+=1
        if tdew_w[tt]==-2:
            ms_mw+=1
print(ms_m,ms_mw,ms_mw/ms_m)