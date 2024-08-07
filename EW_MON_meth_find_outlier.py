#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 23:58:55 2023

@author: cwuhank
"""


import numpy as np
from datetime import date as dt
import netCDF4 as nc
from math import*
import matplotlib.pyplot as plt
tdese=np.load('tdese.npy')
tddp=np.load('tddp_0.npy')
tdew=np.load('tdew_m.npy')
tdlat=np.load('tdlat.npy')
tdlon=np.load('tdlon.npy')
tdtimestr=np.load('tdtimestr.npy',allow_pickle=True)
NAME=np.load('name.npy',allow_pickle=True)
tdyear=np.load('tdyear.npy',allow_pickle=True)
LMIR=np.load('LMIR.npy')
"""
#EW
for tt in range(tdew.size):
    if tdew[tt]==1 and tdlon[tt]>153 and tdlon[tt]<158:
        print(tdese[tt],NAME[tt],tdyear[tt],tt)
        print(tdlat[tt],tdlon[tt],tdtimestr[tt])


#MD
for tt in range(tdew.size):
    if tdew[tt]==-2 and tddp[tt]<500:
        print(tddp[tt],NAME[tt],tdyear[tt],tt)
        print(tdlat[tt],tdlon[tt],tdtimestr[tt])
"""

#MD
for tt in range(tdew.size):
    if tdew[tt]==-2 and LMIR[tt]>=27:
        print(LMIR[tt],NAME[tt],tdyear[tt],tt)
        print(tdlat[tt],tdlon[tt],tdtimestr[tt])