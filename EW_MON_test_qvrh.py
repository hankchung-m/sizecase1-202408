#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 03:02:14 2023

@author: cwuhank
"""


import numpy as np
from datetime import date as dt
#import netCDF4 as nc

def testdt(datestr):
    try:
        dt.fromisoformat(datestr)
    except:
        return False
    else:
        return True

tdqv_5=np.load('tdqv_5.npy')
tdrh_5=np.load('tdrh_5.npy')
tdqv_10=np.load('tdqv_10.npy')
tdrh_10=np.load('tdrh_10.npy')

tdlat_ny=np.load('tdlat_ny.npy')
tdlon_nx=np.load('tdlon_nx.npy')
tdroci=np.load('tdroci.npy')

tdtimestr=np.load('tdtimestr.npy', allow_pickle=True)
tdhour=np.load('tdhour.npy', allow_pickle=True)
name=np.load('name.npy', allow_pickle=True)
LMIR=np.load('LMIR.npy')

tdew=np.zeros(tdtimestr.size)



tdsw=np.zeros(tdtimestr.size)
tdnw=np.zeros(tdtimestr.size)
tdse=np.zeros(tdtimestr.size)
tdese=np.zeros(tdtimestr.size)
tddp=np.zeros(tdtimestr.size)
tdvo=np.zeros(tdtimestr.size)
tdvoe=np.zeros(tdtimestr.size)
tdwc=np.zeros(tdtimestr.size)
tdsh=np.zeros(tdtimestr.size)
tdqv=np.zeros(tdtimestr.size)
tdrh=np.zeros(tdtimestr.size)
for tt in range(tdtimestr.size):
    
    print(tt)
    
    #sort real timestr
    if testdt(tdtimestr[tt]) == False :
        tdsw[tt]=np.nan
        tddp[tt]=np.nan
        tdvo[tt]=np.nan
        tdvoe[tt]=np.nan
        tdsh[tt]=np.nan
        tdqv[tt]=np.nan
        tdrh[tt]=np.nan
        continue
    
    a=tdhour[tt]
    if a != '00' and a != '03' and a != '06' and a != '09' and a != '12' and a != '15' and a != '18' and a != '21':
        tdsw[tt]=np.nan
        tddp[tt]=np.nan
        tdvo[tt]=np.nan
        tdvoe[tt]=np.nan
        tdsh[tt]=np.nan
        tdqv[tt]=np.nan
        tdrh[tt]=np.nan
        continue
    
    tdqv[tt]=(tdqv_10[tt]*4-tdqv_5[tt])/3
    tdrh[tt]=(tdrh_10[tt]*4-tdrh_5[tt])/3
    
if tdqv[tt]==0:
    tdqv[tt]=np.nan
if tdrh[tt]==0:
    tdrh[tt]=np.nan

np.save('tdqv',tdqv)
np.save('tdrh',tdrh)