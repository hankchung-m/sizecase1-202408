#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  4 15:35:45 2022

@author: cwuhank
"""


import numpy as np
from datetime import date as dt
import netCDF4 as nc
from math import*

def Distance(lat1,lng1,lat2,lng2): #https://www.itread01.com/content/1550369361.html
    radlat1=radians(lat1)  
    radlat2=radians(lat2)  
    a=radlat1-radlat2  
    b=radians(lng1)-radians(lng2)  
    s=2*asin(sqrt(pow(sin(a/2),2)+cos(radlat1)*cos(radlat2)*pow(sin(b/2),2)))  
    earth_radius=6378.137  
    s=s*earth_radius  
    if s<0:  
        return -s  
    else:  
        return s

def testdt(datestr):
    try:
        dt.fromisoformat(datestr)
    except:
        return False
    else:
        return True

ilat=np.load('ilat.npy')
ilon=np.load('ilon.npy')
tslat=np.load('tslat.npy')
tslon=np.load('tslon.npy')

ilat_ny=np.load('ilat_ny.npy')
ilon_nx=np.load('ilon_nx.npy')
tslat_ny=np.load('tslat_ny.npy')
tslon_nx=np.load('tslon_nx.npy')

itimestr=np.load('itimestr.npy', allow_pickle=True)
ihour=np.load('ihour.npy', allow_pickle=True)
tstimestr=np.load('tstimestr.npy', allow_pickle=True)
tshour=np.load('tshour.npy', allow_pickle=True)

#!!! 0:itime 1:tstime !!!choose
ii=0
if ii==1:
    itimestr=tstimestr
    ihour=tshour

interxy=8#8degree(+-4degree)

ivws=np.zeros(itimestr.size)
irh=np.zeros(itimestr.size)
for tt in range(itimestr.size):
    
    #print(tt)
    
    #sort real timestr
    if testdt(itimestr[tt]) == False :
        continue
    
    a=ihour[tt]
    if a != '00' and a != '03' and a != '06' and a != '09' and a != '12' and a != '15' and a != '18' and a != '21':
        continue
    
    if tt/10<1 :
        ncfile = nc.Dataset("ERA5_I00%d.nc"%(tt))
    elif tt/100<1 :
        ncfile = nc.Dataset("ERA5_I0%d.nc"%(tt))
    else:
        ncfile = nc.Dataset("ERA5_I%d.nc"%(tt))
        
    XLAT = ncfile["latitude"]
    XLONG = ncfile["longitude"]
    p = ncfile["level"]
    u = ncfile["u"]
    v = ncfile["v"]
    z = ncfile["z"]
    r = ncfile["r"]
    #z=z/9.8
    #print(u.shape)
    
    
    
    # Domain x and y position
    x2=ilon_nx[tt]
    y2=tslat_ny[tt]
    
    #calculate
    point=0
    ushearall=0
    vshearall=0
    rhall=0#20220328
    for i in range(XLONG.size) :
        for j in range(XLAT.size) :
            if Distance(XLAT[y2],XLONG[x2],XLAT[j],XLONG[i])>200 and Distance(XLAT[y2],XLONG[x2],XLAT[j],XLONG[i])<600 :
                point+=1
                ushear=u[0,14,j,i]-u[0,30,j,i]
                vshear=v[0,14,j,i]-v[0,30,j,i]
                ushearall+=ushear
                vshearall+=vshear
                rh=np.average(r[0,17:22,j,i])
                rhall+=rh
                #print(XLAT[j])
                #print(XLONG[i])
                
    if point!=0 :
        ushearavg=ushearall/point
        vshearavg=vshearall/point
        ivws[tt]=sqrt(ushearavg**2+vshearavg**2)
        rhavg=rhall/point
        irh[tt]=rhavg
    else:
        ivws[tt]=np.nan
        irh[tt]=np.nan
    #print(tt)
    #print(itimestr[tt])
    #print(ilat[tt])
    #print(ilon[tt])
    #print(ivws[tt])
    #print(irh[tt])
    
    
np.save('ivws',ivws)
np.save('irh',irh)
    
