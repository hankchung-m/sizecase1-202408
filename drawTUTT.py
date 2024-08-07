#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 20:14:53 2022

@author: cwuhank
"""


import numpy as np
from datetime import date as dt
from math import*
import matplotlib.pyplot as plt
import netCDF4 as nc
import cartopy.crs as ccrs
from composite import composite
from vorticity import vorticity
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

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

def interp(f1,f2,f3,f4,XLONGi,XLONGi1,nlon,XLATj,XLATj1,nlat):#https://blog.csdn.net/qq_45039390/article/details/105633440
    a1=f1*(XLONGi1-nlon)*(XLATj1-nlat)
    a2=f2*(nlon-XLONGi)*(XLATj1-nlat)
    a3=f3*(XLONGi1-nlon)*(nlat-XLATj)
    a4=f4*(nlon-XLONGi)*(nlat-XLATj)
    a=a1+a2+a3+a4
    return a


LMIR=np.load('LMIR.npy')

tdlat=np.load('tdlat.npy')
tdlon=np.load('tdlon.npy')

tstimestr=np.load('tstimestr.npy', allow_pickle=True)
tshour=np.load('tshour.npy', allow_pickle=True)
tdtimestr=np.load('tdtimestr.npy', allow_pickle=True)
tdhour=np.load('tdhour.npy', allow_pickle=True)
tdew=np.load('tdew.npy')
tdutcl=np.load('tdutcl.npy')
tdutcllat=np.load('tdutcllat.npy')
tdutcllon=np.load('tdutcllon.npy')

relat=np.zeros(tdtimestr.size)
relon=np.zeros(tdtimestr.size)

# 製作figure  
fig1 = plt.figure()   

#圖表的設定

"""
ax1 = fig1.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax1.gridlines()
ax1.set_xticks([-2000,-1500,-1000,-500,0,500,1000,1500,2000], crs=ccrs.PlateCarree())
ax1.set_yticks([-2000,-1500,-1000,-500,0,500,1000,1500,2000], crs=ccrs.PlateCarree())  
"""

ax1 = fig1.add_subplot(1, 1, 1, aspect='equal')
ax1.set_xticks([-2000,-1500,-1000,-500,0,500,1000,1500,2000])
ax1.set_yticks([-2000,-1500,-1000,-500,0,500,1000,1500,2000])  
ax1.tick_params(labelsize=8)
ax1.grid(True)

#lon_formatter = LongitudeFormatter()   
#lat_formatter = LatitudeFormatter()
#ax1.xaxis.set_major_formatter(lon_formatter)
#ax1.yaxis.set_major_formatter(lat_formatter)

for tt in range(tdtimestr.size):
    
    #sort real timestr
    if testdt(tdtimestr[tt]) == False :
        continue
    
    a=tdhour[tt]
    if a != '00' and a != '03' and a != '06' and a != '09' and a != '12' and a != '15' and a != '18' and a != '21':
        continue
    
    relat[tt]=tdlat[tt]-tdutcllat[tt]
    relon[tt]=tdlon[tt]-tdutcllon[tt]
    
    redx=Distance(tdutcllat[tt],tdlon[tt],tdutcllat[tt],tdutcllon[tt])
    redy=Distance(tdlat[tt],tdlon[tt],tdutcllat[tt],tdlon[tt])
    if tdlon[tt]<tdutcllon[tt]:
        redx=-redx
    if tdlat[tt]<tdutcllat[tt]:
        redy=-redy
    
    if tdew[tt]==1:
        ax1.scatter(redx,redy,5,'red')
    elif tdew[tt]==-1:
        ax1.scatter(redx,redy,5,'purple')
    elif tdew[tt]==-2:
        ax1.scatter(redx,redy,5,'blue')
    elif tdew[tt]==-1.5:
        ax1.scatter(redx,redy,5,'green')
        
nthe = 360                  # angle for azimuthal mean
inv_the = 0.1                 # angle resolution for azimuthal mean 
for the in range(0,int(nthe/inv_the)):
    rthe=the*inv_the
    x=1700*cos(radians(rthe))
    y=1700*sin(radians(rthe))
    ax1.scatter(x,y,0.2,'black')
plt.scatter(0, 0, marker='x', color='black')
ax1.set_xlim(left=-2000, right=2000)
ax1.set_ylim(bottom=-2000, top=2000)

ax1.set_xlabel('(km)', fontsize=8)
ax1.set_ylabel('(km)', fontsize=8)

plt.savefig('fig/utcl.png', dpi=300)