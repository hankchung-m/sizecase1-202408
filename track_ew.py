#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 17:06:19 2022

@author: cwuhank
"""


import numpy as np
from datetime import date as dt
from math import*
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import netCDF4 as nc
#from wrf import (getvar,interplevel,ALL_TIMES)
from pandas import to_datetime#, DataFrame

def testdt(datestr):
    try:
        dt.fromisoformat(datestr)
    except:
        return False
    else:
        return True

#!!!begin year:2006 (>=1999) , 2006 RMW more universal
start = 3755
end = 4240 #4240

ncfile = nc.Dataset("IBTrACS/IBTrACS.WP.v04r00.nc")

# 讀取變數
#number = ncfile["number"]
#numobs = ncfile["numobs"]
n_time  = ncfile.dimensions['date_time'].size #360
n_storm = ncfile.dimensions['storm'].size #4236
n_storm = end

#1d
year = ncfile.variables["season"][start:n_storm].data.astype('int')

#2d
status = ncfile.variables["usa_status"][start:n_storm].data.astype('str')
wind = ncfile.variables["usa_wind"][start:n_storm].data.astype('float')
lat = ncfile.variables["usa_lat"][start:n_storm].data.astype('float')
lon = ncfile.variables["usa_lon"][start:n_storm].data.astype('float')
#name = ncfile.variables["name"][start:n_storm].data.astype('str')

tdew=np.load('tdew_m.npy')#_test
name=np.load('name.npy')

typeew4=[1,-1,-1.5,-2]
t1=[0,195]
t2=[195,n_storm-start]

title_num=[['a','b'],['c','d']]


for i in range(2):
    #draw https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.streamplot.html
    fig,ax1=plt.subplots(2,2,sharex=True,sharey=True,figsize=(18,10), subplot_kw={'projection': ccrs.PlateCarree()})

    for ew in range(4): 
        typeew=typeew4[ew]#!!! ew=1, md=-2, mc=-1, ms=-1.5 !!!

        ax1[ew//2,ew%2].gridlines()
        ax1[ew//2,ew%2].coastlines(resolution='50m',linewidth=0.5)
            
        ax1[ew//2,ew%2].set_xticks([100,110,120,130,140,150,160,170,180], crs=ccrs.PlateCarree())
        ax1[ew//2,ew%2].set_yticks([0,5,10,15,20,25,30,35,40], crs=ccrs.PlateCarree())  
        lon_formatter = LongitudeFormatter()   
        lat_formatter = LatitudeFormatter()
        ax1[ew//2,ew%2].xaxis.set_major_formatter(lon_formatter)
        ax1[ew//2,ew%2].yaxis.set_major_formatter(lat_formatter)

        for tt in range(t1[i],t2[i]):#!!! (195)2006-2012, (195,n_storm-start)2013-2021 !!!
            tstime=0
            itime=0
            dp=0
            land=0
            alltm=0
            alltdroci=0
            first=0
            
            drawlon=np.zeros(n_time-8)
            drawlat=np.zeros(n_time-8)
                    
            for tm in range(n_time-8):
                
                if wind[tt,tm]>=5 and lon[tt,tm]>100: #wind>=25kt, wind value, (first wind had already >25kt also)
                    drawlon[tm]=lon[tt,tm]
                    drawlat[tm]=lat[tt,tm]
                else:
                    drawlon[tm]=np.nan
                    drawlat[tm]=np.nan
                
                            
                if wind[tt,tm]>=25 and first==0 : #wind>=25kt, wind value, (first wind had already >25kt also)
                    first=1
                    if tdew[tt] == typeew :
                        if typeew == 1:
                            ax1[ew//2,ew%2].scatter(lon[tt,tm],lat[tt,tm],10,'red', marker='o')
                            print(tt,year[tt],name[tt],'EW')
                        if typeew == -2:
                            ax1[ew//2,ew%2].scatter(lon[tt,tm],lat[tt,tm],10,'blue', marker='o')
                            print(tt,year[tt],name[tt],'MD')
                        if typeew == -1:
                            ax1[ew//2,ew%2].scatter(lon[tt,tm],lat[tt,tm],10,'purple', marker='o')
                            print(tt,year[tt],name[tt],'MC')
                        if typeew == -1.5:
                            ax1[ew//2,ew%2].scatter(lon[tt,tm],lat[tt,tm],10,'green', marker='o')
                            print(tt,year[tt],name[tt],'MS')
                        
            
            if tdew[tt] == typeew :
                if typeew == 1:
                    ax1[ew//2,ew%2].plot(drawlon,drawlat,c='red',linewidth=0.5)
                if typeew == -2:
                    ax1[ew//2,ew%2].plot(drawlon,drawlat,c='blue',linewidth=0.5)
                if typeew == -1:
                    ax1[ew//2,ew%2].plot(drawlon,drawlat,c='purple',linewidth=0.5)
                if typeew == -1.5:
                    ax1[ew//2,ew%2].plot(drawlon,drawlat,c='green',linewidth=0.5)
                    
                        
        #!!!圖例，標題等 
        ax1[ew//2,ew%2].set_xlim(left=100, right=180)
        ax1[ew//2,ew%2].set_ylim(bottom=0, top=40)
        ax1[ew//2,ew%2].set_title('('+title_num[ew//2][ew%2]+')',loc='left', fontsize=25)
        
        if typeew == 1:
            ax1[ew//2,ew%2].set_title('EW', fontsize=25)
 
        if typeew == -2:
            ax1[ew//2,ew%2].set_title('MD', fontsize=25)

        if typeew == -1:
            ax1[ew//2,ew%2].set_title('MC', fontsize=25)

        if typeew == -1.5:
            ax1[ew//2,ew%2].set_title('MS', fontsize=25)

    plt.tight_layout()

    if i ==0:
        plt.savefig('fig/track.png', dpi=600)
    if i ==1:
        plt.savefig('fig/track_2013.png', dpi=600)
    plt.show()