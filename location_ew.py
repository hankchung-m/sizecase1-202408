#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 16:06:11 2022

@author: cwuhank
"""


import numpy as np
from datetime import date as dt
from math import*
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

def testdt(datestr):
    try:
        dt.fromisoformat(datestr)
    except:
        return False
    else:
        return True

#!!! 0:itime 1:tstime 2:tdtime !!!choose
#draw https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.streamplot.html
fig,ax1=plt.subplots(3,1,sharex=True,sharey=True,figsize=(10,12), subplot_kw={'projection': ccrs.PlateCarree()})
title_num=['a','b','c']
for ii in range(3):
    #ii=0
    LMIR=np.load('LMIR.npy')
    tsrmw=np.load('tsrmw.npy')
    if ii==2:
        tslat=np.load('ilat.npy')
        tslon=np.load('ilon.npy')
    elif ii==0:
        tslat=np.load('tdlat.npy')
        tslon=np.load('tdlon.npy')
    else:
        tslat=np.load('tslat.npy')
        tslon=np.load('tslon.npy')
    tdew=np.load('tdew_m.npy')
    tstimestr=np.load('tstimestr.npy', allow_pickle=True)
    tdmonth=np.load('tdmonth.npy')
    tshour=np.load('tshour.npy', allow_pickle=True)
    
    gen=['EW','MD','MC','MS']
#for i in range(1):
    
    # 製作figure  
    #fig1 = plt.figure()   

    #圖表的設定
    #ax1 = fig1.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax1[ii].gridlines()
    ax1[ii].coastlines(resolution='50m',linewidth=0.5)
    
    ax1[ii].set_xticks([100,110,120,130,140,150,160,170,180], crs=ccrs.PlateCarree())
    ax1[ii].set_yticks([0,5,10,15,20,25,30], crs=ccrs.PlateCarree())  
    lon_formatter = LongitudeFormatter()   
    lat_formatter = LatitudeFormatter()
    ax1[ii].xaxis.set_major_formatter(lon_formatter)
    ax1[ii].yaxis.set_major_formatter(lat_formatter)
    

    n1=0
    n2=0
    n3=0
    n4=0
    ew=0
    md=0
    mc=0
    ms=0
    ewlon=np.zeros(tstimestr.size)
    ewlat=np.zeros(tstimestr.size)
    mdlon=np.zeros(tstimestr.size)
    mdlat=np.zeros(tstimestr.size)
    mclon=np.zeros(tstimestr.size)
    mclat=np.zeros(tstimestr.size)
    mslon=np.zeros(tstimestr.size)
    mslat=np.zeros(tstimestr.size)

    for tt in range(tstimestr.size):
    
        #sort real timestr
        if testdt(tstimestr[tt]) == False :
            
            continue
    
        a=tshour[tt]
        if a != '00' and a != '03' and a != '06' and a != '09' and a != '12' and a != '15' and a != '18' and a != '21':
            
            continue
        
        
        #for selecting cases
        if tdmonth[tt] != 8 and tdmonth[tt] != 9:
            continue
        
        
        if tslon[tt]<0:
            tslon[tt]=360+tslon[tt]
        
        if tdew[tt] == 1 :
            ew+=1
            ewlon[tt]=tslon[tt]
            ewlat[tt]=tslat[tt]
            
            if n1==1 :
                ax1[ii].scatter(tslon[tt],tslat[tt],10,'red',alpha=0.5)
            else:
                ax1[ii].scatter(tslon[tt],tslat[tt],10,'red',alpha=0.5,label=gen[0])
                n1=1
            
        elif tdew[tt] == -2 :
            md+=1
            mdlon[tt]=tslon[tt]
            mdlat[tt]=tslat[tt]
            
            if n2==1 :
                ax1[ii].scatter(tslon[tt],tslat[tt],10,'blue',alpha=0.5)
            else:
                ax1[ii].scatter(tslon[tt],tslat[tt],10,'blue',alpha=0.5,label=gen[1])
                n2=1
                
        elif tdew[tt] == -1 :
            mc+=1
            mclon[tt]=tslon[tt]
            mclat[tt]=tslat[tt]
            
            if n3==1 :
                ax1[ii].scatter(tslon[tt],tslat[tt],10,'purple',alpha=0.5)
            else:
                ax1[ii].scatter(tslon[tt],tslat[tt],10,'purple',alpha=0.5,label=gen[2])
                n3=1
                
        elif tdew[tt] == -1.5 :
            ms+=1
            mslon[tt]=tslon[tt]
            mslat[tt]=tslat[tt]
            
            if n4==1 :
                ax1[ii].scatter(tslon[tt],tslat[tt],10,'green',alpha=0.5)
            else:
                ax1[ii].scatter(tslon[tt],tslat[tt],10,'green',alpha=0.5,label=gen[3])
                n4=1
        
        ewlon[ewlon[tt]==0]=np.nan
        ewlat[ewlat[tt]==0]=np.nan
        mdlon[mdlon[tt]==0]=np.nan
        mdlat[mdlat[tt]==0]=np.nan
        mclon[mclon[tt]==0]=np.nan
        mclat[mclat[tt]==0]=np.nan
        mslon[mslon[tt]==0]=np.nan
        mslat[mslat[tt]==0]=np.nan
    
    ewwlon=np.nanmean(ewlon)
    ewwlat=np.nanmean(ewlat)
    mddlon=np.nanmean(mdlon)
    mddlat=np.nanmean(mdlat)
    mcclon=np.nanmean(mclon)
    mcclat=np.nanmean(mclat)
    msslon=np.nanmean(mslon)
    msslat=np.nanmean(mslat)
    
    ewwlonstd=np.nanstd(ewlon)
    ewwlatstd=np.nanstd(ewlat)
    mddlonstd=np.nanstd(mdlon)
    mddlatstd=np.nanstd(mdlat)
    mcclonstd=np.nanstd(mclon)
    mcclatstd=np.nanstd(mclat)
    msslonstd=np.nanstd(mslon)
    msslatstd=np.nanstd(mslat)
    """
    ewwlon_25=np.nanpercentile(ewlon,25)
    ewwlat_25=np.nanpercentile(ewlat,25)
    mddlon_25=np.nanpercentile(mdlon,25)
    mddlat_25=np.nanpercentile(mdlat,25)
    mcclon_25=np.nanpercentile(mclon,25)
    mcclat_25=np.nanpercentile(mclat,25)
    msslon_25=np.nanpercentile(mslon,25)
    msslat_25=np.nanpercentile(mslat,25)
    
    ewwlon_50=np.nanpercentile(ewlon,50)
    ewwlat_50=np.nanpercentile(ewlat,50)
    mddlon_50=np.nanpercentile(mdlon,50)
    mddlat_50=np.nanpercentile(mdlat,50)
    mcclon_50=np.nanpercentile(mclon,50)
    mcclat_50=np.nanpercentile(mclat,50)
    msslon_50=np.nanpercentile(mslon,50)
    msslat_50=np.nanpercentile(mslat,50)
    
    ewwlon_75=np.nanpercentile(ewlon,75)
    ewwlat_75=np.nanpercentile(ewlat,75)
    mddlon_75=np.nanpercentile(mdlon,75)
    mddlat_75=np.nanpercentile(mdlat,75)
    mcclon_75=np.nanpercentile(mclon,75)
    mcclat_75=np.nanpercentile(mclat,75)
    msslon_75=np.nanpercentile(mslon,75)
    msslat_75=np.nanpercentile(mslat,75)
    """
    
    """
    ax1.plot([ewwlon_25,ewwlon_75],[ewwlat_50,ewwlat_50],color='red')
    ax1.plot([ewwlon_50,ewwlon_50],[ewwlat_25,ewwlat_75],color='red')
    ax1.plot([mddlon_25,mddlon_75],[mddlat_50,mddlat_50],color='blue')
    ax1.plot([mddlon_50,mddlon_50],[mddlat_25,mddlat_75],color='blue')
    ax1.plot([mcclon_25,mcclon_75],[mcclat_50,mcclat_50],color='purple')
    ax1.plot([mcclon_50,mcclon_50],[mcclat_25,mcclat_75],color='purple')
    ax1.plot([msslon_25,msslon_75],[msslat_50,msslat_50],color='green')
    ax1.plot([msslon_50,msslon_50],[msslat_25,msslat_75],color='green')
    """
    ax1[ii].plot([ewwlon-ewwlonstd,ewwlon+ewwlonstd],[ewwlat,ewwlat],color='red')
    ax1[ii].plot([ewwlon,ewwlon],[ewwlat-ewwlatstd,ewwlat+ewwlatstd],color='red')
    ax1[ii].plot([mddlon-mddlonstd,mddlon+mddlonstd],[mddlat,mddlat],color='blue')
    ax1[ii].plot([mddlon,mddlon],[mddlat-mddlatstd,mddlat+mddlatstd],color='blue')
    ax1[ii].plot([mcclon-mcclonstd,mcclon+mcclonstd],[mcclat,mcclat],color='purple')
    ax1[ii].plot([mcclon,mcclon],[mcclat-mcclatstd,mcclat+mcclatstd],color='purple')
    ax1[ii].plot([msslon-msslonstd,msslon+msslonstd],[msslat,msslat],color='green')
    ax1[ii].plot([msslon,msslon],[msslat-msslatstd,msslat+msslatstd],color='green')
    
    ax1[ii].set_xlim(left=100, right=180)
    ax1[ii].set_ylim(bottom=0, top=30)
    
    #ax1[ii].set_title('location')
    ax1[ii].set_title('('+title_num[ii]+')',loc='left', fontsize=25)
        
    if ii==2:
        ax1[ii].set_title('I_TIME', fontsize=25)
        #plt.savefig('fig/location_I_EW.png', dpi=300)
    elif ii==0:
        ax1[ii].set_title('TD_TIME', fontsize=25)
        #plt.savefig('fig/location_TD_EW.png', dpi=300)
    else:
        ax1[ii].set_title('TS_TIME', fontsize=25)
        #plt.savefig('fig/location_EW.png', dpi=300)

#圖例，標題等
ax1[0].text(102,16,str(ew),color='red', fontsize=15)
ax1[0].text(102,19,str(md),color='blue', fontsize=15)
ax1[0].text(102,25,str(mc),color='purple', fontsize=15)
ax1[0].text(102,22,str(ms),color='green', fontsize=15)
ax1[0].legend(loc='upper right',fontsize=15)#
plt.tight_layout()
plt.savefig('fig/location_EW.png', dpi=600)
plt.show()
    #plt.title('location, level='+str(lev))
    #plt.savefig('fig/location_'+str(lev)+'.png')
    #plt.show()
    
    #ax[0].title('location, level='+str(lev))
    #plt.legend(loc='upper right')
    #ax[0].savefig('fig/location_'+str(lev)+'.png')
    #ax[0].show()