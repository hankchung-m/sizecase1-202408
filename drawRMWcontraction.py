#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 09:19:22 2022

@author: cwuhank
"""


import matplotlib.pyplot as plt
from datetime import date as dt
import numpy as np
#import cartopy.crs as ccrs
import netCDF4 as nc
#from wrf import (getvar,interplevel,ALL_TIMES)
from pandas import to_datetime#, DataFrame
from sklearn.metrics import r2_score
from matplotlib.patches import Ellipse

#!!!begin year:2006 (>=1999) , 2006 RMW more universal
start = 3755
end = 4240 #4240

#thanks to Ken Shiu
def get_str(var, number, time = True):
    return np.reshape(np.array(list(map(''.join, zip(*[iter(var.flatten())]*number)))), (n_storm-start, n_time) if time == True else (n_storm-start))

def reg(x,y):
    num=np.polyfit(x, y, 1)#np.log(x) nomial.polynomial.Polynomial.
    eq=np.poly1d(num)
    cor=r2_score(y,eq(x))
    return num,cor

def testdt(datestr):
    try:
        dt.fromisoformat(datestr)
    except:
        return False
    else:
        return True

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
wind = ncfile.variables["usa_wind"][start:n_storm].data.astype('float')
lat = ncfile.variables["usa_lat"][start:n_storm].data.astype('float')
lon = ncfile.variables["usa_lon"][start:n_storm].data.astype('float')
rmw = ncfile.variables["usa_rmw"][start:n_storm].data.astype('float')
rmw = rmw*1.852 #nmile to km https://www.hackmath.net/en/calculator/conversion-of-length-units?unit1=n+mile&unit2=km&dir=1
landfall = ncfile.variables["landfall"][start:n_storm].data.astype('float')
dist2land = ncfile.variables["dist2land"][start:n_storm].data.astype('float')

name = ncfile.variables["name"][start:n_storm].data.astype('str')
name = get_str(name, 128, time = False)
#print(landfall[1])

#3d
wpid = ncfile.variables["usa_atcf_id"][start:n_storm].data.astype('str')

time = ncfile.variables["iso_time"][start:n_storm].data.astype('str')
time = get_str(time, 19)
hour = np.reshape(to_datetime(time.flatten(), format = '%Y-%m-%d %H:%M:%S').strftime('%H'), (n_storm-start, n_time))
timestr = np.reshape(to_datetime(time.flatten(), format = '%Y-%m-%d %H:%M:%S').strftime('%Y-%m-%d'), (n_storm-start, n_time))
#print(timestr)
#print(hour)

t = np.arange(start,n_storm)
#plt.plot(t,year,"r.")



ivws=np.load('ivws.npy')




#Lifetime Max. Intensification Rate (LMIR), Max. Intensification rmw (irmw), TS rmw (tsrmw), TS to Max. intensification(ts2i)
LMIR=np.zeros((n_storm-start))

tsrmw=np.zeros((n_storm-start))
tdrmw=np.zeros((n_storm-start))
ts2i=np.ones((n_storm-start))
#TD time (tdtimestr,tdhour), TD lat & lon (tdlat,tdlon) for ERA5 (composite)
tdtimestr=np.empty_like(timestr,shape=n_storm-start,dtype=object)
tdhour=np.empty_like(hour,shape=n_storm-start,dtype=object)
tdlat=np.zeros((n_storm-start))
tdlon=np.zeros((n_storm-start))
#TS time (tstimestr,tshour), TS lat & lon (tslat,tslon) for ERA5 (composite)
tstimestr=np.empty_like(timestr,shape=n_storm-start,dtype=object)
tshour=np.empty_like(hour,shape=n_storm-start,dtype=object)
tslat=np.zeros((n_storm-start))
tslon=np.zeros((n_storm-start))
#Max. Intensification time (itimestr,ihour), Max. Intensification lat & lon (ilat,ilon) for ERA5 (shear, SST, OHC)
itimestr=np.empty_like(timestr,shape=n_storm-start,dtype=object)
ihour=np.empty_like(hour,shape=n_storm-start,dtype=object)
ilat=np.zeros((n_storm-start))
ilon=np.zeros((n_storm-start))

#intensity
wind_t24=np.zeros((n_storm-start))
wind_0=np.zeros((n_storm-start))
wind_24=np.zeros((n_storm-start))
wind_48=np.zeros((n_storm-start))
wind_72=np.zeros((n_storm-start))
#RMW
irmw_t24=np.zeros((n_storm-start))
irmw=np.zeros((n_storm-start))
irmw_24=np.zeros((n_storm-start))
irmw_48=np.zeros((n_storm-start))
irmw_72=np.zeros((n_storm-start))
#24
itimestr_24=np.empty_like(timestr,shape=n_storm-start,dtype=object)
ihour_24=np.empty_like(hour,shape=n_storm-start,dtype=object)
ilat_24=np.zeros((n_storm-start))
ilon_24=np.zeros((n_storm-start))

irmw_c=np.zeros((n_storm-start))

for tt in range(n_storm-start):
    tstime=0
    itime=0
    dp=0
    land=0
    for tm in range(n_time-8):
        
        if wind[tt,tm]>=25 : #wind>=25kt, wind value, (first wind had already >25kt also)
            if tdrmw[tt]==0 : #first ts
                tdrmw[tt]=rmw[tt,tm]
                tdtimestr[tt]=timestr[tt,tm]
                tdhour[tt]=hour[tt,tm]
                tdlat[tt]=lat[tt,tm]
                tdlon[tt]=lon[tt,tm]
                
            if  dist2land[tt,tm]<200 :   #first ts tsrmw[tt]!=0 and
                land=2
            if  landfall[tt,tm]==0 :#after landfall
                land=1
                    
        
        if wind[tt,tm]>34 and rmw[tt,tm]>0 : #wind>34kt, wind value, rmw value
            if tsrmw[tt]==0 : #first ts
                tsrmw[tt]=rmw[tt,tm]
                tstimestr[tt]=timestr[tt,tm]
                tshour[tt]=hour[tt,tm]
                tslat[tt]=lat[tt,tm]
                tslon[tt]=lon[tt,tm]
                
                tstime=tm
                    
                
        
        a=hour[tt,tm+8]
        b=hour[tt,tm]
        if a=='00'or a=='03'or a=='06'or a=='09'or a=='12'or a=='15'or a=='18'or a=='21' : #value
            if a==b and wind[tt,tm+8]>34 and wind[tt,tm]>0 and rmw[tt,tm]>0 : #intval 24hour, wind>34kt, wind value, rmw value
                dp=wind[tt,tm+8]-wind[tt,tm] #Intensification Rate
                if dp>LMIR[tt] : #sort Lifetime Max. Intensification Rate
                    LMIR[tt]=dp
                    
                    wind_t24[tt]=wind[tt,tm+8]
                    wind_0[tt]=wind[tt,tm]
                    if wind[tt,tm-8]>0:
                        wind_24[tt]=wind[tt,tm-8]
                    if wind[tt,tm-16]>0:
                        wind_48[tt]=wind[tt,tm-16]
                    if wind[tt,tm-24]>0:
                        wind_72[tt]=wind[tt,tm-24]
                    
                    if rmw[tt,tm+8]>0:
                        irmw_t24[tt]=rmw[tt,tm+8]
                    if rmw[tt,tm]>0:
                        irmw[tt]=rmw[tt,tm]
                    if rmw[tt,tm-8]>0:
                        irmw_24[tt]=rmw[tt,tm-8]
                    if rmw[tt,tm-16]>0:
                        irmw_48[tt]=rmw[tt,tm-16]
                    if rmw[tt,tm-24]>0:
                        irmw_72[tt]=rmw[tt,tm-24]
                        
                    itimestr_24[tt]=timestr[tt,tm-8]
                    ihour_24[tt]=hour[tt,tm-8]
                    if lat[tt,tm-8]>-180:
                        ilat_24[tt]=lat[tt,tm-8]
                    if lon[tt,tm-8]>-180:
                        ilon_24[tt]=lon[tt,tm-8]
                    
                    itimestr[tt]=timestr[tt,tm] #+0 =start of intensifying, +8 =end of intensifying
                    ihour[tt]=hour[tt,tm]
                    ilat[tt]=lat[tt,tm]
                    ilon[tt]=lon[tt,tm]
                    
                    itime=tm
                    
                    
                    if landfall[tt,tm+8]==0 or land==1 or land==2:
                        LMIR[tt]=np.nan
                        
                        wind_t24[tt]=np.nan
                        wind_0[tt]=np.nan
                        wind_24[tt]=np.nan
                        wind_48[tt]=np.nan
                        wind_72[tt]=np.nan
                        irmw_t24[tt]=np.nan
                        irmw[tt]=np.nan
                        irmw_24[tt]=np.nan
                        irmw_48[tt]=np.nan
                        irmw_72[tt]=np.nan
                        itimestr_24[tt]=np.nan
                        ihour_24[tt]=np.nan
                        ilat_24[tt]=np.nan
                        ilon_24[tt]=np.nan
                        
                        itimestr[tt]=np.nan
                        ihour[tt]=np.nan
                        ilat[tt]=np.nan
                        ilon[tt]=np.nan
                        
                        #TS composite should exclude unfavorable factors
                        tsrmw[tt]=np.nan
                        tstimestr[tt]=np.nan
                        tshour[tt]=np.nan
                        tslat[tt]=np.nan
                        tslon[tt]=np.nan
                        #TD composite should exclude unfavorable factors
                        tdrmw[tt]=np.nan
                        tdtimestr[tt]=np.nan
                        tdhour[tt]=np.nan
                        tdlat[tt]=np.nan
                        tdlon[tt]=np.nan
                        
    if tslat[tt]>26 or tslon[tt]<100 or tslon[tt]<121 or ivws[tt]>8 :
        LMIR[tt]=np.nan
        
        wind_t24[tt]=np.nan
        wind_0[tt]=np.nan
        wind_24[tt]=np.nan
        wind_48[tt]=np.nan
        wind_72[tt]=np.nan
        irmw_t24[tt]=np.nan
        irmw[tt]=np.nan
        irmw_24[tt]=np.nan
        irmw_48[tt]=np.nan
        irmw_72[tt]=np.nan
        itimestr_24[tt]=np.nan
        ihour_24[tt]=np.nan
        ilat_24[tt]=np.nan
        ilon_24[tt]=np.nan
        
        itimestr[tt]=np.nan
        ihour[tt]=np.nan
        ilat[tt]=np.nan
        ilon[tt]=np.nan
        
        
        #TS composite should exclude unfavorable factors
        tsrmw[tt]=np.nan
        tstimestr[tt]=np.nan
        tshour[tt]=np.nan
        tslat[tt]=np.nan
        tslon[tt]=np.nan
        #TD composite should exclude unfavorable factors
        tdrmw[tt]=np.nan
        tdtimestr[tt]=np.nan
        tdhour[tt]=np.nan
        tdlat[tt]=np.nan
        tdlon[tt]=np.nan
                    
                    
                    
    ts2i[tt]=(itime-tstime)*3
    

        
    #exclude no data
    if LMIR[tt]==0 or tsrmw[tt]==0 or irmw[tt]==0 or ts2i[tt]==1 or np.isnan(LMIR[tt])==True or ilon_24[tt]==0:
        LMIR[tt]=np.nan
        #print(name[tt],year[tt],time[tt,0],max(wind[tt]))
        
        wind_t24[tt]=np.nan
        wind_0[tt]=np.nan
        wind_24[tt]=np.nan
        wind_48[tt]=np.nan
        wind_72[tt]=np.nan
        irmw_t24[tt]=np.nan
        irmw[tt]=np.nan
        irmw_24[tt]=np.nan
        irmw_48[tt]=np.nan
        irmw_72[tt]=np.nan
        
        ts2i[tt]=np.nan
        
        tsrmw[tt]=np.nan
        tstimestr[tt]=np.nan
        tshour[tt]=np.nan
        tslat[tt]=np.nan
        tslon[tt]=np.nan
        
        tdrmw[tt]=np.nan
        tdtimestr[tt]=np.nan
        tdhour[tt]=np.nan
        tdlat[tt]=np.nan
        tdlon[tt]=np.nan
        
        itimestr_24[tt]=np.nan
        ihour_24[tt]=np.nan
        ilat_24[tt]=np.nan
        ilon_24[tt]=np.nan
        
    irmw_c[tt]=irmw_24[tt]-irmw[tt]

np.save('itimestr_24',itimestr_24)
np.save('ihour_24',ihour_24)
np.save('ilat_24',ilat_24)
np.save('ilon_24',ilon_24)
np.save('irmw_c',irmw_c)

np.save('wind_t24', wind_t24)
np.save('wind_0', wind_0)
np.save('wind_24', wind_24)
np.save('wind_48', wind_48)
np.save('wind_72', wind_72)

np.save('irmw_t24', irmw_t24)
np.save('irmw', irmw)
np.save('irmw_24', irmw_24)
np.save('irmw_48', irmw_48)
np.save('irmw_72', irmw_72)


hhll=['h','l']
rmwind=['rmw_','wind_']
itime_5=['t24','0','24','48','72']
color_5=['blue','purple','red','orange','yellow']

# 製作figure  
for lev in range(30,80,5):#30-75
    fig,ax=plt.subplots(1,2,sharex=False,sharey=True,figsize=(18,8))
    
    n1=0
    n2=0
    h=0
    l=0
    hrmw_t24=np.zeros(tstimestr.size)
    hrmw_0=np.zeros(tstimestr.size)
    hrmw_24=np.zeros(tstimestr.size)
    hrmw_48=np.zeros(tstimestr.size)
    hrmw_72=np.zeros(tstimestr.size)
    hwind_t24=np.zeros(tstimestr.size)
    hwind_0=np.zeros(tstimestr.size)
    hwind_24=np.zeros(tstimestr.size)
    hwind_48=np.zeros(tstimestr.size)
    hwind_72=np.zeros(tstimestr.size)
    lrmw_t24=np.zeros(tstimestr.size)
    lrmw_0=np.zeros(tstimestr.size)
    lrmw_24=np.zeros(tstimestr.size)
    lrmw_48=np.zeros(tstimestr.size)
    lrmw_72=np.zeros(tstimestr.size)
    lwind_t24=np.zeros(tstimestr.size)
    lwind_0=np.zeros(tstimestr.size)
    lwind_24=np.zeros(tstimestr.size)
    lwind_48=np.zeros(tstimestr.size)
    lwind_72=np.zeros(tstimestr.size)
    
    
    levi=int((lev-30)/5)
    for tt in range(tstimestr.size):
        #sort real timestr
        if testdt(tstimestr[tt]) == False :
            
            continue
    
        a=tshour[tt]
        if a != '00' and a != '03' and a != '06' and a != '09' and a != '12' and a != '15' and a != '18' and a != '21':
            
            continue
        
        
        
        if LMIR[tt] >= lev :
            
            h+=1
            hrmw_t24[tt]=irmw_t24[tt]
            hrmw_0[tt]=irmw[tt]
            hrmw_24[tt]=irmw_24[tt]
            hrmw_48[tt]=irmw_48[tt]
            hrmw_72[tt]=irmw_72[tt]
            hwind_t24[tt]=wind_t24[tt]
            hwind_0[tt]=wind_0[tt]
            hwind_24[tt]=wind_24[tt]
            hwind_48[tt]=wind_48[tt]
            hwind_72[tt]=wind_72[tt]
            

            for j in range(len(rmwind)):
                    for k in range(len(itime_5)):
                        exec(hhll[0]+rmwind[j]+itime_5[k]+'['+hhll[0]+rmwind[j]+itime_5[k]+'==0]=np.nan')
                        #exec(hhll[i]+rmwind[j]+itime_5[k]+'['+hhll[i]+rmwind[j]+itime_5[k]+'==-9999]=np.nan')
                        #exec(hhll[i]+rmwind[j]+itime_5[k]+'['+hhll[i]+rmwind[j]+itime_5[k]+'==-18518.1]=np.nan')
        
            
            ax[0].scatter(hrmw_t24[tt],hwind_t24[tt],35,'blue', alpha=0.5, linewidths=0)
            ax[0].scatter(hrmw_0[tt],hwind_0[tt],35,'purple', alpha=0.5, linewidths=0)
            ax[0].scatter(hrmw_24[tt],hwind_24[tt],35,'red', alpha=0.5, linewidths=0)
            ax[0].scatter(hrmw_48[tt],hwind_48[tt],35,'orange', alpha=0.5, linewidths=0)
            ax[0].scatter(hrmw_72[tt],hwind_72[tt],35,'yellow', alpha=0.5, linewidths=0)
            
            
            
        elif LMIR[tt] < lev :
            
            l+=1
            lrmw_t24[tt]=irmw_t24[tt]
            lrmw_0[tt]=irmw[tt]
            lrmw_24[tt]=irmw_24[tt]
            lrmw_48[tt]=irmw_48[tt]
            lrmw_72[tt]=irmw_72[tt]
            lwind_t24[tt]=wind_t24[tt]
            lwind_0[tt]=wind_0[tt]
            lwind_24[tt]=wind_24[tt]
            lwind_48[tt]=wind_48[tt]
            lwind_72[tt]=wind_72[tt]
            

            for j in range(len(rmwind)):
                    for k in range(len(itime_5)):
                        exec(hhll[1]+rmwind[j]+itime_5[k]+'['+hhll[1]+rmwind[j]+itime_5[k]+'==0]=np.nan')
        
            
            if n2==1 :
                ax[1].scatter(lrmw_t24[tt],lwind_t24[tt],35,'blue', alpha=0.5, linewidths=0)
                ax[1].scatter(lrmw_0[tt],lwind_0[tt],35,'purple', alpha=0.5, linewidths=0)
                ax[1].scatter(lrmw_24[tt],lwind_24[tt],35,'red', alpha=0.5, linewidths=0)
                ax[1].scatter(lrmw_48[tt],lwind_48[tt],35,'orange', alpha=0.5, linewidths=0)
                ax[1].scatter(lrmw_72[tt],lwind_72[tt],35,'yellow', alpha=0.5, linewidths=0)
            else:
                ax[1].scatter(lrmw_t24[tt],lwind_t24[tt],35,'blue', alpha=0.5, linewidths=0,label='+24h')
                ax[1].scatter(lrmw_0[tt],lwind_0[tt],35,'purple', alpha=0.5, linewidths=0,label='-0h')
                ax[1].scatter(lrmw_24[tt],lwind_24[tt],35,'red', alpha=0.5, linewidths=0,label='-24h')
                ax[1].scatter(lrmw_48[tt],lwind_48[tt],35,'orange', alpha=0.5, linewidths=0,label='-48h')
                ax[1].scatter(lrmw_72[tt],lwind_72[tt],35,'yellow', alpha=0.5, linewidths=0,label='-72h')
                n2=1
            
    
    for i in range(len(hhll)):
            for j in range(len(rmwind)):
                for k in range(len(itime_5)):
                    exec(hhll[i]+hhll[i]+rmwind[j]+itime_5[k]+'=np.nanmean('+hhll[i]+rmwind[j]+itime_5[k]+')')
    
    for i in range(len(hhll)):
            for j in range(len(rmwind)):
                for k in range(len(itime_5)):
                    exec(hhll[i]+hhll[i]+rmwind[j]+itime_5[k]+'_std=np.nanstd('+hhll[i]+rmwind[j]+itime_5[k]+')*2')
    
    for i in range(len(hhll)):
        for k in range(len(itime_5)):
            exec('ax['+str(i)+'].add_patch(Ellipse(xy = ('+hhll[i]+hhll[i]+rmwind[0]+itime_5[k]+','+hhll[i]+hhll[i]+rmwind[1]+itime_5[k]+
                     '), width = '+hhll[i]+hhll[i]+rmwind[0]+itime_5[k]+'_std'+
                     ', height = '+hhll[i]+hhll[i]+rmwind[1]+itime_5[k]+'_std'+',facecolor="'+color_5[k]+'",alpha=0.1,linewidth=5,edgecolor="'+color_5[k]+'"))')
        exec(hhll[i]+'=str(np.around('+hhll[i]+hhll[i]+rmwind[0]+itime_5[2]+'-'+hhll[i]+hhll[i]+rmwind[0]+itime_5[1]+',decimals=4))')
        ax[i].text(100,150,'RMW contraction(km)=', fontsize=15)
        exec('ax['+str(i)+'].text(100,140,'+hhll[i]+', fontsize=15)')    
        
        
        ax[i].grid(True)
        ax[i].grid(True)
        ax[i].set_xlim(left=0, right=300)
        ax[i].set_ylim(bottom=-0, top=160)
        ax[i].set_xlabel('RMW (km)', fontsize=15)
        ax[i].set_ylabel('Intensity (kt)', fontsize=15)
        ax[i].tick_params(labelsize=15)
        ax[1].legend(prop={'size': 15})
        if i==0:
            ax[i].set_title('RMW vs. Intensity, >='+str(lev)+'kt/24h', fontsize=15)
        else:
            ax[i].set_title('RMW vs. Intensity, <'+str(lev)+'kt/24h', fontsize=15)
        
        
    fig.savefig('fig/RMWcontraction_'+str(lev)+'.png', dpi=300)
    