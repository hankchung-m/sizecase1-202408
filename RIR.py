#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 15:12:50 2021

@author: cwuhank
"""


import matplotlib.pyplot as plt
import numpy as np
#import cartopy.crs as ccrs
import netCDF4 as nc
#from wrf import (getvar,interplevel,ALL_TIMES)
from pandas import to_datetime#, DataFrame
from sklearn.metrics import r2_score
from scipy.optimize import curve_fit
from scipy.stats import pearsonr
from scipy.stats import spearmanr

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

def func(x,a,b):
    return a/x+b

#!!!do you want to calculus VWS, RH, UTCL, or EW (means ALL cases, download ERA5)? yes:1 no:0
forVWS=0

ncfile = nc.Dataset("IBTrACS/IBTrACS.WP.v04r00.nc")

# 讀取變數
#number = ncfile["number"]
#numobs = ncfile["numobs"]
n_time  = ncfile.dimensions['date_time'].size #360
n_storm = ncfile.dimensions['storm'].size #4236 => 4253
#print(n_storm)
n_storm = end


#1d
year = ncfile.variables["season"][start:n_storm].data.astype('int')

#t = np.arange(4235,4238)
#yeart=ncfile.variables["season"][4235:4238].data.astype('int')
#plt.plot(t,yeart,"r.")

#2d
status = ncfile.variables["usa_status"][start:n_storm].data.astype('str')
wind = ncfile.variables["usa_wind"][start:n_storm].data.astype('float')
lat = ncfile.variables["usa_lat"][start:n_storm].data.astype('float')
lon = ncfile.variables["usa_lon"][start:n_storm].data.astype('float')
rmw = ncfile.variables["usa_rmw"][start:n_storm].data.astype('float')
rmw = rmw*1.852 #nmile to km https://www.hackmath.net/en/calculator/conversion-of-length-units?unit1=n+mile&unit2=km&dir=1
roci = ncfile.variables["usa_roci"][start:n_storm].data.astype('float')
roci = roci*1.852 #nmile to km https://www.hackmath.net/en/calculator/conversion-of-length-units?unit1=n+mile&unit2=km&dir=1
landfall = ncfile.variables["landfall"][start:n_storm].data.astype('float')
dist2land = ncfile.variables["dist2land"][start:n_storm].data.astype('float')

name = ncfile.variables["name"][start:n_storm].data.astype('str')
name = get_str(name, 128, time = False)
#print(landfall[1])
#print(name[74],wind[74])

#3d
#wpid = ncfile.variables["usa_atcf_id"][start:n_storm].data.astype('str')

time = ncfile.variables["iso_time"][start:n_storm].data.astype('str')
time = get_str(time, 19)
hour = np.reshape(to_datetime(time.flatten(), format = '%Y-%m-%d %H:%M:%S').strftime('%H'), (n_storm-start, n_time))
month = np.reshape(to_datetime(time.flatten(), format = '%Y-%m-%d %H:%M:%S').strftime('%m'), (n_storm-start, n_time))
year = np.reshape(to_datetime(time.flatten(), format = '%Y-%m-%d %H:%M:%S').strftime('%Y'), (n_storm-start, n_time))
timestr = np.reshape(to_datetime(time.flatten(), format = '%Y-%m-%d %H:%M:%S').strftime('%Y-%m-%d'), (n_storm-start, n_time))
ttimestr = np.reshape(to_datetime(time.flatten(), format = '%Y-%m-%d %H:%M:%S').strftime('%Y%m%d%H'), (n_storm-start, n_time))

#print(timestr)
#print(hour)




if forVWS==0:
    ivws=np.load('ivws.npy')
    irh=np.load('irh.npy')

else:
    ivws=np.zeros((n_storm-start))
    irh=np.zeros((n_storm-start))

#!!!
tdew=np.load('tdew_m.npy')

#Lifetime Max. Intensification Rate (LMIR), Max. Intensification rmw (irmw), TS rmw (tsrmw), TS to Max. intensification(ts2i)
LMIR=np.zeros((n_storm-start))
#irmw=np.zeros((n_storm-start))
tsrmw=np.zeros((n_storm-start))
tdrmw=np.zeros((n_storm-start))
tdroci=np.zeros((n_storm-start))
ts2i=np.ones((n_storm-start))
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

#TD time (tdtimestr,tdhour), TD lat & lon (tdlat,tdlon) for ERA5 (composite)
tdtimestr=np.empty_like(timestr,shape=n_storm-start,dtype=object)
tdtimestr_sat=np.empty_like(timestr,shape=n_storm-start,dtype=object)
tdhour=np.empty_like(hour,shape=n_storm-start,dtype=object)
tdmonth=np.zeros((n_storm-start))
tdyear=np.empty_like(hour,shape=n_storm-start,dtype=object)
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
#TUTT UTCL
tdtttimestr=np.empty_like(ttimestr,shape=n_storm-start,dtype=object)

#enviromental factors
env_lf_i=np.zeros((n_storm-start))
env_lf_s=np.zeros((n_storm-start))
env_dl=np.zeros((n_storm-start))
env_tslat=np.zeros((n_storm-start))
env_tslon=np.zeros((n_storm-start))
env_ivws=np.zeros((n_storm-start))

#shear
l_LMIR=np.zeros((n_storm-start))
l_irmw=np.zeros((n_storm-start))
l_tsrmw=np.zeros((n_storm-start))
m_LMIR=np.zeros((n_storm-start))
m_irmw=np.zeros((n_storm-start))
m_tsrmw=np.zeros((n_storm-start))
h_LMIR=np.zeros((n_storm-start))
h_irmw=np.zeros((n_storm-start))
h_tsrmw=np.zeros((n_storm-start))

#print(status.shape)
#print(status[275,0,0])
#print(wind[474,:])

for tt in range(n_storm-start):
    tstime=0
    itime=0
    dp=0
    land=0
    alltm=0
    alltdroci=0
    tdts=0
    for tm in range(n_time-8):
        if roci[tt,tm]>=500:
            alltdroci+=roci[tt,tm]
            alltm+=1
        
        
        #if wind[tt,tm]>=25 and rmw[tt,tm]<=0 and tdrmw[tt]==0 :
        #    print('no tdrmw',name[tt],year[tt],wind[tt,tm-1],wind[tt,tm])
        if wind[tt,tm]>=25 : #wind>=25kt, wind value, (first wind had already >25kt also)
            if tdmonth[tt]==0 : #first td
                
                tdtimestr[tt]=timestr[tt,tm]
                tdtimestr_sat[tt]=ttimestr[tt,tm]
                tdhour[tt]=hour[tt,tm]
                tdmonth[tt]=float(month[tt,tm])
                tdyear[tt]=year[tt,tm]
                tdlat[tt]=lat[tt,tm]
                tdlon[tt]=lon[tt,tm]
                tdtttimestr[tt]=ttimestr[tt,tm]
                tdroci[tt]=roci[tt,tm]
                
                #if wind[tt,tm]>25:
                #    print('>25',name[tt],year[tt],time[tt,tm-1],wind[tt,tm-1],wind[tt,tm])
                    
            if  dist2land[tt,tm]<200 :   #first ts tsrmw[tt]!=0 and
                land=2
            if  landfall[tt,tm]==0 :#after landfall
                land=1
                
            if tsrmw[tt]==0 and rmw[tt,tm]>0:
                tdrmw[tt]+=rmw[tt,tm]
                tdts+=1
                    
        #if wind[tt,tm]>34 and rmw[tt,tm]<=0 and tsrmw[tt]==0 :
        #    print('no tsrmw',name[tt],year[tt],wind[tt,tm-1],wind[tt,tm])
        if wind[tt,tm]>34 and rmw[tt,tm]>0 : #wind>34kt, wind value, rmw value
            if tsrmw[tt]==0 : #first ts
                #if alltm>0:
                #    tdroci[tt]=alltdroci/alltm
                tsrmw[tt]=rmw[tt,tm]
                tstimestr[tt]=timestr[tt,tm]
                tshour[tt]=hour[tt,tm]
                tslat[tt]=lat[tt,tm]
                tslon[tt]=lon[tt,tm]
                
                tstime=tm
                    
                #if wind[tt,tm-1]>34:
                #    print('>34',name[tt],year[tt],wind[tt,tm-1],wind[tt,tm])
        
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
                    
                    #irmw[tt]=rmw[tt,tm]
                    itimestr[tt]=timestr[tt,tm] #+0 =start of intensifying, +8 =end of intensifying
                    ihour[tt]=hour[tt,tm]
                    ilat[tt]=lat[tt,tm]
                    ilon[tt]=lon[tt,tm]
                    
                    itime=tm
                    
                    #VWS calculating shouldn't exclude unfavorable factors
                    if forVWS==0:
                        if landfall[tt,tm+8]==0 or land==1 or land==2:
                            if landfall[tt,tm+8]==0:
                                env_lf_i[tt]=1
                            if land==1:
                                env_lf_s[tt]=1
                            if land==2:
                                env_dl[tt]=1
                            LMIR[tt]=np.nan
                            #irmw[tt]=np.nan
                            itimestr[tt]=np.nan
                            ihour[tt]=np.nan
                            ilat[tt]=np.nan
                            ilon[tt]=np.nan
                            
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
                            #TS composite should exclude unfavorable factors
                            tsrmw[tt]=np.nan
                            tstimestr[tt]=np.nan
                            tshour[tt]=np.nan
                            tslat[tt]=np.nan
                            tslon[tt]=np.nan
                            #TD composite should exclude unfavorable factors
                            tdroci[tt]=np.nan
                            tdrmw[tt]=np.nan
                            tdtimestr[tt]=np.nan
                            tdtimestr_sat[tt]=np.nan
                            tdhour[tt]=np.nan
                            tdmonth[tt]=np.nan
                            tdyear[tt]=np.nan
                            tdlat[tt]=np.nan
                            tdlon[tt]=np.nan
                            
                            tdtttimestr[tt]=np.nan
                            #if land==1:
                            #    print(name[tt],year[tt],time[tt,0])
        
        #if status[tt,tm,0]=='M':
            #tdew[tt]=-2
            #print(tt,status[tt,tm],name[tt],year[tt],tsrmw[tt],LMIR[tt])
        
        #if status[tt,tm,0]=='E':# and status[tt,tm,1]=='W':
        #    print(tt,status[tt,tm],name[tt],year[tt],tsrmw[tt],LMIR[tt])
        
    
    if  tdlon[tt]<=100 or tdlat[tt]>=30:
        #if  tslat[tt]>30:
        #    print('not in domain',name[tt],year[tt],tslon[tt],tslat[tt])
        LMIR[tt]=np.nan
        #irmw[tt]=np.nan
        itimestr[tt]=np.nan
        ihour[tt]=np.nan
        ilat[tt]=np.nan
        ilon[tt]=np.nan
        
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
        #TS composite should exclude unfavorable factors
        tsrmw[tt]=np.nan
        tstimestr[tt]=np.nan
        tshour[tt]=np.nan
        tslat[tt]=np.nan
        tslon[tt]=np.nan
        #TD composite should exclude unfavorable factors
        tdroci[tt]=np.nan
        tdrmw[tt]=np.nan
        tdtimestr[tt]=np.nan
        tdtimestr_sat[tt]=np.nan
        tdhour[tt]=np.nan
        tdmonth[tt]=np.nan
        tdyear[tt]=np.nan
        tdlat[tt]=np.nan
        tdlon[tt]=np.nan
        
        tdtttimestr[tt]=np.nan
        
    
    
    #VWS calculating shouldn't exclude unfavorable factors
    if forVWS==0:
        #if tslat[tt]>=25 and tslat[tt]<=26 and ivws[tt]<=8:
        #    print(name[tt],year[tt],tsrmw[tt],LMIR[tt])
        if tslat[tt]>26 or ivws[tt]>8 or tslon[tt]<121 :#!!!remember to 26
            if tslat[tt]>26:
                env_tslat[tt]=1
            if tslon[tt]<121:
                env_tslon[tt]=1
            if ivws[tt]>8:
                env_ivws[tt]=1
            LMIR[tt]=np.nan
            #irmw[tt]=np.nan
            itimestr[tt]=np.nan
            ihour[tt]=np.nan
            ilat[tt]=np.nan
            ilon[tt]=np.nan
            
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
            
            #TS composite should exclude unfavorable factors
            tsrmw[tt]=np.nan
            tstimestr[tt]=np.nan
            tshour[tt]=np.nan
            tslat[tt]=np.nan
            tslon[tt]=np.nan
            #TD composite should exclude unfavorable factors
            tdroci[tt]=np.nan
            tdrmw[tt]=np.nan
            tdtimestr[tt]=np.nan
            tdtimestr_sat[tt]=np.nan
            tdhour[tt]=np.nan
            tdmonth[tt]=np.nan
            tdyear[tt]=np.nan
            tdlat[tt]=np.nan
            tdlon[tt]=np.nan
            
            tdtttimestr[tt]=np.nan
                    
        
    tdrmw[tt]=tdrmw[tt]/tdts
    ts2i[tt]=(itime-tstime)*3
    
    #if ts2i[tt]<=-50 and LMIR[tt]!=0 : #just show
        #print(name[tt],year[tt],time[tt,0])
    
    #if tslat[tt]>=25 :
    #    print(name[tt],itimestr[tt],ihour[tt],max(wind[tt]))
        
    #exclude no data
    if LMIR[tt]==0 or tsrmw[tt]==0 or irmw[tt]==0 or ts2i[tt]==1 or np.isnan(LMIR[tt])==True :
        LMIR[tt]=np.nan
        #print(name[tt],year[tt],time[tt,0],max(wind[tt]))
        #irmw[tt]=np.nan
        ts2i[tt]=np.nan
        
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
        
        tsrmw[tt]=np.nan
        tstimestr[tt]=np.nan
        tshour[tt]=np.nan
        tslat[tt]=np.nan
        tslon[tt]=np.nan
        
        tdroci[tt]=np.nan
        tdrmw[tt]=np.nan
        tdtimestr[tt]=np.nan
        tdtimestr_sat[tt]=np.nan
        tdhour[tt]=np.nan
        tdmonth[tt]=np.nan
        tdyear[tt]=np.nan
        tdlat[tt]=np.nan
        tdlon[tt]=np.nan
        
        tdtttimestr[tt]=np.nan
        
        #print('no data',name[tt],year[tt])
        
    #for diff. VWS
    if ivws[tt]<4.5 :
        l_LMIR[tt]=LMIR[tt]
        l_tsrmw[tt]=tsrmw[tt]
        l_irmw[tt]=irmw[tt]
        #if l_LMIR[tt]<=25 and l_tsrmw[tt]<100:
        #    print(name[tt],itimestr[tt],ihour[tt],max(wind[tt]))
    elif ivws[tt]<10 :
        m_LMIR[tt]=LMIR[tt]
        m_tsrmw[tt]=tsrmw[tt]
        m_irmw[tt]=irmw[tt]
    if ivws[tt]>10 :
        h_LMIR[tt]=LMIR[tt]
        h_tsrmw[tt]=tsrmw[tt]
        h_irmw[tt]=irmw[tt]
        
    if l_LMIR[tt]==0 or l_tsrmw[tt]==0 or l_irmw[tt]==0 or np.isnan(l_LMIR[tt])==True :
        l_LMIR[tt]=np.nan
        l_tsrmw[tt]=np.nan
        l_irmw[tt]=np.nan
    if m_LMIR[tt]==0 or m_tsrmw[tt]==0 or m_irmw[tt]==0 or np.isnan(m_LMIR[tt])==True :
        m_LMIR[tt]=np.nan
        m_tsrmw[tt]=np.nan
        m_irmw[tt]=np.nan
    if h_LMIR[tt]==0 or h_tsrmw[tt]==0 or h_irmw[tt]==0 or np.isnan(h_LMIR[tt])==True :
        h_LMIR[tt]=np.nan
        h_tsrmw[tt]=np.nan
        h_irmw[tt]=np.nan
        
    irmw_c[tt]=irmw_24[tt]-irmw[tt]
    
    #if irh[tt]<60 and LMIR[tt]>30:
    #    print(name[tt],year[tt],tsrmw[tt],LMIR[tt],ivws[tt],irh[tt])
        
    #if tdroci[tt]>=500:
    #    print(tt,name[tt],year[tt],tsrmw[tt])
        
    #if tsrmw[tt]>=150 and tdew[tt]>-2:
    #    print(tt,name[tt],year[tt],tsrmw[tt],tdew[tt])
    
    #if tsrmw[tt]<50 and tdew[tt]<=-2:
    #    print('small MD',tt,name[tt],year[tt],tsrmw[tt],tdroci[tt],tdew[tt])
    
    #if tt>400:
    #    print(tt,tdyear[tt],name[tt],LMIR[tt])
        
    #if tslat[tt]<=26 and ilat[tt]>30:
    #    print(tt,tdyear[tt],name[tt],ilat[tt],tslat[tt],tdlat[tt])
        
#plt.plot(t,LMIR,"r.")
#print(tstimestr)
        
# 將變數中的0值替換為np.nan
variables = [wind_t24, wind_0, wind_24, wind_48, wind_72, irmw_t24, irmw, irmw_24, irmw_48, irmw_72]

for var in variables:
    var[var == 0] = np.nan

# npy
np.save('LMIR',LMIR)
np.save('tsrmw',tsrmw)
np.save('tstimestr',tstimestr)
np.save('tshour',tshour)
np.save('tslat',tslat)
np.save('tslon',tslon)
np.save('itimestr',itimestr)
np.save('ihour',ihour)
np.save('ilat',ilat)
np.save('ilon',ilon)
#np.save('irmw',irmw)
np.save('tdroci',tdroci)
np.save('tdrmw',tdrmw)
np.save('tdtimestr',tdtimestr)
np.save('tdtimestr_sat',tdtimestr_sat)
np.save('tdhour',tdhour)
np.save('tdmonth',tdmonth)
np.save('tdyear',tdyear)
np.save('tdlat',tdlat)
np.save('tdlon',tdlon)
np.save('name',name)
np.save('tdtttimestr',tdtttimestr)

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

if forVWS==0:
    np.save('env_lf_i',env_lf_i)
    np.save('env_lf_s',env_lf_s)
    np.save('env_dl',env_dl)
    np.save('env_tslat',env_tslat)
    np.save('env_tslon',env_tslon)
    np.save('env_ivws',env_ivws)

#filt NaN to regress
x=tsrmw[~np.isnan(tsrmw)]
y=LMIR[~np.isnan(LMIR)]
x=1/x#np.log(x)
tsnum,tscor=reg(x,y)
tscorr_s, tsp_s = spearmanr(x,y)
tscorr_p, tsp_p = pearsonr(x,y)
#tsnum,tscor=curve_fit(func, x, y, (700, 10))
tslinex=np.arange(300)
tsliney=tsnum[0]/(tslinex)+tsnum[1]

x=wind_0[~np.isnan(wind_0)]
y=LMIR[~np.isnan(LMIR)]
#x=1/x#np.log(x)
iwnum,iwcor=reg(x,y)
iwcorr_s, iwp_s = spearmanr(x,y)
iwcorr_p, iwp_p = pearsonr(x,y)

x=irmw[~np.isnan(irmw)]
y=LMIR[~np.isnan(LMIR)]
x=1/x#np.log(x)
inum,icor=reg(x,y)
icorr_s, ip_s = spearmanr(x,y)
icorr_p, ip_p = pearsonr(x,y)
#inum,icor=curve_fit(func, x, y, (700, 10))
ilinex=np.arange(300)
iliney=inum[0]/(ilinex)+inum[1]




#draw https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.streamplot.html
fig,ax=plt.subplots(2,2,sharex=False,sharey=False,figsize=(22,15))

# 設定座標軸上數字的大小
for row in ax:
    for axis in row:
        axis.tick_params(labelsize=15)

#散佈圖
#genesis type
ew=0
mc=0
ms=0
mg=0
gen=['ew','mc','ms','md']
xy=['tsrmw','LMIR']
xy3=['tsrmw','LMIR']
for i in range(len(gen)):
    for j in range(len(xy)):
        exec(gen[i]+xy[j]+'=np.zeros(tstimestr.size)')
for tt in range(n_storm-start):
    if tdew[tt]==1:
        for j in range(len(xy)):
            exec(gen[0]+xy[j]+'['+str(tt)+']='+xy[j]+'['+str(tt)+']')
        if ew>=1:
            ax[1,1].scatter(tsrmw[tt], LMIR[tt],100, color='red', marker='s', alpha=0.2, linewidths=0)
            ew+=1
        else:
            ax[1,1].scatter(tsrmw[tt], LMIR[tt],100, color='red', marker='s', alpha=0.2, linewidths=0, label='EW')
            ew=1
            
    elif tdew[tt]==-1:
        for j in range(len(xy)):
            exec(gen[1]+xy[j]+'['+str(tt)+']='+xy[j]+'['+str(tt)+']')
        if mc>=1:
            ax[1,1].scatter(tsrmw[tt], LMIR[tt],100, color='purple', marker='o', alpha=0.2, linewidths=0)
            mc+=1
        else:
            ax[1,1].scatter(tsrmw[tt], LMIR[tt],100, color='purple', marker='o', alpha=0.2, linewidths=0, label='MC')
            mc=1
            
    elif tdew[tt]==-1.5:
        for j in range(len(xy)):
            exec(gen[2]+xy[j]+'['+str(tt)+']='+xy[j]+'['+str(tt)+']')
        if ms>=1:
            ax[1,1].scatter(tsrmw[tt], LMIR[tt],100, color='green', marker='o', alpha=0.2, linewidths=0)
            ms+=1
        else:
            ax[1,1].scatter(tsrmw[tt], LMIR[tt],100, color='green', marker='o', alpha=0.2, linewidths=0, label='MS')
            ms=1
            
    elif tdew[tt]==-2:
        for j in range(len(xy)):
            exec(gen[3]+xy[j]+'['+str(tt)+']='+xy[j]+'['+str(tt)+']')
        if mg>=1:
            ax[1,1].scatter(tsrmw[tt], LMIR[tt],100, color='blue', marker='v', alpha=0.2, linewidths=0)
            mg+=1
        else:
            ax[1,1].scatter(tsrmw[tt], LMIR[tt],100, color='blue', marker='v', alpha=0.2, linewidths=0, label='MD')
            mg=1
            

    else:
        ax[1,1].scatter(tsrmw[tt], LMIR[tt],100, color='black', marker='*', alpha=0.2, linewidths=1)
        


#calculate crossline
for i in range(len(gen)):
    for j in range(len(xy)):
        exec(gen[i]+xy[j]+'['+gen[i]+xy[j]+'==0]=np.nan')
        exec(gen[i]+xy[j]+'mean=np.nanmean('+gen[i]+xy[j]+')')
        exec(gen[i]+xy[j]+'std=np.nanstd('+gen[i]+xy[j]+')')


#draw crossline
for i in range(len(gen)):
    for j in range(len(xy)):
        if j==0:
            if i==0:
                exec('ax[1,1].plot(['+gen[i]+xy[j]+'mean-'+gen[i]+xy[j]+'std,'+gen[i]+xy[j]+'mean+'+gen[i]+xy[j]+'std],['+
                 gen[i]+xy[j+1]+'mean,'+gen[i]+xy[j+1]+'mean],color="red", linewidth=2)')
            if i==1:
                exec('ax[1,1].plot(['+gen[i]+xy[j]+'mean-'+gen[i]+xy[j]+'std,'+gen[i]+xy[j]+'mean+'+gen[i]+xy[j]+'std],['+
                 gen[i]+xy[j+1]+'mean,'+gen[i]+xy[j+1]+'mean],color="purple", linewidth=2)')
            if i==2:
                exec('ax[1,1].plot(['+gen[i]+xy[j]+'mean-'+gen[i]+xy[j]+'std,'+gen[i]+xy[j]+'mean+'+gen[i]+xy[j]+'std],['+
                 gen[i]+xy[j+1]+'mean,'+gen[i]+xy[j+1]+'mean],color="green", linewidth=2)')
            if i==3:
                exec('ax[1,1].plot(['+gen[i]+xy[j]+'mean-'+gen[i]+xy[j]+'std,'+gen[i]+xy[j]+'mean+'+gen[i]+xy[j]+'std],['+
                 gen[i]+xy[j+1]+'mean,'+gen[i]+xy[j+1]+'mean],color="blue", linewidth=2)')
            
        if j==1:
            if i==0:
                exec('ax[1,1].plot(['+gen[i]+xy[j-1]+'mean,'+gen[i]+xy[j-1]+'mean],['+
                 gen[i]+xy[j]+'mean-'+gen[i]+xy[j]+'std,'+gen[i]+xy[j]+'mean+'+gen[i]+xy[j]+'std],color="red", linewidth=2)')
            if i==1:
                exec('ax[1,1].plot(['+gen[i]+xy[j-1]+'mean,'+gen[i]+xy[j-1]+'mean],['+
                 gen[i]+xy[j]+'mean-'+gen[i]+xy[j]+'std,'+gen[i]+xy[j]+'mean+'+gen[i]+xy[j]+'std],color="purple", linewidth=2)')
            if i==2:
                exec('ax[1,1].plot(['+gen[i]+xy[j-1]+'mean,'+gen[i]+xy[j-1]+'mean],['+
                 gen[i]+xy[j]+'mean-'+gen[i]+xy[j]+'std,'+gen[i]+xy[j]+'mean+'+gen[i]+xy[j]+'std],color="green", linewidth=2)')
            if i==3:
                exec('ax[1,1].plot(['+gen[i]+xy[j-1]+'mean,'+gen[i]+xy[j-1]+'mean],['+
                 gen[i]+xy[j]+'mean-'+gen[i]+xy[j]+'std,'+gen[i]+xy[j]+'mean+'+gen[i]+xy[j]+'std],color="blue", linewidth=2)')
            

#ax1.plot(tslinex,tsliney, color='red')

#圖例，標題等
#ax1.text(155,90,'y='+str(np.around(tsnum[0],decimals=4))+'ln(x)+'+str(np.around(tsnum[1],decimals=4)))
#ax1.text(155,85,'correlation='+str(np.around(tscor,decimals=4)))
#ax1.text(230,70.5,str(ew),color='red')
#ax1.text(230,93,str(mc),color='purple')
#ax1.text(230,78,str(mg),color='blue')
#ax1.text(230,85.5,str(ms),color='green')
#ax[1,1].legend()
ax[1,1].grid(True)
ax[1,1].set_xlim(left=0, right=300)
ax[1,1].set_ylim(bottom=-0, top=100)
#ax1.set_xscale("symlog")
#ax1.set_yscale("symlog")
ax[1,1].set_xlabel('TS_RMW (km)', fontsize=25)
#ax[1,1].set_ylabel('LMIR (kt/24h)')
#ax[1,1].set_title('TS_RMW vs. Lifetime Max. Intensification Rate')
#plt.savefig('fig/LMIR_tsrmw_sh.png', dpi=300)
#plt.show()





# 製作figure  
#fig2 = plt.figure()   

#圖表的設定
#ax2 = fig2.add_subplot(1, 1, 1)

#散佈圖
ax[0,0].scatter(irmw, LMIR,100, color='red', alpha=0.3, linewidths=2, edgecolors='black')
ax[0,0].scatter(203, 95,100, color='red', alpha=0.3, linewidths=2, edgecolors='black')
for i in range(3):
    ax[0,0].scatter(203, 89,100, color='red', alpha=0.3, linewidths=2, edgecolors='black')
for i in range(5):
    ax[0,0].scatter(203, 83,100, color='red', alpha=0.3, linewidths=2, edgecolors='black')

#圖例，標題等
ax[0,0].text(210,93,'1  case', fontsize=15)
ax[0,0].text(210,87,'3  cases', fontsize=15)
ax[0,0].text(210,81,'5  cases', fontsize=15)

ax[0,0].plot(ilinex,iliney, color='black')

#圖例，標題等
ax[0,0].text(155,69,'y='+str(np.around(inum[0],decimals=4))+'/x+'+str(np.around(inum[1],decimals=4)), fontsize=15)
#ax2.text(155,85,'correlation='+str(np.around(icor,decimals=4)))
ax[0,0].grid(True)
ax[0,0].set_xlim(left=0, right=300)
ax[0,0].set_ylim(bottom=-0, top=100)
#ax2.set_xscale("symlog")
#ax2.set_yscale("symlog")
ax[0,0].set_xlabel('I_RMW (km)', fontsize=25)
ax[0,0].set_ylabel('$LMIR (kt \\, day^{-1})$', fontsize=25)
#ax2.set_title('I_RMW vs. Lifetime Max. Intensification Rate')
#plt.savefig('fig/LMIR_irmw_sh_noew.png', dpi=300)
#plt.show()




# 製作figure  
#fig4 = plt.figure()   

#圖表的設定
#ax4 = fig4.add_subplot(1, 1, 1)

#散佈圖
ax[1,0].scatter(tsrmw, LMIR,100, color='red', alpha=0.3, linewidths=2, edgecolors='black')
#ax4.plot(ilinex,iliney, color='red')

#圖例，標題等
#ax4.text(155,90,'y='+str(np.around(inum[0],decimals=4))+'/x+'+str(np.around(inum[1],decimals=4)))
#ax2.text(155,85,'correlation='+str(np.around(icor,decimals=4)))
ax[1,0].grid(True)
ax[1,0].set_xlim(left=0, right=300)
ax[1,0].set_ylim(bottom=-0, top=100)
#ax2.set_xscale("symlog")
#ax2.set_yscale("symlog")
ax[1,0].set_xlabel('TS_RMW (km)', fontsize=25)
ax[1,0].set_ylabel('$LMIR (kt \\, day^{-1})$', fontsize=25)
#ax4.set_title('TS_RMW vs. Lifetime Max. Intensification Rate')
#plt.savefig('fig/LMIR_tsrmw_sh_noew.png', dpi=300)
#plt.show() 



# 製作figure  
#fig6 = plt.figure()   

#圖表的設定
#ax6 = fig6.add_subplot(1, 1, 1)

#散佈圖
#genesis type
ew=0
mc=0
ms=0
mg=0
gen=['ew','mc','ms','md']
xy=['irmw','LMIR']
xy3=['irmw','LMIR']
for i in range(len(gen)):
    for j in range(len(xy)):
        exec(gen[i]+xy[j]+'=np.zeros(tstimestr.size)')
for tt in range(n_storm-start):
    if tdew[tt]==1:
        for j in range(len(xy)):
            exec(gen[0]+xy[j]+'['+str(tt)+']='+xy[j]+'['+str(tt)+']')
        if ew>=1:
            ax[0,1].scatter(irmw[tt], LMIR[tt],100, color='red', marker='s', alpha=0.2, linewidths=0)
            ew+=1
        else:
            ax[0,1].scatter(irmw[tt], LMIR[tt],100, color='red', marker='s', alpha=0.2, linewidths=0, label='EW')
            ew=1
            
    elif tdew[tt]==-1:
        for j in range(len(xy)):
            exec(gen[1]+xy[j]+'['+str(tt)+']='+xy[j]+'['+str(tt)+']')
        if mc>=1:
            ax[0,1].scatter(irmw[tt], LMIR[tt],100, color='purple', marker='o', alpha=0.2, linewidths=0)
            mc+=1
        else:
            ax[0,1].scatter(irmw[tt], LMIR[tt],100, color='purple', marker='o', alpha=0.2, linewidths=0, label='MC')
            mc=1
            
    elif tdew[tt]==-1.5:
        for j in range(len(xy)):
            exec(gen[2]+xy[j]+'['+str(tt)+']='+xy[j]+'['+str(tt)+']')
        if ms>=1:
            ax[0,1].scatter(irmw[tt], LMIR[tt],100, color='green', marker='o', alpha=0.2, linewidths=0)
            ms+=1
        else:
            ax[0,1].scatter(irmw[tt], LMIR[tt],100, color='green', marker='o', alpha=0.2, linewidths=0, label='MS')
            ms=1
            
    elif tdew[tt]==-2:
        for j in range(len(xy)):
            exec(gen[3]+xy[j]+'['+str(tt)+']='+xy[j]+'['+str(tt)+']')
        if mg>=1:
            ax[0,1].scatter(irmw[tt], LMIR[tt],100, color='blue', marker='v', alpha=0.2, linewidths=0)
            mg+=1
        else:
            ax[0,1].scatter(irmw[tt], LMIR[tt],100, color='blue', marker='v', alpha=0.2, linewidths=0, label='MD')
            mg=1
            

    else:
        ax[0,1].scatter(irmw[tt], LMIR[tt],100, color='black', marker='*', alpha=0.2, linewidths=1)
        


#calculate crossline
for i in range(len(gen)):
    for j in range(len(xy)):
        exec(gen[i]+xy[j]+'['+gen[i]+xy[j]+'==0]=np.nan')
        exec(gen[i]+xy[j]+'mean=np.nanmean('+gen[i]+xy[j]+')')
        exec(gen[i]+xy[j]+'std=np.nanstd('+gen[i]+xy[j]+')')


#draw crossline
for i in range(len(gen)):
    for j in range(len(xy)):
        if j==0:
            if i==0:
                exec('ax[0,1].plot(['+gen[i]+xy[j]+'mean-'+gen[i]+xy[j]+'std,'+gen[i]+xy[j]+'mean+'+gen[i]+xy[j]+'std],['+
                 gen[i]+xy[j+1]+'mean,'+gen[i]+xy[j+1]+'mean],color="red", linewidth=2)')
            if i==1:
                exec('ax[0,1].plot(['+gen[i]+xy[j]+'mean-'+gen[i]+xy[j]+'std,'+gen[i]+xy[j]+'mean+'+gen[i]+xy[j]+'std],['+
                 gen[i]+xy[j+1]+'mean,'+gen[i]+xy[j+1]+'mean],color="purple", linewidth=2)')
            if i==2:
                exec('ax[0,1].plot(['+gen[i]+xy[j]+'mean-'+gen[i]+xy[j]+'std,'+gen[i]+xy[j]+'mean+'+gen[i]+xy[j]+'std],['+
                 gen[i]+xy[j+1]+'mean,'+gen[i]+xy[j+1]+'mean],color="green", linewidth=2)')
            if i==3:
                exec('ax[0,1].plot(['+gen[i]+xy[j]+'mean-'+gen[i]+xy[j]+'std,'+gen[i]+xy[j]+'mean+'+gen[i]+xy[j]+'std],['+
                 gen[i]+xy[j+1]+'mean,'+gen[i]+xy[j+1]+'mean],color="blue", linewidth=2)')
            
        if j==1:
            if i==0:
                exec('ax[0,1].plot(['+gen[i]+xy[j-1]+'mean,'+gen[i]+xy[j-1]+'mean],['+
                 gen[i]+xy[j]+'mean-'+gen[i]+xy[j]+'std,'+gen[i]+xy[j]+'mean+'+gen[i]+xy[j]+'std],color="red", linewidth=2)')
            if i==1:
                exec('ax[0,1].plot(['+gen[i]+xy[j-1]+'mean,'+gen[i]+xy[j-1]+'mean],['+
                 gen[i]+xy[j]+'mean-'+gen[i]+xy[j]+'std,'+gen[i]+xy[j]+'mean+'+gen[i]+xy[j]+'std],color="purple", linewidth=2)')
            if i==2:
                exec('ax[0,1].plot(['+gen[i]+xy[j-1]+'mean,'+gen[i]+xy[j-1]+'mean],['+
                 gen[i]+xy[j]+'mean-'+gen[i]+xy[j]+'std,'+gen[i]+xy[j]+'mean+'+gen[i]+xy[j]+'std],color="green", linewidth=2)')
            if i==3:
                exec('ax[0,1].plot(['+gen[i]+xy[j-1]+'mean,'+gen[i]+xy[j-1]+'mean],['+
                 gen[i]+xy[j]+'mean-'+gen[i]+xy[j]+'std,'+gen[i]+xy[j]+'mean+'+gen[i]+xy[j]+'std],color="blue", linewidth=2)')
            

#ax1.plot(tslinex,tsliney, color='red')

#圖例，標題等
#ax1.text(155,90,'y='+str(np.around(tsnum[0],decimals=4))+'ln(x)+'+str(np.around(tsnum[1],decimals=4)))
#ax1.text(155,85,'correlation='+str(np.around(tscor,decimals=4)))
#ax1.text(230,70.5,str(ew),color='red')
#ax1.text(230,93,str(mc),color='purple')
#ax1.text(230,78,str(mg),color='blue')
#ax1.text(230,85.5,str(ms),color='green')
ax[0,1].legend(fontsize=15)#labels=['MD', 'MS', 'MC', 'EW'],
ax[0,1].grid(True)
ax[0,1].set_xlim(left=0, right=300)
ax[0,1].set_ylim(bottom=-0, top=100)
#ax1.set_xscale("symlog")
#ax1.set_yscale("symlog")
ax[0,1].set_xlabel('I_RMW (km)', fontsize=25)
#ax6.set_ylabel('LMIR (kt/24h)')
#ax6.set_title('I_RMW vs. Lifetime Max. Intensification Rate')
#plt.savefig('fig/LMIR_irmw_sh.png', dpi=300)

title_num=[['a','b'],['c','d']]
for i in range(2):
    for j in range(2):
        ax[i,j].set_title('('+title_num[i][j]+')',loc='left', fontsize=25)

#fig.subplots_adjust(right=0.95)
fig.subplots_adjust(left=0.087, right=0.97, bottom=0.087, top=0.97)
plt.savefig('fig/LMIR_RMW.png', dpi=600)
plt.show()



fig, ax = plt.subplots(2, 1, sharex=False, sharey=False, figsize=(11, 15))

# 设置座标轴上数字的大小
for axis in ax:
    axis.tick_params(labelsize=15)

# 上方散佈图
ax[0].scatter(irmw, LMIR, 100, color='red', alpha=0.3, linewidths=2, edgecolors='black')
ax[0].scatter(203, 95, 100, color='red', alpha=0.3, linewidths=2, edgecolors='black')
for i in range(3):
    ax[0].scatter(203, 89, 100, color='red', alpha=0.3, linewidths=2, edgecolors='black')
for i in range(5):
    ax[0].scatter(203, 83, 100, color='red', alpha=0.3, linewidths=2, edgecolors='black')

# 图例，标题等
ax[0].text(210, 93, '1  case', fontsize=15)
ax[0].text(210, 87, '3  cases', fontsize=15)
ax[0].text(210, 81, '5  cases', fontsize=15)

ax[0].plot(ilinex, iliney, color='black')

# 图例，标题等
ax[0].text(155, 69, 'y=' + str(np.around(inum[0], decimals=4)) + '/x+' + str(np.around(inum[1], decimals=4)), fontsize=15)
ax[0].grid(True)
ax[0].set_xlim(left=0, right=300)
ax[0].set_ylim(bottom=-0, top=100)
ax[0].set_xlabel('I_RMW (km)', fontsize=25)
ax[0].set_ylabel('$LMIR (kt \\, day^{-1})$', fontsize=25)
ax[0].set_title('(a)', loc='left', fontsize=25)

# 下方散佈图
ax[1].scatter(tsrmw, LMIR, 100, color='red', alpha=0.3, linewidths=2, edgecolors='black')

ax[1].grid(True)
ax[1].set_xlim(left=0, right=300)
ax[1].set_ylim(bottom=-0, top=100)
ax[1].set_xlabel('TS_RMW (km)', fontsize=25)
ax[1].set_ylabel('$LMIR (kt \\, day^{-1})$', fontsize=25)
ax[1].set_title('(b)', loc='left', fontsize=25)

fig.subplots_adjust(left=0.087, right=0.97, bottom=0.087, top=0.97)
plt.savefig('fig/LMIR_RMW_ALL.png', dpi=600)
plt.show()




# 製作figure  
fig3, ax = plt.subplots(1, 1, figsize=(11, 8))

# 设置座标轴上数字的大小
ax.tick_params(labelsize=15)
fig3 = plt.figure()   

#圖表的設定
#ax3 = fig3.add_subplot(1, 1, 1)

#散佈圖
ax.scatter(ts2i, LMIR, color='red', alpha=0.3, linewidths=2, edgecolors='black')

#圖例，標題等
ax.grid(True)
ax.set_xlim(left=-50, right=300)
ax.set_ylim(bottom=-0, top=100)
ax.set_xlabel('TS_TIME to I_TIME (hr)', fontsize=25)
ax.set_ylabel('$LMIR (kt \\, day^{-1})$', fontsize=25)
ax.set_title('TS_TIME to I_TIME (hr) vs. $LMIR (kt \\, day^{-1})$')
plt.savefig('fig/LMIR_ts2i.png', dpi=600)
plt.show() 



# 创建子图
fig, ax = plt.subplots(2, 2, figsize=(22, 15))

# 设置座标轴上数字的大小
for row in ax:
    for axis in row:
        axis.tick_params(labelsize=15)

# 初始化计数器
ew = mc = ms = mg = 0
gen = ['ew', 'mc', 'ms', 'md']
xy = ['tsrmw', 'LMIR']
xy3 = ['tsrmw', 'LMIR']
axes = {'ew': ax[0, 0], 'mc': ax[0, 1], 'ms': ax[1, 0], 'md': ax[1, 1]}
colors = {'ew': 'red', 'mc': 'purple', 'ms': 'green', 'md': 'blue'}
markers = {'ew': 's', 'mc': 'o', 'ms': 'o', 'md': 'v'}
labels = {'ew': 'EW', 'mc': 'MC', 'ms': 'MS', 'md': 'MD'}

# 初始化数组
for i in range(len(gen)):
    for j in range(len(xy)):
        exec(gen[i] + xy[j] + '=np.zeros(tstimestr.size)')

# 绘制散点图
for tt in range(n_storm - start):
    if tdew[tt] == 1:
        for j in range(len(xy)):
            exec(gen[0] + xy[j] + '[' + str(tt) + ']=' + xy[j] + '[' + str(tt) + ']')
        ax[0, 0].scatter(tsrmw[tt], LMIR[tt], 100, color='red', marker='o', alpha=0.3, linewidths=2, edgecolors='black', label=labels['ew'] if ew == 0 else "")
        ew += 1
    elif tdew[tt] == -1:
        for j in range(len(xy)):
            exec(gen[1] + xy[j] + '[' + str(tt) + ']=' + xy[j] + '[' + str(tt) + ']')
        ax[0, 1].scatter(tsrmw[tt], LMIR[tt], 100, color='purple', marker='o', alpha=0.3, linewidths=2, edgecolors='black', label=labels['mc'] if mc == 0 else "")
        mc += 1
    elif tdew[tt] == -1.5:
        for j in range(len(xy)):
            exec(gen[2] + xy[j] + '[' + str(tt) + ']=' + xy[j] + '[' + str(tt) + ']')
        ax[1, 0].scatter(tsrmw[tt], LMIR[tt], 100, color='green', marker='o', alpha=0.3, linewidths=2, edgecolors='black', label=labels['ms'] if ms == 0 else "")
        ms += 1
    elif tdew[tt] == -2:
        for j in range(len(xy)):
            exec(gen[3] + xy[j] + '[' + str(tt) + ']=' + xy[j] + '[' + str(tt) + ']')
        ax[1, 1].scatter(tsrmw[tt], LMIR[tt], 100, color='blue', marker='o', alpha=0.3, linewidths=2, edgecolors='black', label=labels['md'] if mg == 0 else "")
        mg += 1
    else:
        ax[1, 1].scatter(tsrmw[tt], LMIR[tt], 100, color='black', marker='*', alpha=0.3, linewidths=2, edgecolors='black')

# 添加图例
for key, axis in axes.items():
    axis.scatter(203, 95, 100, color=colors[key], alpha=0.3, linewidths=2, edgecolors='black')
    for _ in range(3):
        axis.scatter(203, 89, 100, color=colors[key], alpha=0.3, linewidths=2, edgecolors='black')
    for _ in range(5):
        axis.scatter(203, 83, 100, color=colors[key], alpha=0.3, linewidths=2, edgecolors='black')
    axis.text(210, 93, '1  case', fontsize=15)
    axis.text(210, 87, '3  cases', fontsize=15)
    axis.text(210, 81, '5  cases', fontsize=15)

# 计算均值和标准差
for i in range(len(gen)):
    for j in range(len(xy)):
        exec(gen[i] + xy[j] + '[' + gen[i] + xy[j] + '==0]=np.nan')
        exec(gen[i] + xy[j] + 'mean=np.nanmean(' + gen[i] + xy[j] + ')')
        exec(gen[i] + xy[j] + 'std=np.nanstd(' + gen[i] + xy[j] + ')')

# 绘制水平和垂直均值和标准差线
for i in range(len(gen)):
    # 绘制水平线
    exec(f"ax[{i//2},{i%2}].plot([{gen[i]}{xy[0]}mean-{gen[i]}{xy[0]}std, {gen[i]}{xy[0]}mean+{gen[i]}{xy[0]}std], [{gen[i]}{xy[1]}mean, {gen[i]}{xy[1]}mean], color=colors['{gen[i]}'], linewidth=4)")
    # 绘制垂直线
    exec(f"ax[{i//2},{i%2}].plot([{gen[i]}{xy[0]}mean, {gen[i]}{xy[0]}mean], [{gen[i]}{xy[1]}mean-{gen[i]}{xy[1]}std, {gen[i]}{xy[1]}mean+{gen[i]}{xy[1]}std], color=colors['{gen[i]}'], linewidth=4)")


# 设置子图标题、网格、轴标签等
titles = ['EW', 'MC', 'MS', 'MD']
title_num=[['a','b'],['c','d']]
for i in range(2):
    
    ax[i, 0].set_ylabel('$LMIR (kt \\, day^{-1})$', fontsize=25)
    for j in range(2):
        ax[i, j].grid(True)
        ax[i, j].set_xlim(left=0, right=300)
        ax[i, j].set_ylim(bottom=-0, top=100)
        ax[1, j].set_xlabel('TS_RMW (km)', fontsize=25)
        
        ax[i, j].set_title('('+title_num[i][j]+')', loc='left', fontsize=25)
        ax[i, j].set_title(titles[i*2+j], fontsize=25)
        #ax[i, j].legend(fontsize=15)

fig.subplots_adjust(left=0.087, right=0.97, bottom=0.087, top=0.97)
plt.savefig('fig/LMIR_RMW_TS.png', dpi=600)
plt.show()

"""


#散佈圖
#genesis type
ew=0
mc=0
ms=0
mg=0
gen=['ew','mc','ms','md']
xy=['tsrmw','LMIR']
xy3=['tsrmw','LMIR']
for i in range(len(gen)):
    for j in range(len(xy)):
        exec(gen[i]+xy[j]+'=np.zeros(tstimestr.size)')
for tt in range(n_storm-start):
    if tdew[tt]==1:
        for j in range(len(xy)):
            exec(gen[0]+xy[j]+'['+str(tt)+']='+xy[j]+'['+str(tt)+']')
        if ew>=1:
            ax[1,1].scatter(tsrmw[tt], LMIR[tt],100, color='red', marker='s', alpha=0.2, linewidths=0)
            ew+=1
        else:
            ax[1,1].scatter(tsrmw[tt], LMIR[tt],100, color='red', marker='s', alpha=0.2, linewidths=0, label='EW')
            ew=1
            
    elif tdew[tt]==-1:
        for j in range(len(xy)):
            exec(gen[1]+xy[j]+'['+str(tt)+']='+xy[j]+'['+str(tt)+']')
        if mc>=1:
            ax[1,1].scatter(tsrmw[tt], LMIR[tt],100, color='purple', marker='o', alpha=0.2, linewidths=0)
            mc+=1
        else:
            ax[1,1].scatter(tsrmw[tt], LMIR[tt],100, color='purple', marker='o', alpha=0.2, linewidths=0, label='MC')
            mc=1
            
    elif tdew[tt]==-1.5:
        for j in range(len(xy)):
            exec(gen[2]+xy[j]+'['+str(tt)+']='+xy[j]+'['+str(tt)+']')
        if ms>=1:
            ax[1,1].scatter(tsrmw[tt], LMIR[tt],100, color='green', marker='o', alpha=0.2, linewidths=0)
            ms+=1
        else:
            ax[1,1].scatter(tsrmw[tt], LMIR[tt],100, color='green', marker='o', alpha=0.2, linewidths=0, label='MS')
            ms=1
            
    elif tdew[tt]==-2:
        for j in range(len(xy)):
            exec(gen[3]+xy[j]+'['+str(tt)+']='+xy[j]+'['+str(tt)+']')
        if mg>=1:
            ax[1,1].scatter(tsrmw[tt], LMIR[tt],100, color='blue', marker='v', alpha=0.2, linewidths=0)
            mg+=1
        else:
            ax[1,1].scatter(tsrmw[tt], LMIR[tt],100, color='blue', marker='v', alpha=0.2, linewidths=0, label='MD')
            mg=1
            

    else:
        ax[1,1].scatter(tsrmw[tt], LMIR[tt],100, color='black', marker='*', alpha=0.2, linewidths=1)
        


#calculate crossline
for i in range(len(gen)):
    for j in range(len(xy)):
        exec(gen[i]+xy[j]+'['+gen[i]+xy[j]+'==0]=np.nan')
        exec(gen[i]+xy[j]+'mean=np.nanmean('+gen[i]+xy[j]+')')
        exec(gen[i]+xy[j]+'std=np.nanstd('+gen[i]+xy[j]+')')


#draw crossline
for i in range(len(gen)):
    for j in range(len(xy)):
        if j==0:
            if i==0:
                exec('ax[1,1].plot(['+gen[i]+xy[j]+'mean-'+gen[i]+xy[j]+'std,'+gen[i]+xy[j]+'mean+'+gen[i]+xy[j]+'std],['+
                 gen[i]+xy[j+1]+'mean,'+gen[i]+xy[j+1]+'mean],color="red", linewidth=2)')
            if i==1:
                exec('ax[1,1].plot(['+gen[i]+xy[j]+'mean-'+gen[i]+xy[j]+'std,'+gen[i]+xy[j]+'mean+'+gen[i]+xy[j]+'std],['+
                 gen[i]+xy[j+1]+'mean,'+gen[i]+xy[j+1]+'mean],color="purple", linewidth=2)')
            if i==2:
                exec('ax[1,1].plot(['+gen[i]+xy[j]+'mean-'+gen[i]+xy[j]+'std,'+gen[i]+xy[j]+'mean+'+gen[i]+xy[j]+'std],['+
                 gen[i]+xy[j+1]+'mean,'+gen[i]+xy[j+1]+'mean],color="green", linewidth=2)')
            if i==3:
                exec('ax[1,1].plot(['+gen[i]+xy[j]+'mean-'+gen[i]+xy[j]+'std,'+gen[i]+xy[j]+'mean+'+gen[i]+xy[j]+'std],['+
                 gen[i]+xy[j+1]+'mean,'+gen[i]+xy[j+1]+'mean],color="blue", linewidth=2)')
            






# 調整x坐標為1/x的形式
inv_irmw = 1 / np.array(irmw)
inv_203 = 1 / 203

fig, ax = plt.subplots(1, 1, figsize=(6, 8))

ax.scatter(inv_irmw, LMIR, 100, color='red', alpha=0.2, linewidths=0)
#ax.scatter(inv_203, 95, 100, color='red', alpha=0.2, linewidths=0)
#for i in range(3):
#    ax.scatter(inv_203, 89, 100, color='red', alpha=0.2, linewidths=0)
#for i in range(5):
#    ax.scatter(inv_203, 83, 100, color='red', alpha=0.2, linewidths=0)

# 圖例，標題等
#ax.text(inv_203 + 0.005, 93, '1 case', fontsize=15)
#ax.text(inv_203 + 0.005, 87, '3 cases', fontsize=15)
#ax.text(inv_203 + 0.005, 81, '5 cases', fontsize=15)

# 更新反x值的線性圖
iline_inv_x = 1 / np.array(ilinex)
ax.plot(iline_inv_x, iliney, color='black')

# 圖例，標題等
#ax.text(1/155, 69, 'y=' + str(np.around(inum[0], decimals=4)) + '/x+' + str(np.around(inum[1], decimals=4)), fontsize=15)
ax.grid(True)

# 調整x軸限制以反映1/x的範圍
ax.set_xlim(left=0, right=1/12)
ax.set_ylim(bottom=-0, top=100)

ax.set_xlabel('1 / I_RMW (km)', fontsize=25)
ax.set_ylabel('$LMIR (kt \\, day^{-1})$', fontsize=25)

# 保存圖像
#plt.savefig('inverse_x_plot.png', dpi=600)

# 顯示圖像
plt.show()




#!!!draw VWS!!!

#filt NaN to regress
x=l_tsrmw[~np.isnan(l_tsrmw)]
y=l_LMIR[~np.isnan(l_LMIR)]
x=np.log(x)
l_tsnum,l_tscor=reg(x,y)
l_tslinex=np.arange(300)
l_tsliney=l_tsnum[0]*np.log(l_tslinex)+l_tsnum[1]#np.log

x=l_irmw[~np.isnan(l_irmw)]
y=l_LMIR[~np.isnan(l_LMIR)]
x=np.log(x)
l_inum,l_icor=reg(x,y)
l_ilinex=np.arange(300)
l_iliney=l_inum[0]*np.log(l_ilinex)+l_inum[1]#np.log

x=m_tsrmw[~np.isnan(m_tsrmw)]
y=m_LMIR[~np.isnan(m_LMIR)]
x=np.log(x)
m_tsnum,m_tscor=reg(x,y)
m_tslinex=np.arange(300)
m_tsliney=m_tsnum[0]*np.log(m_tslinex)+m_tsnum[1]#np.log

x=m_irmw[~np.isnan(m_irmw)]
y=m_LMIR[~np.isnan(m_LMIR)]
x=np.log(x)
m_inum,m_icor=reg(x,y)
m_ilinex=np.arange(300)
m_iliney=m_inum[0]*np.log(m_ilinex)+m_inum[1]#np.log

x=h_tsrmw[~np.isnan(h_tsrmw)]
y=h_LMIR[~np.isnan(h_LMIR)]
x=np.log(x)
h_tsnum,h_tscor=reg(x,y)
h_tslinex=np.arange(300)
h_tsliney=h_tsnum[0]*np.log(h_tslinex)+h_tsnum[1]#np.log

x=h_irmw[~np.isnan(h_irmw)]
y=h_LMIR[~np.isnan(h_LMIR)]
x=np.log(x)
h_inum,h_icor=reg(x,y)
h_ilinex=np.arange(300)
h_iliney=h_inum[0]*np.log(h_ilinex)+h_inum[1]#np.log


# 製作figure  
fig1 = plt.figure()   

#圖表的設定
ax1 = fig1.add_subplot(1, 1, 1)

#散佈圖
ax1.scatter(l_tsrmw, l_LMIR, color='purple', alpha=0.2, linewidths=0)
ax1.plot(l_tslinex,l_tsliney, color='purple')
ax1.scatter(m_tsrmw, m_LMIR, color='blue', alpha=0.2, linewidths=0)
ax1.plot(m_tslinex,m_tsliney, color='blue')
ax1.scatter(h_tsrmw, h_LMIR, color='green', alpha=0.2, linewidths=0)
ax1.plot(h_tslinex,h_tsliney, color='green')

#圖例，標題等
#ax1.text(155,90,'y='+str(np.around(tsnum[0],decimals=4))+'ln(x)+'+str(np.around(tsnum[1],decimals=4)))
#ax1.text(155,85,'correlation='+str(np.around(tscor,decimals=4)))
ax1.legend(['VWS:<4.5m/s','VWS:4.5-10m/s','VWS:>10m/s'])
#ax1.legend(['VWS:>10m/s'])
ax1.grid(True)
ax1.set_xlim(left=0, right=300)
ax1.set_ylim(bottom=-0, top=100)
#ax1.set_xscale("symlog")
#ax1.set_yscale("symlog")
ax1.set_xlabel('RMW of TS formation (km)')
ax1.set_ylabel('LMIR (kt/24h)')
ax1.set_title('RMW of TS formation vs. Lifetime Max. Intensification Rate')
plt.savefig('fig/LMIR_tsrmw_shear.png', dpi=300)


# 製作figure  
fig2 = plt.figure()   

#圖表的設定
ax2 = fig2.add_subplot(1, 1, 1)

#散佈圖
ax2.scatter(l_irmw, l_LMIR, color='purple', alpha=0.2, linewidths=0)
ax2.plot(l_ilinex,l_iliney, color='purple')
ax2.scatter(m_irmw, m_LMIR, color='blue', alpha=0.2, linewidths=0)
ax2.plot(m_ilinex,m_iliney, color='blue')
ax2.scatter(h_irmw, h_LMIR, color='green', alpha=0.2, linewidths=0)
ax2.plot(h_ilinex,h_iliney, color='green')

#圖例，標題等
#ax2.text(155,90,'y='+str(np.around(inum[0],decimals=4))+'ln(x)+'+str(np.around(inum[1],decimals=4)))
#ax2.text(155,85,'correlation='+str(np.around(icor,decimals=4)))
ax2.legend(['VWS:<4.5m/s','VWS:4.5-10m/s','VWS:>10m/s'])
#ax2.legend(['VWS:>10m/s'])
ax2.grid(True)
ax2.set_xlim(left=0, right=300)
ax2.set_ylim(bottom=-0, top=100)
#ax2.set_xscale("symlog")
#ax2.set_yscale("symlog")
ax2.set_xlabel('RMW of Max. Intensification (km)')
ax2.set_ylabel('LMIR (kt/24h)')
ax2.set_title('RMW of Max. Intensification vs. Lifetime Max. Intensification Rate')
plt.savefig('fig/LMIR_irmw_shear.png', dpi=300)
plt.show() 
"""
