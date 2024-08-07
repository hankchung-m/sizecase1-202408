#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 08:13:34 2022

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
tdlat=np.load('tdlat.npy')
tdlon=np.load('tdlon.npy')

itimestr=np.load('itimestr.npy', allow_pickle=True)
ihour=np.load('ihour.npy', allow_pickle=True)
tstimestr=np.load('tstimestr.npy', allow_pickle=True)
tshour=np.load('tshour.npy', allow_pickle=True)
tdtimestr=np.load('tdtimestr.npy', allow_pickle=True)
tdhour=np.load('tdhour.npy', allow_pickle=True)

# 0:itime 1:tstime 2:tdtime !!!choose
for ii in range(3):
    if ii==1:
        itimestr=tstimestr
        ihour=tshour
        ilat=tslat
        ilon=tslon
    elif ii==2:
        itimestr=tdtimestr
        ihour=tdhour
        ilat=tdlat
        ilon=tdlon
    
    interxy=8#8degree(+-4degree)
    
    ilon_nx=np.zeros(itimestr.size)
    ilat_ny=np.zeros(itimestr.size)
    tslon_nx=np.zeros(itimestr.size)
    tslat_ny=np.zeros(itimestr.size)
    tdlon_nx=np.zeros(itimestr.size)
    tdlat_ny=np.zeros(itimestr.size)
    
    tsvt=np.zeros((itimestr.size,40)) # !!!(+-1000km)
    tsvr=np.zeros((itimestr.size,40))
    tsvo=np.zeros((itimestr.size,40))
    tsqv=np.zeros((itimestr.size,40))
    tsrh=np.zeros((itimestr.size,40))
    tdrmw_n=np.zeros(itimestr.size)
    tsrmw_n=np.zeros(itimestr.size)
    tsvtm=np.zeros(itimestr.size)
    
    for tt in range(itimestr.size):
        
        #print(tt)
        
        #sort real timestr
        if testdt(itimestr[tt]) == False :
            continue
        
        a=ihour[tt]
        if a != '00' and a != '03' and a != '06' and a != '09' and a != '12' and a != '15' and a != '18' and a != '21':
            continue
        
        if ii==0:
            if tt/10<1 :
                ncfile = nc.Dataset("ERA5_I00%d.nc"%(tt))
            elif tt/100<1 :
                ncfile = nc.Dataset("ERA5_I0%d.nc"%(tt))
            else:
                ncfile = nc.Dataset("ERA5_I%d.nc"%(tt))
        elif ii==1 :
            if tt/10<1 :
                ncfile = nc.Dataset("ERA5_TS00%d.nc"%(tt))
            elif tt/100<1 :
                ncfile = nc.Dataset("ERA5_TS0%d.nc"%(tt))
            else:
                ncfile = nc.Dataset("ERA5_TS%d.nc"%(tt))
        elif ii==2 :
            if tt/10<1 :
                ncfile = nc.Dataset("ERA5_TD00%d.nc"%(tt))
            elif tt/100<1 :
                ncfile = nc.Dataset("ERA5_TD0%d.nc"%(tt))
            else:
                ncfile = nc.Dataset("ERA5_TD%d.nc"%(tt))
                
        XLAT = ncfile["latitude"]
        XLONG = ncfile["longitude"]
        p = ncfile["level"]
        u = ncfile["u"]
        v = ncfile["v"]
        z = ncfile["z"]
        qv = ncfile["q"]
        rh = ncfile["r"]
        vo = ncfile["vo"]
        #z=z/9.8
        
        
        #find low & left
        nz=np.zeros((interxy*4+1,interxy*4+1))
        low=0
        left=0
        for i in range(XLONG.size) :
            for j in range(XLAT.size) :
                lon1=ilon[tt]-(interxy/2)
                lat1=ilat[tt]-(interxy/2)
                # abs(lat1-XLAT[j]):0-0.25/2 , abs(lon1-XLONG[i]):0-0.25/2
                if abs(lat1-XLAT[j])<0.25/2 and abs(lat1-XLAT[j])>=0:
                    low=j
                if abs(lon1-XLONG[i])<0.25/2 and abs(lon1-XLONG[i])>=0:
                    left=i
        """
        for i in range(left,left+(interxy*4)+1,1) :
            for j in range(low,low-(interxy*4)-1,-1) : #j is 0 at top
                if abs(z[0,30,j,i])<7000 :# and sp[0,j,i]>85000:
                    nz[low-j,i-left]=z[0,30,j,i]
                else:
                    nz[low-j,i-left]=0
                    break
        """
        nz=z[0,30,low-(interxy*4)-1:low,left:left+(interxy*4)+1]
        
        #Find center
    
        # step 1
        x1 = np.where(nz==np.nanmin(nz))[1][0]
        y1 = np.where(nz==np.nanmin(nz))[0][0]
        
        # step 2 and 3
        XLAT_tmp=np.zeros((interxy*4+1,interxy*4+1))
        XLONG_tmp=np.zeros((interxy*4+1,interxy*4+1))
        for i in range((interxy*4)+1):
            XLAT_tmp[:,i]=XLAT[low-(interxy*4)-1:low]
            XLONG_tmp[i,:]=XLONG[left:left+(interxy*4)+1]
        SLP_anom=nz-np.nanmean(nz)
        SLP_anom[SLP_anom>0]=0
        SLP_anom=-SLP_anom
        # step 4
        lon=round(np.average(XLONG_tmp,weights=SLP_anom),2)
        lat=round(np.average(XLAT_tmp,weights=SLP_anom),2)
        # Domain x and y position
        x2=np.where(abs(XLONG[:]-lon)<=0.125)[0]
        y2=np.where(abs(XLAT[:]-lat)<=0.125)[0]
        
        
        
        if ii==0:
            ilon_nx[tt]=x2
            ilat_ny[tt]=y2
            dd=Distance(XLAT[y2],XLONG[x2],ilat[tt],ilon[tt])
            if dd>100:
                print(tt,'relocate>100 km,','Distance=',dd)
        elif ii==1:
            tslon_nx[tt]=x2
            tslat_ny[tt]=y2
            dd=Distance(XLAT[y2],XLONG[x2],tslat[tt],tslon[tt])
            if dd>100:
                print(tt,'relocate>100 km,','Distance=',dd)
        elif ii==2:
            tdlon_nx[tt]=x2
            tdlat_ny[tt]=y2
            dd=Distance(XLAT[y2],XLONG[x2],tdlat[tt],tdlon[tt])
            if dd>100:
                print(tt,'relocate>100 km,','Distance=',dd)
            
        if ii!=0:
            # Vr and Vt
            #  建立Δx和Δy矩陣 !!!(+-10 degree)
            distance_x, distance_y = np.meshgrid(np.arange(-40,41),np.arange(-40,41))
            distance=np.zeros_like((distance_x))
            for i in range(-40,41):
                for j in range(-40,41):
                    if i<0:
                        distance_x[j+40,i+40]=-Distance(XLAT[y2],XLONG[x2],XLAT[y2],XLONG[x2]+(i/4))
                    else:
                        distance_x[j+40,i+40]=Distance(XLAT[y2],XLONG[x2],XLAT[y2],XLONG[x2]+(i/4))
                    if j<0:
                        distance_y[j+40,i+40]=-Distance(XLAT[y2],XLONG[x2],XLAT[y2]+(j/4),XLONG[x2])
                    else:
                        distance_y[j+40,i+40]=Distance(XLAT[y2],XLONG[x2],XLAT[y2]+(j/4),XLONG[x2])
                    distance[j+40,i+40]=Distance(XLAT[y2],XLONG[x2],XLAT[y2]+(j/4),XLONG[x2]+(i/4))
            
            #  建立與中心距離的矩陣
            #distance = (distance_x**2+distance_y**2)**0.5
    
            VT=np.zeros((81,81))
            VR=np.zeros((81,81))
            VO=np.zeros((81,81))
            QV=np.zeros((81,81))
            RH=np.zeros((81,81))
            ss=int(y2+40) #j is 0 at top
            nn=int(y2-41)
            ww=int(x2-40)
            ee=int(x2+41)
            outerxy=20    #(+-10 degree)
            for i in range((outerxy*4)+1) :
                for j in range((outerxy*4)+1) : #j is 0 at top
                    erai=ww+i
                    eraj=ss-j
                    VT[j,i] = -u[0,35,eraj,erai]*distance_y[j,i]/distance[j,i]+v[0,35,eraj,erai]*distance_x[j,i]/distance[j,i]
                    VR[j,i] = u[0,35,eraj,erai]*distance_x[j,i]/distance[j,i]+v[0,35,eraj,erai]*distance_y[j,i]/distance[j,i]
                    VO[j,i] = vo[0,35,eraj,erai]
                    QV[j,i] = qv[0,30,eraj,erai]
                    RH[j,i] = rh[0,30,eraj,erai]
                    
            
            #print(VT[45,45])
            #print(RH[45,45])
    
            #  中心處設為0
            VT[40,40] = 0
            VR[40,40] = 0
            
            # Azimuthal mean
            #distance=np.round(distance)
            VT_azi=np.zeros(40)
            VR_azi=np.zeros(40)
            VO_azi=np.zeros(40)
            QV_azi=np.zeros(40)
            RH_azi=np.zeros(40)
            ndistance=np.zeros((40,81,81))
            
            first=0
            for r in range(40):
                for i in range(-40,41):
                    for j in range(-40,41):
                        if distance[j+40,i+40]>=25*(r-0.9) and distance[j+40,i+40]<=25*(r+0.9):
                            ndistance[r,j+40,i+40]=r 
                VT_azi[r] = np.mean(VT[ndistance[r]==r])
                VR_azi[r] = np.mean(VR[ndistance[r]==r])
                VO_azi[r] = np.mean(VO[ndistance[r]==r])
                QV_azi[r] = np.mean(QV[ndistance[r]==r])
                RH_azi[r] = np.mean(RH[ndistance[r]==r])
                #print(ndistance[r])
                if r>1 and VT_azi[r]>0 and first==0: #let first Vt max is rmw
                    inner=VT_azi[r-1]-VT_azi[r-2]
                    outer=VT_azi[r-1]-VT_azi[r]
                    if inner>0 and outer>0:
                        tsrmw_n[tt]=(r-1)*25
                        tsvtm[tt]=VT_azi[r-1]
                        first=1
                if r==40 and first==0: #let rmw>1000km is 1000km
                    if np.max(VT_azi)>0:
                        tsrmw_n[tt]=np.argmax(VT_azi)*25
                        tsvtm[tt]=np.max(VT_azi)
    
                    #tdrmw_n[tt]=r*0.25
                    #print(tt,'RMW>10 degree,','Max RMW=',0.25*np.argmax(VT_azi),'Max Vt=',max(VT_azi))
            
            tsvt[tt]=VT_azi
            tsvr[tt]=VR_azi
            tsvo[tt]=VO_azi
            tsqv[tt]=QV_azi
            tsrh[tt]=RH_azi
            
            #print(tsvt[tt])
            
            if tdrmw_n[tt]==0:
                tdrmw_n[tt]=np.nan
            if tsrmw_n[tt]==0:
                tsrmw_n[tt]=np.nan
                tsvtm[tt]=np.nan
                tsvt[tt]=np.nan
       
            
    if ii==0:
        np.save('ilon_nx',ilon_nx)
        np.save('ilat_ny',ilat_ny)
    elif ii==1:
        np.save('tslon_nx',tslon_nx)
        np.save('tslat_ny',tslat_ny)
        np.save('tsvt_azi',tsvt)
        np.save('tsvr_azi',tsvr)
        np.save('tsvo_azi',tsvo)
        np.save('tsqv_azi',tsqv)
        np.save('tsrh_azi',tsrh)
        np.save('tsvtm',tsvtm)
        np.save('tsrmw_n',tsrmw_n)
    elif ii==2:
        np.save('tdlon_nx',tdlon_nx)
        np.save('tdlat_ny',tdlat_ny)
        np.save('tdvt_azi',tsvt)
        np.save('tdvr_azi',tsvr)
        np.save('tdvo_azi',tsvo)
        np.save('tdqv_azi',tsqv)
        np.save('tdrh_azi',tsrh)
        np.save('tdvtm',tsvtm)
        np.save('tdrmw_n',tsrmw_n)