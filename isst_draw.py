#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 18:09:40 2022

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

#!!!
#6:SST
var=6

#!!!
#1:i
#2:ts
#3:i24
drawt=1

LMIR=np.load('LMIR.npy')
if drawt==1:
    tslat=np.load('ilat.npy')
    tslon=np.load('ilon.npy')
    tstimestr=np.load('tdtimestr.npy', allow_pickle=True)
    tshour=np.load('tdhour.npy', allow_pickle=True)
if drawt==2:
    tslat=np.load('ilat_24.npy')
    tslon=np.load('ilon_24.npy')
    tstimestr=np.load('tstimestr.npy', allow_pickle=True)
    tshour=np.load('tshour.npy', allow_pickle=True)
if drawt==3:
    tslat=np.load('ilat_24.npy')
    tslon=np.load('ilon_24.npy')
    tstimestr=np.load('tdtimestr.npy', allow_pickle=True)
    tshour=np.load('tdhour.npy', allow_pickle=True)


interxy=10 #!!!10degree(+-5degree)
f1=np.zeros((int(tstimestr.size),interxy*4+1,interxy*4+1))
f2=np.zeros((int(tstimestr.size),interxy*4+1,interxy*4+1))
f3=np.zeros((int(tstimestr.size),interxy*4+1,interxy*4+1))
f4=np.zeros((int(tstimestr.size),interxy*4+1,interxy*4+1))
f1com=np.zeros((10,2,interxy*4+1,interxy*4+1)) #(30-75 , h/l , interxy*4+1,interxy*4+1)
fill=np.zeros((10,2,interxy*4+1,interxy*4+1))

title_num=[['a','b'],['c','d'],['e','f'],['g','h']]

for tt in range(tstimestr.size):
    
    #print(tt)
    
    #sort real timestr
    if testdt(tstimestr[tt]) == False :
        continue
    
    a=tshour[tt]
    if a != '00' and a != '03' and a != '06' and a != '09' and a != '12' and a != '15' and a != '18' and a != '21':
        continue
    
    if np.isnan(tslat[tt])==True or np.isnan(tslon[tt])==True:
        continue
    
    
    if drawt==1:
        if tt/10<1 :
            ncfile = nc.Dataset("ERA5_SFC_TD00%d.nc"%(tt))
        elif tt/100<1 :
            ncfile = nc.Dataset("ERA5_SFC_TD0%d.nc"%(tt))
        else:
            ncfile = nc.Dataset("ERA5_SFC_TD%d.nc"%(tt))
    if drawt==2:
        if tt/10<1 :
            ncfile = nc.Dataset("ERA5_SFC_TS00%d.nc"%(tt))
        elif tt/100<1 :
            ncfile = nc.Dataset("ERA5_SFC_TS0%d.nc"%(tt))
        else:
            ncfile = nc.Dataset("ERA5_SFC_TS%d.nc"%(tt))
    if drawt==3:
        if tt/10<1 :
            ncfile = nc.Dataset("ERA5_SFC_TD00%d.nc"%(tt))
        elif tt/100<1 :
            ncfile = nc.Dataset("ERA5_SFC_TD0%d.nc"%(tt))
        else:
            ncfile = nc.Dataset("ERA5_SFC_TD%d.nc"%(tt))
        
    XLAT = ncfile["latitude"]
    XLONG = ncfile["longitude"]
    #p = ncfile["level"]
    #u = ncfile["u"]#(1,37,321,521)
    #v = ncfile["v"]
    #w = ncfile["w"]
    #z = ncfile["z"]
    #r = ncfile["r"]
    #vo = ncfile["vo"]
    #sp = ncfile["sp"]
    sst = ncfile["sst"]
    #print(v.shape)
    
    if var==1:
        f11=sst[0,:,:]
        f22=sst[0,:,:]
        f33=sst[0,:,:]
    if var==2:
        f11=sst[0,:,:]
        f22=sst[0,:,:]
        f33=sst[0,:,:]
    if var==3:
        f11=sst[0,:,:]
        f22=sst[0,:,:]
        f33=sst[0,:,:]
    if var==4:
        f11=sst[0,:,:]
        f22=sst[0,:,:]
        f33=sst[0,:,:]
    if var==5:
        f11=sst[0,:,:]
        f22=sst[0,:,:]
        f33=sst[0,:,:]
    if var==6:
        f11=sst[0,:,:]
        #f22=sst[0,:,:]
        #f33=sst[0,:,:]
        
    #f11=u[0,30,:,:]
    #f22=v[0,30,:,:]
    #f33=z[0,21,:,:]/9.8
    #f11=r[0,25,:,:]
    #f22=w[0,21,:,:]
    #f33=z[0,25,:,:]/9.8
    #f33=vo[0,30,:,:]
    
    
    #find low & left
    low=0
    left=0
    for i in range(XLONG.size) :
        for j in range(XLAT.size) :
            lon1=tslon[tt]-(interxy/2)
            lat1=tslat[tt]-(interxy/2)
            # abs(lat1-XLAT[j]):0-0.25/2 , abs(lon1-XLONG[i]):0-0.25/2
            if abs(lat1-XLAT[j])<0.25/2 and abs(lat1-XLAT[j])>=0:
                low=j
            if abs(lon1-XLONG[i])<0.25/2 and abs(lon1-XLONG[i])>=0:
                left=i
                    
    for lev in range(30,80,5):#30-75
        levi=int((lev-30)/5)
        
        #no longer interp
        for i in range(left,left+(interxy*4)+1,1) :
            for j in range(low,low-(interxy*4)-1,-1) : #j is 0 at top
                
                #nlon=tslon[tt]-(interxy/2)+(0.25*(i-left))
                #nlat=tslat[tt]-(interxy/2)+(0.25*(low-j))
                #b=(XLONG[i+1]-XLONG[i])*(XLAT[j-1]-XLAT[j])
                    
                #sort missing value
                #f22range=400
                if var==1 or var==2:
                    f11range=400
                    #f33range=400
                if var==3:
                    f11range=400
                    #f33range=400
                if var==4 or var==5 or var==6:
                    f11range=400
                    #f33range=400
                if abs(f11[j,i])<f11range and f11[j,i]>-1:
                    
                    if LMIR[tt] >= lev :
                        fill[levi,0,low-j,i-left]+=1
                        f1com[levi,0,low-j,i-left]+=f11[j,i]-273.15
                    elif LMIR[tt] < lev :
                        fill[levi,1,low-j,i-left]+=1
                        f1com[levi,1,low-j,i-left]+=f11[j,i]-273.15
                    
                    #f1[tt,low-j,i-left]=f11[j,i]-273.15
                    #f2[tt,low-j,i-left]=f22[j,i]-273.15
                    #f3[tt,low-j,i-left]=f33[j,i]-273.15
                
                    #f1[tt,low-j,i-left]=interp(f11[j,i],f11[j,i+1],f11[j-1,i],f11[j-1,i+1],XLONG[i],XLONG[i+1],nlon,XLAT[j],XLAT[j-1],nlat)/b
                    #f2[tt,low-j,i-left]=interp(f22[j,i],f22[j,i+1],f22[j-1,i],f22[j-1,i+1],XLONG[i],XLONG[i+1],nlon,XLAT[j],XLAT[j-1],nlat)/b
                    #f3[tt,low-j,i-left]=interp(f33[j,i],f33[j,i+1],f33[j-1,i],f33[j-1,i+1],XLONG[i],XLONG[i+1],nlon,XLAT[j],XLAT[j-1],nlat)/b
                #else:
                #    f1[tt,interxy*2,interxy*2]=0
                #    f2[tt,interxy*2,interxy*2]=0
                #    f3[tt,interxy*2,interxy*2]=0
                    #print(f11[j,i],f22[j,i],f33[j,i])
                #    break
                    
        #f1=f11[low-(interxy*4)-1:low,left:left+(interxy*4)+1]
        #f2=f22[low-(interxy*4)-1:low,left:left+(interxy*4)+1]
        #f3=f33[low-(interxy*4)-1:low,left:left+(interxy*4)+1]
                
        
        #if f1[tt,interxy*2,interxy*2]==0 or np.isnan(f1[tt,interxy*2,interxy*2])==True or f2[tt,interxy*2,interxy*2]==0 or np.isnan(f2[tt,interxy*2,interxy*2])==True or f3[tt,interxy*2,interxy*2]==0 or np.isnan(f3[tt,interxy*2,interxy*2])==True :
        #    f1[tt,:,:]=np.nan
        #    f2[tt,:,:]=np.nan
        #    f3[tt,:,:]=np.nan
            
        #print(f1[tt,:,:])

for lev in range(30,80,5):#30-75
    levi=int((lev-30)/5)
    for i in range(left,left+(interxy*4)+1,1) :
        for j in range(low,low-(interxy*4)-1,-1) : #j is 0 at top
            for hl in range(2):
                f1com[levi,hl,low-j,i-left]=f1com[levi,hl,low-j,i-left]/fill[levi,hl,low-j,i-left]

#f1com=composite(f1,interxy)
#f2com=composite(f2,interxy)
#f3com=composite(f3,interxy)

# npy
np.save('f1com_sst',f1com)
#np.save('f2com',f2com)
#np.save('f3com',f3com)

"""

f1com=np.load('f1com_sst.npy')
#f2com=np.load('f2com.npy')
#f3com=np.load('f3com.npy')

"""

#print(f3com)
if var==1:
    title1='850hPa wind, 500hPa geopotential height, '
    title3='850w_500z_'
if var==2:
    title1='850hPa vorticity, 850hPa wind, '
    title3='850v_850w_'
if var==3:
    title1='700hPa RH, 700hPa geopotential height, '
    title3='700r_700z_'
if var==4:
    title1='200hPa wind, '
    title3='200w_200z_'
if var==5:
    title1='200-850hPa VWS, '
    title3='200_850vws_'
if var==6:
    title1='SST, '
    title3='sst_'

NULAT=np.arange(-int(interxy/2),int(interxy/2)+0.25,0.25)
NULONG=np.arange(-int(interxy/2),int(interxy/2)+0.25,0.25)

#draw https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.streamplot.html
fig,ax=plt.subplots(4,2,sharex=True,sharey=True,figsize=(18,32))
for lev in range(30,70,10):#30-75
    levi=int((lev-30)/10)
    for hl in range(2):
        #levi=int((lev-30)/5)

        ax[levi,hl].axis([-5,5,-5,5])
        ax[levi,hl].tick_params(labelsize=25)
        ax[levi,hl].grid(True)
        ax[levi,hl].set_aspect('equal')

        #plt.axes(projection=ccrs.PlateCarree()).coastlines(resolution='auto', color='k',linewidth=0.5)
        
        if var==6:
            c=ax[levi,hl].contourf(NULONG[:],NULAT[:],f1com[levi,hl,:,:],cmap="YlOrRd",levels=np.arange(28.5,30.1,0.1),extent=(-5,-5,5,5),extend='both')
            e=ax[levi,hl].contour(NULONG[:],NULAT[:],f1com[levi,hl,:,:],levels=np.arange(29,32,0.5),color="b",linewidths=1)
            ax[levi,hl].clabel(e, inline=True, fontsize=25)
        if hl==0:
            ax[levi,hl].set_title(title1+'>='+str(lev)+'kt/24h', fontsize=25)
        else:
            ax[levi,hl].set_title(title1+'<'+str(lev)+'kt/24h', fontsize=25)
        ax[levi,hl].set_title('('+title_num[levi][hl]+')',loc='left', fontsize=25)
        ax[3,hl].set_xlabel('$(^o)$', fontsize=25)
        ax[levi,0].set_ylabel('$(^o)$', fontsize=25)
        

cbar_ax = fig.add_axes([0.82, 0.15, 0.02, 0.7]) #x,y,dx,dy
cbar=fig.colorbar(c, cax=cbar_ax)
cbar.ax.tick_params(labelsize=25)
fig.subplots_adjust(right=0.8)
fig.subplots_adjust(right=0.8)

if drawt==1:
    fig.savefig('fig/td_i_'+title3+'.png', dpi=300)
if drawt==2:
    fig.savefig('fig/ts_'+title3+'.png', dpi=300)
if drawt==3:
    fig.savefig('fig/td_i24_'+title3+'.png', dpi=300)