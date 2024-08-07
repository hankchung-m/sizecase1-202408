#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 17:59:10 2022

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
#1:850hPa wind, 500hPa geopotential height
#2:850hPa vorticity, 850hPa wind
#3:700hPa RH, 700hPa geopotential height
#4:200hPa wind, 200hPa geopotential height
#5:200-850hPa VWS
#6:500hPa RH, 500hPa geopotential height
var=6

#!!!
#1:i
#2:i24
drawt=1

LMIR=np.load('LMIR.npy')

if drawt==1:
    tslat=np.load('ilat.npy')
    tslon=np.load('ilon.npy')
if drawt==2:
    tslat=np.load('ilat_24.npy')
    tslon=np.load('ilon_24.npy')
    
tstimestr=np.load('itimestr.npy', allow_pickle=True)
tshour=np.load('ihour.npy', allow_pickle=True)

if var==3 or var==6:
    interxy=20
else:
    interxy=40 #!!!40degree(+-20degree)
f1=np.zeros((int(tstimestr.size),interxy*4+1,interxy*4+1))
f2=np.zeros((int(tstimestr.size),interxy*4+1,interxy*4+1))
f3=np.zeros((int(tstimestr.size),interxy*4+1,interxy*4+1))
f4=np.zeros((int(tstimestr.size),interxy*4+1,interxy*4+1))

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
    
    
    if tt/10<1 :
        ncfile = nc.Dataset("ERA5_I00%d.nc"%(tt))
    elif tt/100<1 :
        ncfile = nc.Dataset("ERA5_I0%d.nc"%(tt))
    else:
        ncfile = nc.Dataset("ERA5_I%d.nc"%(tt))
        
    XLAT = ncfile["latitude"]
    XLONG = ncfile["longitude"]
    p = ncfile["level"]
    u = ncfile["u"]#(1,37,321,521)
    v = ncfile["v"]
    w = ncfile["w"]
    z = ncfile["z"]
    r = ncfile["r"]
    vo = ncfile["vo"]
    #sp = ncfile["sp"]
    #print(v.shape)
    
    if var==1:
        f11=u[0,30,:,:]
        f22=v[0,30,:,:]
        f33=z[0,21,:,:]/9.8
    if var==2:
        f11=u[0,30,:,:]
        f22=v[0,30,:,:]
        f33=vo[0,30,:,:]
    if var==3:
        f11=r[0,25,:,:]
        f22=w[0,21,:,:]
        f33=z[0,25,:,:]/9.8
    if var==4:
        f11=u[0,14,:,:]
        f22=v[0,14,:,:]
        f33=z[0,14,:,:]/9.8
    if var==5:
        f11=u[0,14,:,:]-u[0,30,:,:]
        f22=v[0,14,:,:]-v[0,30,:,:]
        f33=vo[0,14,:,:]
    if var==6:
        f11=r[0,21,:,:]
        f22=w[0,21,:,:]
        f33=z[0,21,:,:]/9.8
        
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
                
    lack=0
    
    #no longer interp
    for i in range(left,left+(interxy*4)+1,1) :
        for j in range(low,low-(interxy*4)-1,-1) : #j is 0 at top
            
            #nlon=tslon[tt]-(interxy/2)+(0.25*(i-left))
            #nlat=tslat[tt]-(interxy/2)+(0.25*(low-j))
            #b=(XLONG[i+1]-XLONG[i])*(XLAT[j-1]-XLAT[j])
                
            #sort missing value
            f22range=100
            if var==1 or var==2:
                f11range=100
                f33range=7000
            if var==3:
                f11range=200
                f33range=7000
            if var==4 or var==5:
                f11range=100
                f33range=13000
            if var==6:
                f11range=200
                f33range=7000
            if abs(f11[j,i])<f11range and abs(f22[j,i])<f22range and abs(f33[j,i])<f33range and f33[j,i]>-1:
                f1[tt,low-j,i-left]=f11[j,i]
                f2[tt,low-j,i-left]=f22[j,i]
                f3[tt,low-j,i-left]=f33[j,i]
            
                #f1[tt,low-j,i-left]=interp(f11[j,i],f11[j,i+1],f11[j-1,i],f11[j-1,i+1],XLONG[i],XLONG[i+1],nlon,XLAT[j],XLAT[j-1],nlat)/b
                #f2[tt,low-j,i-left]=interp(f22[j,i],f22[j,i+1],f22[j-1,i],f22[j-1,i+1],XLONG[i],XLONG[i+1],nlon,XLAT[j],XLAT[j-1],nlat)/b
                #f3[tt,low-j,i-left]=interp(f33[j,i],f33[j,i+1],f33[j-1,i],f33[j-1,i+1],XLONG[i],XLONG[i+1],nlon,XLAT[j],XLAT[j-1],nlat)/b
            else:
                lack=1
                #print(f11[j,i],f22[j,i],f33[j,i])
                break
        
        if lack==1:
            break
    #f1=f11[low-(interxy*4)-1:low,left:left+(interxy*4)+1]
    #f2=f22[low-(interxy*4)-1:low,left:left+(interxy*4)+1]
    #f3=f33[low-(interxy*4)-1:low,left:left+(interxy*4)+1]
            
    
    #if f1[tt,interxy*2,interxy*2]==0 or np.isnan(f1[tt,interxy*2,interxy*2])==True or f2[tt,interxy*2,interxy*2]==0 or np.isnan(f2[tt,interxy*2,interxy*2])==True or f3[tt,interxy*2,interxy*2]==0 or np.isnan(f3[tt,interxy*2,interxy*2])==True :
    if lack==1:
        f1[tt,:,:]=np.nan
        f2[tt,:,:]=np.nan
        f3[tt,:,:]=np.nan
        
    #print(f1[tt,:,:])

for tt in range(tstimestr.size):
    if f1[tt,interxy*2,interxy*2]==0 or f2[tt,interxy*2,interxy*2]==0 or f3[tt,interxy*2,interxy*2]==0:
        f1[tt,:,:]=np.nan
        f2[tt,:,:]=np.nan
        f3[tt,:,:]=np.nan

f1com=composite(f1,interxy)
#f2com=composite(f2,interxy)
#f3com=composite(f3,interxy)


# npy
np.save('f1com_'+str(var),f1com)
#np.save('f2com',f2com)
#np.save('f3com',f3com)

"""

f1com=np.load('f1com_'+str(var)+'.npy')
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
    title1='700hPa RH, '
    title3='700r_700z_'
if var==4:
    title1='200hPa wind, '
    title3='200w_200z_'
if var==5:
    title1='200-850hPa VWS, '
    title3='200_850vws_'
if var==6:
    title1='500hPa RH, '
    title3='500r_500z_'

NULAT=np.arange(-int(interxy/2),int(interxy/2)+0.25,0.25)
NULONG=np.arange(-int(interxy/2),int(interxy/2)+0.25,0.25)

#draw https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.streamplot.html
fig,ax=plt.subplots(4,2,sharex=True,sharey=True,figsize=(18,32))
for lev in range(30,70,10):#30-75
    levi=int((lev-30)/10)
    for hl in range(2):
        #levi=int((lev-30)/5)
        
        ax[levi,hl].axis([-10,10,-10,10])
        ax[levi,hl].set_xticks(np.arange(-10,15,5))
        ax[levi,hl].set_yticks(np.arange(-10,15,5))
        ax[levi,hl].tick_params(labelsize=25)
        ax[levi,hl].grid(True)
        ax[levi,hl].set_aspect('equal')

        #plt.axes(projection=ccrs.PlateCarree()).coastlines(resolution='auto', color='k',linewidth=0.5)
        
        if var==1 or var==2 or var==4 or var==5:
            f12com=np.zeros((interxy*4+1,interxy*4+1))
            for i in range(interxy*4+1) :
                for j in range(interxy*4+1) :
                    f12com[j,i]=sqrt(f1com[levi,hl,j,i]**2+f2com[levi,hl,j,i]**2)

        if var==1:
            c=ax[hl].contourf(NULONG[:],NULAT[:],f12com[:,:],cmap="Reds",levels=np.arange(0,15,2),extent=(-20,-20,20,20),extend='max')
            d=ax[hl].streamplot(NULONG[:],NULAT[:],f1com[levi,hl,:,:],f2com[hl,:,:],color="b")
            e=ax[hl].contour(NULONG[:],NULAT[:],f3com[hl,:,:],levels=np.arange(5280,6480,15),color="g",linewidths=5)
            ax[hl].clabel(e, inline=True, fontsize=15)
        
        if var==2:
            c=ax[hl].contourf(NULONG[:],NULAT[:],f3com[levi,hl,:,:],cmap="bwr",levels=np.arange(-0.000015,0.000018,0.000003),extent=(-20,-20,20,20),extend='both')
            d=ax[hl].streamplot(NULONG[:],NULAT[:],f1com[levi,hl,:,:],f2com[hl,:,:],color="b")
            e=ax[hl].contour(NULONG[:],NULAT[:],f12com[:,:],levels=np.arange(0,16,4),color="yellow",linewidths=3)
            ax[hl].clabel(e, inline=True, fontsize=15)
        
        if var==3:
            c=ax[levi,hl].contourf(NULONG[:],NULAT[:],f1com[levi,hl,:,:],cmap="Blues",levels=np.arange(50,100,10),extent=(-20,-20,20,20),extend='both')
            #d=ax[hl].contour(NULONG[:],NULAT[:],f2com[hl,:,:],levels=np.arange(-0.3,-0.1,0.1),color="red",linewidths=1)
            #e=ax[hl].contour(NULONG[:],NULAT[:],f3com[levi,hl,:,:],levels=np.arange(2820,3420,30),color="g",linewidths=5)
            #ax[hl].clabel(d, inline=True, fontsize=10)
            #ax[hl].clabel(e, inline=True, fontsize=25)
        
        if var==4:
            c=ax[hl].contourf(NULONG[:],NULAT[:],f12com[:,:],cmap="Reds",levels=np.arange(0,29,4),extent=(-20,-20,20,20),extend='max')
            d=ax[hl].streamplot(NULONG[:],NULAT[:],f1com[levi,hl,:,:],f2com[levi,hl,:,:],color="b")
            #e=ax[hl].contour(NULONG[:],NULAT[:],f3com[hl,:,:],levels=np.arange(11760,12600,60),color="g",linewidths=5)
            #ax[hl].clabel(e, inline=True, fontsize=12)
        
        if var==5:
            c=ax[hl].contourf(NULONG[:],NULAT[:],f12com[:,:],cmap="Reds",levels=np.arange(0,29,4),extent=(-20,-20,20,20),extend='max')
            d=ax[hl].streamplot(NULONG[:],NULAT[:],f1com[levi,hl,:,:],f2com[levi,hl,:,:],color="b")
        
        if var==6:
            c=ax[levi,hl].contourf(NULONG[:],NULAT[:],f1com[levi,hl,:,:],cmap="Blues",levels=np.arange(50,100,10),extent=(-20,-20,20,20),extend='both')
            #d=ax[hl].contour(NULONG[:],NULAT[:],f2com[hl,:,:],levels=np.arange(-0.3,-0.1,0.1),color="red",linewidths=1)
            #e=ax[hl].contour(NULONG[:],NULAT[:],f3com[levi,hl,:,:],levels=np.arange(5280,6480,15),color="g",linewidths=5)
            #ax[hl].clabel(d, inline=True, fontsize=10)
            #ax[hl].clabel(e, inline=True, fontsize=25)
        
        if hl==0:
            ax[levi,hl].set_title('>='+str(lev)+'kt/24h', fontsize=25)
        else:
            ax[levi,hl].set_title('<'+str(lev)+'kt/24h', fontsize=25)
        ax[levi,hl].set_title('('+title_num[levi][hl]+')',loc='left', fontsize=25)
        ax[3,hl].set_xlabel('$(^o)$', fontsize=25)
        ax[levi,0].set_ylabel('$(^o)$', fontsize=25)
        
cbar_ax = fig.add_axes([0.82, 0.15, 0.02, 0.7]) #x,y,dx,dy
cbar=fig.colorbar(c, cax=cbar_ax)
cbar.ax.tick_params(labelsize=25)
fig.subplots_adjust(right=0.8)

if drawt==1:
    fig.savefig('fig/i_'+title3+'.png', dpi=300)
if drawt==2:
    fig.savefig('fig/i24_'+title3+'.png', dpi=300)
