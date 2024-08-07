#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 18:26:20 2022

@author: cwuhank
"""


import numpy as np
from datetime import date as dt
from math import*
import matplotlib.pyplot as plt
import netCDF4 as nc
import cartopy.crs as ccrs
from ew_composite import ew_composite
from vorticity import vorticity
from Car2cylin import Car2cylin
from matplotlib.ticker import ScalarFormatter

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

#!!!do you want to calculate ALL cases? yes:1 no:0 (ideal cases only)
forALL=0

#1:850hPa wind, 500hPa geopotential height
#2:850hPa vorticity, 850hPa wind
#3:700hPa RH, 700hPa geopotential height
#4:200hPa wind, 200hPa geopotential height
#5:200-850hPa VWS
#6:850hPa specific humidity, 975hPa wind
#7:850hPa vorticity, 500hPa vorticity (useless)
#var=6


#1:ew
#2:mon
#ew=2

LMIR=np.load('LMIR.npy')
tdrmw=np.load('tdrmw.npy')
tdlat=np.load('tdlat_ny.npy')
tdlon=np.load('tdlon_nx.npy')
tdtimestr=np.load('tdtimestr.npy', allow_pickle=True)
tdhour=np.load('tdhour.npy', allow_pickle=True)
file_wm=['_w','_m']
title_num=[['a','b'],['c','d']]
#!!!
for wm in range(1,2):
    tdew=np.load('tdew'+file_wm[wm]+'.npy')
    print(file_wm[wm])
    for var in range(1,7):
        print(var)
        #draw https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.streamplot.html
        fig,ax=plt.subplots(2,2,sharex=True,sharey=True,figsize=(18,18))
        
        
        for ew in range(1,3):
    
            if var==6:
                interxy=20
            else:
                interxy=40 #40degree(+-20degree)
            f1=np.zeros((int(tdtimestr.size),interxy*4+1,interxy*4+1))
            f2=np.zeros((int(tdtimestr.size),interxy*4+1,interxy*4+1))
            f3=np.zeros((int(tdtimestr.size),interxy*4+1,interxy*4+1))
            f4=np.zeros((int(tdtimestr.size),interxy*4+1,interxy*4+1))
            
            
            
            for tt in range(tdtimestr.size):
                
                #print(tt)
                
                #sort real timestr
                if testdt(tdtimestr[tt]) == False :
                    f1[tt,:,:]=np.nan
                    f2[tt,:,:]=np.nan
                    f3[tt,:,:]=np.nan
                    continue
                
                a=tdhour[tt]
                if a != '00' and a != '03' and a != '06' and a != '09' and a != '12' and a != '15' and a != '18' and a != '21':
                    f1[tt,:,:]=np.nan
                    f2[tt,:,:]=np.nan
                    f3[tt,:,:]=np.nan
                    continue
            
                #!!!
                if ew==1:
                    type1=1
                    type2=-1
                    notype1=-1.5
                    notype2=-2
                if ew==2:
                    type1=-1.5
                    type2=-2
                    notype1=1
                    notype2=-1
                
                if tdew[tt] == notype1 or tdew[tt] == notype2 :
                    continue
                
                if tt/10<1 :
                    ncfile = nc.Dataset("ERA5_TD00%d.nc"%(tt))
                elif tt/100<1 :
                    ncfile = nc.Dataset("ERA5_TD0%d.nc"%(tt))
                else:
                    ncfile = nc.Dataset("ERA5_TD%d.nc"%(tt))
                    
                XLAT = ncfile["latitude"]
                XLONG = ncfile["longitude"]
                p = ncfile["level"]
                u = ncfile["u"]#(1,37,321,521)
                v = ncfile["v"]
                w = ncfile["w"]
                z = ncfile["z"]
                r = ncfile["r"]
                vo = ncfile["vo"]
                q = ncfile["q"]
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
                    f11=r[0,30,:,:]
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
                    f11=u[0,35,:,:]
                    f22=v[0,35,:,:]
                    f33=q[0,30,:,:]
                if var==7:
                    f11=u[0,30,:,:]
                    f22=vo[0,21,:,:]
                    f33=vo[0,30,:,:]
                    
                #f11=u[0,30,:,:]
                #f22=v[0,30,:,:]
                #f33=z[0,21,:,:]/9.8
                #f33=vo[0,30,:,:]
                
                
                
                #find low & left
                low=0
                left=0
                for i in range(XLONG.size) :
                    for j in range(XLAT.size) :
                        lon1=XLONG[tdlon[tt]]-(interxy/2)
                        lat1=XLAT[tdlat[tt]]-(interxy/2)
                        # abs(lat1-XLAT[j]):0-0.25/2 , abs(lon1-XLONG[i]):0-0.25/2
                        if abs(lat1-XLAT[j])<0.25/2 and abs(lat1-XLAT[j])>=0:
                            low=j
                        if abs(lon1-XLONG[i])<0.25/2 and abs(lon1-XLONG[i])>=0:
                            left=i
                            
                            
                
                #no longer interp
                for i in range(left,left+(interxy*4)+1,1) :
                    for j in range(low,low-(interxy*4)-1,-1) : #j is 0 at top
                        
                        #nlon=tslon[tt]-(interxy/2)+(0.25*(i-left))
                        #nlat=tslat[tt]-(interxy/2)+(0.25*(low-j))
                        #b=(XLONG[i+1]-XLONG[i])*(XLAT[j-1]-XLAT[j])
                            
                        #sort missing value
                        f22range=100
                        if var==1 or var==2 or var==6 or var==7:
                            f11range=100
                            f33range=7000
                        if var==3:
                            f11range=200
                            f33range=4000
                        if var==4 or var==5:
                            f11range=100
                            f33range=13000
                        if abs(f11[j,i])<f11range and abs(f22[j,i])<f22range and abs(f33[j,i])<f33range and f33[j,i]>-1:
                            f1[tt,low-j,i-left]=f11[j,i]
                            f2[tt,low-j,i-left]=f22[j,i]
                            f3[tt,low-j,i-left]=f33[j,i]
                            #f1[tt,low-j,i-left]=interp(f11[j,i],f11[j,i+1],f11[j-1,i],f11[j-1,i+1],XLONG[i],XLONG[i+1],nlon,XLAT[j],XLAT[j-1],nlat)/b
                            #f2[tt,low-j,i-left]=interp(f22[j,i],f22[j,i+1],f22[j-1,i],f22[j-1,i+1],XLONG[i],XLONG[i+1],nlon,XLAT[j],XLAT[j-1],nlat)/b
                            #f3[tt,low-j,i-left]=interp(f33[j,i],f33[j,i+1],f33[j-1,i],f33[j-1,i+1],XLONG[i],XLONG[i+1],nlon,XLAT[j],XLAT[j-1],nlat)/b
                        else:
                            f1[tt,interxy*2,interxy*2]=0
                            f2[tt,interxy*2,interxy*2]=0
                            f3[tt,interxy*2,interxy*2]=0
                            #print(f11[j,i],f22[j,i],f33[j,i])
                            break
                            
                #f1=f11[low-(interxy*4)-1:low,left:left+(interxy*4)+1]
                #f2=f22[low-(interxy*4)-1:low,left:left+(interxy*4)+1]
                #f3=f33[low-(interxy*4)-1:low,left:left+(interxy*4)+1]
                        
                
                if f1[tt,interxy*2,interxy*2]==0 or np.isnan(f1[tt,interxy*2,interxy*2])==True or f2[tt,interxy*2,interxy*2]==0 or np.isnan(f2[tt,interxy*2,interxy*2])==True or f3[tt,interxy*2,interxy*2]==0 or np.isnan(f3[tt,interxy*2,interxy*2])==True :
                    f1[tt,:,:]=np.nan
                    f2[tt,:,:]=np.nan
                    f3[tt,:,:]=np.nan
                    
                #print(f1[tt,:,:])
                
            # npy
            if forALL==1:
                np.save('f1'+file_wm[wm]+'_var'+str(var)+'_ew'+str(ew),f1)
                np.save('f2'+file_wm[wm]+'_var'+str(var)+'_ew'+str(ew),f2)
                np.save('f3'+file_wm[wm]+'_var'+str(var)+'_ew'+str(ew),f3)
            else:
                np.save('f1'+file_wm[wm]+'_var'+str(var)+'_ew'+str(ew)+'_i',f1)
                np.save('f2'+file_wm[wm]+'_var'+str(var)+'_ew'+str(ew)+'_i',f2)
                np.save('f3'+file_wm[wm]+'_var'+str(var)+'_ew'+str(ew)+'_i',f3)
            
            
            
            f1com=ew_composite(f1,interxy,type1,type2,wm)
            f2com=ew_composite(f2,interxy,type1,type2,wm)
            f3com=ew_composite(f3,interxy,type1,type2,wm)
            
            
            # npy
            if forALL==1:
                np.save('f1com'+file_wm[wm]+'_var'+str(var)+'_ew'+str(ew),f1com)
                np.save('f2com'+file_wm[wm]+'_var'+str(var)+'_ew'+str(ew),f2com)
                np.save('f3com'+file_wm[wm]+'_var'+str(var)+'_ew'+str(ew),f3com)
            else:
                np.save('f1com'+file_wm[wm]+'_var'+str(var)+'_ew'+str(ew)+'_i',f1com)
                np.save('f2com'+file_wm[wm]+'_var'+str(var)+'_ew'+str(ew)+'_i',f2com)
                np.save('f3com'+file_wm[wm]+'_var'+str(var)+'_ew'+str(ew)+'_i',f3com)
    
            """
            
            if forALL==1:
                f1com=np.load('f1com'+file_wm[wm]+'_var'+str(var)+'_ew'+str(ew)+'.npy')
                f2com=np.load('f2com'+file_wm[wm]+'_var'+str(var)+'_ew'+str(ew)+'.npy')
                f3com=np.load('f3com'+file_wm[wm]+'_var'+str(var)+'_ew'+str(ew)+'.npy')
            else:
                f1com=np.load('f1com'+file_wm[wm]+'_var'+str(var)+'_ew'+str(ew)+'_i.npy')
                f2com=np.load('f2com'+file_wm[wm]+'_var'+str(var)+'_ew'+str(ew)+'_i.npy')
                f3com=np.load('f3com'+file_wm[wm]+'_var'+str(var)+'_ew'+str(ew)+'_i.npy')
            
            """
    
            #print(f3com)
            
            NULAT=np.arange(-int(interxy/2),int(interxy/2)+0.25,0.25)
            NULONG=np.arange(-int(interxy/2),int(interxy/2)+0.25,0.25)
            
            
            if var==1:
                title1='850hPa wind speed, '
                title3='850w_500z_'
            if var==2:
                title1='850hPa vorticity, '
                title3='850v_850w_'
            if var==3:
                title1='850hPa RH, '
                title3='850r_700z_'
            if var==4:
                title1='200hPa wind speed, '
                title3='200w_200z_'
            if var==5:
                title1='200-850hPa VWS, '
                title3='200_850vws_'
            if var==6:
                title1='850hPa specific humidity, '
                title3='850q_975w_'
            if var==7:
                title1='500hPa vorticity, '
                title3='850v_500v_'
                
            if ew==1:
                title4='ew'
            if ew==2:
                title4='mon'
            
            for hl in range(2):#!!!
                
                    if ew==1 and hl==0:
                        title2='EW'
                    if ew==2 and hl==0:
                        if wm==0:
                            title2='CS'
                        else:
                            title2='MS'
                    if ew==1 and hl==1:
                        title2='MC'
                    if ew==2 and hl==1:
                        if wm==0:
                            title2='MS'
                        else:
                            title2='MD'
            
                    if var==6:
                        ax[ew-1,hl].axis([-10,10,-10,10])
                        ax[ew-1,hl].set_xticks(np.arange(-10,15,5))
                        ax[ew-1,hl].set_yticks(np.arange(-10,15,5))
                    else:
                        ax[ew-1,hl].axis([-20,20,-20,20])
                        ax[ew-1,hl].set_xticks(np.arange(-20,30,10))
                        ax[ew-1,hl].set_yticks(np.arange(-20,30,10))
                    ax[ew-1,hl].tick_params(labelsize=25)
                    ax[ew-1,hl].grid(True)
                    ax[ew-1,hl].set_aspect('equal')
                    #plt.axes(projection=ccrs.PlateCarree()).coastlines(resolution='auto', color='k',linewidth=0.5)
                    
            
                    if var==1 or var==2 or var==4 or var==5 or var==6:
                        f12com=np.zeros((interxy*4+1,interxy*4+1))
                        for i in range(interxy*4+1) :
                            for j in range(interxy*4+1) :
                                f12com[j,i]=sqrt(f1com[hl,j,i]**2+f2com[hl,j,i]**2)
            
                    if var==1:
                        c=ax[ew-1,hl].contourf(NULONG[:],NULAT[:],f12com[:,:],cmap="Reds",levels=np.arange(4,16,2),extent=(-20,-20,20,20),extend='max')
                        d=ax[ew-1,hl].streamplot(NULONG[:],NULAT[:],f1com[hl,:,:],f2com[hl,:,:],color="b")
                        e=ax[ew-1,hl].contour(NULONG[:],NULAT[:],f3com[hl,:,:],levels=np.arange(5880,5895,15),colors="g",linewidths=10)
                        ax[ew-1,hl].clabel(e, inline=True, fontsize=25, fmt="%d")
                        #cbar_ax = fig.add_axes([0.82, 0.15, 0.02, 0.7]) #x,y,dx,dy
                        cbar_ax = fig.add_axes([0.87, 0.087, 0.02, 0.8]) #x,y,dx,dy #!!!x==0.87 y==0.087 dy==0.8
                        cbar=fig.colorbar(c, cax=cbar_ax)
                        cbar.ax.tick_params(labelsize=25)
                        cbar.set_label('$[m \\, s^{-1}]$', fontsize=25)
                    
                    if var==2:
                        c=ax[ew-1,hl].contourf(NULONG[:],NULAT[:],f3com[hl,:,:],cmap="bwr",levels=np.arange(-1.5*0.00001,1.8*0.00001,0.3*0.00001),extent=(-20,-20,20,20),extend='both')
                        d=ax[ew-1,hl].streamplot(NULONG[:],NULAT[:],f1com[hl,:,:],f2com[hl,:,:],color="k")
                        #e=ax[ew-1,hl].contour(NULONG[:],NULAT[:],f12com[:,:],levels=np.arange(8,20,4),colors="c",linewidths=5)
                        #ax[ew-1,hl].clabel(e, inline=True, fontsize=25)
                        #cbar_ax = fig.add_axes([0.82, 0.15, 0.02, 0.7]) #x,y,dx,dy
                        cbar_ax = fig.add_axes([0.87, 0.087, 0.02, 0.8]) #x,y,dx,dy #!!!x==0.87 y==0.087 dy==0.8
                        #cbar=fig.colorbar(c, cax=cbar_ax)#*100000
                        cbar=fig.colorbar(c, cax=cbar_ax, format=ScalarFormatter(useMathText=True))
                        # 調整科學記號文字大小
                        cbar.formatter._exponent = 'e'
                        cbar.formatter.set_powerlimits((0, 0))
                        
                        # 擴大科學記號的文字大小
                        cbar_ax.yaxis.offsetText.set(size=25)
                        
                        # 調整科學記號在水平方向上的位置
                        cbar_ax.yaxis.offsetText.set(x=2.5)
                        cbar.ax.tick_params(labelsize=25)
                        cbar.set_label('$[s^{-1}]$', fontsize=25)
                        
                        #RMW
                        VTcom=np.zeros((interxy*4+1,interxy*4+1))
                                
                        for i in range(interxy*4+1) :
                            for j in range(interxy*4+1) :
                                erai=i-(interxy*2)
                                eraj=j-(interxy*2)
                                distance=sqrt((eraj**2)+(erai**2))
                                VTcom[j,i] = -f1com[hl,j,i]*eraj/distance+f2com[hl,j,i]*erai/distance
                                
                        close_range = 20        # 20 degree square
                        close_range_g = 4*close_range        # 20*4 grids square
                        Rdim = close_range_g-1        # radial grid points for azimuthal mean
                        inv_r = 1                   # radial resolution for azimuthal mean (0.2 grid, 5km)
                        nthe = 360                  # angle for azimuthal mean
                        inv_the = 5                 # angle resolution for azimuthal mean 
                        cf12com=Car2cylin(close_range_g, close_range_g, 1, Rdim, inv_r, nthe, inv_the, VTcom[np.newaxis,:,:])[0]
                        avgrmw=0
                        avg=0
                        for the in range(0,int(nthe/inv_the)):
                            rmwloc=np.where(cf12com[the,:]==np.nanmax(cf12com[the,:]))[0][0]
                            if rmwloc>0:
                                rthe=the*inv_the
                                rmw=(rmwloc)*inv_r
                                avgrmw+=rmw
                                avg+=1
                                #x=rmw/4*cos(radians(rthe))
                                #y=rmw/4*sin(radians(rthe))
                                #ax[ew-1,hl].scatter(x,y,30,'yellow')
                        avgrmw=avgrmw/avg/4*111
                        print(title2,np.around(avgrmw,decimals=0))
                    
                    if var==3:
                        c=ax[ew-1,hl].contourf(NULONG[:],NULAT[:],f1com[hl,:,:],cmap="Blues",levels=np.arange(50,100,10),extent=(-20,-20,20,20),extend='both')
                        #d=ax[ew-1,hl].contour(NULONG[:],NULAT[:],f2com[hl,:,:],levels=np.arange(-0.3,-0.1,0.1),color="red",linewidths=1)
                        e=ax[ew-1,hl].contour(NULONG[:],NULAT[:],f3com[hl,:,:],levels=np.arange(2820,3420,30),color="g",linewidths=5)
                        #ax[hl].clabel(d, inline=True, fontsize=10)
                        ax[ew-1,hl].clabel(e, inline=True, fontsize=25, fmt="%d")
                        #cbar_ax = fig.add_axes([0.82, 0.15, 0.02, 0.7]) #x,y,dx,dy
                        cbar_ax = fig.add_axes([0.87, 0.087, 0.02, 0.8]) #x,y,dx,dy #!!!x==0.87 y==0.087 dy==0.8
                        cbar=fig.colorbar(c, cax=cbar_ax)
                        cbar.ax.tick_params(labelsize=25)
                        cbar.set_label('[%]', fontsize=25)
                    
                    if var==4:
                        c=ax[ew-1,hl].contourf(NULONG[:],NULAT[:],f12com[:,:],cmap="Reds",levels=np.arange(0,29,4),extent=(-20,-20,20,20),extend='max')
                        d=ax[ew-1,hl].streamplot(NULONG[:],NULAT[:],f1com[hl,:,:],f2com[hl,:,:],color="b")
                        #e=ax[ew-1,hl].contour(NULONG[:],NULAT[:],f3com[hl,:,:],levels=np.arange(11760,12600,60),color="g",linewidths=5)
                        #ax[ew-1,hl].clabel(e, inline=True, fontsize=15)
                        cbar_ax = fig.add_axes([0.82, 0.15, 0.02, 0.7]) #x,y,dx,dy
                        cbar=fig.colorbar(c, cax=cbar_ax)
                        cbar.ax.tick_params(labelsize=25)
                        cbar.set_label('$[m \\, s^{-1}]$', fontsize=25)
                    
                    if var==5:
                        c=ax[ew-1,hl].contourf(NULONG[:],NULAT[:],f12com[:,:],cmap="Reds",levels=np.arange(0,29,4),extent=(-20,-20,20,20),extend='max')
                        d=ax[ew-1,hl].streamplot(NULONG[:],NULAT[:],f1com[hl,:,:],f2com[hl,:,:],color="b")
                        cbar_ax = fig.add_axes([0.82, 0.15, 0.02, 0.7]) #x,y,dx,dy
                        cbar=fig.colorbar(c, cax=cbar_ax)
                        cbar.ax.tick_params(labelsize=25)
                        cbar.set_label('$[m \\, s^{-1}]$', fontsize=25)
                        
                        #VWS
                        
                        
                    if var==6:
                        c=ax[ew-1,hl].contourf(NULONG[:],NULAT[:],f3com[hl,:,:],cmap="Blues",levels=np.arange(0.012,0.014,0.0005),extent=(-10,-10,10,10),extend='both')
                        d=ax[ew-1,hl].streamplot(NULONG[:],NULAT[:],f1com[hl,:,:],f2com[hl,:,:],color="k")
                        #e=ax[ew-1,hl].contour(NULONG[:],NULAT[:],f12com[:,:],levels=np.arange(0,16,4),color="y",linewidths=5)
                        #ax[ew-1,hl].clabel(e, inline=True, fontsize=15)
                        #cbar_ax = fig.add_axes([0.82, 0.15, 0.02, 0.7]) #x,y,dx,dy
                        cbar_ax = fig.add_axes([0.87, 0.087, 0.02, 0.8]) #x,y,dx,dy #!!!x==0.87 y==0.087 dy==0.8
                        cbar=fig.colorbar(c, cax=cbar_ax)
                        cbar.ax.tick_params(labelsize=25)
                        cbar.set_label('$[kg \\, kg^{-1}]$', fontsize=25)
                    
                    if var==7:#useless
                        c=ax[ew-1,hl].contourf(NULONG[:],NULAT[:],f2com[hl,:,:]*100000,cmap="bwr",levels=np.arange(-1.5,1.8,0.3),extent=(-20,-20,20,20),extend='both')
                        #d=ax[ew-1,hl].streamplot(NULONG[:],NULAT[:],f1com[hl,:,:],f2com[hl,:,:],color="k")
                        e=ax[ew-1,hl].contour(NULONG[:],NULAT[:],f3com[hl,:,:]*100000,levels=np.arange(1.2,1.5,0.3),colors="c",linewidths=5)
                        ax[ew-1,hl].clabel(e, inline=True, fontsize=25)
                        cbar_ax = fig.add_axes([0.82, 0.15, 0.02, 0.7]) #x,y,dx,dy
                        cbar=fig.colorbar(c, cax=cbar_ax)
                        cbar.ax.tick_params(labelsize=25)
                        cbar.set_label('$(10^-5)(s^-1)$', fontsize=25)
                    
                    ax[ew-1,hl].set_title(title2, fontsize=25)
                    ax[ew-1,hl].set_title('('+title_num[ew-1][hl]+')',loc='left', fontsize=25)
                    ax[1,hl].set_xlabel('$(^o)$', fontsize=25)
                    ax[ew-1,0].set_ylabel('$(^o)$', fontsize=25)
                        
        #fig.subplots_adjust(right=0.8)
        fig.subplots_adjust(left=0.087, right=0.85, bottom=0.01, top=0.99, hspace=0.01)#!!!
            
        if forALL==1:
            fig.savefig('fig/'+title3+file_wm[wm]+'.png', dpi=600)
        else:
            fig.savefig('fig/'+title3+file_wm[wm]+'_i.png', dpi=600)
        plt.show()
        
        