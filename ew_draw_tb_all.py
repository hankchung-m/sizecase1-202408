#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 23:49:19 2023

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
#(3:700hPa RH, 700hPa geopotential height
#(4:200hPa wind, 200hPa geopotential height
#(5:200-850hPa VWS
#6:850hPa specific humidity, 975hPa wind
#7:TB
#var=2
    


LMIR=np.load('LMIR.npy')
tdrmw=np.load('tdrmw.npy')
tdlat=np.load('tdlat.npy')
tdlon=np.load('tdlon.npy')
tdew=np.load('tdew_m.npy')
tdtimestr=np.load('tdtimestr.npy', allow_pickle=True)
tdtimestr_sat=np.load('tdtimestr_sat.npy', allow_pickle=True)
tdhour=np.load('tdhour.npy', allow_pickle=True)
name=np.load('name.npy', allow_pickle=True)
tdyear=np.load('tdyear.npy', allow_pickle=True)

for var in range(7,8):
    if var==7:
        interxy=125
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
            continue
        
        a=tdhour[tt]
        if a != '00' and a != '03' and a != '06' and a != '09' and a != '12' and a != '15' and a != '18' and a != '21':
            continue
        
        
        
        ncfile = nc.Dataset('/Dellwork5/cwutyp/GPM_MERGIR/'+tdyear[tt]+'/merg_'+tdtimestr_sat[tt]+'_4km-pixel.nc4')
            
        XLAT = ncfile["lat"]
        XLONG = ncfile["lon"]
        tb = ncfile["Tb"]
        
        
        f11=tb[0,:,:]
        
        #find low & left
        low=0
        left=0
        for i in range(7144,XLONG.size) :
            if abs(tdlon[tt]-XLONG[i])>1:
                continue
            for j in range(816,XLAT.size) :
                if abs(tdlat[tt]-XLAT[j])>1:
                    continue
                
                # abs(lat1-XLAT[j]):0-0.04/2 , abs(lon1-XLONG[i]):0-0.04/2
                if abs(tdlat[tt]-XLAT[j])<0.04/2 and abs(tdlat[tt]-XLAT[j])>=0:
                    low=j-(interxy*2)
                if abs(tdlon[tt]-XLONG[i])<0.04/2 and abs(tdlon[tt]-XLONG[i])>=0:
                    left=i-(interxy*2)
                    
                if low!=0 and left!=0:
                    break
            if low!=0 and left!=0:
                break
                    
        
        #no longer interp
        for i in range(left,left+int(interxy*4)+1,1) :
            for j in range(low,low+int(interxy*4)+1,1) : #j is 0 at bottom
                    
                #near 180 degree
                if i>=9896:
                    ni=i-9896
                    #sort missing value
                    if f11[j,ni]>=150:
                    
                        f1[tt,j-low,i-left]=f11[j,ni]
                
                else:
                    #sort missing value
                    if f11[j,i]>=150:
                    
                        f1[tt,j-low,i-left]=f11[j,i]
     

        
        #f1[tt,:,:]=np.flip(f1[tt,:,:],axis=0)
        #f2[tt,:,:]=np.flip(f2[tt,:,:],axis=0)
        #f3[tt,:,:]=np.flip(f3[tt,:,:],axis=0)
        """ 
        NULAT=np.arange(-int(interxy/2),int(interxy/2)+0.25,0.25)
        NULONG=np.arange(-int(interxy/2),int(interxy/2)+0.25,0.25)
        
        #draw https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.streamplot.html
        fig,ax=plt.subplots(1,1,sharex=False,sharey=True,figsize=(9,8))
    
        if var==6:
            ax.axis([-10,10,-10,10])
        else:
            ax.axis([-20,20,-20,20])
        ax.tick_params(labelsize=15)
        ax.grid(True)
        #plt.axes(projection=ccrs.PlateCarree()).coastlines(resolution='auto', color='k',linewidth=0.5)
        
    
        if var==1 or var==6:
            f12com=np.zeros((interxy*4+1,interxy*4+1))
            for i in range(interxy*4+1) :
                for j in range(interxy*4+1) :
                    f12com[j,i]=sqrt(f1[tt,j,i]**2+f2[tt,j,i]**2)
    
        if var==1:
            c=ax.contourf(NULONG[:],NULAT[:],f12com[:,:],cmap="Reds",levels=np.arange(0,15,2),extent=(-20,-20,20,20),extend='max')
            d=ax.streamplot(NULONG[:],NULAT[:],f1[tt,:,:],f2[tt,:,:],color="b")
            #e=ax.contour(NULONG[:],NULAT[:],f3[tt,:,:],levels=np.arange(0.012,0.02,0.001),color="g",linewidths=3)
            #ax.clabel(e, inline=True, fontsize=15)
        if var==2:
            c=ax.contourf(NULONG[:],NULAT[:],f3[tt,:,:],cmap="bwr",levels=np.arange(-0.000015,0.000018,0.000003),extent=(-20,-20,20,20),extend='both')
            d=ax.streamplot(NULONG[:],NULAT[:],f1[tt,:,:],f2[tt,:,:],color="k")
            #e=ax.contour(NULONG[:],NULAT[:],f12com[:,:],levels=np.arange(0,16,4),color="y",linewidths=5)
            #ax.clabel(e, inline=True, fontsize=15)
                    
    
        if var==6:
            c=ax.contourf(NULONG[:],NULAT[:],f3[tt,:,:],cmap="Blues",levels=np.arange(0.012,0.015,0.0005),extent=(-10,-10,10,10),extend='both')
            d=ax.streamplot(NULONG[:],NULAT[:],f1[tt,:,:],f2[tt,:,:],color="k")
            e=ax.contour(NULONG[:],NULAT[:],f12com[:,:],levels=np.arange(12,16,4),color="r",linewidths=5)
            ax.clabel(e, inline=True, fontsize=15)
        
        """
        
        
        NULAT=np.arange(-int(interxy*4*4/2),int(interxy*4*4/2)+4,4)
        NULONG=np.arange(-int(interxy*4*4/2),int(interxy*4*4/2)+4,4)
        
        #draw
        fig,ax=plt.subplots(1,1,sharex=False,sharey=True,figsize=(9,8))
        
        title1='Tb, '
        title3='tb_'


        ax.axis([-int(interxy*4*4/2),int(interxy*4*4/2),-int(interxy*4*4/2),int(interxy*4*4/2)])
        ax.set_xticks(np.arange(-int(interxy*4*4/2),int(interxy*4*4/2)+500,500))
        ax.set_yticks(np.arange(-int(interxy*4*4/2),int(interxy*4*4/2)+500,500))
        ax.tick_params(labelsize=25)
        ax.grid(True)
        #plt.axes(projection=ccrs.PlateCarree()).coastlines(resolution='auto', color='k',linewidth=0.5)
        
        c=ax.contourf(NULONG[:],NULAT[:],f1[tt,:,:],cmap="Greys",
                          levels=np.arange(210,290,10),
                          extent=(-int(interxy*4*4/2),-int(interxy*4*4/2),int(interxy*4*4/2),int(interxy*4*4/2)),
                          extend='max')
                

        if tdew[tt]==1:
            plt.title(str(tt)+' '+name[tt]+'('+tdyear[tt]+'),EW', fontsize=15)
        if tdew[tt]==-1:
            plt.title(str(tt)+' '+name[tt]+'('+tdyear[tt]+'),MC', fontsize=15)
        if tdew[tt]==-1.5:
            plt.title(str(tt)+' '+name[tt]+'('+tdyear[tt]+'),MS', fontsize=15)
        if tdew[tt]==-2:
            plt.title(str(tt)+' '+name[tt]+'('+tdyear[tt]+'),MD', fontsize=15)
                

            
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.82, 0.15, 0.02, 0.7]) #x,y,dx,dy
        cbar=fig.colorbar(c, cax=cbar_ax)
        cbar.ax.tick_params(labelsize=25)
        cbar.set_label('(K)', fontsize=25)
        
        
        if var==7:
            if tt/10<1 :
                fig.savefig("fig/850tb/00%d.png"%(tt), dpi=300)
            elif tt/100<1 :
                fig.savefig("fig/850tb/0%d.png"%(tt), dpi=300)
            else:
                fig.savefig("fig/850tb/%d.png"%(tt), dpi=300)
                
        plt.show()