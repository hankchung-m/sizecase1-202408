#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 15:20:42 2022

@author: cwuhank
"""


import numpy as np
from datetime import date as dt
from math import*
import matplotlib.pyplot as plt
import netCDF4 as nc
import cartopy.crs as ccrs
#from ew_composite import ew_composite
#from vorticity import vorticity

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

    

#1:ew
#2:mon
#ew=1

LMIR=np.load('LMIR.npy')
tdrmw=np.load('tdrmw.npy')
tdlat=np.load('tdlat.npy')
tdlon=np.load('tdlon.npy')
tdew=np.load('tdew.npy')
tdtimestr=np.load('tdtimestr.npy', allow_pickle=True)
tdtimestr_sat=np.load('tdtimestr_sat.npy', allow_pickle=True)
tdyear=np.load('tdyear.npy', allow_pickle=True)
tdhour=np.load('tdhour.npy', allow_pickle=True)


interxy=125 #2000km/4km/4(+-1000km/4)

file_wm=['_w','_m']
title_num=[['a','b'],['c','d']]
#!!!
for wm in range(1,2):
    tdew=np.load('tdew'+file_wm[wm]+'.npy')
    print(file_wm[wm])
    
    #draw https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.streamplot.html
    fig,ax=plt.subplots(2,2,sharex=True,sharey=True,figsize=(18,18))
    
    for ew in range(1,3):
        
        
        """
        
        
        ewfill=np.zeros((int(interxy*4+1),int(interxy*4+1)))
        monfill=np.zeros((int(interxy*4+1),int(interxy*4+1)))
        f1com=np.zeros((2,int(interxy*4+1),int(interxy*4+1)))
        
        
        for tt in range(tdtimestr.size):
            
            #print(tt)
            
            #sort real timestr
            if testdt(tdtimestr[tt]) == False :
                continue
            
            a=tdhour[tt]
            if a != '00' and a != '03' and a != '06' and a != '09' and a != '12' and a != '15' and a != '18' and a != '21':
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
                        
                            if tdew[tt] == type1 :
                                ewfill[j-low,i-left]+=1
                                f1com[0,j-low,i-left]+=f11[j,ni]
                            
                            elif tdew[tt] == type2 :
                                monfill[j-low,i-left]+=1
                                f1com[1,j-low,i-left]+=f11[j,ni]
                    
                    else:
                        #sort missing value
                        if f11[j,i]>=150:
                        
                            if tdew[tt] == type1 :
                                ewfill[j-low,i-left]+=1
                                f1com[0,j-low,i-left]+=f11[j,i]
                            
                            elif tdew[tt] == type2 :
                                monfill[j-low,i-left]+=1
                                f1com[1,j-low,i-left]+=f11[j,i]
         
        
        for i in range(left,left+int(interxy*4)+1,1) :
            for j in range(low,low+int(interxy*4)+1,1) :
                f1com[0,j-low,i-left]=f1com[0,j-low,i-left]/ewfill[j-low,i-left]
                f1com[1,j-low,i-left]=f1com[1,j-low,i-left]/monfill[j-low,i-left]
            
        
        # npy
        np.save('f1com'+file_wm[wm]+'_sat_ew'+str(ew),f1com)
        
        """
        
        f1com=np.load('f1com'+file_wm[wm]+'_sat_ew'+str(ew)+'.npy')
        
        
        
        #print(f3com)
        
        NULAT=np.arange(-int(interxy*4*4/2),int(interxy*4*4/2)+4,4)
        NULONG=np.arange(-int(interxy*4*4/2),int(interxy*4*4/2)+4,4)
        
        
        title1='Tb, '
        title3='tb_'
            
        if ew==1:
            title4='ew'
        if ew==2:
            title4='mon'
        
        for hl in range(2):
            
                if ew==1 and hl==0:
                    title2='EW'
                    rmw=251#251
                if ew==2 and hl==0:
                    if wm==0:
                        title2='CS'
                        rmw=307
                    else:
                        title2='MS'
                        rmw=256
                if ew==1 and hl==1:
                    title2='MC'
                    rmw=254#260
                if ew==2 and hl==1:
                    if wm==0:
                        title2='MS'
                        rmw=306
                    else:
                        title2='MD'
                        rmw=386
                
                ax[ew-1,hl].axis([-int(interxy*4*4/2),int(interxy*4*4/2),-int(interxy*4*4/2),int(interxy*4*4/2)])
                ax[ew-1,hl].set_xticks(np.arange(-int(interxy*4*4/2),int(interxy*4*4/2)+500,500))
                ax[ew-1,hl].set_yticks(np.arange(-int(interxy*4*4/2),int(interxy*4*4/2)+500,500))
                ax[ew-1,hl].tick_params(labelsize=25)
                ax[ew-1,hl].grid(True)
                ax[ew-1,hl].set_aspect('equal')
                #plt.axes(projection=ccrs.PlateCarree()).coastlines(resolution='auto', color='k',linewidth=0.5)
                
                c=ax[ew-1,hl].contourf(NULONG[:],NULAT[:],f1com[hl,:,:],cmap="Greys",
                                  levels=np.arange(210,290,10),
                                  extent=(-int(interxy*4*4/2),-int(interxy*4*4/2),int(interxy*4*4/2),int(interxy*4*4/2)),
                                  extend='max')
                
                nthe = 360                  # angle for azimuthal mean
                inv_the = 1                 # angle resolution for azimuthal mean 
                for the in range(0,int(nthe/inv_the)):
                    rthe=the*inv_the
                    x=rmw*cos(radians(rthe))
                    y=rmw*sin(radians(rthe))
                    ax[ew-1,hl].scatter(x,y,1,'red')
                
                ax[ew-1,hl].set_title(title2, fontsize=25)
                ax[ew-1,hl].set_title('('+title_num[ew-1][hl]+')',loc='left', fontsize=25)
                ax[1,hl].set_xlabel('(km)', fontsize=25)
                ax[ew-1,0].set_ylabel('(km)', fontsize=25)
                    
        
        #cbar_ax = fig.add_axes([0.82, 0.15, 0.02, 0.7]) #x,y,dx,dy
        cbar_ax = fig.add_axes([0.883, 0.087, 0.02, 0.8]) #x,y,dx,dy #!!!x==0.883 y==0.087 dy==0.8
        cbar=fig.colorbar(c, cax=cbar_ax)
        cbar.ax.tick_params(labelsize=25)
        cbar.set_label('(K)', fontsize=25)
        #fig.subplots_adjust(right=0.8)
        fig.subplots_adjust(left=0.1, right=0.863, bottom=0.01, top=0.99, hspace=0.01)#!!!
        
        fig.savefig('fig/'+title3+file_wm[wm]+'.png', dpi=600)

"""

#---------azi_mean---------


distance_x, distance_y = np.meshgrid(np.arange(-interxy*8,interxy*8+1,4),np.arange(-interxy*8,interxy*8+1,4))
distance = ((distance_x**2+distance_y**2)**0.5)/10
distance=10*np.round(distance)

TB1_azi=np.zeros(100)
TB2_azi=np.zeros(100)
TB1=f1com[0,:,:]
TB2=f1com[1,:,:]
print(TB1.shape)
for r in range(100):
    TB1_azi[r] = np.mean(TB1[distance==r*10])
    TB2_azi[r] = np.mean(TB2[distance==r*10])

deg=np.arange(0,1000,10)

#draw
fig,ax=plt.subplots(1,2,sharex=False,sharey=True,figsize=(18,8))

title1='azi_Tb, '
title3='azi_tb_'
    
if ew==1:
    title4='ew'
if ew==2:
    title4='mon'

for hl in range(2):
    
        ax[hl].axis([0,1000,173,293])
        
        
        if ew==1 and hl==0:
            ax[hl].plot(deg,TB1_azi,'-',color='red')
            title2='ew'
        if ew==2 and hl==0:
            ax[hl].plot(deg,TB1_azi,'-',color='purple')
            title2='mc'
        if ew==1 and hl==1:
            ax[hl].plot(deg,TB2_azi,'-',color='blue')
            title2='md'
        if ew==2 and hl==1:
            ax[hl].plot(deg,TB2_azi,'-',color='green')
            title2='ms'


        ax[hl].set_title(title1+title2)
            

fig.savefig('fig/'+title3+title4+'.png')
"""
print('tb done')