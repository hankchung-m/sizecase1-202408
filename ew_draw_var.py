#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 01:04:30 2022

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

def ew_var(field,interxy,type1,type2,type1num,type2num,wm):
    
    file_wm=['_w','_m']
    tdew=np.load('tdew'+file_wm[wm]+'.npy')
    tdtimestr=np.load('tdtimestr.npy', allow_pickle=True)

    
    fcom1=np.zeros((int(type1num),interxy*4+1,interxy*4+1)) #(ew/mon , interxy*4+1,interxy*4+1)
    fcom2=np.zeros((int(type2num),interxy*4+1,interxy*4+1))
    fcom=np.zeros((2,interxy*4+1,interxy*4+1))
    ew=0
    mon=0
        
    for tt in range(tdtimestr.size):
        if np.isnan(field[tt,interxy*2,interxy*2])==False :
            if tdew[tt] == type1 :
                fcom1[ew,:,:]=field[tt,:,:]
                ew+=1
            elif tdew[tt] == type2 :
                fcom2[mon,:,:]=field[tt,:,:]
                mon+=1
    
    for i in range((interxy*4)+1) :
        for j in range((interxy*4)+1) :
            fcom[0,j,i]=np.var(fcom1[:,j,i])
            fcom[1,j,i]=np.var(fcom2[:,j,i])
        
        
    return fcom


#!!!do you want to calculate ALL cases? yes:1 no:0 (ideal cases only)
forALL=0

#1:850hPa wind, 500hPa geopotential height
#2:850hPa vorticity, 850hPa wind
#3:700hPa RH, 700hPa geopotential height
#4:200hPa wind, 200hPa geopotential height
#5:200-850hPa VWS
#6:850hPa specific humidity, 975hPa wind
#7:850hPa wind, 500hPa geopotential height
#var=6
    

#1:ew
#2:mon
#ew=2

LMIR=np.load('LMIR.npy')
tdrmw=np.load('tdrmw.npy')
tdlat=np.load('tdlat_ny.npy')
tdlon=np.load('tdlon_nx.npy')
#tdew=np.load('tdew.npy')
tdtimestr=np.load('tdtimestr.npy', allow_pickle=True)
tdhour=np.load('tdhour.npy', allow_pickle=True)

#!!!
file_wm=['_w','_m']
title_num=[['a','b'],['c','d']]

for wm in range(1,2):
    tdew=np.load('tdew'+file_wm[wm]+'.npy')
    print(file_wm[wm])
    
    for var in range(8,9):
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
            type1num=0
            type2num=0
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
                
                if tdew[tt] == type1:
                    type1num+=1
                    
                if tdew[tt] == type2:
                    type2num+=1
                
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
                
                if var==7:
                    f11=u[0,30,:,:]
                    f22=v[0,30,:,:]
                    f33=z[0,21,:,:]/9.8
                if var==2:
                    f11=u[0,30,:,:]
                    f22=v[0,30,:,:]
                    f33=vo[0,30,:,:]
                if var==8:
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
                        if var==7 or var==2 or var==6:
                            f11range=100
                            f33range=7000
                        if var==3 or var==8:
                            f11range=200
                            f33range=4000
                        if var==4 or var==5:
                            f11range=100
                            f33range=13000
                        if abs(f11[j,i])<f11range and abs(f22[j,i])<f22range and abs(f33[j,i])<f33range and f33[j,i]>-1:
                            f1[tt,low-j,i-left]=f11[j,i]
                            f2[tt,low-j,i-left]=f22[j,i]
                            f3[tt,low-j,i-left]=f33[j,i]
                            #if tdew[tt]==type1:
                            #    type1f1[tt,low-j,i-left]=f11[j,i]
                            #if tdew[tt]==type2:
                            #    type2f1[tt,low-j,i-left]=f11[j,i]
                            
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
                np.save('f1'+file_wm[wm]+'_var_var'+str(var)+'_ew'+str(ew),f1)
                np.save('f2'+file_wm[wm]+'_var_var'+str(var)+'_ew'+str(ew),f2)
                np.save('f3'+file_wm[wm]+'_var_var'+str(var)+'_ew'+str(ew),f3)
            else:
                np.save('f1'+file_wm[wm]+'_var_var'+str(var)+'_ew'+str(ew)+'_i',f1)
                np.save('f2'+file_wm[wm]+'_var_var'+str(var)+'_ew'+str(ew)+'_i',f2)
                np.save('f3'+file_wm[wm]+'_var_var'+str(var)+'_ew'+str(ew)+'_i',f3)
            
            f1com=ew_var(f1,interxy,type1,type2,type1num,type2num,wm)
            f2com=ew_var(f2,interxy,type1,type2,type1num,type2num,wm)
            f3com=ew_var(f3,interxy,type1,type2,type1num,type2num,wm)
            
            
            # npy
            if forALL==1:
                np.save('f1com'+file_wm[wm]+'_var_var'+str(var)+'_ew'+str(ew),f1com)
                np.save('f2com'+file_wm[wm]+'_var_var'+str(var)+'_ew'+str(ew),f2com)
                np.save('f3com'+file_wm[wm]+'_var_var'+str(var)+'_ew'+str(ew),f3com)
            else:
                np.save('f1com'+file_wm[wm]+'_var_var'+str(var)+'_ew'+str(ew)+'_i',f1com)
                np.save('f2com'+file_wm[wm]+'_var_var'+str(var)+'_ew'+str(ew)+'_i',f2com)
                np.save('f3com'+file_wm[wm]+'_var_var'+str(var)+'_ew'+str(ew)+'_i',f3com)
    
            """
            
            if forALL==1:
                f1com=np.load('f1com'+file_wm[wm]+'_var_var'+str(var)+'_ew'+str(ew)+'.npy')
                f2com=np.load('f2com'+file_wm[wm]+'_var_var'+str(var)+'_ew'+str(ew)+'.npy')
                f3com=np.load('f3com'+file_wm[wm]+'_var_var'+str(var)+'_ew'+str(ew)+'.npy')
                #f1com=np.load('f1com'+file_wm[wm]+'_var_var'+str(var)+'ew_'+str(ew)+'.npy')
                #f2com=np.load('f2com'+file_wm[wm]+'_var_var'+str(var)+'ew_'+str(ew)+'.npy')
                #f3com=np.load('f3com'+file_wm[wm]+'_var_var'+str(var)+'ew_'+str(ew)+'.npy')
                
            else:
                f1com=np.load('f1com'+file_wm[wm]+'_var_var'+str(var)+'_ew'+str(ew)+'_i.npy')
                f2com=np.load('f2com'+file_wm[wm]+'_var_var'+str(var)+'_ew'+str(ew)+'_i.npy')
                f3com=np.load('f3com'+file_wm[wm]+'_var_var'+str(var)+'_ew'+str(ew)+'_i.npy')
            
            """
    
            
            if var==6:
                var_o=6
                f3com_o=np.load('f3com'+file_wm[wm]+'_var'+str(var_o)+'_ew'+str(ew)+'_i.npy')
            if var==7:
                var_o=2
                f1com_o=np.load('f1com'+file_wm[wm]+'_var'+str(var_o)+'_ew'+str(ew)+'_i.npy')
                f2com_o=np.load('f2com'+file_wm[wm]+'_var'+str(var_o)+'_ew'+str(ew)+'_i.npy')
            if var==8:
                var_o=3
                f1com_o=np.load('f1com'+file_wm[wm]+'_var'+str(var_o)+'_ew'+str(ew)+'_i.npy')
            #print(f3com)
            
            NULAT=np.arange(-int(interxy/2),int(interxy/2)+0.25,0.25)
            NULONG=np.arange(-int(interxy/2),int(interxy/2)+0.25,0.25)
            
            
            if var==7:
                title1='850hPa wind VAR, '
                title3='850w_500z_'
            if var==2:
                title1='850hPa vorticity, 850hPa wind, '
                title3='850v_850w_'
            if var==8:
                title1='850hPa RH, 700hPa geopotential height, '
                title3='700r_700z_'
            if var==4:
                title1='200hPa wind, '
                title3='200w_200z_'
            if var==5:
                title1='200-850hPa VWS, '
                title3='200_850vws_'
            if var==6:
                title1='850hPa specific humidity VAR, '
                title3='850q_975w_'
                
            if ew==1:
                title4='ew'
            if ew==2:
                title4='mon'
            
            for hl in range(2):
                
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
                    else:
                        ax[ew-1,hl].axis([-20,20,-20,20])
                    ax[ew-1,hl].tick_params(labelsize=25)
                    ax[ew-1,hl].grid(True)
                    ax[ew-1,hl].set_aspect('equal')
                    #plt.axes(projection=ccrs.PlateCarree()).coastlines(resolution='auto', color='k',linewidth=0.5)
                    
            
                    #if var==7 or var==2 or var==4 or var==5 or var==6:
                    #    f12com=np.zeros((interxy*4+1,interxy*4+1))
                    #    for i in range(interxy*4+1) :
                    #        for j in range(interxy*4+1) :
                    #            f12com[j,i]=sqrt(f1com[hl,j,i]**2+f2com[hl,j,i]**2)
            
                    if var==7:
                        c=ax[ew-1,hl].contourf(NULONG[:],NULAT[:],f1com[hl,:,:],cmap="Reds",levels=np.arange(5,50,5),extent=(-20,-20,20,20),extend='max')
                        d=ax[ew-1,hl].streamplot(NULONG[:],NULAT[:],f1com_o[hl,:,:],f2com_o[hl,:,:],color="b")
                        #d=ax[ew-1,hl].streamplot(NULONG[:],NULAT[:],f1com[hl,:,:],f2com[hl,:,:],color="b")
                        #e=ax[ew-1,hl].contour(NULONG[:],NULAT[:],f2com[hl,:,:],levels=np.arange(20,100,10),colors="g",linewidths=5)
                        e=ax[ew-1,hl].contour(NULONG[:],NULAT[:],f2com[hl,:,:],levels=np.arange(20,70,10),cmap='YlGn',linewidths=5)
                        ax[ew-1,hl].clabel(e, inline=True, fontsize=25, fmt="%d")
                    
                    if var==2:
                        c=ax[ew-1,hl].contourf(NULONG[:],NULAT[:],f3com[hl,:,:],cmap="bwr",levels=np.arange(-0.000015,0.000018,0.000003),extent=(-20,-20,20,20),extend='both')
                        d=ax[ew-1,hl].streamplot(NULONG[:],NULAT[:],f1com[hl,:,:],f2com[hl,:,:],color="k")
                        #e=ax[hl].contour(NULONG[:],NULAT[:],f12com[:,:],levels=np.arange(0,16,4),color="y",linewidths=5)
                        #ax[hl].clabel(e, inline=True, fontsize=15)
                    
                    if var==8:
                        c=ax[ew-1,hl].contourf(NULONG[:],NULAT[:],f1com[hl,:,:],cmap="Blues",levels=np.arange(50,100,10),extent=(-20,-20,20,20),extend='both')
                        #d=ax[hl].contour(NULONG[:],NULAT[:],f2com[hl,:,:],levels=np.arange(-0.3,-0.1,0.1),color="red",linewidths=1)
                        e=ax[ew-1,hl].contour(NULONG[:],NULAT[:],f1com_o[hl,:,:],levels=np.arange(80,100,10),color="g",linewidths=5)
                        #ax[hl].clabel(d, inline=True, fontsize=10)
                        ax[ew-1,hl].clabel(e, inline=True, fontsize=15)
                    
                    if var==4:
                        #c=ax[hl].contourf(NULONG[:],NULAT[:],f12com[:,:],cmap="Reds",levels=np.arange(0,29,4),extent=(-20,-20,20,20),extend='max')
                        d=ax[hl].streamplot(NULONG[:],NULAT[:],f1com[hl,:,:],f2com[hl,:,:],color="b")
                        #e=ax[hl].contour(NULONG[:],NULAT[:],f3com[hl,:,:],levels=np.arange(11760,12600,60),color="g",linewidths=5)
                        #ax[hl].clabel(e, inline=True, fontsize=15)
                    
                    if var==5:
                        #c=ax[hl].contourf(NULONG[:],NULAT[:],f12com[:,:],cmap="Reds",levels=np.arange(0,29,4),extent=(-20,-20,20,20),extend='max')
                        d=ax[hl].streamplot(NULONG[:],NULAT[:],f1com[hl,:,:],f2com[hl,:,:],color="b")
                        
                    if var==6:
                        c=ax[ew-1,hl].contour(NULONG[:],NULAT[:],f3com_o[hl,:,:],levels=np.arange(0.013,0.063,0.005),color="k",linewidths=5)
                        c=ax[ew-1,hl].contourf(NULONG[:],NULAT[:],f3com[hl,:,:],cmap="Blues",levels=np.arange(0.0000005,0.0000065,0.0000005),extent=(-10,-10,10,10),extend='max')
                        #d=ax[ew-1,hl].streamplot(NULONG[:],NULAT[:],f1com[hl,:,:],f2com[hl,:,:],color="k")
                        #e=ax[ew-1,hl].contour(NULONG[:],NULAT[:],f12com[:,:],levels=np.arange(0,16,4),color="y",linewidths=5)
                        #ax[ew-1,hl].clabel(e, inline=True, fontsize=25)
                    
                    ax[ew-1,hl].set_title(title2, fontsize=25)
                    ax[ew-1,hl].set_title('('+title_num[ew-1][hl]+')',loc='left', fontsize=25)
                    ax[1,hl].set_xlabel('$(^o)$', fontsize=25)
                    ax[ew-1,0].set_ylabel('$(^o)$', fontsize=25)
                        
        
        #cbar_ax = fig.add_axes([0.82, 0.15, 0.02, 0.7]) #x,y,dx,dy
        cbar_ax = fig.add_axes([0.87, 0.087, 0.02, 0.8]) #x,y,dx,dy #!!!x==0.87 y==0.087 dy==0.8
        if var==6:
            cbar=fig.colorbar(c, cax=cbar_ax, format=ScalarFormatter(useMathText=True))
            # 調整科學記號文字大小
            cbar.formatter._exponent = 'e'
            cbar.formatter.set_powerlimits((0, 0))
            
            # 擴大科學記號的文字大小
            cbar_ax.yaxis.offsetText.set(size=25)
            
            # 調整科學記號在水平方向上的位置
            cbar_ax.yaxis.offsetText.set(x=2.5)
        else:
            cbar=fig.colorbar(c, cax=cbar_ax)
        cbar.ax.tick_params(labelsize=25)

        #fig.subplots_adjust(right=0.8)
        fig.subplots_adjust(left=0.087, right=0.85, bottom=0.01, top=0.99, hspace=0.01)#!!!
        
        if forALL==1:
            fig.savefig('fig/var_'+title3+file_wm[wm]+'.png', dpi=600)
        else:
            fig.savefig('fig/var_'+title3+file_wm[wm]+'_i.png', dpi=600)

    
    