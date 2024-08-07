#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  9 11:58:43 2022

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
import scipy.stats as stats

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

def j_composite(field,interxy,type1,type2,type1num,type2num,tdew):
    

    #tdew=np.load('tdew.npy')
    #tdtimestr=np.load('tdtimestr.npy', allow_pickle=True)

    
    fcom=np.zeros((2,interxy*4+1,interxy*4+1)) #(ew/mon , interxy*4+1,interxy*4+1)
    ew=0
    mon=0
        
    for tt in range(281,431):
        if np.isnan(field[tt,interxy*2,interxy*2])==False :
            
            if tdew[tt] == type1 :
                ew+=1
                fcom[0,:,:]+=field[tt,:,:]
            elif tdew[tt] == type2 :
                mon+=1
                fcom[1,:,:]+=field[tt,:,:]
    
    fcom[0,:,:]=fcom[0,:,:]/type1num
    fcom[1,:,:]=fcom[1,:,:]/type2num
        
        
    return fcom

#!!!
#1:850hPa wind, 500hPa geopotential height
#2:850hPa specific humidity
#3:850hPa RH, 700hPa geopotential height
#4:200hPa wind, 200hPa geopotential height
#5:200-850hPa VWS
var=1
    
#!!!
#1:ew
#2:mon
ew=1


LMIR=np.load('LMIR.npy')
tdrmw=np.load('tdrmw.npy')
tdlat=np.load('tdlat_ny.npy')
tdlon=np.load('tdlon_nx.npy')
tdew=np.load('tdew.npy')
tdtimestr=np.load('tdtimestr.npy', allow_pickle=True)
tdhour=np.load('tdhour.npy', allow_pickle=True)
tdyear=np.load('tdyear.npy', allow_pickle=True)

title_num=[['a','b'],['c','d']]

#draw https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.streamplot.html
fig,ax=plt.subplots(2,2,sharex=False,sharey=False,figsize=(18,18))

for var in range(1,3):
    if var==2:
        interxy=20
    else:
        interxy=40 #40degree(+-20degree)
    f1=np.zeros((int(tdtimestr.size),interxy*4+1,interxy*4+1))
    f2=np.zeros((int(tdtimestr.size),interxy*4+1,interxy*4+1))
    f3=np.zeros((int(tdtimestr.size),interxy*4+1,interxy*4+1))
    f4=np.zeros((int(tdtimestr.size),interxy*4+1,interxy*4+1))
    
    """
    
    #tdew=np.zeros(tdtimestr.size)
    j_md=0
    j_else=0
    
    for tt in range(281,431):#259,454
        
        #print(tt)
        
        #sort real timestr
        if testdt(tdtimestr[tt]) == False :
            continue
        
        a=tdhour[tt]
        if a != '00' and a != '03' and a != '06' and a != '09' and a != '12' and a != '15' and a != '18' and a != '21':
            continue
        
        if ew==1:
            type1=-1.5#-1
            type2=-2#-1.5
            notype1=1#1
            notype2=-1#-2
        
        if tdew[tt] == type2:
            tdew[tt] = notype2
            print(tt)
        
        if tt==281 or tt==283 or tt==297 or tt==339 or tt==368 or tt==399\
             or tt==400 or tt==404 or tt==405 or tt==423 or tt==430:
                 tdew[tt] = -2
                 j_md+=1
        else:
            tdew[tt] = -1.5
            j_else+=1
        
        
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
        #vo = ncfile["vo"]
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
            f33=q[0,30,:,:]
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
                if var==1 or var==2:
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
        
    f1com=j_composite(f1,interxy,type1,type2,j_else,j_md,tdew)
    f2com=j_composite(f2,interxy,type1,type2,j_else,j_md,tdew)
    f3com=j_composite(f3,interxy,type1,type2,j_else,j_md,tdew)
    
    #ANOVA P-value
    
    others=np.zeros((int(j_else),interxy*4+1,interxy*4+1))
    md=np.zeros((int(j_md),interxy*4+1,interxy*4+1))
    
    oth=0
    mon=0
    for tt in range(281,431):
        if np.isnan(f1[tt,interxy*2,interxy*2])==False :
            
            if tdew[tt] == type1 :
                if var==1:
                    others[oth,:,:]=f1[tt,:,:]
                if var==2:
                    others[oth,:,:]=f3[tt,:,:]
                oth+=1
            elif tdew[tt] == type2 :
                if var==1:
                    md[mon,:,:]=f1[tt,:,:]
                if var==2:
                    md[mon,:,:]=f3[tt,:,:]
                mon+=1
    
    fpcom=np.zeros((interxy*4+1,interxy*4+1))
    for i in range(interxy*4+1) :
        for j in range(interxy*4+1) :
            st,fpcom[i,j]=stats.f_oneway(others[:,i,j], md[:,i,j])
    
    
    # npy
    np.save('f1com_J_var'+str(var),f1com)
    np.save('f2com_J_var'+str(var),f2com)
    np.save('f3com_J_var'+str(var),f3com)
    np.save('fpcom_J_var'+str(var),fpcom)
    
    """
    
    f1com=np.load('f1com_J_var'+str(var)+'.npy')
    f2com=np.load('f2com_J_var'+str(var)+'.npy')
    f3com=np.load('f3com_J_var'+str(var)+'.npy')
    fpcom=np.load('fpcom_J_var'+str(var)+'.npy')
    #f1com[1,:,:]=f1com[1,:,:]*61/11
    #f2com[1,:,:]=f2com[1,:,:]*61/11
    
    
    
    #print(f3com)
    
    NULAT=np.arange(-int(interxy/2),int(interxy/2)+0.25,0.25)
    NULONG=np.arange(-int(interxy/2),int(interxy/2)+0.25,0.25)
    
    
    
    if var==1:
        title1=''
        title3='850w_'
    if var==2:
        title1=''
        title3='850q_850w_'
    if var==3:
        title1=''
        title3='850r_700z_'
    if var==4:
        title1='200hPa wind, '
        title3='200w_200z_'
    if var==5:
        title1='200-850hPa VWS, '
        title3='200_850vws_'
        
    if ew==1:
        title4='JTWC_MD'
    if ew==2:
        title4='mon'
    
    for hl in range(2):
        
            if ew==1 and hl==1:
                title2='others'
            if ew==2 and hl==0:
                title2='MC'
            if ew==1 and hl==0:
                title2='JTWC_MD'
            if ew==2 and hl==1:
                title2='MS'
    
            if var==2:
                ax[var-1,hl].axis([-10,10,-10,10])
                ax[var-1,hl].set_xticks(np.arange(-10,12.5,2.5))
                ax[var-1,hl].set_yticks(np.arange(-10,12.5,2.5))
            else:
                ax[var-1,hl].axis([-20,20,-20,20])
                ax[var-1,hl].set_xticks(np.arange(-20,25,5))
                ax[var-1,hl].set_yticks(np.arange(-20,25,5))
            ax[var-1,hl].tick_params(labelsize=15)
            ax[var-1,hl].grid(True)
            ax[var-1,hl].set_aspect('equal')
            #plt.axes(projection=ccrs.PlateCarree()).coastlines(resolution='auto', color='k',linewidth=0.5)
            
            #inverse (fku jj)
            lh=-hl+1
    
            if var==1 or var==2 or var==4 or var==5:
                f12com=np.zeros((interxy*4+1,interxy*4+1))
                for i in range(interxy*4+1) :
                    for j in range(interxy*4+1) :
                        f12com[j,i]=sqrt(f1com[lh,j,i]**2+f2com[lh,j,i]**2)
    
            if var==1:
                c=ax[var-1,hl].contourf(NULONG[:],NULAT[:],f12com[:,:],cmap="Reds",levels=np.arange(4,16,2),extent=(-20,-20,20,20),extend='max')
                d=ax[var-1,hl].streamplot(NULONG[:],NULAT[:],f1com[lh,:,:],f2com[lh,:,:],color="b")
                e=ax[var-1,hl].contour(NULONG[:],NULAT[:],fpcom[:,:],levels=np.arange(0.05,0.06,0.01),colors=["#FFFF00"],linewidths=5)
                #ax[var-1,hl].clabel(e, inline=True, fontsize=12)
                cbar_ax = fig.add_axes([0.87, 0.57, 0.02, 0.35]) #x,y,dx,dy #!!!x==0.87 y==0.57
                cbar1=fig.colorbar(c, cax=cbar_ax)
                cbar1.ax.tick_params(labelsize=25)
                cbar1.set_label('$[m \\, s^{-1}]$', fontsize=25)
            
            if var==2:
                c=ax[var-1,hl].contourf(NULONG[:],NULAT[:],f3com[lh,:,:],cmap="Blues",levels=np.arange(0.012,0.014,0.0005),extent=(-10,-10,10,10),extend='both')
                d=ax[var-1,hl].streamplot(NULONG[:],NULAT[:],f1com[lh,:,:],f2com[lh,:,:],color="k")
                e=ax[var-1,hl].contour(NULONG[:],NULAT[:],fpcom[:,:],levels=np.arange(0.05,0.06,0.01),colors=["#FFFF00"],linewidths=5)
                #e=ax[var-1,hl].contour(NULONG[:],NULAT[:],f12com[:,:],levels=np.arange(0,16,4),color="y",linewidths=5)
                #ax[var-1,hl].clabel(e, inline=True, fontsize=15)
                cbar_ax = fig.add_axes([0.87, 0.087, 0.02, 0.35]) #x,y,dx,dy #!!!x==0.87 y==0.087
                cbar2=fig.colorbar(c, cax=cbar_ax)
                cbar2.ax.tick_params(labelsize=25)
                cbar2.set_label('$[kg \\, kg^{-1}]$', fontsize=25)
            
            if var==3:
                c=ax[var-1,hl].contourf(NULONG[:],NULAT[:],f1com[lh,:,:],cmap="Blues",levels=np.arange(50,100,10),extent=(-20,-20,20,20),extend='both')
                #d=ax[var-1,hl].contour(NULONG[:],NULAT[:],f2com[hl,:,:],levels=np.arange(-0.3,-0.1,0.1),color="red",linewidths=1)
                #e=ax[var-1,hl].contour(NULONG[:],NULAT[:],f3com[hl,:,:],levels=np.arange(2820,3420,30),color="g",linewidths=5)
                #ax[var-1,hl].clabel(d, inline=True, fontsize=10)
                #ax[var-1,hl].clabel(e, inline=True, fontsize=12)
                cbar_ax = fig.add_axes([0.87, 0.15, 0.02, 0.7]) #x,y,dx,dy
                cbar=fig.colorbar(c, cax=cbar_ax)
                cbar.ax.tick_params(labelsize=25)
                cbar.set_label('(%)', fontsize=25)
            
            if var==4:
                c=ax[var-1,hl].contourf(NULONG[:],NULAT[:],f12com[:,:],cmap="Reds",levels=np.arange(0,29,4),extent=(-20,-20,20,20),extend='max')
                d=ax[var-1,hl].streamplot(NULONG[:],NULAT[:],f1com[lh,:,:],f2com[lh,:,:],color="b")
                #e=ax[var-1,hl].contour(NULONG[:],NULAT[:],f3com[hl,:,:],levels=np.arange(11760,12600,60),color="g",linewidths=5)
                #ax[var-1,hl].clabel(e, inline=True, fontsize=12)
            
            if var==5:
                c=ax[var-1,hl].contourf(NULONG[:],NULAT[:],f12com[:,:],cmap="Reds",levels=np.arange(0,29,4),extent=(-20,-20,20,20),extend='max')
                d=ax[var-1,hl].streamplot(NULONG[:],NULAT[:],f1com[hl,:,:],f2com[hl,:,:],color="b")
            
            ax[0,hl].set_title(title1+title2, fontsize=25)
            ax[var-1,hl].set_title('('+title_num[var-1][hl]+')',loc='left', fontsize=25)
            ax[1,hl].set_xlabel('$(^o)$', fontsize=25)
            ax[var-1,0].set_ylabel('$(^o)$', fontsize=25)
                

#fig.subplots_adjust(right=0.8)
fig.subplots_adjust(left=0.087, right=0.85, bottom=0.01, top=0.99, hspace=0.01)#!!!
    
fig.savefig('fig/'+title4+'.png', dpi=600)