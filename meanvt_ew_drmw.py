#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 14:28:45 2022

@author: cwuhank
"""


import numpy as np
from datetime import date as dt
from math import*
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

def testdt(datestr):
    try:
        dt.fromisoformat(datestr)
    except:
        return False
    else:
        return True

def diff(fl2,fl1,fr1,fr2):
    a=4/3*((fr1-fl1)/50)
    b=1/3*((fr2-fl2)/100)
    dfdt=a-b
    return dfdt

def diffedge(fl,f):
    dfdt=(f-fl)/25
    return dfdt

def diff2(fl2,fl1,f,fr1,fr2):
    dfdt1=(f-fl2)/50
    dfdt2=(fr2-f)/50
    ddfdtt=(dfdt2-dfdt1)/50
    return ddfdtt

def diffedge2(fl2,fl1,f):
    dfdt1=(fl1-fl2)/50
    dfdt2=(f-fl1)/50
    ddfdtt=(dfdt2-dfdt1)/50
    return ddfdtt

LMIR=np.load('LMIR.npy')
tstimestr=np.load('tdtimestr.npy', allow_pickle=True)
tshour=np.load('tdhour.npy', allow_pickle=True)
tsyear=np.load('tdyear.npy', allow_pickle=True)
tsmonth=np.load('tdmonth.npy')
tsvt=np.load('tdvt_azi.npy')
tsvr=np.load('tdvr_azi.npy')
tsvo=np.load('tdvo_azi.npy')
tsrh=np.load('tdrh_azi.npy')
tdvtm=np.load('tdvtm.npy')
tdrmw_n=np.load('tdrmw_n.npy')
tdew=np.load('tdew_m.npy')
name=np.load('name.npy')

#!!!
#1:vt1000
#2:vt
#3:vo
#4:rh
#5:qv
#6:vtt
var=6


#1:ew
#2:mon
#ew=1
#for ew in range(1,3):
if var==1:
    tdraw=np.load('tdvt_azi.npy')
    title1='Vt'
    title1_1=' (m/s)'
    title3='tdvt1000_'
    rnum=40
    deg=np.arange(0,1000,25)
    
if var==2:
    tdraw=np.load('tdvt_azi.npy')[:,1:21]
    title1='Vt'
    title1_1=' (m/s)'
    title3='tdvt_'
    rnum=20
    deg=np.arange(25,525,25)
    
if var==3:
    tdraw=np.load('tdvo_azi.npy')[:,1:21]
    title1='Vorticity'
    title1_1=' (s-1)'
    title3='tdvo_'
    rnum=20
    deg=np.arange(25,525,25)
    
if var==5:
    tdraw=np.load('tdqv_azi.npy')[:,1:21]
    title1='Specific Humidity'
    title1_1=' (g/kg)'
    title3='tdqv_'
    rnum=20
    deg=np.arange(25,525,25)
    
if var==6:
    vt=np.load('tdvt_azi.npy')[:,1:21]
    vr=np.load('tdvr_azi.npy')[:,1:21]
    vo=np.load('tdvo_azi.npy')[:,1:21]
    vtt=(-1)*vr*vo
    title1='Vt tendency'
    title1_1=' (ms-1d-1)'
    title3='tdvtt_'
    rnum=20
    deg=np.arange(25,525,25)
    


tdraw_azi_ew=np.zeros((2,tstimestr.size,rnum))
tdraw_azi_q=np.zeros((2,3,rnum))

fig,ax=plt.subplots(1,1,sharex=False,sharey=True,figsize=(9,8))
    
ax.axis([0,500,0,500])
ax.tick_params(labelsize=15)
ax.grid(True)
        
n1=0
n2=0
n3=0
n4=0
ew=0
md=0
mc=0
ms=0
ewrmw=np.zeros(tstimestr.size)
ewdrmw=np.zeros(tstimestr.size)
mcrmw=np.zeros(tstimestr.size)
mcdrmw=np.zeros(tstimestr.size)
msrmw=np.zeros(tstimestr.size)
msdrmw=np.zeros(tstimestr.size)
mdrmw=np.zeros(tstimestr.size)
mddrmw=np.zeros(tstimestr.size)
gen=['EW','MD','MC','MS']
for tt in range(tstimestr.size):

    #sort real timestr
    if testdt(tstimestr[tt]) == False :
        ewrmw[tt]=np.nan
        ewdrmw[tt]=np.nan
        mdrmw[tt]=np.nan
        mddrmw[tt]=np.nan
        continue

    a=tshour[tt]
    if a != '00' and a != '03' and a != '06' and a != '09' and a != '12' and a != '15' and a != '18' and a != '21':
        ewrmw[tt]=np.nan
        ewdrmw[tt]=np.nan
        mdrmw[tt]=np.nan
        mddrmw[tt]=np.nan
        continue

    #print(np.where(vt[tt]==np.nanmax(vt[tt]))[0])
    rmwloc=np.where(vt[tt]==np.nanmax(vt[tt]))[0][0]
    if rmwloc>0:
        rmw=(rmwloc)*25
    else:
        rmw=np.nan
        print(tt,name[tt],tdew[tt])
    #rmw=tdrmw_n[tt]
    #print(rmw)
    if rmwloc<=17:
        a=diff(vtt[tt,rmwloc-2], vtt[tt,rmwloc-1], vtt[tt,rmwloc+1], vtt[tt,rmwloc+2])
        b=diff2(vt[tt,rmwloc-2], vt[tt,rmwloc-1], vt[tt,rmwloc], vt[tt,rmwloc+1], vt[tt,rmwloc+2])
    else:
        a=diffedge(vtt[tt,rmwloc-1], vtt[tt,rmwloc])
        b=diffedge2(vt[tt,rmwloc-2], vt[tt,rmwloc-1], vt[tt,rmwloc])
    drmw=a/b*86400
    
    if tdew[tt] == 1 :
        ew+=1
        ewrmw[tt]=rmw
        ewdrmw[tt]=drmw
        
        if n1==1 :
            ax.scatter(ewrmw[tt],ewdrmw[tt],50,'red',alpha=0.5)
        else:
            ax.scatter(ewrmw[tt],ewdrmw[tt],50,'red',alpha=0.5,label=gen[0])
            n1=1
        
    else:
        ewrmw[tt]=np.nan
        ewdrmw[tt]=np.nan
    """
    if tdew[tt] == -1 :
        mc+=1
        mcrmw[tt]=rmw
        mcdrmw[tt]=drmw
        
        if n2==1 :
            ax.scatter(mcrmw[tt],mcdrmw[tt],50,'purple',alpha=0.5)
        else:
            ax.scatter(mcrmw[tt],mcdrmw[tt],50,'purple',alpha=0.5,label=gen[2])
            n2=1
    
    if tdew[tt] == -1.5 :
        ms+=1
        msrmw[tt]=rmw
        msdrmw[tt]=drmw
        
        if n3==1 :
            ax.scatter(msrmw[tt],msdrmw[tt],50,'green',alpha=0.5)
        else:
            ax.scatter(msrmw[tt],msdrmw[tt],50,'green',alpha=0.5,label=gen[3])
            n3=1
    """
    if tdew[tt] == -2 :
        md+=1
        mdrmw[tt]=rmw
        mddrmw[tt]=drmw
        
        if n4==1 :
            ax.scatter(mdrmw[tt],mddrmw[tt],50,'blue',alpha=0.5)
        else:
            ax.scatter(mdrmw[tt],mddrmw[tt],50,'blue',alpha=0.5,label=gen[1])
            n4=1
    else:
        mdrmw[tt]=np.nan
        mddrmw[tt]=np.nan
    

ewq2=np.nanmean(ewrmw)
ewq2h=450#np.nanmean(ewdrmw)  
ewq3=ewq2+np.nanstd(ewrmw)
ewq1=ewq2-np.nanstd(ewrmw)
ewm=np.nanpercentile(ewrmw,50)

mdq2=np.nanmean(mdrmw)
mdq2h=400#np.nanmean(mddrmw)
mdq3=mdq2+np.nanstd(mdrmw)
mdq1=mdq2-np.nanstd(mdrmw)
mdm=np.nanpercentile(mdrmw,50)

"""
ax.plot([ewq1,ewq3],[ewq2h,ewq2h],color="red",linewidth=5)
ax.scatter(ewq2,ewq2h,100,'red')
ax.plot([mdq1,mdq3],[mdq2h,mdq2h],color="blue",linewidth=5)
ax.scatter(mdq2,mdq2h,100,'blue')
                
#linestyle=[':','--',':']
    
ax.legend(prop={'size': 15})#loc='upper right'
#ax.set_title(title1, fontsize=15)
ax.set_xlabel('RMW (km)', fontsize=15)
ax.set_ylabel('dRMW (km d-1)', fontsize=15)
fig.subplots_adjust(right=0.8)
fig.savefig('fig/rmw_drmw.png', dpi=300)
"""

#corelation of median & mean
#RH:0, qv:1
varcom=0

f1com_var3_1=np.load('f1com_m_var3_ew1_i.npy')
##f2com_var1_1=np.load('f2com_m_var1_ew1_i.npy')
f3com_var6_1=np.load('f3com_m_var6_ew1_i.npy')
##f3com_var2_1=np.load('f3com_m_var2_ew1_i.npy')
f1_var3_1=np.load('f1_m_var3_ew1_i.npy')
##f2_var1_1=np.load('f2_m_var1_ew1_i.npy')
f3_var6_1=np.load('f3_m_var6_ew1_i.npy')
##f3_var2_1=np.load('f3_m_var2_ew1_i.npy')

f1com_var3_2=np.load('f1com_m_var3_ew2_i.npy')
f3com_var6_2=np.load('f3com_m_var6_ew2_i.npy')
f1_var3_2=np.load('f1_m_var3_ew2_i.npy')
f3_var6_2=np.load('f3_m_var6_ew2_i.npy')

f3com_var_var6_1=np.load('f3com_m_var_var6_ew1_i.npy')
f3com_var_var6_2=np.load('f3com_m_var_var6_ew2_i.npy')

f1com_var_var8_1=np.load('f1com_m_var_var8_ew1_i.npy')
f1com_var_var8_2=np.load('f1com_m_var_var8_ew2_i.npy')

#check Quasi-idealized background field (2011-2019, Jul-Aug)
for tt in range(tstimestr.size):#119,tstimestr.size
    
    if ewrmw[tt]==ewm and (tsmonth[tt]>=8 and tsmonth[tt]<=9):
        #wind, sh
        #com_var1=np.hstack((f1com_var1[0,40:121,40:121],f2com_var1[0,40:121,40:121],f3com_var6[0,:,:]*1000,f3com_var6[0,:,:]*1000))#,f2com_var1[0,40:121,40:121],f3com_var6[0,:,:]
        #tt_var1=np.hstack((f1_var1[tt,40:121,40:121],f2_var1[tt,40:121,40:121],f3_var6[tt,:,:]*1000,f3_var6[tt,:,:]*1000))#,f2_var1[tt,40:121,40:121],f3_var6[tt,:,:]
        if varcom==0:
            com_var1=np.hstack(f1com_var3_1[0,40:121,40:121])#*100000)
            tt_var1=np.hstack(f1_var3_1[tt,40:121,40:121])#*100000)
            std_var1=np.hstack(f1com_var_var8_1[0,40:121,40:121])#*100000)
        else:
            com_var1=np.hstack(f3com_var6_1[0,:,:]*100000)
            tt_var1=np.hstack(f3_var6_1[tt,:,:]*100000)
            std_var1=np.hstack(f3com_var_var6_1[0,:,:]*100000)
        err=(tt_var1-com_var1)/np.sqrt(std_var1)
        square_err=np.power(err,2)
        mean_square_err=np.mean(square_err)
        if varcom==0:
            print('EW',tt,tsyear[tt],name[tt],mean_square_err)
        else:
            print('EW',tt,tsyear[tt],name[tt],np.round(mean_square_err))
        
    if mdrmw[tt]==mdm and (tsmonth[tt]>=8 and tsmonth[tt]<=9):
        #wind, sh
        #com_var1=np.hstack((f1com_var1[1,40:121,40:121],f2com_var1[1,40:121,40:121],f3com_var6[1,:,:]*1000,f3com_var6[1,:,:]*1000))#,f2com_var1[1,40:121,40:121],f3com_var6[1,:,:]
        #tt_var1=np.hstack((f1_var1[tt,40:121,40:121],f2_var1[tt,40:121,40:121],f3_var6[tt,:,:]*1000,f3_var6[tt,:,:]*1000))#,f2_var1[tt,40:121,40:121],f3_var6[tt,:,:]
        if varcom==0:
            com_var1=np.hstack(f1com_var3_2[1,40:121,40:121])#*100000)
            tt_var1=np.hstack(f1_var3_2[tt,40:121,40:121])#*100000)
            std_var1=np.hstack(f1com_var_var8_2[1,40:121,40:121])#*100000)
        else:
            com_var1=np.hstack(f3com_var6_2[1,:,:]*100000)
            tt_var1=np.hstack(f3_var6_2[tt,:,:]*100000)
            std_var1=np.hstack(f3com_var_var6_2[1,:,:]*100000)
        err=(tt_var1-com_var1)/np.sqrt(std_var1)
        square_err=np.power(err,2)
        mean_square_err=np.mean(square_err)
        if varcom==0:
            print('MD',tt,tsyear[tt],name[tt],mean_square_err)
        else:
            print('MD',tt,tsyear[tt],name[tt],np.round(mean_square_err))