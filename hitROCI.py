#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 21:44:45 2022

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

def func(x,a,b):
    return a/x+b

tddp=np.load('tddp_0.npy')
#tddp_1_2=np.load('tddp_1_2.npy')
tdroci=np.load('tdroci.npy')
tdlon=np.load('tdlon.npy')
tddp=tddp/2
#tddp_1_2=tddp_1_2/2

num=0
for tt in range(tddp.size):
    if tddp[tt]>0 and tddp[tt]<1500 and tdroci[tt]>0 and tdroci[tt]<1500:# and tdlon[tt]<121 and tdlon[tt]>100:
        num+=1

print(num) #357 JTWC's problem
x = np.zeros(num)
y = np.zeros(num)
mad = np.zeros(num)
bia = np.zeros(num)
uu=0
for tt in range(tddp.size):
    if tddp[tt]>0 and tddp[tt]<1500 and tdroci[tt]>0 and tdroci[tt]<1500:# and tdlon[tt]<121 and tdlon[tt]>100:
        x[uu]=tdroci[tt]
        y[uu]=tddp[tt]
        mad[uu]=abs(tddp[tt]-tdroci[tt])
        bia[uu]=tddp[tt]-tdroci[tt]
        #if tddp[tt]!=tddp_1_2[tt]:
        #    print(tt,tddp[tt],tddp_1_2[tt])
        uu+=1
        if uu==num:
            break

#x=tdroci[~np.isnan(tdroci)]
#y=tddp[~np.isnan(tddp)]
#x=1/x#np.log(x)
#inum,icor=reg(x,y)
r = np.corrcoef(x, y)
a, b = np.polyfit(x, y, 1)
#inum,icor=curve_fit(func, x, y, (700, 20))
#ilinex=np.arange(1500)
#iliney=inum[0]/(ilinex)+inum[1]

# 製作figure  
fig = plt.figure()   

#圖表的設定
ax = fig.add_subplot(1, 1, 1, aspect='equal')

#散佈圖
ax.scatter(x, y, color='red', alpha=0.8, linewidths=0, s=5)
plt.plot(x, a*x+b)
#x.plot(ilinex,iliney, color='red')
ax.text(250,1400,'correlation='+str(np.around(r[0,1],decimals=4)))
ax.text(250,1300,'mean abs. diff.='+str(np.around(np.mean(mad),decimals=4)))
ax.text(250,1200,'biase='+str(np.around(np.mean(bia),decimals=4)))

#圖例，標題等
#ax.text(155,90,'y='+str(np.around(inum[0],decimals=4))+'/x+'+str(np.around(inum[1],decimals=4)))
#ax2.text(155,85,'correlation='+str(np.around(icor,decimals=4)))
ax.set_xticks([250,500,750,1000,1250,1500])
ax.set_yticks([250,500,750,1000,1250,1500])  
#ax.tick_params(labelsize=15)
ax.grid(True)
ax.set_xlim(left=0, right=1500)
ax.set_ylim(bottom=0, top=1500)
ax.set_xlabel('JTWC ROCI (km)')
ax.set_ylabel('ERA5 ROCI (km)')
plt.savefig('fig/ROCI_JTWC_ERA5.png', dpi=300)
#plt.savefig('fig/ROCI_JTWC_ERA5_SCS.png')
plt.show() 