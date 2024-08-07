#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 20:44:15 2022

@author: cwuhank
"""


import numpy as np
import matplotlib.pyplot as plt

tdew=np.load('tdew_m.npy')#!!!_test
LMIR=np.load('LMIR.npy')
tsrmw=np.load('tsrmw.npy')
tstimestr=np.load('tstimestr.npy', allow_pickle=True)

colors = ['red','purple','green','blue']
labels=['EW','MC','MS','MD']
ewtypes=[1,-1,-1.5,-2]

#!!!
#1:LMIR
#2:tsrmw
var=1

#draw https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.streamplot.html
fig,ax=plt.subplots(2,1,sharex=True,sharey=False,figsize=(7,10))
title_num=['a','b']

for var in range(2):
    
    if var==1:
        drawvar=LMIR
        title1='LMIR '
        title2='$ (kt \\, day^{-1})$'
        
    if var==0:
        drawvar=tsrmw
        title1='TS_RMW'
        title2=' (km)'
    
    num=0
    ewnum=np.zeros(4)
    for tt in range(tstimestr.size):
        if np.isnan(drawvar[tt])==False : 
            num+=1
            for i in range(4):
                if tdew[tt]==ewtypes[i]:
                    ewnum[i]+=1
    
    print(num) #186 375
    draw1=np.zeros(int(ewnum[0]))
    draw2=np.zeros(int(ewnum[1]))
    draw3=np.zeros(int(ewnum[2]))
    draw4=np.zeros(int(ewnum[3]))
    uu1=0
    uu2=0
    uu3=0
    uu4=0
    for tt in range(tstimestr.size):
        if np.isnan(drawvar[tt])==False :
                if tdew[tt]==1:
                    draw1[uu1]=drawvar[tt]
                    uu1+=1
                if tdew[tt]==-1:
                    draw2[uu2]=drawvar[tt]
                    uu2+=1
                if tdew[tt]==-1.5:
                    draw3[uu3]=drawvar[tt]
                    uu3+=1
                if tdew[tt]==-2:
                    draw4[uu4]=drawvar[tt]
                    uu4+=1
    
    draw=[draw1,draw2,draw3,draw4]
    #fig, ax = plt.subplots()
    bplot=ax[var].boxplot(draw,patch_artist=True,labels=labels)
    ax[var].set_ylabel(title1+title2)
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)
    ax[var].yaxis.grid(True)
    ax[var].set_title('('+title_num[var]+')',loc='left', fontsize=25)
    plt.savefig('fig/boxplot.png', dpi=600)