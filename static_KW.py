#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 12:19:04 2022

@author: cwuhank
"""
#https://statistics-using-python.blogspot.com/2019/08/kruskal-wallis-h-testnon-parametric.html

import numpy as np
import scipy.stats as stats

tdew=np.load('tdew.npy')
tdlon=np.load('tdlon.npy')
tdlat=np.load('tdlat.npy')
tddp=np.load('tddp.npy')
tdese=np.load('tdese.npy')
tdsw=np.load('tdsw.npy')
tdrmw_n=np.load('tdrmw_n.npy')
tsrmw=np.load('tsrmw.npy')
irmw=np.load('irmw.npy')
LMIR=np.load('LMIR.npy')

ewtype=[1,-1,-1.5,-2]

for var in range(2):
    
    if var==0:
        tdraw=tsrmw
        print('tsrmw')
    if var==1:
        tdraw=LMIR
        print('LMIR')
    idealized=1#!!!
    
    for e in range(4):
        num=0
        for tt in range(tddp.size):
            if idealized==1:
                if tddp[tt]>0 and tddp[tt]<3000 and tdlon[tt]>=121 and tdew[tt]==ewtype[e]:
                    num+=1
            if idealized==0:
                if tddp[tt]>0 and tddp[tt]<3000 and tdew[tt]==ewtype[e]:
                    num+=1
    
        #print(num) #350 JTWC's problem
        
        x = np.zeros(num)
    
        uu=0
        for tt in range(tddp.size):
            if idealized==1:
                if tddp[tt]>0 and tddp[tt]<3000 and tdlon[tt]>=121 and tdew[tt]==ewtype[e]:
                    x[uu]=tdraw[tt]
                    
                    uu+=1
                    if uu==num:
                        break
            if idealized==0:
                if tddp[tt]>0 and tddp[tt]<3000 and tdew[tt]==ewtype[e]:
                    x[uu]=tdraw[tt]
                    
                    uu+=1
                    if uu==num:
                        break
        
        if e==0:
            x0=x
        if e==1:
            x1=x
        if e==2:
            x2=x
        if e==3:
            x3=x
        x4=tdraw[~np.isnan(tdraw)]
    
    #print('EW:',stats.shapiro(x0))
    #print('MC:',stats.shapiro(x1))
    #print('MS:',stats.shapiro(x2))
    #print('MD:',stats.shapiro(x3))
    #print(stats.levene(x0, x1, x2, x3, center = 'mean'))
    print(stats.kruskal(x0, x1, x2, x3))
    ST,P=stats.kruskal(x0, x1, x2, x3)
    if P<0.05:
        print('different')
    else:
        print('NOT different')
        
    
    print('------------')
    print('check normal distribution (col. 2>0.05)')
    print('ALL',stats.shapiro(x4))
    print('EW',stats.shapiro(x0))
    print('MC',stats.shapiro(x1))
    print('MS',stats.shapiro(x2))
    print('MD',stats.shapiro(x3))
    print('check equal variances (p>0.05)')
    print(stats.levene(x0, x1, x2, x3, center='mean'))
    print('ANOVA: (p<0.05)')
    print(stats.f_oneway(x2, x3))
    print('------------------------')