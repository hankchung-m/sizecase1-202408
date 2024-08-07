#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 16:48:52 2022

@author: cwuhank
"""
#https://statistics-using-python.blogspot.com/2019/08/kruskal-wallis-h-testnon-parametric.html

import numpy as np
import scipy.stats as stats
import scikit_posthocs as sp
import pandas as pd

tdew=np.load('tdew_m.npy')#!!!_test
tdlon=np.load('tdlon.npy')
tddp=np.load('tddp.npy')
tdese=np.load('tdese.npy')
tdsw=np.load('tdsw.npy')
tdrmw_n=np.load('tdrmw_n.npy')
tsrmw=np.load('tsrmw.npy')
LMIR=np.load('LMIR.npy')

ewtype=[1,-1,-1.5,-2]

for var in range(6):
    
    if var==0:
        tdraw=tddp
        print('tddp')
    if var==1:
        tdraw=tdese
        print('tdese')
    if var==2:
        tdraw=tdsw
        print('tdsw')
    if var==3:
        tdraw=tdrmw_n
        print('tdrmw_n')
    if var==4:
        tdraw=tsrmw
        print('tsrmw')
    if var==5:
        tdraw=LMIR
        print('LMIR')


    num=0
    for tt in range(tddp.size):
            if tddp[tt]>0 and tddp[tt]<3000 and tdlon[tt]>=121 and tdew[tt]==ewtype[3]:
                num+=1
    
    #print(num) #350 JTWC's problem
        
    x = np.zeros(num)
    
    uu=0
    for tt in range(tddp.size):
            if tddp[tt]>0 and tddp[tt]<3000 and tdlon[tt]>=121 and tdew[tt]==ewtype[3]:
                x[uu]=tdraw[tt]
                
                uu+=1
                if uu==num:
                    break
        
    x3=x
    
    
    num=0
    for tt in range(tddp.size):
            if tddp[tt]>0 and tddp[tt]<3000 and tdlon[tt]>=121 and tdew[tt]!=0:
                num+=1
    
    #print(num) #350 JTWC's problem
        
    x = np.zeros(num)
    y = np.zeros(num)
    
    uu=0
    for tt in range(tddp.size):
            if tddp[tt]>0 and tddp[tt]<3000 and tdlon[tt]>=121 and tdew[tt]!=0:
                x[uu]=tdraw[tt]
                y[uu]=tdew[tt]
                
                uu+=1
                if uu==num:
                    break
        
    x4=x
    y4=y
    
    dat = {'raw':x4,'Type':y4}
    df = pd.DataFrame(dat)
    print(sp.posthoc_conover(df, val_col = 'raw', group_col = 'Type', p_adjust = 'holm'))
    print('------------------------')
    #print('EW+MC+MS:',stats.shapiro(x4))
    #print('MD:',stats.shapiro(x3))
    #print(stats.levene(x4, x3, center = 'mean'))
    #print(stats.kruskal(x4, x3))
    #ST,P=stats.kruskal(x4, x3)
    #if P<0.05:
    #    print('different')
    #else:
    #    print('NOT different')
