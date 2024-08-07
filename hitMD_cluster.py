#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 20:22:11 2022

@author: cwuhank
"""


import numpy as np

tdew=np.load('tdew_m.npy')#_test
tddp=np.load('tddp.npy')
tdroci=np.load('tdroci.npy')
tdtimestr=np.load('tdtimestr.npy', allow_pickle=True)
name=np.load('name.npy')

allnum=0
num=0
hit=0
mis=0
for tt in range(tdtimestr.size):
    if tt>280 and tt<431 and tdew[tt]!=0:#280 431 258 455
        allnum+=1
        if tdew[tt]==-2 : #!!!
            print(tt,name[tt])
            num+=1
            
            if tt==281 or tt==283 or tt==297 or tt==339 or tt==368 or tt==399\
                 or tt==400 or tt==404 or tt==405 or tt==423 or tt==430:
                     hit+=1
                     
        else:
            if tt==281 or tt==283 or tt==297 or tt==339 or tt==368 or tt==399\
                 or tt==400 or tt==404 or tt==405 or tt==423 or tt==430:
                     mis+=1
        
pre=hit/num
rec=hit/11
f1=2*pre*rec/(pre+rec)
print('num:',num)
print('hit:',hit)
print('mis:',mis)
print('f1:',f1)

"""

for tt in range(tdtimestr.size):
    if tt>=195 and tdroci[tt]>=500:
        print(tt,name[tt])
"""