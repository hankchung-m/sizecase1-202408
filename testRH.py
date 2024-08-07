#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 18:12:46 2024

@author: cwuhank
"""

import numpy as np
from datetime import date as dt

def testdt(datestr):
    try:
        dt.fromisoformat(datestr)
    except:
        return False
    else:
        return True

irh=np.load('irh.npy')
itimestr=np.load('itimestr.npy', allow_pickle=True)
ihour=np.load('ihour.npy', allow_pickle=True)
tstimestr=np.load('tstimestr.npy', allow_pickle=True)
tshour=np.load('tshour.npy', allow_pickle=True)

env_lf_i=np.load('env_lf_i.npy')
env_lf_s=np.load('env_lf_s.npy')
env_dl=np.load('env_dl.npy')
env_tslat=np.load('env_tslat.npy')
env_tslon=np.load('env_tslon.npy')
env_ivws=np.load('env_ivws.npy')


#!!! 0:itime 1:tstime !!!choose
ii=0
if ii==1:
    itimestr=tstimestr
    ihour=tshour

interxy=8#8degree(+-4degree)

numl=0
num=0
for tt in range(itimestr.size):
    #sort real timestr
    if testdt(itimestr[tt]) == False :
        continue
    
    a=ihour[tt]
    if a != '00' and a != '03' and a != '06' and a != '09' and a != '12' and a != '15' and a != '18' and a != '21':
        continue
    
    if env_lf_i[tt]==0 and env_lf_s[tt]==0 and env_dl[tt]==0 \
        and env_tslat[tt]==0 and env_tslon[tt]==0 and env_ivws[tt]==0:
            
        if irh[tt]<50 and irh[tt]!=0:
            print(tt,irh[tt])
            numl+=1
        if irh[tt]!=0:
            num+=1
print(num)
print(numl)
print(numl/num)