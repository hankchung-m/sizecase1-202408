#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 17:45:17 2024

@author: cwuhank
"""

import numpy as np
import scipy.stats as stats
import scikit_posthocs as sp
import pandas as pd
from datetime import date as dt

def testdt(datestr):
    try:
        dt.fromisoformat(datestr)
    except:
        return False
    else:
        return True
itimestr=np.load('itimestr.npy', allow_pickle=True)
ihour=np.load('ihour.npy', allow_pickle=True)
tstimestr=np.load('tstimestr.npy', allow_pickle=True)
tshour=np.load('tshour.npy', allow_pickle=True)
tdew=np.load('tdew_m.npy')#!!!_test
tdlon=np.load('tdlon.npy')
tddp=np.load('tddp.npy')
tdese=np.load('tdese.npy')
tdsw=np.load('tdsw.npy')
tdrmw_n=np.load('tdrmw_n.npy')
tsrmw=np.load('tsrmw.npy')
LMIR=np.load('LMIR.npy')
env_lf_i=np.load('env_lf_i.npy')
env_lf_s=np.load('env_lf_s.npy')
env_dl=np.load('env_dl.npy')
env_tslat=np.load('env_tslat.npy')
env_tslon=np.load('env_tslon.npy')
env_ivws=np.load('env_ivws.npy')
#print(np.nanmean(LMIR[tsrmw>100]))
ii=0
jj=0
for tt in range(tddp.size):
    if tsrmw[tt]<=100:
        ii+=1
        if LMIR[tt]>=30:
            jj+=1
            
print(ii,jj,jj/ii)



ew=[1,-1,-1.5,-2]
ewname=['ew','mc','ms','md']
for eeww in range(4):
    ii=0
    jj=0
    for tt in range(itimestr.size):
        if testdt(itimestr[tt]) == False :
            continue
        
        a=ihour[tt]
        if a != '00' and a != '03' and a != '06' and a != '09' and a != '12' and a != '15' and a != '18' and a != '21':
            continue
        if env_lf_i[tt]==0 and env_lf_s[tt]==0 and env_dl[tt]==0 \
            and env_tslat[tt]==0 and env_tslon[tt]==0 and env_ivws[tt]==0:
            if tdew[tt]==ew[eeww]:
                ii+=1
                if LMIR[tt]>=30:
                    jj+=1
                
    print(ewname[eeww],ii,jj,jj/ii)