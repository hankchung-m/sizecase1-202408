#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 21:34:50 2022

@author: cwuhank
"""


import numpy as np

tdew=np.load('tdew_m.npy')
LMIR=np.load('LMIR.npy')
tstimestr=np.load('tstimestr.npy', allow_pickle=True)
env_lf_i=np.load('env_lf_i.npy')
env_lf_s=np.load('env_lf_s.npy')
env_dl=np.load('env_dl.npy')
env_tslat=np.load('env_tslat.npy')
env_tslon=np.load('env_tslon.npy')
env_ivws=np.load('env_ivws.npy')

labels=['ew','mc','ms','md']
ewtypes=[1,-1,-1.5,-2]
ewnum=np.zeros(4)
ewnum_f=np.zeros(4)
uf=['lf_i','lf_s','dl','tslat','tslon','ivws']
ewnum_lf_i=np.zeros(4)
ewnum_lf_s=np.zeros(4)
ewnum_dl=np.zeros(4)
ewnum_tslat=np.zeros(4)
ewnum_tslon=np.zeros(4)
ewnum_ivws=np.zeros(4)

for e in range(4):
    for tt in range(tstimestr.size):
        if tdew[tt]==ewtypes[e]:
            ewnum[e]+=1
            if env_lf_i[tt]==0 and env_lf_s[tt]==0 and env_dl[tt]==0 \
                and env_tslat[tt]==0 and env_tslon[tt]==0 and env_ivws[tt]==0:
                ewnum_f[e]+=1
            
            if env_lf_i[tt]==1:
                ewnum_lf_i[e]+=1
            if env_lf_s[tt]==1:
                ewnum_lf_s[e]+=1
            if env_dl[tt]==1:
                ewnum_dl[e]+=1
            if env_tslat[tt]==1:
                ewnum_tslat[e]+=1
            if env_tslon[tt]==1:
                ewnum_tslon[e]+=1
            if env_ivws[tt]==1:
                ewnum_ivws[e]+=1
                
    ewnum_f[e]=ewnum_f[e]/ewnum[e]
    for u in range(6):
        exec('ewnum_'+uf[u]+'[e]=ewnum_'+uf[u]+'[e]/ewnum[e]')

print('type')
print(labels)
print('number')
print(ewnum)
print('favorable')
print(ewnum_f)
print('unfavorable: lf_i')
print(ewnum_lf_i)
print('stucture change: lf_s')
print(ewnum_lf_s)
print('stucture change: dl')
print(ewnum_dl)
print('unfavorable: tslat')
print(ewnum_tslat)
print('unfavorable: tslon')
print(ewnum_tslon)
print('unfavorable: ivws')
print(ewnum_ivws)