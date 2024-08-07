#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 18:17:46 2022

@author: cwuhank
"""


import numpy as np
from datetime import date as dt
from math import*

def ew_composite(field,interxy,type1,type2,wm):
    
    file_wm=['_w','_m']
    tdew=np.load('tdew'+file_wm[wm]+'.npy')
    tdtimestr=np.load('tdtimestr.npy', allow_pickle=True)

    
    fcom=np.zeros((2,interxy*4+1,interxy*4+1)) #(ew/mon , interxy*4+1,interxy*4+1)
    ew=0
    mon=0
        
    for tt in range(tdtimestr.size):
        if np.isnan(field[tt,interxy*2,interxy*2])==False :
            if tdew[tt] == type1 :
                ew+=1
                fcom[0,:,:]+=field[tt,:,:]
            elif tdew[tt] == type2 :
                mon+=1
                fcom[1,:,:]+=field[tt,:,:]
    
    fcom[0,:,:]=fcom[0,:,:]/ew
    fcom[1,:,:]=fcom[1,:,:]/mon
        
        
    return fcom