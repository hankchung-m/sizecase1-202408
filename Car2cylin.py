#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 29 19:44:55 2022

@author: cwuhank
"""


import numpy as np
from math import*

# define common function
def Car2cylin_wei(the, rr, i_ctr_tmp, j_ctr_tmp, var):
        """  Transform variable from cartision to cylindrical coordinate
             (using 4 points weighting method, computing whole column directly)
        """
        #print(var.shape)
        angle = 2*pi*the/360
        i_float = rr*np.cos(angle)
        j_float = rr*np.sin(angle)
        i_int = int(i_float+i_ctr_tmp)  # round down
        j_int = int(j_float+j_ctr_tmp)  # round down
        di=(i_float+i_ctr_tmp)-i_int
        dj=(j_float+j_ctr_tmp)-j_int
        if di*dj==0:
            w1 = 1
            w2 = 0
            w3 = 0
            w4 = 0
        else:
            w1 = 1/(di*dj)
            w2 = 1/((1-di)*dj)
            w3 = 1/(di*(1-dj))
            w4 = 1/((1-di)*(1-dj))   
        var_cylin = (var[:,j_int,i_int]*w1 + var[:,j_int,i_int+1]*w2 + var[:,j_int+1,i_int]*w3 + var[:,j_int+1,i_int+1]*w4)/(w1+w2+w3+w4)
        return var_cylin
    
def Car2cylin(tcx, tcy, Zdim, Rdim, inv_r, nthe, inv_the, var):
        """  Transform var matrix from cartision to cylindrical coordinate
        """
        var_cylin = np.zeros((Zdim,int(nthe/inv_the),int(Rdim/inv_r)+1))
        for rr in range(0, int(Rdim/inv_r)+1):
            rrr=rr*inv_r
            for the in range(0,int(nthe/inv_the)):
                rthe=the*inv_the
                var_cylin[:,the,rr] = Car2cylin_wei(rthe, rrr, tcx, tcy, var) 
        return var_cylin