#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 17:46:09 2022

@author: cwuhank
"""


import numpy as np
import datetime
from math import*

def Distance(lat1,lng1,lat2,lng2): #https://www.itread01.com/content/1550369361.html
    radlat1=radians(lat1)  
    radlat2=radians(lat2)  
    a=radlat1-radlat2  
    b=radians(lng1)-radians(lng2)  
    s=2*asin(sqrt(pow(sin(a/2),2)+cos(radlat1)*cos(radlat2)*pow(sin(b/2),2)))  
    earth_radius=6378.137  
    s=s*earth_radius  
    if s<0:  
        return -s  
    else:  
        return s

def testdt(datestr):
    try:
        datetime.date.fromisoformat(datestr)
    except:
        return False
    else:
        return True

tslat=np.load('tslat.npy')
tslon=np.load('tslon.npy')
tdlat=np.load('tdlat.npy')
tdlon=np.load('tdlon.npy')
tdew=np.load('tdew.npy')

tstimestr=np.load('tstimestr.npy', allow_pickle=True)
tshour=np.load('tshour.npy', allow_pickle=True)
tdtimestr=np.load('tdtimestr.npy', allow_pickle=True)
tdhour=np.load('tdhour.npy', allow_pickle=True)
tdtttimestr=np.load('tdtttimestr.npy', allow_pickle=True)


name=np.load('name.npy', allow_pickle=True)

utcltimestr=np.loadtxt('UTCL_bst_2000_2021.txt', delimiter=' ', usecols=(0))
utcllat=np.loadtxt('UTCL_bst_2000_2021.txt', delimiter=' ', usecols=(1))
utcllon=np.loadtxt('UTCL_bst_2000_2021.txt', delimiter=' ', usecols=(2))

tdutcl=np.zeros(tdtimestr.size)
#tdutcl_dx=np.zeros(tdtimestr.size)
#tdutcl_dy=np.zeros(tdtimestr.size)
tdutcllat=np.zeros(tdtimestr.size)
tdutcllon=np.zeros(tdtimestr.size)

pretdutcl=np.zeros(tdtimestr.size)
pretdutcllat=np.zeros(tdtimestr.size)
pretdutcllon=np.zeros(tdtimestr.size)
preutcl=np.zeros(tdtimestr.size)
fromutcl=np.zeros(tdtimestr.size)

for tt in range(tdtimestr.size):
    tdutcl[tt]=99999
    pretdutcl[tt]=99999
    #print(tt)
    
    #sort real timestr
    if testdt(tdtimestr[tt]) == False :
        continue
    
    a=tdhour[tt]
    if a != '00' and a != '03' and a != '06' and a != '09' and a != '12' and a != '15' and a != '18' and a != '21':
        continue
    
    dd=99999
    predd=99999
    for x in range(utcltimestr.size):
        
        #>29day continue
        tdtttimestrs = str(tdtttimestr[tt])
        utcltimestrs = str(int(utcltimestr[x]))
        utcltimestrs_check = str(int(utcltimestr[x-4]))
        tutttimestrstr = datetime.datetime.strptime(tdtttimestrs,"%Y%m%d%H")
        utcltimestrstr = datetime.datetime.strptime(utcltimestrs,"%Y%m%d%H")
        utcltimestrstr_check = datetime.datetime.strptime(utcltimestrs_check,"%Y%m%d%H")
        
        check_pre_x=utcltimestrstr_check-utcltimestrstr
        
        intv=utcltimestrstr-tutttimestrstr
        if str(intv) =='29 days, 0:00:00':
            continue
        
        
        #find +-3hr, calculate Distance
        if tdtttimestrs==utcltimestrs:
            dd=Distance(float(tdlat[tt]),float(tdlon[tt]),float(utcllat[x]),float(utcllon[x]))
            if dd<tdutcl[tt] and dd<1700:
                tdutcl[tt]=dd
                tdutcllat[tt]=float(utcllat[x])
                tdutcllon[tt]=float(utcllon[x])
            
        elif str(intv) =='3:00:00' or str(intv) =='-1 day, 21:00:00':
            if tdutcl[tt]==99999:
                dd=Distance(float(tdlat[tt]),float(tdlon[tt]),float(utcllat[x]),float(utcllon[x]))
                if dd<tdutcl[tt] and dd<1700:
                    tdutcl[tt]=dd
                    tdutcllat[tt]=float(utcllat[x])
                    tdutcllon[tt]=float(utcllon[x])
            
        
    
        #find precede TD
        precede1=(str(intv) =='-1 day, 18:00:00' or str(intv) =='-1 day, 15:00:00' or 
                 str(intv) =='-1 day, 12:00:00')  
        precede2=(str(intv) =='-1 day, 9:00:00' or
                 str(intv) =='-1 day, 6:00:00' or str(intv) =='-1 day, 3:00:00' or str(intv) =='-1 day')
        
        precede3=(str(intv) =='-2 day, 21:00:00' or str(intv) =='-2 day, 18:00:00' or str(intv) =='-2 day, 15:00:00' or 
                 str(intv) =='-2 day, 12:00:00')  
        precede4=(str(intv) =='-2 day, 9:00:00' or
                 str(intv) =='-2 day, 6:00:00' or str(intv) =='-2 day, 3:00:00' or str(intv) =='-2 day')
        
        precede5=(str(intv) =='-3 day, 21:00:00' or str(intv) =='-3 day, 18:00:00' or str(intv) =='-3 day, 15:00:00' or 
                 str(intv) =='-3 day, 12:00:00') 
        precede6=(str(intv) =='-3 day, 9:00:00' or
                 str(intv) =='-3 day, 6:00:00' or str(intv) =='-3 day, 3:00:00' or str(intv) =='-3 day')
        
        real_pre_x=(str(check_pre_x) =='-1 day')
        
        #predd=Distance(float(tdlat[tt]),float(tdlon[tt]),float(utcllat[x]),float(utcllon[x]))
        if (precede1 or precede2 or precede3 or precede4 or precede5 or precede6) and real_pre_x:# and predd<1700:
            
            #predict location
            if precede1:
                pre_x=1/4
            if precede2:
                pre_x=3/4
            if precede3:
                pre_x=5/4
            if precede4:
                pre_x=7/4
            if precede5:
                pre_x=9/4
            if precede6:
                pre_x=11/4
            prelat=pre_x*(float(utcllat[x])-float(utcllat[x-4]))+float(utcllat[x])
            prelon=pre_x*(float(utcllon[x])-float(utcllon[x-4]))+float(utcllat[x])
            movdd=Distance(float(utcllat[x]),float(utcllon[x]),prelat,prelon)
            predd=Distance(float(tdlat[tt]),float(tdlon[tt]),prelat,prelon)
            if predd<tdutcl[tt] and predd<1700:
                pretdutcl[tt]=predd
                pretdutcllat[tt]=prelat
                pretdutcllon[tt]=prelon
                
                #genesis from utcl = 1
                if predd<movdd/2:
                    fromutcl[tt]=1
        
    
    
    #replace
    if tdutcl[tt]==99999:
        tdutcl[tt]=pretdutcl[tt]
        tdutcllat[tt]=pretdutcllat[tt]
        tdutcllon[tt]=pretdutcllon[tt]
        preutcl[tt]=1
        
        
    
    """
    #check different
    if pretdutcl[tt]!=99999 and abs(tdutcllon[tt]-pretdutcllon[tt])>3:
        print('UTCL',tdutcllon[tt],tdutcllat[tt])
        print('preUTCL',pretdutcllon[tt],pretdutcllat[tt])
        #precede is near
        if tdutcl[tt]>pretdutcl[tt]:
            tdutcl[tt]=pretdutcl[tt]
            tdutcllat[tt]=pretdutcllat[tt]
            tdutcllon[tt]=pretdutcllon[tt]
    """
    
    if tdutcl[tt]!=99999:
        print(tt)
        print(name[tt])
        print('EW',tdew[tt])
        print(tdtimestr[tt])
        print('TD',tdlon[tt],tdlat[tt])
        print(tdutcl[tt])
        print('UTCL',tdutcllon[tt]-tdlon[tt],tdutcllat[tt]-tdlat[tt])
        print('genesis pre utcl',preutcl[tt])
        print('genesis from utcl',fromutcl[tt])
        print('----------')
        #print('preUTCL',pretdutcllon[tt],pretdutcllat[tt])
    
np.save('tdutcl',tdutcl)
#np.save('tdutcl_dx',tdutcl_dx)
#np.save('tdutcl_dy',tdutcl_dy)
np.save('tdutcllon',tdutcllon)
np.save('tdutcllat',tdutcllat)
np.save('preutcl',preutcl)
np.save('fromutcl',fromutcl)