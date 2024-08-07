#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  9 06:01:28 2022

@author: cwuhank
"""


import numpy as np
from datetime import date as dt
import netCDF4 as nc
from math import*
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

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

def DDistance(lat1,lng1,lat2,lng2): #https://www.itread01.com/content/1550369361.html
    s=sqrt(((lat1-lat2)**2)+((lng1-lng2)**2))
    return s

def testdt(datestr):
    try:
        dt.fromisoformat(datestr)
    except:
        return False
    else:
        return True

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

def mddp(tt,phi,ctrx_phi,ctry_phi,XLONG,XLAT):
    
    close_range = 14        # 14 degree square (Weber et al. 2014: 1500km)
    close_range_g = 4*close_range        # 14*4 grids square
    Rdim = close_range_g-1        # radial grid points for azimuthal mean
    inv_r = 0.2                   # radial resolution for azimuthal mean (0.2 grid, 5km) (Weber et al. 2014: 0.5km)
    nthe = 360                  # angle for azimuthal mean
    inv_the = 0.5                 # angle resolution for azimuthal mean (Weber et al. 2014: 0.625)
    cri_phi = 1000000               # criteria for UTCL
    
    #find min_phi
    min_phi=9999
    small_close_range = 4        # 4 degree square (Weber et al. 2014: 1500km)
    small_close_range_g = 4*small_close_range        # 20*4 grids square
    small_Rdim = small_close_range_g-1        # radial grid points for azimuthal mean
    rrxy=small_close_range_g#4 degree (Weber et al. 2014:500km)
    for xx in range(int(ctrx_phi-rrxy),int(ctrx_phi+rrxy)+1):
        for yy in range(int(ctry_phi-rrxy),int(ctry_phi+rrxy)+1):
            if (phi[yy,xx]<=min_phi):
                min_phi=phi[yy,xx]
                #if (rr>=min_rr):
                minx_phi=xx
                miny_phi=yy
                
    analysis_clon=XLONG[minx_phi]
    analysis_clat=XLAT[miny_phi]
    #print(analysis_clon,analysis_clat)
    #print(XLONG[ctrx_phi],XLAT[ctry_phi])
                
    """
    #choose phi_close to transfer coordinate (Weber et al. 2014)
    phi_close = phi[int(miny_phi-small_close_range_g):int(miny_phi+small_close_range_g)+1,\
                                int(minx_phi-small_close_range_g):int(minx_phi+small_close_range_g)+1]
    
    #transfer phi_close (Weber et al. 2014)
    phi_clo_cylin = Car2cylin(small_close_range_g, small_close_range_g, 1, small_Rdim, inv_r, nthe, inv_the, phi_close[np.newaxis,:,:])[0]

    
    #find ICI center (Weber et al. 2014)
    icilon=np.zeros(int(nthe/inv_the))
    icilat=np.zeros(int(nthe/inv_the))
    for the in range(0,int(nthe/inv_the)):
        rthe=the*inv_the

        for rr in range(2,int(small_Rdim/inv_r)-1):
            rrr=(rr-1)*inv_r

            #smoothing
            a=(phi_clo_cylin[the,rr-2]+phi_clo_cylin[the,rr-1]+phi_clo_cylin[the,rr]+phi_clo_cylin[the,rr+1]+phi_clo_cylin[the,rr+2])/5
            if (a>int(min_phi)+3):
                icilon[the]=analysis_clon+rrr/4*cos(radians(rthe))
                icilat[the]=analysis_clat+rrr/4*sin(radians(rthe))
                break
    
    mass_clon=np.mean(icilon)
    mass_clat=np.mean(icilat)
    minx_phi=np.where(abs(XLONG[:]-mass_clon)<=0.125)[0]
    miny_phi=np.where(abs(XLAT[:]-mass_clat)<=0.125)[0]
    """
    
    #again!!! (Weber et al. 2014)
    
    #choose phi_close to transfer coordinate
    phi_close = phi[int(miny_phi-close_range_g):int(miny_phi+close_range_g)+1,\
                                int(minx_phi-close_range_g):int(minx_phi+close_range_g)+1]
    
    #transfer phi_close
    phi_clo_cylin = Car2cylin(close_range_g, close_range_g, 1, Rdim, inv_r, nthe, inv_the, phi_close[np.newaxis,:,:])[0]
    
    #find POCI
    for everyphi in range(int(min_phi)+3,int(min_phi)+363,8): #!!!+1,101,1hPa / +3,363,8gpm
        lowphi=everyphi
        highphi=everyphi+8#!!!+1,+8
        case=np.zeros(int(nthe/inv_the))
        
        for the in range(0,int(nthe/inv_the)):
            rthe=the*inv_the
            lowpass=0
            highpass=0
            countswitch=0
            lowrr=0
            highrr=0
            Rstart=2
            """
            #for TD in land
            if lsm_clo_cylin[the,2]>0:
                Rstart=8#2 degree start
            else:
                Rstart=2#0.5 degree start
            """
            for rr in range(Rstart,int(Rdim/inv_r)-1):
                rrr=rr*inv_r
                """
                #except land effect: dist. covered by land being longer then 5 degree
                if lsm_clo_cylin[the,rr]>0:
                    land+=1
                if land>=5*4:
                    case[the]=0
                    break
                """
                #!!!smoothing (NOT in Weber et al. 2014)
                a=(phi_clo_cylin[the,rr-2]+phi_clo_cylin[the,rr-1]+phi_clo_cylin[the,rr]+phi_clo_cylin[the,rr+1]+phi_clo_cylin[the,rr+2])/5
                    
                #count lower phi contour num.
                if countswitch==0 and a>lowphi:
                    """
                    #except land effect: isobar in land
                    if lsm_clo_cylin[the,rr]>0:
                            case[the]=0
                            break
                    """
                    lowpass+=1
                    countswitch=1
                    lowrr=rrr
                
                
                #two TC seldom near 5 degree
                if countswitch==1 and a<lowphi and rrr>=5*4:#(Weber et al. 2014: 500km)
                    """
                    #except land effect: isobar in land
                    if lsm_clo_cylin[the,rr]>0:
                            case[the]=0
                            break
                    """
                    lowpass+=1
                    countswitch=0
                        
                #count higher phi contour num.
                if a>highphi:
                    """
                    #except land effect: isobar in land
                    if lsm_clo_cylin[the,rr]>0:
                            case[the]=0
                            break
                    """
                    highpass=1
                    highrr=rrr
                
                #else:
                    #case[the]=0
                    #break
                #    land=1
                
                #case1: both lower higher are circle
                #from geostrophic wind(Vg)= 5m/s => Vg=(dphi/dn)/f => dn=(dphi/Vg)/f=9*9.8/5/(5*10**-5)
                #1hPa vs. 326km / 9gpm vs. 353km
                #from large-scale tropical motion => dphi/L=10**-4
                #1hPa vs. 816.3km / 9gpm vs. 882km
                #avg: 3.83 degree / 4.14 degree
                #from gradient wind(V)= 5.25m/s => VV/R+fV=(dphi/dn) => dn=dphi/(fV+(VV/R))
                #dn=9.8*9/(0.00005*5.25+(5.25*5.25/(rr*111000)))
                #ndn=((dn/1000)+882)/2/111
                if lowpass==1 and highpass==1:
                    case[the]=1
                    cri_phi_alter=lowphi
                    break
                
                #case2: lower are two circle, higher wrap lower (NOT in Weber et al. 2014)
                if lowpass>=3 and highpass==1:
                    case[the]=2
                    break
                
                #case3: both lower higher are circle, but higher is far away from lower
                #if lowpass==1 and highpass==1:
                #    case[the]=3
                    #print('pgf too low')
                #    break
                
                #case4: higher doesn't wrap lower
                if lowpass>=1 and highpass==0:
                    case[the]=4
                    
                
                
                #case5: land effect
                #if land==1:
                #    case[the]=-1
                #    break
                
            #case2 or 3 or 4
            if case[the]==2 or case[the]==3 or case[the]==4:
                cri_phi=lowphi
                dp_case=case[the]
                if case[the]==3:
                    print(rthe,rrr,case[the],'pgf too low',highrr,lowrr)
                else:
                    print(rthe,rrr,case[the],cri_phi)
                    #print(phi_clo_cylin[the+1,140-10:140+10])
                break
            
        #case2 or 3 or 4 find cri_phi
        if cri_phi==lowphi:
            break

    #print(cri_phi)

    #if this can't find POCI, pick edge phi
    if cri_phi==1000000:
        cri_phi=cri_phi_alter
        dp_case=1
        
    
    option=3#!!!#(Weber et al. 2014: option3(0))
    #(Weber et al. 2014: option1(1) or option3(0))
    if option==1 or option==3:
        roci_mass_cdiff=1000000
    while roci_mass_cdiff>97:
        #find ROCI
        ddlon=np.zeros(int(nthe/inv_the))
        ddlat=np.zeros(int(nthe/inv_the))
        cri=np.zeros(int(nthe/inv_the))
        for the in range(0,int(nthe/inv_the)):
            rthe=the*inv_the
    
            for rr in range(2,int(Rdim/inv_r)-1):
                rrr=(rr-1)*inv_r
    
                #smoothing
                a=(phi_clo_cylin[the,rr-2]+phi_clo_cylin[the,rr-1]+phi_clo_cylin[the,rr]+phi_clo_cylin[the,rr+1]+phi_clo_cylin[the,rr+2])/5
                if (a>cri_phi):
                    ddlon[the]=analysis_clon+rrr/4*cos(radians(rthe))
                    ddlat[the]=analysis_clat+rrr/4*sin(radians(rthe))
                    cri[the]=Distance(analysis_clat,analysis_clon,ddlat[the],ddlon[the])
                    #if cri[the]>1500:
                    #    print(cri[the],ddlon,ddlat)
                    break
        
        mass_clon=np.mean(ddlon)
        mass_clat=np.mean(ddlat)
        if option==1:
            roci_mass_cdiff=Distance(analysis_clat,analysis_clon,mass_clat,mass_clon)
            if roci_mass_cdiff>97:
                cri_phi=cri_phi-8
        if option==3:
            roci_mass_cdiff=0
        
        
    """
    #(Weber et al. 2014: option1(2))
    ellipticity=0.5
    while ellipticity>0.4:
        #find ROCI
        ddlon=np.zeros(int(nthe/inv_the))
        ddlat=np.zeros(int(nthe/inv_the))
        cri=np.zeros(int(nthe/inv_the))
        for the in range(0,int(nthe/inv_the)):
            rthe=the*inv_the
    
            for rr in range(2,int(Rdim/inv_r)-1):
                rrr=(rr-1)*inv_r
    
                #smoothing
                a=(phi_clo_cylin[the,rr-2]+phi_clo_cylin[the,rr-1]+phi_clo_cylin[the,rr]+phi_clo_cylin[the,rr+1]+phi_clo_cylin[the,rr+2])/5
                if (a>cri_phi):
                    ddlon[the]=analysis_clon+rrr/4*cos(radians(rthe))
                    ddlat[the]=analysis_clat+rrr/4*sin(radians(rthe))
                    cri[the]=Distance(analysis_clat,analysis_clon,ddlat[the],ddlon[the])
                    #if cri[the]>1500:
                    #    print(cri[the],ddlon,ddlat)
                    break
        
        crim=np.sort(cri)
        max10=crim[0:20]#max10degree
        min10=crim[int(nthe/inv_the)-20:int(nthe/inv_the)]#min10degree
        ar=np.mean(max10)
        br=np.mean(min10)
        #ar=np.max(cri)
        #br=np.min(cri)
        #print(ar,br)
        ellipticity=(ar-br)/ar
        if ellipticity>0.4:
            cri_phi=cri_phi-8
    """
    
    #print(cri[2])
    
    #find DOCI
    num=0
    allcri=0
    for i in range(int(nthe/inv_the)):
        if np.isnan(cri[i])==False:#234cases no nan
            num+=1
            allcri+=cri[i]


    if num!=0:
        dp=(allcri/num)*2#=radius*2=diameter
    else:
        dp=np.nan
        print('!!!!!!not found',ctrx_phi)#ctrx_phi=164 means lon=121
    
    return dp,dp_case,cri_phi



tdlat_ny=np.load('tdlat_ny.npy')
tdlon_nx=np.load('tdlon_nx.npy')
tdroci=np.load('tdroci.npy')

tdtimestr=np.load('tdtimestr.npy', allow_pickle=True)
tdhour=np.load('tdhour.npy', allow_pickle=True)
name=np.load('name.npy', allow_pickle=True)
tdyear=np.load('tdyear.npy', allow_pickle=True)
LMIR=np.load('LMIR.npy')

tdew=np.zeros(tdtimestr.size)



tdsw=np.zeros(tdtimestr.size)
tdnw=np.zeros(tdtimestr.size)
tdse=np.zeros(tdtimestr.size)
tdese=np.zeros(tdtimestr.size)
tddp=np.zeros(tdtimestr.size)
tddp_case=np.zeros(tdtimestr.size)
tdpp=np.zeros(tdtimestr.size)
tdwc=np.zeros(tdtimestr.size)
tdsh=np.zeros(tdtimestr.size)
for tt in range(tdtimestr.size):
    
    print(tt)
    
    #sort real timestr
    if testdt(tdtimestr[tt]) == False :
        tdsw[tt]=np.nan
        tddp[tt]=np.nan
        tddp_case[tt]=np.nan
        tdsh[tt]=np.nan
        continue
    
    a=tdhour[tt]
    if a != '00' and a != '03' and a != '06' and a != '09' and a != '12' and a != '15' and a != '18' and a != '21':
        tdsw[tt]=np.nan
        tddp[tt]=np.nan
        tddp_case[tt]=np.nan
        tdsh[tt]=np.nan
        continue

    
    if tt/10<1 :
        ncfile = nc.Dataset("ERA5_TD00%d.nc"%(tt))
        #ncfiles = nc.Dataset("ERA5_SFC_TD00%d.nc"%(tt))
    elif tt/100<1 :
        ncfile = nc.Dataset("ERA5_TD0%d.nc"%(tt))
        #ncfiles = nc.Dataset("ERA5_SFC_TD0%d.nc"%(tt))
    else:
        ncfile = nc.Dataset("ERA5_TD%d.nc"%(tt))
        #ncfiles = nc.Dataset("ERA5_SFC_TD%d.nc"%(tt))
        
    XLAT = ncfile["latitude"]
    XLONG = ncfile["longitude"]

    z = ncfile["z"]
    
    """
    # 製作figure  
    #fig1 = plt.figure()   

    #圖表的設定
    #ax1 = fig1.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    #ax1.gridlines()
    #ax1.coastlines(resolution='50m',linewidth=0.5)
    #ax1.set_xticks([XLONG[tdlon_nx[tt]]-20,XLONG[tdlon_nx[tt]]-15,XLONG[tdlon_nx[tt]]-10,XLONG[tdlon_nx[tt]]-5,XLONG[tdlon_nx[tt]],XLONG[tdlon_nx[tt]]+5,XLONG[tdlon_nx[tt]]+10,XLONG[tdlon_nx[tt]]+15,XLONG[tdlon_nx[tt]]+20], crs=ccrs.PlateCarree())
    #ax1.set_yticks([XLAT[tdlat_ny[tt]]-20,XLAT[tdlat_ny[tt]]-15,XLAT[tdlat_ny[tt]]-10,XLAT[tdlat_ny[tt]]-5,XLAT[tdlat_ny[tt]],XLAT[tdlat_ny[tt]]+5,XLAT[tdlat_ny[tt]]+10,XLAT[tdlat_ny[tt]]+15,XLAT[tdlat_ny[tt]]+20], crs=ccrs.PlateCarree())  
    #lon_formatter = LongitudeFormatter()   
    #lat_formatter = LatitudeFormatter()
    #ax1.xaxis.set_major_formatter(lon_formatter)
    #ax1.yaxis.set_major_formatter(lat_formatter)
    
    #plot z_draw
    z_draw = z[0,30,int(tdlat_ny[tt]-4*20):int(tdlat_ny[tt]+4*20)+1,\
               int(tdlon_nx[tt]-4*20):int(tdlon_nx[tt]+4*20)+1]/9.8
    xy=np.arange(int(4*20*2)+1)
    CS=plt.contour(xy,xy, np.flip(z_draw[:,:],axis=0),levels=np.arange(900,2700,3))
    plt.clabel(CS, inline=1, fontsize=10)
    
    #圖例，標題等 
    #ax1.set_xlim(left=XLONG[tdlon_nx[tt]]-20, right=XLONG[tdlon_nx[tt]]+20)
    #ax1.set_ylim(bottom=XLAT[tdlat_ny[tt]]-20, top=XLAT[tdlat_ny[tt]]+20)
    #ax1.set_title(name[tt]+'('+tdyear[tt]+')')
    plt.title(name[tt]+'('+tdyear[tt]+')')
    if tt/10<1 :
        plt.savefig("fig/roci/00%d.png"%(tt), dpi=300)
    elif tt/100<1 :
        plt.savefig("fig/roci/0%d.png"%(tt), dpi=300)
    else:
        plt.savefig("fig/roci/%d.png"%(tt), dpi=300)
    plt.show()
    """

    #DP choose MSLP, or choose 850hPa to prevent outside moutain
    #print(XLONG[tdlon_nx[tt]],XLAT[tdlat_ny[tt]])
    if XLAT[tdlat_ny[tt]]<=30:
        dp,dp_case,pp=mddp(tt,z[0,30,:,:]/9.8,tdlon_nx[tt],tdlat_ny[tt],XLONG,XLAT)
        #dp=md.mddp(lsm[0,:,:],mslp[0,:,:]/100,tdlon_nx[tt],tdlat_ny[tt])
    print(dp)

    tddp[tt]=dp
    tddp_case[tt]=dp_case
    tdpp[tt]=pp

    if tddp[tt]==0:
        tddp[tt]=np.nan
    if tddp_case[tt]==0:
        tddp_case[tt]=np.nan
    if tdpp[tt]==0:
        tdpp[tt]=np.nan

#!!!
#option1(1),option1(12),option3(0)
np.save('tddp_0',tddp)
np.save('tddp_case_0',tddp_case)
np.save('tdpp_0',tdpp)