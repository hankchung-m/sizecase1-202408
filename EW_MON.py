# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 23:33:17 2022

@author: hank0
"""

import numpy as np
from datetime import date as dt
import netCDF4 as nc
from math import*
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import md
import gap

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

def naive_sharding(ds, k):
    """
    Create cluster centroids using deterministic naive sharding algorithm.
    
    Parameters
    ----------
    ds : numpy array
        The dataset to be used for centroid initialization.
    k : int
        The desired number of clusters for which centroids are required.
    Returns
    -------
    centroids : numpy array
        Collection of k centroids as a numpy array.
    """
    
    n = np.shape(ds)[1]
    m = np.shape(ds)[0]
    centroids = np.mat(np.zeros((k,n)))

    # Sum all elements of each row, add as col to original dataset, sort
    composite = np.mat(np.sum(ds, axis=1))
    ds = np.append(composite.T, ds, axis=1)
    ds.sort(axis=0)

    # Step value for dataset sharding
    step = floor(m/k)

    # Vectorize mean ufunc for numpy array
    vfunc = np.vectorize(_get_mean)

    # Divide matrix rows equally by k-1 (so that there are k matrix shards)
    # Sum columns of shards, get means; these columnar means are centroids
    for j in range(k):
        if j == k-1:
            centroids[j:] = vfunc(np.sum(ds[j*step:,1:], axis=0), step)
        else:
            centroids[j:] = vfunc(np.sum(ds[j*step:(j+1)*step,1:], axis=0), step)

    return centroids

def _get_mean(sums, step):
    """
    Vectorizable ufunc for getting means of summed shard columns.
    
    Parameters
    ----------
    sums : float
        The summed shard columns.
    step : int
        The number of instances per shard.
    Returns
    -------
    sums/step (means) : numpy array
        The means of the shard columns.
    """

    return sums/step

tdlat_ny=np.load('tdlat_ny.npy')
tdlon_nx=np.load('tdlon_nx.npy')
tdroci=np.load('tdroci.npy')

tdtimestr=np.load('tdtimestr.npy', allow_pickle=True)
tdhour=np.load('tdhour.npy', allow_pickle=True)
name=np.load('name.npy', allow_pickle=True)
LMIR=np.load('LMIR.npy')

tdew=np.zeros(tdtimestr.size)

tdsws=np.zeros(tdtimestr.size)
tdses=np.zeros(tdtimestr.size)
tdeses=np.zeros(tdtimestr.size)

tdwsw=np.zeros(tdtimestr.size)
tdnn=np.zeros(tdtimestr.size)

tdsw=np.zeros(tdtimestr.size)
tdnw=np.zeros(tdtimestr.size)
tdne=np.zeros(tdtimestr.size)
tdse=np.zeros(tdtimestr.size)
tdese=np.zeros(tdtimestr.size)
tddp=np.zeros(tdtimestr.size)
tdvo=np.zeros(tdtimestr.size)
tdvoe=np.zeros(tdtimestr.size)
tdwc=np.zeros(tdtimestr.size)
tdsh=np.zeros(tdtimestr.size)
tdqv=np.zeros(tdtimestr.size)
tdrh=np.zeros(tdtimestr.size)
for tt in range(tdtimestr.size):
    
    #print(tt)
    
    #sort real timestr
    if testdt(tdtimestr[tt]) == False :
        tdsw[tt]=np.nan
        tddp[tt]=np.nan
        tdvo[tt]=np.nan
        tdvoe[tt]=np.nan
        tdsh[tt]=np.nan
        tdqv[tt]=np.nan
        tdrh[tt]=np.nan
        continue
    
    a=tdhour[tt]
    if a != '00' and a != '03' and a != '06' and a != '09' and a != '12' and a != '15' and a != '18' and a != '21':
        tdsw[tt]=np.nan
        tddp[tt]=np.nan
        tdvo[tt]=np.nan
        tdvoe[tt]=np.nan
        tdsh[tt]=np.nan
        tdqv[tt]=np.nan
        tdrh[tt]=np.nan
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
    
    #p = ncfile["level"]
    u = ncfile["u"]
    v = ncfile["v"]
    vo = ncfile["vo"]
    z = ncfile["z"]
    """
    qv = ncfile["q"]
    rh = ncfile["r"]
    """
    #mslp = ncfiles["msl"]
    #lsm = ncfiles["lsm"]
    #print(u.shape)
    
    pointNN=0
    pointWSW=0
    
    pointNE=0
    pointNW=0
    pointSW=0
    pointSE=0
    pointESE=0
    pointVO=0
    pointVOE=0
    pointWC=0
    pointWCall=0
    pointQV=0
    pointRH=0
    #pointSH=0
    
    uNNall=0
    uWSWall=0
    
    uNEall=0
    uNWall=0
    uSWall=0
    uSEall=0
    uESEall=0
    wSHall=0
    mdp=0
    QVall=0
    RHall=0
    
    for i in range(XLONG.size) :
        for j in range(XLAT.size) :
            
            #ESES
            if abs(XLAT[j]-XLAT[tdlat_ny[tt]])<=10 and abs(XLONG[i]-XLONG[tdlon_nx[tt]])<=10:
                if XLAT[j]<XLAT[tdlat_ny[tt]]-5 and XLONG[i]>XLONG[tdlon_nx[tt]]+5:
                    pointESE+=1
                    uESE=u[0,30,j,i]
                    uESEall+=uESE
                    
            if abs(XLAT[j]-XLAT[tdlat_ny[tt]])<=10 and abs(XLONG[i]-XLONG[tdlon_nx[tt]])<=5:
            #SES
                if XLAT[j]<XLAT[tdlat_ny[tt]]-5 and XLONG[i]>XLONG[tdlon_nx[tt]]:
                    pointSE+=1
                    uSE=u[0,30,j,i]
                    uSEall+=uSE
                    
            #SWS
                if XLAT[j]<XLAT[tdlat_ny[tt]]-5 and XLONG[i]<XLONG[tdlon_nx[tt]]:
                    pointSW+=1
                    uSW=u[0,30,j,i]
                    uSWall+=uSW
            
            """
            
            #NN
            if abs(XLAT[j]-XLAT[tdlat_ny[tt]])<7.5 and abs(XLONG[i]-XLONG[tdlon_nx[tt]])<=2.5:
                if XLAT[j]>XLAT[tdlat_ny[tt]]+2.5:
                    pointNN+=1
                    uNN=u[0,30,j,i]
                    uNNall+=uNN
                    
            #WSW
            if abs(XLAT[j]-XLAT[tdlat_ny[tt]])<=10 and abs(XLONG[i]-XLONG[tdlon_nx[tt]])<=10:
                if XLAT[j]<XLAT[tdlat_ny[tt]]-5 and XLONG[i]<XLONG[tdlon_nx[tt]]-5:
                    pointWSW+=1
                    uWSW=u[0,30,j,i]
                    uWSWall+=uWSW
            
            
            
            #ESE
            if abs(XLAT[j]-XLAT[tdlat_ny[tt]])<=5 and abs(XLONG[i]-XLONG[tdlon_nx[tt]])<=10:
                if XLAT[j]<XLAT[tdlat_ny[tt]] and XLONG[i]>XLONG[tdlon_nx[tt]]+5:
                    pointESE+=1
                    uESE=u[0,30,j,i]
                    uESEall+=uESE
            
            if abs(XLAT[j]-XLAT[tdlat_ny[tt]])<=5 and abs(XLONG[i]-XLONG[tdlon_nx[tt]])<=5:
                
                #SW
                if XLAT[j]<XLAT[tdlat_ny[tt]] and XLONG[i]<XLONG[tdlon_nx[tt]]:
                    pointSW+=1
                    uSW=u[0,30,j,i]
                    uSWall+=uSW
                    
                #NW
                if XLAT[j]>XLAT[tdlat_ny[tt]] and XLONG[i]<XLONG[tdlon_nx[tt]]:
                    pointNW+=1
                    uNW=u[0,30,j,i]
                    uNWall+=uNW
                    
                #NE
                if XLAT[j]>XLAT[tdlat_ny[tt]] and XLONG[i]>XLONG[tdlon_nx[tt]]:
                    pointNE+=1
                    uNE=u[0,30,j,i]
                    uNEall+=uNE
                    
                #SE
                if XLAT[j]<XLAT[tdlat_ny[tt]] and XLONG[i]>XLONG[tdlon_nx[tt]]:
                    pointSE+=1
                    uSE=u[0,30,j,i]
                    uSEall+=uSE
                    
                #SH
                #if XLONG[i]==XLONG[tdlon_nx[tt]]:
                #    pointSH+=1
                #    SH=(u[0,30,j,i]-u[0,30,j-1,i])/111190 #https://m4.hhlink.com/%E7%B6%93%E7%B7%AF%E5%BA%A6
                #    weight=3-(abs(XLAT[j]-XLAT[tdlat_ny[tt]])/1.25)
                #    wSH=weight*SH
                #    wSHall+=wSH
                    
            #if abs(XLAT[j]-XLAT[tdlat_ny[tt]])<=10 and abs(XLONG[i]-XLONG[tdlon_nx[tt]])<=20:
                
                #VO
            #    if XLONG[i]<=XLONG[tdlon_nx[tt]] and vo[0,30,j,i]>0.000004:
            #        pointVO+=1
                    
                #VOE
            #    if XLONG[i]>XLONG[tdlon_nx[tt]] and vo[0,30,j,i]>0.000004:
            #        pointVOE+=1
            
            #WC
            #if DDistance(XLAT[tdlat_ny[tt]],XLONG[tdlon_nx[tt]],XLAT[j],XLONG[i])<=2:
            #    pointWCall+=1
            #    ws=sqrt(u[0,30,j,i]**2+v[0,30,j,i]**2)
            #    if ws<=4:
            #        pointWC+=1
                    
            
            
            #QV
            if abs(XLAT[j]-XLAT[tdlat_ny[tt]])<=5 and abs(XLONG[i]-XLONG[tdlon_nx[tt]])<=5 and abs(XLAT[j]-XLAT[tdlat_ny[tt]])>2.5 and abs(XLONG[i]-XLONG[tdlon_nx[tt]])>2.5:
                pointQV+=1
                QVall+=qv[0,30,j,i]
            
            
            #RH
            if abs(XLAT[j]-XLAT[tdlat_ny[tt]])<=5 and abs(XLONG[i]-XLONG[tdlon_nx[tt]])<=5 and abs(XLAT[j]-XLAT[tdlat_ny[tt]])>2.5 and abs(XLONG[i]-XLONG[tdlon_nx[tt]])>2.5:
                pointRH+=1
                RHall+=rh[0,30,j,i]
               
            """
            
            #DP
            #if DDistance(XLAT[tdlat_ny[tt]],XLONG[tdlon_nx[tt]],XLAT[j],XLONG[i])<7.75 and DDistance(XLAT[tdlat_ny[tt]],XLONG[tdlon_nx[tt]],XLAT[j],XLONG[i])>7.25 :
            #    dp= z[0,30,j,i]-z[0,30, tdlat_ny[tt], tdlon_nx[tt]]
            #    if dp>mdp :
            #        mdp=dp
                    
    
    #DP choose MSLP, or choose 850hPa to prevent outside moutain
    #if XLAT[tdlat_ny[tt]]<=30:
    #    dp=md.mddp(z[0,30,:,:]/9.8,tdlon_nx[tt],tdlat_ny[tt])
        #dp=md.mddp(lsm[0,:,:],mslp[0,:,:]/100,tdlon_nx[tt],tdlat_ny[tt])
    #print(dp)
    
    if pointSE!=0 :
        uSEavg=uSEall/pointSE
    if pointESE!=0 :
        uESEavg=uESEall/pointESE
    if pointSW!=0 :
        uSWavg=uSWall/pointSW
    """
    if pointNN!=0 :
        uNNavg=uNNall/pointNN
    if pointWSW!=0 :
        uWSWavg=uWSWall/pointWSW
    
    if pointNE!=0 :
        uNEavg=uNEall/pointNE
    if pointNW!=0 :
        uNWavg=uNWall/pointNW
    #if pointWCall!=0 :
        #wSHavg=wSHall/pointSH
    
    if pointQV!=0 :
        QVavg=QVall/pointQV
    else:
        QVavg=np.nan
        
    if pointRH!=0 :
        RHavg=RHall/pointRH
    else:
        RHavg=np.nan
        
    tdqv[tt]=QVavg
    tdrh[tt]=RHavg
    """
    
    #if uSWavg==0:
    #    uSWavg=np.nan
    #if wSHavg==0:
    #    wSHavg=np.nan
    

    
    tdsws[tt]=uSWavg
    
    tdses[tt]=uSEavg
    
    tdeses[tt]=uESEavg
    
    """
    
    tdnn[tt]=uNNavg
    
    tdwsw[tt]=uWSWavg    
    
    
    
    pointALL=(20*4+1)**2
        
    tdsw[tt]=uSWavg
    
    tdnw[tt]=uNWavg
    
    tdne[tt]=uNEavg
    
    tdse[tt]=uSEavg
    
    tdese[tt]=uESEavg

    #tddp[tt]=dp
    
    tdvo[tt]=pointVO/pointALL*100
    
    tdvoe[tt]=pointVOE/pointALL*100
    
    #tdwc[tt]=pointWC/pointWCall*100
    
    #tdsh[tt]=wSHavg
    """
    
    if tdnn[tt]==0:
        tdnn[tt]=np.nan
    if tdwsw[tt]==0:
        tdwsw[tt]=np.nan
    
    if tdsw[tt]==0:
        tdsw[tt]=np.nan
    if tdnw[tt]==0:
        tdnw[tt]=np.nan
    if tdne[tt]==0:
        tdne[tt]=np.nan
    if tdse[tt]==0:
        tdse[tt]=np.nan
    if tdese[tt]==0:
        tdese[tt]=np.nan
    if tdvo[tt]==0:
        tdvo[tt]=np.nan
    if tdvoe[tt]==0:
        tdvoe[tt]=np.nan
    #if tdwc[tt]==0:
    #    tdwc[tt]=np.nan
    #if tdsh[tt]==0:
    #    tdsh[tt]=np.nan
    if tdqv[tt]==0:
        tdqv[tt]=np.nan
    if tdrh[tt]==0:
        tdrh[tt]=np.nan

np.save('tdsws',tdsws)
np.save('tdses',tdses)
np.save('tdeses',tdeses)

"""

np.save('tdnn',tdnn)
np.save('tdwsw',tdwsw)



np.save('tdsw',tdsw)
np.save('tdnw',tdnw)
np.save('tdne',tdne)
np.save('tdse',tdse)
np.save('tdese',tdese)
#np.save('tdwc',tdwc)



np.save('tdqv',tdqv)
np.save('tdrh',tdrh) 
#np.save('tdsh',tdsh)
#np.save('tdvo',tdvo)
#np.save('tdvoe',tdvoe)

"""

#SW histogram
#plt.hist(tdsw, bins=60, density=False, edgecolor='#E6E6E6', color = 'purple', cumulative = False)
#plt.legend()
#plt.xlabel('avg. zonal wind in SW(m/s)')
#plt.savefig('fig/sw.png')
#plt.show()

#DP histogram
#plt.hist(tddp, bins=60, density=False, edgecolor='#E6E6E6', color = 'purple', cumulative = False)
#plt.legend()
#plt.xlabel('diameter of outermost closed isobar(degree)')
#plt.savefig('fig/dp.png')
#plt.show()

#VO histogram
#plt.hist(tdvo, bins=60, density=False, edgecolor='#E6E6E6', color = 'purple', cumulative = False)
#plt.legend()
#plt.xlabel('percentage of vor.>0.000008s-1')
#plt.savefig('fig/vo.png')
#plt.show()

#VOE histogram
#plt.hist(tdvoe, bins=60, density=False, edgecolor='#E6E6E6', color = 'purple', cumulative = False)
#plt.legend()
#plt.xlabel('percentage of vor.>0.000008s-1')
#plt.savefig('fig/voe.png')
#plt.show()

#WC histogram
#plt.hist(tdvo, bins=60, density=False, edgecolor='#E6E6E6', color = 'purple', cumulative = False)
#plt.legend()
#plt.xlabel('num. of grid points where ws.<4ms-1')
#plt.savefig('fig/wc.png')
#plt.show()

#SH histogram
#plt.hist(tdsh, bins=60, density=False, edgecolor='#E6E6E6', color = 'purple', cumulative = False)
#plt.legend()
#plt.xlabel('weighted average meridional shear')
#plt.savefig('fig/sh.png')
#plt.show()

"""
    
tdsw=np.load('tdsw.npy')
tdnw=np.load('tdnw.npy')
tdse=np.load('tdse.npy')
tdese=np.load('tdese.npy')
tddp=np.load('tddp.npy')
tdvo=np.load('tdvo.npy')
tdvoe=np.load('tdvoe.npy')
#tdsh=np.load('tdsh.npy')
tdwc=np.load('tdwc.npy')

"""

for tt in range(tdtimestr.size):
    if tdsw[tt]>-7 and tdsw[tt]<23 and tddp[tt]>=8 and tdwc[tt]>=10: # and tddp[tt]>0 and tddp[tt]<30 and tdwc[tt]>0 and tdwc[tt]<160 and tdvo[tt]>50 and tdvo[tt]<1200
        tdew[tt]=-2

"""

num=0
for tt in range(tdtimestr.size):
    if tdsw[tt]>-7 and tdsw[tt]<23 and tddp[tt]>0 and tddp[tt]<30: # and tdwc[tt]>=0 and tdwc[tt]<200 and tdvo[tt]>50 and tdvo[tt]<1200
        num+=1
        
    #if tddp[tt]>=500:
    #    tdew[tt]=-2
#num=201 #more than 10yrs better
print(num) #368
X = np.zeros((num,4))
NX = np.zeros((num,4))
DX = np.zeros((num,2))
tdnum = np.zeros(num)
uu=0
for tt in range(tdtimestr.size):
    if tdsw[tt]>-7 and tdsw[tt]<23 and tddp[tt]>0 and tddp[tt]<30:
        X[uu,1]=tdsw[tt]#+tdse[tt]
        X[uu,2]=tddp[tt]
        X[uu,0]=tdese[tt]#+tdsw[tt]
        DX[uu,1]=tdvo[tt]-tdvoe[tt]
        X[uu,3]=tdse[tt]
        #X[uu,3]=tdnw[tt]
        #X[uu,4]=tdwc[tt]
        #X[uu,3]=abs(tdnw[tt])+abs(tdsw[tt])
        tdnum[uu]=tt
        uu+=1
        if uu==num:
            break

NX[:,0]=( X[:,0] - np.mean(X[:,0]) ) / ( np.std(X[:,0]) )
NX[:,1]=( X[:,1] - np.mean(X[:,1]) ) / ( np.std(X[:,1]) )
NX[:,2]=( X[:,2] - np.mean(X[:,2]) ) / ( np.std(X[:,2]) )
NX[:,3]=( X[:,3] - np.mean(X[:,3]) ) / ( np.std(X[:,3]) )
#NX[:,4]=( X[:,4] - np.mean(X[:,4]) ) / ( np.std(X[:,4]) )
#NX[:,5]=( X[:,5] - np.mean(X[:,5]) ) / ( np.std(X[:,5]) )
#print(X)

#https://www.kdnuggets.com/2017/03/naive-sharding-centroid-initialization-method.html
#centroids_scaled = naive_sharding(X, 4)
#print(centroids_scaled)
#iiii=np.zeros((4,4))
#iiii[0]=[0,400,500,400]
#iiii[1]=[6,600,500,600]
#iiii[2]=[12,600,500,800]
#iiii[3]=[12,800,1000,800]
#print(iiii)

#https://www.analyticsvidhya.com/blog/2021/05/k-mean-getting-the-optimal-number-of-clusters/
#https://stackoverflow.com/questions/19197715/scikit-learn-k-means-elbow-criterion
silhouette_avg = []
sse = {}
for i in range(3,14):
    kmeans_fit = KMeans(n_clusters = i).fit(NX)
    sse[i] = kmeans_fit.inertia_
    silhouette_avg.append(silhouette_score(NX, kmeans_fit.labels_))

plt.figure()
plt.plot(range(3,14), list(sse.values()))
plt.xlabel("Number of cluster")
plt.ylabel("SSE")
plt.show()

ks,gaps=gap.gap(NX,ks=range(3, 14))
plt.plot(range(3,14), gaps)
plt.xlabel("Number of cluster")
plt.ylabel("Gap")
plt.savefig('fig/n_clusters_gap.png')
plt.show()

plt.plot(range(3,14), silhouette_avg)
plt.xlabel("Number of cluster")
plt.ylabel("Silhouette Score")
plt.savefig('fig/n_clusters_silhouette_score.png')
plt.show()

n_clusters=4

kmeans = KMeans(n_clusters=n_clusters)
kmeans.fit(NX)
Z = kmeans.predict(NX)
centroids = kmeans.cluster_centers_
#score = silhouette_score(NX,kmeans.fit(NX).labels)
#print(score)
#plt.scatter(tdsw.T[0], tdsw.T[1], c=new_dy, cmap=plt.cm.Set1)
colors = ['red','blue','purple','green']

for n in range(n_clusters):
    ys = X[:,0][ Z==n ]
    #if min(ys)!=-999:
    xs = X[:,1][ Z==n ]
    if centroids.T[0][n]>0.9: # ==max(centroids.T[0][:])  and centroids.T[1][n]>0.4
            #print(n,min(ys))
            plt.scatter(xs, ys, s=10, color=colors[3])
            for uu in range(tdnum.size):
                if Z[uu]==n:
                    tdew[int(tdnum[uu])]=-1.5 #mons
            
    elif centroids.T[0][n]<-0.63:
            #print(n,max(ys))
            plt.scatter(xs, ys, s=10, color=colors[0])
            for uu in range(tdnum.size):
                if Z[uu]==n:
                    tdew[int(tdnum[uu])]=1 #EW
                    
    elif centroids.T[0][n]>-0.63 and centroids.T[0][n]<-0.25:
            #print(n,max(ys))
            plt.scatter(xs, ys, s=10, color=colors[2])
            for uu in range(tdnum.size):
                if Z[uu]==n:
                    tdew[int(tdnum[uu])]=-1 #monc

    else:
            plt.scatter(xs, ys, s=10, color=colors[1])
            for uu in range(tdnum.size):
                if Z[uu]==n:
                    tdew[int(tdnum[uu])]=-2 #mond

    #plt.scatter(centroids.T[0][n], centroids.T[0][n], marker='x', color='black')
    print(n,centroids.T[0][n])
    print(n,xs.shape)
    


plt.ylabel('avg. zonal wind in ESE(m/s)')
plt.xlabel('avg. zonal wind in SW(m/s)')
plt.subplots_adjust(left=0.25,right=0.75)
plt.show()



plt.savefig('fig/ews1.png')



#---------mon---------

for n in range(n_clusters):
    ys = X[:,2][ Z==n ]
    #if min(ys)!=-999:
    xs = X[:,1][ Z==n ]
    if centroids.T[0][n]>0.9: #==max(centroids.T[0][:])   and centroids.T[1][n]>0.4
            #print(n,min(ys))
            plt.scatter(xs, ys, s=10, color=colors[3])

    elif centroids.T[0][n]<-0.63:
            #print(n,max(ys))
            plt.scatter(xs, ys, s=10, color=colors[0])
            
    elif centroids.T[0][n]>-0.63 and centroids.T[0][n]<-0.25:
            #print(n,max(ys))
            plt.scatter(xs, ys, s=10, color=colors[2])


    else:
            plt.scatter(xs, ys, s=10, color=colors[1])

    


plt.ylabel('diameter of outermost closed isobar(degree)')
plt.xlabel('avg. zonal wind in SW(m/s)')
plt.subplots_adjust(left=0.25,right=0.75)
plt.show()



plt.savefig('fig/ews2.png')



for n in range(n_clusters):
    ys = X[:,2][ Z==n ]
    #if min(ys)!=-999:
    xs = X[:,0][ Z==n ]
    if centroids.T[0][n]>0.9: #==max(centroids.T[0][:])   and centroids.T[1][n]>0.4
            #print(n,min(ys))
            plt.scatter(xs, ys, s=10, color=colors[3])

    elif centroids.T[0][n]<-0.63:
            #print(n,max(ys))
            plt.scatter(xs, ys, s=10, color=colors[0])
            
    elif centroids.T[0][n]>-0.63 and centroids.T[0][n]<-0.25:
            #print(n,max(ys))
            plt.scatter(xs, ys, s=10, color=colors[2])


    else:
            plt.scatter(xs, ys, s=10, color=colors[1])

    


plt.ylabel('diameter of outermost closed isobar(degree)')
plt.xlabel('avg. zonal wind in ESE(m/s)')
plt.subplots_adjust(left=0.25,right=0.75)
plt.show()




plt.savefig('fig/ews3.png')

np.save('tdew',tdew)
"""

"""

for tt in range(tdtimestr.size):
    if tdew[tt]!=3:
        tdvo[tt]=np.nan
        tdsh[tt]=np.nan


X = np.zeros((480,2))
X[:,0]=tdvo
X[:,1]=tdsh
kmeans = KMeans(n_clusters=2)
kmeans.fit(X)
Z = kmeans.predict(X)
centroids = kmeans.cluster_centers_
#plt.scatter(tdsw.T[0], tdsw.T[1], c=new_dy, cmap=plt.cm.Set1)
colors = ['purple','blue']

if centroids.T[0][0]>centroids.T[0][1]: #VO>5
    plt.scatter(tdsh[ Z==0 ], tdvo[ Z==0 ], s=10, color=colors[1])
    tdew[Z==0]=2 #mong
    plt.scatter(tdsh[ Z==1 ], tdvo[ Z==1 ], s=10, color=colors[0])
    tdew[Z==1]=1 #monsc
else:
    plt.scatter(tdsh[ Z==1 ], tdvo[ Z==1 ], s=10, color=colors[1])
    tdew[Z==1]=2 #mong
    plt.scatter(tdsh[ Z==0 ], tdvo[ Z==0 ], s=10, color=colors[0])
    tdew[Z==0]=1 #monsc
    
plt.scatter(centroids.T[0][n], centroids.T[0][n], marker='x', color='black')
#print(n,centroids.T[0][n])
#print(n,xs.shape)

plt.savefig('fig/ews2.png')
plt.show()
"""




