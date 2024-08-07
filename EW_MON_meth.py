#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  8 22:34:40 2022

@author: cwuhank
"""


import numpy as np
from datetime import date as dt
import netCDF4 as nc
from math import*
import matplotlib.pyplot as plt
#!!!
test=0
should_save=1
should_save_result=0
#!!!
method=1
if method==1:
    from sklearn.cluster import KMeans as cluster
if method==2 or method==0:
    from sklearn.cluster import AgglomerativeClustering as cluster
if method==3:
    from sklearn.mixture import GaussianMixture as cluster
from sklearn.metrics import silhouette_score
import gap
from scipy import stats
from classifiability import kmeanscluster

tslat=np.load('tslat.npy')
tdlon=np.load('tdlon.npy')

tdlat_ny=np.load('tdlat_ny.npy')
tdlon_nx=np.load('tdlon_nx.npy')
tdroci=np.load('tdroci.npy')

tdtimestr=np.load('tdtimestr.npy', allow_pickle=True)
tdhour=np.load('tdhour.npy', allow_pickle=True)
name=np.load('name.npy', allow_pickle=True)
LMIR=np.load('LMIR.npy')

tdew=np.zeros(tdtimestr.size)

tdsws=np.load('tdsws.npy')
tdses=np.load('tdses.npy')
tdeses=np.load('tdeses.npy')
tdnn=np.load('tdnn.npy')
tdwsw=np.load('tdwsw.npy')
tdsw=np.load('tdsw.npy')
tdnw=np.load('tdnw.npy')
tdne=np.load('tdne.npy')
tdse=np.load('tdse.npy')
tdese=np.load('tdese.npy')
tddp=np.load('tddp_0.npy')
tdvo=np.load('tdvo.npy')
tdvoe=np.load('tdvoe.npy')
#tdsh=np.load('tdsh.npy')
tdwc=np.load('tdwc.npy')
tdtb=np.load('tdtb.npy')
tdqv=np.load('tdqv.npy')
tdqv_5_10=np.load('tdqv_5_10.npy')
tdrh=np.load('tdrh.npy')

"""

for tt in range(tdtimestr.size):
    if tdsw[tt]>-7 and tdsw[tt]<23 and tddp[tt]>=8 and tdwc[tt]>=10: # and tddp[tt]>0 and tddp[tt]<30 and tdwc[tt]>0 and tdwc[tt]<160 and tdvo[tt]>50 and tdvo[tt]<1200
        tdew[tt]=-2

"""

fig1,ax1=plt.subplots(1,1,sharex=False,sharey=True,figsize=(8,8))
ax1.grid(True)
fig2,ax2=plt.subplots(3,2,sharex=False,sharey=False,figsize=(20,30))
#ax2.grid(True)

title_wm=['Method-W','Method-MW']
title_num=[['a','b'],['c','d'],['e','f']]
file_wm=['_w','_m']

for wm in range(2):
    num=0
    for tt in range(tdtimestr.size):
        if tdsw[tt]>-7 and tdsw[tt]<23 and tddp[tt]>0 and tddp[tt]<3000 and tdlon[tt]>0:# and tdtb[tt]>0: # and tdwc[tt]>=0 and tdwc[tt]<200 and tdvo[tt]>50 and tdvo[tt]<1200
            num+=1
            
        #if tddp[tt]>=500:
        #    tdew[tt]=-2
    #num=201 #more than 10yrs better
    print(num) #376
    if wm==0:
        var=3
    else:
        var=5#!!!
    X = np.zeros((num,var))
    NX = np.zeros((num,var))
    DX = np.zeros((num,2))
    tdnum = np.zeros(num)
    uu=0
    for tt in range(tdtimestr.size):
        if tdsw[tt]>-7 and tdsw[tt]<23 and tddp[tt]>0 and tddp[tt]<3000 and tdlon[tt]>0:# and tdtb[tt]>0:
    
            X[uu,1]=tdsw[tt]#(tdsw[tt]+tdsws[tt])/2#tdsw[tt]#
            #X[uu,5]=tddp[tt]#tdqv_5_10[tt]
            X[uu,2]=tdse[tt]#(tdse[tt]+tdses[tt])/2#tdse[tt]#
            X[uu,0]=tdese[tt]#(tdese[tt]+tdeses[tt])/2#tdese[tt]#
            if wm==1:
                X[uu,3]=tdqv[tt]
                X[uu,4]=tdnn[tt]
            DX[uu,0]=tdqv[tt]
            DX[uu,1]=tdnn[tt]
            #DX[uu,1]=tdvo[tt]-tdvoe[tt]
            #DX[uu,0]=(tdese[tt]+tdse[tt]+tdsw[tt])/3
            #X[uu,3]=tdtb[tt]
            
            #X[uu,5]=tdnw[tt]
            #X[uu,3]=abs(tdnw[tt])+abs(tdsw[tt])
            tdnum[uu]=tt
            uu+=1
            if uu==num:
                break
    
    NX[:,0]=( X[:,0] - np.mean(X[:,0]) ) / ( np.std(X[:,0]) )
    NX[:,1]=( X[:,1] - np.mean(X[:,1]) ) / ( np.std(X[:,1]) )
    NX[:,2]=( X[:,2] - np.mean(X[:,2]) ) / ( np.std(X[:,2]) )
    if wm==1:
        NX[:,3]=( X[:,3] - np.mean(X[:,3]) ) / ( np.std(X[:,3]) )
        NX[:,4]=( X[:,4] - np.mean(X[:,4]) ) / ( np.std(X[:,4]) )
    
    
    #NX[:,5]=( X[:,5] - np.mean(X[:,5]) ) / ( np.std(X[:,5]) )
    #print('ese',np.mean(X[:,0]),np.std(X[:,0]))
    #print('qv',np.mean(X[:,1]),np.std(X[:,1]))
    #print('dp',np.mean(X[:,2]),np.std(X[:,2]))
    #print('sw',np.mean(X[:,3]),np.std(X[:,3]))
    #print('se',np.mean(X[:,4]),np.std(X[:,4]))
    #print('qv_outer',np.mean(X[:,5]),np.std(X[:,5]))
    #print(X)
    
    
    """
    
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
    P=np.zeros(12)
    for i in range(2,14):
        if method!=3:
            kmeans_fit = cluster(n_clusters = i).fit(NX)
            sse[i] = kmeans_fit.inertia_
            silhouette_avg.append(silhouette_score(NX, kmeans_fit.labels_))
        else:
            gmm = cluster(n_components=i).fit(NX)
            labels=gmm.predict(NX)
            silhouette_avg.append(silhouette_score(NX, labels))
        if method==1:
            kmeans = cluster(n_clusters=i)
            kmeans.fit(NX)
            Z = kmeans.predict(NX)
            for e in range(i):
                num=0
                x = X[:,0][ Z==e ]
                num=x.shape[0]
                if e==0:
                    x0=x
                if e==1:
                    x1=x
                if e==2:
                    x2=x
                if e==3:
                    x3=x
                if e==4:
                    x4=x
                if e==5:
                    x5=x
                if e==6:
                    x6=x
                if e==7:
                    x7=x
                if e==8:
                    x8=x
                if e==9:
                    x9=x
                if e==10:
                    x10=x
                if e==11:
                    x11=x
                if e==12:
                    x12=x
            if i==2:
                ST,P[i-2]=stats.kruskal(x0, x1)
            if i==3:
                ST,P[i-2]=stats.kruskal(x0, x1, x2)
            if i==4:
                ST,P[i-2]=stats.kruskal(x0, x1, x2, x3)
            if i==5:
                ST,P[i-2]=stats.kruskal(x0, x1, x2, x3, x4)
            if i==6:
                ST,P[i-2]=stats.kruskal(x0, x1, x2, x3, x4, x5)
            if i==7:
                ST,P[i-2]=stats.kruskal(x0, x1, x2, x3, x4, x5, x6)
            if i==8:
                ST,P[i-2]=stats.kruskal(x0, x1, x2, x3, x4, x5, x6, x7)
            if i==9:
                ST,P[i-2]=stats.kruskal(x0, x1, x2, x3, x4, x5, x6, x7, x8)
            if i==10:
                ST,P[i-2]=stats.kruskal(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9)
            if i==11:
                ST,P[i-2]=stats.kruskal(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
            if i==12:
                ST,P[i-2]=stats.kruskal(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11)
            if i==13:
                ST,P[i-2]=stats.kruskal(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12)
            
    
    
    plt.figure()
    plt.plot(range(2,14), list(sse.values()))
    plt.xlabel("Number of clusters")
    plt.ylabel("SSE")
    plt.show()
    
    #ks,gaps=gap.gap(NX,ks=range(2, 14))
    #plt.plot(range(2,14), gaps)
    #plt.xlabel("Number of clusters")
    #plt.ylabel("Gap")
    #plt.savefig('fig/n_clusters_gap.png')
    #plt.show()
    
    plt.plot(range(2,14), silhouette_avg)
    plt.xlabel("Number of clusters")
    
    """
    
    if method==1:
        ci=kmeanscluster(NX)
        a=ax1.plot(range(3,14), ci[1:,0],label=title_wm[wm])
        ax1.scatter(np.arange(3, 14, 1), ci[1:,0], marker='o', s=50)

        #plt.plot(range(2,14), P)
        #plt.xlabel("Number of clusters")
        #plt.ylabel("P-value")
        #if test==0:
        #    plt.savefig('fig/n_clusters_P.png')
        #plt.show()
    
    
    n_clusters=4#!!!
    n_iter=100
    
    
    if method==1:
        result=np.zeros((n_iter,n_clusters))
        for i in range(n_iter):
            kmeans = cluster(n_clusters=n_clusters)
            kmeans.fit(NX)
            Z = kmeans.predict(NX)
            centroids = kmeans.cluster_centers_
            for n in range(n_clusters):
                xs = X[:,0][ Z==n ]
                result[i,n]=xs.shape[0]
            result[i,0]=result[i,0]*result[i,1]*result[i,2]*result[i,3]#*result[i,4]*result[i,5] #!!!
        #print(stats.mode(result[:,0]))
        maxresult=stats.mode(result[:,0]).mode.item()
        print(maxresult)
        #if test==0:
        newresult=0
        while newresult!=maxresult:
            kmeans = cluster(n_clusters=n_clusters)
            kmeans.fit(NX)
            Z = kmeans.predict(NX)
            centroids = kmeans.cluster_centers_
            for n in range(n_clusters):
                xs = X[:,0][ Z==n ]
                result[i,n]=xs.shape[0]
            newresult=result[i,0]*result[i,1]*result[i,2]*result[i,3]#*result[i,4]*result[i,5] #!!!
            
    
    if method==2 or method==0:
        ml = cluster(n_clusters=n_clusters)
        Z=ml.fit_predict(NX)
        #Z = ml.predict(NX)
        #centroids = ml.cluster_centers_
        
    if method==3:
        gmm = cluster(n_components=n_clusters).fit(NX)
        Z=gmm.predict(NX)
    
    colors = ['red','blue','purple','green']
    
    if method==0:
        Z=np.zeros(num)
        for uu in range(tdnum.size):
            if X[uu,2]>=1000:
                Z[uu]=-2
                tdew[int(tdnum[uu])]=-2
                plt.scatter(X[uu,0], X[uu,2], s=10, color=colors[1])
            else:
                if X[uu,0]>=-3 and X[uu,0]<=2:
                    Z[uu]=-1
                    tdew[int(tdnum[uu])]=-1
                    plt.scatter(X[uu,0], X[uu,2], s=10, color=colors[2])
                elif X[uu,0]>2:
                    Z[uu]=-1.5
                    tdew[int(tdnum[uu])]=-1.5
                    plt.scatter(X[uu,0], X[uu,2], s=10, color=colors[3])
                else:
                    Z[uu]=1
                    tdew[int(tdnum[uu])]=1
                    plt.scatter(X[uu,0], X[uu,2], s=10, color=colors[0])
        plt.ylabel('diameter of outermost closed isobar(degree)')
        plt.xlabel('avg. zonal wind in ESE(m/s)')
        plt.subplots_adjust(left=0.25,right=0.75)
        plt.show()
    
    if method==1:
        for draw in range(3):
            for n in range(n_clusters):
                if draw==0:
                    ys = X[:,1][ Z==n ]
                    xs = X[:,0][ Z==n ]
                if draw==1:
                    ys = DX[:,0][ Z==n ]
                    xs = X[:,0][ Z==n ]
                if draw==2:
                    ys = DX[:,1][ Z==n ]
                    xs = X[:,2][ Z==n ]
                
                #print('QV',centroids.T[3][n])
                #print('ESE',centroids.T[0][n])
                
                if test==0:
                    """
                    if centroids.T[2][n]>0.5:
                        plt.scatter(xs, ys, s=10, color=colors[1])
                        for uu in range(tdnum.size):
                            if Z[uu]==n:
                                tdew[int(tdnum[uu])]=-2
                    """
                    #else:
                    if centroids.T[0][n]<-0.66:
                        ax2[draw,wm].scatter(xs, ys, s=100, color=colors[0])
                        for uu in range(tdnum.size):
                            if Z[uu]==n:
                                tdew[int(tdnum[uu])]=1
                    elif centroids.T[0][n]>0.33 and centroids.T[0][n]<1:
                        if wm==0:
                            ax2[draw,wm].scatter(xs, ys, s=100, color=colors[3])
                            for uu in range(tdnum.size):
                                if Z[uu]==n:
                                    tdew[int(tdnum[uu])]=-1.5
                        else:
                            ax2[draw,wm].scatter(xs, ys, s=100, color=colors[1])
                            for uu in range(tdnum.size):
                                if Z[uu]==n:
                                    tdew[int(tdnum[uu])]=-2
                    elif centroids.T[0][n]>1:
                        if wm==0:
                            ax2[draw,wm].scatter(xs, ys, s=100, color=colors[1])
                            for uu in range(tdnum.size):
                                if Z[uu]==n:
                                    tdew[int(tdnum[uu])]=-2
                        else:
                            ax2[draw,wm].scatter(xs, ys, s=100, color=colors[3])
                            for uu in range(tdnum.size):
                                if Z[uu]==n:
                                    tdew[int(tdnum[uu])]=-1.5
                    else:
                        ax2[draw,wm].scatter(xs, ys, s=100, color=colors[2])
                        for uu in range(tdnum.size):
                            if Z[uu]==n:
                                tdew[int(tdnum[uu])]=-1
                
                if test==1:
                    if n==0:
                        plt.scatter(xs, ys, s=10, color='red')
                        for uu in range(tdnum.size):
                            if Z[uu]==n:
                                tdew[int(tdnum[uu])]=1
                        
                    if n==1:
                        plt.scatter(xs, ys, s=10, color='purple')
                        for uu in range(tdnum.size):
                            if Z[uu]==n:
                                tdew[int(tdnum[uu])]=-1
                    if n==2:
                        plt.scatter(xs, ys, s=10, color='green')
                        for uu in range(tdnum.size):
                            if Z[uu]==n:
                                tdew[int(tdnum[uu])]=-1.5
                        
                    if n==3:
                        plt.scatter(xs, ys, s=10, color='blue')
                        for uu in range(tdnum.size):
                            if Z[uu]==n:
                                tdew[int(tdnum[uu])]=-2
                    if n==4:
                        plt.scatter(xs, ys, s=10, color='black')
                        for uu in range(tdnum.size):
                            if Z[uu]==n:
                                tdew[int(tdnum[uu])]=2
                        
                    if n==5:
                        plt.scatter(xs, ys, s=10, color='yellow')
                        for uu in range(tdnum.size):
                            if Z[uu]==n:
                                tdew[int(tdnum[uu])]=3
                
                
                
                #print('ese',n,centroids.T[0][n])
                #print('qv',n,centroids.T[1][n])
                #print('dp',n,centroids.T[2][n])
                #print('sw',n,centroids.T[3][n])
                #print('se',n,centroids.T[4][n])
                #print('qv_outer',n,centroids.T[5][n])
                print(n,xs.shape)
            
            if draw==0:
                ax2[draw,wm].set_title(title_wm[wm], fontsize=40)
            ax2[draw,wm].set_title('('+title_num[draw][wm]+')',loc='left', fontsize=40)
            ax2[draw,wm].tick_params(labelsize=25)
            if draw==0:
                ax2[draw,wm].set_ylabel('u_SW (m/s)', fontsize=25)
                ax2[draw,wm].set_xlabel('u_ESE (m/s)', fontsize=25)
            elif draw==1:
                ax2[draw,wm].set_ylabel('q_inner (kg/kg)', fontsize=25)
                ax2[draw,wm].set_xlabel('u_ESE (m/s)', fontsize=25)
            else:
                ax2[draw,wm].set_ylabel('u_NN (m/s)', fontsize=25)
                ax2[draw,wm].set_xlabel('u_SE (m/s)', fontsize=25)

            #ax1.legend(prop={'size': 25})
            #plt.tight_layout()

        #plt.subplots_adjust(left=0.25,right=0.75)

        #plt.show()
        
        
    if method==2:
        
        for uu in range(tdnum.size):
            if Z[uu]==3:
                tdew[int(tdnum[uu])]=-2
            if Z[uu]==1:
                tdew[int(tdnum[uu])]=-1.5
            if Z[uu]==2:
                tdew[int(tdnum[uu])]=-1
            if Z[uu]==0:
                tdew[int(tdnum[uu])]=1
            if Z[uu]==4:
                tdew[int(tdnum[uu])]=2
        
        plt.scatter(X[:,0],X[:,2], s=10, c=Z, cmap='tab10')
        plt.ylabel('diameter of outermost closed isobar(degree)')
        plt.xlabel('avg. zonal wind in ESE(m/s)')
        plt.subplots_adjust(left=0.25,right=0.75)
        if test==0:
            plt.savefig('fig/ews_ese_doci_m2.png', dpi=300)
        plt.show()
        
    if method==3:
        for uu in range(tdnum.size):
            if Z[uu]==3:
                tdew[int(tdnum[uu])]=-2
            if Z[uu]==1:
                tdew[int(tdnum[uu])]=-1.5
            if Z[uu]==2:
                tdew[int(tdnum[uu])]=-1
            if Z[uu]==0:
                tdew[int(tdnum[uu])]=1
        plt.scatter(X[:,0],X[:,2], s=10, c=Z, cmap='tab10')
        plt.ylabel('diameter of outermost closed isobar(degree)')
        plt.xlabel('avg. zonal wind in ESE(m/s)')
        plt.subplots_adjust(left=0.25,right=0.75)
        plt.show()
    
    if test==0 and should_save==1 and should_save_result==1:
        np.save('tdew'+file_wm[wm],tdew)
    if test==1:
        np.save('tdew_test',tdew)
        

ax1.set_xticks(np.arange(3,14,1))
ax1.tick_params(labelsize=15)
ax1.set_xlabel("Number of clusters", fontsize=15)
ax1.set_ylabel("CI", fontsize=15)
ax1.legend(prop={'size': 25})
plt.tight_layout()
#fig1.subplots_adjust(right=0.8)
#if test==0:
fig1.savefig('fig/n_clusters_ci.png', dpi=600)
plt.show()

plt.tight_layout()
fig2.subplots_adjust(top=0.97,bottom=0.03,right=0.99,left=0.1)
if test==0 and should_save==1:
    fig2.savefig('fig/ews.png', dpi=600)