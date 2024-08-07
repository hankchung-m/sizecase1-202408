#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 11:39:21 2022

@author: cwuhank
"""


import numpy as np
import pandas as pd
from datetime import date as dt
from math import*
import matplotlib.pyplot as plt
    
tdtimestr=np.load('tdtimestr.npy', allow_pickle=True)
tdmonth=np.load('tdmonth.npy')
tdew=np.load('tdew_m.npy')
#print(tdmonth[0])
#tdmonth[tdew!=-2]=np.nan
tdmonth_ew = np.zeros((tdew.size,4))
tdmonth_str=np.empty_like(tdmonth_ew,shape=(tdew.size,4),dtype=object)

#dtm = lambda x: int(x[5:7])
#tdmonth=list(map(dtm, tdtimestr))

# 定義月份的名稱
month_str = [dt(1900, i, 1).strftime('%b') for i in range(1, 13)]

# 定義 tdew 類型和顏色
tdew_types = [-2, -1.5, -1, 1]
colors = ['blue', 'green', 'purple', 'red']

# 創建一個 Pandas DataFrame，以便繪製兩張子圖
plotdata_counts = pd.DataFrame(columns=tdew_types)
plotdata_percentages = pd.DataFrame(columns=tdew_types)

# 計算 counts 和 percentages
for month in range(1, 13):
    counts = [np.sum((tdew == tdew_type) & (tdmonth == month)) for tdew_type in tdew_types]
    total_counts = np.sum(counts)
    percentages = [count * 100 / total_counts for count in counts]
    plotdata_counts.loc[month] = counts
    plotdata_percentages.loc[month] = percentages

fig = plt.figure(figsize=(10, 8))
# 繪製上方的子圖 (數量)
ax1 = fig.add_subplot(211)
ax1 = plotdata_counts.plot(kind='bar', stacked=False, color=colors, ax=ax1)
ax1.set_xticks(range(12))
ax1.set_xticklabels(month_str)
ax1.set_ylabel('Count')
ax1.set_title('(a)',loc='left')
ax1.get_legend().remove()
ax1.legend(labels=['MD', 'MS', 'MC', 'EW'])

#plt.savefig('fig/season_EW.png', dpi=300)
#plt.show()

# 繪製下方的子圖 (百分比)
ax2 = fig.add_subplot(212)
ax2 = plotdata_percentages.plot(kind='bar', stacked=True, color=colors, ax=ax2)
ax2.set_xticks(range(12))
ax2.set_xticklabels(month_str)
#ax2.set_xlabel('Month')
ax2.set_ylabel('Percentage')
ax2.set_title('(b)',loc='left')
ax2.set_ylim([0, 100])
ax2.get_legend().remove()

plt.tight_layout()

plt.savefig('fig/season_EW_relative.png', dpi=600)
plt.show()




"""

bins=np.arange(1,14)
month_str=[dt(1900,i,1).strftime('%b') for i in bins[:-1]]

ew=np.zeros(12)
mc=np.zeros(12)
ms=np.zeros(12)
md=np.zeros(12)

colors=['red','purple','green','blue']
labels=['EW','MC','MS','MD']

bins = np.arange(1, 14)
month_str = [dt(1900, i, 1).strftime('%b') for i in bins[:-1]]
month_str = [dt(1900, i, 1).strftime('%b') for i in range(1, 13)]
month_s = [i for i in range(1, 13)]

tdew_types = [-2, -1.5, -1, 1]
colors = ['blue', 'green', 'purple', 'red']

# 創建一個 Pandas DataFrame，以便繪製兩張子圖
plotdata_counts = pd.DataFrame(columns=tdew_types, index=month_s)
plotdata_percentages = pd.DataFrame(columns=tdew_types, index=month_s)

for month in range(1, 13):
    counts = [np.sum((tdew == tdew_type) & (tdmonth == month)) for tdew_type in tdew_types]
    total_counts = np.sum(counts)
    percentages = [count * 100 / total_counts for count in counts]
    plotdata_counts.loc[month] = counts
    plotdata_percentages.loc[month] = percentages

# 繪製兩個子圖
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(10, 8))

# 繪製上方的子圖 (數量)
plotdata_counts.plot(kind='bar', stacked=False, color=colors, ax=ax1)
ax1.set_xticklabels(month_str)  # 設定 x 軸標籤
ax1.set_xlabel('')
ax1.set_ylabel('Count')
ax1.set_title('TDew Types Count by Month')

# 繪製下方的子圖 (百分比)
plotdata_percentages.plot(kind='bar', stacked=True, color=colors, ax=ax2)
ax2.set_xticklabels(month_str)  # 設定 x 軸標籤
ax2.set_xlabel('Month')
ax2.set_ylabel('Percentage')
ax2.set_ylim([0, 100])

# 刪除圖例
ax1.get_legend().remove()
ax2.get_legend().remove()



plt.tight_layout()

"""



"""

#fig,ax=plt.subplots(2,1,sharex=False,sharey=True,figsize=(10,15))

for tt in range(tdew.size):
    if tdew[tt]==1:
        for m in range(12):
            if tdmonth[tt]==m+1:
                ew[m]+=1
    if tdew[tt]==-1:
        for m in range(12):
            if tdmonth[tt]==m+1:
                mc[m]+=1
    if tdew[tt]==-1.5:
        for m in range(12):
            if tdmonth[tt]==m+1:
                ms[m]+=1
    if tdew[tt]==-2:
        for m in range(12):
            if tdmonth[tt]==m+1:
                md[m]+=1
#tdmonth_ew[tdmonth_ew==0]=np.nan

plotdata = pd.DataFrame({
    "EW":ew,
    "MC":mc,
    "MS":ms,
    "MD":md
    }, index=month_str
)

stacked_data = plotdata.apply(lambda x: x*100/sum(x), axis=1)

for tt in range(tdew.size):
    if tdew[tt]==1:
        tdmonth_ew[tt,0]=tdmonth[tt]
    if tdew[tt]==-1:
        tdmonth_ew[tt,1]=tdmonth[tt]
    if tdew[tt]==-1.5:
        tdmonth_ew[tt,2]=tdmonth[tt]
    if tdew[tt]==-2:
        tdmonth_ew[tt,3]=tdmonth[tt]
tdmonth_ew[tdmonth_ew==0]=np.nan



ax = plt.subplot(212)

stacked_data.plot(kind="bar", stacked=True, color = colors, ax=ax, legend=False)
#plt.savefig('fig/season_EW_relative.png', dpi=300)

#fig, ax = plt.subplots()
ax = plt.subplot(211)
ax.hist(tdmonth_ew, bins=bins, density=False, edgecolor='#E6E6E6', color = colors, label=labels, cumulative = False, align='left')
ax.legend()

ax.set_xticks(bins[:-1])
ax.set_xticklabels([dt(1900,i,1).strftime('%b') for i in bins[:-1]] )

"""


#plt.savefig('fig/season_EW.png', dpi=600)
#plt.show()
