#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 16:24:16 2024

@author: cwuhank
"""

import matplotlib.pyplot as plt
import numpy as np
#import cartopy.crs as ccrs
from datetime import date as dt

from sklearn.metrics import r2_score
from scipy.optimize import curve_fit
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from matplotlib.patches import Ellipse

def testdt(datestr):
    try:
        dt.fromisoformat(datestr)
    except:
        return False
    else:
        return True

wind_t24 = np.load('wind_t24.npy')
wind_0 = np.load('wind_0.npy')
wind_24 = np.load('wind_24.npy')
wind_48 = np.load('wind_48.npy')
wind_72 = np.load('wind_72.npy')

irmw_t24 = np.load('irmw_t24.npy')
irmw = np.load('irmw.npy')
irmw_24 = np.load('irmw_24.npy')
irmw_48 = np.load('irmw_48.npy')
irmw_72 = np.load('irmw_72.npy')

tdew=np.load('tdew_m.npy')

tstimestr=np.load('tstimestr.npy', allow_pickle=True)
tshour=np.load('tshour.npy', allow_pickle=True)

#hhll=['h','l']
rmwind=['rmw_','wind_']
itime_5=['t24','0','24','48']#,'72']
color_5=['blue','purple','red','orange']#,'yellow']
itime=['72','48','24','0','t24']

# 製作figure  
# 创建子图
fig, ax = plt.subplots(2, 2, figsize=(22, 22))
#ax2 = ax
# 设置座标轴上数字的大小
for row in ax:
    for axis in row:
        axis.tick_params(labelsize=15)


axes = {'ew': ax[0, 0], 'mc': ax[0, 1], 'ms': ax[1, 0], 'md': ax[1, 1]}

titles = ['EW', 'MC', 'MS', 'MD']
title_num = [['a', 'b'], ['c', 'd']]

# 初始化計數器
ew = 0
mc = 0
ms = 0
mg = 0

xy = ['irmw_t24', 'irmw', 'irmw_24', 'irmw_48', 'irmw_72', 'wind_t24', 'wind_0', 'wind_24', 'wind_48', 'wind_72']
gen = ['ew_', 'mc_', 'ms_', 'md_']
time=[-72,-48,-24,0,24]

for i in range(len(gen)):
    for j in range(len(xy)):
        exec(gen[i] + xy[j] + '=np.zeros(tstimestr.size)')

# 初始化圖例標記變量
legend_flag = {'ew': 0, 'mc': 0, 'ms': 0, 'md': 0}


"""

axes2 = {'ew': ax[0, 0].twinx(), 'mc': ax[0, 1].twinx(), 'ms': ax[1, 0].twinx(), 'md': ax[1, 1].twinx()}

for tt in range(tstimestr.size):
    irmw_case=[irmw_72[tt],irmw_48[tt],irmw_24[tt],irmw[tt],irmw_t24[tt]]
    wind_case=[wind_72[tt],wind_48[tt],wind_24[tt],wind_0[tt],wind_t24[tt]]
    
    if testdt(tstimestr[tt]) == False:
        continue

    a = tshour[tt]
    if a not in ['00', '03', '06', '09', '12', '15', '18', '21']:
        continue

    if tdew[tt] == 1:
        for j in range(len(xy)):
            exec(gen[0] + xy[j] + '[' + str(tt) + ']=' + xy[j] + '[' + str(tt) + ']')
        axes['ew'].plot(time,irmw_case, color='blue', alpha=0.5, label='RMW')
        axes2['ew'].plot(time,wind_case, color='orange', alpha=0.5, label='Intensity')
        ew += 1

    elif tdew[tt] == -1:
        for j in range(len(xy)):
            exec(gen[1] + xy[j] + '[' + str(tt) + ']=' + xy[j] + '[' + str(tt) + ']')
        axes['mc'].plot(time,irmw_case, color='blue', alpha=0.5, label='RMW')
        axes2['mc'].plot(time,wind_case, color='orange', alpha=0.5, label='Intensity')
        mc += 1
    elif tdew[tt] == -1.5:
        for j in range(len(xy)):
            exec(gen[2] + xy[j] + '[' + str(tt) + ']=' + xy[j] + '[' + str(tt) + ']')
        axes['ms'].plot(time,irmw_case, color='blue', alpha=0.5, label='RMW')
        axes2['ms'].plot(time,wind_case, color='orange', alpha=0.5, label='Intensity')
        ms += 1
    elif tdew[tt] == -2:
        for j in range(len(xy)):
            exec(gen[3] + xy[j] + '[' + str(tt) + ']=' + xy[j] + '[' + str(tt) + ']')
        axes['md'].plot(time,irmw_case, color='blue', alpha=0.5, label='RMW')
        axes2['md'].plot(time,wind_case, color='orange', alpha=0.5, label='Intensity')
        mg += 1

#From Claude
# 将0值替换为np.nan，计算平均值和标准差
results = {}
for g in gen:
    for var in xy:
        var_name = f"{g}{var}"
        exec(f"{var_name}[{var_name} == 0] = np.nan")
        mean = eval(f"np.nanmean({var_name})")
        std = eval(f"np.nanstd({var_name})")
        results[f"{var_name}_mean"] = mean
        results[f"{var_name}_std"] = std





# 绘制椭圆图
for i, g in enumerate(gen):
    row = i // 2
    col = i % 2
    irmw_mean = np.zeros(5)
    irmw_std = np.zeros(5)
    wind_mean = np.zeros(5)
    wind_std = np.zeros(5)
    for k, t in enumerate(itime):
        irmw_var = f"irmw_{t}" if t != '0' else "irmw"
        wind_var = f"wind_{t}"
        
        irmw_mean[k] = results[f"{g}{irmw_var}_mean"]
        irmw_std[k] = results[f"{g}{irmw_var}_std"]
        wind_mean[k] = results[f"{g}{wind_var}_mean"]
        wind_std[k] = results[f"{g}{wind_var}_std"]
    print(g,t,'rmw',irmw_mean)
    print(g,t,'wind',wind_mean)

        
        # 添加圆心
        #ax[row, col].scatter(irmw_mean, wind_mean, color='black', s=100, zorder=5)
    
    # 计算RMW收缩
    rmw_contraction = results[f"{g}irmw_mean"] - results[f"{g}irmw_24_mean"]
    ax[row, col].text(-72, 275, 'RMW contraction (km)=', fontsize=15)
    ax[row, col].text(-72, 250, f'{rmw_contraction:.4f}', fontsize=15)
    #print(wind_mean[4]-wind_mean[3])

    ax[row, col].plot(time,irmw_mean, marker='o', markersize=20,linewidth=10, color='blue', label='RMW')
    if i==0:
        loc='ew'
    if i==1:
        loc='mc'
    if i==2:
        loc='ms'
    if i==3:
        loc='md'
    axes2[loc].plot(time,wind_mean, marker='o', markersize=20,linewidth=10, color='orange', label='Intensity')


# 设置子图标题和样式
axes2['mc'].set_ylabel('Intensity (kt)', color='orange', fontsize=25)
axes2['md'].set_ylabel('Intensity (kt)', color='orange', fontsize=25)
for i in range(2):
    ax[i, 0].set_ylabel('RMW (km)', color='blue', fontsize=25)
    #ax[i, 1].twinx().set_ylabel('Intensity (kt)', color='orange', fontsize=25)
    for j in range(2):
        if i==0 and j==0:
            loc='ew'
        if i==0 and j==1:
            loc='mc'
        if i==1 and j==0:
            loc='ms'
        if i==1 and j==1:
            loc='md'
        ax[i, j].set_xticks(time)
        ax[i, j].grid(True)
        ax[i, j].set_xlim(left=-78, right=30)
        ax[i, j].set_ylim(bottom=0, top=300)
        axes2[loc].set_ylim(bottom=0, top=160)
        ax[1, j].set_xlabel('Time (h)', fontsize=25)
        ax[i, j].tick_params(labelsize=15)
        ax[i, j].tick_params(axis='y',labelcolor='blue')
        axes2[loc].tick_params(labelsize=15)
        axes2[loc].tick_params(axis='y',labelcolor='orange')
        ax[i, j].set_title('(' + title_num[i][j] + ')', loc='left', fontsize=25)
        ax[i, j].set_title(titles[i*2+j], fontsize=25)
        #ax[i, j].legend(prop={'size': 15})

fig.subplots_adjust(left=0.087, right=0.97, bottom=0.087, top=0.97)
#fig.savefig('fig/RMWevolution_EW.png', dpi=600)

"""


for tt in range(tstimestr.size):
    if testdt(tstimestr[tt]) == False:
        continue

    a = tshour[tt]
    if a not in ['00', '03', '06', '09', '12', '15', '18', '21']:
        continue

    if tdew[tt] == 1:
        for j in range(len(xy)):
            exec(gen[0] + xy[j] + '[' + str(tt) + ']=' + xy[j] + '[' + str(tt) + ']')
        axes['ew'].scatter(irmw_t24[tt], wind_t24[tt], 100, 'blue', alpha=0.3, linewidths=0, edgecolor='black')
        axes['ew'].scatter(irmw[tt], wind_0[tt], 100, 'purple', alpha=0.3, linewidths=0, edgecolor='black')
        axes['ew'].scatter(irmw_24[tt], wind_24[tt], 100, 'red', alpha=0.3, linewidths=0, edgecolor='black')
        axes['ew'].scatter(irmw_48[tt], wind_48[tt], 100, 'orange', alpha=0.3, linewidths=0, edgecolor='black')
        #axes['ew'].scatter(irmw_72[tt], wind_72[tt], 100, 'yellow', alpha=0.3, linewidths=0, edgecolor='black')
        ew += 1
    elif tdew[tt] == -1:
        for j in range(len(xy)):
            exec(gen[1] + xy[j] + '[' + str(tt) + ']=' + xy[j] + '[' + str(tt) + ']')
        if legend_flag['mc'] == 1:
            axes['mc'].scatter(irmw_t24[tt], wind_t24[tt], 100, 'blue', alpha=0.3, linewidths=0, edgecolor='black')
            axes['mc'].scatter(irmw[tt], wind_0[tt], 100, 'purple', alpha=0.3, linewidths=0, edgecolor='black')
            axes['mc'].scatter(irmw_24[tt], wind_24[tt], 100, 'red', alpha=0.3, linewidths=0, edgecolor='black')
            axes['mc'].scatter(irmw_48[tt], wind_48[tt], 100, 'orange', alpha=0.3, linewidths=0, edgecolor='black')
            #axes['mc'].scatter(irmw_72[tt], wind_72[tt], 100, 'yellow', alpha=0.3, linewidths=0, edgecolor='black')
        else:
            axes['mc'].scatter(irmw_t24[tt], wind_t24[tt], 100, 'blue', alpha=0.3, linewidths=0, edgecolor='black', label='+24h')
            axes['mc'].scatter(irmw[tt], wind_0[tt], 100, 'purple', alpha=0.3, linewidths=0, edgecolor='black', label='-0h')
            axes['mc'].scatter(irmw_24[tt], wind_24[tt], 100, 'red', alpha=0.3, linewidths=0, edgecolor='black', label='-24h')
            axes['mc'].scatter(irmw_48[tt], wind_48[tt], 100, 'orange', alpha=0.3, linewidths=0, edgecolor='black', label='-48h')
            #axes['mc'].scatter(irmw_72[tt], wind_72[tt], 100, 'yellow', alpha=0.3, linewidths=0, edgecolor='black', label='-72h')
            legend_flag['mc'] = 1
        mc += 1
    elif tdew[tt] == -1.5:
        for j in range(len(xy)):
            exec(gen[2] + xy[j] + '[' + str(tt) + ']=' + xy[j] + '[' + str(tt) + ']')
        axes['ms'].scatter(irmw_t24[tt], wind_t24[tt], 100, 'blue', alpha=0.3, linewidths=0, edgecolor='black')
        axes['ms'].scatter(irmw[tt], wind_0[tt], 100, 'purple', alpha=0.3, linewidths=0, edgecolor='black')
        axes['ms'].scatter(irmw_24[tt], wind_24[tt], 100, 'red', alpha=0.3, linewidths=0, edgecolor='black')
        axes['ms'].scatter(irmw_48[tt], wind_48[tt], 100, 'orange', alpha=0.3, linewidths=0, edgecolor='black')
        #axes['ms'].scatter(irmw_72[tt], wind_72[tt], 100, 'yellow', alpha=0.3, linewidths=0, edgecolor='black')
        ms += 1
    elif tdew[tt] == -2:
        for j in range(len(xy)):
            exec(gen[3] + xy[j] + '[' + str(tt) + ']=' + xy[j] + '[' + str(tt) + ']')
        axes['md'].scatter(irmw_t24[tt], wind_t24[tt], 100, 'blue', alpha=0.3, linewidths=0, edgecolor='black')
        axes['md'].scatter(irmw[tt], wind_0[tt], 100, 'purple', alpha=0.3, linewidths=0, edgecolor='black')
        axes['md'].scatter(irmw_24[tt], wind_24[tt], 100, 'red', alpha=0.3, linewidths=0, edgecolor='black')
        axes['md'].scatter(irmw_48[tt], wind_48[tt], 100, 'orange', alpha=0.3, linewidths=0, edgecolor='black')
        #axes['md'].scatter(irmw_72[tt], wind_72[tt], 100, 'yellow', alpha=0.3, linewidths=0, edgecolor='black')
        mg += 1

#From Claude
# 将0值替换为np.nan，计算平均值和标准差
results = {}
for g in gen:
    for var in xy:
        var_name = f"{g}{var}"
        exec(f"{var_name}[{var_name} == 0] = np.nan")
        mean = eval(f"np.nanmean({var_name})")
        std = eval(f"np.nanstd({var_name})")
        results[f"{var_name}_mean"] = mean
        results[f"{var_name}_std"] = std

# 绘制椭圆图
for i, g in enumerate(gen):
    row = i // 2
    col = i % 2
    for k, t in enumerate(itime_5):
        irmw_var = f"irmw_{t}" if t != '0' else "irmw"
        wind_var = f"wind_{t}"
        
        irmw_mean = results[f"{g}{irmw_var}_mean"]
        irmw_std = results[f"{g}{irmw_var}_std"]
        wind_mean = results[f"{g}{wind_var}_mean"]
        wind_std = results[f"{g}{wind_var}_std"]
        print(g,t,'rmw',irmw_mean)
        print(g,t,'wind',wind_mean)
        
        ellipse = Ellipse(xy=(irmw_mean, wind_mean), 
                          width=irmw_std*2, 
                          height=wind_std*2,
                          facecolor=color_5[k],
                          alpha=0.5,
                          linewidth=5,
                          edgecolor=color_5[k])
        ax[row, col].add_patch(ellipse)
        
        # 添加圆心
        ax[row, col].scatter(irmw_mean, wind_mean, color='black', s=100, zorder=5)
    
    # 计算RMW收缩
    rmw_contraction = results[f"{g}irmw_mean"] - results[f"{g}irmw_24_mean"]
    ax[row, col].text(100, 150, 'RMW contraction (km)=', fontsize=15)
    ax[row, col].text(100, 140, f'{rmw_contraction:.4f}', fontsize=15)



# 设置子图标题和样式
for i in range(2):
    ax[i, 0].set_ylabel('Intensity (kt)', fontsize=25)
    for j in range(2):
        ax[i, j].grid(True)
        ax[i, j].set_xlim(left=0, right=300)
        ax[i, j].set_ylim(bottom=0, top=160)
        ax[1, j].set_xlabel('RMW (km)', fontsize=25)
        ax[i, j].tick_params(labelsize=15)
        ax[i, j].set_title('(' + title_num[i][j] + ')', loc='left', fontsize=25)
        ax[i, j].set_title(titles[i*2+j], fontsize=25)
        ax[i, j].legend(prop={'size': 15})

fig.subplots_adjust(left=0.087, right=0.97, bottom=0.087, top=0.97)
fig.savefig('fig/RMWcontraction_EW.png', dpi=600)
