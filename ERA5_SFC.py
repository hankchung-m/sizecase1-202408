import threading
import cdsapi
import numpy as np
import subprocess
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
tdtimestr=np.load('tdtimestr.npy', allow_pickle=True)
tdhour=np.load('tdhour.npy', allow_pickle=True)

#!!! 0:itime 1:tstime 2:tdtime !!!choose
ii=2
if ii==1:
    itimestr=tstimestr
    ihour=tshour
elif ii==2:
    itimestr=tdtimestr
    ihour=tdhour

for tt in range(itimestr.size):
    
    print(tt)
    
    #sort real timestr
    if testdt(itimestr[tt]) == False :
        continue
    
    a=ihour[tt]
    if a != '00' and a != '03' and a != '06' and a != '09' and a != '12' and a != '15' and a != '18' and a != '21':
        continue
    
    start = itimestr[tt] #-1d?
    end = itimestr[tt] #+1d?
    area = '50/80/-30/-150' #N/W/S/E
    
    timearray = np.arange(np.datetime64(start), np.datetime64(end)+1)
    date = str(timearray[0])
    year,month,day = date[0:4],date[5:7],date[8:10]

    c = cdsapi.Client()
    r = c.retrieve('reanalysis-era5-single-levels',
    {
        'variable':['surface_pressure','mean_sea_level_pressure','land_sea_mask','sea_surface_temperature','t2m','d2m','skt','u10','v10'],
        'product_type':'reanalysis',
        'date': [str(e) for e in timearray],

        'area': [area],
        'time': [ihour[tt]+':00'],
        'format':'netcdf'
    })
    
    filestr = ' '.join(str(e)[5:7]+str(e)[8:10] for e in timearray)
    print(filestr)
    
    if ii==0 :
            if tt/10<1 :
                r.download('ERA5_SFC_I00%d.nc'%(tt))
            elif tt/100<1 :
                r.download('ERA5_SFC_I0%d.nc'%(tt))
            else:
                r.download('ERA5_SFC_I%d.nc'%(tt))
        
    elif ii==1 :
            if tt/10<1 :
                r.download('ERA5_SFC_TS00%d.nc'%(tt))
            elif tt/100<1 :
                r.download('ERA5_SFC_TS0%d.nc'%(tt))
            else:
                r.download('ERA5_SFC_TS%d.nc'%(tt))
        
    elif ii==2 :
            if tt/10<1 :
                r.download('ERA5_SFC_TD00%d.nc'%(tt))
            elif tt/100<1 :
                r.download('ERA5_SFC_TD0%d.nc'%(tt))
            else:
                r.download('ERA5_SFC_TD%d.nc'%(tt))

"""
'variable':['134','151','151189','151214','172','sea_surface_temperature'],
        'product_type':'reanalysis',
        'date': [str(e) for e in timearray],
        'area': [area],
        'time': [ihour[tt]+':00'],
        'format':'netcdf'
        
        'year': [date[0:4]],
        'month':[date[5:7]],
        'day' : [date[8:10]],
"""