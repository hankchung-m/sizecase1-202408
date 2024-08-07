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
ii=0
if ii==1:
    itimestr=tstimestr
    ihour=tshour
elif ii==2:
    itimestr=tdtimestr
    ihour=tdhour

for tt in range(446,itimestr.size):
    
    print(tt)
    
    #sort real timestr
    if testdt(itimestr[tt]) == False :
        continue
    
    a=ihour[tt]
    if a != '00' and a != '03' and a != '06' and a != '09' and a != '12' and a != '15' and a != '18' and a != '21':
        continue
    
    c = cdsapi.Client()

    def job(date,area):

        year,month,day = date[0:4],date[5:7],date[8:10]

        c.retrieve('reanalysis-era5-pressure-levels',
                   dict(variable=['129','130','131','132','133','135','138','157'],
                        pressure_level=[
                        '1','2','3','5','7','10','20','30','50','70', #0-9
                        '100','125','150','175',
                        '200','225','250', #14-16
                        '300','350','400','450', #17-20
                        '500','550','600','650', #21-24
                        '700','750','775','800','825','850','875', #25-31
                        '900','925','950','975','1000'], #32-36
                        product_type='reanalysis',
                        year=[year],
                        month=[month],
                        day=[day],
                        area = [area],
                        time=[ihour[tt]+':00'], format='netcdf'),'%s%s' %(month,day))
             
    start = itimestr[tt] #-1d?
    end = itimestr[tt] #+1d?
    area = '50/80/-30/-150' #N/W/S/E

    timearray = np.arange(np.datetime64(start), np.datetime64(end)+1)
    filestr = ' '.join(str(e)[5:7]+str(e)[8:10] for e in timearray)
    print(filestr)
    threads = []

    for i in range(len(timearray)):
        threads.append(threading.Thread(target = job, args = (str(timearray[i]),area)))
        threads[i].start()
  
    for s in threads:
        s.join()
        
        if ii==0 :
            if tt/10<1 :
                subprocess.Popen('cat %s > ERA5_I00%d.nc'%(filestr,tt), shell=True)
            elif tt/100<1 :
                subprocess.Popen('cat %s > ERA5_I0%d.nc'%(filestr,tt), shell=True)
            else:
                subprocess.Popen('cat %s > ERA5_I%d.nc'%(filestr,tt), shell=True)
        
        elif ii==1 :
            if tt/10<1 :
                subprocess.Popen('cat %s > ERA5_TS00%d.nc'%(filestr,tt), shell=True)
            elif tt/100<1 :
                subprocess.Popen('cat %s > ERA5_TS0%d.nc'%(filestr,tt), shell=True)
            else:
                subprocess.Popen('cat %s > ERA5_TS%d.nc'%(filestr,tt), shell=True)
        
        elif ii==2 :
            if tt/10<1 :
                subprocess.Popen('cat %s > ERA5_TD00%d.nc'%(filestr,tt), shell=True)
            elif tt/100<1 :
                subprocess.Popen('cat %s > ERA5_TD0%d.nc'%(filestr,tt), shell=True)
            else:
                subprocess.Popen('cat %s > ERA5_TD%d.nc'%(filestr,tt), shell=True)