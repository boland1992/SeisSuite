# -*- coding: utf-8 -*-
"""
Created on Thu May  7 09:19:18 2015

@author: boland
"""

import pyproj
import csv
import io
import numpy as np
stat_names = ["ADE", "DNL", "MRAT", "BRAD", "CLIF", "HODL", "HOLS", "KRAN", 
              "LOCU", "MELU", "MRDN", "NARR", "PEGM", "PEGU", "S88U", 
              "SOMU", "WALM", "OUTU"]

from pylab import *
def info_csv():
    info = io.open('/home/boland/Dropbox/University/UniMelb/AGOS/PROGRAMS/ANT/Versions/07.05.2015/tools/stations_info.csv', "rb") 
    info = csv.reader(info)#import in csv file format
    info.next() #get rid of titles
    info = [row for row in info] #create array of data
    stat_names = [row[0] for row in info]
    stat_names = filter(None, stat_names) # fastest
    lats = [row[1] for row in info]
    lats = filter(None, lats) # fastest
    lons = [row[2] for row in info] 
    lons = filter(None, lons) # fastest


    return stat_names, lats, lons
    

stat_names, lats, lons = info_csv()

print(stat_names, lats, lons)

def max_period(dist):
    c = 2.87#set estimated average seismic shear wave velocity for region in km/s 

    max_per = dist / c    
    return max_per

def distance(stat_names, lats, lons):
    """
    Uses pyproj module to calculate the great circle distance between two lat 
    lon points. imports either individual lat lon co-ordinates or a list 
    of lat lon coordinates. Returns distance in km betwen stations
    """
    g = pyproj.Geod(ellps='WGS84')
    
    data = []

    for stat1, lat1, lon1 in zip(stat_names, lats, lons):
        
        for stat2, lat2, lon2 in zip(stat_names, lats, lons):
            
                (az12, az21, dist) = g.inv(lon1, lat1, lon2, lat2)
                dist = dist/1e3
                max_per = max_period(dist)
                if not dist == 0:               
                    data.append([stat1, stat2, dist, max_per])

    return data
    
data = distance(stat_names, lats, lons)



for i in data: print(i)

per = [info[3] for info in data]
avg = np.mean(per)
stnd = np.std(per)
Q75, Q25 = np.percentile(per, [75 ,25])

IQR = Q75 - Q25
IQR = np.subtract(*np.percentile(per, [75, 25])) #calculate interquartile range

print("\nPeriod range for Gippsland is between %0.2f and %0.2f seconds"%(min(per),max(per)))
print("\naverage of %0.2f seconds" % (avg))
print("\nstandard deviation of %0.2f seconds" %(stnd))
print("\nInterquartile range of %0.2f seconds between %0.2f and %0.2f seconds" %(IQR, Q25, Q75)) 





# notched plot
figure(1)
title("Box Plot of maximum periods resolvable between Gippsland stations using\
 ANT. \nResolved using equation from Bensen et al. (2007) with a representative \
Gippsland shear wave velocity of 2.87km/s")
xlabel(
"Period range for Gippsland stations is between %0.2f and %0.2f seconds\
\nmean of %0.2f seconds\
\nstandard deviation of %0.2f seconds\
\nInterquartile range of %0.2f seconds between %0.2f and %0.2f seconds"
%(min(per),max(per), avg, stnd, IQR, Q25, Q75)
)

ylabel("Wave Period (s)")
boxplot(per)


show()



