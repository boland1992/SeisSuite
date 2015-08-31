# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 12:28:32 2015

@author: boland
"""
from math import sqrt, atan2, radians,degrees, cos, tan, sin, pi, atan, asin, acos
import matplotlib.pyplot as plt
import numpy as np
import datetime

lat1 = -33.
lon1 = -71.6
lat2 = 31.4
lon2 = 121.8
npts = 1e6





#arc = acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon_diff))

#d = arc*R

    


def waypoint_init(lat1,lon1,lat2,lon2,npts):
    
    R = 6371

    lat1 = radians(lat1)
    lon1 = radians(lon1)
    lat2 = radians(lat2)
    lon2 = radians(lon2)
    
    lon_diff = lon2-lon1

    alpha1 = atan2(sin(lon_diff),(cos(lat1)*tan(lat2)-sin(lat1)*cos(lon_diff)))
    #alpha2 = atan2(sin(lon_diff),(-cos(lat2)*tan(lat1)+sin(lat2)*cos(lon_diff)))
    sigma12 = acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon_diff))
    sigma01 = atan2(tan(lat1), cos(alpha1))
    #sigma02 = sigma01+sigma12
    alpha0 = asin(sin(alpha1)*cos(lat1))
    sigma01 = atan2(tan(lat1), cos(alpha1))
    lon01 = atan2(sin(alpha0)*sin(sigma01), cos(sigma01))
    lon0 = lon1 - lon01
    
    d = sigma12*R
    all_d = np.linspace(0,d,npts)  
    
    def waypoints(dist, sigma01=sigma01, alpha0=alpha0, lon0=lon0):
        sigma = sigma01 + dist/R
        lat = degrees(asin(cos(alpha0)*sin(sigma)))
        lon = degrees(atan2(sin(alpha0)*sin(sigma), cos(sigma))) + degrees(lon0)
        #alpha = atan2(tan(alpha0),cos(sigma))
        
        return [lon, lat]
    
    coords = map(waypoints,all_d)

    return np.vstack(coords)





t0=datetime.datetime.now()
coords = waypoint_init(lat1,lon1,lat2,lon2,npts)
t1=datetime.datetime.now()

print "waypoints", t1-t0
plt.figure()
plt.scatter(coords[:,0], coords[:,1])
plt.show()