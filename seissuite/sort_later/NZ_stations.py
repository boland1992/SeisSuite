# -*- coding: utf-8 -*-
"""
Created on Fri Jul  17 15:38:50 2015

@author: boland
"""

from obspy.fdsn import Client
from obspy import UTCDateTime
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
import pickle
client = Client("GEONET")

starttime = UTCDateTime("2014-01-01")
endtime = UTCDateTime("2015-01-01")

inventory = client.get_stations(network="*", 
                                station="*",
                                loc='*',
                                channel="*Z",
                                starttime=starttime, 
                                endtime=endtime,
                                level="response")



# save all response plots


#inventory[0].plot_response(min_freq=1E-4, 
#                           channel="BHZ",  
#                           location="10",
#                           outfile=None)


#help(inventory[0][0])
# goal: to populate a list of stations with appropriate seismic noise frequency 
# response ranges. 

#for net in inventory: 
#    for sta in net:
#        print sta.code
#        channels = sta.channels
#        for channel in channels:
#            print channel.code 


# list of acceptible channels for ambient noise studies
acceptible_channels = ['BHZ', 'MHZ', 'LHZ', 'VHZ', 'UHZ', 
                       'BNZ', 'MNZ', 'LNZ', 'VNZ', 'UNZ']


#print inventory.get_contents()['channels']

#print inventory.get_contents().keys()

#for inv in inventory: 

#    try:
#        inv.plot_response(min_freq=1E-4, 
#                          channel="BHZ",  
#                          location="10",
#                          outfile=None)        
        
            

def get_latlon(inv, check_channels=False):
    """
    Function to return latitude and longitude coordinates of all stations in
    an obspy inventory class object. 
    """

    lats = []
    lons = []
    labels = []
    for net in inv:
        for sta in net:
            if sta.latitude is None or sta.longitude is None:
                msg = ("Station '%s' does not have latitude/longitude "
                       "information and will not be plotted." % label)
                warnings.warn(msg)
                continue
            
            # perform another loop to check if the channels for the station contain
            # any of the acceptible channels for ambient noise tomography. 
            if check_channels:
                channels = sta.channels
                channel_list = []
                for channel in channels:
                    channel_list.append(channel.code)             

                if any(item in acceptible_channels for item in channel_list):
                    label_ = "   " + ".".join((net.code, sta.code))
                    lats.append(sta.latitude)
                    lons.append(sta.longitude)
                    labels.append(label_)                    

                
            else:
                label_ = "   " + ".".join((net.code, sta.code))
                lats.append(sta.latitude)
                lons.append(sta.longitude)
                labels.append(label_)
            

            
    return np.column_stack((lons, lats))



coords_withcheck = get_latlon(inventory, check_channels=True)
coords_withoutcheck = get_latlon(inventory, check_channels=False)



# set boundaries
bbox = [100, 179, -50, -30]

def remove_coords(coordinates, bbox):
    """
    Function that removes coordinates from outside of a specified bbox. 
    
    coordinates: (2,N) numpy array, or python array or list
    bbox: [xmin, xmax, ymin, ymax]
    """
    xmin, xmax, ymin, ymax = bbox[0], bbox[1], bbox[2], bbox[3] 
    
    
    100, 179, -50, -30

    #convert to python list
    coords = list(coordinates)
    for i, coord in enumerate(coords):
        if coord[0] < xmin or coord[0] > xmax or \
           coord[1] < ymin or coord[1] > ymax:
            
            del coords[i]    

    return np.asarray(coords)




        
coords_withoutcheck = remove_coords(coords_withoutcheck, bbox)


fig1 = plt.figure(figsize=(15,15))
plt.title('Locations of All Available NZ Geonet Seismic Stations')
plt.ylabel('Latitude (Degrees)')
plt.xlabel('Longitude (Degrees)')
plt.scatter(coords_withoutcheck[:,0], coords_withoutcheck[:,1])
fig1.savefig('NZ_Geonet_Stations.svg', format='SVG')
plt.clf


        
coords_withcheck = remove_coords(coords_withcheck, bbox)

fig2 = plt.figure(figsize=(15,15))
plt.title('Locations of All NZ Geonet Seismic Stations \n \
with Seismic Noise Range Capable Instrument Responses')
plt.ylabel('Latitude (Degrees)')
plt.xlabel('Longitude (Degrees)')
plt.scatter(coords_withcheck[:,0], coords_withcheck[:,1])
fig2.savefig('NZ_Filtered_Geonet_Stations.svg', format='SVG')
plt.clf

# mainland New Zealand Geonet 2014 operational station locations
NZ_COORDS = coords_withcheck 

with open(u'NZ_COORDS.pickle', 'wb') as f:
    pickle.dump(NZ_COORDS, f, protocol=2)