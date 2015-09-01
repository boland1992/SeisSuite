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
import os
from scipy.optimize import fsolve
import pylab


# set range of periods that seismic noise gives a resolveable signal:
period_range = [1.,40.]
global freq_range
freq_range = [1./max(period_range), 1./min(period_range)]


global acceptible_channels
acceptible_channels = ['BHZ', 'MHZ', 'LHZ', 'VHZ', 'UHZ']
                       #'BNZ', 'MNZ', 'LNZ', 'VNZ', 'UNZ']

outfolder = '/storage/ANT/NZ Station Responses'
# create list of all possible FDSN clients that work under obspy. 
client_list = (u'BGR', u'ETH', u'GEONET', u'GFZ', u'INGV',
               u'IPGP', u'IRIS', u'KOERI', u'LMU', u'NCEDC', 
               u'NEIP', u'NERIES', u'ODC', u'ORFEUS', u'RESIF',
               u'SCEDC', u'USGS', u'USP')
                
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

for net in inventory: 
    print net
    for sta in net:
        print sta

quit()

# save all response plots
#inventory[0].plot_response(min_freq=1E-4, 
#                           channel="BHZ",  
#                           location="10",
#                           outfile=None)


#help(inventory[0][0])
# goal: to populate a list of stations with appropriate seismic noise frequency 
# response ranges. 

def find_sample(reponse):
    """
    Function that can find the sampling rate for a given station.
    """

    for stage in reponse.response_stages[::-1]:
        if (stage.decimation_input_sample_rate is not None and
                stage.decimation_factor is not None):
            sampling_rate = (stage.decimation_input_sample_rate /
                             stage.decimation_factor)
            break
    else:
        msg = ("Failed to autodetect sampling rate of channel from "
               "response stages. Please manually specify parameter "
               "`sampling_rate`")
        raise Exception(msg)

    return sampling_rate

def get_response(min_freq, response, sampling_rate):

    t_samp = 1.0 / sampling_rate
    nyquist = sampling_rate / 2.0
    nfft = sampling_rate / min_freq

    cpx_response, freq = response.get_evalresp_response(
        t_samp=t_samp, nfft=nfft)
    
    return cpx_response, freq 

def response_window(cpx_response, freq, tolerance=0.7):
    """
    Function that can evaluate the response of a given seismic instrument and
    return a frequency "window" for which the instrument is most effective. 
    The lower the tolerance value (must be float between 0 and 1), the larger but
    less accurate the frequency window will be.
    """
    #make sure that the gain response array is a numpy array
    cpx_response = np.asarray(cpx_response)
    # first find maximum gain response in cpx_reponse
    max_gain = np.max(cpx_response)
    
    gain_tol = max_gain * tolerance

    arr2 = np.column_stack((freq, abs(cpx_response)))
    
    
    # find indices of cpx_reponse where the grain is above the tolerance
    gain_above = np.argwhere(cpx_response >= gain_tol)
   
    lower_index, upper_index = gain_above[0], gain_above[-1]
        
    arr3 = arr2[lower_index:upper_index]
    
    window = np.vstack((arr3[0], arr3[-1]))
    #plt.figure()
    #plt.plot(freq, abs(cpx_response))
    #plt.plot(arr3[:,0], arr3[:,1], c='r')
    #plt.scatter(window[:,0], window[:,1], c='g', s=30)
    #plt.show()
    return window



def freq_check(freq_range, freq_window):
    """
    Function to return True if any of frequencies in the frequency range
    found using the response_window function are contained within the
    freq_range set in the initial variables of this programme. 
    """
    boolean = False    
    
    if any(np.min(freq_range) < freq < np.max(freq_range) \
                                 for freq in freq_window):
        boolean = True

    return boolean

def response_plots(inventory, outfolder, acceptible_channels):

    min_freq = 1e-4
    for net in inventory: 
        for sta in net:
            #print sta.code
            channels = sta.channels
            for channel in channels:
                if str(channel.code) in acceptible_channels:
                    resp = channel.response
                    sample_rate = find_sample(resp)
                    cpx_response, freq = get_response(min_freq, resp, sample_rate)
                    window = response_window(cpx_response, freq)
                    #plt.figure()
                    #plt.loglog(freq, abs(cpx_response))
                    #plt.plot(window[:,0], window[:,1],'r')
                    #plt.scatter(freq_range, [ np.max(cpx_response), 
                    #                     np.max(cpx_response)], c='g', s=35)
                    #plt.show()
                    
                
                    outname = '{}.{}.{}.{}.svg'.format(str(net.code),
                                                       str(sta.code), 
                                                       str(channel.location_code),
                                                       str(channel.code))
                    print outname

                    outfile = os.path.join(outfolder,outname)
                    resp.plot(min_freq, outfile = outfile)
                    freq_window = window[:,0]
                    print freq_check(freq_range, freq_window)

def resp_in_window(inventory, freq_range, acceptible_channels):
    """
    Function to return a list of all station codes that whose
    frequency response window contains any frequencies in the frequency
    range specified for the given study e.g. 0.025-1Hz for current ambient
    noise studies (2015).
    """
    min_freq = 1e-4
    chan_codes = []
    for net in inventory: 
        for sta in net:
            #print sta.code
            channels = sta.channels
            for channel in channels:
                if str(channel.code) in acceptible_channels:
                    resp = channel.response
                    sample_rate = find_sample(resp)
                    cpx_response, freq = get_response(min_freq, 
                                                      resp, 
                                                      sample_rate)
                    window = response_window(cpx_response, freq)
                    #plt.figure()
                    #plt.loglog(freq, abs(cpx_response))
                    #plt.plot(window[:,0], window[:,1],'r')
                    #plt.scatter(freq_range, [ np.max(cpx_response), 
                    #                     np.max(cpx_response)], c='g', s=35)
                    #plt.show()
                    
                    freq_window = window[:,0]
                    check = freq_check(freq_range, freq_window)

                    chan_code = '{}.{}.{}.{}'.format(str(net.code),
                                                     str(sta.code), 
                                                     str(channel.location_code),
                                                     str(channel.code))
                                                     
                    if chan_code not in chan_codes and check:
                        chan_codes.append(chan_code)
                    
    
    return chan_codes


#help(inventory[0][0])

#resp = inventory[0][0][0].response
#print resp            


# list of acceptible channels for ambient noise studies


#print inventory.get_contents()['channels']

#print inventory.get_contents().keys()

#for inv in inventory: 

#    try:
#        inv.plot_response(min_freq=1E-4, 
#                          channel="BHZ",  
#                          location="10",
#                          outfile=None)        
        
            

def get_latlon(inv, check_channels=False, check_codes=False):
    """
    Function to return latitude and longitude coordinates of all stations in
    an obspy inventory class object. 
    """
        
    lats = []
    lons = []
    for net in inv:
        for sta in net:
            label_ = "   " + ".".join((net.code, sta.code))

            if sta.latitude is None or sta.longitude is None:
                msg = ("Station '%s' does not have latitude/longitude "
                       "information and will not be plotted." % label_)
                print msg
                continue
            
            for channel in sta.channels:
            # perform another loop to check if the channels for the station contain
            # any of the acceptible channels for ambient noise tomography. 
                if check_channels:
                    channels = sta.channels
                    channel_list = []
                    for channel in channels:
                        channel_list.append(channel.code)             

                    if any(item in acceptible_channels for item in channel_list):
                        lats.append(sta.latitude)
                        lons.append(sta.longitude)
                    
                elif check_codes: 
                    chan_codes = resp_in_window(inv, 
                                                freq_range, 
                                                acceptible_channels)
                                            
                    chan_code = '{}.{}.{}.{}'.format(str(net.code),
                                                     str(sta.code), 
                                                     str(channel.location_code),
                                                     str(channel.code))    
                
                
                    if chan_code in chan_codes:
                        lats.append(sta.latitude)
                        lons.append(sta.longitude)

                elif check_codes and check_channels:
                    chan_codes = resp_in_window(inv, 
                                                freq_range, 
                                                acceptible_channels)
                    chan_code = '{}.{}.{}.{}'.format(str(net.code),
                                                     str(sta.code), 
                                                     str(channel.location_code),
                                                     str(channel.code))    
                
                    channels = sta.channels
                    channel_list = []
                
                    for channel in channels:
                        channel_list.append(channel.code)  
                    
                    if chan_code in chan_codes and \
                    any(item in acceptible_channels for item in channel_list):
                        lats.append(sta.latitude)
                        lons.append(sta.longitude)
                    
                else:
                    lats.append(sta.latitude)
                    lons.append(sta.longitude)
            

            
    return np.column_stack((lons, lats))


coords_original = get_latlon(inventory, check_channels=False)
coords_checkchannels = get_latlon(inventory, check_channels=True)
coords_checkfreq = get_latlon(inventory, check_codes=True)
coords_combcheck = get_latlon(inventory, check_channels=True, check_codes=True)

# set boundaries
bbox = [100, 179, -50, -30]

def remove_coords(coordinates, bbox):
    """
    Function that removes coordinates from outside of a specified bbox. 
    
    coordinates: (2,N) numpy array, or python array or list
    bbox: [xmin, xmax, ymin, ymax]
    """
    xmin, xmax, ymin, ymax = bbox[0], bbox[1], bbox[2], bbox[3] 
    

    #convert to python list
    coords = list(coordinates)
    for i, coord in enumerate(coords):
        if coord[0] < xmin or coord[0] > xmax or \
           coord[1] < ymin or coord[1] > ymax:
            
            del coords[i]    

    return np.asarray(coords)




        
coords_original = remove_coords(coords_original, bbox)


fig1 = plt.figure(figsize=(15,15))
plt.title('Locations of All Available NZ Geonet Seismic Stations')
plt.ylabel('Latitude (Degrees)')
plt.xlabel('Longitude (Degrees)')
plt.scatter(coords_original[:,0], coords_original[:,1])
fig1.savefig('NZ_Geonet_Stations.svg', format='SVG')
plt.clf


        
coords_checkchannels = remove_coords(coords_checkchannels, bbox)

fig2 = plt.figure(figsize=(15,15))
plt.title('Locations of All NZ Geonet Seismic Stations \n \
with Chosen List of Channel Names')
plt.ylabel('Latitude (Degrees)')
plt.xlabel('Longitude (Degrees)')
plt.scatter(coords_checkchannels[:,0], coords_checkchannels[:,1])
fig2.savefig('NZ_coords_checkchannels_Geonet_Stations.svg', format='SVG')
plt.clf




coords_checkfreq = remove_coords(coords_checkfreq, bbox)

fig3 = plt.figure(figsize=(15,15))
plt.title('Locations of All NZ Geonet Seismic Stations \n \
within Top 30% of Instrument Frequency Response Range')
plt.ylabel('Latitude (Degrees)')
plt.xlabel('Longitude (Degrees)')
plt.scatter(coords_checkfreq[:,0], coords_checkfreq[:,1])
fig3.savefig('NZ_coords_checkfreq_Geonet_Stations.svg', format='SVG')
plt.clf



coords_combcheck = remove_coords(coords_combcheck, bbox)

fig4 = plt.figure(figsize=(15,15))
plt.title('Locations of All NZ Geonet Seismic Stations \n \
with Combined Channel List and Frequency Response Range Checks')
plt.ylabel('Latitude (Degrees)')
plt.xlabel('Longitude (Degrees)')
plt.scatter(coords_combcheck[:,0], coords_combcheck[:,1])
fig4.savefig('NZ_coords_combcheck_Geonet_Stations.svg', format='SVG')
plt.clf



# mainland New Zealand Geonet 2014 operational station locations
NZ_COORDS = coords_combcheck 

with open(u'NZ_COORDS.pickle', 'wb') as f:
    pickle.dump(NZ_COORDS, f, protocol=2)