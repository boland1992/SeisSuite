# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 09:03:13 2015

@author: boland
"""

from xml.etree.ElementTree import parse
import datetime
import numpy as np
import itertools

xml_path = 'NZ_GEONET.xml'

t0 = datetime.datetime.now()
XML = parse(xml_path).getroot()
t1 = datetime.datetime.now()


print "The xml file took:", t1-t0

#help(XML)


#for network in XML.findall('Network'):
#    print network
instr_list = []

#for instr in XML.iter('{http://www.fdsn.org/xml/station/1}Type'):
#    if instr.text not in instr_list:
#        instr_list.append(instr.text)


# create a nested dictionary with the following structure:
# Level 1: Networks
# Level 2: Stations
# Level 3: Channels, Location


METADATA = {}

for net in XML.findall('{http://www.fdsn.org/xml/station/1}Network'):
    network = net.get('code')
    
    # check to see if the network code already exists inside the METADATA dict
    if network not in METADATA.keys():
        METADATA[network] = {}
    for stat in net.findall('{http://www.fdsn.org/xml/station/1}Station'):
        station = stat.get('code')
        if station not in METADATA[network].keys():
            METADATA[network][station] = {}
        
        
        latitude = float(stat[0].text)
        longitude = float(stat[1].text)
        elevation = float(stat[2].text)
        
        location = {'latitude': latitude, 
                    'longitude': longitude,
                    'elevation': elevation}
        
        METADATA[network][station]['location'] = location
        

        
        #for chan in stat.findall('{http://www.fdsn.org/xml/station/1}Channel'):
        chan_list = []
        sensor_list = []
        # find all channels for a given station - 
        for chan in stat.iter('{http://www.fdsn.org/xml/station/1}Channel'):
            chan_name = chan.get('code')
            
            if chan_name not in chan_list:
                chan_list.append(chan_name)
                
            
        for sensor in stat.iter('{http://www.fdsn.org/xml/station/1}Type'):
            if sensor.text not in sensor_list:
                sensor_list.append(sensor.text)
        
        METADATA[network][station]['channels'] = chan_list
        
        
        
        
        METADATA[network][station]['instruments'] = sensor_list
        
for stat in METADATA['NZ'].keys():
    print METADATA['NZ'][stat]['instruments']      
 