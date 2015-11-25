# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 14:08:48 2015

@author: iese
"""
from obspy.station.channel import Channel
from obspy.station.station import Station
from obspy.station.network import Network
from obspy.station.inventory import Inventory
from obspy.core import UTCDateTime
import pickle 
import os
import sqlite3 as lite
import itertools as it
from pandas import read_csv
import glob


network_code = 'XX'
chan_list = ['DHZ', 'DHN', 'DHE']
site_description = 'UNAM local network'
starttime = UTCDateTime('2015-01-01')
endtime = UTCDateTime('2015-12-01')


# import CONFIG class initalised in ./configs/tmp_config.pickle
config_pickle = 'configs/tmp_config.pickle'
f = open(name=config_pickle, mode='rb')
CONFIG = pickle.load(f)
f.close()
    
# import variables from initialised CONFIG class.
DATABASE_DIR = CONFIG.DATABASE_DIR



# call on database to find station information from raw data headers
database_path = os.path.join(DATABASE_DIR, 'timeline.db')
csv_file = glob.glob(os.path.join(DATABASE_DIR, '*.csv'))[0]


csv_info = read_csv(csv_file)

station_list = list(csv_info['station'])
lats_list = list(csv_info['latitude'])
lons_list = list(csv_info['longitude'])
elev_list = list(csv_info['elevation'])

# constructuct coordinates from station information
stat_dict = {}

for i in range(0, len(station_list)):
    stat_code = station_list[i]    
    lon = lons_list[i]
    lat = lats_list[i]
    elev = elev_list[i]
    stat_dict[stat_code] = (lon, lat, elev)


print 'Initialising SQL database ... '
conn = lite.connect(database_path)
c = conn.cursor()
print 'SQL database cursor initialised.'

# mapping old station names to new station names
UNAM_statmap = {'C0200':'AM01', 'C0230':'AM03', 'C0120':'AM04', 'C01E0':'AM05', 
                'C0180':'AM06', 'C01B0':'AM10', 'C0220':'AM12', 'C00D0':'AM13', 
                'C00E0':'AM15', 'C00A0':'AM16', 'C0240':'AM17', 'C0130':'AM15'}

# find unique station names
unique_stations = list(it.chain(*list(c.execute('''SELECT DISTINCT station 
                                                   FROM timeline'''))))



change_stats = False

statlist = []
for code in unique_stations:  
    
    net, stat, chan = code.split('.')
    
    if change_stats:
        try:
            # use station map dictionary from above to generate new names
            new_stat = UNAM_statmap[stat]
            statlist.append((net, new_stat, chan))
    
        except Exception as error:
            print "station cannot be mapped correctly: {}".format(error)
    
    else:
        
        statlist = unique_stations


#statlist = [str(c).split('.')[1] for c in statlist]

print "statlist: ", statlist


station_list = []


class Site:
    
    def __init__(self):
        self.name = "UNAM"
        self.description = "UNAM"
        self.town = "UNAM"
        self.county = "UNAM"
        self.region = "UNAM"
        self.country = "Mexico"

site = Site()

for stat_info in statlist:
    net, stat, chan = str(stat_info).split('.')
    
    if stat in stat_dict.keys():
        lon, lat, elev = stat_dict[stat]
    
        # create obspy.station.channel opject for each channel
        channels_list = []
        for chan in chan_list:
        
            if 'Z' in chan:
                loc_code = '11'
            elif 'N' in chan:
                loc_code = '12'
            elif 'E' in chan:
                loc_code = '13'
            
            depth = 0.0

            channel = Channel(chan, loc_code, lat, lon,
                          elev, depth)

            channels_list.append(channel)
               
        station = Station(stat, lat, lon, elev, channels=channels_list, 
                      site=site, vault=None, 
                      geology=None, equipments=None,
                      operators=None, 
                      creation_date=starttime, termination_date=endtime,
                      total_number_of_channels=None,
                      selected_number_of_channels=None, description=None,
                      comments=None, start_date=starttime, end_date=starttime,
                      restricted_status=None, alternate_code=None,
                      historical_code=None, data_availability=None)  
    
    station_list.append(station)


network = Network(network_code, stations=station_list, total_number_of_stations=None, 
                  selected_number_of_stations=None, description=None, 
                  comments=None, start_date=None, end_date=None, 
                  restricted_status=None, alternate_code=None, 
                  historical_code=None, data_availability=None)

net_list = [network]
                  
inv = Inventory(net_list, 'IESE Ltd.')
print inv
inv.write('{}_network.xml'.format(network_code), 'STATIONXML')