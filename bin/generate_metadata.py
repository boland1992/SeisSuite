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
site_description = 'Wairakei local network'
starttime = UTCDateTime('2016-01-01')
endtime = UTCDateTime('2016-01-10')


# import CONFIG class initalised in ./configs/tmp_config.pickle
config_pickle = 'configs/tmp_config.pickle'
f = open(name=config_pickle, mode='rb')
CONFIG = pickle.load(f)
f.close()
    
# import variables from initialised CONFIG class.
DATABASE_DIR = CONFIG.DATABASE_DIR
STATIONXML_DIR = CONFIG.STATIONXML_DIR



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
statmap = {'A324':'W01', 'A320':'W02', 'AA09':'W03', 'A333':'W04', 
           'A971':'W05', 'A356':'W06', 'A316':'W25', '9994':'W08', 
           'B15A':'W26','B2FA':'W10', 'B205':'W11', 'A276':'W19', 
           'A956':'W24', 'B205':'W28', 'ALRZ':'W80', 'HRRZ':'W82', 
           'GRRZ':'W81', 'RITZ':'W97', 'WATZ':'W98', 'WHTZ':'W99',
           'KRVZ':'W83', 'OTVZ':'W84', 'PRRZ':'W85', 'WPRZ':'W86', 
           'HSRZ':'W87', 'MRHZ':'W88', 'ARAZ':'W90', 'HATZ':'W91', 
           'HITZ':'W92', 'KATZ':'W93','KUTZ':'W94','POIZ':'W95', 
           'RATZ':'W96' }

# find unique station names
unique_stations = list(it.chain(*list(c.execute('''SELECT DISTINCT station 
                                                   FROM timeline'''))))


print unique_stations


change_stats = False

statlist = []
for code in unique_stations:  
    
    net, stat, chan = code.split('.')
    
    if change_stats:
        try:
            # use station map dictionary from above to generate new names
            new_stat = statmap[stat]
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
        self.name = "WAIRAKEI"
        self.description = "WAIRAKEI"
        self.town = "WAIRAKEI"
        self.county = "WAIRAKEI"
        self.region = " WAIRAKEI"
        self.country = "NZ"

site = Site()

for stat_info in statlist:
    net, stat, chan = str(stat_info).split('.')
    
    if stat in stat_dict.keys():
        lon, lat, elev = stat_dict[stat]
    
        # create obspy.station.channel opject for each channel
        channels_list = []
        for chan in chan_list:
        
            if 'Z' in chan:
                loc_code = ''
            elif 'N' in chan:
                loc_code = ''
            elif 'E' in chan:
                loc_code = ''
            
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

print network_code

network = Network(network_code, stations=station_list, total_number_of_stations=None, 
                  selected_number_of_stations=None, description=None, 
                  comments=None, start_date=None, end_date=None, 
                  restricted_status=None, alternate_code=None, 
                  historical_code=None, data_availability=None)

net_list = [network]
                  
inv = Inventory(net_list, 'IESE Ltd.')
print inv
inv.write(os.path.join(STATIONXML_DIR, '{}_network.xml').format(network_code), 
          'STATIONXML')