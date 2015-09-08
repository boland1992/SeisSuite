# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 10:23:55 2015

@author: boland

Goal of this code is to achieve a list of all intrumentation make and models
with combined frequency response ranges for each channel in each instrument
to make a database of which instruments a seismologist should use for a 
given period, or frequency range. 
"""

import numpy as np
from obspy import read_inventory
import os
#from obspy.station.response import Response
import sys
from xml.etree.ElementTree import parse
import matplotlib.pyplot as plt
import sqlite3 as lite
import itertools as it

try:
    import cPickle as pickle
except:
    import pickle
    print "Caution, database code may run slow due to cPickle failed import"

# import CONFIG class initalised in ./configs/tmp_config.pickle
config_pickle = 'configs/tmp_config.pickle'
f = open(name=config_pickle, mode='rb')
CONFIG = pickle.load(f)
f.close()
    
# import variables from initialised CONFIG class.
AUTOMATE = CONFIG.AUTOMATE
MSEED_DIR = CONFIG.MSEED_DIR
DATABASE_DIR = CONFIG.DATABASE_DIR
DATALESS_DIR = CONFIG.DATALESS_DIR
STATIONXML_DIR = CONFIG.STATIONXML_DIR

RESP_CHECK = CONFIG.RESP_CHECK
RESP_FREQS = CONFIG.RESP_FREQS
RESP_TOL = CONFIG.RESP_TOL
RESP_EFFECT = CONFIG.RESP_EFFECT


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
    #nyquist = sampling_rate / 2.0
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

    return window

def freq_check(freq_range, freq_window):
    """
    Function to return True if any of frequencies in the frequency range
    found using the response_window function are contained within the
    freq_range set in the initial variables of this programme. 
    """
    condition = False    
    
    if any(np.min(freq_range) < freq < np.max(freq_range) \
                                 for freq in freq_window):
        condition = True

    return condition

def windows_in_inv(inventory):
    min_freq = 1e-4
    for net in inventory: 
        for sta in net:
            #print sta.code
            channels = sta.channels
            for channel in channels:
                resp = channel.response
                sample_rate = find_sample(resp)
                cpx_response, freq = get_response(min_freq, resp, sample_rate)
                window = response_window(cpx_response, freq)                    
                return window[:,0]
                

def paths_sortsize(paths):
    """
    Function used in order to sort file paths by file size. 
    
    paths - this input MUST either be an input array or list of absolute
            paths to the files that are required to be sorted by size. 
    """
    # Initialise a new tuple.
    path_list = []
    # Re-populate list with filename, size tuples
    for i, path in enumerate(paths):
        inst_list = [path, float(os.path.getsize(path))]
        path_list.append(inst_list)

    return np.asarray(sorted(path_list,key=lambda x: x[1]))[:,0]

    
    
def paths(folder_path=None, extension='xml'):
    """
    Function that returns a list of desired absolute paths called abs_paths
    of files that contains a given extension e.g. .txt should be entered as
    folder_path, txt. This function will run recursively through and find
    any and all files within this folder with that extension!
    """

    abs_paths = []
    
    for root, dirs, files in os.walk(folder_path):       
        for f in files:           
            fullpath = os.path.join(root, f)           
            if os.path.splitext(fullpath)[1] == '.{}'.format(extension):               
                abs_paths.append(fullpath)       
    return abs_paths

def window_overlap(freq_range, freq_window):
    """
    Function to return the percentage of overlap between two windows 
    containing the min and max of a frequency range. 
    Example input and output:
    freq_range = [1, 10]
    freq_window = [1, 8]
    the output would be 0.8 for this case
    """
    min_window, max_window = np.min(freq_window), np.max(freq_window)
    min_range, max_range = np.min(freq_range), np.max(freq_range)
    range_diff = abs(max_range - min_range)
    
    #try:
        #case one: both min and max freq_window in freq_range
    if all(min_range <= freq <= max_range for freq in freq_window):
        overlap = abs(max_window - min_window) / range_diff
        return overlap

        #case two: only max freq_window in freq_range
    elif min_range <= max_window <= max_range:
        overlap = abs(max_window - min_range) / range_diff
        return overlap

        #case three: only min freq_window in freq_range
    elif min_range <= min_window <= max_range:
        overlap = abs(max_range - min_window) / range_diff
        return overlap

def process_response(code):
    """
    This function takes a trace in the format of: net_stat_chan and
    feeds that into the response.db SQL database table called 'resp_windows'
    and either sets to remove or keep said trace. The default is to keep
    the trace if no code is found. This is the 'benefit of the doubt' 
    protocol. 
    """
    #stats = trace.stats
    # create code to check with SQL database
    #code = '{}_{}_{}'.format(stats.network, stats.station, 
    #                         stats.channel)
                             
    SQL_path = os.path.join(DATABASE_DIR, 'response.db')

    conn = lite.connect(SQL_path)
    c = conn.cursor()            
    overlap = None
    # test all stations 
    c.execute('SELECT station FROM resp_windows')
    station_list = list(it.chain(*list(c.fetchall())))
    for code in station_list:
       c.execute('SELECT min_window, max_window FROM resp_windows WHERE station=(?)', (code,))
       # ensure that if there are multiple entries for the same station code
       # then only choose the first one as all should be the same.     
       window = c.fetchall()[0]

       if freq_check(RESP_FREQS, window):
           overlap = window_overlap(RESP_FREQS, window),
    
    # commit changes
    conn.commit()
    # close database
    conn.close()     

    return overlap

class Instrument:
    """
    Class created to perform fast operations to extract information from 
    a parsed station XML file. This file format is generally needed to be
    taken from an XML file downloaded via obspy.FDSN. The key modules required
    are:
    
    - from xml.etree.ElementTree import parse
    - from obspy import read_inventory
    
    xml_path is the absolute or relative path to the station XML file
    that is to be processed!
    """
    
    def __init__(self, xml_path):
        
        self.xml_path = xml_path
        self.inventory = None
        self.XML = None
        self.inv_dict = {}
        self.networks = None
        self.stations = None
        self.channels = None
        self.output_dict = {}
        self.window = None

    def import_XML(self):
        """
        Import station XML from XML file absolute or relative path. The 
        function then returns two outputs to self.inventory and self.XML
        """
        
        xml_path = self.xml_path
        
        if self.XML or self.inventory is not None: 
            try:
                 self.XML  = parse(xml_path).getroot()
                 self.inventory = read_inventory(xml_path, format="STATIONXML")
            except Exception as error:
                #print error
                raise Exception('There was a problem with importing \n\
                                 the XML file.')
                
        self.XML = parse(xml_path).getroot()
        self.inventory = read_inventory(xml_path, format="STATIONXML")
        
        #return self.XML, self.inventory


    def index_inv(self):
        """
        Function that returns a dictionary whose keys are the network, station,
        and then channel names in that leveled order, and can be called upon
        to find the index in an inventory object of said station. 
        
        Example for use:
        
        inv_index[network][station][channel] = (1,3,4)
        This means that when calling this channel's reponse from the inventory 
        object you must use chan = inventory[1][3][4].
        """
        
        if self.XML or self.inventory is None:
            try:
                self.import_XML()
            except Exception as error:
                #print error
                raise Exception('Please run Instrument.import_XML with a valid\
 input xml file path')
            

        inventory = self.inventory
        inv_dict = self.inv_dict
        for i, net in enumerate(inventory):
            for j, sta in enumerate(net):
                channels = sta.channels
                for k, channel in enumerate(channels):
                    index_code = '{}.{}.{}'.format(net.code, 
                                                   sta.code, 
                                                   channel.code)
                    inv_dict[index_code] = (i, j, k)
                    print index_code,
        self.inv_dict = inv_dict
        

        
    def find_sample(self, response):
        """
        Function that can find the sampling rate for a given station.
        """

        for stage in response.response_stages[::-1]:
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

    def get_response(self, min_freq, response, sampling_rate):

        t_samp = 1.0 / sampling_rate
        #nyquist = sampling_rate / 2.0
        nfft = sampling_rate / min_freq

        cpx_response, freq = response.get_evalresp_response(
            t_samp=t_samp, nfft=nfft)
    
        return cpx_response, freq 

    def response_window(self, cpx_response, freq, tolerance=0.7):
        """
        Function that can evaluate the response of a given seismic instrument 
        and return a frequency "window" for which the instrument is most 
        effective. The lower the tolerance value (must be float between 0 
        and 1), the larger but less accurate the frequency window will be.
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
        self.window = np.vstack((arr3[0], arr3[-1]))
        print self.window
        return  self.window
        
    def plot_window(self, cpx_response, freq, code, tolerance, window=None):
        """
        Function that plots the frequency response of an individual instrument
        channel, and adds two red dots to signify where the reponse window is. 
        """
        if window is None:  
            window = self.window
        
        # initialise figure
        fig = plt.figure()
        plt.loglog(freq, abs(cpx_response))
        plt.scatter(window[:,0], window[:,1], c='r') 
        plt.title('{} Instrument Frequency Response with {}% tolerance'.format(code, tolerance*100))
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Gain (Counts/m/s)')
        plt.show()
        # clear figure
        fig.clf()
        
    def channel_operations(self, channel):
        """
        Function used in order to extract both instrument response and 
        instrument sensor make and model information from channel 
        information. 
        """
        
        min_freq = 1e-4
        inventory = self.inventory
        output_dict = self.output_dict
        
        try:
            chan = channel[0].get('code')
            index_code = '{}.{}'.format(channel[1], chan)
            #print index_code
            indices = self.inv_dict[index_code]

            #get reponse
            resp = inventory[indices[0]][indices[1]][indices[2]].response
            #calculate frequency response window
            sample_rate = self.find_sample(resp)

            cpx_response, freq = self.get_response(min_freq, resp, sample_rate)
            
            tolerance = 0.7
            window = list(self.response_window(cpx_response, freq, 
                                               tolerance=tolerance)[:,0])
            #print window
            
            show = True
            
            if show:
                self.plot_window(cpx_response, freq, index_code, tolerance)
                
            #find the instrument type for the channel
            for sensor in channel[0].findall\
            ('{http://www.fdsn.org/xml/station/1}Sensor'):
                for descript in sensor.findall\
            ('{http://www.fdsn.org/xml/station/1}Description'):
                    if descript.text not in output_dict.keys():
                        output_dict[descript.text] = {}
                        output_dict[descript.text][chan] = window
                    else:
                        output_dict[descript.text][chan] = window
                          
        except Exception as error:
            a=5
            #print error
            
        self.output_dict = output_dict            
            
    def station_operations(self, station):
        """
        Function used in order to extract both instrument response and instrument
        sensor make and model information from station information. 
        """
        
        network_code = station[1]
        station_code = station[0].get('code')
        combo_code = '{}.{}'.format(network_code, station_code)        
        
        #print combo_code
        
        channels = station[0].findall(\
        '{http://www.fdsn.org/xml/station/1}Channel')
        
        channels = [[channel, combo_code] for channel in channels]        

        map(self.channel_operations, channels)
        
    def network_operations(self, network):
        """
        Function used in order to extract both instrument response and instrument
        sensor make and model information from network information. 
        """
        
        #print type(network)
        try:
            stations = network.findall(\
            '{http://www.fdsn.org/xml/station/1}Station')

            
        except Exception as error:
                #print error
                try:
                    self.inv_dict()
                    self.XML_operations()
                    
                except Exception as error:
                    #print error
                    raise Exception('Please run Instrument.import_XML and \n\
                                     Instrument.index_inv with a valid input\n\
                                     xml file path')

        
        stations = [[station, network.get('code')] for station in stations]        
        
        map(self.station_operations, stations)


    def XML_operations(self):
        """
        Function used in order to extract both instrument response and instrument
        sensor make and model information from network information. 
        """
        
        try:
            XML = self.XML
            networks = XML.findall(\
            '{http://www.fdsn.org/xml/station/1}Network')
            
        except Exception as error:
                #print error
                try:
                    self.inv_dict()
                    XML = self.XML
                    networks = XML.findall(\
                    '{http://www.fdsn.org/xml/station/1}Network')
                
                except Exception as error:
                    #print error
                    raise Exception('Please run Instrument.import_XML and \n\
                                     Instrument.index_inv with a valid input\n\
                                     xml file path')
            
        map(self.network_operations, networks)
        
def fast_instrument(xml_path):
    """
    Function used for rapid mapping of all xml_paths to the Instrument
    class and its functions!
    """        
    INSTRUMENT = Instrument(xml_path)
    #INSTRUMENT.import_XML()
    INSTRUMENT.index_inv()
    INSTRUMENT.XML_operations()
    #print INSTRUMENT.output_dict

if __name__ == "__main__":
    
        
    # set sys.argv
    args = sys.argv
    if len(args) < 2:
        raise Exception('Please input a folder path to scan for station \
XML files!') 
    folder_path = args[1]
    #folder_path =  '/storage/ANT/METADATA/METADATA'
    #extension = 'xml'

    # get xml files!
    abs_paths = paths(folder_path=folder_path)
    abs_paths = paths_sortsize(abs_paths)
    
    for xml_path in abs_paths:
        print '\n\nScanning file: ', xml_path
        print 'Processing channels ... '
        INSTRUMENT = Instrument(xml_path)
        INSTRUMENT.index_inv()
        INSTRUMENT.XML_operations()
    #map(fast_instrument, abs_paths)    


