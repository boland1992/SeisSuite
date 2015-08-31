# -*- coding: utf-8 -*-
"""
Created on Thu Jul 24 09:36:19 2015

@author: boland

Goal of this code is to achieve a list of all intrumentation make and models
with combined frequency response ranges for each channel in each instrument
to make a database of which instruments a seismologist should use for a 
given period, or frequency range. 
"""
import numpy as np
import pickle
from xml.etree.ElementTree import parse
import datetime
from obspy import read_inventory
import os
from lxml import etree

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
                


def call_response(inventory, inv_dict):
    """
    Function that allows the user to call the response of a given station
    in a given network for a given channel.
    """
    a=5


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
    
    
    
def paths(folder_path, extension):
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

folder_path =  '/storage/ANT/METADATA/METADATA'
extension = 'xml'

# get files and sort by file size!
abs_paths = paths(folder_path, extension)
# sort files by size
abs_paths = paths_sortsize(abs_paths)


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
                print error
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
                print error
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
        window = np.vstack((arr3[0], arr3[-1]))

        return window

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
            print index_code
            indices = self.inv_dict[index_code]

            #get reponse
            resp = inventory[indices[0]][indices[1]][indices[2]].response
            #calculate frequency response window
            sample_rate = self.find_sample(resp)

            cpx_response, freq = self.get_response(min_freq, resp, sample_rate)
            
            window = list(self.response_window(cpx_response, freq)[:,0])
            print window
                   
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
            print error
            
        self.output_dict = output_dict            
            
    def station_operations(self, station):
        """
        Function used in order to extract both instrument response and instrument
        sensor make and model information from station information. 
        """
        
        network_code = station[1]
        station_code = station[0].get('code')
        combo_code = '{}.{}'.format(network_code, station_code)        
        
        print combo_code
        
        channels = station[0].findall(\
        '{http://www.fdsn.org/xml/station/1}Channel')
        
        channels = [[channel, combo_code] for channel in channels]        

        map(self.channel_operations, channels)
        
    def network_operations(self, network):
        """
        Function used in order to extract both instrument response and instrument
        sensor make and model information from network information. 
        """
        
        print type(network)
        try:
            stations = network.findall(\
            '{http://www.fdsn.org/xml/station/1}Station')

            
        except Exception as error:
                print error
                try:
                    self.inv_dict()
                    self.XML_operations()
                    
                except Exception as error:
                    print error
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
                print error
                try:
                    self.inv_dict()
                    XML = self.XML
                    networks = XML.findall(\
                    '{http://www.fdsn.org/xml/station/1}Network')
                
                except Exception as error:
                    print error
                    raise Exception('Please run Instrument.import_XML and \n\
                                     Instrument.index_inv with a valid input\n\
                                     xml file path')
                


        
        map(self.network_operations, networks)
        
# choose smallest size path xml file to test Instrument class
import xml.etree.cElementTree as ET

total_instr_list = []

for xml_path in abs_paths:
    
    loop_t0 = datetime.datetime.now()
    root = ET.iterparse(xml_path)
    instr_elements = [element[1] for element in root \
                      if 'Sensor' in element[1].tag]
    
    instr_list = np.unique([j.text for i in  instr_elements for j in i])
    
    for instr in instr_list: total_instr_list.append(instr)
        
    total_instr_list = list(np.unique(total_instr_list))
    
    outfile = 'total_instrument_list.pickle'

    with open(outfile, 'wb') as f:
        pickle.dump(total_instr_list, f, protocol=2)
        
    loop_t1 = datetime.datetime.now()
    print 'That loop took: ', loop_t1-loop_t0

quit()


for element in root: 
    if 'Sensor' in element[1].tag:
        for child in element[1]:
            print child.text




# initialise Instrument class
INSTRUMENT = Instrument(xml_path)
#INSTRUMENT.import_XML()
INSTRUMENT.index_inv()
INSTRUMENT.XML_operations()

print INSTRUMENT.output_dict
quit()
# create a nested dictionary with the following structure:
# Level 1: Instrument
# Level 2: Channels 
# Level 3: Response Window (Hz)



#for net in XML.findall('{http://www.fdsn.org/xml/station/1}Network'):



outbasename =  os.path.basename(xml_path)
outname = os.path.splitext(outbasename)[0]
outfile = '{}.pickle'.format(outname)
outpath = os.path.join(os.path.dirname(xml_path), outfile)
with open(outpath, 'wb') as f:
    pickle.dump(INSTRUMENT, f, protocol=2)
                
t_net1  = datetime.datetime.now()
print "The total instrument dictionary loop took: ", t_net1-t_net0

                
