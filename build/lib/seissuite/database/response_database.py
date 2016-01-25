# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 10:23:55 2015

@author: boland

Goal of this code is to achieve a table of frequency response information 
for each with combined frequency response ranges for each channel in each 
instrument to make a database of which instruments a seismologist 
should use for a  given period, or frequency range. 
"""

import numpy as np
from obspy import read_inventory
import os
import itertools as it
from seissuite.response.resp import (freq_check, window_overlap)
#from obspy.station.response import Response
#import sys
#from xml.etree.ElementTree import parse
#import matplotlib.pyplot as plt
import sqlite3 as lite
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

import itertools 
from obspy.xseed import Parser
from obspy.signal.invsim import pazToFreqResp


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
    
    if len(path_list) > 0:
        return np.asarray(sorted(path_list, key=lambda x: x[1]))[:,0]
    else:
        raise Exception('There is no data to process, please place one or \
more MSEED files in the MSEED_DIR')
    
    
    
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
    
    def __init__(self, metadata_path, SQL_path):
        
        self.SQL_path = SQL_path
        self.metadata_path = metadata_path
        self.extension = os.path.basename(metadata_path).split('.')[-1]
        self.inventory = None
        self.inv_dict = {}
        self.networks = None
        self.stations = None
        self.channels = None
        self.window = None

    def output_SQL(self, code, gains, freqs):
        """
        Function to connect to and output response information to SQL database.
        """
        
        #print 'creating SQL database ... '
        conn = lite.connect(self.SQL_path)
        c = conn.cursor()
        #print 'SQL database cursor initialised'
        
        code = code.replace('.', '_')

        c.execute('CREATE TABLE IF NOT EXISTS {} (gain real, frequency real)'\
        .format(code))
        #c.execute("SELECT name FROM sqlite_master WHERE type='table';")
        #print 'Created frequency response table for {}'.format(code)
        
            
        response_rows = []
        # response lists are too large to store in their entirety. 
        # Reduce 100 fold
        for gain, freq in zip(gains, freqs):
            response_rows.append((gain, freq))

        response_rows = tuple(response_rows)
        c.executemany('INSERT INTO {} VALUES (?,?)'.format(code),response_rows)
        
#        for row in c.execute('SELECT * FROM {}'.format(code)):
#            print row
#        c.execute("SELECT name FROM sqlite_master WHERE type='table';")
#        print c.fetchall()
 
       # commit changes
        conn.commit()
        # close database
        conn.close() 
        
    def output_resp(self, code, gains, freqs):
        """
        Creates a table called resp_windows. This table
        will contain three columns. The first will be the instrument
        code e.g. AU_BENZ_BHZ. The second column will contain the
        minimum effective frequency as calculated by taking the
        response curve, and finding where the line of RESP_EFFECT*gain_max 
        intersect. This will be two points, the min and max of the effective
        frequency window. min_window will be the second column and max
        window will be the third column. This table will be saved to the
        response.db database. 
        """
        
        #print 'creating SQL database ... '
        conn = lite.connect(self.SQL_path)
        c = conn.cursor()
        #print 'SQL database cursor initialised'
        
        code = code.replace('.', '_')
        
        c.execute('''CREATE TABLE IF NOT EXISTS resp_windows (network text, 
                  station text, channel text, min_window real, 
                  max_window real, overlap real)''')

        # check if the net_stat_chan code already exists in the table!
        c.execute('SELECT station FROM resp_windows')
        # create list of codes already in the database
        stations_list = list(it.chain(*list(c.fetchall())))

        # insert information for codes that do not already exist in the table!
        if unicode(code) not in stations_list:
            print code,
        
            # calculate the response window
            window = self.response_window(gains, freqs, tolerance=RESP_EFFECT)

                                     
            overlap = window_overlap(RESP_FREQS, window)
     
            
            net, stat, chan = code.split('_')
            window_tuple = (code, window[0][0], 
                            window[1][0], overlap)
                            
            #print window_tuple
        
            c.execute('INSERT INTO resp_windows VALUES (?,?,?,?,?,?)', 
                      window_tuple)
            #for row in c.execute('SELECT * FROM resp_windows'):
            #    print row
            
        # commit changes
        conn.commit()
        # close database
        conn.close() 

        
        
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

        gain, freq = response.get_evalresp_response(t_samp=t_samp, nfft=nfft)
    
        return np.real(gain)[0::100], freq[0::100]

    def response_window(self, cpx_response, freq, tolerance=RESP_EFFECT):
        """
        Function that can evaluate the response of a given seismic instrument 
        and return a frequency "window" for which the instrument is most 
        effective. The lower the tolerance value (must be float between 0 
        and 1), the larger but less accurate the frequency window will be.
        """

        #print "cpx_resp in window: ", cpx_response
        #print "freqs in window: ", freq        

        #make sure that the gain response array is a numpy array
        cpx_response = np.asarray(cpx_response)
        
        #print "cpx_resp as array: ", cpx_response

        # first find maximum gain response in cpx_reponse
        max_gain = np.max(cpx_response)        
        #print "max gain: ", max_gain        

        gain_tol = max_gain * tolerance
        #print "gain_tol: ", gain_tol        

        arr2 = np.column_stack((freq, abs(cpx_response)))   
        #print "arr2: ", arr2        

        # find indices of cpx_reponse where the grain is above the tolerance
        gain_above = np.argwhere(cpx_response >= gain_tol)
        #print "gain_above: ", gain_above        

        lower_index, upper_index = gain_above[0], gain_above[-1]   
        #print "lower_index: ", lower_index        
        #print "upper_index: ", upper_index        
        
        arr3 = arr2[lower_index:upper_index]
        #print "arr3: ", arr3        
        self.window = np.vstack((arr3[0], arr3[-1]))

        return  self.window
        
    def import_xml(self):
        """
        Import station XML from XML file absolute or relative path. The 
        function then returns two outputs to self.inventory and self.XML
        """
        metadata_path = self.metadata_path
        
        self.inventory = read_inventory(metadata_path, format='STATIONXML')
        
    def import_dataless(self):
        """
        Import station XML from XML file absolute or relative path. The 
        function then returns two outputs to self.inventory and self.XML
        """
        metadata_path = self.metadata_path
        
        return Parser(metadata_path)
    
        
    def xml_resp(self):
        """
        Function used in order to extract both instrument response and 
        instrument sensor make and model information from channel 
        information. 
        """
        
        if self.inventory is None:
            self.import_xml()
            
        min_freq = 1e-4
        inventory = self.inventory
        
        for network in inventory:
            for station in network:
                for channel in station:
                    try:
                        code = '{}_{}_{}'.format(network.code, 
                                                     station.code, 
                                                     channel.code)
                        resp = channel.response
                        #calculate frequency response window
                        sample_rate = self.find_sample(resp)
                        cpx_resp, freqs = self.get_response(min_freq, resp, 
                                                         sample_rate)
                                                         
                        # reduce the number of values saved to database
                        cpx_resp, freqs  = cpx_resp[::100], freqs[::100]  
                        # make sure no gains or frequencies are negative                          
                        cpx_resp, freqs = (np.abs(cpx_resp), np.abs(freqs))
                        self.output_SQL(code, cpx_resp, freqs)
                        self.output_resp(code, np.real(cpx_resp), 
                                         np.real(freqs))

                    except:
                        #print error
                        continue
        
    def dataless_resp(self):
        """
        Function used in order to extract both instrument response and 
        instrument sensor make and model information from channel 
        information. 
        """
        
        sp = self.import_dataless()
        
        # get station information
        min_freq = 1e-4
        
        inventory = sp.getInventory()
        channels = inventory['channels']
        for channel in channels:
            try:
                code = channel["channel_id"]
                net, stat, loc, chan = code.split('.')
                # check for and resolve channel naming problems from IRIS
                if chan == 'MHE':
                    chan = 'LHE'
                if chan == 'MHZ':
                    chan = 'LHZ'                    
                if chan == 'MHN':
                    chan = 'LHN'                    
                    
                sample_rate = channel['sampling_rate']
                data = sp.getPAZ(code)
                poles = data['poles']
                zeros = data['zeros']

             
                t_samp = 1.0 / sample_rate
                #nyquist = sampling_rate / 2.0
                nfft = sample_rate / min_freq
                
                cpx_resp, freqs = pazToFreqResp(poles, zeros, 1,
                                                t_samp, nfft, freq=True)
                                                
                #print "cpx_resp before: ", cpx_resp
                #print "freqs before: ", freqs
                
                # reduce the number of values saved to database; factor of 100
                cpx_resp, freqs  =(np.real(cpx_resp[::100]), 
                                   np.real(freqs[::100]))
                cpx_resp, freqs  =(np.abs(cpx_resp), 
                                   np.abs(freqs))
                #print "cpx_resp after: ", cpx_resp
                #print "freqs after: ", freqs                       
                out_code = '{}_{}_{}'.format(net,stat,chan)
                #print out_code,

                self.output_SQL(out_code, cpx_resp, freqs)
                self.output_resp(out_code, cpx_resp, freqs)
            except:
                a = 5
                
    def dataless_or_xml(self):
        """
        Function to determine whether or not this class should run either
        xml or dataless file functions. 
        """
        if self.extension in ['xml', 'XML', 'Xml']:
            # run xml functions
            self.xml_resp()
        elif self.extension in ['dataless', 'DATALESS', 'Dataless']:
            # run dataless functions    
            self.dataless_resp()
                    
                    
# create database if it doesn't exist already, if it does, stop the programme.
database_name = os.path.join(DATABASE_DIR, 'response.db')

if not AUTOMATE:
    if os.path.exists(database_name):
        yeses = ['y','Y','yes','Yes','YES']    
        nos = ['n','N','no','No','NO']    
    
        condition = False
        while condition is False:
    
            answer = raw_input('Would you like to remove the existing database\
 and start afresh? (y/n): ')
    
            if answer in yeses:
                os.remove(database_name)
                condition = True
            
            elif answer in nos:
                raise Exception("The SQL database {} already exists, \
quitting programme.".format(database_name))
                condition = True

            else:
                print "The input answer must be of yes or no format."
                condition = False
                
else:
    os.remove(database_name)
    
    
xml_paths = paths(folder_path=STATIONXML_DIR, extension='xml')
dataless_paths = paths(folder_path=DATALESS_DIR, extension='dataless')
# combine both xml and dataless paths
abs_paths = list(itertools.chain(*[xml_paths, dataless_paths]))    
# sort by size
abs_paths = paths_sortsize(abs_paths)
for abs_path in abs_paths:
    print '\nScanning file: ', os.path.basename(abs_path)
    INVENTORY = Instrument(abs_path, database_name)
    INVENTORY.dataless_or_xml()