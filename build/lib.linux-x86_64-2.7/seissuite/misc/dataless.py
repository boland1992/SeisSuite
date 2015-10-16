# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 15:39:42 2015

@author: boland
"""

from obspy.xseed import Parser
import numpy as np
import os

class Dataless:
    """
    Class with representative functions to process input dataless SEED files. 
    """    
    
    def __init__(self, dataless_path=None):
        
        self.metadata_path = dataless_path
        
    def dSEED_XML(self, metadata_path=None):
        """
        Function that imports a given dataless SEED format file from
        its either absolute or relative file path, to a station XML
        file format. This XML location and file name will be x.dataless
        convereted into x.XML and the locations of both files will be the
        same!
        """
        if metadata_path is None: 
            metadata_path = self.metadata_path
                
        dataless_basename =  os.path.basename(metadata_path)
        xml_name = os.path.splitext(dataless_basename)[0]
        xml_path = '{}.xml'.format(xml_name)
        sp = Parser(metadata_path)
        sp.writeXSEED(xml_path) 
            

    def locs_from_dataless(self, metadata_path=None):
        """
        Function that returns a numpy (2,N) shaped array of the longitude
        latitude coordinates (in degree, decimal) from a dataless SEED file. 
        """
        if metadata_path is None: 
            metadata_path = self.metadata_path
            
        sp = Parser(metadata_path)

        metadata = Parser.getInventory(sp)

        lats = np.asarray([float(i['latitude']) for i 
                           in metadata['channels']])
                           
        lons = np.asarray([float(i['longitude']) for i 
                           in metadata['channels']])

        elev = np.asarray([float(i['elevation_in_m']) for i 
                           in metadata['channels']])
        
        return np.column_stack((lons, lats, elev))
        
    def stats_from_dataless(self, metadata_path=None):
        """
        Function that returns a (1,N) shaped array of the station names
        from a dataless SEED file. 
        """
        if metadata_path is None: 
            metadata_path = self.metadata_path
            
        sp = Parser(metadata_path)

        metadata = Parser.getInventory(sp)
        stats = np.asarray([stat['station_id'] for 
                           stat in metadata['stations']])
        return stats
