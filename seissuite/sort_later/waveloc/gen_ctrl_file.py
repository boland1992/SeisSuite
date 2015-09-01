# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 08:56:57 2015

@author: boland

The following script is used to generate a control file for 
the programme NonLinLoc. The format is such:

Source description (multiple sources can be specified)

GTSRCE  stat_name loc_type x_srce  y_srce   z_srce   elev

Examples:
GTSRCE  STA   XYZ  	27.25  -67.78  0.0  1.242
GTSRCE  CALF  LATLON   43.753  6.922  0.0  1.242
GTSRCE  JOU  LATLONDM  43 38.00 N  05 39.52 E   0.0   0.300

For each control file, the location types should be consistent. Elevation is
measured in km. 

For Example:
GTSRCE 101         LATLON  49.9520 12.5110 0.0 0.677
GTSRCE 102         LATLON  49.9660 12.5440 0.0 0.595
GTSRCE 103         LATLON  49.9780 12.5860 0.0 0.548
GTSRCE 104         LATLON  49.9910 12.6120 0.0 0.531
GTSRCE 105         LATLON  50.0070 12.6490 0.0 0.733
"""

import sys; sys.path.append('../..')
import os
from obspy import read
from classes.dataless import Dataless
import numpy as np

# set input parameters. 
# this script can take either station XML, dataless SEED metadata or MSEED 
# set one import path to the type of file that you want to import for your
# metadata 

in_path = 'S.BHZ.01.2014.696174.dataless'

use_dataless = True

# set output file location and name
outfolder = os.getcwd()
outfile = 'station_info'
outpath = os.path.join(outfolder, outfile)

if use_dataless:
    obj_dataless = Dataless(in_path)
    coords = obj_dataless.locs_from_dataless()
    stats = obj_dataless.stats_from_dataless()
    
    #combine stations IDs with location data
    info = np.column_stack((stats, coords))

if os.path.exists(outpath):
    # generate a searchable object 
    search = True
    search_list = []
else:
    search = False







# now construct the control file line syntax and output to outpath
with open(outpath, "a+") as f:
    
    if search: 
        # if the path already exists, then make sure to only append new lines!
        for line in f:
            search_list.append(line)
            
    for row in info:
        line_str = 'GTSRCE %s LATLON %0.4f %0.4f 0.0 %0.3f\n'\
                                               %(str(row[0].split('.')[-1]), 
                                               float(row[2]),
                                               float(row[1]),
                                               float(row[3])/1e3)
        if search:                                       
            if not line_str in search_list:
                f.writelines(line_str)
                
                
                
