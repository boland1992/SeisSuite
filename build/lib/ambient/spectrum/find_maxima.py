# -*- coding: utf-8 -*-
"""
Created on Fri July 6 11:04:03 2015

@author: boland
"""

import os

import datetime
import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt

from scipy import signal
from obspy import read
from scipy.signal import argrelextrema
from info_dataless import locs_from_dataless

from matplotlib.colors import LogNorm
import pickle 
import fiona
from shapely import geometry
from shapely.geometry import asPolygon, Polygon
from math import sqrt, radians, cos, sin, asin
from info_dataless import locs_from_dataless
from descartes.patch import PolygonPatch
from matplotlib.colors import LogNorm
from scipy.spatial import ConvexHull
from scipy.cluster.vq import kmeans
from shapely.affinity import scale
from matplotlib.path import Path
#------------------------------------------------------------------------------
# CLASSES
#------------------------------------------------------------------------------

class InShape:
    """
    Class defined in order to define a shapefile boundary AND quickly check
    if a given set of coordinates is contained within it. This class uses
    the shapely module.
    """    
    
    def __init__(self, input_shape, coords=0.):
        #initialise boundary shapefile location string input
        self.boundary = input_shape
        #initialise coords shape input        
        self.dots = coords
        #initialise boundary polygon 
        self.polygon = 0.
        #initialise output coordinates that are contained within the polygon 
        self.output = 0.
    
    def shape_poly(self):
        with fiona.open(self.boundary) as fiona_collection:
            # In this case, we'll assume the shapefile only has one later
            shapefile_record = fiona_collection.next()
            # Use Shapely to create the polygon
            self.polygon = geometry.asShape( shapefile_record['geometry'] )
            return self.polygon
            
    def point_check(self, coord):
        """
        Function that takes a single (2,1) shape input, converts the points
        into a shapely.geometry.Point object and then checks if the coord
        is contained within the shapefile.
        """
        self.polygon = self.shape_poly()
        point = geometry.Point(coord[0], coord[1]) 
        if self.polygon.contains(point):
            return coord
            
    def shape_bounds(self):
        """
        Function that returns the bounding box coordinates xmin,xmax,ymin,ymax
        """
        self.polygon = self.shape_poly()
        return self.polygon.bounds
        
    def shape_buffer(self, shape=None, size=1., res=1):
        """
        Function that returns a new polygon of the larger buffered points. 
        Can import polygon into function if desired. Default is 
        self.shape_poly()
        """
        if shape is None:
            self.polygon = self.shape_poly()
        
        return asPolygon(self.polygon.buffer(size, resolution=res)\
        .exterior)
        
    def extract_poly_coords(self, poly):
        if poly.type == 'Polygon':
            exterior_coords = poly.exterior.coords[:]
        elif poly.type == 'MultiPolygon':
            exterior_coords = []
            for part in poly:
                epc = np.asarray(self.extract_poly_coords(part)) # Recursive call
                exterior_coords.append(epc)
        else:
            raise ValueError('Unhandled geometry type: ' + repr(poly.type))
        
        return np.vstack(exterior_coords)
            
    def external_coords(self, shape=None, buff=None, size=1., res=1):
        """
        Function that returns the external coords of a buffered shapely 
        polygon. Note that shape variable input
        MUST be a shapely Polygon object.
        """
        if shape is not None and buff is not None:
            poly = self.shape_buffer(shape=shape, size=size, res=res)
        elif shape is not None: 
            poly = shape 
        else:
            poly = self.shape_poly()
        
        exterior_coords = self.extract_poly_coords(poly)
        
        return  exterior_coords


#------------------------------------------------------------------------------
# IMPORT PATHS TO MSEED FILES
#------------------------------------------------------------------------------
def spectrum(tr):
    wave = tr.data #this is how to extract a data array from a mseed file
    fs = tr.stats.sampling_rate
    #hour = str(hour).zfill(2) #create correct format for eqstring
    f, Pxx_spec = signal.welch(wave, fs, 'flattop', nperseg=1024, scaling='spectrum')
                    #plt.semilogy(f, np.sqrt(Pxx_spec))
    
    if len(f) >= 256:
        column = np.column_stack((f[:255], np.abs(np.sqrt(Pxx_spec)[:255])))   
        return column
    else:
        return 0.

#    x = np.linspace(0, 10, 1000)
#    f_interp = interp1d(np.sqrt(Pxx_spec),f, kind='cubic')
    #x.reverse()
    #y.reverse()
#    print f_interp(x)
    #f,np.sqrt(Pxx_spec),'o',
#    plt.figure()
#    plt.plot(x,f_interp(x),'-' )
#    plt.show()
    


def paths_sort(path):
    """
    Function defined for customised sorting of the abs_paths list
    and will be used in conjunction with the sorted() built in python
    function in order to produce file paths in chronological order.
    """
    base_name = os.path.basename(path)
    
    stat_name = base_name.split('.')[0]    

    date = base_name.split('.')[1]
    
    try:
        date = datetime.datetime.strptime(date, '%Y-%m-%d')
        
        return date, stat_name
    except Exception as e:
        a=4
        
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

    abs_paths = sorted(abs_paths, key=paths_sort)
       
    return abs_paths
    
    
# import background shapefile location
    
shape_path = "/home/boland/Dropbox/University/UniMelb\
/AGOS/PROGRAMS/ANT/Versions/26.04.2015/shapefiles/aus.shp"

# generate shape object

# Generate InShape class
SHAPE = InShape(shape_path)
# Create shapely polygon from imported shapefile 
UNIQUE_SHAPE = SHAPE.shape_poly()


dataless_path = 'ALL_AUSTRALIA.870093.dataless'
stat_locs = locs_from_dataless(dataless_path)


#folder_path = '/storage/ANT/INPUT/DATA/AU-2014'
folder_path = '/storage/ANT/INPUT/DATA/AU-2014'
extension = 'mseed'

paths_list = paths(folder_path, extension)

t0_total = datetime.datetime.now()
figs_counter = 0


#fig1 = plt.figure(figsize=(15,10))
#ax1 = fig1.add_subplot(111)    
#ax1.set_title("Seismic Waveform Power Density Spectrum\n{}".format('S | 2014'))
#ax1.set_xlabel('Frequency (Hz)')
#ax1.set_ylabel('Power Density Spectrum (V RMS)')
#ax1.set_xlim([0,4])
#ax1.grid(True, axis='both', color='gray')
#ax1.set_autoscaley_on(True)
#ax1.set_yscale('log')

# initialise dictionary to hold all maxima information for a given station
# this will be used to return a new dictionary of the average maxima
# for each station over the course of a year.

maxima_dict0 = {}
maxima_dict1 = {}
a=5
for s in paths_list[:2]:

    try:

        split_path = s.split('/')
        stat_info = split_path[-1][:-6]
        net = stat_info.split('.')[0]
        stat = stat_info.split('.')[1]
        net_stat = '{}.{}'.format(net,stat)
        year = split_path[-2].split('-')[0]
        
        t0 = datetime.datetime.now()
        st = read(s)
        t1 = datetime.datetime.now()
        
        if a == 5: # net == 'S':

            print "time taken to import one month mseed was: ", t1-t0
            # set up loop for all traces within each imported stream.
            t0 = datetime.datetime.now()
            pool = mp.Pool()
            spectra = pool.map(spectrum, st[:])
            pool.close()
            pool.join()    
            t1 = datetime.datetime.now()
            print "time taken to calculate monthly spectra: ", t1-t0
        

            # Caclulate weighted average spectrum for this station for this month
            spectra = np.asarray(spectra)
            search = np.where(spectra==0.)
            spectra = np.delete(spectra, search)    
            spectra = np.average(spectra, axis=0)
    
    
            X, Y = spectra[:,0], spectra[:,1]
            extrema_indices = argrelextrema(Y, np.greater)[0]
            maxima_X = X[extrema_indices]
            maxima_Y = Y[extrema_indices]

            local_extrema = np.column_stack((maxima_X, maxima_Y))
            # sort local maxima            
            local_extrema = local_extrema[local_extrema[:, 1].argsort()]
            local_extrema = local_extrema[::-1]
            
            # retrieve the top two maxima from the PDS plot for use on 
            # noise map. 
            max0, max1 = local_extrema[0], local_extrema[1]
            maxes = [max0,max1]



            if not net_stat in maxima_dict0.keys():
                maxima_dict0[net_stat] = []
                
            if net_stat in maxima_dict0.keys():
                #if not len(maxima_dict[stat]) >= 1:
                maxima_dict0[net_stat].append(max0)
            
            
            if net_stat not in maxima_dict1.keys():
                maxima_dict1[net_stat] = []

            if net_stat in maxima_dict1.keys():     
                maxima_dict1[net_stat].append(max1)

            #smooth_Y = np.convolve(X,Y)
            #smooth_X = np.linspace(np.min(X), np.max(X),len(smooth_Y))
            #plt.plot(smooth_X, smooth_Y, c='b', alpha=0.8)
            #plt.plot(X, Y, c='k', alpha=0.5)
            #plt.scatter(maxima_X, maxima_Y, c='r', s=30)
            #plt.show()
            #plt.clf()
    except:
        a=5
    
#plt.figure()
#stack and find average values for all of the above for each station
#for key in maxima_dict0.keys():
#    stat_locs[key]
    
#    maxima_dict0[key] = np.asarray(maxima_dict0[key])
    
#    plt.scatter(maxima_dict0[key][:,0],maxima_dict0[key][:,1], c='b', s=10)    
    
#    maxima_dict0[key] = np.average(maxima_dict0[key], axis=0)
    
#    plt.scatter(maxima_dict0[key][0],maxima_dict0[key][1], c='r', s=30)    

#    print maxima_dict0[key] 

#for key in maxima_dict1.keys():
#    maxima_dict1[key] = np.asarray(maxima_dict1[key])
    
#    plt.scatter(maxima_dict1[key][:,0],maxima_dict1[key][:,1], c='b', s=10)    
    
#    maxima_dict1[key] = np.average(maxima_dict1[key], axis=0)
    
#    plt.scatter(maxima_dict1[key][0],maxima_dict1[key][1], c='r', s=30)  

#plt.show()


noise_info0 = []


#stack and find average values for all of the above for each station
for key in maxima_dict0.keys():
    
    maxima_dict0[key] = np.asarray(maxima_dict0[key])
    
    maxima_dict0[key] = np.average(maxima_dict0[key], axis=0)
        
    noise_info0.append([stat_locs[key][0],
                       stat_locs[key][1], 
                       maxima_dict0[key][1]])
    

noise_info0 = np.asarray(noise_info0)
# dump  noise_info1

with open('noise_info0.pickle', 'wb') as f:
        pickle.dump(noise_info0, f, protocol=2)

fig = plt.figure(figsize=(15,10), dpi=1000)
plt.title('Average Seismic Noise First Peak Maximum PDS\n S Network | 2014')
plt.xlabel('Longitude (degrees)')
plt.ylabel('Latitude (degrees)')
cm = plt.cm.get_cmap('RdYlBu')
cmin, cmax = np.min(noise_info0[:,2]), np.max(noise_info0[:,2])
sc = plt.scatter(noise_info0[:,0], noise_info0[:,1], c=noise_info0[:,2], 
                 norm=LogNorm(vmin=100, vmax=3e4),  s=35, cmap=cm)                 
col = plt.colorbar(sc)
col.ax.set_ylabel('Maximum Power Density Spectrum (V RMS)')      
fig.savefig('station_pds_maxima/Peak1 PDS Average Maxima 2014.svg', format='SVG')





quit()

noise_info1 = []
#stack and find average values for all of the above for each station
for key in maxima_dict1.keys():
    
    maxima_dict1[key] = np.asarray(maxima_dict1[key])
    
    maxima_dict1[key] = np.average(maxima_dict1[key], axis=0)
        
    noise_info1.append([stat_locs[key][0],
                        stat_locs[key][1], 
                        maxima_dict1[key][1]])
    

noise_info1 = np.asarray(noise_info1)


# dump  noise_info1

with open('noise_info0.pickle', 'wb') as f:
        pickle.dump(noise_info0, f, protocol=2)

fig1 = plt.figure(figsize=(15,10), dpi=1000)
plt.title('Average Seismic Noise Second Peak Maximum PDS\n S Network | 2014')
plt.xlabel('Longitude (degrees)')
plt.ylabel('Latitude (degrees)')
cm = plt.cm.get_cmap('RdYlBu')
cmin, cmax = np.min(noise_info1[:,2]), np.max(noise_info1[:,2])
sc = plt.scatter(noise_info1[:,0], noise_info1[:,1], c=noise_info1[:,2], 
                 norm=LogNorm(vmin=100, vmax=3e4),  s=35, cmap=cm)

col = plt.colorbar(sc)
col.ax.set_ylabel('Maximum Power Density Spectrum (V RMS)')


if shape_path is not None and UNIQUE_SHAPE is not None:
    patch = PolygonPatch(UNIQUE_SHAPE, facecolor='white',\
    edgecolor='k', zorder=1)
    ax = fig.add_subplot(111)
    ax.add_patch(patch)
    
fig1.savefig('station_pds_maxima/Peak2 PDS Average Maxima 2014.svg', format='SVG')




with open('noise_info1.pickle', 'wb') as f:
        pickle.dump(noise_info1, f, protocol=2)
