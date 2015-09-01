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
from scipy import interpolate
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

from scipy.interpolate import griddata
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

# set plotting limits for shapefile boundaries

lonmin, latmin, lonmax, latmax = SHAPE.shape_bounds()

print lonmin, latmin, lonmax, latmax
#lonmin, lonmax, latmin, latmax = SHAPE.plot_lims()


dataless_path = 'ALL_AUSTRALIA.870093.dataless'
stat_locs = locs_from_dataless(dataless_path)


#folder_path = '/storage/ANT/INPUT/DATA/AU-2014'
folder_path = '/storage/ANT/INPUT/DATA/AU-2014'
extension = 'mseed'

paths_list = paths(folder_path, extension)

t0_total = datetime.datetime.now()
figs_counter = 0


pickle_file =  '/storage/ANT/spectral_density/station_pds_maxima/\
S Network 2014/noise_info0_SNetwork2014.pickle'

f = open(name=pickle_file, mode='rb')
noise_info0 = pickle.load(f)
f.close()

# dump  noise_info1

fig = plt.figure(figsize=(15,10), dpi=1000)
plt.title('Average Seismic Noise First Peak Maximum PDS\n S Network | 2014')
plt.xlabel('Longitude (degrees)')
plt.ylabel('Latitude (degrees)')


patch = PolygonPatch(UNIQUE_SHAPE, facecolor='white',\
                     edgecolor='k', zorder=1)
ax = fig.add_subplot(111)
ax.add_patch(patch)


#create 5000 Random points distributed within the circle radius 100
x, y = noise_info0[:,0], noise_info0[:,1]
points = np.column_stack((x,y))

xmin, xmax = np.min(x), np.max(x)
ymin, ymax = np.min(y), np.max(y)


values = noise_info0[:,2]

#now we create a grid of values, interpolated from our random sample above
y = np.linspace(ymin, ymax, 100)
x = np.linspace(xmin, xmax, 100)
gridx, gridy = np.meshgrid(x, y)
heat_field = griddata(points, values, (gridx, gridy), method='cubic',fill_value=0)

print heat_field


heat_field = np.where(heat_field < 0, 1, heat_field)
heat_field = np.ma.masked_where(heat_field==0,heat_field)  


plt.pcolor(gridx, gridy, heat_field,
           cmap='rainbow',alpha=0.5, norm=LogNorm(vmin=100, vmax=3e4),  
           zorder=2)
           
plt.scatter(noise_info0[:,0], noise_info0[:,1], c=noise_info0[:,2], 
            norm=LogNorm(vmin=100, vmax=3e4),  s=35, cmap='rainbow', zorder=3)           

#cmin, cmax = np.min(noise_info0[:,2]), np.max(noise_info0[:,2])
#sc = plt.scatter(noise_info0[:,0], noise_info0[:,1], c=noise_info0[:,2], 
#                 norm=LogNorm(vmin=100, vmax=3e4),  s=50, cmap=cm, zorder=2)                 

col = plt.colorbar()
col.ax.set_ylabel('Maximum Power Density Spectrum (V RMS)')      



ax.set_xlim(lonmin-0.05*abs(lonmax-lonmin), \
            lonmax+0.05*abs(lonmax-lonmin))
ax.set_ylim(latmin-0.05*abs(latmax-latmin), \
            latmax+0.05*abs(latmax-latmin))
            
            
fig.savefig('station_pds_maxima/check1.svg', format='SVG')
