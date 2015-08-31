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
import shapefile
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
import itertools
from scipy.interpolate import griddata
import random
from sklearn.cluster import DBSCAN

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

    
class InPoly:
    """
    Class defined in order to define a shapefile boundary AND quickly check
    if a given set of coordinates is contained within it. The class uses
    the matplotlib Path class. 
    """    
    
    def __init__(self, input_shape, coords=0.):
        #initialise boundary shapefile location string input
        self.boundary = input_shape
        #initialise coords shape input        
        self.dots = coords
        #initialise boundary polygon 
        self.polygon = 0.
        #initialise boundary polygon nodes
        self.nodes = 0.
        #initialise output coordinates that are contained within the polygon 
        self.output = 0.
        
    def poly_nodes(self):
        """
        Function that returns the nodes of a shapefile as a (n,2) array.
        """
        sf = shapefile.Reader(self.boundary)
        poly = sf.shapes()[0]
        #find polygon nodes lat lons
        self.nodes = np.asarray(poly.points)            
        return self.nodes
        
    def points_from_path(self, poly):
        """
        Function that returns nodes from matplotlib Path object. 
        """
        return poly.vertices
    
    def shapefile_poly(self):
        """
        Function that imports a shapefile location path and returns
        a matplotlib Path object representing this shape.
        """        
        self.nodes = self.poly_nodes()
        #convert to a matplotlib path class!
        self.polygon = Path(self.nodes)
        return self.polygon
        
    def node_poly(self, nodes):
        """
        Function creates a matplotlib Path object from input nodes.
        """
        #convert to a matplotlib path class!
        polygon = Path(nodes)
        return polygon
        
    def points_in_shapefile_poly(self):
        """
        Function that takes a single (2,1) coordinate input, and uses the 
        contains() function in class matplotlib Path to check if point is
        in the polygon. 
        """
        self.polygon = self.shapefile_poly()
        points_in = self.polygon.contains_points(self.dots)
        self.output = self.dots[points_in == True]
        return np.asarray(self.output)
        
    def points_in(self, points, poly=None, IN=True, indices=False):
        """
        Function that takes a many (2,N) points, and uses the 
        contains() function in class matplotlib Path to check if point is
        in the polygon. If IN=True then the function will return points inside
        the matplotlib Path object, else if IN=False then the function will
        return the points outside the matplotlib Path object. 
        """
        
        if poly is None: 
            poly = self.shapefile_poly()
        
        points_test = poly.contains_points(points)
        
        if indices:
            return points_test
        else:
            output = points[points_test == IN]
            return np.asarray(output)
        
    def bounds_poly(self, nodes=None):
        """
        Function that returns boundaries of a shapefile polygon.
        """
        if nodes is None: 
            nodes = self.poly_nodes()
            
        xmin, xmax = np.min(nodes[:,0]), np.max(nodes[:,0])
        ymin, ymax = np.min(nodes[:,1]), np.max(nodes[:,1])
        return xmin, xmax, ymin, ymax
        
    def poly_from_shape(self, shape=None, size=1., res=1):
        """
        Function that returns a matplotlib Path object from 
        buffered shape points. if shape != None then the shape input
        MUST be of type shapely polygon. 
        """
        SHAPE = InShape(self.boundary)
        if shape is None:
            # Generates shape object from shape_file input
            shape = SHAPE
            return self.node_poly(shape.external_coords(size=size, res=res))
        else:
            return self.node_poly(SHAPE.external_coords(shape=shape))
    
    def rand_poly(self, poly=None, N=1e4, IN=True):
        """
        Function that takes an input matplotlib Path object (or the default)
        and generates N random points within the bounding box around it. 
        Then M unknown points are returned that ARE contained within the
        Path object. This is done for speed. If IN=True then the function
        will return points inside the matplotlib Path object, else if 
        IN=False then the function will return the points outside the 
        matplotlib Path object. 
        """
        if poly is None: 
            #poly = self.shapefile_poly()
            xmin, xmax, ymin, ymax = self.bounds_poly()
        else:   
            nodes = self.points_from_path(poly)
            xmin, xmax, ymin, ymax = self.bounds_poly(nodes=nodes)
            
        X = abs(xmax - xmin) * np.random.rand(N,1) + xmin
        Y = abs(ymax - ymin) * np.random.rand(N,1) + ymin
        many_points = np.column_stack((X,Y))
        many_points = self.points_in(many_points, poly=poly, IN=IN)
        
        return many_points
        
    def rand_shape(self, shape=None, N=1e4, IN=True):
        """
        Function that takes an input shapely Polygon object (or the default)
        and generates N random points within the bounding box around it. 
        Then M unknown points are returned that ARE contained within the
        Polygon object. This is done for speed. If IN=True then the function 
        will return points inside the matplotlib Path object, else 
        if IN=False then the function will return the points outside
        the matplotlib Path object. 
        """
        if shape is None:
            # Generates shape object from shape_file input
            INSHAPE = InShape(self.boundary)
            shape = self.node_poly(INSHAPE.external_coords())
            xmin, xmax, ymin, ymax = INSHAPE.shape_bounds()
        poly = self.node_poly(SHAPE.external_coords(shape=shape))
        points = self.rand_poly(poly=poly, N=N, IN=IN)
        return points
        
class Geodesic:
    """
    Class defined in order to create to process points, distances and
    other related geodesic calculations and functions 
    """    
    
    def __init__(self, period_range=[1, 40], km_point=20., max_dist=2e3):
        # initialise period_range as [1,40] default for ambient noise
        self.per_range = period_range
        self.km = km_point
        self.max_dist = max_dist
        
    def remove_distance(self, period_range, max_dist=None):
        """
        Function that returns a given possible resolvable ambient noise
        structure distance range, given the maximum period range
        availabe to the study. The distance returned is in km.
        
        Maximum distance default can be reassigned based on the cut-off found
        by your time-lag plots for your study! 
        """
        if max_dist is None: 
            max_dist = self.max_dist
        
        if type(period_range) == list:
            min_dist = min(period_range) * 9
            return [min_dist, max_dist]
        
        elif type(period_range) == int or float:
            return [period_range*9, max_dist]

    def haversine(self, lon1, lat1, lon2, lat2, R=6371):
        """
        Calculate the great circle distance between two points 
        on the earth (specified in decimal degrees). R is radius of
        spherical earth. Default is 6371km.
        """
        # convert decimal degrees to radians 
        lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
        # haversine formula 
        dlon, dlat = lon2 - lon1, lat2 - lat1
        a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
        c = 2 * asin(sqrt(a)) 
        km = R * c
        return km

    def fast_geodesic(self, lon1, lat1, lon2, lat2, npts):
        """
        Returns a list of *npts* points along the geodesic between
        (and including) *coord1* and *coord2*, in an array of
        shape (*npts*, 2).
        @rtype: L{ndarray}
        """
        if npts < 2:
            raise Exception('nb of points must be at least 2')

        path = wgs84.npts(lon1=lon1, lat1=lat1,
                          lon2=lon2, lat2=lat2,
                          npts=npts-2)
                      
        return np.array([[lon1,lat1]] + path + [[lon2,lat2]])

    def paths_calc(self, path_info, km_points=None, per_lims=None):
        """
        Function that returns an array of coordinates equidistant along 
        a great cricle path between two lat-lon coordinates if these points
        lay within a certain distance range ... otherwise the points return
        only a set of zeros the same size as the array. Default is 1.0km 
        distance per point.
        """
        if per_lims is None:
            # if no new default for period limits is defined, then set the
            # limit to the default.
            per_lims = self.per_range
            
        if km_points is None: 
            km_points = self.km
        
        lon1, lat1, lon2, lat2 = path_info[0], \
        path_info[1], path_info[2], path_info[3]    
        
        # interpoint distance <= 1 km, and nb of points >= 100
        dist = self.haversine(lon1, lat1, lon2, lat2)    
        npts = max(int((np.ceil(dist) + 1) / km_points), 100)    
        path = self.fast_geodesic(lon1, lat1, lon2, lat2, npts)
        dist_range = self.remove_distance(per_lims)
        
        if min(dist_range) < dist < max(dist_range):
        #remove the closest points along this line that fall below the distance 
        #find the index of the first point that is above this distance away!
            pts_km = npts / float((np.ceil(dist) + 1)) #this gives pts/km
            #remove all points below this index in the paths list
            dist_index = pts_km * min(dist_range) 
            path = path[dist_index:]
            return path    
        else:
            return np.zeros_like(path)
        
    def fast_paths(self, coord_list):
        """
        Function that takes many point coordinate combinations and quickly
        passes them through the paths_calc function. coord_list MUST be
        of the shape (4, N) whereby each coordinate combination is in a 
        (4,1) row [lon1,lat1,lon2,lat2]. 
        """
        return map(self.paths_calc, coord_list)
        
    def combine_paths(self, paths):
        """
        Function that takes many paths (should be array of same length as 
        number of stations). This is automatically generated by parallelising
        the fast_paths function above. 
        
        The output array should only contain unique, no repeating paths 
        and should be of the shape (2,N) where N is a large number of coords.
        """
        
        #create a flattened numpy array of size 2xN from the paths created! 
        paths = list(itertools.chain(*paths))
        paths = np.asarray(list(itertools.chain\
                                    (*paths)))
    
        #keep all but the repeated coordinates by keeping only unique whole rows!
        b = np.ascontiguousarray(paths).view(np.dtype\
        ((np.void, paths.dtype.itemsize * \
        paths.shape[1])))
        _, idx = np.unique(b, return_index=True)
        
        paths = np.unique(b).view(paths.dtype)\
        .reshape(-1, paths.shape[1])
        
        return paths
        
    def remove_zeros(self, paths):
        """
        Function that processes the flattened path output from combine_paths
        and removes the zero paths created by paths_calc. Remove zeroes 
        from paths to ensure all paths that were NOT in the distance threshold 
        are removed from the path density calculation! 
        """
    
        path_lons, path_lats = paths[:,0], paths[:,1]
                                               
        FIND_ZERO1 = np.where(paths[:,0]==0)[0]
        FIND_ZERO2 = np.where(paths[:,1]==0)[0]
        if len(FIND_ZERO1) != 0 and len(FIND_ZERO2) != 0:
            path_lons = np.delete(path_lons, FIND_ZERO1)
            path_lats = np.delete(path_lats, FIND_ZERO2)
        
        return np.column_stack((path_lons, path_lats))

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
    

GEODESIC = Geodesic()

    
# import background shapefile location
shape_path = "/home/boland/Dropbox/University/UniMelb\
/AGOS/PROGRAMS/ANT/Versions/26.04.2015/shapefiles/aus.shp"
INPOLY =  InPoly(shape_path)




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


pickle_file0 =  '/storage/ANT/spectral_density/station_pds_maxima/\
AUSTRALIA 2014/noiseinfo_comb.pickle'

pickle_file0 = '/storage/ANT/spectral_density/station_pds_maxima/AUSTRALIA 2014/first_peak_dict_australia_2014.pickle'

pickle_file0 = '/storage/ANT/spectral_density/noise_info0.pickle'

comb_noise = '/storage/ANT/spectral_density/station_pds_maxima/total_noise_combination.pickle'


f = open(name=comb_noise, mode='rb')
noise_info0 = pickle.load(f)
f.close()





# sort the noise
noise_info0 = np.asarray(noise_info0)
#noise_info0 = noise_info0[np.argsort(noise_info0[:, 1])]
# Combine AU with S info

print len(noise_info0)
                   
# find outliers
def reject_outliers(data, m=0.5):
    return data[abs(data - np.mean(data)) < m * np.std(data)] 
outliers = reject_outliers(noise_info0[:,2])

# remove outliers
noise_info0 = np.asarray([info for info in noise_info0 \
                          if info[2] in outliers])
                              
                              
# filter coordinates that are too close together. 

min_dist = 1. #degrees

coords = np.column_stack((noise_info0[:,0], noise_info0[:,1]))


# next remove points outside of the given poly if applicable
coord_indices = INPOLY.points_in(coords, indices=True)
noise_info0 =  noise_info0[coord_indices == True]
print noise_info0


coords = np.column_stack((noise_info0[:,0], noise_info0[:,1]))

coord_combs = np.asarray(list(itertools.combinations(coords,2)))

print len(coord_combs)

def coord_combinations(coord_combs):
    lon1, lat1 = coord_combs[0][0], coord_combs[0][1]
    lon2, lat2 = coord_combs[1][0], coord_combs[1][1]

    return [coord_combs, GEODESIC.haversine(lon1, lat1, lon2, lat2)]
    
t0 = datetime.datetime.now()
pool = mp.Pool()    
comb_dists = pool.map(coord_combinations, coord_combs)
pool.close()
pool.join()
t1 = datetime.datetime.now()
print t1-t0


comb_dists = np.asarray(comb_dists)

# sort by distance
comb_dists = comb_dists[np.argsort(comb_dists[:, 1])]

# find where the distances are less than the min_dist
find_min = np.where(comb_dists[:,1]>min_dist)[0]
# remove points where the distances are less than the min_dist
comb_dists = np.delete(comb_dists, find_min, axis=0)

remaining_coords = comb_dists[:,0]

# get unique coordinates from remaining coords

#paths = list(itertools.chain(*paths))
remaining_coords = np.asarray(list(itertools.chain\
(*remaining_coords)))


#keep all but the repeated coordinates by keeping only unique whole rows!
b = np.ascontiguousarray(remaining_coords).view(np.dtype\
     ((np.void, remaining_coords.dtype.itemsize * \
     remaining_coords.shape[1])))
_, idx = np.unique(b, return_index=True)
        
remaining_coords = np.unique(b).view(remaining_coords.dtype)\
        .reshape(-1, remaining_coords.shape[1])
 
 


 
#scan for all points that are within a degree radius of one another! 
db = DBSCAN(eps=min_dist).fit(coords)
core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
core_samples_mask[db.core_sample_indices_] = True
labels = db.labels_
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

# Black removed and is used for noise instead.
unique_labels = set(labels)

clusters = []
cluster_keep = []

for k in unique_labels:
    if k != -1:
        class_member_mask = (labels == k)        
        cluster = coords[class_member_mask & core_samples_mask]
        #xy = coords[class_member_mask & ~core_samples_mask]
    # Select only 1 random point from each cluster to keep. Remove all others!          
        clusters.append(cluster)
        cluster_keep.append(cluster[random.randint(0,len(cluster)-1)])
        
        
cluster_keep = np.asarray(cluster_keep)

# flatten clusters array 
clusters = np.asarray(list(itertools.chain(*clusters)))

# remove all points in clusters from the overall coords array
coords = np.asarray([coord for coord in coords if coord not in clusters])

# place single representative point from cluster back into overall coord list
coords = np.append(coords, cluster_keep, axis=0)

print len(noise_info0)

# remove cluster coordinates from noise_info0
noise_info0 = np.asarray([info for info in noise_info0 \
                          if info[0] in coords[:,0]])
           


fig = plt.figure(figsize=(15,10), dpi=1000)
plt.title('Average Seismic Noise First Peak Maximum PDS\n Australian Network | 2014')
plt.xlabel('Longitude (degrees)')
plt.ylabel('Latitude (degrees)')

print "number of station points: ", len(noise_info0)
patch = PolygonPatch(UNIQUE_SHAPE, facecolor='white',\
                     edgecolor='k', zorder=1)
ax = fig.add_subplot(111)
ax.add_patch(patch)


x, y = noise_info0[:,0], noise_info0[:,1]
points = np.column_stack((x,y))

xmin, xmax = np.min(x), np.max(x)
ymin, ymax = np.min(y), np.max(y)


values = noise_info0[:,2]

#now we create a grid of values, interpolated from our random sample above
y = np.linspace(ymin, ymax, 200)
x = np.linspace(xmin, xmax, 200)
gridx, gridy = np.meshgrid(x, y)
heat_field = griddata(points, values, (gridx, gridy), 
                      method='cubic',fill_value=0)


#heat_field = np.where(heat_field < 0, 1, heat_field)
heat_field = np.ma.masked_where(heat_field==0,heat_field)  

print gridx

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
            
            
fig.savefig('station_pds_maxima/noise_map_all.svg', 
            format='SVG')

