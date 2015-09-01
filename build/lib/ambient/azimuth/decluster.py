# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 08:44:50 2015

@author: boland
"""
#------------------------------------------------------------------------------
# MODULES
#------------------------------------------------------------------------------
import os
import fiona
import pickle
import pyproj
import datetime
import itertools
import shapefile
import numpy as np
import datetime as dt
import pointshape as ps
from shapely import geometry
import multiprocessing as mp
import matplotlib.pyplot as plt
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
# VARIABLES
#------------------------------------------------------------------------------

verbose = False
#Enter path to boundary shape file.
shape_boundary = False
#shape_path = "/home/boland/Dropbox/University/UniMelb\
#/AGOS/PROGRAMS/ANT/Versions/26.04.2015/shapefiles/aus.shp"
dataless = True
# Enter path to dataless file
#dataless_path = 'ALL_AUSTRALIA.870093.dataless'
#dataless_path = 'UOM.dataless'
# Enter number new stations desired.
N = 3
# Enter km spacing between path density points.
km_points = 20.0
# Reference elipsoid to calculate distance.
wgs84 = pyproj.Geod(ellps='WGS84')
# Enter number of bins for 2D Histogram density calculation. 
nbins = 220
# Enter estimated average shear wave velocity. 3kms-1 is the default!
velocity = 3.0
# Define your ambient noise period range OR individual period in seconds.
global period_range
period_range = [1,40]

dataless_path = 'east-timor/timor.dataless'
dataless_path = '/storage/ANT/spectral_density/USARRAY/full_USARRAY.dataless'


coords = locs_from_dataless(dataless_path)


#shape_path = "/home/boland/Dropbox/University/UniMelb\
#/AGOS/PROGRAMS/ANT/Versions/26.04.2015/shapefiles/aus.shp"

shape_path = 'east-timor/TLS_adm0.shp'


coords = locs_from_dataless(dataless_path)

#shape_path = "/home/boland/Dropbox/University/UniMelb\
#/AGOS/PROGRAMS/ANT/Versions/26.04.2015/shapefiles/aus.shp"


t0 = dt.datetime.now()

# Generate InShape class
SHAPE = InShape(shape_path)
# Create shapely polygon from imported shapefile 
UNIQUE_SHAPE = SHAPE.shape_poly()
print type(UNIQUE_SHAPE)
# Generate InPoly class
INPOLY = InPoly(shape_path)
# Create matplotlib Path object from imported shapefile
#outer_shape = UNIQUE_SHAPE.buffer(1.,resolution=1)
#inner_shape = UNIQUE_SHAPE.buffer(-8,resolution=1)

#outer_poly = INPOLY.poly_from_shape(shape=outer_shape)
#inner_poly = INPOLY.poly_from_shape(shape=inner_shape)
#many_points = INPOLY.rand_poly(poly=outer_poly, N=1e4)

# Scale smaller shape to fit inside larger shape. 
#SMALL_SHAPE = scale(UNIQUE_SHAPE, xfact=0.3, yfact=0.3)
#points_in_small_shape = INPOLY.rand_shape(shape=SMALL_SHAPE, IN=False)
# Generate matplotlib Path object for the small scalled polygon 
#small_poly = INPOLY.node_poly(SHAPE.external_coords(shape=SMALL_SHAPE))
# Remove points that are outside the buffered_poly
#outer_poly_points = INPOLY.points_in(many_points, poly=outer_poly)

# Remove points that are inside the small_poly
#inner_poly_points = np.asarray(INPOLY.points_in(outer_poly_points, 
#                                                poly=inner_poly,
#                                                IN=False))

#cluster_points = np.asarray(kmeans(inner_poly_points, 130)[0])


#plt.figure()
#plt.scatter(inner_poly_points[:,0], inner_poly_points[:,1], c='b')
#plt.scatter(cluster_points[:,0], cluster_points[:,1], c='orange', s=35)
#plt.show()

#-----------------------------------------------------------------------------
# INITIALISE CLASS STATES
#-----------------------------------------------------------------------------
GEODESIC = Geodesic()
COORDS = Coordinates()
INPOLY = InPoly(shape_path)
POLY_NODES = INPOLY.poly_nodes()
#-----------------------------------------------------------------------------
# GENERATE SECOND SET OF VARIABLES AND STATES
#-----------------------------------------------------------------------------
ideal_path = 'ideal_coordinates.pickle'
#if no paths have been done before, start afresh!
#if dataless:
#    coords = locs_from_dataless(dataless_path)
#    original_coords = coords
#elif os.path.exists(ideal_path):
#    f = open(name=ideal_path, mode='rb')
#    coords = pickle.load(f)
#    f.close()

    
# decluster the points to desired specifications. 
coords = COORDS.decluster(inputs=coords, degree_dist=0.5)


lonmin, lonmax = np.floor(min(coords[:,0])), np.ceil(max(coords[:,0]))
latmin, latmax = np.floor(min(coords[:,1])), np.ceil(max(coords[:,1]))
print lonmin,lonmax,latmin,latmax


plt.figure()
plt.scatter(coords[:,0], coords[:,1])
plt.show()



kappa = [np.vstack([[coord1[0],coord1[1],coord2[0],coord2[1]]\
                    for coord2 in coords]) for coord1 in coords]                        


def spread_paths(coord_list):
    return GEODESIC.fast_paths(coord_list)
    
t0 = datetime.datetime.now()
pool = mp.Pool()    
paths = pool.map(spread_paths, kappa)
pool.close()
pool.join()
t1 = datetime.datetime.now()
print t1-t0


counter, counter2 = 0, 0
#cd Desktop/Link\ to\ SIMULATIONS/Network_Tracks/smarter_model/
grad_ideal, grad_check1, grad_check2, H_avg1, H_avg2 = 0, 0, 0, 0, 0
SHAPE = (1,1)
perc_high = 0.01
low_counter = 0
random_counter = 0
#new_coord = 0
infinite_counter = 0
find_it = []
check_coord = None
use_old_path = False
searches_per_point = 3
factor = 0.05
cluster = False
N_cluster_points = False

while infinite_counter <= 1:
    t0 = datetime.datetime.now()
    
    #----------------------------------------------------------------------
    # Generate N new point coordinates
    #----------------------------------------------------------------------
    if cluster:
        new_coords = N_cluster_points
    else:
        new_coords = ps.points_in_shape(shape_path, N)
        
    coords = np.append(coords, new_coords, axis=0)

    coord_set = [np.vstack([[coord1[0],coord1[1],coord2[0],coord2[1]]\
                 for coord2 in coords]) for coord1 in coords]

    t0 = datetime.datetime.now()
    pool = mp.Pool()    
    paths = pool.map(spread_paths, coord_set)
    pool.close()
    pool.join()
    t1 = datetime.datetime.now()
    print "time to generate new paths", t1-t0
    
    # Append new set of paths now that old set has been deleted.
    
    #create a flattened numpy array of size 2xN from the paths created! 
    paths1 = GEODESIC.combine_paths(paths)

    paths = list(paths)

    paths1 = GEODESIC.remove_zeros(paths1)

    DENSITY = Density(paths=paths1)

    H, xedges, yedges = DENSITY.hist2d(paths=paths1)
    grad = DENSITY.hgrad(H=H)
    
    H_avg1 = np.average(H)
    grad_check1 = np.std(grad)
    
    H_masked = DENSITY.transform_h(H=H)
    grad = DENSITY.transform_grad(grad=grad)

    
    search = np.where(H<0.1*np.average(H))
    Hmaxx, Hmaxy =  search[1], search[0]    
    Hmaxx = (lonmax-lonmin)/(nbins) * Hmaxx + lonmin
    Hmaxy = (latmax-latmin)/(nbins) * Hmaxy + latmin
    # Make sure all low density coordinates ARE within shapefile!
    low_density_coords = ps.paths_in_shape(np.column_stack((Hmaxx, Hmaxy)))
    
    #N_cluster_points = kmeans(low_density_coords, N)[0]
    
    
    density_coords = DENSITY.select_points()
    # make sure that your density coords are within the boundary shape        
    density_coords = INPOLY.points_in(density_coords)
    #cluster = True
    
    if counter == 0:
        grad_ideal = 1e6
        avg_ideal = 0  
    

    if grad_check1 < grad_ideal and avg_ideal < H_avg1:     

        with open(u'ideal_coordinates.pickle', 'wb') as f:
            print "\nExporting new ideal coordinates."
            pickle.dump(coords, f, protocol=2)
        
        DENSITY.plot_field()#nodes=POLY_NODES)#SHAPE=UNIQUE_SHAPE)

        grad_ideal = grad_check1
        avg_ideal = H_avg1

    coords = COORDS.del_N(N=N, inputs=coords)
    paths = COORDS.del_N(N=N, inputs=paths)
    paths=list(paths)

    counter+=1
    t1 = datetime.datetime.now()
    print "That loop took: ", t1-t0