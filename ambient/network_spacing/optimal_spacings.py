# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 13:30:21 2015

@author: boland

CODE DESCRIPTION:
The following python script is used to import an irregular polygon, in the
form of an imported shapefile (.shp) and then determine where N means-squared
evenly distributed (distance-wise) points are within this polygon.
"""

import fiona
from shapely import geometry
import numpy as np
from scipy.cluster.vq import kmeans
import multiprocessing as mp

#input_shape = "/home/boland/Dropbox/University/UniMelb/AGOS/PROGRAMS/ANT/Versions/26.04.2015/shapefiles/aus.shp"
#input_shape = '/home/boland/Downloads/NZ Shapefiles/clipped_NZ.shp'
#quit()

def shape_(input_shape):
    
    with fiona.open(input_shape) as fiona_collection:
    
        # In this case, we'll assume the shapefile only has one record/layer (e.g., the shapefile
        # is just for the borders of a single country, etc.).
        shapefile_record = fiona_collection.next()
    
        # Use Shapely to create the polygon
        shape = geometry.asShape( shapefile_record['geometry'] )
    
    return shape


def many_points(shape):
    
    """
    Funtion that returns lat-lon coodinates of many random points within 
    a shapefile shape e.g. boundaries of a state or country.
    """    
    minx, miny, maxx, maxy = shape.bounds
    
    #print minx; print miny; print maxx; print maxy
    bounding_box = geometry.box(minx, miny, maxx, maxy)

    #generate random points within bounding box! 
    n = 1e4
    
    X = abs(maxx - minx) * np.random.rand(n,1) + minx
    
    Y = abs(maxy - miny) * np.random.rand(n,1) + miny
        
    
    coords = np.column_stack((X,Y))
    
    coord_points = []

    pool = mp.Pool(8)
    coord_points = pool.map(point_check, coords)
    pool.close()
    pool.join()
    
    #remove None values from above numpy array! 
    X = [i[0] for i in coord_points if i != None]
    Y = [i[1] for i in coord_points if i != None]
    
    #convert np array to kmeans function friendly nx2 matrix!     
    coord_points = np.column_stack((X,Y))
    
    #output is a nx2 vector stacked matrix called coord_points, 
    #coord_points[:,0] are lons and coord_points[:,1] are lats

    return coord_points

def cluster_points(coord_points):
    """
    Function that returns k which is an nx2 matrix of lon-lat vector columns
    containing the optimal cluster centroid spacings within a large set of random
    numbers e.g. those produced by the many_points() function above!
    """

    k = kmeans(coord_points, 100)
    
    return k
    
