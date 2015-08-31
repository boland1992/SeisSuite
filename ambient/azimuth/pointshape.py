# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 12:28:32 2015

@author: boland
"""

import sys
sys.path.append("/home/boland/Anaconda/lib/python2.7/site-packages")
import fiona
import shapefile
from shapely import geometry
import numpy as np
import matplotlib.pyplot as plt
import pyproj
import datetime
from matplotlib.path import Path
#---------------------------------------------
#DEFINE INPUT PARAMETERS
#---------------------------------------------
#enter shapefile absolute or relative path name as string if optimal = True
shape_path = "/home/boland/Dropbox/University/UniMelb/AGOS/PROGRAMS/ANT/Versions/26.04.2015/shapefiles/aus.shp"
N = 130

def shape_(input_shape):
    
    with fiona.open(input_shape) as fiona_collection:
    
        # In this case, we'll assume the shapefile only has one record/layer (e.g., the shapefile
        # is just for the borders of a single country, etc.).
        shapefile_record = fiona_collection.next()
    
        # Use Shapely to create the polygon
        shape = geometry.asShape( shapefile_record['geometry'] )
    
    return shape


def points_in_shape(shape_path, N):

    shape = shape_(shape_path)
    
    minx, miny, maxx, maxy = shape.bounds
    
    #print minx; print miny; print maxx; print maxy
    #bounding_box = geometry.box(minx, miny, maxx, maxy)

    #generate random points within bounding box! 
    
    N_total = 130**2


    sf = shapefile.Reader(shape_path)
    shape = sf.shapes()[0]

    #find polygon nodes lat lons
    verticies = shape.points
    
    #convert to a matplotlib path class!
    polygon = Path(verticies)

    #points_in_shape = polygon.contains_points(coords)
    #coords = coords[points_in_shape == True][0:N-1]

    X = abs(maxx - minx) * np.random.rand(N_total,1) + minx
    
    Y = abs(maxy - miny) * np.random.rand(N_total,1) + miny

    coords = np.column_stack((X,Y))
    
    points_in_shape = polygon.contains_points(coords)
    
    coords = coords[points_in_shape == True][0:N]
    
    #indices = range(0,len(coords))
    #coords = np.column_stack((coords,indices))

    return coords



def paths_in_shape(paths):

    shape = shape_(shape_path)
    
    minx, miny, maxx, maxy = shape.bounds
    
    #print minx; print miny; print maxx; print maxy
    #bounding_box = geometry.box(minx, miny, maxx, maxy)

    #generate random points within bounding box! 
    
    sf = shapefile.Reader(shape_path)
    
    shape = sf.shapes()[0]

    #find polygon nodes lat lons
    verticies = shape.points
    
    #convert to a matplotlib path class!
    polygon = Path(verticies)
    points_in_shape = polygon.contains_points(paths)
    
    points_in_shape = paths[points_in_shape == True]
    
    return points_in_shape