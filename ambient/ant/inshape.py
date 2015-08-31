# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 08:42:50 2015

@author: boland
"""
#------------------------------------------------------------------------------
# MODULES
#------------------------------------------------------------------------------
import fiona
import numpy as np
from shapely import geometry
from shapely.geometry import asPolygon

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
            self.polygon = geometry.asShape(shapefile_record['geometry'] )
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