# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 08:45:03 2015

@author: boland
"""

import shapefile
import numpy as np
from matplotlib.path import Path
from inshape import InShape

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
        
    def points_in(self, points, poly=None, IN=True):
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
        poly = self.node_poly(shape.external_coords(shape=shape))
        points = self.rand_poly(poly=poly, N=N, IN=IN)
        return points