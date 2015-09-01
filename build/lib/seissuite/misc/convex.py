# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 08:48:56 2015

@author: boland
"""

import matplotlib.pylot as plt
from scipy.spatial import ConvexHull
import numpy as np

class Convex:
    """
    CLASS CURRENTLY NOT WORKING!
    Class defined in order to create a convex hull around an array of points
    and then perform functions on them e.g. produce K random points inside,
    find N equidistant points within etc.
    """    
    def __init__(self, points):
        # initialise points of interest        
        self.dots = points
        # initialise polygon for potential convex hull with matplotlib Path
        self.polygon = 0.
        # initialise output points for points inside convex hull!
        self.output = 0.

    def convex_hull(self, point_set):
        """
        Function to produce a convex hull object that surrounds 
        a set of points. The input must be of the type Nx2 matrix/
        numpy array or equivalent. new_point shape is (2,1)
        """
        return ConvexHull(point_set)
            
    def poly_hull(self):
        """
        Function that generates a matplotlib Path object from convex hull 
        nodes. 
        """
        hull = self.convex_hull(self.dots)
        X, Y = self.dots[hull.vertices,0], self.dots[hull.vertices,1]        
        self.polygon = Path(np.column_stack((X, Y)))
        return self.polygon
        
    def in_poly_hull(self, point_set):
        """
        Function that quickly returns (2,N) array from (2,M) array of
        input points such that M >= N and N points are contained within
        the self.polygon polygon. 
        """
        self.polygon =  self.poly_hull()
        points_in = self.polygon.contains_points(point_set)
        self.output = point_set[points_in == True]
        return self.output
    
    def plot_hull(self, show_points=False):
        """
        Function that plots the boundaries of a convex hull using 
        matplotlib.pyplot. Input hull must be of type:
        scipy.spatial.qhull.ConvexHull 
            
        points input must be of the original coordinates.
        """
        hull = self.convex_hull(self.dots)
        plt.figure()
        for simplex in hull.simplices:
            plt.plot(self.dots[simplex,0], \
            self.dots[simplex,1], 'k-')
        if show_points:
            plt.scatter(self.dots[:,0], \
            self.dots[:,1], s=10,c='g')
            plt.scatter(self.dots[:,0], \
            self.dots[:,1], s=30,c='orange')
            plt.show()
    
    def rand_hull(hull, points, K):
        "Generate K new random points contained within a convex hull"

        minx, maxx = np.min(points[:,0]), np.max(points[:,0])  
        miny, maxy = np.min(points[:,1]), np.max(points[:,1])
        X = abs(maxx - minx) * np.random.rand(10*K**2,1) + minx
        Y = abs(maxy - miny) * np.random.rand(10*K**2,1) + miny
        return np.column_stack((X,Y))
