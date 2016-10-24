# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 08:44:50 2015

@author: boland

CODE DESCRIPTION: 
The following python script searches for M new random points atop N set station
points. The tests performed have to do with point density distribution of
points representing all combinations of great-circlel lines that could 
be physically possible between seismic stations. An extension is to select
low point density points as a new cluster to search for new station points. 
"""
#------------------------------------------------------------------------------
# MODULES
#------------------------------------------------------------------------------
import os
#import fiona
import pysal

import pickle
import pyproj
import datetime
import itertools
import shapefile
import numpy as np
import datetime as dt
import multiprocessing as mp
import matplotlib.pyplot as plt
from math import sqrt, atan2, asin, degrees, radians, tan, sin, cos
from shapely.geometry import asPolygon, Polygon
from descartes.patch import PolygonPatch
from scipy.spatial import ConvexHull
from scipy.cluster.vq import kmeans
from shapely.affinity import scale
from matplotlib.path import Path
from shapely import geometry


#------------------------------------------------------------------------------
# VARIABLES
#------------------------------------------------------------------------------

# Reference elipsoid to calculate distance.
wgs84 = pyproj.Geod(ellps='WGS84')

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
        #with fiona.open(self.boundary) as fiona_collection:
            # In this case, we'll assume the shapefile only has one later
       #     shapefile_record = fiona_collection.next()
            # Use Shapely to create the polygon
       #     self.polygon = geometry.asShape( shapefile_record['geometry'] )
       #     return self.polygon
        # Now, open the shapefile using pysal's FileIO
        shps = pysal.open(self.boundary , 'r')
        poly = shps.next()
        self.polygon = geometry.asShape(poly)
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
        new_coords = np.column_stack((X,Y))
        return new_coords

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
        
    def fast_paths(self, coord_list, km_points=None, per_lims=None):
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


class Coordinates:
    """
    Class defined in order to perform to latitude, longitude coordinates
    operations.
    """
    def __init__(self, input_list=None, N=None):
        # initialise input list of a (2,N) numpy array
        self.input_list = input_list
        # initialise import number
        self.N = N

    def del_N(self, N=None, inputs=None):
        """
        Function that deletes the last N coordinates from a list of coordinate
        """
        if N is None:
            if self.N is not None:
                N = self.N
            elif self.N is None:
                raise "There are no number input. Please enter a desired\
                number of points to remove from the input_list!" 
        if inputs is None:
            if self.input_list is not None:
                inputs = self.input_list
            elif self.input_list is None:
                raise "There are no list input. Please enter a desired\
                list of points to remove N number of points from the end!" 
        if not type(inputs) == 'list':
            inputs = list(inputs)
        del inputs[-N:]
        return np.asarray(inputs)
        

    def decluster(self,  inputs=None, degree_dist=1., verbose=False):
        """
        Function that deletes points that are too close together
        given a set degree range and returns only one point to represent
        that cluster. Default is one degree distance. Inputs must be (2,N)
        lon-lat coordinate arrays/lists.
        """
        from sklearn.cluster import DBSCAN
        import random 
        
        if inputs is None:
            if self.input_list is not None:
                inputs = self.input_list
            elif self.input_list is None:
                raise "There are no list input. "
                
        
        #scan for all points that are within a degree radius of one another! 
        db = DBSCAN(eps=degree_dist).fit(inputs)
        core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True
        labels = db.labels_
        unique_labels = set(labels)

        clusters = []
        cluster_keep = []
        

        for k in unique_labels:
            if k != -1:
                class_member_mask = (labels == k)        
                cluster = inputs[class_member_mask & core_samples_mask]

    # Select only 1 random point from each cluster to keep. Remove all others!          
                clusters.append(cluster)
                cluster_keep.append(cluster[random.randint(0,len(cluster)-1)])
        
        
        cluster_keep = np.asarray(cluster_keep)
        # flatten clusters array 
        clusters = np.asarray(list(itertools.chain(*clusters)))
        
        # remove all points in clusters from the overall coords array
        inputs = np.asarray([point for point in inputs if 
                             point not in clusters])
                             
        if verbose:                     
            print "clusters array shape: ", clusters.shape
            print "inputs array shape: ", inputs.shape
            print "cluster_keep array shape: ",cluster_keep.shape
        
        if len(cluster_keep) > 0:
            output_coords = np.append(inputs, cluster_keep, axis=0)
        # place single representative point from cluster into coord list
        else:
            output_coords = inputs
        
        return output_coords
        
class Density:
    """
    Class defined to perform to density field operations e.g. 2d histogram, 
    gradient fields and averages and standard deviations thereof. 
    """
    def __init__(self, paths=None, nbins=200):
        # initialise path points for shape (2,N) to calculate point density
        self.paths = paths
        # initialise the number of bins per axis for 2d histogram
        self.nbins = nbins
        # initialise the density calculation
        self.H = 0.
        self.H_masked = 0.
        # initialise the density gradient field calculation
        self.grad = 0.
        self.grad_masked = 0.
        # initialise x and y tridiagonal coordinate matrices
        self.xedges = 0.
        self.yedges = 0.
        
    def hist2d(self, paths=None):
        """
        Function that calculates the 2D histogram and generates H, xedges, 
        and yedges. 
        """
        if paths is None:
            paths = self.paths
            
        self.H, self.xedges, self.yedges = np.histogram2d(paths[:,0],
                                           paths[:,1],
                                           bins=self.nbins)
        return self.H, self.xedges, self.yedges
            
    def hgrad(self, H=None):
        """
        Function that calculates the 2D histogram and generates H, xedges, 
        and yedges. 
        """
        #, xedges=None, yedges=None
        if H is None: 
            H, xedges, yedges = self.hist2d()
            
        self.grad = np.abs(np.asarray(np.gradient(H)[0]))
        
        return self.grad
            
    def transform_h(self, H=None):
        """
        Function that rotates, flips and masks the H density field 
        in order for it to be plotted etc.
        """
        if H is None:
            if all( [self.H != 0., self.xedges != 0, self.yedges != 0] ):
                H, xedges, yedges = self.H, self.xedges, self.yedges
            else:
                H, xedges, yedges = self.hist2d()
            
        H = np.rot90(H)
        H = np.flipud(H)
        self.H_masked = np.ma.masked_where(H==0,H)  
        return self.H_masked
        
    def transform_grad(self, grad=None):
        """
        Function that rotates, flips and masks the H density gradient field 
        in order for it to be plotted etc.
        """
        if grad is None: 
            grad, xedges, yedges = self.hgrad()
        grad = np.rot90(grad)
        grad = np.flipud(grad)
        self.grad_masked = np.ma.masked_where(grad==0,grad)
        return self.grad_masked
    

    def plot_lims(self, paths=None):
        if paths is None:
            try:
                lons, lats = self.paths[:,0], self.paths[:,1]
            except Exception as error:
                raise error
        else:
            try:
                lons, lats = paths[:,0], paths[:,1]
            except Exception as error:
                raise error
            
        return np.min(lons), np.max(lons), np.min(lats), np.max(lats)
    
    def select_points(self, perc=0.1, high=None):
        """
        Function that returns the lat-lon coordinates of points below a 
        certain density. This is taken by H < perc*np.average(H) OR if
        high=True or is not None, then high density points are chosen from
        H >perc*np.average(H). Perc=0.1 by default.
        """
        H = self.H
        if high is None:
            search = np.where(H<perc*np.average(self.H))
        else:
            search = np.where(H>perc*np.average(self.H))
         
        xmin, xmax =  np.min(self.xedges), np.max(self.xedges)
        ymin, ymax =  np.min(self.yedges), np.max(self.yedges) 

        Hdensx, Hdensy =  search[1], search[0]    
        Hdensx = (xmax-xmin)/(self.nbins) * Hdensx + xmin
        Hdensy = (ymax-ymin)/(self.nbins) * Hdensy + ymin
        return np.column_stack((Hdensx, Hdensy))

    
    def plot_field(self, grad=False, SHAPE=None, swell=0.00):
        
        lonmin, lonmax, latmin, latmax = self.plot_lims()
        
        fig = plt.figure(figsize=(15,10), dpi=100)
        
        plt.xlabel('longitude (degrees)')
        plt.ylabel('latitude (degrees)')
        plt.xlim(lonmin-swell*abs(lonmax-lonmin),\
                 lonmax+swell*abs(lonmax-lonmin))
        plt.ylim(latmin-swell*abs(latmax-latmin),\
                 latmax+swell*abs(latmax-latmin))
    
        if not grad:
            plt.title("Path Density Distribution")
            if self.H_masked is not 0.:
                H_masked = self.H_masked
            else:                    
                H_masked = self.transform_h()
            plt.pcolor(self.xedges, self.yedges, H_masked, norm=LogNorm(\
            vmin=np.min(H_masked), vmax=np.max(H_masked)), cmap='rainbow',\
            alpha=0.6, zorder = 3)

            col = plt.colorbar()
            col.ax.set_ylabel('Points Per Bin')
        elif grad:
            plt.title("Gradient Path Density Distribution")
            
            if self.grad_masked is not 0.:
                grad_masked = self.grad_masked
            else:           
                raise Exception("grad_masked has not yet been defined. please\
                run the necessary functions e.g. transform_grad before plotting")
                
            plt.pcolor(self.xedges, self.yedges, grad_masked, norm=LogNorm(\
            vmin=np.min(grad), vmax=np.max(grad)), cmap='rainbow',\
            alpha=0.6, zorder = 3)
            col = plt.colorbar()
            col.ax.set_ylabel('Gradient Points Per Bin')      
        else:
            raise Exception("Either you have not chosen to plot anything OR\n\
                   both H and grad are inputed and the function doesn't\
                   know what to do.")

        if SHAPE is not None:
            patch = PolygonPatch(SHAPE, facecolor='white',\
            edgecolor='k', zorder=1)
            ax = fig.add_subplot(111)
            ax.add_patch(patch)
            ax.set_xlim(lonmin-0.05*abs(lonmax-lonmin), \
                        lonmax+0.05*abs(lonmax-lonmin))
            ax.set_ylim(latmin-0.05*abs(latmax-latmin), \
                        latmax+0.05*abs(latmax-latmin))
                            
        #plt.scatter(new_coords[:,0], new_coords[:,1],c='r', s=30)
        fig.savefig("plot_density.png")
        fig.clf()


# The functions below are used to calculate the 

def latitude(dist, sigma01, alpha0, lon0):
    sigma = sigma01 + dist#/R
    lat = degrees(asin(cos(alpha0)*sin(sigma)))
    #alpha = atan2(tan(alpha0),cos(sigma))
    return lat
    
def longitude(dist, sigma01, alpha0, lon0):
    sigma = sigma01 + dist#/R
    lon = degrees(atan2(sin(alpha0)*sin(sigma), cos(sigma))) + degrees(lon0)
    #alpha = atan2(tan(alpha0),cos(sigma))
    return lon

vlat_func = np.vectorize(latitude)
vlon_func = np.vectorize(longitude)

def waypoint_init(path_info, km=20):
    R = 6371
    lon1, lat1, lon2, lat2, dist = radians(path_info[0]), \
    radians(path_info[1]), radians(path_info[2]), \
    radians(path_info[3]), radians(path_info[4])
    #lon1, lat1, lon2, lat2, dist = \
    #map(radians, [path_info[0],path_info[1],path_info[2],
        #path_info[3],path_info[4]])
    lon_diff = lon2-lon1
    alpha1 = atan2(sin(lon_diff),(cos(lat1)*tan(lat2)-sin(lat1)*cos(lon_diff)))
       #alpha2 = atan2(sin(lon_diff),(-cos(lat2)*tan(lat1)+sin(lat2)*cos(lon_diff)))
    #try:
        #sigma12 = acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon_diff))
    #except:
        #return
    sigma01, alpha0 = atan2(tan(lat1), cos(alpha1)), asin(sin(alpha1)*cos(lat1))
    #sigma02 = sigma01+sigma12
    lon01 = atan2(sin(alpha0)*sin(sigma01), cos(sigma01))
    lon0 = lon1 - lon01    
    npts = max(int((np.ceil(dist) + 1)/km), 100)
    all_d = np.linspace(0,dist,npts)/R  
    lons, lats = vlon_func(all_d, sigma01, alpha0, lon0), vlat_func(all_d, sigma01, alpha0, lon0)   
    
    return np.column_stack((lons, lats))
