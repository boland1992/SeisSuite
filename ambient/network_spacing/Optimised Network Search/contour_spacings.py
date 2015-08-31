# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 12:28:32 2015

@author: boland
"""

import sys
sys.path.append('/home/boland/Anaconda/lib/python2.7/site-packages')
import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.vq import kmeans
import multiprocessing as mp
import pyproj
import os
import itertools
import datetime
import pointshape as ps
from math import sqrt, radians, cos, sin, asin
import random
from shapely import geometry
import fiona






verbose = False
shape_path = "/home/boland/Dropbox/University/UniMelb/AGOS/PROGRAMS/ANT/\
Versions/26.04.2015/shapefiles/aus.shp"

N = 130
#enter km spacing between path density points
km_points = 100.0

# reference elipsoid to calculate distance
wgs84 = pyproj.Geod(ellps='WGS84')
nbins = 200

#enter estimated average shear wave velocity. 3kms-1 is the default!
velocity = 3.0
#define your ambient noise period range OR individual period in seconds
global period_range
period_range = [1,40]

#-----------------------------------------------------------------------------
#OLD SHAPEFILE FUNCTIONS
#-----------------------------------------------------------------------------

def shape_(input_shape):
    
    with fiona.open(input_shape) as fiona_collection:
    
        # In this case, we'll assume the shapefile only has one record/layer (e.g., the shapefile
        # is just for the borders of a single country, etc.).
        shapefile_record = fiona_collection.next()
    
        # Use Shapely to create the polygon
        shapes = geometry.asShape( shapefile_record['geometry'] )
        #x = [p.coords.xy[0] for p in points]
        #y = [p.coords.xy[1] for p in points]
        #plt.figure(figsize=(10,10))
        #plt.plot(x,y,'o', color='#f16824')
        #plt.show()
        return shapes

shape = shape_(shape_path)

def point_check(coord, shape=shape):#, coord_points=coord_points):
    
    point = geometry.Point(coord[0], coord[1]) # longitude, latitude
    if shape.contains(point):
    #X = np.delete(X, i); Y = np.delete(Y, i)
        #coord_points.append(coord)

        return coord   

def many_points(shape, coords):
    
    """
    Funtion that returns lat-lon coodinates of many random points within 
    a shapefile shape e.g. boundaries of a state or country.
    """    
    minx, miny, maxx, maxy = shape.bounds
    
    X, Y = coords[:,0], coords[:,1]
    
    coords = np.column_stack((X,Y))
    
    coord_points = []

    pool = mp.Pool()
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
#-----------------------------------------------------------------------------

def haversine(coordinates):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    lon1, lat1, lon2, lat2= coordinates[0],coordinates[1],\
    coordinates[2],coordinates[3]

    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    
    # haversine formula 
    dlon, dlat = lon2 - lon1, lat2 - lat1
    
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    km = 6367 * c
    
    return km
    
def haversine2(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    
    # haversine formula 
    dlon, dlat = lon2 - lon1, lat2 - lat1
    
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    km = 6367 * c
    
    return km

def geodesic(coord1, coord2, npts):
    """
    Returns a list of *npts* points along the geodesic between
    (and including) *coord1* and *coord2*, in an array of
    shape (*npts*, 2).
    @rtype: L{ndarray}
    """
    if npts < 2:
        raise Exception('nb of points must be at least 2')

    path = wgs84.npts(lon1=coord1[0], lat1=coord1[1],
                      lon2=coord2[0], lat2=coord2[1],
                      npts=npts-2)
    return np.array([coord1] + path + [coord2])
    
def new_geodesic(lon1,lat1,lon2,lat2, npts):
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

def cluster_points(coord_points, N):
    """
    Function that returns k which is an nx2 matrix of lon-lat vector columns
    containing the optimal cluster centroid spacings within a large set of random
    numbers e.g. those produced by the many_points() function above!
    """

    k = kmeans(coord_points, N)
    return k[0]
    
def remove_distance(period_range, max_dist = 2000):
    """
    Function that returns a given possible resolvable ambient noise
    tomography structure distance range, given the maximum period range
    availabe to the study. The distance returned is in km.
    
    Maximum distance default can be reassigned based on the cut-off found
    by your time-lag plots for your study!
    """
    
    if type(period_range) == list:
        min_dist = min(period_range) * 9
        return [min_dist, max_dist]
        
    elif type(period_range) == int or float:
        return [period_range*9, max_dist]

global dist_range

dist_range = remove_distance(period_range, max_dist = 2000)

def paths_func(path_info, km=km_points):

    lon1, lat1, lon2, lat2 = path_info[0], \
    path_info[1], path_info[2], path_info[3]    
    #lon1, lat1, lon2, lat2, dist = path_info[0], \
    #path_info[1], path_info[2], path_info[3],    \
    #path_info[4]  
  
    # interpoint distance <= 1 km, and nb of points >= 100
         
    dist = haversine2(lon1, lat1, lon2, lat2)
    
    npts = max(int((np.ceil(dist) + 1)/km), 100)
    
    path = new_geodesic(lon1,lat1,lon2,lat2, npts)

    lons, lats = [lon1 for i in range(0,len(path))], \
    [lat1 for i in range(0,len(path))]
    
    path = np.column_stack((path,lons,lats))
    
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
        
def points_and_distances(path_info):
    lon1, lat1, lon2, lat2 = path_info[0], \
    path_info[1], path_info[2], path_info[3]
    dist = haversine2(lon1, lat1, lon2, lat2)
    
    return [dist, lon2, lat2]



def HIST2D(nbins,paths, grad=False):
    H, xedges, yedges = np.histogram2d(paths[:,0],paths[:,1],bins=nbins)
    #name = "path_density_2Dhist.png"
    if grad:
        H = np.abs(np.asarray(np.gradient(H)[0]))#name = "path_density_2Dhist_grad.png"
    # H needs to be rotated and flipped
    H = np.rot90(H)
    H = np.flipud(H)
     # Mask zeros
    Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero  
    return Hmasked

    
t_total0 = datetime.datetime.now()


t0 = datetime.datetime.now()

ideal_path = 'ideal_coordinates.pickle'
#if no paths have been done before, start afresh!
if not os.path.exists(ideal_path):
    M = 1e5
    many_points = ps.points_in_shape(shape_path, M)
    coords = cluster_points(many_points,N)
    
#else import already processed coordinates if the program has already done so.
else:
    f = open(name=ideal_path, mode='rb')
    coords = pickle.load(f)
    f.close()

#generate N kmeans cluster points from massive M number of randomly distributed
#points inside the shape file. 

lonmin = np.floor(min(coords[:,0]))
lonmax = np.ceil(max(coords[:,0]))
latmin = np.floor(min(coords[:,1]))
latmax = np.ceil(max(coords[:,1]))

print lonmin,lonmax,latmin,latmax


#coords1 = [coord1 for coord1 in coords for coord2 in coords]
                     
#coords2 = [coord2 for coord1 in coords for coord2 in coords]

#columns = np.column_stack((coords1, coords2))


kappa = [np.vstack([[coord1[0],coord1[1],coord2[0],coord2[1]]\
                    for coord2 in coords]) for coord1 in coords]


def fast_dists(nets):
    dist_points = map(points_and_distances, nets)
    #print dist_points
    #find the point closest to the high density points in the array    
    #minimum_distance =  np.min(np.vstack(dist_points)[:,0])
    #print(len(dist_points))
    #find_it = np.where(dist_points==minimum_distance)[0][0]
    
    #return dist_points[find_it]
    return np.vstack(dist_points)



def spread_paths(nets):
    return map(paths_func, nets)
    
t0 = datetime.datetime.now()
pool = mp.Pool()    
paths = pool.map(spread_paths, kappa)
pool.close()
pool.join()
t1 = datetime.datetime.now()
print t1-t0
#paths = list(paths)

counter = 0

#cd Desktop/Link\ to\ SIMULATIONS/Network_Tracks/smarter_model/
grad_ideal, grad_check1, grad_check2, H_avg1, H_avg2 = 0, 0, 0, 0, 0
SHAPE = (1,1)
counter2 = 0
perc_high = 0.01

#counter of how many times the points 
#have been chosen from the lowest path density spots
low_counter = 0

#counter of how many times the points 
#have been chosen from the random spots.
random_counter = 0
new_coord = 0

infinite_counter = 0
find_it = []
check_coord = None
while infinite_counter <= 1:
    
    print "check coord", check_coord
    t0 = datetime.datetime.now()
    
    #the following while loop is a work around fix to a:
    #new paths shape: (130, 100, 4) rather than being (130,) like it should be!
    while SHAPE != (130,):
        #if counter2 >= len(paths)-1:
        #    counter2 = 0
        
        #cycle through paths 
        #----------------------------------------------------------------------
        #old_path = paths[counter2]
        #del paths[counter2]
        #old_coord = [old_path[0][0][0],old_path[0][0][1]]
        #itemindex = np.where(coords==old_coord)[0][0]
        
        #coords = list(coords)
        #find index of array in nested array to remove!
        #del coords[itemindex]
        #----------------------------------------------------------------------

        #or random selection of paths?!
        #----------------------------------------------------------------------
        if check_coord is None:

            while len(find_it) == 0:
            
                #remove a random set of paths associated with a single one of the N coordinates 
                rand_int = random.randint(0,len(paths)-1)
            
                old_path = paths[rand_int]
                #figure out which old coordinate to remove from the coordinates list
                old_coord = [old_path[0][0][0],old_path[0][0][1]]
                #print "old coord:", old_coord
                #NEED TO REMOVE OLD POINT FROM COORDS!
                #find index of array in nested array to remove!
                find_it = np.where(coords==old_coord)[0]
            
        else:
            #find high density coords
            while len(find_it) == 0:

                old_coord = check_coord
                find_it = np.where(coords==old_coord)[0]

                for paths_index, path in enumerate(paths):
                    for dots in path:
                        find_paths = np.where(dots==old_coord[0])[0]
                        if len(find_paths) != 0:
                            rand_int = False
                            break

                
        coord_index = find_it[0]
            
        print "coord index to remove is:", coord_index
        coords = list(coords)
        del coords[coord_index]

        #find index of array in nested array to remove!
        coords = np.asarray(coords)
        
        new_coord_first = new_coord

        #----------------------------------------------------------------------
        #generate new point coordinate
        if  counter <= 1:
            new_coord = ps.points_in_shape(shape_path, 1)[0]
        else:
            new_coord = new_coord
        
        #place new coordinate in old set of coordinates
        coords = np.append(coords, [new_coord], axis=0)

        #generate new array of points in conjunction with the new randomly generated point!
        new_coord_set = np.vstack([[new_coord[0],new_coord[1],coord1[0],\
                                coord1[1]] for coord1 in coords])
    
        #generate new random point in place of all 'popped' points!
        new_paths = map(paths_func, new_coord_set)
        
        SHAPE = np.asarray(new_paths).shape
        
        if not SHAPE == (130,):
            #remove substitute back the old coordinate for the new coordinate!
            coords = list(coords)
            #find index of array in nested array to remove!
            del coords[-1]
            coords = np.asarray(coords)
            #place new coordinate in old set of coordinates
            coords = np.append(coords, [old_coord], axis=0)
    

    if not rand_int:
        print "path index to remove is:", coord_index
        del paths[paths_index]

    else:
        print "paths index to remove is:", coord_index
        del paths[rand_int]
        
    SHAPE = (1,1)

    #place new_paths in original path set!  
    #paths = np.insert(paths, [1], [new_paths], axis=0)
    
    paths = np.append(paths, [new_paths], axis=0)

    #paths = paths.append(new_paths)
    #paths = np.concatenate((paths, [new_paths]), axis=0)
    
    #paths = np.append(paths, new_paths, axis=0)
    
    #create a flattened numpy array of size 2xN from the paths created! 

    paths_density_check = list(itertools.chain(*paths))
    paths_density_check = np.asarray(list(itertools.chain(*paths_density_check)))
    
    #keep all but the repeated coordinates by keeping only unique whole rows!
    #method is slowed without the b contiguous array
    b = np.ascontiguousarray(paths_density_check).view(np.dtype\
    ((np.void, paths_density_check.dtype.itemsize * \
    paths_density_check.shape[1])))
    _, idx = np.unique(b, return_index=True)
        
    paths_density_check = np.unique(b).view(paths_density_check.dtype)\
    .reshape(-1, paths_density_check.shape[1])

    #plt.figure()
    #plt.scatter(paths_density_check[:,0],paths_density_check[:,1])
    #plt.savefig('{}.png'.format(counter))
    
    #remove 3rd and 4th columns
    #paths_density_check = np.column_stack((paths_density_check[:,0],
    #                                   paths_density_check[:,1]))
           
                 
    #remove all path points that lay outside the shape file polygon
    #paths_density_check = ps.paths_in_shape(paths_density_check)
    paths = list(paths)

    #-------------------------------------------------------------------------
    #-------------------------------------------------------------------------
    #remove zeroes from path_density_check to ensure all paths that
    #were NOT in the distance threshold are removed from the path density
    #calculation!
    path_density_lons,path_density_lats = paths_density_check[:,0], \
                                          paths_density_check[:,1]
    
    FIND_ZERO1 = np.where(paths_density_check[:,0]==0)[0]
    FIND_ZERO2 = np.where(paths_density_check[:,1]==0)[0]
    if len(FIND_ZERO1) != 0 and len(FIND_ZERO2) != 0:
        path_density_lons = np.delete(path_density_lons, FIND_ZERO1)
        path_density_lats = np.delete(path_density_lats, FIND_ZERO2)
    #-------------------------------------------------------------------------
    #-------------------------------------------------------------------------

    # Estimate the 2D histogram
    H, xedges, yedges = np.histogram2d(path_density_lons,
                                       path_density_lats,
                                       bins=nbins)
    #edges_new = ps.paths_in_shape(np.column_stack((xedges,yedges)))
    
    GRAD = np.abs(np.asarray(np.gradient(H)[0]))
    
    # H needs to be rotated and flipped
    H = np.rot90(H)
    GRAD = np.rot90(GRAD)
    
    H = np.flipud(H)
    GRAD = np.flipud(GRAD)
    
    # Mask zeros
    H = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero  
    GRAD = np.ma.masked_where(GRAD==0,GRAD) # Mask pixels with a value of zero  

    H_avg1 = np.average(H)
    grad_avg = np.average(GRAD)
    grad_check1 = np.std(GRAD)
        #find the set of coordinates in the original coordinate set that are the 
    #closest to the highest density anomalies. Select these to be moved!
    search = np.where(H > 0.50 * np.max(H))
    #scale these points with respect to the lat-lon limits!
    Hmaxx, Hmaxy =  search[1], search[0]    
    Hmaxx = (lonmax-lonmin)/(nbins) * Hmaxx + lonmin
    Hmaxy = (latmax-latmin)/(nbins) * Hmaxy + latmin

    #make sure all low density coordinates ARE within shapefile!
    highest_density_coords = ps.paths_in_shape(np.column_stack((Hmaxx, Hmaxy)))
    
    rand_indicator = random.randint(1,10)
    
    if 0 < rand_indicator <= 9: 
    #half the time move the coordinates to low density locations.
    
        WHERE = np.where(H < perc_high*H_avg1)
        #scale these points with respect to the lat-lon limits!
    
        Hminx, Hminy =  WHERE[1], WHERE[0]    
        Hminx = (lonmax-lonmin)/(nbins) * Hminx + lonmin
        Hminy = (latmax-latmin)/(nbins) * Hminy + latmin

        #make sure all low density coordinates ARE within shapefile!
        
        low_density_coords = ps.paths_in_shape(np.column_stack((Hminx, Hminy)))
        
        #low_density_coords = many_points(shape, np.column_stack((Hminx, Hminy)))
        
        if len(low_density_coords) == 0:
            new_coord = ps.points_in_shape(shape_path, 1)[0]
            #increase percentage of search if no new low density points are created!
            perc_high +=0.05
        
        elif len(low_density_coords) == 1:
            new_coord = low_density_coords[0]
            perc_high +=0.05

        else:
            new_coord = low_density_coords[random.randint(0,len(low_density_coords)-1)]

    elif 5 < rand_indicator <= 10: 
    #half the time move coordinates to random locations.
        new_coord = ps.points_in_shape(shape_path, 1)[0]

    if counter == 0:
        grad_ideal = 1e6
        avg_ideal = 0


    
    fig1 = plt.figure()
    #plt.scatter(highest_density_coords[:,0],\
    #highest_density_coords[:,1],c='orange', s=20)
    plt.scatter(coords[:,0], coords[:,1], s=5)
    plt.scatter(new_coord[0], new_coord[1],c='r', s=30)
    plt.scatter(old_coord[0], old_coord[1],c='y', s=30)
    
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    levels = (np.average(H) + 0.15*(np.max(H)-np.average(H)),0)
    cset = plt.contour(H, levels, origin='lower',colors=['red','black'],\
    linewidths=(0.5, 0.5),extent=extent)

    #convert high density contour lines into a set of shapely polygons!


    
    if check_coord != None:
        plt.scatter(check_coord[0], check_coord[1],c='g', s=40)
    
    fig1.savefig('coords.png')
    fig1.clf()
    
    
    


    if grad_check1 < grad_ideal:# and avg_ideal < H_avg1:     
        #counter >= 1 and
        #dump the coordinates!
    
        #print grad_check1, grad_ideal
        #print avg_ideal, H_avg1

        with open(u'ideal_coordinates.pickle', 'wb') as f:
            print "\nExporting new ideal coordinates."
            pickle.dump(coords, f, protocol=2)
        
        grad_ideal = grad_check1
        avg_ideal = H_avg1
        # find indices of pixels where H==HMIN
        #HMATMIN = np.ma.masked_where(H>HMIN,H)
        #only select coordinates where the density is 10% of the average or below!
        fig2 = plt.figure()
        plt.pcolormesh(xedges,yedges,H)
        plt.xlabel('longitude (degrees)')
        plt.ylabel('latitude (degrees)')
        plt.xlim(lonmin-0.05*abs(lonmax-lonmin), lonmax+0.05*abs(lonmax-lonmin))
        plt.ylim(latmin-0.05*abs(latmax-latmin),latmax+0.05*abs(latmax-latmin))
        col = plt.colorbar()
        col.ax.set_ylabel('Counts')
            
        plt.scatter(highest_density_coords[:,0],\
        highest_density_coords[:,1],c='orange', s=10)
        plt.scatter(new_coord[0], new_coord[1],c='r', s=30)
        plt.scatter(old_coord[0], old_coord[1],c='g', s=30)
        fig2.savefig("min_density{}.png".format(counter))
        fig2.clf()

        #print(len(paths))
        #print(len(KEEP_PATHS))
            
    else:
        #RESET!
        #remove new coordinate and replace with old coordinate
        coords = list(coords)
        del coords[-1]
        coords = np.asarray(coords)
        #place new coordinate in old set of coordinates
        coords = np.append(coords, [old_coord], axis=0)
        
        #remove new path and replace it with the old set!
        paths = list(paths)
        del paths[-1]
        paths = list(np.append(paths, [old_path], axis=0))
        
    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
    #calculate the combinations of coordinates from real points, to high points! 
    high_point_coord_combs = [np.vstack([[coord1[0],coord1[1],coord2[0],\
    coord2[1]] for coord1 in highest_density_coords for coord2 in coords])]
    
    #print "high point coordinate combitions:", high_point_coord_combs
    pool = mp.Pool()    
    high_point_coords = pool.map(fast_dists, high_point_coord_combs)
    pool.close()
    pool.join()
    
    high_point_coords =  np.vstack(np.asarray(\
                                   list(itertools.chain(*high_point_coords))))
    
    
#    t0 = datetime.datetime.now()
#    high_point_coords1 = np.sort(high_point_coords.view('i8,i8,i8'),\
#                                 order=['f1'], axis=0).view(np.int)
#    t1 = datetime.datetime.now()
#    print "view sorting time: ", t1-t0
    
#    t0 = datetime.datetime.now()
#    high_point_coords1 = high_point_coords[np.lexsort(np.fliplr(high_point_coords).T)]
#    t1 = datetime.datetime.now()
#    print "flip sorting time: ", t1-t0
    
    #sort values of distance from coordinates to high density points
#    t0 = datetime.datetime.now()
    high_point_coords = high_point_coords[high_point_coords[:,0].argsort()]
#   t1 = datetime.datetime.now()
#    print "fancy index sorting time: ", t1-t0
    
    #SELECT ONLY 20% OF THESE LOWEST DISTANCE COORDINATES FROM HIGH DENSITY POINTS! 
    point_index = int(0.2*len(high_point_coords))
    high_point_coords = high_point_coords[:point_index]
    
    #create unique list of coordinates from the above high_point_coords and remove distances
    high_point_coords =  np.column_stack((high_point_coords[:,1], 
                                          high_point_coords[:,2]))
    b = np.ascontiguousarray(high_point_coords).view(np.dtype\
    ((np.void, high_point_coords.dtype.itemsize * \
    high_point_coords.shape[1])))
    _, idx = np.unique(b, return_index=True)
        
    high_point_coords = np.unique(b).view(high_point_coords.dtype)\
    .reshape(-1, high_point_coords.shape[1])
        
    #SET THIS AS THE OLD COORDINATE TO REMOVE!
    check_coord = high_point_coords[random.randint(0,len(high_point_coords))]
    
    if check_coord in coords:
        old_coord = check_coord
    else:
        check_coord = None
    
    #cd Desktop/Link\ to\ SIMULATIONS/Network_Tracks/smarter_model/
    
    #-------------------------------------------------------------------------
    #RANDOMLY SELECT COORD FROM HIGH DENSITY COORDS
    #check_coord = high_point_coords[random.randint(0, len(high_point_coords))]
    #if check_coord in coords:
    #    new_coord = check_coord
    #print new_coord
    #-------------------------------------------------------------------------

    #select high density coordinate!
    #point_to_check = high_point_coords[0][1:]
    #if len(high_point_coords) > 0 \
    #and -180 <= point_to_check[0] <= 180\
    #and -90 <= point_to_check[1] <= 90:
    #    new_coord = point_to_check
    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------


    find_it = []
    counter+=1
    counter2+=1

    t1 = datetime.datetime.now()
    print t1-t0    