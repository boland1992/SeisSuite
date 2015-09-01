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
from math import sqrt, atan2, radians,degrees, cos, tan, sin, asin
import random

shape_path = "/home/boland/Dropbox/University/UniMelb/AGOS/PROGRAMS/ANT/\
Versions/26.04.2015/shapefiles/aus.shp"

N = 130
#enter km spacing between path density points
km_points = 20.0

# reference elipsoid to calculate distance
wgs84 = pyproj.Geod(ellps='WGS84')
nbins = 200

#enter estimated average shear wave velocity. 3kms-1 is the default!
velocity = 3.0
#define your ambient noise period range OR individual period in seconds
global period_range
period_range = 40

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
    
    #print("still going strong\n")
    length = len(path)
    lons = [lon1 for i in range(0,length)]
    lats = [lat1 for i in range(0,length)]
    
    path = np.column_stack((path,lons,lats))
    if min(dist_range) < dist < max(dist_range):

        return path
    else:
        return np.zeros_like(path)

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

def spread_paths(nets):

    #pool = mp.Pool()    
    #paths = pool.map(new_paths, nets)
    #pool.close()
    #pool.join()
    paths = map(paths_func, nets)
    #create a flattened numpy array of size 2xN from the paths created! 
    #paths = np.asarray(list(itertools.chain(*paths)))

    #keep all but the repeated coordinates by keeping only unique whole rows!
    #method is slowed without the b contiguous array
    #b = np.ascontiguousarray(paths).view(np.dtype((np.void, paths.dtype.itemsize * paths.shape[1])))
    #_, idx = np.unique(b, return_index=True)
    
    #paths = np.unique(b).view(paths.dtype).reshape(-1, paths.shape[1])

    #plt.figure()
    #plt.scatter(paths[:,0],paths[:,1])
    #name =  uuid.uuid4()
    #plt.savefig('{}.png'.format(name))
    
    return paths
    
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
while infinite_counter <= 1:
    
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
        #print(counter2)
        #----------------------------------------------------------------------
        
        #or random selection of paths?!
        #----------------------------------------------------------------------
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
        
        itemindex = find_it[0]
        print itemindex
        coords = list(coords)
        #find index of array in nested array to remove!
        del coords[itemindex]
        coords = np.asarray(coords)
        
        new_coord_first = new_coord

        #----------------------------------------------------------------------
        #generate new point coordinate
        if not counter >= 1:
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
            #remove substitude back the old coordinate for the new coordinate!
            coords = list(coords)
            #find index of array in nested array to remove!
            del coords[-1]
            coords = np.asarray(coords)
            #place new coordinate in old set of coordinates
            coords = np.append(coords, [old_coord], axis=0)
    
            
        #print "new paths shape:", SHAPE
        
        #paths = np.asarray(paths)
        #if np.asarray(new_paths).shape != (130,):
        #    print("This one's trouble")
        #    print np.asarray(new_paths).shape
        #    new_paths =  np.asarray(new_paths[0]).reshape(130,)
    
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
    
    grad_check1 = np.std(GRAD)
    
    rand_indicator = random.randint(1,10)
    
    if 0 < rand_indicator <= 5: 
    #half the time move the coordinates to low density locations.
    
        WHERE = np.where(H < perc_high*H_avg1)
        #scale these points with respect to the lat-lon limits!
    
        Hminx, Hminy =  WHERE[1], WHERE[0]    
        Hminx = (lonmax-lonmin)/(nbins) * Hminx + lonmin
        Hminy = (latmax-latmin)/(nbins) * Hminy + latmin

        #make sure all low density coordinates ARE within shapefile!
        low_density_coords = ps.paths_in_shape(np.column_stack((Hminx, Hminy)))
        
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

    plt.figure(4)
    plt.scatter(coords[:,0], coords[:,1])
    plt.scatter(new_coord[0],new_coord[1],c='r')
    plt.scatter(old_coord[0],old_coord[1],c='y')
    plt.savefig('coords.png')
    plt.clf()

    if grad_check1 < grad_ideal and avg_ideal < H_avg1:     
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
        fig = plt.figure()
        plt.pcolormesh(xedges,yedges,H)
        plt.xlabel('longitude (degrees)')
        plt.ylabel('latitude (degrees)')
        plt.xlim(lonmin-0.05*abs(lonmax-lonmin), lonmax+0.05*abs(lonmax-lonmin))
        plt.ylim(latmin-0.05*abs(latmax-latmin),latmax+0.05*abs(latmax-latmin))
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('Counts')
        try:
            plt.scatter(low_density_coords[:,0], 
                        low_density_coords[:,1], 
                        color='red')
        except:
            a = 5 
        fig.savefig("min_density.png".format(counter))
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
        

        # find indices of pixels where H==HMIN
        #HMATMIN = np.ma.masked_where(H>HMIN,H)
        #only select coordinates where the density is 10% of the average or below!
        fig = plt.figure()
        plt.pcolormesh(xedges,yedges,H)
        plt.xlabel('longitude (degrees)')
        plt.ylabel('latitude (degrees)')
        plt.xlim(lonmin-0.05*abs(lonmax-lonmin), lonmax+0.05*abs(lonmax-lonmin))
        plt.ylim(latmin-0.05*abs(latmax-latmin),latmax+0.05*abs(latmax-latmin))
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('Counts')
        try:
            plt.scatter(low_density_coords[:,0], 
                        low_density_coords[:,1], 
                        color='red')
        except:
            a = 5 
        fig.savefig("min_density.png".format(counter))
        
    
    
#plt.scatter(Hminx, Hminy, color='yellow')
    #grad_check2 = grad_check1
    #H_avg2 = H_avg1

    #print(counter)
    find_it = []
    counter+=1
    counter2+=1

    t1 = datetime.datetime.now()
    print t1-t0    