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
from seissuite.spacing.search_station import (InShape, 
                                              InPoly, 
                                              Geodesic, 
                                              Coordinates, 
                                              Density)

import os
import pickle
import pyproj
import datetime
import numpy as np
import datetime as dt
#import pointshape as ps
import multiprocessing as mp
from scipy.cluster.vq import kmeans
from seissuite.misc.dataless import Dataless
import matplotlib.pyplot as plt
import itertools as it
from matplotlib.colors import LogNorm





#------------------------------------------------------------------------------
# VARIABLES
#------------------------------------------------------------------------------

verbose = True
#Enter path to boundary shape file.
shape_boundary = True

dataless = False


# Enter number new stations desired.
n_stations = 6
# Enter km spacing between path density points.
km_points = 0.5
# Reference elipsoid to calculate distance.
wgs84 = pyproj.Geod(ellps='WGS84')
# Enter number of bins for 2D Histogram density calculation. 
nbins = 75
# Enter estimated average shear wave velocity. 3kms-1 is the default!
velocity = 3.0
# Define your ambient noise period range OR individual period in seconds.
global period_range
period_range = [1,40]

# Set box length size in metres (square field)
grid_length = 2000 # (metres)

box_length = grid_length / nbins

# Enter path to dataless file
dataless_path = 'ALL_AUSTRALIA.870093.dataless'

#coords = Dataless.locs_from_dataless(dataless_path)

shape_path = '/home/iese/10kmR.shp' 


t0 = dt.datetime.now()

# Generate InShape class
SHAPE = InShape(shape_path)
# Create shapely polygon from imported shapefile 
UNIQUE_SHAPE = SHAPE.shape_poly()
#print type(UNIQUE_SHAPE)
# Generate InPoly class
INPOLY = InPoly(shape_path)
# Create matplotlib Path object from imported shapefile
#outer_shape = UNIQUE_SHAPE.buffer(1.,resolution=1)
#inner_shape = UNIQUE_SHAPE.buffer(-8,resolution=1)

#outer_poly = INPOLY.poly_from_shape(shape=outer_shape)
#inner_poly = INPOLY.poly_from_shape(shape=inner_shape)

poly = INPOLY.poly_from_shape(shape=UNIQUE_SHAPE)

many_points = INPOLY.rand_poly(poly=poly, N=1e5)

# Scale smaller shape to fit inside larger shape. 
#SMALL_SHAPE = scale(UNIQUE_SHAPE, xfact=0.3, yfact=0.3)
#points_in_small_shape = INPOLY.rand_shape(shape=SMALL_SHAPE, IN=False)
# Generate matplotlib Path object for the small scalled polygon 
#small_poly = INPOLY.node_poly(SHAPE.external_coords(shape=SMALL_SHAPE))
# Remove points that are outside the buffered_poly
#outer_poly_points = INPOLY.points_in(many_points, poly=outer_poly)

# Remove points that are inside the small_poly
poly_points = np.asarray(INPOLY.points_in(many_points, 
                                                poly=poly,
                                                IN=True))

cluster_points = np.asarray(kmeans(poly_points, n_stations)[0])


#plt.figure()
#plt.scatter(poly_points[:,0], poly_points[:,1], c='b')
#plt.scatter(cluster_points[:,0], cluster_points[:,1], c='orange', s=35)
#plt.show()

#-----------------------------------------------------------------------------
# GENERATE SECOND SET OF VARIABLES AND STATES
#-----------------------------------------------------------------------------
ideal_path = 'ideal_coordinates.pickle'
#if no paths have been done before, start afresh!
#if dataless:
#    coords = Dataless.locs_from_dataless(dataless_path)
#    original_coords = coords
#elif os.path.exists(ideal_path):
#    f = open(name=ideal_path, mode='rb')
#    coords = pickle.load(f)
#    f.close()

GEOD = Geodesic()


coords = cluster_points

lons, lats = coords[:,0], coords[:,1]

dists = []
for lon1, lat1 in zip(lons, lats):
    for lon2, lat2 in zip(lons, lats):
       dists.append(GEOD.haversine(lon1, lat1, lon2, lat2))
       
       
for i in dists: 
    print i        
print "Min. Interstation Distance: ", np.min(dists)
print "Avg. Interstation Distance: ", np.average(dists)
print "Max. Interstation Distance: ", np.max(dists)


plt.figure()
plt.scatter(lons, lats)
plt.show()

# dump statistically optimised station spacings to .csv file
np.savetxt('/home/iese/Desktop/optimised_spacings.csv', coords, delimiter=",")




#-----------------------------------------------------------------------------
# POINT DENSITY ESTIMATION
#-----------------------------------------------------------------------------


lonmin, lonmax = np.floor(min(coords[:,0])), np.ceil(max(coords[:,0]))
latmin, latmax = np.floor(min(coords[:,1])), np.ceil(max(coords[:,1]))
print lonmin, lonmax, latmin, latmax

kappa = [np.vstack([[coord1[0],coord1[1],coord2[0],coord2[1]]\
                    for coord2 in coords]) for coord1 in coords]
             
             
GEODESIC = Geodesic(km_point=0.01)

def fast_geodesic(lon1, lat1, lon2, lat2, npts):
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
    return np.array(path)


npts = 100

total_points = []

for i in kappa:
    paths = GEODESIC.fast_paths(i)

    for j in i: 

        lon1, lat1, lon2, lat2 = j
        path = fast_geodesic(lon1, lat1, lon2, lat2, npts)
        
    #    GEODESIC.fast_paths(path)
        
        total_points.append(path)
        
total_points = list(it.chain(*total_points))
total_points = np.array(total_points)

total_points = np.asarray(INPOLY.points_in(total_points, 
                                           poly=poly,
                                           IN=True))


plt.figure()
plt.scatter(total_points[:,0], total_points[:,1])
plt.scatter(coords[:,0], coords[:,1], c='orange')

plt.show()


DENSITY = Density(paths=total_points, nbins=nbins)

H, xedges, yedges = DENSITY.hist2d(paths=total_points)

#histogram_GIS = np.column_stack((H, xedges, yedges))


print H.shape, xedges.shape, yedges.shape


coords = np.array([[x, y] for x in xedges[:-1] for y in yedges[:-1]])


H = np.rot90(H)
H = np.flipud(H)
#H = np.rot90(H)
#H = np.rot90(H)



#plt.figure()
#plt.scatter(coords[:,0], coords[:,1])
#plt.scatter(coords[0][0], coords[0][1], c='orange')
#plt.scatter(coords[-1][0], coords[-1][1], c='red')
#plt.show()
#plt.clf()


#grad = DENSITY.hgrad(H=H)
    
#H_avg1 = np.average(H)
#grad_check1 = np.std(grad)
    
H_masked = DENSITY.transform_h(H=H)
#grad = DENSITY.transform_grad(grad=grad)
    
DENSITY.plot_field(SHAPE=UNIQUE_SHAPE)


#-----------------------------------------------------------------------------
# RESOLUTION ESTIMATION
#-----------------------------------------------------------------------------



# define the point-density to metre resolution conversion function
def res_conv(point_density, length=box_length, k2=100.0):
    """
    Function to convert a point in a point density matrix into an approx. 
    metre resolution value. 
    """
    
    if point_density < k2 / 20:
        return 0
        
    else: 
        k = length
        x = point_density
        
        alpha = k2 * np.log(k - 1)     
    
        return np.exp(alpha / x) + k - 1 

    

H_unpacked = np.array(list(it.chain(*H)))

H_resolution = np.array(map(res_conv, H_unpacked))



GIS_output = np.column_stack((coords[:,0],coords[:,1], 
                              H_unpacked, H_resolution))
         

H_res = H_resolution.reshape(H.shape)


plt.pcolor(xedges, yedges, H_res, norm=LogNorm(\
    vmin=np.min(H_masked), vmax=np.max(H_masked)), cmap='rainbow',\
    alpha=0.6, zorder = 3)


col = plt.colorbar()
col.ax.set_ylabel('Estimated Resolution (m)')

plt.show()


np.savetxt("/home/iese/Documents/Ben/OLSUS_GIS/resolution_matrix.csv", 
           GIS_output, delimiter=",")


print GIS_output

quit()


print kappa

def spread_paths(coord_list):
    return GEODESIC.fast_paths(coord_list)

paths = map(spread_paths, kappa)




print paths
quit()

t0 = datetime.datetime.now()
pool = mp.Pool()    
paths = pool.map(spread_paths, kappa)
pool.close()
pool.join()

t1 = datetime.datetime.now()
print t1-t0


print paths
#create a flattened numpy array of size 2xN from the paths created! 
paths1 = GEODESIC.combine_paths(paths)

paths = list(paths)

paths1 = GEODESIC.remove_zeros(paths1)

plt.figure()
plt.scatter(paths1[:,0], paths1[:,1])
plt.show()


DENSITY = Density(paths=paths1)

H, xedges, yedges = DENSITY.hist2d(paths=paths1)
grad = DENSITY.hgrad(H=H)
    
H_avg1 = np.average(H)
grad_check1 = np.std(grad)
    
H_masked = DENSITY.transform_h(H=H)
grad = DENSITY.transform_grad(grad=grad)
    
DENSITY.plot_field(SHAPE=UNIQUE_SHAPE)


quit()
COORDS = Coordinates()



for i in [0]:
    t0 = datetime.datetime.now()
    
    #----------------------------------------------------------------------
    # Generate N new point coordinates
    #----------------------------------------------------------------------
    #new_coords = N_cluster_points
   
#    if cluster:
 #       new_coords = N_cluster_points
#    else:
    #    new_coords = ps.points_in_shape(shape_path, N)
        
    #coords = np.append(coords, new_coords, axis=0)

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


#    plt.figure()
#    plt.plot(paths1[:,0], paths1[:,1])
#    plt.show()


    DENSITY = Density(paths=paths1)

    H, xedges, yedges = DENSITY.hist2d(paths=paths1)
    grad = DENSITY.hgrad(H=H)
    
    H_avg1 = np.average(H)
    grad_check1 = np.std(grad)
    
    H_masked = DENSITY.transform_h(H=H)
    grad = DENSITY.transform_grad(grad=grad)
    
    DENSITY.plot_field(SHAPE=UNIQUE_SHAPE)

    quit()    
    #search = np.where(H<0.1*np.average(H))
    #Hmaxx, Hmaxy =  search[1], search[0]    
    #Hmaxx = (lonmax-lonmin)/(nbins) * Hmaxx + lonmin
    #Hmaxy = (latmax-latmin)/(nbins) * Hmaxy + latmin
    # Make sure all low density coordinates ARE within shapefile!
    #low_density_coords = ps.paths_in_shape(np.column_stack((Hmaxx, Hmaxy)))
    
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
        
        DENSITY.plot_field(SHAPE=UNIQUE_SHAPE)

        grad_ideal = grad_check1
        avg_ideal = H_avg1

    coords = COORDS.del_N(N=n_stations, inputs=coords)
    paths = COORDS.del_N(N=n_stations, inputs=paths)
    paths=list(paths)

    counter+=1
    t1 = datetime.datetime.now()
    print "That loop took: ", t1-t0
