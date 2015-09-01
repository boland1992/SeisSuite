#!/usr/bin/env python

# combining density estimation and delaunay interpolation for confidence-weighted value mapping
# Dan Stowell, April 2013

import numpy as np
from numpy import random
#from math import exp, log

from scipy import stats, mgrid, c_, reshape, rot90
import matplotlib.delaunay
#import matplotlib.tri as tri
import matplotlib.delaunay.interpolate
from matplotlib.colors import LogNorm

import matplotlib.pyplot as plt
import matplotlib.cm as cm
#from colorsys import hls_to_rgb
import pickle

pickle_file =  '/storage/ANT/spectral_density/station_pds_maxima/\
S Network 2014/noise_info0_SNetwork2014.pickle'

f = open(name=pickle_file, mode='rb')
data = pickle.load(f)
f.close()


#############################
# user settings
n = 100
gridsize = 100
fontsize = 'xx-small'

#############################
# first generate some random [x,y,z] data -- random locations but closest to the middle, and random z-values
# we will add some correlation to the z-values
data[:,2] += data[:,1]
data[:,2] += data[:,0]
# scale the z-values to 0--1 for convenience
zmin = np.min(data[:,2])
zmax = np.max(data[:,2])

xmin = np.min(data[:,0])
xmax = np.max(data[:,0])
ymin = np.min(data[:,1])
ymax = np.max(data[:,1])
zmin = np.min(data[:,2])
zmax = np.max(data[:,2])



##################################################
# plot it simply
plt.figure()
##################################################
# now make a KDE of it and plot that

kdeX, kdeY = mgrid[xmin:xmax:gridsize*1j, ymin:ymax:gridsize*1j]
positions = c_[kdeX.ravel(), kdeY.ravel()]

values = c_[data[:,0], data[:,1]]
kernel = stats.kde.gaussian_kde(values.T)
kdeZ = reshape(kernel(positions.T).T, kdeX.T.shape)


##################################################
# now make a delaunay triangulation of it and plot that

tt = matplotlib.delaunay.triangulate.Triangulation(data[:,0], data[:,1])
print xmin, xmax, ymin, ymax
print gridsize
extrap = tt.nn_extrapolator(data[:,2])
interped = extrap[xmin:xmax:gridsize*1j, ymin:ymax:gridsize*1j]

##################################################
# now combine delaunay with KDE

colours = np.zeros((gridsize, gridsize, 4))
kdeZmin = np.min(kdeZ)
kdeZmax = np.max(kdeZ)
confdepth = 0.45
for x in range(gridsize):
	for y in range(gridsize):
		conf = (kdeZ[x,y] - kdeZmin) / (kdeZmax - kdeZmin)
		val  = min(1., max(0., interped[x,y]))
		colour = list(cm.rainbow(val))
		# now fade it out to white according to conf
		for index in [0,1,2]:
			colour[index] = (colour[index] * conf) + (1.0 * (1. -conf))
		colours[x,y,:] = colour
		#colours[x,y,:] = np.hstack((hls_to_rgb(val, 0.5 + confdepth - (confdepth * conf), 1.0), 1.0))
		#colours[x,y,:] = [conf, conf, 1.0-conf, val]
  
  
print colours
plt.imshow(rot90(colours), cmap=cm.rainbow, norm=LogNorm(\
            vmin=zmin, vmax=zmax))
plt.title("interpolated & confidence-shaded")

plt.ylim([ymin,ymax])
plt.xlim([xmin,xmax])
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)

############################################
plt.savefig("plot_heati_simple.svg", format='SVG')
