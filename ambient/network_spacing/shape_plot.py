# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 14:53:27 2015

@author: boland

CODE DESCRIPTION:
The following python script is used to demonstrate how to import a shapefile
with the extension .shp into python using the fiona module, and then plot
that file using matplotlib.pyplot.plot()
"""



#http://stackoverflow.com/questions/2770356/\
#extract-points-within-a-shape-from-a-raster

import fiona
import shapefile
import random
import matplotlib.pyplot as plt
import numpy as np


# Source shapefile - can be any polygon
shpFilePath = "/home/boland/Dropbox/University/UniMelb/AGOS/PROGRAMS\
/ANT/Versions/26.04.2015/shapefiles/aus.shp"
input_shape = '/home/boland/Downloads/NZ Shapefiles/clipped_NZ.shp'

listx=[]
listy=[]

r = shapefile.Reader(input_shape)

for sr in r.shapeRecords():
    for xNew,yNew in sr.shape.points:
        listx.append(xNew)
        listy.append(yNew)
        
plt.plot(listx,listy)
plt.show()
