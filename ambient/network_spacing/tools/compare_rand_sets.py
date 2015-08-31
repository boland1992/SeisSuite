# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 18:59:35 2015

@author: boland

CODE DESCRIPTION: 
The following python script is used to test how to geospatially compare two 
2D arrays.
"""

from scipy import spatial
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp

n = 10000
minx,maxx,miny,maxy = -180,180,-90,90

#tree = spatial.KDTree(zip(x.ravel(), y.ravel())) 
#tree.data

def rand_coords(n):
    X = abs(maxx - minx) * np.random.rand(n,1) + minx
    Y = abs(maxy - miny) * np.random.rand(n,1) + miny
    return np.column_stack((X.ravel(),Y.ravel()))
    
static_coords = rand_coords(n)
first_set = rand_coords(n)
SUM = np.sum(np.abs((np.abs(static_coords)-np.abs(first_set))))
SUM =101.
while SUM < 100:
    new_set = rand_coords(n)
    SUM = np.sum(np.abs((np.abs(static_coords)-np.abs(new_set))))

TREE = spatial.KDTree(rand_coords(n))

#find nearest neighbour distances for each of the given coordinates pairs
for i in range(0,1000):
    QUERY = TREE.query(rand_coords(n))[0]
    comp_stat = np.sum(QUERY); print comp_stat

