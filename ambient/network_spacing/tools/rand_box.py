# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 15:39:27 2015

@author: boland

CODE DESCRIPTION: 
The following python script is used to generate n random points within a set
bounding box and plot thenm using matplotlib.pyplot.scatter()
"""

import numpy as np
from scipy.cluster.vq import kmeans
import matplotlib.pyplot as plt

n = 100000
counts = 0.
no = 1e6
opac = 1.
variance_lowest = 10
counts = 0
plt.figure()

while counts < 50:
 #   x = np.random.rand(n) #* 10
 #   x = [[i] for i in x]
    
 #   y = np.random.rand(n)
 #   y = [[i] for i in y]
    
    #generate grid of random points within a 1x1 box    
    X = np.random.rand(n,2)

    k = kmeans(X, 20)
    print type(X)
    counts+=1



    #plt.scatter(X[:,0],X[:,1])
    plt.scatter(k[0][:,0], k[0][:,1], alpha = opac)
    plt.xlim(0,1.)
    plt.ylim(0,1.) 
    
    opac *= 0.9

plt.show()

#plt.figure()
#
#plt.xlim(0,1.)
#plt.ylim(0,1.) 
#plt.show()

