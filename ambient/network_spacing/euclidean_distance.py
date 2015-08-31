# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 09:07:28 2015

@author: boland

CODE DESCRIPTION: 
The following python script is used to find the euclidean distances between
(2xN) vector of 2D points. 
"""

import numpy as np
import matplotlib.pyplot as plt
import datetime
from scipy.spatial.distance import pdist

t0 = datetime.datetime.now()
 
 
 
n = 50
counts = 0.
no = 1e6

variance_lowest = 10


while counts < no: 
 
 #   x = np.random.rand(n) #* 10
 #   x = [[i] for i in x]
    
 #   y = np.random.rand(n)
 #   y = [[i] for i in y]
    
    #generate grid of random points within a 1x1 box    
    X = np.random.rand(n,2)

    #find all distances between all points within the 1x1 box        
    distances = pdist(X, 'euclidean')
    
    point_distances = np.split(distances, n-1)
    
    #calculate the mean of the nearest 4 point distances for each point.
    point_distances = [np.mean(i[:3]) for i in point_distances]
    
    #minimise variation between these means to find best station spacing!
    
    #calculate the variance of the distances between all points. Lowest variation wins!
    #standard = np.std(distances)
    
    variance = np.var(point_distances)

    
    if variance < variance_lowest:
        
        ideal_distribution = X
        variance_lowest = variance
    
        
    counts += 1


plt.figure(1)
plt.scatter(ideal_distribution[:,0],ideal_distribution[:,1])
plt.xlim(0,1.)
plt.ylim(0,1.) 
plt.show()


t1 = datetime.datetime.now()


print("\nThe time taken to perform {} tasks was: {}".format(no, t1-t0))



#a = np.matrix([ [1,x[0][0]], [1,x[1][0]], [1,x[2][0]] ])

#print(a)

#b = np.matrix([ [x[0][1]], [x[1][1]], [x[2][1]] ])

#print(b)


#yy = (a.T * a).I * a.T * b

#xx = np.linspace(1, 10, 50)

#y = np.array(yy[0] + yy[1] * xx)





#plt.figure(1)
#plt.plot(xx, y.T, color='r')
#plt.scatter([x[0][0], x[1][0], x[2][0] ], [x[0][1], x[1][1], x[2][1] ]) 
#plt.show()