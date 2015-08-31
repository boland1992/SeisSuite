# -*- coding: utf-8 -*-
"""
Created on Wed May  6 12:39:14 2015

@author: boland
"""

import matplotlib.pyplot as plt 

def shift1(xcorr_len):
    """
    Function that produces an estimatated xcorr shift_length maximum given
    the length of the traces being cross-correlated. 
    """    
    
    a = 2047000.0/1023.0
    b = 2488320000000.0/341.0
    
    shift_len = (a*xcorr_len**2 - b) / (xcorr_len**2)
    return shift_len
    
    
    
def shift2(xcorr_len):
    
    a = 10.0 / 837.0
    b = 3.0e4 / 31.0
    
    shift_len = a * xcorr_len + b
    
    
    return shift_len




y = []
counts = 0
for i in range(1, 86400 + 1):
    y.append(shift2(i))
    if shift2(i) < 0:
        counts += 1

print(counts)

plt.figure(1)

plt.plot(y)
plt.show()

