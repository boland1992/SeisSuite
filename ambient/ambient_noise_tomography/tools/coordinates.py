# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 08:53:30 2015

@author: boland
"""

import numpy as np

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