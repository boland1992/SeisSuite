# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 10:54:58 2015

@author: boland
"""

import itertools as it
import numpy as np

# the following options will remain the same throughout all processing configs
TDD = True
RESP_REMOVE = True
DOWNSAMPLE = True
COMPLETENESS = True
TIME_NOMALISATION = True
SPEC_WHITENING = True

# the following will vary depending on the combination of config options

HIGHAMP_REMOVE = True
RESP_CHECK = True
BANDPASS = True
                 
config_options = ['HIGHAMP_REMOVE', 'RESP_CHECK', 'BANDPASS']

zeros = np.zeros_like(config_options)

options_zeros = list(it.chain(*[config_options, zeros]))


options = list(it.combinations(options_zeros, len(config_options)))

total_options = []

for opt in options:
    if opt not in total_options:
        total_options.append(opt)


# create a dictionary of the options available
for opt in total_options:
    print opt
    
    
print "Number of processing configurations: ", len(total_options)