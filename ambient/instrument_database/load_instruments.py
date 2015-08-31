# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 18:08:51 2015

@author: boland
"""

import pickle

pickle_path = 'total_instrument_list.pickle'

f = open(name=pickle_path, mode='rb')
instrument_list = pickle.load(f)
f.close()

for instrument in instrument_list: print instrument
print 'Number of possible instruments present:', len(instrument_list)