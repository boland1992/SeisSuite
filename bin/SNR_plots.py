# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 23:34:00 2015

@author: boland
"""

import pickle
import os
# import CONFIG class initalised in ./configs/tmp_config.pickle
pickle_file = '/storage/MASTERS/CONFIGURATIONS/TEST/OUTPUT/CROSS/16.10.2015-22:48:34/XCORR-STACK_01.01.2014-31.01.2014_datalesspaz.pickle'

config_file = '/home/boland/seissuite/bin/configs/TEST.cnf'

f = open(name=pickle_file, mode='rb')
xc = pickle.load(f)
f.close()
OUTFILE = os.getcwd()
#xc.plot_SNR(plot_type='all', outfile=OUTFILE, 
#                    config=os.path.basename(config_file))

print xc.pairs()