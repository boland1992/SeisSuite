# -*- coding: utf-8 -*-
"""
Created on Thu May  7 11:35:27 2015

@author: boland
"""
import numpy as np
import itertools
import datetime
import os

t0 = datetime.datetime.now()
periods1 = np.arange(4,31,1)
periods2 = np.arange(5,32,1)

periods = np.column_stack((periods1,periods2))


from pysismo import pscrosscorr

#PICKLE_PATH = '/storage/ANT/PROGRAMS/ANT_OUTPUT/OUTPUT/CROSS/06.05.2015-15:53:28/XCORR-STACK_01.01.2014-31.12.2014_datalesspaz.pickle'
PICKLE_PATH = '/storage/ANT/OUTPUT/CROSS/xcorr_2014-2014_datalesspaz_AU-S-2014-FULL.pickle'

xc = pscrosscorr.load_pickled_xcorr(PICKLE_PATH)

# optimizing time-scale: max time = max distance / vmin (vmin = 2.5 km/s)
maxdist = max([xc[s1][s2].dist() for s1, s2 in xc.pairs()])
maxt = min(1500., maxdist/2.5)


output_folder = '/storage/ANT/spectral_density/period_distance_plots'

for period in periods:     
    
    file_name = 'distance_plot_periodband_{}-{}s.png'.format(period[0], period[1])
    absfile_name = os.path.join(output_folder, file_name)
    
    
    #plot distance plot of cross-correlations with individual period combos
    xc.plot(plot_type='period_distance', xlim=(-maxt, maxt), 
    outfile=absfile_name, showplot=False,periodband=period)

    print 'period-band {}-{}s complete'.format(period[0], period[1])





t1 = datetime.datetime.now()
print "total time to compute all possible {} combinations for period distance \
plot is: ".format(len(periods)), t1-t0