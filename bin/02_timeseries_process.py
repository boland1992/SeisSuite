#!/usr/bin/python -u
"""
[Advice: run this script using python with unbuffered output:
`python -u 01_timeseries_process.py`]

This script reads seismic waveform data from a set of stations, and
calculates the cross-correlations between all pairs of stations. The
data (in miniseed format) must be located in folder *MSEED_DIR*. The
stations information (coordinates, instrument response) can be read
from dataless seed files (if *USE_DATALESSPAZ* = True) located in
folder *DATALESS_DIR*, and/or stationXML files (if *USE_STATIONXML* =
True) located in folder *STATIONXML_DIR*. Note that two different
stations MUST HAVE DIFFERENT NAMES, even if they do not belong to
the same network. Also, one given station cannot have several
sets of coordinates: if so, it will be skipped.

In the current version of the program, miniseed files MUST be
organized inside their directory as:
<year>-<month>/<network>.<station>.<channel>.mseed, e.g.:
1988-10/BL.JFOB.BHZ.mseed
So, there is one sub-directory per month, and inside it, one miniseed
file  per month and per station.

The implemented algorithm follows the lines of Bensen et al.,
"Processing seismic ambient noise data to obtain reliable broad-band
surface wave dispersion measurements", Geophys. J. Int. (2007).

The procedure consists in stacking daily cross-correlations between
pairs of stations, from *FIRSTDAY* to *LASTDAY* and, in each given day,
rejecting stations whose data fill is < *MINFILL*. Define a subset of
stations to cross-correlate in *CROSSCORR_STATIONS_SUBSET* (or let it
empty to cross-correlate all stations). Define a list of locations to
skip in *CROSSCORR_SKIPLOCS*, if any. The cross-correlations are
calculated between -/+ *CROSSCORR_TMAX* seconds.

Several pre-processing steps are applied to the daily seismic waveform
data, before the daily cross-correlation is calculated and stacked:

(1) removal of the instrument response, the mean and the trend;

(2) band-pass filter between *PERIODMIN* and *PERIODMAX* sec

(3) down-sampling to sampling step = *PERIOD_RESAMPLE* sec

(4) time-normalization:

    - if *ONEBIT_NORM* = False, normalization of the signal by its
      (smoothed) absolute amplitude in the earthquake period band,
      defined as *PERIODMIN_EARTHQUAKE* - *PERIODMIN_EARTHQUAKE* sec.
      The smoothing window is *PERIODMAX_EARTHQUAKE* / 2;

    - if *ONEBIT_NORM* = False, one-bit normalization, wherein
      only the sign of the signal is kept (+1 or -1);

(5) spectral whitening of the Fourier amplitude spectrum: the Fourier
    amplitude spectrum of the signal is divided by a smoothed version
    of itself. The smoonthing window is *WINDOW_FREQ*.

Note that all the parameters mentioned above are defined in the
configuration file.

When all the cross-correlations are calculated, the script exports
several files in dir *CROSS*

"""

from seissuite.ant import (pscrosscorr,
                           psstation, 
                           pspreprocess, 
                           pserrors, 
                           psstationSQL)
import os
import warnings
import datetime as dt
import itertools as it
import pickle
import obspy.signal.cross_correlation
import time
import glob
import matplotlib.pyplot as plt
#import numpy as np
# turn on multiprocessing to get one merged trace per station?
# to preprocess trace? to stack cross-correlations?
MULTIPROCESSING = {'merge trace': True,
                   'process trace': True,
                   'cross-corr': True}
# how many concurrent processes? (set None to let multiprocessing module decide)
NB_PROCESSES = None
if any(MULTIPROCESSING.values()):
    import multiprocessing as mp

# ====================================================
# parsing configuration file to import some parameters
# ====================================================

from seissuite.ant.psconfig import (MSEED_DIR, DATABASE_DIR, DATALESS_DIR, 
                                    STATIONXML_DIR, CROSSCORR_DIR,
                                    USE_DATALESSPAZ, USE_STATIONXML, 
                                    CROSSCORR_STATIONS_SUBSET, 
                                    CROSSCORR_SKIPLOCS, FIRSTDAY, LASTDAY, 
                                    MINFILL, FREQMIN, FREQMAX, CORNERS, 
                                    ZEROPHASE, PERIOD_RESAMPLE, ONEBIT_NORM, 
                                    FREQMIN_EARTHQUAKE, FREQMAX_EARTHQUAKE, 
                                    WINDOW_TIME, WINDOW_FREQ,
                                    XCORR_INTERVAL, CROSSCORR_TMAX)
                                    
print "\nProcessing parameters:"
print "- dir of miniseed data: " + MSEED_DIR
print "- dir of dataless seed data: " + DATALESS_DIR
print "- dir of stationXML data: " + STATIONXML_DIR
print "- output dir: " + CROSSCORR_DIR
print "- cross-correlation length (mins): " + str(XCORR_INTERVAL)
print "- cross-correlation maximum time interval (s): " + str(CROSSCORR_TMAX)


print "- band-pass: {:.1f}-{:.1f} s".format(1.0 / FREQMAX, 1.0 / FREQMIN)
if ONEBIT_NORM:
    print "- normalization in time-domain: one-bit normalization"
else:
    s = ("- normalisation in time-domain: "
         "running normalisation in earthquake band ({:.1f}-{:.1f} s)")
    print s.format(1.0 / FREQMAX_EARTHQUAKE, 1.0 / FREQMIN_EARTHQUAKE)
fmt = '%d/%m/%Y'
s = "- cross-correlation will be stacked between {}-{}"
print s.format(FIRSTDAY.strftime(fmt), LASTDAY.strftime(fmt))
subset = CROSSCORR_STATIONS_SUBSET
if subset:
    print "  for stations: {}".format(', '.join(subset))
print

    
# Initializing collection of cross-correlations
xc = pscrosscorr.CrossCorrelationCollection()

#create a metadata list, may need dictionary based on how much info required         
metadata = [] 

#ask if system has crashed or stopped before another process was finished?

print "\nScanning for partial pickle cross-correlation files ... "

#maybe create pause statement for interesting load menu.

# loading cross-correlations (looking for *.part.pickle files in folders in
#in dir *CROSSCORR_DIR*)
folder_list = sorted(glob.glob(os.path.join(CROSSCORR_DIR, '*')))

pickle_list = []
index = 0 #creating index for call 
for folder in folder_list:
    #check to see if there are any pickle files in the xcorr time folder 
    if len(glob.glob(os.path.join(folder, '*.part.pickle'))) < 1:
        #print("There are no .pickle files in this folder. Skipping ...")
        continue
    else:
        #append name of pickle file path location string to pickle_list 
        pickle_list.append(glob.glob(os.path.join(folder, \
        '*.part.pickle'))[0])

if len(pickle_list) < 1: 
    print("\nThere are no partial pickle files to begin again from.")
    print("\nThe program will start from the beginning")
    res = ""
else:
    print "\nPlease choose a file to begin again from, or a combination thereof."
    print "Else hit enter to continue anew"
    #print combinations of partial pickle files available
    print '\n0 - All except backups (*~)'    
    print '\n'.join('{} - {}'.format(i + 1, os.path.basename(f))
        for i, f in enumerate(folder_list))
            
            
    #change folder_list to pickle_list if this gives problems
    res = False#raw_input('\n')

#IF LIST INDEX OUT OF RANGE START PROGRAM ALSO    

#if beginning again, reset time-series intervals to the where the selected 
# .part.pickle file left off! 

if not res:
    # ========================================
    #set output file name as normal
    # ========================================
    time_string = str(time.strftime("%d.%m.%Y") + "-" + time.strftime("%X"))
    responsefrom = []
    if USE_DATALESSPAZ:
        responsefrom.append('datalesspaz')
    if USE_STATIONXML:
        responsefrom.append('xmlresponse')
    OUTBASENAME_PARTS = [
    'XCORR-STACK',
    '-'.join(s for s in CROSSCORR_STATIONS_SUBSET) \
    if CROSSCORR_STATIONS_SUBSET else None,
    '{}-{}'.format(FIRSTDAY.strftime("%d.%m.%Y"),
                   LASTDAY.strftime("%d.%m.%Y")),
    '1bitnorm' if ONEBIT_NORM else None,
    '+'.join(responsefrom)
    ]
    
    OUTFILESNAME = '_'.join(p for p in OUTBASENAME_PARTS if p)

    OUTFILESPATH = os.path.join(CROSSCORR_DIR, time_string, OUTFILESNAME)

    OUTFOLDERS = os.path.join(CROSSCORR_DIR, 
                              time_string,  
                              'XCORR_PLOTS')
                              
    #create unique folder in CROSS output folder named by the present time.
    if not os.path.exists(OUTFOLDERS):\
    os.makedirs(OUTFOLDERS)
    
    METADATA_PATH = '{}metadata.pickle'.format(OUTFILESPATH.\
              replace(os.path.basename(OUTFILESPATH), ""))

else:
    # ========================================
    #reset time as previous time, reset output paths as previous path name
    #reset cross-correlation dictionaries 
    # ========================================
    
    PART_PICKLE = pickle_list[int(res)-1]
    OUTFILESPATH = PART_PICKLE[:-12]
     
    # re-initialising .part.pickle collection of cross-correlations
    xc = pscrosscorr.load_pickled_xcorr(PART_PICKLE)
    #read in metadata to find latest time slot. Then assign this to FIRSTDAY
    METADATA_PATH = '{}metadata.pickle'.format(OUTFILESPATH.\
              replace(os.path.basename(OUTFILESPATH), ""))
  
    metadata = pscrosscorr.load_pickled_xcorr(METADATA_PATH)

    #re-assign FIRSTDAY variable to where the data was cut off
    FIRSTDAY = metadata[len(metadata) - 1] +  \
    dt.timedelta(minutes=XCORR_INTERVAL)




# ============
# Main program
# ============

# Reading inventories in dataless seed and/or StationXML files
dataless_inventories = []
if USE_DATALESSPAZ:
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        dataless_inventories = psstation.get_dataless_inventories(DATALESS_DIR,
                                                                  verbose=True)

xml_inventories = []
if USE_STATIONXML:
    xml_inventories = psstation.get_stationxml_inventories(STATIONXML_DIR,
                                                           verbose=True)

# Getting list of stations
#stations, subdir_len = psstation.get_stations(mseed_dir=MSEED_DIR,
#                                  xml_inventories=xml_inventories,
#                                  dataless_inventories=dataless_inventories,
#                                  startday=FIRSTDAY,
#                                  endday=LASTDAY,
#                                  verbose=False)

#connect SQL database
SQL_db = os.path.join(DATABASE_DIR, 'timeline.db')

stations, subdir_len = psstationSQL.get_stationsSQL(SQL_db, 
                                  xml_inventories=xml_inventories,
                                  dataless_inventories=dataless_inventories,
                                  startday=FIRSTDAY,
                                  endday=LASTDAY,
                                  verbose=False)


# Loop on time interval
 #number of time steps
N = int(((LASTDAY - FIRSTDAY).days + 1)*60*24 / XCORR_INTERVAL)

dates = [FIRSTDAY + dt.timedelta(minutes=i) for i in \
         [j*XCORR_INTERVAL for j in range(N)]]


#begin = raw_input("\nPress enter to begin the program ")

# initialise preprocess class: METHOD - Bensen et al. (2007)
Preprocess = pspreprocess.Preprocess(FREQMIN, FREQMAX,
                                     FREQMIN_EARTHQUAKE,
                                     FREQMAX_EARTHQUAKE,
                                     CORNERS, ZEROPHASE,
                                     PERIOD_RESAMPLE,
                                     WINDOW_TIME,
                                     WINDOW_FREQ,
                                     ONEBIT_NORM)
                       
#loop on time-series. Date now represents XCORR_INTERVAL long time intervals
counter = 0 

for date in dates:
    #there may be an overlap in metadata times. Maybe add xcorr interval to this?
    metadata.append(date)
    #have the program restart from the beginning of each stack! not each month
    #this gives more redundancy!
    #(allows to restart after a crash from that date)

    with open(u'{}.part.pickle'.format(OUTFILESPATH), 'wb') as f:
        print "\nExporting cross-correlations calculated until now to: " + f.name
        pickle.dump(xc, f, protocol=2)

        
    #also create a metadata dump file for use only if the program needs to be restarted
    #use replace() to get rid of basename to create file named metadata.pickle in 
    #correct path
    with open(METADATA_PATH, 'wb') as f:
        print "\nExporting re-start metadata of time-series calculated until \
now to: " + f.name
        pickle.dump(metadata, f, protocol=2)

    print "\nProcessing data for date {} with a {} minute cross-correlation\
 time-interval between times: {} and {}".format(date.date(), \
    int(XCORR_INTERVAL)  , date.time(), \
    (date + dt.timedelta(minutes=XCORR_INTERVAL)).time())
    
    iterate_stations = sorted(sta for sta in stations)

    # subset if stations (if provided)
    if CROSSCORR_STATIONS_SUBSET:
        iterate_stations = [sta for sta in iterate_stations
                          if sta.name in CROSSCORR_STATIONS_SUBSET]
                              
    # =================
    # processing traces
    # =================                          

    # =============================================================
    # preparing functions that get one merged trace per station
    # and pre-process trace, ready to be parallelized (if required)
    # =============================================================

    def get_merged_trace(station):
        """
        Preparing func that returns one trace from selected station,
        at current date. Function is ready to be parallelized.
        """

        try:
            trace = Preprocess.get_merged_trace(station=station,
                                                 date=date,
                                                 xcorr_interval=XCORR_INTERVAL,
                                                 skiplocs=CROSSCORR_SKIPLOCS,
                                                 minfill=MINFILL)
    
            errmsg = None
        except pserrors.CannotPreprocess as err:
            # cannot preprocess if no trace or daily fill < *minfill*
            trace = None
            errmsg = '{}: skipping'.format(err)
        except Exception as err:
            # unhandled exception!
            trace = None
            errmsg = 'Unhandled error: {}'.format(err)

        if errmsg:
            # printing error message
            print '{}.{} [{}] '.format(station.network, station.name, errmsg),

        return trace

    def preprocessed_trace((trace, response)):
        """
        Preparing func that returns processed trace: processing includes
        removal of instrumental response, band-pass filtering, demeaning,
        detrending, downsampling, time-normalization and spectral whitening
        (see pscrosscorr.preprocess_trace()'s doc)

        Function is ready to be parallelized.
        """
        

        if not trace or response is False:
            return
        

        try:
            Preprocess.preprocess_trace(trace=trace, paz=response)
            msg = 'ok'
            print '{}.{} [{}] '.format(trace.stats.network, 
                                        trace.stats.station, 
                                        msg),

        except pserrors.CannotPreprocess as err:
            # cannot preprocess if no instrument response was found,
            # trace data are not consistent etc. (see function's doc)
            trace = None
            print(err)

            print 'skipping'

        except Exception as err:
            # unhandled exception!
            trace = None
            print(err)
            print 'skipping'

        # printing output (error or ok) message

        # although processing is performed in-place, trace is returned
        # in order to get it back after multi-processing
        return trace

    # ====================================
    # getting one merged trace per station
    # ====================================
    
    t0 = dt.datetime.now()
    if MULTIPROCESSING['merge trace']:
        # multiprocessing turned on: one process per station
        pool = mp.Pool(None)
        traces = pool.map(get_merged_trace, iterate_stations)
        pool.close()
        pool.join()
    else:
        # multiprocessing turned off: processing stations one after another
        traces = [get_merged_trace(s) for s in iterate_stations]

    # =====================================================
    # getting or attaching instrumental response
    # (parallelization is difficult because of inventories)
    # =====================================================


    responses = []
    for tr in traces:

        
        if not tr:
            responses.append(None)
            continue

        # responses elements can be (1) dict of PAZ if response found in
        # dataless inventory, (2) None if response found in StationXML
        # inventory (directly attached to trace) or (3) False if no
        # response found

        try:
            response = Preprocess.get_or_attach_response(
                trace=tr,
                dataless_inventories=dataless_inventories,
                xml_inventories=xml_inventories)
            errmsg = None
        except pserrors.CannotPreprocess as err:
            # response not found
            response = False
            errmsg = '{}: skipping'.format(err)
        except Exception as err:
            # unhandled exception!
            response = False
            errmsg = 'Unhandled error: {}'.format(err)

        responses.append(response)
        if errmsg:
            # printing error message
            print '{}.{} [{}] '.format(tr.stats.network, 
                                       tr.stats.station, 
                                       errmsg),

    # =================
    # processing traces
    # =================


        
    if MULTIPROCESSING['process trace']:
        # multiprocessing turned on: one process per station
        pool = mp.Pool(NB_PROCESSES)
        traces = pool.map(preprocessed_trace, zip(traces, responses))
        pool.close()
        pool.join()
    else:
        # multiprocessing turned off: processing stations one after another
        traces = [preprocessed_trace((tr, res)) for tr, 
                  res in zip(traces, responses)]

    # setting up dict of current date's traces, {station: trace}
    tracedict = {s.name: trace for s, trace in zip(iterate_stations, 
                                                   traces) if trace}

    delta = (dt.datetime.now() - t0).total_seconds()
    print "\nProcessed stations in {:.1f} seconds".format(delta)
    
    # create tmp folder for tracedict
    if not os.path.exists('tmp'): os.makedirs('tmp')   
    #dump the time interval's pre-processed items in tracedict to a pickle
    with open('tmp/preprocessed_tracedict.pickle', 'wb') as f:
        print "\nExporting pre-processed traces of time-series to: " + f.name
        pickle.dump(tracedict, f, protocol=2)
    
                                   
    # import tracedict from output pickle produced with preprocess_total
    tracedict_pickle = 'tmp/preprocessed_tracedict.pickle'
    f = open(name=tracedict_pickle, mode='rb')
    tracedict = pickle.load(f)
    f.close()
    
    # remove preprocessed tracedict pickle file
    if os.path.isfile(tracedict_pickle): os.remove(tracedict_pickle)
        
        
    # ======================================================
    # stacking cross-correlations of the current time-series
    # ======================================================

    
    if len(tracedict) < 2:
        print "No cross-correlation for this day"
        continue

    t0 = dt.datetime.now()
    xcorrdict = {}
    if MULTIPROCESSING['cross-corr']:
        # if multiprocessing is turned on, we pre-calculate cross-correlation
        # arrays between pairs of stations (one process per pair) and feed
        # them to xc.add() (which won't have to recalculate them)
        print "Pre-calculating cross-correlation arrays"

        def xcorr_func(pair):
            """
            Preparing func that returns cross-correlation array
            beween two traces
            """
            (s1, tr1), (s2, tr2) = pair
            print '{}-{} '.format(s1, s2),
            shift = int(CROSSCORR_TMAX / PERIOD_RESAMPLE)
            xcorr = obspy.signal.cross_correlation.xcorr(
                tr1, tr2, shift_len=shift, full_xcorr=True)[2]
            return xcorr

        pairs = list(it.combinations(sorted(tracedict.items()), 2))
        pool = mp.Pool(NB_PROCESSES)
        xcorrs = pool.map(xcorr_func, pairs)
        pool.close()
        pool.join()
        xcorrdict = {(s1, s2): xcorr for ((s1, _), (s2, _)), 
                     xcorr in zip(pairs, xcorrs)}
        print
        
    #print "Stacking cross-correlations"
    xc.add(tracedict=tracedict,
           stations=stations,
           xcorr_tmax=CROSSCORR_TMAX,
           xcorrdict=xcorrdict,
           date=date,
           verbose=not MULTIPROCESSING['cross-corr'])
    


#==============================================================================    

    print "Calculated and stacked cross-correlations in \
{:.1f} seconds".format(delta)
           

# exporting cross-correlations
if not xc.pairs():
    print "No cross-correlation could be calculated: nothing to export!"
else:
    # exporting to binary and ascii files
    xc.export(outprefix=OUTFILESPATH, stations=stations, verbose=True)

    # exporting to png file
    print "Exporting cross-correlations to file: {}.png".format(OUTFILESPATH)
    # optimizing time-scale: max time = max distance / vmin (vmin = 2.5 km/s)
    maxdist = max([xc[s1][s2].dist() for s1, s2 in xc.pairs()])
    maxt = min(CROSSCORR_TMAX, maxdist / 2.5)
    
    #plot distance plot of cross-correlations
    xc.plot(plot_type='distance', xlim=(-maxt, maxt), 
            outfile=os.path.join(OUTFOLDERS, OUTFILESNAME)\
            + '.png', showplot=False)
    #plot individual cross-correlations
    #change their names!
    xc.plot(plot_type='classic', xlim=(-maxt, maxt), 
            outfile=os.path.join(OUTFOLDERS, OUTFILESNAME)\
            + '.png', showplot=False)

# removing file containing periodical exports of cross-corrs
try:
    os.remove(u'{}.part.pickle'.format(OUTFILESPATH))
except:
    pass