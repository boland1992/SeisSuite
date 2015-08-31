# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 13:44:40 2015

@author: boland
"""

import obspy.core
import os
import glob
import pickle
import datetime as dt
from obspy import read
import obspy.xseed
import obspy.signal
from obspy.core import Trace, Stream, UTCDateTime
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np


start = UTCDateTime("2014-01-25")
end = UTCDateTime("2014-01-26") 

MSEED_DIR = '/storage/ANT/INPUT/DATA/AU-2014/2014-01/AU.ARMA.BHZ.mseed'

DATALESS_DIR = '/storage/ANT/INPUT/DATALESS/AU.BHZ.02.2014.821983.dataless'

st = read(MSEED_DIR, startime=start, endtime=end)
st.merge()
tr = st[0]


def spectrum(tr):
    wave = tr.data #this is how to extract a data array from a mseed file
    fs = tr.stats.sampling_rate
    #hour = str(hour).zfill(2) #create correct format for eqstring
    f, Pxx_spec = signal.welch(wave, fs, 'flattop', nperseg=1024, scaling='spectrum')
                    #plt.semilogy(f, np.sqrt(Pxx_spec))
    
    if len(f) >= 256:
        column = np.column_stack((f[:255], np.abs(np.sqrt(Pxx_spec)[:255])))   
        return column
    else:
        return 0.

column = spectrum(tr)


plt.figure()
plt.semilogy(column[:,0], column[:,1])
plt.show()




def get_dataless_inventories(dataless_dir=DATALESS_DIR, verbose=False):
    """
    Reads inventories in all dataless seed (*.dataless) and
    pickle (*.pickle) files of specified dir

    @type dataless_dir: unicode or str
    @type verbose: bool
    @rtype: list of L{obspy.xseed.parser.Parser}
    """
    inventories = []

    # list of *.dataless files
    flist = [DATALESS_DIR]#glob.glob(pathname=os.path.join(dataless_dir, "*.dataless"))

    if verbose:
        if flist:
            print "Reading inventory in dataless seed file:",
        else:
            s = u"Could not find any dalatess seed file (*.dataless) in dir: {}!"
            print s.format(dataless_dir)

    for f in flist:
        if verbose:
            print os.path.basename(f),
        inv = obspy.xseed.Parser(f)
        inventories.append(inv)

    # list of *.pickle files
    flist = glob.glob(pathname=os.path.join(dataless_dir, "*.pickle"))

    if flist and verbose:
        print "\nReading inventory in pickle file:",

    for f in flist:
        if verbose:
            print os.path.basename(f),
        f = open(f, 'rb')
        inventories.extend(pickle.load(f))
        f.close()

    if flist and verbose:
        print

    return inventories

dataless_inventories = get_dataless_inventories(dataless_dir=DATALESS_DIR, 
                                                verbose=True)
                                                                                                
                    

def get_paz(channelid, t, inventories):
    """
    Gets PAZ from list of dataless (or pickled dict) inventories
    @type channelid: str
    @type t: L{UTCDateTime}
    @type inventories: list of L{obspy.xseed.parser.Parser} or dict
    @rtype: dict
    """

    for inv in inventories:
        try:
            if hasattr(inv, 'getPAZ'):
                paz = inv.getPAZ(channelid, t)
            else:
                assert channelid == inv['channelid']
                assert not inv['startdate'] or t >= inv['startdate']
                assert not inv['enddate'] or t <= inv['enddate']
                paz = inv['paz']
        except Exception as error:
            print error
            continue
        else:
            return paz
    else:
        # no matching paz found
        raise Exception("no matching paz found")



def get_or_attach_response(trace, dataless_inventories=(), xml_inventories=()):
    """
    Returns or attach instrumental response, from dataless seed inventories
    (as returned by psstation.get_dataless_inventories()) and/or StationXML
    inventories (as returned by psstation.get_stationxml_inventories()).
    If a response if found in a dataless inventory, then a dict of poles
    and zeros is returned. If a response is found in a StationXML
    inventory, then it is directly attached to the trace and nothing is
    returned.

    Raises CannotPreprocess exception if no instrumental response is found.

    @type trace: L{Trace}
    @param dataless_inventories: inventories from dataless seed files (as returned by
                                 psstation.get_dataless_inventories())
    @type dataless_inventories: list of L{obspy.xseed.parser.Parser}
    @param xml_inventories: inventories from StationXML files (as returned by
                            psstation.get_stationxml_inventories())
    @type xml_inventories: list of L{obspy.station.inventory.Inventory}
    """
    t1 = dt.datetime.now()
    # looking for instrument response...
    try:
        # ...first in dataless seed inventories
        paz = get_paz(channelid=trace.id,
                      t=trace.stats.starttime,
                      inventories=dataless_inventories)
        return paz
    except:
        # ... then in dataless seed inventories, replacing 'BHZ' with 'HHZ'
        # in trace's id (trick to make code work with Diogo's data)
        try:
            paz = get_paz(channelid=trace.id.replace('BHZ', 'HHZ'),
                          t=trace.stats.starttime,
                          inventories=dataless_inventories)
            return paz
        except Exception as error:
            print error

    
    delta = (dt.datetime.now() - t1).total_seconds()
    print "\nProcessed response attachment in {:.1f} seconds".format(delta)


paz = get_or_attach_response(tr, dataless_inventories=dataless_inventories)
tr = tr.simulate(paz_remove=paz,
                 paz_simulate=obspy.signal.cornFreq2Paz(0.01),
                 remove_sensitivity=True,
                 simulate_sensitivity=True,
                 nfft_pow2=True)
                 

column = spectrum(tr)


plt.figure()
plt.semilogy(column[:,0], column[:,1])
plt.show()
quit()
st = Stream(traces=[tr])
st.plot()                 


