"""
Module that parses global parameters from a configuration
file at first import, to make them available to the other
parts of the program.
"""

import ConfigParser
import os
import glob
import json
import datetime as dt
import numpy as np

# cPickle is faster than pickle but doesn't always import correctly. 
try:
    import cPickle as pickle
except:
    import pickle
    

cnf_path = os.path.join(os.getcwd(), 'configs')
print cnf_path
def cnf_list(cnf_path, ext='cnf'):
    return glob.glob(os.path.join(cnf_path, u'*.{}'.format(ext)))


def create_config_list(basedir=os.path.join(os.getcwd(), 'configs'), 
                       ext='cnf', verbose=True):
    """
    Function that runs the individual .cnf file if only one exists, 
    otherwise checks for AUTOMATE option in the .cnf and if any of the
    .cnf files have AUTOMATE = True then both will run in sequence 
    (in named order of .cnf files), finally run an user input if there
    is more than one .cnf file and ALL state AUTOMATE = False
    
    This function will return a list of only the .cnf files to be processed
    """
    
    config_files = cnf_list(basedir, ext='cnf')
    
    condition = False
    config = ConfigParser.ConfigParser()
    
    for config_file in config_files:
        print 'config_file: ', config_file
        config.read(config_file)
        AUTOMATE = config.getboolean('automate', 'AUTOMATE')
        if AUTOMATE:
            condition = True
            
    if condition or len(config_files) == 1:
        # if automate is set in ANY of the control files, the the whole list
        # will be returned!
        config_list = config_files
        
    elif len(config_files) == 1:
        raise Exception('There are no configuration files in the \
./configs/ folder')

    else:
        # run user input if automation is False for ALL files.
        print "Please select a configuration file:"
        for i, f in enumerate(config_files, start=1):
            print "{} - {}".format(i, os.path.basename(f))
        res = int(raw_input(''))
        config_file = config_files[res - 1]
        config_list = [config_file]

    return config_list
    
    
def shift(xcorr_len):
    
    a = 10.0 / 837.0
    b = 3.0e4 / 31.0
    
    shift_len = int(a * xcorr_len + b)
    
    return shift_len



class Config:
    """
    The following class contains all the required configuration information
    required to process the whole set of applications!
    """
    
    def __init__(self, config_file):
    
        # initialise config file parser object
        self.config = ConfigParser.ConfigParser()
        self.config.read(config_file)
        
        # initialise ALL configuration variables on definition of class!
        # -----
        # paths
        # -----
        
        self.AUTOMATE = self.config.getboolean('automate', 'AUTOMATE')    
        
        self.FOLDER = self.config.get('paths', 'FOLDER')
        
        #TIMELINE_DB = config.get('paths', 'TIMELINE_DB')
        #RESPONSE_DB = config.get('paths', 'RESPONSE_DB')
        
        
        # input dirs
        if self.FOLDER == "DEFAULT":
            #fold = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
            fold = os.getcwd()
            self.MSEED_DIR = "{}/INPUT/DATA".format(fold)
            self.STATIONXML_DIR = "{}/INPUT/XML".format(fold)
            self.DATALESS_DIR = "{}/INPUT/DATALESS".format(fold)
            self.DATABASE_DIR = "{}/INPUT/DATABASES".format(fold)
            # output dirs
            self.CROSSCORR_DIR = "{}/OUTPUT/CROSS".format(fold)
            self.FTAN_DIR = "{}/OUTPUT/FTAN".format(fold)
            self.TOMO_DIR = "{}/OUTPUT/TOMO".format(fold)
            self.DEPTHMODELS_DIR = "{}/OUTPUT/DEPTH".format(fold)
        
        else:
            self.MSEED_DIR = "{}/INPUT/DATA".format(self.FOLDER)
            self.STATIONXML_DIR = "{}/INPUT/XML".format(self.FOLDER)
            self.DATALESS_DIR = "{}/INPUT/DATALESS".format(self.FOLDER)
            self.DATABASE_DIR = "{}/INPUT/DATABASES".format(self.FOLDER)
        
            # output dirs
            self.CROSSCORR_DIR = "{}/OUTPUT/CROSS".format(self.FOLDER)
            self.FTAN_DIR = "{}/OUTPUT/FTAN".format(self.FOLDER)
            self.TOMO_DIR = "{}/OUTPUT/TOMO".format(self.FOLDER)
            self.DEPTHMODELS_DIR = "{}/OUTPUT/DEPTH".format(self.FOLDER)  
            
        # dir of the Computer Programs in Seismology (can be None)
        self.COMPUTER_PROGRAMS_IN_SEISMOLOGY_DIR = \
        self.config.get('paths', 'COMPUTER_PROGRAMS_IN_SEISMOLOGY_DIR')
        #===========
        # processing
        #===========
        
        #set the individual preprocessing techniques that you want performed on your analysis. Each
        # must be set either True or False to work. Any other options with give an error
        self.MAX_DISTANCE = self.config.get('processing', 'MAX_DISTANCE')
        self.TDD = self.config.getboolean('processing', 'TDD')	 

        #self.EVENT_REMOVE = self.config.getboolean('processing', 
        #                                           'EVENT_REMOVE')
        self.RESP_REMOVE = self.config.getboolean('processing', 
                                                  'RESP_REMOVE')    
        self.HIGHAMP_REMOVE = self.config.getboolean('processing', 
                                                     'HIGHAMP_REMOVE')
        self.RESP_CHECK = self.config.getboolean('processing', 
                                                 'RESP_CHECK')
                                                 
        self.RESP_RANGE = json.loads(self.config.get('processing', 
                                                     'RESP_RANGE'))  
                                                     
        self.RESP_FREQS = [1./max(self.RESP_RANGE), 1./min(self.RESP_RANGE)]       
        
        self.RESP_TOL =  self.config.getfloat('processing', 'RESP_TOL')
        self.RESP_EFFECT = self.config.getfloat('processing', 'RESP_EFFECT')
        
        self.BANDPASS = self.config.getboolean('processing', 'BANDPASS')           
        self.DOWNSAMPLE = self.config.getboolean('processing', 'DOWNSAMPLE')         
        self.COMPLETENESS = self.config.getboolean('processing', 
                                                   'COMPLETENESS')       
        self.TIME_NOMALISATION = self.config.getboolean('processing', 
                                                        'TIME_NOMALISATION')  
        self.SPEC_WHITENING = self.config.getboolean('processing', 
                                                     'SPEC_WHITENING')   
        
        self.FULL_COMB = self.config.getboolean('processing', 
                                                     'FULL_COMB')   
        # ---------------
        # maps parameters
        # ---------------
        
        # paths to shapefiles (coasts, tectonic provinces and labels)
        self.COAST_SHP = self.config.get('maps', 'COAST_SHP')
        self.TECTO_SHP = self.config.get('maps', 'TECTO_SHP')
        self.TECTO_LABELS = self.config.get('maps', 'TECTO_LABELS')
        
        # colors of tectonic provinces
        self.TECTO_COLORS = json.loads(self.config.get('maps', 'TECTO_COLORS'))
        
        # bounding boxes
        self.BBOX_LARGE = json.loads(self.config.get('maps', 'BBOX_LARGE'))
        self.BBOX_SMALL = json.loads(self.config.get('maps', 'BBOX_SMALL'))
        
        # --------------------------------------
        # cross-correlation / spectra parameters
        # --------------------------------------
        
        
        self.RANDOM_STACK = self.config.getboolean('cross-correlation', 
                                                   'RANDOM_STACK')
                                                 
                                                 
        # use dataless files or stationXML files to remove instrument response?
        self.USE_DATALESSPAZ = self.config.getboolean('cross-correlation', 
                                                 'USE_DATALESSPAZ')
        self.USE_STATIONXML = self.config.getboolean('cross-correlation', 
                                                'USE_STATIONXML')
        
        # subset of stations to cross-correlate
        self.CROSSCORR_STATIONS_SUBSET = \
        self.config.get('cross-correlation', 'CROSSCORR_STATIONS_SUBSET')
        
        self.CROSSCORR_STATIONS_SUBSET = \
        json.loads(self.CROSSCORR_STATIONS_SUBSET)
        
        # locations to skip
        self.CROSSCORR_SKIPLOCS = \
        json.loads(self.config.get('cross-correlation', 'CROSSCORR_SKIPLOCS'))
        
        #GET RID OF .day() for FIRSTDAY and LASTDAY variables. 
        #This is is to allow for the interval to be a datetime object rather than just
        #a date object!
        
        # first and last day, minimum data fill per day
        self.FIRSTDAY = self.config.get('cross-correlation', 'FIRSTDAY')
        self.FIRSTDAY = dt.datetime.strptime(self.FIRSTDAY, '%d/%m/%Y')
        self.LASTDAY = self.config.get('cross-correlation', 'LASTDAY')
        self.LASTDAY = dt.datetime.strptime(self.LASTDAY, '%d/%m/%Y')
        self.MINFILL = self.config.getfloat('cross-correlation', 'MINFILL')
        
        # band-pass parameters
        self.PERIODMIN = self.config.getfloat('cross-correlation', 
                                              'PERIODMIN')
        self.PERIODMAX = self.config.getfloat('cross-correlation', 
                                              'PERIODMAX')
        self.FREQMIN = 1.0 / self.PERIODMAX
        self.FREQMAX = 1.0 / self.PERIODMIN
        self.CORNERS = self.config.getint('cross-correlation', 'CORNERS')
        self.ZEROPHASE = self.config.getboolean('cross-correlation', 
                                                'ZEROPHASE')
        # resample period (to decimate traces, after band-pass)
        self.PERIOD_RESAMPLE = self.config.getfloat('cross-correlation', 
                                               'PERIOD_RESAMPLE')
        
        # Time-normalization parameters:
        self.ONEBIT_NORM = self.config.getboolean('cross-correlation', 
                                                  'ONEBIT_NORM')
        # earthquakes period bands
        self.PERIODMIN_EARTHQUAKE = self.config.getfloat('cross-correlation', 
                                                    'PERIODMIN_EARTHQUAKE')
        self.PERIODMAX_EARTHQUAKE = self.config.getfloat('cross-correlation', 
                                                    'PERIODMAX_EARTHQUAKE')
        
        self.FREQMIN_EARTHQUAKE = 1.0 / self.PERIODMAX_EARTHQUAKE
        self.FREQMAX_EARTHQUAKE = 1.0 / self.PERIODMIN_EARTHQUAKE
        # time window (s) to smooth data in earthquake band
        # and calculate time-norm weights
        self.WINDOW_TIME = 0.5 * self.PERIODMAX_EARTHQUAKE
        
        # frequency window (Hz) to smooth ampl spectrum
        # and calculate spect withening weights
        self.WINDOW_FREQ = self.config.getfloat('cross-correlation', 
                                                'WINDOW_FREQ')
        
        #45 minute long xcorr time intervals for superior SNR
        self.XCORR_INTERVAL = self.config.getfloat('cross-correlation', 
                                              'XCORR_INTERVAL') 
        
        # Max time window (s) for cross-correlation
        #set this to automatic based on function shift() above.
        try:
            self.CROSSCORR_TMAX = self.config.getfloat('cross-correlation', 
                                                  'CROSSCORR_TMAX') 
            
        except:
            self.CROSSCORR_TMAX = shift(self.XCORR_INTERVAL * 60.0)
        
        self.PLOT_CLASSIC = self.config.getboolean('cross-correlation', 
                                                  'PLOT_CLASSIC')
                                                  
        self.PLOT_DISTANCE = self.config.getboolean('cross-correlation', 
                                                  'PLOT_DISTANCE')
        
        
        # ---------------
        # FTAN parameters
        # ---------------
        
        # default period bands, used to:
        # - plot cross-correlation by period bands, in plot_FTAN(), 
        # - plot_by_period_bands()
        # - plot spectral SNR, in plot_spectral_SNR()
        # - estimate min spectral SNR, in FTANs()
        self.PERIOD_BANDS = json.loads(self.config.get('FTAN', 'PERIOD_BANDS'))
        
        # default parameters to define the signal and noise windows used to
        # estimate the SNR:
        # - the signal window is defined according to a min and a max velocity as:
        #   dist/vmax < t < dist/vmin
        # - the noise window has a fixed size and starts after a fixed trailing
        #   time from the end of the signal window
        
        self.SIGNAL_WINDOW_VMIN = self.config.getfloat('FTAN', 
                                                       'SIGNAL_WINDOW_VMIN')
        self.SIGNAL_WINDOW_VMAX = self.config.getfloat('FTAN', 
                                                       'SIGNAL_WINDOW_VMAX')
        self.SIGNAL2NOISE_TRAIL = self.config.getfloat('FTAN', 
                                                       'SIGNAL2NOISE_TRAIL')
        self.NOISE_WINDOW_SIZE = self.config.getfloat('FTAN', 
                                                      'NOISE_WINDOW_SIZE')
        
        # smoothing parameter of FTAN analysis
        self.FTAN_ALPHA = self.config.getfloat('FTAN', 'FTAN_ALPHA')
        
        # periods and velocities of FTAN analysis
        self.RAWFTAN_PERIODS_STARTSTOPSTEP = \
        self.config.get('FTAN', 'RAWFTAN_PERIODS_STARTSTOPSTEP')
        self.RAWFTAN_PERIODS_STARTSTOPSTEP = \
        json.loads(self.RAWFTAN_PERIODS_STARTSTOPSTEP)
        
        self.RAWFTAN_PERIODS = np.arange(*self.RAWFTAN_PERIODS_STARTSTOPSTEP)
        
        self.CLEANFTAN_PERIODS_STARTSTOPSTEP = \
        self.config.get('FTAN', 'CLEANFTAN_PERIODS_STARTSTOPSTEP')
        
        self.CLEANFTAN_PERIODS_STARTSTOPSTEP = \
        json.loads(self.CLEANFTAN_PERIODS_STARTSTOPSTEP)
        
        self.CLEANFTAN_PERIODS = \
        np.arange(*self.CLEANFTAN_PERIODS_STARTSTOPSTEP)
        
        self.FTAN_VELOCITIES_STARTSTOPSTEP = \
        self.config.get('FTAN', 'FTAN_VELOCITIES_STARTSTOPSTEP')
        self.FTAN_VELOCITIES_STARTSTOPSTEP = \
        json.loads(self.FTAN_VELOCITIES_STARTSTOPSTEP)
        self.FTAN_VELOCITIES = np.arange(*self.FTAN_VELOCITIES_STARTSTOPSTEP)
        self.FTAN_VELOCITIES_STEP = self.FTAN_VELOCITIES_STARTSTOPSTEP[2]
        
        # relative strength of the smoothing term in the penalty function that
        # the dispersion curve seeks to minimize
        self.STRENGTH_SMOOTHING = self.config.getfloat('FTAN', 
                                                       'STRENGTH_SMOOTHING')
        
        # replace nominal frequancy (i.e., center frequency of Gaussian filters)
        # with instantaneous frequency (i.e., dphi/dt(t=arrival time) with phi the
        # phase of the filtered analytic signal), in the FTAN and dispersion curves?
        # See Bensen et al. (2007) for technical details.
        self.USE_INSTANTANEOUS_FREQ = \
        self.config.getboolean('FTAN', 'USE_INSTANTANEOUS_FREQ')
        
        # if the instantaneous frequency (or period) is used, we need to discard bad
        # values from instantaneous periods. So:
        # - instantaneous periods whose relative difference with respect to
        #   nominal period is greater than ``MAX_RELDIFF_INST_NOMINAL_PERIOD``
        #   are discarded,
        # - instantaneous periods lower than ``MIN_INST_PERIOD`` are discarded,
        # - instantaneous periods whose relative difference with respect to the
        #   running median is greater than ``MAX_RELDIFF_INST_MEDIAN_PERIOD`` are
        #   discarded; the running median is calculated over
        #   ``HALFWINDOW_MEDIAN_PERIOD`` points to the right and to the left
        #   of each period.
        
        self.MAX_RELDIFF_INST_NOMINAL_PERIOD = \
        self.config.getfloat('FTAN', 'MAX_RELDIFF_INST_NOMINAL_PERIOD')
        self.MIN_INST_PERIOD = self.config.getfloat('FTAN', 'MIN_INST_PERIOD')
        self.HALFWINDOW_MEDIAN_PERIOD = \
        self.config.getint('FTAN', 'HALFWINDOW_MEDIAN_PERIOD')
        self.MAX_RELDIFF_INST_MEDIAN_PERIOD = \
        self.config.getfloat('FTAN', 'MAX_RELDIFF_INST_MEDIAN_PERIOD')
        
        # --------------------------------
        # Tomographic inversion parameters
        # --------------------------------
        
        # Default parameters related to the velocity selection criteria
        
        # min spectral SNR to retain velocity
        self.MINSPECTSNR = self.config.getfloat('tomography', 'MINSPECTSNR')
        # min spectral SNR to retain velocity if no std dev
        self.MINSPECTSNR_NOSDEV = self.config.getfloat('tomography', 
                                                       'MINSPECTSNR_NOSDEV')
        # max sdt dev (km/s) to retain velocity
        self.MAXSDEV = self.config.getfloat('tomography', 'MAXSDEV')
        # min nb of trimesters to estimate std dev
        self.MINNBTRIMESTER = self.config.getint('tomography', 
                                                 'MINNBTRIMESTER')
        # max period = *MAXPERIOD_FACTOR* * pair distance
        self.MAXPERIOD_FACTOR = self.config.getfloat('tomography', 
                                                     'MAXPERIOD_FACTOR')
        
        # Default internode spacing of grid
        self.LONSTEP = self.config.getfloat('tomography', 'LONSTEP')
        self.LATSTEP = self.config.getfloat('tomography', 'LATSTEP')
        
        # Default correlation length of the smoothing kernel:
        # S(r,r') = exp[-|r-r'|**2 / (2 * correlation_length**2)]
        self.CORRELATION_LENGTH = self.config.getfloat('tomography', 
                                                  'CORRELATION_LENGTH')
        
        # Default strength of the spatial smoothing term (alpha) and the
        # weighted norm penalisation term (beta) in the penalty function
        self.ALPHA = self.config.getfloat('tomography', 'ALPHA')
        self.BETA = self.config.getfloat('tomography', 'BETA')
        
        # Default parameter in the damping factor of the norm penalization term,
        # such that the norm is weighted by exp(- lambda_*path_density)
        # With a value of 0.15, penalisation strong when path density < ~20
        # With a value of 0.30, penalisation strong when path density < ~10
        self.LAMBDA = self.config.getfloat('tomography', 'LAMBDA')
        

def run_config(config_file):
    """
    Function that inialises an instance of the Config class with a given 
    input config file path, and then saves a temporary pickle dump 
    of that class for reference from other scripts! This class dump
    will contain all global variables inside the Config class! 
    """
    # set tmp_config.pickle to be saved in the same config file directory
    root_config = os.path.dirname(config_file)
    pickle_path = os.path.join(root_config, 'tmp_config.pickle')
    
    # initialise the config class with the desired config file
    CONFIG = Config(config_file)
    
    with open(pickle_path, 'wb') as f:
        print '\nSaving temporary configuration variables ... '
        pickle.dump(CONFIG, f, protocol=2)
        
        
def remove_config(config_file):
    """
    Function that removes only the 'tmp_config.pickle' file in order to 
    both end the current processing configuration loop, and allow the
    next one to begin fresh. 
    """
    # set tmp_config.pickle to be saved in the same config file directory
    root_config = os.path.dirname(config_file)
    pickle_path = os.path.join(root_config, 'tmp_config.pickle')
    os.remove(pickle_path)
