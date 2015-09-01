import logging
from waveloc.options import WavelocOptions
from waveloc.SDS_processing import do_SDS_processing_setup_and_run
from waveloc.migration import do_migration_setup_and_run
from waveloc.locations_trigger import do_locations_trigger_setup_and_run
from waveloc.plot_locations2 import do_plotting_setup_and_run
import waveloc.hdf5_grids as hdf5 
import glob
import os
from obspy import UTCDateTime
import datetime as dt 



#==============================================================================
# INITIAL SET-UP
#==============================================================================

# please set the estimated date-time extents of your data! 
# (these will be split into hour chucks) FORMAT MUST REMAIN THE SAME!
# please choose a time bracket that has exact hour increments! otherwise
# residuals at the end will be unaccounted for. 
starttime =  "01/01/2014"
endtime =  "01/02/2014"


starttime = dt.datetime.strptime(starttime, '%d/%m/%Y')
endtime = dt.datetime.strptime(endtime, '%d/%m/%Y')

#endtime = starttime + datetime.timedelta(hours=1)

#print starttime
#print endtime

# create loop of all applicable hours to process!
# get the number of hours of all applicable processes

n_hours = int((endtime - starttime).days) * 24

# Loop on hour interval
times = [starttime + dt.timedelta(hours=1)*hour for hour in range(n_hours)]


def basenames(folder_path, extension):
    """
    Function that returns a list of desired absolute paths called abs_paths
    of files that contains a given extension e.g. .txt should be entered as
    folder_path, txt. This function will run recursively through and find
    any and all files within this folder with that extension!
    """
    nets = []
    stats = []
    chans = []
    for root, dirs, files in os.walk(folder_path):
        for f in files:
            #fullpath = os.path.join(root, f)
            
            if os.path.splitext(f)[1] == '.{}'.format(extension):
                comps = f.split('.')
                if comps[0] not in nets:
                    nets.append(comps[0])
                if comps[1] not in stats:
                    stats.append(comps[1])
                if comps[2] not in chans:
                    chans.append(comps[2])                

                    
    return nets, stats, chans
    

def comb_strings(strings):
    emptry_str = ''
    for i in strings:
        emptry_str += i + ', '
        
    return emptry_str[-1]

logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s : %(asctime)s : %(message)s')

# set up default parameters
wo = WavelocOptions()
wo.opdict['base_path'] = '/home/boland/Documents/waveloc_examples' # os.getcwd()
#print wo.opdict['base_path']
# set base path to $WAVELOC_PATH
wo.verify_base_path()

# steps to automating the generation of the necessary grids for waveloc!

#1. generate NLL_control file from raw station metadata
#2. run checks that this control file won't fail
#3. run the NLL programme Vel2Grid

#./Vel2Grid '/home/boland/Desktop/Link to SEIS_SUITE/ambient-v0.1.1/ambient/nonlinloc/control_files/NLL_control_Vel2Grid.in'
#./Grid2Time '/home/boland/Desktop/Link to SEIS_SUITE/ambient-v0.1.1/ambient/nonlinloc/control_files/NLL_control_Vel2Grid.in'
#4. run the NLL programme Grid2Time
#5. Link waveloc to this dir with the outputs of both Vel2Grid and Grid2Time
#6. Convert all .time.hdr and .time.buf files to .hdf5 for use with waveloc


# find all .hdr and .buf files from /lib/ using glob module
bufs_list = glob.glob('lib/*.time.buf')

if len(bufs_list) > 0:
    # convert .buf and .hdr files to .hdf5 grid time files!
    for name in glob.glob('lib/*.time.buf'):
        nll_name = name[:-4]
        h5_name = nll_name + '.hdf5'
        hdf5.nll2hdf5(nll_name, h5_name)
    # remove the .time.buf and .time.hdr files!
    os.system('rm -r lib/*.time.buf'); os.system('rm -r lib/*.time.hdr')


print 'a1'
print times[:1]
for time in times[:1]:
    ##########################################
    # set waveloc options
    ##########################################
    print 'a2'

    wo.opdict['time'] = True
    wo.opdict['verbose'] = False

    wo.opdict['test_datadir'] = '/home/boland/Documents/waveloc_examples/test_data/raw_data' #'/storage/ANT/PROGRAMS/ANT_OUTPUT/INPUT/DATA'

    # create list of all raw file names that are being processed in test_datadir
    extension = 'MSEED'
    nets, stats, chans = basenames(wo.opdict['test_datadir'], extension)

    nets_str = comb_strings(nets)
    stats_str = comb_strings(stats)
    chans_str = comb_strings(chans)

    wo.opdict['datadir'] = 'EXAMPLE'
    wo.opdict['outdir'] = 'EXAMPLE_fullRes'

    # set up automatic station list and network list populator from raw data!
    print nets_str
    wo.opdict['net_list'] = nets_str
    wo.opdict['sta_list'] = stats_str
    wo.opdict['comp_list'] = chans_str
    
    # set up max and min times populator from raw data!
    #wo.opdict['starttime'] =  str(UTCDateTime(time))
    #wo.opdict['endtime'] =  str(UTCDateTime(time + dt.timedelta(hours=1)))
    
    wo.opdict['starttime'] = '2010-10-14T00:14:00.00'
    wo.opdict['endtime'] =  '2010-10-17T00:14:00.00'
    
    wo.opdict['time_grid'] = 'Slow_len.100m.P'
    wo.opdict['search_grid'] = 'grid.Taisne.search.hdr'
    
    #automate the creation of UTM from lat-lon station txt document
    wo.opdict['stations'] = 'station_info'

    wo.opdict['resample'] = False
    wo.opdict['fs'] = None

    wo.opdict['c1'] = 4.0
    wo.opdict['c2'] = 10.0

    wo.opdict['kwin'] = 4
    wo.opdict['krec'] = False
    wo.opdict['kderiv'] = True

    wo.opdict['data_length'] = 300
    wo.opdict['data_overlap'] = 20
    
    wo.opdict['dataglob'] = '*filt.mseed'
    wo.opdict['kurtglob'] = '*kurt.mseed'
    wo.opdict['gradglob'] = '*grad.mseed'
    
    wo.opdict['load_ttimes_buf'] = True

    wo.opdict['loclevel'] = 5000.0
    wo.opdict['snr_limit'] = 10.0
    wo.opdict['sn_time'] = 10.0
    wo.opdict['n_kurt_min'] = 4

    wo.opdict['plot_tbefore'] = 4
    wo.opdict['plot_tafter'] = 6
    wo.opdict['plot_otime_window'] = 2

    ##########################################
    # end of option setting - start processing
    ##########################################

    wo.verify_SDS_processing_options()
    do_SDS_processing_setup_and_run(wo.opdict)

    wo.verify_migration_options()
    do_migration_setup_and_run(wo.opdict)

    # do trigger location
    wo.verify_location_options()
    do_locations_trigger_setup_and_run(wo.opdict)

    # This will do plotting of grids and stacks for locations
    wo.verify_plotting_options()
    do_plotting_setup_and_run(wo.opdict, plot_wfm=True, plot_grid=True)