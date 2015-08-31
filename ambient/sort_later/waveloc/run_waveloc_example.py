import logging
from waveloc.options import WavelocOptions
from waveloc.SDS_processing import do_SDS_processing_setup_and_run
from waveloc.migration import do_migration_setup_and_run
from waveloc.locations_trigger import do_locations_trigger_setup_and_run
from waveloc.plot_locations2 import do_plotting_setup_and_run
import os


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
    
def comb_strings(str_list):
    
    out_str = ''    
    for i in str_list: 
        out_str += i + ','
        
    return out_str[:-1] #gets rid of the last comma 
    
    

logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s : %(asctime)s : %(message)s')

# set up default parameters
wo = WavelocOptions()
wo.opdict['base_path'] = os.getcwd()
#print wo.opdict['base_path']
# set base path to $WAVELOC_PATH
wo.verify_base_path()

##########################################
# set waveloc options
##########################################




wo.opdict['time'] = True
wo.opdict['verbose'] = False

wo.opdict['test_datadir'] = '/storage/ANT/PROGRAMS/ANT_OUTPUT/INPUT/DATA'

# create list of all raw file names that are being processed in test_datadir
extension = 'mseed'
nets, stats, chans = basenames(wo.opdict['test_datadir'], extension)

nets_str = comb_strings(nets)
stats_str = comb_strings(stats)
chans_str = comb_strings(chans)

wo.opdict['datadir'] = 'EXAMPLE'
wo.opdict['outdir'] = 'EXAMPLE_fullRes'

# set up automatic station list and network list populator from raw data!

wo.opdict['net_list'] = nets_str
wo.opdict['sta_list'] = stats_str
wo.opdict['comp_list'] = chans_str


# set up max and min times populator from raw data!

wo.opdict['starttime'] = "2014-01-01T00:00:00.0Z"
wo.opdict['endtime'] = "2014-02-01T00:00:00.0Z"

wo.opdict['time_grid'] = 'Slow_len.100m.P'
wo.opdict['search_grid'] = 'grid.Taisne.search.hdr'
wo.opdict['stations'] = 'coord_stations_test'

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
