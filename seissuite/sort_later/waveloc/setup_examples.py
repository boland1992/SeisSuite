#!/usr/bin/env python
import os
import glob
import logging
from waveloc.make_SDS_data_links import make_SDS_data_links

def paths(folder_path, extension):
    """
    Function that returns a list of desired absolute paths called abs_paths
    of files that contains a given extension e.g. .txt should be entered as
    folder_path, txt. This function will run recursively through and find
    any and all files within this folder with that extension!
    """
    abs_paths = []
    for root, dirs, files in os.walk(folder_path):

        for f in files:
            fullpath = os.path.join(root, f)
            if os.path.splitext(fullpath)[1] == '.{}'.format(extension):
                abs_paths.append(fullpath)

    dirs = []
    for path in abs_paths:
        
        sub_dir = path.split('/')[-2]
        if sub_dir not in dirs:
            dirs.append(sub_dir)
        
    return dirs, abs_paths


test_download = False

logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s : %(asctime)s : %(message)s')

base_path = os.getcwd()
# get basic information
#base_path = os.getenv('WAVELOC_PATH')
#test_data_dir = os.path.join(base_path, 'test_data', 'raw_data')
test_data_dir = '/storage/ANT/PROGRAMS/ANT_OUTPUT/INPUT/DATA/2014-01'

dirs, paths = paths(test_data_dir, 'mseed')

if not os.path.exists('test_data'): 
    os.mkdir('test_data')
    for directory in dirs:
        #make a directory for each month for the previous format! 
        os.mkdir('test_data/raw_{}'.format(directory))
        
else:
    raise Exception('test_data dir already exists. will not overwrite!')

#quit()

def symlink_dest(src):
    """
    Function to create symbolic links from the ambient noise tomography 
    program input data format file structure, to SDS where waveloc 
    can use it. This happens without copying or cutting of the original
    data. 
    """
    
    #dest_basename = '{}.{}.00.{}.MSEED'.format(net, stat, chan) 

#os.symlink()


data_dir = os.path.join(base_path, 'data', 'EXAMPLE')
lib_dir = os.path.join(base_path, 'lib')
if not os.path.exists(data_dir):
    os.makedirs(data_dir)
if not os.path.exists(lib_dir):
    os.makedirs(lib_dir)

# make the data links
make_SDS_data_links(test_data_dir, '*mseed', data_dir)

# make link for test grid file etc
test_files = ['coord_stations_test', 'grid.Taisne.search.hdr']
for tfile in test_files:
    try:
        os.symlink(os.path.join(base_path, 'test_data', tfile),
                   os.path.join(base_path, 'lib', tfile))
        logging.info("Linked %s" % tfile)
    except OSError:
        logging.info("File %s already linked" % 'tfile')
        logging.info("Removing old %s" % tfile)
        os.remove(os.path.join(base_path, 'lib', tfile))
        os.symlink(os.path.join(base_path, 'test_data', tfile),
                   os.path.join(base_path, 'lib', tfile))
        logging.info("Linked %s" % tfile)


    
# make links for PDF time grids
test_files = glob.glob(os.path.join(base_path, 'test_data',
                                    'time_grids', 'Slow*'))
                                    
                                    
#if test_files == [] or (not os.path.exists('test_data') and test_download):
#    logging.error('Download \
#    https://github.com/downloads/amaggi/waveloc/test_data.tgz and \
#    unpack it in the %s directory, then re-run' % (base_path))
    
#    if test_download:
        # download the necessary packages for nonlinloc
#        os.system('wget -O test_data.tar.gz \
#https://github.com/downloads/amaggi/waveloc/test_data.tgz')
        # unpack the tar.gz ball into a folder in the same location as the download file
#        os.system('tar -xvzf test_data.tar.gz ')

        # remove tar ball once unpacked
#        os.system('rm -r test_data.tar.gz')




for tfile in test_files:
    try:
        os.symlink(os.path.join(base_path, 'test_data', 'time_grids',
                                os.path.basename(tfile)),
                   os.path.join(base_path, 'lib', os.path.basename(tfile)))
        logging.info("Linked %s" % tfile)
    except OSError:
        pass
