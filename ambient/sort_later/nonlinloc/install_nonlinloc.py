# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 13:30:37 2015

@author: boland
"""

import os

#maybe write an initial part to check if the correct directories already exist?

# create tmp folder for storage of the latest tar.gz file 


# download the necessary packages for nonlinloc
os.system('wget -O nonlinloc.tar.gz http://alomax.free.fr/nlloc/\
soft6.00/tar/NLL6.00_src.tgz')

# unpack the tar.gz ball into a folder in the same location as the download file
os.system('tar -xvzf nonlinloc.tar.gz')

# remove tar ball once unpacked
os.system('rm -r nonlinloc.tar.gz')

# compile all necessary C applications for nonlinloc 
os.system('cd src && make -R all')

print "\nFinished compiling and installing Nonlinloc"
print "\nHAPPY HUNTING!"



#newer than local copy python  wget ‐‐continue ‐‐timestamping wordpress.org/latest.zip

