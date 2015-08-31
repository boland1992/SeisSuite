# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 13:30:37 2015

@author: boland
"""

import os

#maybe write an initial part to check if the correct directories already exist?

# download the necessary packages for waveloc
os.system('wget -O waveloc.tar.gz \
https://s3.amazonaws.com/waveloc/waveloc-0.2.3.tar.gz')
# unpack the tar.gz ball into a folder in the same location as the download file
os.system('tar -xvzf waveloc.tar.gz ')

# remove tar ball once unpacked
os.system('rm -r waveloc.tar.gz')

# compile all necessary C applications for nonlinloc 
os.system('cd waveloc-0.2.3 && python setup.py install')

print "\nFinished installing waveloc"

#print "\nNow running examples"
#os.system('cd examples && python setup_examples.py')


