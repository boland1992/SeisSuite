# -*- coding: utf-8 -*-
"""
Created on Sun May  3 13:14:51 2015

@author: boland
"""
import os
import sqlite3 as lite

# ====================================================
# SETUP INITIAL FILE STRUCTURE
# ====================================================


FOLDER = '/home/iese/Documents/Ben/WAIRAKEI'
DATABASE_DIR = os.path.join(FOLDER, 'DATABASES')


#try:
#from seissuite.ant.psconfig import FOLDER, DATABASE_DIR

#except Exception as error:
#    print error
#    print "If there was an import error, please install Anaconda\n\
#http://continuum.io/downloads"

if FOLDER == 'DEFAULT':
    FOLDER = os.getcwd()    
    
print FOLDER

if os.path.exists(FOLDER) and FOLDER != os.getcwd():
     raise Exception("\nThis folder path is already in use. Please go back to config\
 file and choose another path.\n")
 
#folder is cwd if set to false in config file
if FOLDER == "False":
    #FOLDER is the parent directory to the psconfig directory
    FOLDER = os.getcwd()

#first check if FOLDER_PATH is a true path
if not os.path.isabs(FOLDER):
    print("\nWhat you have entered into the FOLDER variable in the config\
 file is not a true path. Please go back and enter a true path\n")    
    quit()

else: 
    #create file structure
    INPUT = "{}/INPUT".format(FOLDER)
    OUTPUT = "{}/OUTPUT".format(FOLDER)
    input_dirs = ["DATA", "DATALESS", "XML", "DATABASES"]
    output_dirs = ["CROSS", "FTAN", "TOMO", "DEPTH"]

    
    #create input file if not already created
    if not os.path.exists(INPUT):\
    os.makedirs(INPUT)
    
    for i in input_dirs:
        if not os.path.exists("{}/{}".format(INPUT, i)):\
        os.makedirs("{}/{}".format(INPUT, i))    
    
    #create output file if not already created
    if not os.path.exists(OUTPUT):\
    os.makedirs(OUTPUT)
    
    for i in output_dirs:
        if not os.path.exists("{}/{}".format(OUTPUT, i)):\
        os.makedirs("{}/{}".format(OUTPUT, i))   





