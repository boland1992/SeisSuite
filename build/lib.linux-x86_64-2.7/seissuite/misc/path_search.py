# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 11:57:14 2015

@author: boland
"""

import os
import datetime

#SCANNING FUNCTIONS 
def paths_sort(path):
    """
    Function defined for customised sorting of the abs_paths list
    and will be used in conjunction with the sorted() built in python
    function in order to produce file paths in chronological order.
    """
    base_name = os.path.basename(path)
    stat_name, date= base_name.split('.')[0], base_name.split('.')[1]     
    try:
        date = datetime.datetime.strptime(date, '%Y-%m-%d')
        return date, stat_name
    except Exception as e:
        print(e)
        
def paths(folder_path, extension, sort=False):
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
    if sort:
        abs_paths = sorted(abs_paths, key=paths_sort)
        
    return abs_paths
    
