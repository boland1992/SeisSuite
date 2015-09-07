# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 21:46:00 2015

@author: boland
"""
import os
from seissuite.ant.psconfig import cnf_list, Config

cnf_path = os.path.join(os.getcwd(), 'configs')
config_list = cnf_list(cnf_path)
print config_list

# initialise Config class
CONFIG = Config(config_list[0])
print CONFIG
