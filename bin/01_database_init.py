# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 15:35:46 2015

@author: boland
"""
from seissuite.ant.psconfig import DATABASE_DIR
import os
import sqlite3 as lite
#if no two SQL databases exist, then create them! 
TIMELINE_DB = os.path.join(DATABASE_DIR, 'timeline.db')
RESP_DB = os.path.join(DATABASE_DIR, 'response.db')

print TIMELINE_DB
print RESP_DB

lite.connect(TIMELINE_DB); lite.connect(RESP_DB)

from seissuite.database import create_database, response_database

    