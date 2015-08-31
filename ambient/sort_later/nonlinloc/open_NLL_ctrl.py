# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 12:34:36 2015

@author: boland
"""

fname = 'NLL_control.in'

var_lines = []
with open(fname, 'rb') as f:
    for line in f:
        if line[0] is not '#' and line != '\n':
            # use the -2 index as a simple way to remove the '\n' at the end 
            # of every line!
            var_lines.append(line[:-1])
        
            


with open('new_control', 'a+') as f2:
    for line in var_lines: 
        print line
        f2.writelines(line + '\n')            
