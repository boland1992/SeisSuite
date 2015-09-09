#!/usr/bin/env python
"""
@author: Benjamin Boland
"""

from distutils.core import setup

setup(name='seissuite',
	version='0.1.3',
	description='Python Tools for Ambient Noise Seismology',
	author='Benjamin Boland',
	author_email='bolandb@student.unimelb.edu.au',
	url='tbc',
	packages=['seissuite',
                  'seissuite/misc',
                  'seissuite/azimuth',
                  'seissuite/ant', 
                  'seissuite/spacing',
                  'seissuite/spectrum',
                  'seissuite/trigger',
                  'seissuite/database',
                  'seissuite/response'
                  ],
      long_description=open('README.txt').read(),
	#package_dir = {'waveloc' : 'PyProgs', 'waveloc_examples' : 'examples'},
	#scripts=['scripts/grid2hdf5', 'scripts/pyms2sac', 'scripts/make_SDS'],
        requires=[
		'numpy(>=1.6.1)',
		'obspy.core(>=0.7.1)',
		'h5py(>=2.0.0)',
           'scipy(>=0.13.3)',
           'matplotlib(>=1.3.1)',
           'pyshp(>=1.1.7)',
           'pyproj(>=1.8.9)', 
           'pyPdf(>=1.13)'
	],
	classifiers=[
		'Environment :: Console',
		'Intended Audience :: Science/Research',
		'License :: OSI Approved :: GNU General Public License (GPL)',
		'Programming Language :: Python :: 2.7',
		'Topic :: Scientific/Engineering',
	], 
)
