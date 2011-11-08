#!/usr/bin/env python

import sys

from distutils.core import setup
from distutils.core import Extension
from distutils.command.build import build
from distutils.command.install_data import install_data

# Supress warning that distutils generates for the install_requires option
import warnings
warnings.simplefilter('ignore', UserWarning, lineno =236)

# check dependencies
if not hasattr(sys, 'version_info') or sys.version_info < (2,5,0,'final'):
    raise SystemExit,  \
        "Dependency error: WebLogo requires Python 2.5 or later."
 
 
from weblogolib import __version__

def main() :     
    long_description = open("README.txt").read()

      
    setup(
        name             =  "weblogo",
        version          =  __version__,
        description      = "WebLogo3 : Sequence Logos Redrawn",
        long_description  = long_description,
        maintainer       = "Gavin Crooks",
        maintainer_email = "gec@threeplusone.com",
        url              = "http://code.google.com/p/weblogo/",    
        
        download_url     = 'http://weblogo.googlecode.com/files/weblogo-%s.tar.gz' % __version__ ,
        classifiers      =[
            'Development Status :: 5 - Production/Stable',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: BSD License',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Programming Language :: Python',
            'Natural Language :: English',
            'Operating System :: OS Independent',
            'Topic :: Software Development :: Libraries',
            'Topic :: Software Development :: Libraries :: Python Modules',
            ],
        
        scripts = [ 'weblogo', 'transformseq'],
	    packages  = [ 
	        'corebio', 
	        'corebio.db', 
	        'corebio.secstruc',
	        'corebio.seq_io', 
	        'corebio.seq_io._nexus', 
	        'corebio.ssearch_io',         
	        'corebio.utils',
			'weblogolib',
	        ],

		package_data={
		'weblogolib': ['htdocs/*.*','htdocs/img/*.*','htdocs/examples/*.*','template.eps'],
		'corebio': ['data/*.*']
			},
        requires=['numpy'],        
        
    )


if __name__ == '__main__' :
    main()
    
 
