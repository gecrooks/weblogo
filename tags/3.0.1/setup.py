#!/usr/bin/env python

import sys

from distutils.core import setup
from distutils.core import Extension
from distutils.command.build import build
from distutils.command.install_data import install_data

# Supress warning that distutils generates for the install_requires option
import warnings
warnings.simplefilter('ignore', UserWarning, lineno =236)

# check dependancies
if not hasattr(sys, 'version_info') or sys.version_info < (2,3,0,'final'):
    raise SystemExit,  \
        "Dependancy error: WebLogo requires Python 2.3 or later."
 
 
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
        
        scripts = [ 'weblogo', ],
        packages  = [ 'weblogolib',],
        data_files = ['weblogolib/htdocs/*.*','weblogolib/template.eps'],
        install_requires=['numpy', 'corebio'],        
        
        cmdclass= {"install_data" : _install_data},
    )


# Python 2.3 compatability 
# Rework the install_data command to act like the package_data distutils
# command included with python 2.4.
# Adapted from biopython, which was adapted from mxtexttools
class _install_data(install_data):
    def finalize_options(self):
        if self.install_dir is None:
            installobj = self.distribution.get_command_obj('install')
            # Use install_lib rather than install_platlib because we are
            # currently a pure python distribution (No c extensions.)
            self.install_dir = installobj.install_lib 
            #print installobj.install_lib 
        install_data.finalize_options(self)

    def run (self):
        import glob
        import os
        if not self.dry_run:
            self.mkpath(self.install_dir)
        data_files = self.get_inputs()
        for entry in data_files:
            if type(entry) is not type(""):
                raise ValueError("data_files must be strings")
            # Unix- to platform-convention conversion
            entry = os.sep.join(entry.split("/"))
            filenames = glob.glob(entry)
            for filename in filenames:
                dst = os.path.join(self.install_dir, filename)
                dstdir = os.path.split(dst)[0]
                if not self.dry_run:
                    self.mkpath(dstdir)
                    outfile = self.copy_file(filename, dst)[0]
                else:
                    outfile = dst
                self.outfiles.append(outfile)
  
if __name__ == '__main__' :
    main()
    
 