#!/usr/bin/env python

import sys
from setuptools import setup, find_packages

# Suppress warning that distutils generates for the install_requires option
#import warnings

#warnings.simplefilter('ignore', UserWarning, lineno=236)


from weblogolib import __version__


def main():
    long_description = open("README.txt").read()

    setup(
            name="weblogo",
            version=__version__,
            description="WebLogo3 : Sequence Logos Redrawn",
            long_description=long_description,
            license = 'BSD',
            maintainer="Gavin Crooks",
            maintainer_email="gec@threeplusone.com",
            url="https://github.com/WebLogo/weblogo",
            download_url='https://github.com/WebLogo/weblogo/archive/%s.zip' % __version__,
            classifiers=[
                'Development Status :: 5 - Production/Stable',
                'Intended Audience :: Science/Research',
                'License :: OSI Approved :: BSD License',
                'Topic :: Scientific/Engineering :: Bio-Informatics',
                'Programming Language :: Python',
                'Natural Language :: English',
                'Operating System :: OS Independent',
                'Topic :: Software Development :: Libraries',
                'Topic :: Software Development :: Libraries :: Python Modules',
                
                # Specify the Python versions you support here. In particular, ensure
                # that you indicate whether you support Python 2, Python 3 or both.
                'Programming Language :: Python :: 2',
                'Programming Language :: Python :: 2.6',
                'Programming Language :: Python :: 2.7',
                'Programming Language :: Python :: 3',
                'Programming Language :: Python :: 3.2',
                'Programming Language :: Python :: 3.3',
                'Programming Language :: Python :: 3.4',                
                'Programming Language :: Python :: 3.5',                
                
            ],

            scripts=['weblogo', 'transformseq'],
            
            packages=[
                'corebio',
                'corebio.seq_io',
                'corebio.seq_io._nexus',
                'corebio.utils',
                'weblogolib',
            ],

            package_data={
                'weblogolib': ['htdocs/*.*', 'htdocs/img/*.*', 'htdocs/examples/*.*', 'template.eps'],
                'corebio': ['data/*.*']
            },
            
            install_requires=['numpy'],

    )


if __name__ == '__main__':
    main()
