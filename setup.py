#!/usr/bin/env python

from setuptools import setup


long_description = """
WebLogo (https://github.com/WebLogo/weblogo) is a tool for creating sequence
logos from biological sequence alignments.  It can be run on the command line,
as a standalone webserver, as a CGI webapp, or as a python library.

The main WebLogo webserver is located at http://weblogo.threeplusone.com

Please consult the manual for installation instructions and more information:
./weblogo/htdocs/manual.html

(Also located at http://weblogo.threeplusone.com/manual.html.)
"""

setup(
    name="weblogo",
    python_requires='>=3.6',

    install_requires=['numpy', 'scipy', 'setuptools'],

    use_scm_version=True,
    setup_requires=['setuptools_scm'],


    description="WebLogo3 : Sequence Logos Redrawn",
    long_description=long_description,
    license='BSD',
    maintainer="Gavin Crooks",
    maintainer_email="gec@threeplusone.com",
    url="https://github.com/WebLogo/weblogo",

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
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',

    ],

    packages=[
        'weblogo',
        'weblogo.seq_io',
        'weblogo.seq_io._nexus',
        'weblogo.utils',
    ],

    package_data={
        'weblogo': ['htdocs/*.*', 'htdocs/img/*.*', 'htdocs/examples/*.*',
                    'template.eps', 'data/*.*']
    },



    entry_points={
        'console_scripts': [
            'weblogo = weblogo._cli:main',
            'transformseq = weblogo._transformseq:main'
        ],
    },
)
