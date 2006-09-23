===============================================================================

                               WebLogo

                        http://code.google.com/p/weblogo/
===============================================================================

WebLogo is an application for generating sequence logos, graphical
representations of sequence conservation within amino acid or nucleic acid 
multiple sequence alignments. WebLogo makes extensive use of the CoreBio toolkit. 

WebLogo can be run on the command line, as a standalone webserver, as a CGI 
webapp, or it can be used as a python library.

For help on the command line interface run
    ./weblogo.py --help

To build a simple logo run
    ./weblogo.py  < cap.fa > logo0.eps
    
To run as a standalone webserver at localhost:8080 
    ./weblogo.py --server


------------
Installation
------------
Download the latest distribution from the project home page, http://code.google.com/p/weblogo/	
    > curl http://weblogo.googlecode.com/svn/dist/weblogo.tar.gz >weblogo.tar.gz    
    > tar -zxf weblogo.tar.gz
    > cd weblogo-*

WebLogo requires python version 2.3 (or above). The following command will quickly determine if python is available. 

    > python -V
    Python 2.4.2


Download python : http://www.python.org/download/

To install WebLogo from the main WebLogo directory:
    > python setup.py --help install           # Show installation options
    > python setup.py install                  # Install code.
    > python setup.py install --home ~/local   # Install code into alternative location



-----------------------------
Distribution and Modification
-----------------------------

This package is distributed under the new BSD Open Source License. 
Please see the LICENSE.txt file for details on copyright and licensing.


