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

Please consult the manual for more information:
    http://bespoke.lbl.gov/weblogo/manual.html
(Also located in the weblogo/weblogo_htdocs subdirectory.)

-----------------------------
Distribution and Modification
-----------------------------

This package is distributed under the new BSD Open Source License. 
Please see the LICENSE.txt file for details on copyright and licensing.


