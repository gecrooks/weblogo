
WebLogo (http://code.google.com/p/weblogo/) is a tool for creating sequence 
logos from biological sequence alignments.  It can be run on the command line,
as a standalone webserver, as a CGI webapp, or as a python library.

The main WebLogo webserver is located at http://bespoke.lbl.gov/weblogo/

Please consult the manual for installation instructions and more information:
(Also located in the weblogo/weblogo_htdocs subdirectory.)

    http://bespoke.lbl.gov/weblogo/manual.html

For help on the command line interface run
    ./weblogo.py --help

To build a simple logo run
    ./weblogo.py  < cap.fa > logo0.eps
    
To run as a standalone webserver at localhost:8080 
    ./weblogo.py --server

To create a logo in python code:
    >>> from weblogolib import *
    >>> fin = open('cap.fa')
    >>> seqs = read_seq_data(fin) 
    >>> data = LogoData.from_seqs(seqs)
    >>> options = LogoOptions()
    >>> options.title = "A Logo Title"
    >>> format = LogoFormat(data, options)
    >>> fout = open('cap.eps', 'w') 
    >>> eps_formatter( data, format, fout)


-- Distribution and Modification --
This package is distributed under the new BSD Open Source License. 
Please see the LICENSE.txt file for details on copyright and licensing.
The WebLogo source code can be downloaded from http://code.google.com/p/weblogo/
WebLogo requires Python 2.3, 2.4 or 2.5, the corebio python toolkit for computational biology (http://code.google.com/p/corebio), and the python array 
package 'numpy' (http://www.scipy.org/Download)
