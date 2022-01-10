.. WebLogo documentation master file, created by
   sphinx-quickstart on Sat Nov  3 20:56:21 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

WebLogo 3: Sequence Logos Redrawn
===================================

WebLogo is software designed to make the generation of sequence logos easy and painless. A sequence logo is a graphical representation of an amino acid or nucleic acid multiple sequence alignment. Each logo consists of stacks of symbols, one stack for each position in the sequence. The overall height of the stack indicates the sequence conservation at that position, while the height of symbols within the stack indicates the relative frequency of each amino or nucleic acid at that position. In general, a sequence logo provides a richer and more precise description of, for example, a binding site, than would a consensus sequence.

WebLogo features a web interface (http://weblogo.threeplusone.com), and a command line interface provides more options
and control (http://weblogo.threeplusone.com/manual.html#CLI). These pages document the API.

The main WebLogo webserver is located at http://weblogo.threeplusone.com

Please consult the manual for installation instructions and more information:
(Also located in the weblogolib/htdocs subdirectory.)
http://weblogo.threeplusone.com/manual.html

For help on the command line interface run

.. code-block:: console

    weblogo --help

To build a simple logo run

.. code-block:: console

    weblogo  < cap.fa > logo0.eps

To run as a standalone webserver at localhost:8080

.. code-block:: console

    weblogo --serve


Distribution and Modification
-----------------------------

This package is distributed under the MIT Open Source License.
Please see the LICENSE.txt file for details on copyright and licensing.
The WebLogo source code can be downloaded from
https://github.com/WebLogo/weblogo

WebLogo requires Python 3.6 or 3.7. Generating logos in PDF or bitmap graphics formats
require that the ghostscript
program 'gs' be installed. Scalable Vector Graphics (SVG) format also requires
the program 'pdf2svg'.


.. contents:: :local:
.. currentmodule:: weblogo

.. toctree::
 	:maxdepth: 3
   	:caption: Contents:

   	intro
   	seq
   	seq_io
   	logo
  	formatting


* :ref:`genindex`
* :ref:`modindex`

