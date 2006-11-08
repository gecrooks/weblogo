#!/usr/bin/env python

# -------------------------------- WebLogo --------------------------------

#  Copyright (c) 2003-2004 The Regents of the University of California.
#  Copyright (c) 2005 Gavin E. Crooks
#  Copyright (c) 2006, The Regents of the University of California, through 
#  Lawrence Berkeley National Laboratory (subject to receipt of any required
#  approvals from the U.S. Dept. of Energy).  All rights reserved.

#  This software is distributed under the new BSD Open Source License.
#  <http://www.opensource.org/licenses/bsd-license.html>
#
#  Redistribution and use in source and binary forms, with or without 
#  modification, are permitted provided that the following conditions are met: 
#
#  (1) Redistributions of source code must retain the above copyright notice, 
#  this list of conditions and the following disclaimer. 
#
#  (2) Redistributions in binary form must reproduce the above copyright 
#  notice, this list of conditions and the following disclaimer in the 
#  documentation and or other materials provided with the distribution. 
#
#  (3) Neither the name of the University of California, Lawrence Berkeley 
#  National Laboratory, U.S. Dept. of Energy nor the names of its contributors 
#  may be used to endorse or promote products derived from this software 
#  without specific prior written permission. 
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
#  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
#  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
#  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
#  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
#  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
#  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
#  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
#  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
#  POSSIBILITY OF SUCH DAMAGE. 

""" WebLogo is a tool for creating sequence logos from biological sequence
alignments.  It can be run on the command line, as a standalone webserver, as a
CGI webapp, or as a python library.

For help on the command line interface run
    ./weblogo.py --help

To build a simple logo run
    ./weblogo.py  < cap.fa > logo0.eps
    
To run as a standalone webserver at localhost:8080 
    ./weblogo.py --serve

To create a logo in python code:
    >>> seqs = seq_io.read( fin) 
    >>> data = LogoData(data)
    >>> options = LogoOptions()
    >>> options.title = "A Logo Title"
    >>> format = LogoFormat(data, options) 
    >>> logo = format( data, format)

WebLogo makes extensive use of the corebio python toolkit for computational biology.  (http://code.google.com/p/corebio)

"""

#TODO: Split out Motif class
#TODO: Add LogoValueError exception class rather than current tuple hack.
#TODO: Refactor to use Ghostscript class.


import sys
import copy
import os
from datetime import datetime
from StringIO import StringIO

# We have our own private copy of corebio stored as weblogo._corebio.
# First we import our copy, then we mangle sys.modules so that package
# has correct name
import weblogo._corebio
sys.modules['corebio'] = sys.modules['weblogo._corebio']

from  corebio.data import rna_letters, dna_letters, amino_acid_letters
import random

# python2.3 compatability
from corebio._future import Template
from corebio._future.subprocess import *

from math import log, sqrt

# Avoid 'from numarray import *' since numarray has lots of names defined
from numarray import array, asarray, Float32, Float64, ones, zeros, Int32,all,any, shape
import numarray as na

from corebio.utils.deoptparse import DeOptionParser
from optparse import OptionGroup

from weblogo.color import *
from weblogo.colorscheme import *
from corebio.seq import Alphabet, Seq, SeqList
from corebio import seq_io
from corebio.utils import isfloat, find_command
from corebio.moremath import *
from corebio.data import amino_acid_composition
from corebio.seq import unambiguous_rna_alphabet, unambiguous_dna_alphabet, unambiguous_protein_alphabet


# ------ META DATA ------

description  = "Create sequence logos from biological sequence alignments." 

__version__ = "3.0b9"

# These keywords are subsituted by subversion.
# The date and revision will only  tell the truth after a branch or tag,
# since different files in trunk will have been changed at different times
release_date ="$Date: 2006-09-20 21:55:15 -0700 (Wed, 20 Sep 2006) $".split()[1]
release_build = "$Revision: 152 $".split()[1]
release_description = "WebLogo %s (%s)" % (__version__,  release_date)





def cgi(htdocs_directory) :
    import weblogo._cgi
    weblogo._cgi.main(htdocs_directory)
        
class Ghostscript(object) :
    formats = ('png', 'pdf', 'jpeg')
    
    def __init__(self, path=None) :
        command = find_command('gs', path=path) 
        self.command = command
        
    
    def version(self) :
        args = [self.command, '--version']
        try :
            p = Popen(args, stdout=PIPE)
            (out,err) = p.communicate() 
        except OSError :
            raise RuntimeError("Cannot communicate with ghostscript.")  
        return out.strip()
       
    def convert(self, format,fin, fout,  width,  height, resolution=300) :
        device_map = { 'png':'png16m',  'pdf':'pdfwrite', 'jpeg':'jpeg'}
       
        try :
            device = device_map[format]
        except KeyError:
            raise ValueError("Unsupported format.")
        
        args = [self.command, 
            "-sDEVICE=%s" % device, 
            "-dPDFSETTINGS=/printer",
            "-q",   # Quite: Do not dump messages to stdout.
            "-sOutputFile=-", # Stdout
            "-dDEVICEWIDTHPOINTS=%s" % str(width),  
            "-dDEVICEHEIGHTPOINTS=%s" % str(height),  
            "-dSAFER",  # For added security
            "-dNOPAUSE",]
            
        if device != 'pdf' :
            args.append("-r%s" % str(resolution) ) 
            if resolution < 300 : # Antialias if resolution is Less than 300 DPI
                args.append("-dGraphicsAlphaBits=4")
                args.append("-dTextAlphaBits=4")
                args.append("-dAlignToPixels=0")
        
        args.append("-")  # Read from stdin. Must be last argument.
        
        error_msg = "Unrecoverable error : Ghostscript conversion failed " \
                    "(Invalid postscript?). %s" % " ".join(args) 

        try :
            p = Popen(args, stdin=fin, stdout = fout) 
            (out,err) = p.communicate() 
        except OSError :
            raise RuntimeError(error_msg)
    
        if p.returncode != 0 : raise RuntimeError(error_msg)
# end class Ghostscript


aa_composition = [ amino_acid_composition[_k] for _k in unambiguous_protein_alphabet]



# ------  DATA ------

classic = ColorScheme([
    ColorGroup("G",  "orange" ),
    ColorGroup("TU", "red"),
    ColorGroup("C",  "blue"),
    ColorGroup("A",  "green")
    ] )
    
std_color_schemes = {"auto":            None,   # Depends on sequence type
                     "monochrome":      monochrome,
                     "base pairing":    base_pairing,
                     "classic":         classic,
                     "hydrophobicity" : hydrophobicity,  
                     "chemistry" :      chemistry,
                     "charge" :         charge,
                     }#

default_color_schemes = {
    unambiguous_protein_alphabet: hydrophobicity, 
    unambiguous_rna_alphabet: base_pairing,
    unambiguous_dna_alphabet: base_pairing 
}


std_units = {
    "bits"    : 1./log(2),
    "nats"    : 1.,
    "digits"  : 1./log(10),
    "kT"      : 1.,
    "kJ/mol"  : 8.314472 *298.15 /1000.,
    "kcal/mol": 1.987 *298.15  /1000.,
    "probability" : None,
}

class LogoSize(object) :
    def __init__(self, stack_width, stack_height) :
         self.stack_width = stack_width
         self.stack_height = stack_height


# The base stack width is set equal to 9pt Courier. 
# (Courier has a width equal to 3/5 of the point size.)
# Check that can get 80 characters in journal page @small
# 40 chacaters in a journal column
std_sizes = {
    "small" : LogoSize( stack_width = 5.4, stack_height = 5.4*1*5), 
    "medium" : LogoSize( stack_width = 5.4*2, stack_height = 5.4*2*5),
    "large"  : LogoSize( stack_width = 5.4*3, stack_height = 5.4*3*5),    
}
            

            
std_alphabets = {
    'protein': unambiguous_protein_alphabet, 
    'rna': unambiguous_rna_alphabet,
    'dna': unambiguous_dna_alphabet}

std_percentCG = {
    'H. sapiens'    : 40.,
    'E. coli'       : 50.5,
    'S. cerevisiae' : 38.,
    'C. elegans'    : 36.,
    'D. melanogaster': 43.,
    'M. musculus'   :  42.,
    'T. thermophilus' : 69.4,
}

# Thermus thermophilus: Henne A, Bruggemann H, Raasch C, Wiezer A, Hartsch T,
# Liesegang H, Johann A, Lienard T, Gohl O, Martinez-Arias R, Jacobi C, 
# Starkuviene V, Schlenczeck S, Dencker S, Huber R, Klenk HP, Kramer W, 
# Merkl R, Gottschalk G, Fritz HJ: The genome sequence of the extreme 
# thermophile Thermus thermophilus.
# Nat Biotechnol 2004, 22:547-53
            

  
class LogoOptions(object) :
    """ A container for all logo formating options. Not all of these
    are directly accesible through the CLI or web interfaces.    
    """
    # TODO: Make sure all of these options are actually used.
    # TODO: Plot type, units and the default yaxis label
    # TODO: Add docstring documentation for all of these options.
    def __init__(self) :
        self.creator_text = release_description,

        self.logo_title = ""
        self.logo_label = ""
        self.stacks_per_line = 40
        #self.total_stacks = None #FIXME
        
        self.unit_name = "bits"
     
        self.show_yaxis = True
        self.yaxis_label = None         # default depends on other settings. See LogoFormat
        self.yaxis_tic_interval = 1.
        self.yaxis_minor_tic_ratio = 5
        self.yaxis_scale = None

        self.show_xaxis = True
        self.xaxis_label = ""
        self.xaxis_tic_interval =1
        self.rotate_numbers = False
        self.number_interval = 5            
        self.show_ends = False
          
        self.show_fineprint = True
        self.fineprint =  "WebLogo "+__version__ 
    
        self.show_boxes = False
        self.shrink_fraction = 0.5
        self.outline_linewidth  = 0.2           # FIXME: Used?
  
        self.show_errorbars = True
        self.errorbar_fraction = 0.0
        self.errorbar_width_fraction = 0.25      
        self.errorbar_gray = 0.5
 
        self.resolution  = 96.     # Dots per inch
           
        self.default_color = Color.by_name("black")    
        self.color_scheme = None
        self.show_color_key = False
        
        self.alphabet = None    
                
        self.debug = False
        
        self.logo_margin = 2
        self.stroke_width = 0.5
        self.tic_length = 5
        
        self.size = std_sizes["medium"]        
        
        self.stack_margin = 0.5
        self.pad_right = False  

        self.small_fontsize = 6
        self.fontsize = 10
        self.title_fontsize = 12
        self.number_fontsize = 8

        self.text_font    = "ArialMT"
        self.logo_font    = "Arial-BoldMT"
        self.title_font   = "ArialMT"

        self.first_index = 1        
        self.logo_start = None      
        self.logo_end=None          

        # Scale width of characters proportional to gaps
        self.scale_width = True

# End class LogoOptions

        

class LogoFormat(LogoOptions) :     
    """ Specifies the format of the logo. Requires a LogoData and LogoOptions 
    objects.
    
    >>> data = LogoData(data)
    >>> options = LogoOptions()
    >>> options.title = "A Logo Title"
    >>> format = LogoFormat(data, options) 
    """
    
    def __init__(self, data, options= LogoOptions() ) :
        LogoOptions.__init__(self)
        self.__dict__.update(options.__dict__)

        self.data = data
        
        self.alphabet = data.seqs.alphabet
         
        self.show_title = False
        self.show_xaxis_label = False
        self.yaxis_minor_tic_interval = None
        self.lines_per_logo       = None
        self.char_width       = None
        self.line_margin_left = None
        self.line_margin_right    = None
        self.line_margin_bottom = None
        self.line_margin_top = None
        self.title_height = None
        self.xaxis_label_height = None
        self.line_height = None
        self.line_width = None
        self.logo_height = None
        self.logo_width = None
        self.creation_date = None
        self.end_type = None

        if self.stacks_per_line< 1 :
            raise ValueError (("stacks_per_line", "Stacks per line should be greater than zero.") )
        if self.size.stack_height<=0.0 : 
            raise ValueError (("stack_height", "Stack height must be greater than zero.") )
        if self.small_fontsize <= 0 or self.fontsize <=0 or self.title_fontsize<=0 :
            raise ValueError("Font sizes must be positive.")
        if self.errorbar_fraction<0.0 or self.errorbar_fraction>1.0 :
            raise ValueError(
                "The visible fraction of the error bar must be between zero and one.")
        if self.yaxis_tic_interval<=0.0 :
            raise ValueError(
                ('yaxis_tic_interval',"The yaxis tic interval cannot be negative.") )
        if self.size.stack_width <= 0.0 :
            raise ValueError("The width of a stack should be a positive number.")
        if self.yaxis_minor_tic_interval and self.yaxis_minor_tic_interval<=0.0 : 
            raise ValueError("Distances cannot be negative.")
        if self.xaxis_tic_interval<=0 :
            raise ValueError("Tic interval must be greater than zero.")
        if self.number_interval<=0 :
            raise ValueError("Invalid interval between numbers.")
        if self.shrink_fraction<0.0 or self.shrink_fraction>1.0 :
            raise ValueError("Invalid shrink fraction.")
        if self.stack_margin<=0.0 : 
            raise ValueError("Invalid stack margin."  )
        if self.logo_margin<=0.0 : 
            raise ValueError("Invalid logo margin."  )
        if self.stroke_width<=0.0 : 
            raise ValueError("Invalid stroke width.")  
        if self.tic_length<=0.0 : 
            raise ValueError("Invalid tic length.") 

        # FIXME: More validation
         
        # Inclusive upper and lower bounds
        # FIXME: Validate here. Move from eps_formatter        
        if self.logo_start is  None: self.logo_start = self.first_index
        if self.logo_end is  None : self.logo_end = len(data.seqs[0]) + self.first_index -1 
        self.total_stacks = self.logo_end - self.logo_start +1

        if self.logo_start - self.first_index <0 :
            raise ValueError(
                ('logo_range',"Logo range extends before start of available sequence.") )
        if self.logo_end - self.first_index  >= len(data.seqs[0])  : 
            raise ValueError(
                ('logo_range',"Logo range extends beyond end of available sequence.") )
    
        if self.logo_title      : self.show_title = True
        if not self.fineprint   : self.show_fineprint = False
        if self.xaxis_label     : self.show_xaxis_label = True

        if self.yaxis_label is None : 
            self.yaxis_label = self.unit_name
        if self.yaxis_label : 
            self.show_yaxis_label = True
        else :
            self.show_yaxis_label = False
            self.show_ends = False
              
        if not self.yaxis_scale :
            conversion_factor = std_units[self.unit_name]
            if conversion_factor :
                self.yaxis_scale=log(len(self.alphabet))*conversion_factor
            else :
                self.yaxis_scale=1.0    # probability units

        if self.yaxis_scale<=0.0 : 
            raise ValueError(('yaxis_scale', "Invalid yaxis scale"))

        if self.yaxis_tic_interval >= self.yaxis_scale:
            self.yaxis_tic_interval /= 2. 

        self.yaxis_minor_tic_interval \
            = float(self.yaxis_tic_interval)/self.yaxis_minor_tic_ratio
                      
        if self.color_scheme is None :
            if self.alphabet in default_color_schemes :
                self.color_scheme = default_color_schemes[self.alphabet]
            else :
                self.color_scheme = monochrome

        self.lines_per_logo = 1+ ( (self.total_stacks-1) / self.stacks_per_line)
    
        if self.lines_per_logo==1 and not self.pad_right:
            self.stacks_per_line = min(self.stacks_per_line, self.total_stacks)

        self.char_width = self.size.stack_width - 2* self.stack_margin
    
    
        if self.show_yaxis :
            self.line_margin_left = self.fontsize * 3.0
        else :
            self.line_margin_left = 0

        if self.show_ends :
            self.line_margin_right = self.fontsize *1.5 
        else :
            self.line_margin_right = self.fontsize             

        if self.show_xaxis :
            if self.rotate_numbers :
                self.line_margin_bottom = self.number_fontsize *2.5
            else:
                self.line_margin_bottom = self.number_fontsize *1.5
        else :
            self.line_margin_bottom = 4

        self.line_margin_top = 4
    
        if self.show_title :
            self.title_height = self.title_fontsize 
        else :
            self.title_height = 0

        self.xaxis_label_height =0.
        if self.show_xaxis_label :
            self.xaxis_label_height += self.fontsize
        if self.show_fineprint :
            self.xaxis_label_height += self.small_fontsize

        self.line_height = (self.size.stack_height + self.line_margin_top +    
                            self.line_margin_bottom )
        self.line_width  = (self.size.stack_width*self.stacks_per_line + 
                            self.line_margin_left + self.line_margin_right )

        self.logo_height = 2*self.logo_margin + self.title_height \
                + self.xaxis_label_height + self.line_height*self.lines_per_logo
        self.logo_width = 2*self.logo_margin + self.line_width 


        self.creation_date = datetime.now().isoformat(' ') 

        end_type = '-'
        end_types = {
            unambiguous_protein_alphabet: 'p', 
            unambiguous_rna_alphabet: '-',
            unambiguous_dna_alphabet: 'd' 
        }
        if self.show_ends and self.alphabet in end_types:
            end_type = end_types[self.alphabet]
        self.end_type = end_type
    # End __init__
# End class LogoFormat



# ------ Logo Formaters ------
# Each formatter is a function f(LogoData, LogoFormat, output file).
# that draws a represntation of the logo into the given file.
# The main graphical formatter is eps_formatter. A mapping 'formatters'
# containing all available formatters is located after the formatter
# definitions. 

def pdf_formatter(data, format, fout) :
    """ Generate a logo in PDF format."""
    device = "pdfwrite"
    args = ["gs", 
        "-sDEVICE=%s" % device, 
        "-dPDFSETTINGS=/printer",
        "-q",   # Quite: Do not dump messages to stdout.
        "-sOutputFile=-",
        "-dDEVICEWIDTHPOINTS=%s" % str(format.logo_width),  
        "-dDEVICEHEIGHTPOINTS=%s" % str(format.logo_height),  
        "-dSAFER",  # For added security
        "-dNOPAUSE",
        "-", # Read from stdin. Must be last argument.
        ]

    f = StringIO()
    eps_formatter(data, format, f)
    eps = f.getvalue()

    p = Popen(args, stdin=PIPE,stdout = PIPE) 
    (out,err) = p.communicate(eps) 
    if p.returncode != 0 or len(out)==0:
        raise RuntimeError("Unrecoverable seqlogo error : Ghostscript conversion failed (Invalid postscript?). Sorry.") 
    print >>fout, out

def _bitmap_formatter(data, format, fout, device) :

    args = ["gs", 
        "-sDEVICE=%s" % device, 
        "-dPDFSETTINGS=/printer",
        "-q",
        "-r%s" % str(format.resolution), 
        "-sOutputFile=-",
        "-dDEVICEWIDTHPOINTS=%s" % str(format.logo_width),  
        "-dDEVICEHEIGHTPOINTS=%s" % str(format.logo_height),  
        "-dSAFER",
        "-dNOPAUSE",]
        
    if format.resolution<300 : # Antialias if resolution is Less than 300 DPI
        args.append("-dGraphicsAlphaBits=4")
        args.append("-dTextAlphaBits=4")
        args.append("-dAlignToPixels=0")
        
    args.append("-")      

    f = StringIO()
    eps_formatter(data, format, f)
    eps = f.getvalue()

    p = Popen(args, stdin=PIPE,stdout = PIPE) 
    try :
        (out,err) = p.communicate(eps) 
    except OSError, err :
        raise RuntimeError("Unrecoverable seqlogo error : Ghostscript conversion failed (Invalid postscript?). Sorry. %s" % args )
    
    if p.returncode != 0 or len(out)==0:
        raise RuntimeError("","Unrecoverable seqlogo error : Ghostscript conversion failed (Invalid postscript?). Sorry." )


    print >>fout, out


def jpeg_formatter(data, format, fout) : 
    """ Generate a logo in JPEG format."""
    _bitmap_formatter(data, format, fout, device="jpeg")


def png_formatter(data, format, fout) : 
    """ Generate a logo in PNG format."""
    _bitmap_formatter(data, format, fout, device="png16m")


def png_print_formatter(data, format, fout) : 
    """ Generate a logo in PNG format with print quality (600 DPI) resolution."""
    format.resolution = 600
    _bitmap_formatter(data, format, fout, device="png16m")


def txt_formatter( logodata, format, fout) :
    """ Create a text representation of the logo data. 
    Status: Alpha
    """
    # TODO: write me
    print >>fout, logodata
    seq_from = format.logo_start- format.first_index
    seq_to = format.logo_end - format.first_index +1
    #from corebio.distributions import Dirichlet 
    #from numarray import asarray, Float64
    
    # seq_index : zero based index into sequence data
    # logo_index : User visible coordinate, first_index based
    # stack_index : zero based index of visible stacks
    for seq_index in range(seq_from, seq_to) :
        logo_index = seq_index + format.first_index 
        #stack_index = seq_index - seq_from
        
    
        print >>fout, '>', logo_index, 
    
        print >>fout, 'E:', logodata.entropy[seq_index], ' A:', logodata.counts[seq_index].sum(),
        
        for c in logodata.counts[seq_index] :
            print >>fout, c,
        print >>fout, ''         

     #   print >>fout, Dirichlet( asarray(logodata.counts[seq_index], Float64)+0.25).mean_entropy()
   

    
def eps_formatter( logodata, format, fout) :
    """ Generate a logo in Encapsulated Postscript (EPS)"""
    
    subsitutions = {}
    from_format =[
        "creation_date",    "logo_width",           "logo_height",      
        "lines_per_logo",   "line_width",           "line_height",
        "line_margin_right","line_margin_left",     "line_margin_bottom",
        "line_margin_top",  "title_height",         "xaxis_label_height",
        "creator_text",     "logo_title",           "logo_margin",
        "stroke_width",     "tic_length",           
        "stacks_per_line",  "stack_margin",
        "yaxis_label",      "yaxis_tic_interval",   "yaxis_minor_tic_interval",
        "xaxis_label",      "xaxis_tic_interval",   "number_interval",
        "fineprint",        "shrink_fraction",      "errorbar_fraction",
        "errorbar_width_fraction",
        "errorbar_gray",    "small_fontsize",       "fontsize",
        "title_fontsize",   "number_fontsize",      "text_font",
        "logo_font",        "title_font",           "outline_linewidth",
        "logo_label",       "yaxis_scale",          "end_type",
        "debug",            "show_title",           "show_xaxis",
        "show_xaxis_label", "show_yaxis",           "show_yaxis_label",
        "show_boxes",       "show_errorbars",       "show_fineprint",
        "rotate_numbers",   "show_ends",            
        ]
   
    for s in from_format :
        subsitutions[s] = getattr(format,s)


    from_format_size = ["stack_height", "stack_width"]
    for s in from_format_size :
        subsitutions[s] = getattr(format.size,s)

    subsitutions["shrink"] = str(format.show_boxes).lower()


    # --------- COLORS --------------
    def format_color(color):
        return  " ".join( ("[",str(color.red) , str(color.green), 
            str(color.blue), "]"))  

    subsitutions["default_color"] = format_color(format.default_color)

    colors = []  
    for group in format.color_scheme.groups :
        cf = format_color(group.color)
        for s in group.symbols :
            colors.append( "  ("+s+") " + cf )
    subsitutions["color_dict"] = "\n".join(colors)
        
    data = []
    
    # Unit conversion. 'None' for probability units
    conv_factor = std_units[format.unit_name]
    
    data.append("StartLine")
    

    seq_from = format.logo_start- format.first_index
    seq_to = format.logo_end - format.first_index +1

    # seq_index : zero based index into sequence data
    # logo_index : User visible coordinate, first_index based
    # stack_index : zero based index of visible stacks
    for seq_index in range(seq_from, seq_to) :
        logo_index = seq_index + format.first_index 
        stack_index = seq_index - seq_from
        
        if stack_index!=0 and (stack_index % format.stacks_per_line) ==0 :
            data.append("")
            data.append("EndLine")
            data.append("StartLine")
            data.append("")
        
        if logo_index % format.number_interval == 0 : 
            data.append("(%d) StartStack" % logo_index)
        else :            
            data.append("() StartStack" )

        if conv_factor: 
            stack_height = logodata.entropy[seq_index] * std_units[format.unit_name]
        else :
            stack_height = 1.0 # Probability

        if logodata.entropy_interval is not None and conv_factor:
            # Draw Error bars
            low, high = logodata.entropy_interval[seq_index]
            center = logodata.entropy[seq_index]


            down = (center - low) * conv_factor
            up   = (high - center) * conv_factor
            data.append(" %f %f %f DrawErrorbarFirst" % (down, up, stack_height) )
        
        s = zip(logodata.counts[seq_index], logodata.seqs.alphabet)
        def mycmp( c1, c2 ) :
            # Sort by frequency. If equal frequency then reverse alphabetic
            if c1[0] == c2[0] : return cmp(c2[1], c1[1])
            return cmp(c1[0], c2[0])
        
        s.sort(mycmp)

        C = float(sum(logodata.counts[seq_index])) 
        if C > 0.0 :
            fraction_width = 1.0
            if format.scale_width :
                fraction_width = C / len(logodata.seqs)
            # print >>sys.stderr, fraction_width
            for c in s:
                data.append(" %f %f (%s) ShowSymbol" % (fraction_width, c[0]*stack_height/C, c[1]) )

        # Draw error bar on top of logo. Replaced by DrawErrorbarFirst above.
        #if logodata.entropy_interval is not None and conv_factor:
        #    low, high = logodata.entropy_interval[seq_index]
        #    center = logodata.entropy[seq_index]


        #    down = (center - low) * conv_factor
        #    up   = (high - center) * conv_factor
        #    data.append(" %f %f DrawErrorbar" % (down, up) )
            
        data.append("EndStack")
        data.append("")
               
    data.append("EndLine")
    subsitutions["logo_data"] = "\n".join(data)  


    # Create and output logo
    logo = Template(_template).substitute(subsitutions)
    print >>fout, logo
 

# map between output format names and logo  
formatters = {
    'eps': eps_formatter, 
    'pdf': pdf_formatter,
    'png': png_formatter,
    'png_print' : png_print_formatter,
    'jpeg'  : jpeg_formatter,
    #'txt' : txt_formatter,          # FIXME: Implement
    }     
    
default_formatter = eps_formatter





def parse_prior(composition, alphabet, weight=None) :
    """ Parse a description of the expected monomer distribution of a sequence.
    
    Valid compositions:
    
    - None or 'none' :      No composition sepecified 
    - 'auto' or 'automatic': Use the typical average distribution
                            for proteins and an equiprobable distribution for
                            everything else.    
    - 'equiprobable' :      All monomers have the same probability.
    - a percentage, e.g. '45%' or a fraction '0.45':
                            The fraction of CG bases for nucleotide alphabets
    - a species name, e.g. 'E. coli', 'H. sapiens' :
                            Use the average CG percentage for the specie's      
                            genome.
    - An explicit distribution,  e.g. {'A':10, 'C':40, 'G':40, 'T':10}
    """
    if composition is None: return None
    composition = composition.strip()
    comp = composition.lower()
    
    if comp == 'none': return None
    
    if weight is None: weight = float(len(alphabet))
    
    if comp == 'equiprobable' :
        prior = weight * equiprobable_distribution(len(alphabet)) 
    elif comp == 'auto' or comp == 'automatic':
        if alphabet == unambiguous_protein_alphabet :
            prior =  weight * asarray(aa_composition, Float64)
        else :
            prior = weight * equiprobable_distribution(len(alphabet)) 
    
    elif comp in std_percentCG :
        prior = weight * base_distribution(std_percentCG[comp])

    elif comp[-1] == '%' :
        prior = weight * base_distribution( float(comp[:-1]))

    elif isfloat(comp) :
        prior = weight * base_distribution( float(comp)*100. )

    elif composition[0] == '{' and composition[-1] == '}' : 
        explicit = composition[1: -1]
        explicit = explicit.replace(',',' ').replace("'", ' ').replace('"',' ').replace(':', ' ').split()
        
        if len(explicit) != len(alphabet)*2 :
            print explicit
            raise ValueError("Explicit prior does not match length of alphabet")
        prior = - ones(len(alphabet), Float64) 
        try :
            for r in range(len(explicit)/2) :
                letter = explicit[r*2]
                index = alphabet.ord(letter)
                value = float(explicit[r*2 +1])
                prior[index] = value
        except ValueError :
            raise ValueError("Cannot parse explicit composition")
    
        if any(prior==-1.) :
            raise ValueError("Explicit prior does not match alphabet") 
        
    else : 
        raise ValueError("NO SUCH COMPOSITION")
    
    if len(prior) != len(alphabet) :
        raise ValueError("The sequence alphabet and composition are incompatible.")
    return prior

    
def base_distribution(percentCG) :
    A = (1. - (percentCG/100.))/2.
    C = (percentCG/100.)/2.
    G = (percentCG/100.)/2.
    T = (1. - (percentCG/100))/2.
    return asarray((A,C,G,T), Float64)   

def equiprobable_distribution( length) :
    return ones( shape = (length), type=Float64) /length   
  



def read_seq_data(fin, input_parser=seq_io.read, alphabet=None, cased=False, max_file_size=0,):

    # If max_file_size is set, or if fin==stdin (which is non-seekable), we
    # read the data and replace fin with a StringIO object. 
    if(max_file_size>0) :
        data = fin.read(max_file_size)
        more_data = fin.read(2)
        if more_data != "" :
            raise IOError("File too large %d"  % max_file_size) #FIXME: Better message
        fin = StringIO(data)
    elif fin == sys.stdin:
        fin = StringIO(fin.read())
        
    seqs = input_parser(fin)

    if seqs is None or len(seqs) ==0 :
        raise ValueError("Please provide a multiple sequence alignment")
    
    if cased :
        # Case is significant. Do not count lower case letters.
        for i,s in enumerate(seqs) :
            seqs[i] = s.mask()

    # Add alphabet to seqs.
    if alphabet :
        seqs.alphabet = alphabet 
    else :
        seqs.alphabet = which_alphabet(seqs)

    return seqs

 
#def guess_alphabet(seqs) :    
#    total_count = 0
#    amino_count = 0
#    rna_count   = 0
#    dna_count   = 0 
#    
#    for s in seqs:
#        for c in s:
#            total_count +=1
#            if c in amino_acid_letters: amino_count +=1
#            if c in dna_letters : dna_count +=1
#            if c in rna_letters: rna_count += 1
#    #print total_count, amino_count, rna_count, dna_count
#                
#    best = max(amino_count, dna_count, rna_count)
#    if best == dna_count : return unambiguous_dna_alphabet
#    if best == rna_count : return unambiguous_rna_alphabet
#    return unambiguous_protein_alphabet


#TODO: Move to seq_io? 
#   Would have to check that no symbol outside of full alphabet?
def which_alphabet(seqs) :
    """ Returns the most appropriate unambiguous protien, rna or dna alphabet
    for a Seq or SeqList.
    """
    alphabets = (unambiguous_protein_alphabet,
                unambiguous_rna_alphabet,
                unambiguous_dna_alphabet,
            )
    # Heuristic
    # Count occurances of each letter. Downweight longer alphabet.
    score = [1.0*asarray(seqs.tally(a)).sum()/sqrt(len(a)) for a in alphabets]
    #print score
    best = argmax(score) # Ties go to last in list.
    return alphabets[best]
    

    
    
class LogoData(object) :
    """The data needed to generate a sequence logo.
       
    - seqs  -- The actual sequence data. A corebio.seq.SeqList object 
    - counts  -- An array of character counts
    - entropy -- The relative entropy of each column
    - entropy_interval -- entropy confidence interval
     """
    # FIXME: Clean this up
    def __init__(self, seqs, prior= None):

        # --- VALIDATE DATA ---
        # check that at least one sequence of length at least 1 long
        if len(seqs)==0 or len(seqs[0]) ==0:
            raise ValueError("No sequence data found.")   
    
        # Check sequence lengths    
        seq_length = len(seqs[0])
        for i,s in enumerate(seqs) :
            #print i,s, len(s)
            if seq_length != len(s) :
                raise ValueError(( 'sequences',
                  "Sequence number %d differs in length from the previous sequences" % (i+1) ) )

        # FIXME: Check seqs.alphabet?


        # --- Build Logo ---
        
        if prior is not None:
            prior = array(prior, Float64)
        
        self.seqs = seqs
        
                
        counts = asarray(seqs.tally())
        #print repr(seqs.alphabet)
        #print seqs
        #print seqs.tally()
        #sys.exit()
        
        entropy_interval = None
        
        if prior is None :
            R = log(len(seqs.alphabet))
            ent = zeros( shape= seq_length, type=Float64)
            for i in range (0, seq_length) :
                C = sum(counts[i]) 
                #FIXME: fixup corebio.moremath.entropy()?
                if C == 0 :
                    ent[i] = 0.0
                else :
                    ent[i] = R - entropy(counts[i])
        else :
            ent = zeros( shape= seq_length, type=Float64)
            entropy_interval = zeros( shape=(seq_length,2) , type=Float64)
        
            R = log(len(seqs.alphabet))
            
            
            for i in range (0, seq_length) :
                alpha = array(counts[i] , type=Float64)
                alpha += prior
                
                posterior = Dirichlet(alpha)
                ent[i] = posterior.mean_relative_entropy(prior/sum(prior)) 
                entropy_interval[i][0], entropy_interval[i][1] = \
                    posterior.interval_relative_entropy(prior/sum(prior), 0.95) 
 
        
        self.seqs = seqs
        self.counts = counts

        self.entropy = ent
        self.entropy_interval = entropy_interval



# ====================== Main: Parse Command line =============================
def main(): 
    """WebLogo command line interface """
    
    # ------ Parse Command line ------
    parser = _build_option_parser()                
    (opts, args) = parser.parse_args(sys.argv[1:])
    if args : parser.error("Unparsable arguments: %s " % args)
    
    if opts.serve:
        httpd_serve_forever(opts.port) # Never returns?
        sys.exit(0) 
  
            
    # ------ Create Logo ------
    try:
        data = _build_logodata(opts)
        format = _build_logoformat(data, opts)
        
        formatter = opts.formatter
        formatter(data, format, opts.fout)

    except ValueError, err :
        print >>sys.stderr, 'Error:', err
        sys.exit(2)
    except KeyboardInterrupt, err:
        sys.exit(0)
# End main()            
        

def httpd_serve_forever(port=8080) :
    """ Start a webserver on a local port."""
    import BaseHTTPServer
    import CGIHTTPServer 
    
    class __HTTPRequestHandler(CGIHTTPServer.CGIHTTPRequestHandler):
        def is_cgi(self) :
            if self.path == "/create.cgi": 
                self.cgi_info = '', 'create.cgi'
                return True
            return False
    
    # Add current directory to PYTHONPATH. This is
    # so that we can run the standalone server
    # without having to run the install script.      
    pythonpath = os.getenv("PYTHONPATH", '')
    pythonpath += ":" + os.path.abspath(sys.path[0])
    os.environ["PYTHONPATH"] = pythonpath
    

    # TODO: Add command line option to override default webpage directory
    htdocs = os.path.join(os.path.dirname(__file__), "weblogo_htdocs") 
    os.chdir(htdocs) 

    
    HandlerClass = __HTTPRequestHandler
    ServerClass = BaseHTTPServer.HTTPServer
    httpd = ServerClass(('', port), HandlerClass)
    print "Serving HTTP on localhost:%d ..." % port
    
    try :
        httpd.serve_forever()
    except KeyboardInterrupt:
        sys.exit(0)
# end httpd_serve_forever()  
    
    

def _build_logodata(options) :
    seqs = read_seq_data(options.fin, 
        options.input_parser.read,
        alphabet=options.alphabet)      

    # FIXME: Error handling
    prior = parse_prior( options.composition,seqs.alphabet, options.weight)
    data = LogoData(seqs, prior)

    return data
     
             
def _build_logoformat( logodata, opts) :
    """ Extract and process relevant option values and return a 
    LogoFormat object.""" 

    args = {}  
    direct_from_opts = [
        "stacks_per_line", 
        "logo_title",
        "yaxis_label", 
        "show_xaxis",
        "show_yaxis",        
        "xaxis_label", 
        "show_ends",
        "fineprint",  
        "show_errorbars", 
        "show_boxes",  
        "yaxis_tic_interval", 
        "resolution",      
        "alphabet",
        "debug",
        "show_ends",
        "default_color",
        #"show_color_key",
        "color_scheme",
        "unit_name",
        "logo_label",
        "yaxis_scale",
        "first_index",
        "logo_start",
        "logo_end",
        "scale_width", 
        ]
  
    for k in direct_from_opts:
        args[k] = opts.__dict__[k]

    logo_size = copy.copy(opts.__dict__['logo_size'])
    size_from_opts = ["stack_width", "stack_height"]
    for k in size_from_opts :
        length = getattr(opts, k)
        if length : setattr( logo_size, k, length )
    args["size"] = logo_size    


    if opts.colors:
        color_scheme = ColorScheme()
        for color, symbols, desc in opts.colors:
            try :
                #c = Color.from_string(color)
                color_scheme.groups.append( ColorGroup(symbols, color, desc)  )
            except ValueError : 
                raise ValueError(
                     "error: option --color: invalid value: '%s'" % color )          
                 
        args["color_scheme"] = color_scheme
    
    logooptions = LogoOptions() 
    logooptions.__dict__.update(args) # FIXME: inelegant

    
    theformat =  LogoFormat(logodata, logooptions )
    return theformat

    
   
    
   
    
# ========================== OPTIONS ==========================
def _build_option_parser() :
    defaults = LogoOptions()
    parser = DeOptionParser(usage="%prog [options]  < sequence_data.fa > sequence_logo.eps",
        description = description,
        version     = __version__ ,
        add_verbose_options = False
        )    
        
    io_grp = OptionGroup(parser, "Input/Output Options",)
    data_grp = OptionGroup(parser, "Logo Data Options",)
    format_grp = OptionGroup(parser, "Logo Format Options", 
        "These options control the format and display of the logo.")
    color_grp = OptionGroup(parser, "Color Options", 
        "Colors can be specified using CSS2 syntax. e.g. 'red', '#FF0000', etc.")
    advanced_grp = OptionGroup(parser, "Advanced Format Options", 
        "These options provide fine control over the display of the logo. ")
    server_grp = OptionGroup(parser, "WebLogo Server",
        "Run a standalone webserver on a local port.")


    parser.add_option_group(io_grp)
    parser.add_option_group(data_grp)
    parser.add_option_group(format_grp)  
    parser.add_option_group(color_grp)
    parser.add_option_group(advanced_grp)
    parser.add_option_group(server_grp)
    
    # ========================== IO OPTIONS ==========================
    


    io_grp.add_option( "-f", "--fin",
        dest="fin",
        action="store",
        type="file_in",
        default=sys.stdin,
        help="Sequence input file (default: stdin)",
        metavar="FILENAME")

    io_grp.add_option("", "--fin-format", 
        dest="input_parser",
        action="store", type ="dict",
        default = seq_io,
        choices = seq_io.format_names(),
        help="Multiple sequence alignment format: (%s)" % 
           ', '.join([ f.names[0] for f in seq_io.formats]),
        metavar="FORMAT")

    io_grp.add_option("-o", "--fout", dest="fout",
        type="file_out",
        default=sys.stdout,
        help="Output file (default: stdout)",
        metavar="FILENAME")

    io_grp.add_option( "-F", "--format",
        dest="formatter",
        action="store",
        type="dict",
        choices = formatters,
        metavar= "FORMAT",
        help="Format of output: eps (default), png, pdf, gif, txt",
        default = default_formatter)
        #TODO Implement all these formatters

    # ========================== Data OPTIONS ==========================


 
    data_grp.add_option( "-A", "--sequence-type",
        dest="alphabet",
        action="store",
        type="dict",
        choices = std_alphabets,
        help="The type of sequence data: 'protein', 'rna' or 'dna'.",
        metavar="TYPE")
                       
    data_grp.add_option( "-a", "--alphabet",
        dest="alphabet",
        action="store",
        help="The set of symbols to count, e.g. 'AGTC'. "
             "All characters not in the alphabet are ignored. "
             "If neither the alphabet nor sequence-type are specified then weblogo will examine the input data and make an educated guess. "
             "See also --sequence-type, --cased" )                       
    
    # FIXME Add test?
    data_grp.add_option( "", "--ignore-lower-case",
        dest="cased",
        action="store",
        type = "boolean",
        default=False,
        metavar = "YES/NO",
        help="Disregard lower case letters and only count upper case letters in sequences?"
       )
       
    # FIXME: Check that these are all valid units
    data_grp.add_option( "-U", "--units",
        dest="unit_name",
        action="store",
        choices = std_units.keys(),    
        type="choice",
        default = defaults.unit_name,
        help="A unit of entropy ('bits' (default), 'nats'), or a unit of free energy ('kT', 'kJ/mol', 'kcal/mol'), or 'frequency' for probabilities",
        metavar = "NUMBER") 


    data_grp.add_option( "", "--composition",
        dest="composition",
        action="store",
        type="string",
        default = "auto",
        help="The expected composition of the sequences: 'auto' (default), 'equiprobable', 'none' (Do not perform any compositional adjustment), a CG percentage, a species name (e.g. 'E. coli', 'H. sapiens'), or an explicit distribution (e.g. {'A':10, 'C':40, 'G':40, 'T':10}). The automatic option uses a typical distribution for proteins and equiprobable distribution for everything else. ",
        metavar="COMP.")

    data_grp.add_option( "", "--weight",
        dest="weight",
        action="store",
        type="float",
        default = None,
        help="The weight of prior data.  Default: total pseudocounts equal to the number of monomer types.",
        metavar="NUMBER")


    # FIXME: Better name needed. Coordinate?
    data_grp.add_option( "-i", "--first-index",
        dest="first_index",
        action="store",
        type="int",
        default = 1,
        help="Index of first position in sequence data (default: 1)",
        metavar="INDEX")


    # FIXME: Upper, Lower should be format options?
    # Need refactoring anyways.
    data_grp.add_option( "-l", "--lower",
        dest="logo_start",
        action="store",
        type="int",
        help="Lower bound of sequence to display",
        metavar="INDEX")
    
    data_grp.add_option( "-u", "--upper",
        dest="logo_end",
        action="store",
        type="int",
        help="Upper bound of sequence to display",
        metavar="INDEX")

    # ========================== FORMAT OPTIONS ==========================

    #TODO: YAxis Interval. Default maximum.
    
    format_grp.add_option( "-s", "--size",
        dest="logo_size",
        action="store",
        type ="dict",
        choices = std_sizes,
        metavar = "LOGOSIZE",
        default = defaults.size,
        help="Specify a standard logo size (small, medium (default), large)" )
            
    format_grp.add_option( "-n", "--stacks-per-line",
        dest="stacks_per_line",
        action="store",
        type="int",
        help="Maximum number of logo stacks per logo line. (default: %default)",
        default = defaults.stacks_per_line,
        metavar="COUNT")

    format_grp.add_option( "-t", "--title",
        dest="logo_title",
        action="store",
        type="string",
        help="Logo title text.",
        default = defaults.logo_title,
        metavar="TEXT")

    format_grp.add_option( "", "--label",
        dest="logo_label",
        action="store",
        type="string",
        help="A figure label, e.g. '2a'",
        default = defaults.logo_label,
        metavar="TEXT")

    format_grp.add_option( "-X", "--show-xaxis",
        action="store",
        type = "boolean",
        default= defaults.show_xaxis,
        metavar = "YES/NO",
        help="Display sequence numbers along x-axis? (default: %default)")
                       
    format_grp.add_option( "-x", "--xlabel",
        dest="xaxis_label",
        action="store",
        type="string",
        default = defaults.xaxis_label,
        help="X-axis label",
        metavar="TEXT")

    format_grp.add_option( "-S", "--yaxis",
        dest="yaxis_scale",
        action="store",
        type="float",
        help="Height of yaxis in units. ",
        metavar = "UNIT") 

    format_grp.add_option( "-Y", "--show-yaxis",
        action="store",
        type = "boolean",
        dest = "show_yaxis",
        default= defaults.show_yaxis,
        metavar = "YES/NO",
        help="Display entropy scale along y-axis? (default: %default)")

    format_grp.add_option( "-y", "--ylabel",
        dest="yaxis_label",
        action="store",
        type="string",
        help="Y-axis label  (default depends on plot type and units)",
        metavar="TEXT")

    format_grp.add_option( "-E", "--show-ends",
        action="store",
        type = "boolean",
        default= defaults.show_ends,
        metavar = "YES/NO",
        help="Label the ends of the sequence? (default: %default)")
        
    format_grp.add_option( "-P", "--fineprint",
        dest="fineprint",
        action="store",
        type="string",
        default= defaults.fineprint,
        help="The fine print (default: weblogo version)",
        metavar="TEXT")

    format_grp.add_option( "", "--ticmarks",
        dest="yaxis_tic_interval",
        action="store",
        type="float",
        default= defaults.yaxis_tic_interval,
        help="Distance between ticmarks (default: %default)",
        metavar = "NUMBER")

        
    format_grp.add_option( "", "--errorbars",
        dest = "show_errorbars",
        action="store",
        type = "boolean",
        default= defaults.show_errorbars,
        metavar = "YES/NO",
        help="Display error bars? (default: %default)")
     
       
        
    # ========================== Color OPTIONS ==========================
    # TODO: Future Feature
    # color_grp.add_option( "-K", "--color-key",
    #    dest= "show_color_key",
    #    action="store",
    #    type = "boolean",
    #    default= defaults.show_color_key,
    #    metavar = "YES/NO",
    #    help="Display a color key (default: %default)")
        
        
    color_scheme_choices = std_color_schemes.keys()    
    color_scheme_choices.sort()    
    color_grp.add_option( "-c", "--color-scheme",
        dest="color_scheme",
        action="store",
        type ="dict",
        choices = std_color_schemes,
        metavar = "SCHEME",
        default = None, # Auto
        help="Specify a standard color scheme (%s)" % \
            ", ".join(color_scheme_choices) )
            
    color_grp.add_option( "-C", "--color",
        dest="colors",
        action="append",
        metavar="COLOR SYMBOLS DESCRIPTION ",
        nargs = 3,
        default=[],
        help="Specify symbol colors, e.g. --color black AG 'Purine' --color red TC 'Pyrimidine' ")    

    color_grp.add_option( "", "--default-color",
        dest="default_color",
        action="store",
        metavar="COLOR",
        default= defaults.default_color,
        help="Symbol color if not otherwise specified.")


    # ========================== Advanced options =========================   
                
    advanced_grp.add_option( "-W", "--stack-width",
        dest="stack_width",
        action="store",
        type="float",
        default= None,
        help="Width of a logo stack (default: %s)"% defaults.size.stack_width,
        metavar="POINTS" )

    advanced_grp.add_option( "-H", "--stack-height",
        dest="stack_height",
        action="store",
        type="float",
        default= None,
        help="Height of a logo stack (default: %s)"%defaults.size.stack_height,
        metavar="POINTS" )    

    advanced_grp.add_option( "", "--box",
        dest="show_boxes",
        action="store",
        type = "boolean",
        default=False, #FIXME
        metavar = "YES/NO",
        help="Draw boxes around symbols? (default: no)")

    # TODO: Add this information to manual.
    # Low resolution bitmaps (DPI<300) are antialiased.
    advanced_grp.add_option( "", "--resolution",
        dest="resolution",
        action="store",
        type="float",
        default=96,
        help="Bitmap resolution in dots per inch (DPI).  (default: 96 DPI, except png_print, 600 DPI)",
        metavar="DPI")  

    advanced_grp.add_option( "", "--scale-width",
        dest="scale_width",
        action="store",
        type = "boolean",
        default= True, 
        metavar = "YES/NO",
        help="Scale the visible stack width by the fraction of symbols in the column?  (i.e. columns with many gaps of unknowns are narrow.)  (default: yes)")
   
    advanced_grp.add_option( "", "--debug",
        action="store",
        type = "boolean",
        default= defaults.debug,
        metavar = "YES/NO",
        help="Output additional diagnostic information. (default: %default)")


    # ========================== Server options =========================   
    server_grp.add_option( "", "--serve",
        dest="serve",
        action="store_true",
        default= False,
        help="Start a standalone WebLogo server for creating sequence logos.")    
    
    server_grp.add_option( "", "--port",
        dest="port",
        action="store",
        type="int",
        default= 8080,
        help="Listen to this local port. (Default: %default)",
        metavar="PORT")

    # FIXME: Implement
   # parser.add_option( "", "--max",
   #     dest="max_file_size",
   #     action="store",
   #     type="int",
   #     help="Maximum size of input sequence file. (default: no limit)",
   #     metavar="BYTES")

    # FIXME: Implement
    #server_grp.add_option( "", "--htdocs",
    #    dest="htdocs",
    #    action="store",
    #    type="string",
    #    default= "weblogo_htdocs",
    #    help="Directory containing weblogo's web-pages. Default is weblogo_htdocs inside current directory."
   #     )


    return parser
    
    # END _build_option_parser




# The postscript template used by eps_formattter 
_template = r"""%!PS-Adobe-3.0 EPSF-3.0
%%Title:        Sequence Logo: ${logo_title}
%%Creator:      ${creator_text}
%%CreationDate: ${creation_date}
%%BoundingBox:  0  0  ${logo_width}  ${logo_height} 
%%Pages: 0
%%DocumentFonts: 
%%EndComments


% ---- VARIABLES ----

/True   true def
/False  false def

/debug              ${debug}        def

/logo_height        ${logo_height}  def
/logo_width         ${logo_width}  def
/logo_title         (${logo_title}) def
/show_title         ${show_title} def

/logo_margin        ${logo_margin} def
/xaxis_label_height ${xaxis_label_height} def
/title_height       ${title_height} def
/stroke_width       ${stroke_width} def
/tic_length         ${tic_length} def

/lines_per_logo     ${lines_per_logo} def
/line_width         ${line_width} def
/line_height        ${line_height} def
/line_margin_left   ${line_margin_left} def
/line_margin_right  ${line_margin_right} def
/line_margin_bottom ${line_margin_bottom} def
/line_margin_top    ${line_margin_top} def

/stack_width         ${stack_width} def
/stack_height        ${stack_height} def
/stacks_per_line     ${stacks_per_line} def
/stack_margin        ${stack_margin} def

/show_yaxis             ${show_yaxis} def      
/show_yaxis_label       ${show_yaxis_label} def
/yaxis_label            (${yaxis_label}) def
/yaxis_scale          ${yaxis_scale} def              % height in units 
/yaxis_tic_interval     ${yaxis_tic_interval} def           % in units
/yaxis_minor_tic_interval ${yaxis_minor_tic_interval} def   % in units

/show_xaxis_label       ${show_xaxis_label} def             % True or False
/show_xaxis             ${show_xaxis} def                   % True or False
/xaxis_label            (${xaxis_label}) def
/xaxis_tic_interval     ${xaxis_tic_interval} def
/rotate_numbers         ${rotate_numbers} def               % True or False
/number_interval        ${number_interval} def
/show_ends              ${show_ends} def          
/end_type               (${end_type}) def          % d: DNA, p: PROTEIN, -: none

/show_fineprint         ${show_fineprint} def
/fineprint              (${fineprint}) def
/logo_label             (${logo_label}) def

/show_boxes             ${show_boxes} def    % True or False
/shrink                 ${shrink} def    % True or False
/shrink_fraction        ${shrink_fraction} def             
/outline_linewidth      ${outline_linewidth} def    

/show_errorbars         ${show_errorbars} def      % True or False
/errorbar_fraction      ${errorbar_fraction} def
/errorbar_width_fraction  ${errorbar_width_fraction} def
/errorbar_gray          ${errorbar_gray} def

/fontsize               ${fontsize} def
/small_fontsize         ${small_fontsize} def
/title_fontsize         ${title_fontsize} def
/number_fontsize        ${number_fontsize} def



/default_color ${default_color} def 
/color_dict << 
${color_dict}
>> def



% ---- DERIVED PARAMETERS ----

/char_width stack_width 2 stack_margin mul sub def
/char_width2 char_width 2 div def
/char_width4 char_width 4 div def

% movements to place 5'/N and 3'/C symbols
/leftEndDeltaX  fontsize neg         def
/leftEndDeltaY  fontsize 1.25 mul neg def
/rightEndDeltaX fontsize 0.25 mul     def
/rightEndDeltaY leftEndDeltaY        def


% ---- PROCEDURES ----


/SetTitleFont {/${title_font} findfont title_fontsize scalefont setfont} bind def
/SetLogoFont  {/${logo_font} findfont char_width  scalefont setfont} bind def
/SetStringFont{/${text_font} findfont fontsize scalefont setfont} bind def
/SetPrimeFont {/Symbol findfont fontsize scalefont setfont} bind def
/SetSmallFont {/${text_font} findfont small_fontsize scalefont setfont} bind def
/SetNumberFont {/${text_font} findfont number_fontsize scalefont setfont} bind def

/DrawBox { % width height 
    /hh exch def
    /ww exch def
    gsave
        0.2 setlinewidth
        %0.5 setgray
        
        %0 0 moveto 
        hh 0 rlineto
        0 ww rlineto
        hh neg 0 rlineto
        0 ww neg rlineto
        stroke
    grestore
} bind def


/StartLogo { 
  %save 
  gsave 

  
  debug { 
    logo_margin logo_margin moveto
    logo_height logo_margin 2 mul sub
    logo_width logo_margin 2 mul sub
    DrawBox } if
    
  show_title { DrawTitle } if
  show_xaxis_label { DrawXaxisLable } if
  show_fineprint { DrawFineprint } if
  DrawLogoLabel
  
  
  MoveToFirstLine
} bind def


/DrawLogoLabel {
  gsave 
    SetTitleFont
    
    logo_margin
    logo_height title_fontsize sub logo_margin sub
    moveto
    
    debug { title_fontsize logo_label stringwidth pop DrawBox } if
    0 title_fontsize 4 div rmoveto % Move up to baseline (approximatly)
    logo_label show
    
  grestore
} bind def

/DrawTitle {
  gsave 
    SetTitleFont
    
    logo_width 2 div logo_title stringwidth pop 2 div sub 
    logo_height title_fontsize sub logo_margin sub
    moveto
    
    debug { title_fontsize logo_title stringwidth pop DrawBox } if
    
    0 title_fontsize 4 div rmoveto % Move up to baseline (approximatly)
    logo_title show
    
  grestore
} bind def

/DrawXaxisLable {
  % Print X-axis label, bottom center
  gsave
    SetStringFont

    logo_width 2 div xaxis_label stringwidth pop 2 div sub
    xaxis_label_height logo_margin add fontsize sub
    moveto
    %fontsize 3 div
    
    debug { fontsize xaxis_label stringwidth pop DrawBox } if

    xaxis_label show

  grestore
} bind def


/DrawFineprint {
    gsave
    
    SetSmallFont

    logo_width fineprint stringwidth pop sub 
        logo_margin sub line_margin_right sub
    logo_margin 
    moveto
    
    debug { small_fontsize fineprint stringwidth pop DrawBox } if
        
    fineprint show
    grestore
} bind def

/MoveToFirstLine {
    logo_margin 
    logo_height logo_margin sub title_height sub line_height sub
    moveto
} bind def

/EndLogo { 
  grestore 
  %showpage 
  %restore 
} bind def


/StartLine{ 
    gsave 
  
    % Draw outer box
    debug { line_height line_width DrawBox } if  

    % Move to lower left corner of content area
    line_margin_left line_margin_bottom rmoveto

    % Draw inner content box
    debug { 
        line_height line_margin_bottom sub line_margin_top sub
        line_width line_margin_left sub line_margin_right sub
        DrawBox
    } if  
  
    show_yaxis { DrawYaxis } if
    show_xaxis { DrawLeftEnd } if
    
} bind def

/EndLine{ 
    show_xaxis { DrawRightEnd } if
    grestore 
    0 line_height neg rmoveto 
} bind def


/DrawYaxis {
  gsave    
    stack_margin neg 0 translate
    DrawYaxisBar
    DrawYaxisLabel
  grestore
} bind def


/DrawYaxisBar { 
  gsave  
    stack_margin neg 0 rmoveto
    
    SetNumberFont
    stroke_width setlinewidth
    
    /str 10 string def % string to hold number  
    /smallgap stack_margin  def

    % Draw first tic and bar
    gsave    
      tic_length neg 0 rmoveto 
      tic_length 0 rlineto 
      0 stack_height rlineto
      stroke
    grestore

    % Draw the tics
    % initial increment limit proc for
    0 yaxis_tic_interval yaxis_scale abs 
    {/loopnumber exch def

      % convert the number coming from the loop to a string
      % and find its width
      loopnumber 10 str cvrs
      /stringnumber exch def % string representing the number

      stringnumber stringwidth pop
      /numberwidth exch def % width of number to show

      /halfnumberheight
         stringnumber CharBoxHeight 2 div
      def

      gsave
        numberwidth % move back width of number
        neg loopnumber stack_height yaxis_scale div mul % shift on y axis
        halfnumberheight sub % down half the digit
        rmoveto % move back the width of the string

        tic_length neg smallgap sub % Move back a bit more  
        0 rmoveto % move back the width of the tic  

        stringnumber show
        smallgap 0 rmoveto % Make a small gap  

        % now show the tic mark
        0 halfnumberheight rmoveto % shift up again
        tic_length 0 rlineto
        stroke
     grestore
    } for
    
    % Draw the minor tics
    % initial increment limit proc for
    0 yaxis_minor_tic_interval yaxis_scale abs 
    {/loopnumber2 exch def
      gsave
        0
        loopnumber2 stack_height yaxis_scale div mul 
        rmoveto 

        tic_length 2 div neg 0 rlineto
        stroke
     grestore
    } for    
    
  grestore
} bind def

/DrawYaxisLabel {
  gsave
    SetStringFont

    % How far we move left depends on the size of
    % the tic labels.
    /str 10 string def % string to hold number  
    yaxis_scale yaxis_tic_interval div cvi yaxis_tic_interval mul 
    str cvs stringwidth pop
    tic_length 1.25 mul  add neg  

    stack_height
    yaxis_label stringwidth pop
    sub 2 div

    rmoveto
    90 rotate
    
    yaxis_label show
  grestore
} bind def


%Take a single character and return the bounding box
/CharBox { % <char> CharBox <lx> <ly> <ux> <uy>
  gsave
    newpath
    0 0 moveto
    % take the character off the stack and use it here:
    true charpath 
    flattenpath 
    pathbbox % compute bounding box of 1 pt. char => lx ly ux uy
    % the path is here, but toss it away ...
  grestore
} bind def


% The height of a characters bounding box
/CharBoxHeight { % <char> CharBoxHeight <num>
  CharBox
  exch pop sub neg exch pop
} bind def


% The width of a characters bounding box
/CharBoxWidth { % <char> CharBoxHeight <num>
  CharBox
  pop exch pop sub neg 
} bind def


/DrawLeftEnd {
  gsave
    SetStringFont
    leftEndDeltaX leftEndDeltaY rmoveto
    
    show_ends {
        debug { leftEndDeltaY neg leftEndDeltaX neg DrawBox } if
        end_type (d) eq {(5) show DrawPrime} if
        end_type (p) eq {(N) show} if
    } if
  grestore
} bind def

/DrawRightEnd { 
  gsave
    SetStringFont
    rightEndDeltaX rightEndDeltaY rmoveto
    
    show_ends {
        debug { rightEndDeltaY neg leftEndDeltaX neg  DrawBox } if
        end_type (d) eq {(3) show DrawPrime} if
        end_type (p) eq {(C) show} if
    } if
  grestore
} bind def

/DrawPrime {
  gsave
    SetPrimeFont
    (\242) show 
  grestore
} bind def


/StartStack {  % <stackNumber> startstack
  show_xaxis {DrawNumber}{pop} ifelse
  gsave
  debug { stack_height stack_width DrawBox } if
    
} bind def

/EndStack {
  grestore
  stack_width 0 rmoveto
} bind def


/DrawNumber { % number MakeNumber
    /n exch def
    
    
  gsave
    %0 stack_margin neg rmoveto
    stroke_width setlinewidth
    stack_width  0 rlineto
    stack_width 2 div neg 0 rmoveto
    
    n () eq 
    {  0 tic_length 4 div neg rlineto  } 
    { 0 tic_length 2 div neg rlineto } 
    ifelse

    stroke
  grestore   



  gsave
    n
    SetNumberFont
    stack_width 2 div tic_length 2 div neg rmoveto 

    rotate_numbers {
        90 rotate     
        dup stringwidth pop neg % find the length of the number
        stack_margin sub        % Move down a bit
        (0) CharBoxHeight 2 div neg % left half height of numbers
        rmoveto 
        show
    } {
        dup stringwidth pop neg 2 div number_fontsize neg  rmoveto
        show
    } ifelse
    
    

  grestore
} bind def



% Draw a character whose height is proportional to symbol bits
/ShowSymbol{ % interval character ShowSymbol
    /char exch def
    /interval exch def
    /fraction_width exch def
    
    /char_height 
       interval yaxis_scale div stack_height mul
       stack_margin sub
       dup 
       % if char_height is negative or very small replace with zero 
       % BUG FIX: This used to be '0.0 gt' but it seems that DrawHeight
       % has a finite, non-zero minimum, which results in a rangecheck error
       0.001 gt {}{pop 0.0} ifelse 
    def 

    char_height 0.0 gt {
        show_boxes {
            gsave
                /ww char_height stack_margin add def
                /hh stack_width def
                stroke_width setlinewidth
                hh 0 rlineto
                0 ww rlineto
                hh neg 0 rlineto
                0 ww neg rlineto
                stroke
            grestore     
        } if

        gsave
            stack_margin stack_margin rmoveto
            debug { char_height char_width DrawBox } if
            1 fraction_width sub char_width mul 2 div  0 rmoveto
            fraction_width char_width mul char_height char DrawChar
        grestore
        
    } if
    0 interval yaxis_scale div stack_height mul rmoveto
} bind def


/DrawChar { % <width> <height> <char> ShowChar
    /tc exch def    % The character
    /ysize exch def % the y size of the character
    /xsize exch def % the x size of the character
    /xmulfactor 1 def 
    /ymulfactor 1 def
    
    gsave
        SetLogoFont    
        tc SetColor

        % IReplacementHack
        % Deal with the lack of bars on the letter 'I' in Arial and Helvetica
        % by replacing with 'I' from Courier.
        tc (I) eq {
            /Courier findfont char_width  scalefont setfont    
        } if


        shrink {
            xsize 1 shrink_fraction sub 2 div mul
            ysize 1 shrink_fraction sub 2 div mul rmoveto 
            shrink_fraction shrink_fraction scale
        } if
        
        % Calculate the font scaling factors
        % Loop twice to catch small correction due to first scaling
        2 {
            gsave
            xmulfactor ymulfactor scale
      
            ysize % desired size of character in points
            tc CharBoxHeight 
            dup 0.0 ne {
                div % factor by which to scale up the character
                /ymulfactor exch def
            } {pop pop} ifelse

            xsize % desired size of character in points
            tc CharBoxWidth  
            dup 0.0 ne {
                div % factor by which to scale up the character
                /xmulfactor exch def
            } {pop pop} ifelse
            grestore
        } repeat


  
        % Draw the character
        xmulfactor ymulfactor scale
        % Move lower left corner of character to start point
        tc CharBox pop pop % llx lly : Lower left corner
        exch neg exch neg
        rmoveto
        
        tc show    

    grestore
} bind def

/SetColor{ % <char> SetColor
  dup color_dict exch known {
    color_dict exch get aload pop setrgbcolor
  } {
    pop
    default_color aload pop setrgbcolor
  } ifelse 
} bind def


/DrawErrorbar{ % interval_down interval_up DrawErrorbar
    
    gsave
    /points_per_unit stack_height yaxis_scale div def 
    /height_up   exch points_per_unit mul def 
    /height_down exch points_per_unit mul def
    
    show_errorbars {
    
    stroke_width setlinewidth
    errorbar_gray setgray     
    stack_width 2 div 0 rmoveto
    
    /errorbar_width char_width errorbar_width_fraction mul def
    /errorbar_width2 errorbar_width 2 div def
    
    gsave 
        0 height_down neg rmoveto
        errorbar_width2 neg 0 rlineto
        errorbar_width 0 rlineto
        errorbar_width2 neg 0 rlineto
        0 height_down errorbar_fraction mul rlineto 
        stroke
    grestore

    gsave 
        0 height_up  rmoveto
        errorbar_width2 neg 0 rlineto
        errorbar_width 0 rlineto
        errorbar_width2 neg 0 rlineto
        0 height_up neg errorbar_fraction mul rlineto 
        stroke
    grestore
    } if
           
    grestore
    
} bind def

/DrawErrorbarFirst{ % interval_down interval_up center DrawErrorbarFirst
    gsave
    /points_per_unit stack_height yaxis_scale div def 
    /center   exch points_per_unit mul def 

    0 center rmoveto
    DrawErrorbar
    grestore
} bind def

%%EndProlog

%%Page: 1 1

% Example Data
%StartLogo
%    StartLine
%        (1) StartStack
%            1.2 (C) ShowSymbol
%            2.2 (I) ShowSymbol
%            0.5 0.5 DrawErrorbar
%        EndStack
%        (2) StartStack
%            0.5 (I) ShowSymbol
%            0.9 (L) ShowSymbol
%            1.0 (G) ShowSymbol
%            
%            0.5 0.5 DrawErrorbar        
%        EndStack
%        (234) StartStack
%        EndStack
%        (235) StartStack
%        EndStack
%    EndLine
%EndLogo

StartLogo

${logo_data}

EndLogo


%%EOF

"""

##############################################################

class Dirichlet(object) :
    """The Dirichlet probability distribution. The Dirichlet is a continuous 
    multivariate probability distribution across non-negative unit length
    vectors. In other words, the Dirichlet is a probability distribution of 
    probability distributions. It is conjugate to the multinomial
    distribution and is widely used in Bayesian statistics.
    
    The Dirichlet probability distribution of order K-1 is 

     p(theta_1,...,theta_K) d theta_1 ... d theta_K = 
        (1/Z) prod_i=1,K theta_i^{alpha_i - 1} delta(1 -sum_i=1,K theta_i)

    The normalization factor Z can be expressed in terms of gamma functions:

      Z = {prod_i=1,K Gamma(alpha_i)} / {Gamma( sum_i=1,K alpha_i)}  

    The K constants, alpha_1,...,alpha_K, must be positive. The K parameters, 
    theta_1,...,theta_K are nonnegative and sum to 1.
    
    Status:
        Alpha
    """
    __slots__ = 'alpha', '_total', '_mean', 
    
    
    

    def __init__(self, alpha) :
        """
        Args:
            - alpha  -- The parameters of the Dirichlet prior distribution.
                        A vector of non-negative real numbers.  
        """
        # TODO: Check that alphas are positive
        #TODO : what if alpha's not one dimensional?
        #from numarray import asarray, Float64
        self.alpha = asarray(alpha, Float64)
        
        self._total = sum(alpha)
        self._mean = None
        
        
    def sample(self) :
        """Return a randomly generated probability vector.
        
        Random samples are generated by sampling K values from gamma
        distributions with parameters a=\alpha_i, b=1, and renormalizing. 
    
        Ref:
            A.M. Law, W.D. Kelton, Simulation Modeling and Analysis (1991).
        Authors:
            Gavin E. Crooks <gec@compbio.berkeley.edu> (2002)
        """
        #from numarray import asarray, Float64, zeros
        alpha = self.alpha
        K = len(alpha)
        theta = zeros( (K,), Float64)

        for k in range(K):
            theta[k] = random.gammavariate(alpha[k], 1.0) 
        theta /= sum(theta)

        return theta

    def mean(self) :
        if  self._mean ==None:
            self._mean = self.alpha / self._total
        return self._mean
    
    def covariance(self) :
        #from numarray import asarray, Float64, zeros
        alpha = self.alpha
        A = sum(alpha)
        #A2 = A * A
        K = len(alpha)
        cv = zeros( (K,K), Float64) 
        
        for i in range(K) :
            cv[i,i] = alpha[i] * (1. - alpha[i]/A) / (A * (A+1.) )
        
        for i in range(K) :
            for j in range(i+1,K) :
                v = - alpha[i] * alpha[j] / (A * A * (A+1.) )
                cv[i,j] = v
                cv[j,i] = v
        return cv
        
    def mean_x(self, x) :
        
        x = asarray(x, Float64)
        if shape(x) != shape(self.alpha) :
            raise ValueError("Argument must be same dimension as Dirichlet")
        # TODO: Check tha shape(x) == shape(alpha)
        return sum( x * self.mean()) 

    def variance_x(self, x) :
        #from numarray import asarray, Float64, zeros
        x = asarray(x, Float64)
        # x needs to be the same shape as alpha
        cv = self.covariance()
        
        var = na.matrixmultiply(na.matrixmultiply(na.transpose( x), cv), x)
        return var


    def mean_entropy(self) :
        """Calculate the average entropy of probabilities sampled
        from this Dirichlet distribution. 
        
        Returns:
            The average entropy.
            
        Ref:
            Wolpert & Wolf, PRE 53:6841-6854 (1996) Theorem 7
            (Warning: this paper contains typos.)
        Status:
            Alpha
        Authors:
            GEC 2005
    
        """
        # TODO: Optimize
        alpha = self.alpha
        A = float(sum(alpha))
        ent = 0.0
        for a in alpha:
            if a>0 : ent += - 1.0 * a * digamma( 1.0+a) # FIXME: Check
        ent /= A
        ent += digamma(A+1.0)
        return ent    



    def variance_entropy(self):
        """Calculate the variance of the Dirichlet entropy. 

        Ref:
            Wolpert & Wolf, PRE 53:6841-6854 (1996) Theorem 8
            (Warning: this paper contains typos.)
        """
        alpha = self.alpha
        A = float(sum(alpha))
        A2 = A * (A+1)
        L = len(alpha)
        
        dg1 = zeros( (L) , Float64)
        dg2 = zeros( (L) , Float64)
        tg2 = zeros( (L) , Float64)        
        
        for i in range(L) :
            dg1[i] = digamma(alpha[i] + 1.0)
            dg2[i] = digamma(alpha[i] + 2.0)
            tg2[i] = trigamma(alpha[i] + 2.0)

        dg_Ap2 = digamma( A+2. )
        tg_Ap2 = trigamma( A+2. )
        
        mean = self.mean_entropy()
        var = 0.0
    
        for i in range(L) :
            for j in range(L) :
                if i != j :
                    var += (
                        ( dg1[i] - dg_Ap2 ) * (dg1[j] - dg_Ap2 ) - tg_Ap2 
                        ) * (alpha[i] * alpha[j] ) / A2
                else : 
                    var += (
                        ( dg2[i] - dg_Ap2 ) **2 + ( tg2[i] - tg_Ap2 )
                        ) * ( alpha[i] * (alpha[i]+1.) ) / A2

        var -= mean**2
        return var
        
        
        
    def mean_relative_entropy(self, pvec) :
        ln_p = na.log(pvec)
        return - self.mean_x(ln_p) - self.mean_entropy() 
    
    
    def variance_relative_entropy(self, pvec) :
        ln_p = na.log(pvec)
        return self.variance_x(ln_p) + self.variance_entropy()
    
    
    def interval_relative_entropy(self, pvec, frac) :
        mean = self.mean_relative_entropy(pvec) 
        variance = self.variance_relative_entropy(pvec) 
        # If the variance is small, use the standard 95% 
        # confidence interval: mean +/- 1.96 * sd
        if variance< 0.1 :
            sd = sqrt(variance)
            return max(0.0, mean - sd*1.96), mean + sd*1.96
        
        g = Gamma.from_mean_variance(mean, variance)
        low_limit = g.inverse_cdf( (1.-frac)/2.)
        high_limit = g.inverse_cdf( 1. - (1.-frac)/2. )
        
        return low_limit, high_limit


# Standard python voodoo for CLI
if __name__ == "__main__":
    ## Code Profiling. Uncomment these lines
    #import hotshot, hotshot.stats
    #prof = hotshot.Profile("stones.prof")
    #prof.runcall(main)
    #prof.close()
    #stats = hotshot.stats.load("stones.prof")
    #stats.strip_dirs()
    #stats.sort_stats('cumulative', 'calls')
    #stats.print_stats(40)
    #sys.exit()

    main()
    





