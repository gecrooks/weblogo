# -------------------------------- WebLogo --------------------------------

#  Copyright (c) 2003-2004 The Regents of the University of California.
#  Copyright (c) 2005 Gavin E. Crooks
#  Copyright (c) 2006-2015, The Regents of the University of California, through
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

"""
"""

import os
import sys

from typing.io import TextIO
from datetime import datetime
from math import log, sqrt
import shutil
import tempfile
from io import StringIO
from urllib.request import urlopen, Request
from urllib.parse import urlparse, urlunparse

# Avoid 'from numpy import *' since numpy has lots of names defined
from numpy import array, asarray, float64, ones, zeros, any
import numpy as np

from .color import Color
from .colorscheme import (ColorScheme, SymbolColor, monochrome, hydrophobicity, base_pairing,
                          chemistry, charge)
from .logomath import Dirichlet


from . import seq_io
from .data import amino_acid_composition
from scipy.stats import entropy
from .seq import (Alphabet, unambiguous_dna_alphabet,
                  unambiguous_rna_alphabet, unambiguous_protein_alphabet, SeqList)
from .utils import (isfloat, ArgumentError, stdrepr)

from . import __version__

# Shorten development version string of the form weblogo-3.6.1.dev43+g64d9f12.d20190304
if __version__.find('+') != -1:
    __version__ = __version__[:__version__.find('+')]


# from .logo_formatter import (GhostscriptAPI, pdf_formatter, jpeg_formatter, png_formatter,
#                             png_print_formatter,
#                             txt_formatter, eps_formatter, formatters, default_formatter)
# ------ META DATA ------


# __all__ = ['LogoOptions',
#            'description',
#            '__version__',
#            'LogoFormat',
#            'LogoData',
#            'GhostscriptAPI',
#            'std_color_schemes',
#            'default_color_schemes',
#            'classic',
#            'std_units',
#            'std_sizes',
#            'std_alphabets',
#            'std_percentCG',
#            'pdf_formatter',
#            'jpeg_formatter',
#            'png_formatter',
#            'png_print_formatter',
#            'txt_formatter',
#            'eps_formatter',
#            'formatters',
#            'default_formatter',
#            'base_distribution',
#            'equiprobable_distribution',
#            'read_seq_data',
#            'Color',
#            'ColorScheme',
#            'parse_prior',
#            'release_description',
#            'description'
#            ]

description = "Create sequence logos from biological sequence alignments."

release_description = "WebLogo %s" % (__version__)


def cgi(htdocs_directory):  # pragma: no cover
    import weblogo._cgi
    weblogo._cgi.main(htdocs_directory)


aa_composition = [amino_acid_composition[_k] for _k in unambiguous_protein_alphabet]

# ------  DATA ------

classic = ColorScheme([
    SymbolColor("G", "orange"),
    SymbolColor("TU", "red"),
    SymbolColor("C", "blue"),
    SymbolColor("A", "green")
])

std_color_schemes = {"auto": None,  # Depends on sequence type
                     "monochrome": monochrome,
                     "base pairing": base_pairing,
                     "classic": classic,
                     "hydrophobicity": hydrophobicity,
                     "chemistry": chemistry,
                     "charge": charge,
                     }

default_color_schemes = {
    unambiguous_protein_alphabet: hydrophobicity,
    unambiguous_rna_alphabet: base_pairing,
    unambiguous_dna_alphabet: base_pairing
}

std_units = {
    "bits": 1. / log(2),
    "nats": 1.,
    "digits": 1. / log(10),
    "kT": 1.,
    "kJ/mol": 8.314472 * 298.15 / 1000.,
    "kcal/mol": 1.987 * 298.15 / 1000.,
    "probability": None,
}

# The base stack width is set equal to 9pt Courier.
# (Courier has a width equal to 3/5 of the point size.)
# Check that can get 80 characters in journal page @small
# 40 characters in a journal column
std_sizes = {
    "small": 5.4,
    "medium": 5.4 * 2,
    "large": 5.4 * 3
}

std_alphabets = {
    'protein': unambiguous_protein_alphabet,
    'rna': unambiguous_rna_alphabet,
    'dna': unambiguous_dna_alphabet}

std_percentCG = {
    'H. sapiens': 40.,
    'E. coli': 50.5,
    'S. cerevisiae': 38.,
    'C. elegans': 36.,
    'D. melanogaster': 43.,
    'M. musculus': 42.,
    'T. thermophilus': 69.4,
}


# Thermus thermophilus: Henne A, Bruggemann H, Raasch C, Wiezer A, Hartsch T,
# Liesegang H, Johann A, Lienard T, Gohl O, Martinez-Arias R, Jacobi C,
# Starkuviene V, Schlenczeck S, Dencker S, Huber R, Klenk HP, Kramer W,
# Merkl R, Gottschalk G, Fritz HJ: The genome sequence of the extreme
# thermophile Thermus thermophilus.
# Nat Biotechnol 2004, 22:547-53


class LogoOptions(object):
    """ A container for all logo formatting options. Not all of these
    are directly accessible through the CLI or web interfaces.

    To display LogoOption defaults::

    >>> from weblogo import *
    >>> LogoOptions()

    All physical lengths are measured in points. (72 points per inch, 28.3 points per cm)

    Args:
        creator_text:           Embedded as comment in figures.
        logo_title:             Creates title for the sequence logo
        logo_label:             An optional figure label, added to the top left (e.g. '(a)').
        unit_name:              See std_units for options. (Default 'bits')
        yaxis_label:            Defaults to unit_name
        xaxis_label:            Add a label to the x-axis, or hide x-axis altogether.
        fineprint:              Defaults to WebLogo name and version

        show_yaxis:              Display entropy scale along y-axis (default: True)
        show_xaxis:             Display sequence numbers along x-axis (default: True)
        show_ends:              Display label at the ends of the sequence (default: False)
        show_fineprint:          Toggle display of the WebLogo version information in the lower
                                      right corner. Optional, but we appreciate the acknowledgment.
        show_errorbars:          Draw errorbars (default: False)
        show_boxes:              Draw boxes around stack characters (default: True)
        debug:                   Draw extra graphics debugging information.
        rotate_numbers:          Draw xaxis numbers with vertical orientation?
        scale_width:             boolean, scale width of characters proportional to ungaps
        pad_right:               Make a single line logo the same width as multiline logos
                                      (default: False)

        stacks_per_line:          Maximum number of logo stacks per logo line. (Default: 40)
        yaxis_tic_interval:       Distance between ticmarks on y-axis(default: 1.0)
        yaxis_minor_tic_ratio:    Distance between minor tic ratio
        yaxis_scale:              Sets height of the y-axis in designated units
        xaxis_tic_interval:       Distance between ticmarks on x-axis(default: 1.0)
        number_interval:          Distance between ticmarks (default: 1.0)

        shrink_fraction:          Proportional shrinkage of characters if show_boxes is true.

        errorbar_fraction:        Sets error bars display proportion
        errorbar_width_fraction:  Sets error bars display
        errorbar_gray:            Sets error bars' gray scale percentage (default .75)

        resolution:               Dots per inch (default: 96). Used for bitmapped output
                                       formats

        default_color:            Symbol color if not otherwise specified
        color_scheme:             A custom color scheme can be specified using CSS2 (Cascading
                                       Style Sheet) syntax.
                                       E.g. 'red', '#F00', '#FF0000', 'rgb(255, 0, 0)',
                                       'rgb(100%, 0%, 0%)' or 'hsl(0, 100%, 50%)' for the color red.

        stack_width:              Scale the visible stack width by the fraction of symbols in
                                       the column (I.e. columns with many gaps of unknowns are
                                       narrow.)  (Default: yes)
        stack_aspect_ratio:       Ratio of stack height to width (default: 5)

        logo_margin:              Default: 2 pts
        stroke_width:             Default: 0.5 pts
        tic_length:               Default: 5 pts
        stack_margin:             Default: 0.5 pts

        small_fontsize:           Small text font size in points
        fontsize:                 Regular text font size in points
        title_fontsize:           Title text font size in points
        number_fontsize:          Font size for axis-numbers, in points.

        text_font:                Select font for labels
        logo_font:                Select font for Logo
        title_font:               Select font for Logo's title

        first_index:              Index of first position in sequence data
        logo_start:               Lower bound of sequence to display
        logo_end:                 Upper bound of sequence to display

    """

    def __init__(self, **kwargs):
        """ Create a new LogoOptions instance.

        >>> logooptions = LogoOptions(logo_title = "Some Title String")
        >>> logooptions.show_yaxis = False
        >>> repr(logooptions)

        """

        self.alphabet = None

        self.creator_text = release_description

        self.logo_title = ""
        self.logo_label = ""
        self.stacks_per_line = 40

        self.unit_name = "bits"

        self.show_yaxis = True
        # yaxis_lable default depends on other settings. See LogoFormat
        self.yaxis_label = None
        self.yaxis_tic_interval = 1.
        self.yaxis_minor_tic_ratio = 5
        self.yaxis_scale = None

        self.show_xaxis = True
        self.xaxis_label = ""
        self.xaxis_tic_interval = 1
        self.rotate_numbers = False
        self.number_interval = 5
        self.show_ends = False
        self.annotate = None

        self.show_fineprint = True
        self.fineprint = "WebLogo " + __version__

        self.show_boxes = False
        self.shrink_fraction = 0.5

        self.show_errorbars = True
        self.errorbar_fraction = 0.90
        self.errorbar_width_fraction = 0.25
        self.errorbar_gray = 0.75

        self.resolution = 96.  # Dots per inch

        self.default_color = Color.by_name("black")
        self.color_scheme = None
        # self.show_color_key = False # NOT yet implemented

        self.debug = False

        self.logo_margin = 2
        self.stroke_width = 0.5
        self.tic_length = 5

        self.stack_width = std_sizes["medium"]
        self.stack_aspect_ratio = 5

        self.stack_margin = 0.5
        self.pad_right = False

        self.small_fontsize = 6
        self.fontsize = 10
        self.title_fontsize = 12
        self.number_fontsize = 8

        self.text_font = "ArialMT"
        self.logo_font = "Arial-BoldMT"
        self.title_font = "ArialMT"

        self.first_index = 1
        self.logo_start = None
        self.logo_end = None
        self.scale_width = True

        self.reverse_stacks = True  # If true, draw stacks with largest letters on top.

        for k, v in kwargs.items():
            setattr(self, k, v)

    def __repr__(self):
        attributes = list(vars(self).keys())
        attributes.sort()
        return stdrepr(self, attributes)


# End class LogoOptions


class LogoFormat(LogoOptions):
    """ Specifies the format of the logo. Requires LogoData and LogoOptions
    objects.

    >>> logodata = LogoData.from_seqs(seqs)
    >>> logooptions = LogoOptions()
    >>> logooptions.title = "A Logo Title"
    >>> format = LogoFormat(logodata, logooptions)

    Raises:
        ArgumentError: if arguments are invalid.
    """

    def __init__(self, logodata, logooptions=None):
        """ Create a new LogoFormat instance.

        """

        if logooptions is not None:
            self.__dict__.update(logooptions.__dict__)

        self.alphabet = logodata.alphabet
        self.seqlen = logodata.length

        # Derived parameters.
        self.show_title = False
        self.show_xaxis_label = False
        self.yaxis_minor_tic_interval = None
        self.lines_per_logo = None
        self.char_width = None  # Maximum character width. Stack width minus margins.
        self.line_margin_left = None
        self.line_margin_right = None
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

        self.stack_height = self.stack_width * self.stack_aspect_ratio

        # Attribute to test, test, error message
        arg_conditions = (
            ("stacks_per_line", lambda x: x > 0, "Stacks per line must be positive."),
            ("stack_width", lambda x: x > 0.0, "Stack width must be greater than zero."),
            ("stack_aspect_ratio", lambda x: x > 0,
                "Stack aspect ratio must be greater than zero."),
            ("fontsize", lambda x: x > 0, "Font sizes must be positive."),
            ("small_fontsize", lambda x: x > 0, "Font sizes must be positive."),
            ("title_fontsize", lambda x: x > 0, "Font sizes must be positive."),
            ("errorbar_fraction", lambda x: x >= 0.0 and x <= 1.0,
             "The visible fraction of the error bar must be between zero and one."),
            ("yaxis_tic_interval", lambda x: x >= 0.0,
                "The yaxis tic interval cannot be negative."),
            ("yaxis_minor_tic_interval", lambda x: not (x and x < 0.0),
                "Distances cannot be negative."),
            ("xaxis_tic_interval", lambda x: x > 0.0, "Tic interval must be greater than zero."),
            ("number_interval", lambda x: x > 0.0, "Invalid interval between numbers."),
            ("shrink_fraction", lambda x: x >= 0.0 and x <= 1.0, "Invalid shrink fraction."),
            ("stack_margin", lambda x: x > 0.0, "Invalid stack margin."),
            ("logo_margin", lambda x: x > 0.0, "Invalid logo margin."),
            ("stroke_width", lambda x: x > 0.0, "Invalid stroke width."),
            ("tic_length", lambda x: x > 0.0, "Invalid tic length."),
        )

        # Run arguments tests. The second, attribute argument to the ArgumentError is
        # used by the UI to provide user feedback.
        # FIXME: More validation
        for test in arg_conditions:
            if not test[1](getattr(self, test[0])):
                raise ArgumentError(test[2], test[0])

        # Inclusive upper and lower bounds
        # FIXME: Validate here. Move from eps_formatter
        if self.logo_start is None:
            self.logo_start = self.first_index

        if self.logo_end is None:
            self.logo_end = self.seqlen + self.first_index - 1

        self.total_stacks = self.logo_end - self.logo_start + 1
        if self.total_stacks <= 0:
            raise ArgumentError("Logo must contain at least one stack", 'logo_end')

        if self.logo_start - self.first_index < 0:
            raise ArgumentError(
                    "Logo range extends before start of available sequence.",
                    'logo_range')

        if self.logo_end - self.first_index >= self.seqlen:
            raise ArgumentError(
                    "Logo range extends beyond end of available sequence.",
                    'logo_range')

        if self.logo_title:
            self.show_title = True
        if not self.fineprint:
            self.show_fineprint = False
        if self.xaxis_label:
            self.show_xaxis_label = True

        if self.yaxis_label is None:
            self.yaxis_label = self.unit_name

        if self.yaxis_label:
            self.show_yaxis_label = True
        else:
            self.show_yaxis_label = False
            self.show_ends = False

        if not self.yaxis_scale:
            conversion_factor = std_units[self.unit_name]
            if conversion_factor:
                if self.alphabet is None:
                    raise ArgumentError("Need an alphabet", 'alphabet')

                self.yaxis_scale = log(len(self.alphabet)) * conversion_factor
            else:
                self.yaxis_scale = 1.0  # probability units

        if self.yaxis_scale <= 0.0:
            raise ArgumentError("Invalid yaxis scale", 'yaxis_scale', )

        if self.yaxis_tic_interval >= self.yaxis_scale:
            self.yaxis_tic_interval /= 2.

        self.yaxis_minor_tic_interval \
            = float(self.yaxis_tic_interval) / self.yaxis_minor_tic_ratio

        if self.color_scheme is None:
            if self.alphabet in default_color_schemes:
                self.color_scheme = default_color_schemes[self.alphabet]
            else:
                self.color_scheme = monochrome

        self.lines_per_logo = 1 + ((self.total_stacks - 1) // self.stacks_per_line)

        if self.lines_per_logo == 1 and not self.pad_right:
            self.stacks_per_line = min(self.stacks_per_line, self.total_stacks)

        self.char_width = self.stack_width - 2 * self.stack_margin

        if self.show_yaxis:
            self.line_margin_left = self.fontsize * 3.0
        else:
            if self.show_ends and self.show_xaxis:
                self.line_margin_left = self.fontsize * 1.5
            else:
                self.line_margin_left = 4

        if self.show_ends and self.show_xaxis:
            self.line_margin_right = self.fontsize * 1.5
        else:
            self.line_margin_right = 4

        if self.show_xaxis:
            if self.rotate_numbers:
                self.line_margin_bottom = self.number_fontsize * 2.5
            else:
                self.line_margin_bottom = self.number_fontsize * 1.5
        else:
            self.line_margin_bottom = 4

        self.line_margin_top = 4

        if self.show_title:
            self.title_height = self.title_fontsize
        else:
            self.title_height = 0

        self.xaxis_label_height = 0.
        if self.show_xaxis_label:
            self.xaxis_label_height += self.fontsize
        if self.show_fineprint:
            if len(self.fineprint) != 0:
                self.xaxis_label_height += self.small_fontsize

        self.line_height = (self.stack_height + self.line_margin_top +
                            self.line_margin_bottom)
        self.line_width = (self.stack_width * self.stacks_per_line +
                           self.line_margin_left + self.line_margin_right)

        self.logo_height = int(2 * self.logo_margin + self.title_height
                               + self.xaxis_label_height + self.line_height * self.lines_per_logo)
        self.logo_width = int(2 * self.logo_margin + self.line_width)

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

        if self.annotate is None:
            self.annotate = []
            for i in range(self.seqlen):
                index = i + self.first_index
                if index % self.number_interval == 0:
                    self.annotate.append("%d" % index)
                else:
                    self.annotate.append("")

        if len(self.annotate) != self.seqlen:
            raise ArgumentError(
                    "Annotations must be same length as sequences.",
                    'annotate')

    # End __init__


# End class LogoFormat


def parse_prior(composition, alphabet, weight=None):
    """ Parse a description of the expected monomer distribution of a sequence.

    Valid compositions:

    * None or 'none'
        No composition sepecified
    * 'auto' or 'automatic'
        Use the typical average distribution
        for proteins and an equiprobable distribution for
        everything else.
    * 'equiprobable'
        All monomers have the same probability.
    * a percentage, e.g. '45%' or a fraction '0.45'
        The fraction of CG bases for nucleotide alphabets
    * a species name, e.g. 'E. coli', 'H. sapiens',
        Use the average CG percentage for the species's genome.
    * An explicit distribution
        e.g. {'A':10, 'C':40, 'G':40, 'T':10}
    """
    if composition is None:
        return None
    comp = composition.strip()

    if comp.lower() == 'none':
        return None

    if weight is None and alphabet is not None:
        weight = sqrt(float(len(alphabet)))

    if weight < 0:
        raise ValueError("Weight cannot be negative.")

    if comp.lower() == 'equiprobable':
        prior = weight * equiprobable_distribution(len(alphabet))
    elif comp.lower() == 'auto' or comp.lower() == 'automatic':
        if alphabet == unambiguous_protein_alphabet:
            prior = weight * asarray(aa_composition, float64)
        else:
            prior = weight * equiprobable_distribution(len(alphabet))

    elif comp in std_percentCG:
        prior = weight * base_distribution(std_percentCG[comp])

    elif comp[-1] == '%':
        prior = weight * base_distribution(float(comp[:-1]))

    elif isfloat(comp):
        prior = weight * base_distribution(float(comp) * 100.)

    elif composition[0] == '{' and composition[-1] == '}':
        explicit = composition[1: -1]
        explicit = explicit.replace(',', ' ').replace("'", ' ').replace('"', ' ') \
            .replace(':', ' ').split()

        if len(explicit) != len(alphabet) * 2:
            raise ValueError("Explicit prior does not match length of alphabet")
        prior = - ones(len(alphabet), float64)
        try:
            for r in range(len(explicit) // 2):
                letter = explicit[r * 2]
                index = alphabet.ord(letter)
                value = float(explicit[r * 2 + 1])
                prior[index] = value
        except ValueError:
            raise ValueError("Cannot parse explicit composition")

        if any(prior == -1.):
            raise ValueError("Explicit prior does not match alphabet")  # pragma: no cover
        prior /= sum(prior)
        prior *= weight

    else:
        raise ValueError("Unknown or malformed composition: %s" % composition)

    if len(prior) != len(alphabet):
        raise ValueError(
                "The sequence alphabet and composition are incompatible.")  # pragma: no cover
    return prior


def base_distribution(percentCG):
    A = (1. - (percentCG / 100.)) / 2.
    C = (percentCG / 100.) / 2.
    G = (percentCG / 100.) / 2.
    T = (1. - (percentCG / 100)) / 2.
    return asarray((A, C, G, T), float64)


def equiprobable_distribution(length: int) -> np.ndarray:
    return ones((length), float64) / length


def _seq_formats():
    """ Return a dictionary mapping between the names of formats for the sequence data
    and the corresponing parsers.
    """
    # Add position weight matrix formats to input parsers by hand
    fin_choices = dict(seq_io.format_names())
    fin_choices['transfac'] = 'transfac'
    del fin_choices['plain']
    return fin_choices


def _seq_names():
    """ Returns a list of the names of accepted sequence data formats."""
    fin_names = [f.names[0] for f in seq_io.formats]
    fin_names.remove('plain')
    fin_names.append('transfac')
    return fin_names


def read_seq_data(fin,
                  input_parser=seq_io.read,
                  alphabet=None,
                  ignore_lower_case=False,
                  max_file_size=0):
    """Read sequence data from the input stream and return a seqs object.

    The environment variable WEBLOGO_MAX_FILE_SIZE overides the max_file_size argument.
    Used to limit the load on the WebLogo webserver.
    """

    max_file_size = int(os.environ.get("WEBLOGO_MAX_FILE_SIZE", max_file_size))

    # If max_file_size is set, or if fin==stdin (which is non-seekable), we
    # read the data and replace fin with a StringIO object.
    if (max_file_size > 0):
        data = fin.read(max_file_size)
        more_data = fin.read(2)
        if more_data != "":
            raise IOError("File exceeds maximum allowed size: %d bytes" % max_file_size)
        fin = StringIO(data)
    elif fin == sys.stdin:
        fin = StringIO(fin.read())

    fin.seek(0)
    seqs = input_parser(fin)

    if seqs is None or len(seqs) == 0:
        raise ValueError("Please provide a multiple sequence alignment")

    if ignore_lower_case:
        # Case is significant. Do not count lower case letters.
        for i, s in enumerate(seqs):
            seqs[i] = s.mask()

    # Add alphabet to seqs.
    if alphabet:
        seqs.alphabet = Alphabet(alphabet)
    else:
        seqs.alphabet = Alphabet.which(seqs)
    return seqs


class LogoData(object):
    """The data needed to generate a sequence logo.

    Args:
        alphabet: The set of symbols to count.
                   See also --sequence-type, --ignore-lower-case
        length:   All sequences must be the same length, else WebLogo will return an error
        counts:   An array of character counts
        entropy:  The relative entropy of each column
        entropy_interval: entropy confidence interval
     """

    def __init__(self, length=None, alphabet=None, counts=None,
                 entropy=None, entropy_interval=None, weight=None):
        """Creates a new LogoData object"""
        self.length = length
        self.alphabet = alphabet
        self.counts = counts
        self.entropy = entropy
        self.entropy_interval = entropy_interval
        self.weight = weight

    @classmethod
    def from_counts(cls,
                    alphabet: Alphabet,
                    counts: np.ndarray,
                    prior: np.ndarray = None) -> 'LogoData':
        """Build a LogoData object from counts."""
        # Counts is a Motif object?
        # counts = counts.array

        seq_length, A = counts.shape

        if prior is not None:
            prior = array(prior, float64)

        if prior is None or sum(prior) == 0.0:
            R = log(A)
            ent = zeros(seq_length, float64)
            entropy_interval = None
            for i in range(0, seq_length):
                C = sum(counts[i])
                # FIXME: fixup .moremath.entropy()?
                if C == 0:
                    ent[i] = 0.0
                else:
                    ent[i] = R - entropy(counts[i])
        else:
            ent = zeros(seq_length, float64)
            entropy_interval = zeros((seq_length, 2), float64)

            R = log(A)

            for i in range(0, seq_length):
                alpha = array(counts[i], float64)
                alpha += prior

                posterior = Dirichlet(alpha)
                ent[i] = posterior.mean_relative_entropy(prior / sum(prior))
                entropy_interval[i][0], entropy_interval[i][1] = \
                    posterior.interval_relative_entropy(prior / sum(prior), 0.95)

        weight = array(np.sum(counts, axis=1), float)
        max_weight = max(weight)
        if max_weight == 0.0:
            raise ValueError('No counts.')
        weight /= max_weight

        return cls(seq_length, alphabet, counts, ent, entropy_interval, weight)

    @classmethod
    def from_seqs(cls, seqs: SeqList, prior: np.ndarray = None) -> 'LogoData':
        """Build a LogoData object from a SeqList, a list of sequences."""
        # --- VALIDATE DATA ---
        # check that at least one sequence of length at least 1 long
        if len(seqs) == 0 or len(seqs[0]) == 0:
            raise ValueError("No sequence data found.")

        # Check sequence lengths
        seq_length = len(seqs[0])
        for i, s in enumerate(seqs):
            # print(i, s, len(s))
            # TODO: Redundant? Should be checked in SeqList?
            if seq_length != len(s):
                raise ArgumentError(
                        "Sequence number %d differs in length from the previous sequences"
                        % (i + 1), 'sequences')

        # FIXME: Check seqs.alphabet?

        counts = seqs.profile()
        return cls.from_counts(seqs.alphabet, counts, prior)

    def __str__(self) -> str:
        out = StringIO()
        print('## LogoData', file=out)
        print('# First column is position number, counting from zero', file=out)
        print('# Subsequent columns are raw symbol counts', file=out)
        print('# Entropy is mean entropy measured in nats.', file=out)
        print('# Low and High are the 95% confidence limits.', file=out)
        print('# Weight is the fraction of non-gap symbols in the column.', file=out)
        print('#\t', file=out)
        # Show column names
        print('#', end='\t', file=out)
        for a in self.alphabet:
            print(a, end=' \t', file=out)
        print('Entropy\tLow\tHigh\tWeight', file=out)

        # Write the data table
        for i in range(self.length):
            print(i + 1, end=' \t', file=out)
            for c in self.counts[i]:
                print(c, end=' \t', file=out)
            print("%6.4f" % self.entropy[i], end=' \t', file=out)
            if self.entropy_interval is not None:
                print("%6.4f" % self.entropy_interval[i][0], end=' \t',
                      file=out)
                print("%6.4f" % self.entropy_interval[i][1], end=' \t',
                      file=out)
            else:
                print('\t', '\t', end='', file=out)
            if self.weight is not None:
                print("%6.4f" % self.weight[i], end='', file=out)
            print('', file=out)
        print('# End LogoData', file=out)

        return out.getvalue()


def _from_URL_fileopen(target_url: str) -> TextIO:  # pragma: no cover
    """opens files from a remote URL location"""

    # parsing url in component parts
    (scheme, net_location, path, param, query, frag) = urlparse(target_url)

    # checks if string is URL link
    if scheme != "http" and scheme != "https" and scheme != "ftp":
        raise ValueError("Cannot open url: %s", target_url)

    # checks for dropbox link
    if net_location == 'www.dropbox.com':
        # changes dropbox http link into download link
        if query == "dl=0":
            query2 = "dl=1"

        # rebuild download URL, with new query2 variable
        target_url = urlunparse((scheme, net_location, path, param, query2, ""))

    # checks for google drive link
    if net_location == 'drive.google.com':
        # link configuration for direct download instead of html frame
        google_directdl_frag = "https://docs.google.com/uc?export=download&id="

        # pull file id
        (scheme, net_location, path_raw, param, query, frag) = urlparse(target_url)
        id_file = path_raw.split('/')[3]

        # rebuild URL for direct download
        target_url = google_directdl_frag + id_file

    # save url to temporary file
    req = Request(target_url)
    res = urlopen(req)
    temp = tempfile.TemporaryFile()
    shutil.copyfileobj(res, temp)
    temp.seek(0)
    return temp
