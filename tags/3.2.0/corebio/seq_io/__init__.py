
#  Copyright (c) 2005 Gavin E. Crooks <gec@threeplusone.com>
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





""" Sequence file reading and writing.

Biological sequence data is stored and transmitted using a wide variety of
different file formats. This package provides convenient methods to read and
write several of these file fomats.

CoreBio is often capable of guessing the correct file type, either from the
file extension or the structure of the file:
>>> import corebio.seq_io
>>> afile = open("test_corebio/data/cap.fa")
>>> seqs = corebio.seq_io.read(afile)

Alternatively, each sequence file type has a separate module named FILETYPE_io
(e.g. fasta_io, clustal_io).
>>> import corebio.seq_io.fasta_io
>>> afile = open("test_corebio/data/cap.fa")
>>> seqs = corebio.seq_io.fasta_io.read( afile )

Sequence data can also be written back to files:
>>> fout = open("out.fa", "w")
>>> corebio.seq_io.fasta_io.write( fout, seqs )


Supported File Formats
----------------------

Module              Name            Extension  read write features   
---------------------------------------------------------------------------
array_io            array, flatfile             yes  yes    none
clustal_io          clustalw        aln         yes  yes
fasta_io            fasta, Pearson  fa          yes  yes    none
genbank_io          genbank         gb          yes         
intelligenetics_io  intelligenetics ig          yes  yes
msf_io              msf             msf         yes
nbrf_io             nbrf, pir       pir         yes
nexus_io            nexus           nexus       yes
phylip_io           phylip          phy         yes
plain_io            plain, raw                  yes  yes    none
table_io            table           tbl         yes  yes    none

Each IO module defines one or more of the following functions and variables:

read(afile, alphabet=None) 
    Read a file of sequence data and return a SeqList, a collection
    of Seq's (Alphabetic strings) and features.

read_seq(afile, alphabet=None)
    Read a single sequence from a file.

iter_seq(afile, alphabet =None) 
    Iterate over the sequences in a file. 
    
index(afile, alphabet = None)
    Instead of loading all of the sequences into memory, scan the file and
    return an index map that will load sequences on demand. Typically not
    implemented for formats with interleaved sequences.

write(afile, seqlist)
    Write a collection of sequences to the specifed file.

write_seq(afile, seq)
    Write one sequence to the file. Only implemented for non-interleaved, 
    headerless formats, such as fasta and plain.

example
    A string containing a short example of the file format

names
    A list of synonyms for the file format. E.g. for fasta_io, ( 'fasta',    
    'pearson', 'fa'). The first entry is the preferred format name.

extensions
    A list of file name extensions used for this file format. e.g. 
    fasta_io.extensions is ('fa', 'fasta', 'fast', 'seq', 'fsa', 'fst', 'nt',
    'aa','fna','mpfa').  The preferred or standard extension is first in the 
    list.


Attributes :
- formats -- Available seq_io format parsers
- format_names -- A map between format names and format parsers.
- format_extensions -- A map between filename extensions and parsers.

"""

# Dev. References :
#
#    - http://iubio.bio.indiana.edu/soft/molbio/readseq/java/Readseq2-help.html
#    - http://www.ebi.ac.uk/help/formats_frame.html
#    - http://www.cmbi.kun.nl/bioinf/tools/crab_pir.html
#    - http://bioperl.org/HOWTOs/html/SeqIO.html
#    - http://emboss.sourceforge.net/docs/themes/SequenceFormats.html
#    - http://www.cse.ucsc.edu/research/compbio/a2m-desc.html (a2m)
#    - http://www.genomatix.de/online_help/help/sequence_formats.html

from corebio.seq import *

import clustal_io
import fasta_io
import msf_io
import nbrf_io
import nexus_io
import plain_io
import phylip_io
#import null_io
import stockholm_io
import intelligenetics_io
import table_io
import array_io
import genbank_io

__all__ = [
    'clustal_io', 
    'fasta_io',
    'msf_io',
    'nbrf_io',
    'nexus_io',
    'plain_io',
    'phylip_io',
    'null_io',
    'stockholm_io',
    'intelligenetics_io',
    'table_io',
    'array_io',
    'genbank_io',
    'read',
    'formats',
    'format_names',
    'format_extensions',
    ]

formats = ( clustal_io, fasta_io, plain_io, msf_io, genbank_io,nbrf_io, nexus_io, phylip_io, stockholm_io, intelligenetics_io, table_io, array_io)
"""Available seq_io formats"""


def format_names() :   
    """Return a map between format names and format modules"""
    global formats
    fnames = {}
    for f in formats :
        for name in f.names :
            assert name not in fnames # Insanity check
            fnames[name] = f  
    return fnames

def format_extensions() :   
    """Return a map between filename extensions and sequence file types"""
    global formats
    fext = {}
    for f in formats :
        for ext in f.extensions :
            assert ext not in fext # Insanity check
            fext[ext] = f  
    return fext

    
# seq_io._parsers is an ordered list of sequence parsers that are tried, in 
# turn, on files of unknown format. Each parser must raise an exception when
# fed a format further down the list.
#
# The general trend is most common to least common file format. However, 
# 'nbrf_io' is before 'fasta_io' because nbrf looks like fasta with extras, and
# 'array_io' is last, since it is very general.
_parsers = (nbrf_io, fasta_io, clustal_io, phylip_io, genbank_io, stockholm_io, msf_io, nexus_io, table_io, array_io)

 
def _get_parsers(fin) :
    global _parsers
    
    fnames = format_names()
    fext = format_extensions()
    parsers = list(_parsers)
    best_guess = parsers[0]
    
    # If a filename is supplied use the extension to guess the format.
    if hasattr(fin, "name") and '.' in fin.name :
        extension = fin.name.split('.')[-1]
        if extension in  fnames:
            best_guess = fnames[extension]
        elif extension in fext :
            best_guess = fext[extension]
        
    if best_guess in parsers :
        parsers.remove(best_guess)
    parsers.insert(0,best_guess)

    return parsers

 
    
def read(fin, alphabet=None) :
    """ Read a sequence file and attempt to guess its format. 
    First the filename extension (if available) is used to infer the format.
    If that fails, then we attempt to parse the file using several common   
    formats.
    
    returns :
        SeqList
    raises :
        ValueError - If the file cannot be parsed.
        ValueError - Sequence do not conform to the alphabet.
    """

    alphabet = Alphabet(alphabet)
    parsers =  _get_parsers(fin)
    
    for p in _get_parsers(fin) :
        try:    
            return p.read(fin, alphabet)
        except ValueError:
            pass
        fin.seek(0)             # FIXME. Non seakable stdin? 
            
    names = ", ".join([ p.names[0] for p in parsers])
    raise ValueError("Cannot parse sequence file: Tried %s " % names)






