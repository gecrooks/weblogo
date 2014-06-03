#!/usr/bin/env python
 
#  Copyright (c) 2005 Gavin E. Crooks <gec@threeplusone.com>
#
#  This software is distributed under the MIT Open Source License.
#  <http://www.opensource.org/licenses/mit-license.html>
#
#  Permission is hereby granted, free of charge, to any person obtaining a 
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
#  THE SOFTWARE.
#

"""Read and write sequence information in FASTA format.
    
This is a very common format for unannotated biological sequence data, 
accepted by many multiple sequence alignment programs. Each sequence 
consists of a single-line description, followed by lines of sequence data. 
The first character of the description line is a greater-than (">") symbol 
in the first column. The first word of the description is often the name or 
ID of the sequence. Fasta files containing multiple sequences have one 
sequence listed right after another. 


Example Fasta File ::

>Lamprey GLOBIN V - SEA LAMPREY
PIVDTGSVA-P------------------LSAAEKTKIRSAWAPVYSTY---ETSGVDILVKFFTSTPAAQEFFPKFKGL
TT-----ADQLKKSA---DVRWHA-ERIINAVNDAVASMDDTEKMS--MKL-RDLSGKH----AKSFQV-----DPQYFK
VLAAVI-AD-TVAAGD--AGFEKLMSM------I---CILLR----S-----A-----Y------------
>Hagfish GLOBIN III - ATLANTIC HAGFISH
PITDHGQPP-T------------------LSEGDKKAIRESWPQIYKNF---EQNSLAVLLEFLKKFPKAQDSFPKFSAK
KS-------HLEQDP---AVKLQA-EVIINAVNHTIGLMDKEAAMK--KYL-KDLSTKH----STEFQV-----NPDMFK
ELSAVF-VS-TMG-GK--AAYEKLFSI------I---ATLLR----S-----T-----YDA----------
>Frog HEMOGLOBIN BETA CHAIN - EDIBLE FROG
----------GS-----------------------DLVSGFWGKV--DA---HKIGGEALARLLVVYPWTQRYFTTFGNL
GSADAIC-----HNA---KVLAHG-EKVLAAIGEGLKHPENLKAHY--AKL-SEYHSNK----LHVDPANFRLLGNVFIT
VLARHF-QH-EFTPELQ-HALEAHFCA------V---GDALA----K-----A-----YH-----------


"""
from __future__ import absolute_import, division, print_function

import re
from ..utils import *
from ..seq import *
from . import *


names = ( 'fasta', 'pearson', 'fa')
extensions = ('fa', 'fasta', 'fast', 'seq', 'fsa', 'fst', 'nt', 'aa','fna','mpfa', 'faa', 'fnn','mfasta','tfa', 'pfa') 


example = """
>Lamprey GLOBIN V - SEA LAMPREY
PIVDTGSVA-P------------------LSAAEKTKIRSAWAPVYSTY---ETSGVDILVKFFTSTPAAQEFFPKFKGL
TT-----ADQLKKSA---DVRWHA-ERIINAVNDAVASMDDTEKMS--MKL-RDLSGKH----AKSFQV-----DPQYFK
VLAAVI-AD-TVAAGD--AGFEKLMSM------I---CILLR----S-----A-----Y------------

>Hagfish GLOBIN III - ATLANTIC HAGFISH
PITDHGQPP-T------------------LSEGDKKAIRESWPQIYKNF---EQNSLAVLLEFLKKFPKAQDSFPKFSAK
KS-------HLEQDP---AVKLQA-EVIINAVNHTIGLMDKEAAMK--KYL-KDLSTKH----STEFQV-----NPDMFK
ELSAVF-VS-TMG-GK--AAYEKLFSI------I---ATLLR----S-----T-----YDA----------

>Frog HEMOGLOBIN BETA CHAIN - EDIBLE FROG
----------GS-----------------------DLVSGFWGKV--DA---HKIGGEALARLLVVYPWTQRYFTTFGNL
GSADAIC-----HNA---KVLAHG-EKVLAAIGEGLKHPENLKAHY--AKL-SEYHSNK----LHVDPANFRLLGNVFIT
VLARHF-QH-EFTPELQ-HALEAHFCA------V---GDALA----K-----A-----YH-----------

"""    


def read(fin, alphabet=None): 
    """Read and parse a fasta file. 

    Args:
        fin -- A stream or file to read
        alphabet -- The expected alphabet of the data, if given
    Returns: 
        SeqList -- A list of sequences
    Raises: 
        ValueError -- If the file is unparsable
    """         
    seqs = [ s for s in iterseq(fin, alphabet)]
    name = names[0]
    if hasattr(fin, "name") : name = fin.name    
    return SeqList(seqs, name=name)


def readseq(fin, alphabet=None) :
    """Read one sequence from the file, starting 
    from the current file position."""
    return next(iterseq(fin, alphabet))
    
     
def iterseq(fin, alphabet=None):
    """ Parse a fasta file and generate sequences.
    
    Args:
        fin -- A stream or file to read
        alphabet -- The expected alphabet of the data, if given    
    Yeilds: 
        Seq -- One alphabetic sequence at a time.
    Raises: 
        ValueError -- If the file is unparsable
    """
    alphabet = Alphabet(alphabet)

    seqs = []
    comments = []   # FIXME: comments before first sequence are lost.
    header = None
    header_lineno = -1
    
    def build_seq(seqs,alphabet, header, header_lineno,comments) :
        try :
            name = header.split(' ',1)[0]
            if comments :
                header += '\n' + '\n'.join(comments)
            s = Seq( "".join(seqs), alphabet, name=name, description=header)
        except ValueError:
             raise ValueError(
                "Parse failed with sequence starting at line %d: "
                "Character not in alphabet: %s" % (header_lineno, alphabet) )
        return s

    for lineno, line in enumerate(fin) :
        line = line.strip()
        if line == '' : continue
        if line.startswith('>') :
            if header is not None :
                yield  build_seq(seqs,alphabet, header, header_lineno, comments)
                header = None
                seqs = []
            header = line[1:]
            header_lineno = lineno
            comments = []
        elif line.startswith(';') : 
            # Optional (and unusual) comment line
            comments.append(line[1:])           
        else :
            if header is None :
                raise ValueError (
                    "Parse failed on line %d: sequence before header"  
                    % (lineno) )
            seqs.append(line)    

    if not seqs: return
    yield build_seq(seqs,alphabet, header, header_lineno, comments)

     
def write(fout, seqs): 
    """Write a fasta file. 

    Args:
        fout -- A writable stream.
        seqs  -- A list of Seq's
    """ 
    if seqs.description :
        for line in seqs.description.splitlines():
            print(';' + line, file=fout)
    for s in seqs:
        writeseq(fout, s)


def writeseq(afile, seq):
    """ Write a single sequence in fasta format.

    Args:
        afile -- A writable stream.
        seq  -- A Seq instance
    """
    header = seq.description or seq.name or ''
    # We prepend '>' to the first header line
    # Additional lines start with ';' to indicate comment lines
    if header:
        header = header.splitlines()
        print('>' + header[0], file=afile)
        if len(header) > 1:
            for h in header[1:]:
                print(';' + h, file=afile)
    else:
        print('>', file=afile)
    L = len(seq)
    line_length = 80
    for n in range(1 + L // line_length) :
        print(seq[n * line_length : (n+1) * line_length], file=afile)
    print(file=afile)


def index(afile, alphabet=None) :
    """Return a FileIndex for the fasta file. Sequences can be retrieved
    by item number or name.
    """
    def parser( afile) :
        return readseq(afile, alphabet)
    
    key = re.compile(r"^>\s*(\S*)")
    def linekey( line):
        k = key.search(line)
        if k is None : return None
        return k.group(1)
        
    return FileIndex(afile, linekey, parser)
