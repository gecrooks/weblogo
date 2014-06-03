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

"""Read and write raw, unformatted sequence data. The whole file is read 
in as a sequence.  Whitespace is removed.


--- Example Plain/Raw/Text File ---

--------------------------LENSTSPYDYGENESD-------FSDSPPCPQDF
--------------------------LENLEDLF-WELDRLD------NYNDTSLVENH-
--------------------------MSNITDPQMWDFDDLN-------FTGMPPADEDY
-----------------------------------YTSDN---------YSGSGDYDSNK
-SL-------NFDRTFLPALYSLLFLLGLLGNGAVAAVLLSQRTALSSTDTFLLHLAVAD
--LC-PATMASFKAVFVPVAYSLIFLLGVIGNVLVLVILERHRQTRSSTETFLFHLAVAD
-SPC-MLETETLNKYVVIIAYALVFLLSLLGNSLVMLVILYSRVGRSVTDVYLLNLALAD
-EPC-RDENVHFNRIFLPTIYFIIFLTGIVGNGLVILVMGYQKKLRSMTDKYRLHLSVAD
"""
from __future__ import absolute_import, print_function

from ..seq import *
from ..utils import remove_whitespace

example = """
--------------------------LENSTSPYDYGENESD-------FSDSPPCPQDF
--------------------------LENLEDLF-WELDRLD------NYNDTSLVENH-
--------------------------MSNITDPQMWDFDDLN-------FTGMPPADEDY
-----------------------------------YTSDN---------YSGSGDYDSNK
-SL-------NFDRTFLPALYSLLFLLGLLGNGAVAAVLLSQRTALSSTDTFLLHLAVAD
--LC-PATMASFKAVFVPVAYSLIFLLGVIGNVLVLVILERHRQTRSSTETFLFHLAVAD
-SPC-MLETETLNKYVVIIAYALVFLLSLLGNSLVMLVILYSRVGRSVTDVYLLNLALAD
-EPC-RDENVHFNRIFLPTIYFIIFLTGIV
"""

names = ("plain","raw")
extensions = ()

def read(fin, alphabet=None): 
    """Read a file of raw sequence data. 

    Args:
        fin -- A stream or file to read
        alphabet -- The expected alphabet of the data, if given
    Returns: 
        SeqList -- A list of sequences
    Raises: 
        ValueError -- If the file is unparsable
    """         
    seqs = [ s for s in iterseq(fin, alphabet)]
    return SeqList(seqs)


def iterseq(fin, alphabet=None) :
    """ Read the sequence data and yield one (and only one) sequence.

    Args:
        fin -- A stream or file to read
        alphabet -- The expected alphabet of the data, if given    
    Yields: 
        Seq -- One alphabetic sequence at a time.
    Raises: 
        ValueError -- If the file is unparsable
    """

    alphabet = Alphabet(alphabet)
    lines = []    
    for linenum, line in enumerate(fin) :
        if line.isspace(): continue # Blank line
        line = line.strip() 
        
        
        if line[0] == '>' : # probable a fasta file. Fail.
                raise ValueError(
                    "Parse Error on input line: %d " % (linenum) )
        line = remove_whitespace(line)
        
        if not alphabet.alphabetic(line) :
                raise ValueError(
                    "Character on line: %d not in alphabet: %s : %s" % \
                     (linenum, alphabet, line) )
        lines.append(line)

    yield Seq(''.join(lines), alphabet)



def write(afile, seqs): 
    """Write raw sequence data, one line per sequence.

    arguments:
        afile -- A writable stream.
        seqs  -- A list of Seq's
    """         
    for s in seqs :
        writeseq(afile, s)


def writeseq(afile, seq):
    """ Write a single sequence in raw format.

    arguments:
        afile -- A writable stream.
        seq  -- A Seq instance
    """
    print(seq, file=afile)
