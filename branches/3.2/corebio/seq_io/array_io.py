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

"""Read and write a rectangular array of sequence data.
    
One sequence per line and nothing else. Each line must contain the same number
of characters. Blank lines and white space are ignored.
 
--- Example Array ---

--------------------------LENSTSPYDYGENESD-------FSDSPPCPQDF
--------------------------LENLEDLF-WELDRLD------NYNDTSLVENH-
--------------------------MSNITDPQMWDFDDLN-------FTGMPPADEDY
-----------------------------------YTSDN---------YSGSGDYDSNK
-SL-------NFDRTFLPALYSLLFLLGLLGNGAVAAVLLSQRTALSSTDTFLLHLAVAD
--LC-PATMASFKAVFVPVAYSLIFLLGVIGNVLVLVILERHRQTRSSTETFLFHLAVAD
-SPC-MLETETLNKYVVIIAYALVFLLSLLGNSLVMLVILYSRVGRSVTDVYLLNLALAD
-EPC-RDENVHFNRIFLPTIYFIIFLTGIVGNGLVILVMGYQKKLRSMTDKYRLHLSVAD
"""

from corebio.seq import *
from corebio.utils import *

example = """
--------------------------LENSTSPYDYGENESD-------FSDSPPCPQDF
--------------------------LENLEDLF-WELDRLD------NYNDTSLVENH-
--------------------------MSNITDPQMWDFDDLN-------FTGMPPADEDY
-----------------------------------YTSDN---------YSGSGDYDSNK
-SL-------NFDRTFLPALYSLLFLLGLLGNGAVAAVLLSQRTALSSTDTFLLHLAVAD
--LC-PATMASFKAVFVPVAYSLIFLLGVIGNVLVLVILERHRQTRSSTETFLFHLAVAD
-SPC-MLETETLNKYVVIIAYALVFLLSLLGNSLVMLVILYSRVGRSVTDVYLLNLALAD
-EPC-RDENVHFNRIFLPTIYFIIFLTGIVGNGLVILVMGYQKKLRSMTDKYRLHLSVAD
"""

names = ("array",'flatfile')
extensions = ()

def read(fin, alphabet=None): 
    """Read a file of raw sequence alignment data. 

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
    """ Read one line of sequence data and yield the sequence.

    Args:
        fin -- A stream or file to read
        alphabet -- The expected alphabet of the data, if given    
    Yeilds: 
        Seq -- One alphabetic sequence at a time.
    Raises: 
        ValueError -- If the file is unparsable
    """

    alphabet = Alphabet(alphabet)
    line_length = 0
    
    for linenum, line in enumerate(fin) :
        if line.isspace(): continue # Blank line
        line = line.strip() 

        if line[0] == '>' : # probably a fasta file. Fail.
            raise ValueError(
                "Parse Error on input line: %d " % (linenum) )
        
        line = remove_whitespace(line)
        
        if not alphabet.alphabetic(line) :
            raise ValueError(
                "Character on line: %d not in alphabet: %s : %s" % \
                     (linenum, alphabet, line) )
        
        if line_length and line_length != len(line) :
            raise ValueError("Line %d has an incommensurate length." % linenum)
        line_length = len(line)
        
        yield Seq(line, alphabet)


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
    print >>afile, seq

            
