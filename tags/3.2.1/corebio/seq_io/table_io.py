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

"""Read and write sequence information in tab-delimited format.

This very simple format has two columns per line. The first column is a sequence name, the second column is the sequence itself. The columns are separated by a single tab ("\\t") character.

"""
from corebio.utils import *
from corebio.seq import *
from corebio.seq_io import *


names = ( 'table', 'tab')
extensions = ('tbl')


example = """
EC0001	MKRISTTITTTITITTGNGAG
EC0002	MRVLKFGGTSVANAERFLRVADILESNARQGQVATVLSAPAKITNHLVAM
EC0003	MVKVYAPASSANMSVGFDVLGAAVTPVDGALLGDVVTVEAAETFSLNNLG
EC0004	MKLYNLKDHNEQVSFAQAVTQGLGKNQGLFFPHDLPEFSLTEIDEMLKLD
EC0005	MKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGH
EC0006	MLILISPAKTLDYQSPLTTTRYTLPELLDNSQQLIHEARKLTPPQISTLM
EC0007	MPDFFSFINSVLWGSVMIYLLFGAGCWFTFRTGFVQFRYIRQFGKSLKNS
EC0008	MTDKLTSLRQYTTVVADTGDIAAMKLYQPQDATTNPSLILNAAQIPEYRK
EC0009	MNTLRIGLVSISDRASSGVYQDKGIPALEEWLTSALTTPFELETRLIPDE
EC0010	MGNTKLANPAPLGLMGFGMTTILLNLHNVGYFALDGIILAMGIFYGGIAQ
"""    




def read(fin, alphabet=None): 
    """Read and parse file. 

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

    
def iterseq(fin, alphabet=None):
    """ Parse a file and generate sequences.
    
    Args:
        fin -- A stream or file to read
        alphabet -- The expected alphabet of the data, if given    
    Yeilds: 
        Seq -- One alphabetic sequence at a time.
    Raises: 
        ValueError -- If the file is unparsable
    """
    alphabet = Alphabet(alphabet)

    for lineno, line in enumerate(fin) :
        line = line.strip()
        if line == '' : continue

        columns = line.split('\t')
        if len(columns) !=2 :
            raise ValueError( "Parse failed on line %d: did not find two "
             "columns separated by a tab."  % (lineno) )        
        yield Seq(columns[1], alphabet=alphabet, name=columns[0])
     
     
def write(fout, seqs): 
    """Write a two column, tab-delimited file. 

    Args:
        fout -- A writable stream.
        seqs  -- A list of Seq's
    """      
    for s in seqs : writeseq(fout, s)

    
def writeseq(fout, seq):
    """ Write a single sequence in fasta format.

    Args:
        afile -- A writable stream.
        seq  -- A Seq instance
    """
    
    name = seq.name or ''
    print >>fout, name, '\t', seq

    
    
    
    
     


    
    
    
