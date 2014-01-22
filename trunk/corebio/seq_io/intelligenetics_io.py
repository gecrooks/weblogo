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

"""Read and write sequence information in IntelliGenetics format.
    
A sequence file in IG format can contain several sequences, each consisting of a
number of comment lines that must begin with a semicolon (";"), a line with the 
sequence name and the sequence itself terminated with the termination character 
'1' for linear or '2' for circular sequences. The termination character is
defacto optional.

--- Example IG File ---

;H.sapiens fau mRNA, 518 bases
HSFAU
ttcctctttctcgactccatcttcgcggtagctgggaccgccgttcagtc
actcttaagtcttttgtaattctggctttctctaataaaaaagccactta
gttcagtcaaaaaaaaaa1
;H.sapiens fau 1 gene, 2016 bases
HSFAU1
ctaccattttccctctcgattctatatgtacactcgggacaagttctcct
gatcgaaaacggcaaaactaaggccccaagtaggaatgccttagttttcg
gggttaacaatgattaacactgagcctcacacccacgcgatgccctcagc
tcctcgctcagcgctctcaccaacagccgtagcccgcagccccgctggac
accggttctccatccccgcagcgtagcccggaacatggtagctgccatct
ttacctgctacgccagccttctgtgcgcgcaactgtctggtcccgcccc2

"""
from __future__ import absolute_import, division, print_function

from ..utils import *
from ..seq import *
from . import *


names = ( 'intelligenetics', 'ig', 'stanford', )
extensions = ('ig') 


example = """
;H.sapiens fau mRNA, 518 bases
HSFAU
ttcctctttctcgactccatcttcgcggtagctgggaccgccgttcagtc
actcttaagtcttttgtaattctggctttctctaataaaaaagccactta
gttcagtcaaaaaaaaaa1
;H.sapiens fau 1 gene, 2016 bases
HSFAU1
ctaccattttccctctcgattctatatgtacactcgggacaagttctcct
gatcgaaaacggcaaaactaaggccccaagtaggaatgccttagttttcg
gggttaacaatgattaacactgagcctcacacccacgcgatgccctcagc
tcctcgctcagcgctctcaccaacagccgtagcccgcagccccgctggac
accggttctccatccccgcagcgtagcccggaacatggtagctgccatct
ttacctgctacgccagccttctgtgcgcgcaactgtctggtcccgcccc2
"""    




def read(fin, alphabet=None): 
    """Read and parse an IG file. 

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
    """ Parse an IG file and generate sequences.
    
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
    header = []   
    start_lineno = -1
    name = None
    
    def build_seq(seqs,alphabet, name, comments, lineno) :
        try :
            desc = '\n'.join(comments)
            s = Seq( "".join(seqs), alphabet, name=name, description=desc)
        except ValueError :
             raise ValueError(
                "Parse failed with sequence starting at line %d: "
                "Character not in alphabet: %s" % (lineno, alphabet) )
        return s

    for lineno, line in enumerate(fin) :
        line = line.strip()
        if line == '' : continue
        if line.startswith(';') :
            if seqs :
                # end of sequence
                yield  build_seq(seqs,alphabet, name, header, start_lineno)
                header = []
                seqs = []
                name = None
            header.append(line[1:])
            start_lineno = lineno
        elif not name :  
            name = line
        elif line[-1] == '1' or line[-1]=='2': 
            # End of sequence
            seqs.append(remove_whitespace(line[0:-1]))
            yield  build_seq(seqs,alphabet, name, header, start_lineno)
            header = []
            seqs = []       
            name = None
        else:
            seqs.append( remove_whitespace(line))    
    
    if seqs :
        yield build_seq(seqs,alphabet, name, header, start_lineno)
    return
    
    
    


def write(fout, seqs): 
    """Write an IG file. 

    Args:
        fout -- A writable stream.
        seqs  -- A list of Seq's
    Raises:
        ValueError -- If a sequence is missing a name
    """        
    for s in seqs :
        writeseq(fout, s)


def writeseq(fout, seq):
    """ Write a single sequence in IG format.

    Args:
        afile -- A writable stream.
        seq  -- A Seq instance
    Raises:
        ValueError -- If a sequence is missing a name        
    """
    desc = seq.description or ''
    # We prepend ';' to each line
    for h in desc.splitlines():
        print(';' + h, file=fout)
    if not seq.name:
        raise ValueError(
            "Write failed with missing sequence name: %s"% str(seq) )
    print(seq.name, file=fout)
    L = len(seq)
    line_length = 80
    for n in range(1 + L // line_length):
        print(seq[n * line_length : (n+1) * line_length], file=fout)
    print(file=fout)
