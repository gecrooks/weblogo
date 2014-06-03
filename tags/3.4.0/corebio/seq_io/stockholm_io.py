 
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

"""Read a multiple sequence alignment in STOCKHOLM format.

This file format is used by PFAM and HMMER. At present, all annotation
information is ignored.

See: 
    - http://www.cgb.ki.se/cgb/groups/sonnhammer/Stockholm.html
    - HMMER manual

"""
from __future__ import absolute_import, print_function

import re

from ..utils import *
from ..seq import *
from . import *



example = """
# STOCKHOLM 1.0
#=GF ID CBS
#=GF AC PF00571
#=GF DE CBS domain
#=GF AU Bateman A
#=GF CC CBS domains are small intracellular modules mostly found  
#=GF CC in 2 or four copies within a protein. 
#=GF SQ 67
#=GS O31698/18-71 AC O31698
#=GS O83071/192-246 AC O83071
#=GS O83071/259-312 AC O83071
#=GS O31698/88-139 AC O31698
#=GS O31698/88-139 OS Bacillus subtilis
O83071/192-246          MTCRAQLIAVPRASSLAE..AIACAQKM....RVSRVPVYERS
#=GR O83071/192-246 SA  999887756453524252..55152525....36463774777
O83071/259-312          MQHVSAPVFVFECTRLAY..VQHKLRAH....SRAVAIVLDEY
#=GR O83071/259-312 SS  CCCCCHHHHHHHHHHHHH..EEEEEEEE....EEEEEEEEEEE
O31698/18-71            MIEADKVAHVQVGNNLEH..ALLVLTKT....GYTAIPVLDPS
#=GR O31698/18-71 SS    CCCHHHHHHHHHHHHHHH..EEEEEEEE....EEEEEEEEHHH
O31698/88-139           EVMLTDIPRLHINDPIMK..GFGMVINN......GFVCVENDE
#=GR O31698/88-139 SS   CCCCCCCHHHHHHHHHHH..HEEEEEEE....EEEEEEEEEEH
#=GC SS_cons            CCCCCHHHHHHHHHHHHH..EEEEEEEE....EEEEEEEEEEH
O31699/88-139           EVMLTDIPRLHINDPIMK..GFGMVINN......GFVCVENDE
#=GR O31699/88-139 AS   ________________*__________________________
#=GR_O31699/88-139_IN   ____________1______________2__________0____
//
"""



names = ("stockholm", "pfam",)
extensions = ('sth', 'stockholm', 'align')


header_line = re.compile(r'#\s+STOCKHOLM\s+1.\d\s+$')

def iterseq(fin, alphabet=None):
    """Iterate over the sequences in the file."""
    # Default implementation
    return iter(read(fin, alphabet) )


def read(fin, alphabet=None) :  
    alphabet = Alphabet(alphabet)      
    seq_ids = []
    seqs = []
    block_count = 0
    
    
    for token in _scan(fin):
        if token.typeof== "begin_block":
            block_count = 0
        elif token.typeof == "seq_id":
            if len(seqs) <= block_count :
                seq_ids.append(token.data)
                seqs.append([])
        elif token.typeof == "seq":
            if not alphabet.alphabetic(token.data) :
                raise ValueError(
                    "Character on line: %d not in alphabet: %s : %s" % (
                    token.lineno, alphabet, token.data) )
            seqs[block_count].append(token.data)
            block_count +=1

          
    seqs = [ Seq("".join(s), alphabet, name= i) for s,i in zip(seqs,seq_ids)]
    return SeqList(seqs)


def _scan( fin ):

    header, body, block = range(3)
    
    yield Token("begin")
    state = header
    for L, line in enumerate(fin):
        

        if state==header :
            if line.isspace() : continue
            m = header_line.match(line)
            state = body
            if m is not None :
                # print("header: ", m.group())
                yield Token("header", m.group() )
                continue
            else :
                 raise ValueError("Parse error on line: %d" % L)
        
        
        if state == body :
            if line.isspace() : continue
            yield Token("begin_block")
            state = block
            # fall through to block
        
        if state ==  block:        
            if line.isspace() :
                yield Token("end_block")
                state = body
                continue
            if line.strip() == '//' :
                yield Token("end_block")
                return
            
            
            if line[0] =='#' :  # Comment or annotation line 
                continue
   
            name_seq = line.split(None,1) # Split into two parts at first whitespace 
            if len(name_seq) != 2 : 
                raise ValueError("Parse error on line: %d" % L)
 
            
            yield Token("seq_id", name_seq[0].strip() )
            yield Token("seq", name_seq[1].strip() )
            continue

        # END state blocks. If I ever get here something has gone terrible wrong
        raise RuntimeError()

