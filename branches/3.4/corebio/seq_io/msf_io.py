#!/usr/bin/env python

#  Copyright (c) 2005 Clare Gollnick <cgollnick@berkeley.edu>
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


"""Read sequence information in MSF format.
    
This is a file format for biological sequence data. The sequences are interweaved and each line is labeled with the sequence name. The MSF format can be identified in one or more of the following ways: 
1. The word PileUp on the first line (optional)
2. the word !!AA_MULTIPLE_ALIGNMENT or !!NA_MULTIPLE_ALIGNMENT at the start of the file (optional)
3. the word MSF on the first line of the file, and the characters ".." at the end of this line (optional)
4. A header containing sequence information followed by a line with the characters "//"
"""
from __future__ import absolute_import

import re

from ..seq import *
from ..utils import *
from . import *


example= """

 PileUp


MSF: 64 Type: P Check: 767 ..

 Name: Cow              Len:   100  Check: 3761  Weight:  1.00
 Name: Carp             Len:   100  Check: 1550  Weight:  1.00
 Name: Chicken          Len:   100  Check: 2397  Weight:  1.00
 Name: Human            Len:   100  Check: 9021  Weight:  1.00
 Name: Loach            Len:   100  Check:  984  Weight:  1.00
 Name: Mouse            Len:   100  Check: 2993  Weight:  1.00


//

                                                            
    Cow  MAYPMQLGFQ DATSPIMEEL LHFHDHTLMI VFLISSLVLY IISLMLTTKL 
   Carp  MAHPTQLGFK DAAMPVMEEL LHFHDHALMI VLLISTLVLY IITAMVSTKL 
Chicken  MANHSQLGFQ DASSPIMEEL VEFHDHALMV ALAICSLVLY LLTLMLMEKL 
  Human  MAHAAQVGLQ DATSPIMEEL ITFHDHALMI IFLICFLVLY ALFLTLTTKL 
  Loach  MAHPTQLGFQ DAASPVMEEL LHFHDHALMI VFLISALVLY VIITTVSTKL 
  Mouse  MAYPFQLGLQ DATSPIMEEL MNFHDHTLMI VFLISSLVLY IISLMLTTKL 


                                                       
    Cow  THTSTMDAQE VETIWTILPA IILILIALPS LRILYMMDEI NNPSLTVKTM 
   Carp  TNKYILDSQE IEIVWTILPA VILVLIALPS LRILYLMDEI NDPHLTIKAM 
Chicken  S.SNTVDAQE VELIWTILPA IVLVLLALPS LQILYMMDEI DEPDLTLKAI 
  Human  TNTNISDAQE METVWTILPA IILVLIALPS LRILYMTDEV NDPSLTIKSI 
  Loach  TNMYILDSQE IEIVWTVLPA LILILIALPS LRILYLMDEI NDPHLTIKAM 
  Mouse  THTSTMDAQE VETIWTILPA VILIMIALPS LRILYMMDEI NNPVLTVKTM 
 
   """

names = ('msf', 'gcg-msf', 'gcg', 'PileUp')
extensions = ('msf')

end_header=re.compile(r'(//)(\s*)$')
seq_line=re.compile(r'\s*(\S+)\s+([\S\s.?]+)$')

def iterseq(fin, alphabet=None):
    """Iterate over the sequences in the file."""
    # Default implementation
    return iter(read(fin, alphabet) )


        
def read(fin, alphabet=None):
    alphabet =Alphabet(alphabet)
    seq_ids=[]
    seqs=[]
    block_count=0

    for token in _line_is(fin):
        if token.typeof=="begin_block":
                block_count=0
            
        elif token.typeof == "seq_id":
            if len(seqs)<= block_count:
                seq_ids.append(token.data)
                seqs.append([])
        elif token.typeof=="seq":
            if not alphabet.alphabetic(token.data):
                raise ValueError(
                    "Character on line: %d not in alphabet: %s : %s" % (
                    token.lineno, alphabet, token.data) ) 
            seqs[block_count].append(token.data)
            block_count +=1
    if seq_ids==[]:
            raise ValueError("Parse error, possible wrong format")
    seqs = [ Seq("".join(s), alphabet, name= i) for s,i in zip(seqs,seq_ids)]
    return SeqList(seqs)
    
def _line_is(fin):
    header, body, block = range(3) 
    yield Token("begin")
    state=header
    for L, line in enumerate(fin):
        if state==header:
            if line.isspace():continue
            m=end_header.match(line)
            if m is not None:
                yield Token("end_header")
                state=body
                continue            
            else: continue
            
        if state==body:
            if line.isspace():continue
            yield Token("begin_block")
            state=block
            #skips to a block of sequences
            
        if state==block:
            if line.isspace():
                yield Token("end_block") 
                state=body 
                continue
            m=seq_line.match(line)
            if m is None:
                raise ValueError("Parse error on line: %d" % L)
            if m.group(1).isdigit() and m.group(2).strip().isdigit():
                continue
            yield Token("seq_id",m.group(1).strip() )
            data=m.group(2)
            data="".join((data.split()))
            yield Token("seq",data.strip() )

            
            
            

                
