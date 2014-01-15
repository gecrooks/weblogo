#!/usr/bin/env python

# Copyright 2005 Gavin E. Crooks <gec@threeplusone.com>
# Copyright 2005-2006 The Regents of the University of California.
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

"""Read the sequence data from a nexus file.

This IO code only gives read access to the sequence data.

Reference:
'NEXUS: An extensible file format for systematic information'
Maddison, Swofford, Maddison. 1997. Syst. Biol. 46(4):590-621
"""
from __future__ import absolute_import

from ..seq import Seq, SeqList, Alphabet
from ._nexus import Nexus, safename


names = ( 'nexus', 'paup')
extensions = ('nex', 'nexus', 'paup', 'nxs')

def iterseq(fin, alphabet=None):
    """Iterate over the sequences in the file."""
    # Default implementation
    return iter(read(fin, alphabet) )


def read(fin, alphabet=None):          
    """ Extract sequence data from a nexus file."""
    n = Nexus(fin)
    
    seqs = []
    for taxon in n.taxlabels:   
        name = safename(taxon)
        r = n.matrix[taxon]
        if alphabet is None  :
            s = Seq(r, name = name, alphabet=r.alphabet)
        else :
            s = Seq(r, name = name, alphabet=alphabet )
        seqs.append(s)

    if len(seqs) == 0 :
        # Something went terrible wrong.
        raise ValueError("Cannot parse file")
        
    return SeqList(seqs)

    
    
    
    
