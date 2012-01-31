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


from corebio import *
from corebio.seq import *

import corebio.seq_io 
from corebio.seq_io import genbank_io

from test_corebio import *

from StringIO import StringIO

import unittest
import time

def examples():
    return (
        testdata_stream('genbank/NT_019265.gb'),            
        testdata_stream('genbank/cox2.gb'),  
        testdata_stream('genbank/iro.gb'),                  
        testdata_stream('genbank/pri1.gb'),
        testdata_stream('genbank/dbsource_wrap.gb'),        
        testdata_stream('genbank/noref.gb'),                
        testdata_stream('genbank/protein_refseq.gb'),
        testdata_stream('genbank/cor6_6.gb'),               
        testdata_stream('genbank/origin_line.gb'),
        
        #These files are too large to include in the distribution
        #testdata_stream('genbank/arab1.gb'),
        #testdata_stream('genbank/NC_005213.gbk'),
        #testdata_stream('genbank/NC_003888.gbk'),
    )
        



class test_genbank_io(unittest.TestCase) :
    
    # Useful for debugging
    #def test_scan(self) :
    #    for f in examples():
    #        for t in genbank_io._scan(f):
    #            print t
    #        print
    #        print 
            
    def test_parse(self) :
        for f in examples():
            #print f.name
            seqs = genbank_io.read(f)
            #print seqs
                        
    def test_read(self) :
        f = testdata_stream("genbank/cox2.gb")
        seqs = genbank_io.read(f)
        #print seqs

        self.assertEquals(len(seqs), 5)
        self.assertEquals(len(seqs[1]), 210)
        
        f.seek(0)
        seqs = seq_io.read(f)
        self.assertEquals(len(seqs), 5)
        self.assertEquals(len(seqs[1]), 210)
        

         
        f = testdata_stream('genbank/NT_019265.gb') 
        seqs = genbank_io.read(f)
        self.assertEquals(len(seqs), 1)
        self.assertEquals(len(seqs[0]), 0)
        


        


if __name__ == '__main__':
    unittest.main()
