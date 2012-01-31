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
from corebio.seq_io import *
from test_corebio import *

from StringIO import StringIO

import unittest
import time

class test_table_io(unittest.TestCase) :

    def test_read(self) :
        f = StringIO(table_io.example)
        seqs = table_io.read(f)
        self.assertEquals(len(seqs), 10)
        self.assertEquals(seqs[2].name, "EC0003")
        self.assertEquals(len(seqs[1]), 50)
  
    def test_read_fail(self) :
        f = StringIO(plain_io.example)
        # Wrong alphabet
        self.failUnlessRaises(ValueError, table_io.read, f)
   
    def test_write_seq(self) :
        f = StringIO(table_io.example)
        seqs = table_io.read(f)
        
        fout = StringIO()
        table_io.write(fout,seqs)
        
        fout.seek(0)
        seqs2 = table_io.read(fout)
        
        self.assertEquals(seqs, seqs2)
   
   
   
   
          
        
if __name__ == '__main__':
    unittest.main()
