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
from corebio.seq_io import array_io

from StringIO import StringIO

from test_corebio import *
import unittest

class test_array_io(unittest.TestCase) :

    def test_read_example(self) :
        f = StringIO(array_io.example)
        seqs = array_io.read(f)
        #print seqs
        self.assertEquals(len(seqs), 8)
        self.assertEquals(seqs[0].name, None)
        self.assertEquals(len(seqs[1]), 60)
  
    def test_write_seq(self) :
        f = StringIO(array_io.example)
        seqs = array_io.read(f)
        
        fout = StringIO()
        array_io.write(fout,seqs)
        
        fout.seek(0)
        seqs2 = array_io.read(fout)
        
        self.assertEquals(seqs, seqs2)  

    def test_fail(self) :
        # Lengths differ
        example = """
-SPC-MLETETLNKYVVIIAYALVFLLSLLGNSLVMLVILYSRVGRSVTDVYLLNLALAD
-EPC-RDEN
"""
        f = StringIO(example)
        self.failUnlessRaises(ValueError, array_io.read, f)

        
        


             
if __name__ == '__main__':
    unittest.main()
