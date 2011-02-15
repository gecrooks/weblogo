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



from test_corebio import *
from corebio import *
from corebio.seq import *
from corebio.seq_io import *


from StringIO import StringIO

import unittest

class test_clustal_parser(unittest.TestCase) :

    def test_parse_clustal(self) :
        f = testdata_stream("clustal.aln")
        seqs = clustal_io.read(f)
        self.assertEquals(len(seqs), 7)
        self.assertEquals(seqs[1].name, "CATH_HUMAN")
        self.assertEquals(len(seqs[1]), 395)
 
    def test_parse_clustal2_newline(self) :
        # Bug regession test. Clustal barfed on windows line endings, sometimes
        f = testdata_stream("clustalw2.aln")
        s = f.read()

        import re
        s = re.sub("\n", "\r\n", s) #Change to windows line endings
        
        seqs = clustal_io.read(StringIO(s))  
  
    def test_parse_headerless(self):
        f = testdata_stream("clustal_headerless.aln")
        seqs = clustal_io.read(f)
        self.assertEquals(len(seqs), 21)
        self.assertEquals(seqs[2].name, "O16386_CAEEL")
        self.assertEquals(len(seqs[1]), 137)
        
        

    """ Wrong alphabet should throw a parsing error """
    def test_parse_error(self) :
        f = testdata_stream("clustal.aln")
        self.failUnlessRaises(ValueError, 
            clustal_io.read, f, nucleic_alphabet )

    def test_parse_clustal181(self) :
        f = testdata_stream("clustal181.aln")
        seqs = clustal_io.read(f)

    def test_parse_clustal_glualign(self) :
        f = testdata_stream("clustal_glualign.aln")
        seqs = clustal_io.read(f, nucleic_alphabet)

    def test_parse_clustalw182(self) :
        f = testdata_stream("clustalw182.aln")
        seqs = clustal_io.read(f, protein_alphabet)

    def test_parse_fasta_fail(self) :
        # should fail with parse error
        f = StringIO(fasta_io.example)
        self.failUnlessRaises(ValueError, 
            clustal_io.read, f , protein_alphabet )
        self.failUnlessRaises(ValueError, 
            clustal_io.read, f )
                    
    def test_parse_fasta_fail2(self) :
        # should fail with parse error
        f = testdata_stream("globin.fa")
        self.failUnlessRaises(ValueError, 
            clustal_io.read, f )

    def test_parse_clustal_example(self) :
        f = StringIO(clustal_io.example)
        seqs = clustal_io.read(f)
               
    def test_write(self) :
        f = StringIO(clustal_io.example)
        seqs = clustal_io.read(f)
        
        fout = StringIO()
        clustal_io.write(fout, seqs)
        
        fout.seek(0)
        seqs2 = clustal_io.read(fout)
        
        self.assertEquals(seqs, seqs2)

  
    def test_parse_table_fail(self) :
        # should fail with parse error
        f = StringIO(table_io.example)
       
        self.failUnlessRaises(ValueError, 
            clustal_io.read, f  )



             
if __name__ == '__main__':
    unittest.main()
