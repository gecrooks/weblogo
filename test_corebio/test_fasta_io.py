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

example_with_optional_comments="""
>SEQUENCE_1
;comment line 1 (optional)
MTEITAAMVKELRESTGAGMMDCKNALSETNGDFDKAVQLLREKGLGKAAKKADRLAAEG
LVSVKVSDDFTIAAMRPSYLSYEDLDMTFVENEYKALVAELEKENEERRRLKDPNKPEHK
IPQFASRKQLSDAILKEAEEKIKEELKAQGKPEKIWDNIIPGKMNSFIADNSQLDSKLTL
MGQFYVMDDKKTVEQVIAEKEKEFGGKIKIVEFICFEVGEGLEKKTEDFAAEVAAQL
>SEQUENCE_2
;comment line 1 (optional)
;comment line 2 (optional)
SATVSEINSETDFVAKNDQFIALTKDTTAHIQSNSLQSVEELHSSTINGVKFEEYLKSQI
ATIGENLVVRRFATLKAGANGVVNGYIHTNGRVGVVIAAACDSAEVASKSRDLLRQICMH
"""

example3="""
>
AAAGTG
>
AAAGCG
>
TGCCCT
>
TGCCTT
"""

example4="""
>
AAAGTG
>
AAAGCG

TGCCCT
>
TGCCTT
"""






class test_fasta_io(unittest.TestCase) :

    def test_read(self) :
        f = StringIO(fasta_io.example)
        seqs = fasta_io.read(f)
        #print seqs
        self.assertEquals(len(seqs), 3)
        self.assertEquals(seqs[0].description, "Lamprey GLOBIN V - SEA LAMPREY")
        self.assertEquals(seqs[0].name, "Lamprey")
        self.assertEquals(len(seqs[1]), 231)

#    def test_read_long(self):
#        f = testdata_stream("NC_000913.ffn")
#        count = 0
#        start = time.time()
#        for s in seq_io.fasta_io.read_seq(f):
#            count +=1
#        end = time.time()
#        t = end-start
#        
#        self.assertEquals(count, 4243)
#        
#        # Timing is 3s 1.67 GHz G4
#        # print t
  
    def test_read_fail(self) :
        f = StringIO(fasta_io.example)
        # Wrong alphabet
        self.failUnlessRaises(ValueError, fasta_io.read, f, nucleic_alphabet)
          
    def test_parse_globin(self) :
        #f = open_resource(__file__, "test_data","globin.fa")
        f = testdata_stream("globin.fa")
        seqs = fasta_io.read(f,protein_alphabet)
        self.assertEquals(len(seqs), 56)


    def test_parse_clustal_fail(self) :
        # should fail with parse error
        f = StringIO(clustal_io.example)
        self.failUnlessRaises(ValueError, 
            fasta_io.read, f , protein_alphabet )
   
    def test_parse_plain_fail(self) :
        # should fail with parse error
        f = StringIO(plain_io.example)
        self.failUnlessRaises(ValueError, 
            fasta_io.read, f  )
   
   
    def test_write_seq(self) :
        f = StringIO(fasta_io.example)
        seqs = fasta_io.read(f)
        fout = StringIO()
        fasta_io.write(fout,seqs)
        
        fout.seek(0)
        seqs2 = fasta_io.read(fout)
        
        self.assertEquals(seqs, seqs2)
    
    def test_write_with_header(self) :
        f = StringIO(fasta_io.example)
        seqs = fasta_io.read(f)
        seqs.description = 'A description\nMore description'
        fout = StringIO()
        fasta_io.write(fout,seqs)


   
    def test_read_comments(self)   :
        f = StringIO(example_with_optional_comments)
        seqs = fasta_io.read(f)             
        self.assertEquals(len(seqs), 2)
        self.assertEquals( seqs[1].startswith("SATVSEI"), True)
        self.assertEquals( seqs[1].description.splitlines()[1] , 
                ("comment line 1 (optional)"))
        
    def test_write_comments(self)   :
        f = StringIO(example_with_optional_comments)
        seqs = fasta_io.read(f)
        fout = StringIO()
        fasta_io.write(fout, seqs)
        fout.seek(0)
        seqs2 = fasta_io.read(fout)
        self.assertEquals(seqs, seqs2)
        
        self.assertEquals(seqs[1].description, seqs2[1].description)
   
    def test_read_headerless(self) :
        # This example has blank headers.
        f = StringIO(example3)
        seqs = fasta_io.read(f)             
        self.assertEquals(len(seqs), 4)
        #print seqs
        
        fout = StringIO()
        fasta_io.write(fout, seqs)


    def test_index(self) :
        f = StringIO(fasta_io.example)
        idx = fasta_io.index(f)
        #print idx._key_dict
        self.assertEquals(len(idx), 3)
        self.assertEquals(idx[0].description, "Lamprey GLOBIN V - SEA LAMPREY")
        self.assertEquals(idx[0].name, "Lamprey")
        self.assertEquals(idx['Lamprey'].name, "Lamprey")
        self.assertEquals(len(idx['Hagfish']), 231)

    def test_read_empty(self) :
        f = StringIO()
        seqs = fasta_io.read(f) 
        assert len(seqs)==0

    def test_isaligned(self):
        seqs = fasta_io.read(StringIO())
        assert seqs.isaligned()
        seqs = fasta_io.read(StringIO(fasta_io.example))
        assert seqs.isaligned()
        seqs = fasta_io.read(StringIO(example4))
        assert not seqs.isaligned()

    def test_read_with_blank_line(self):
        f=StringIO(example4)
        seqs = fasta_io.read(f)
        assert not seqs.isaligned()
        self.assertEquals(len(seqs), 3)
        
        
        
if __name__ == '__main__':
    unittest.main()
