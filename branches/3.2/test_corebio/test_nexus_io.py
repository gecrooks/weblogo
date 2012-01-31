#!/usr/bin/env python


from corebio import *
from corebio.seq import *
from corebio.seq_io import *
from corebio.utils import *
from test_corebio import *

from StringIO import StringIO

import unittest

class test_fasta_io(unittest.TestCase) :

    def test_read(self) :
        f = testdata_stream("nexus/protein.nex")  
        seqs = nexus_io.read(f)
        #print seqs
        self.assertEquals(len(seqs), 10)
        self.assertEquals(seqs[0].name, "Cow")
        self.assertEquals(len(seqs[1]), 234)
        self.assertEquals( str(seqs[0][0:10]), 'MAYPMQLGFQ')
        
  
   
    def test_parse_StringIO(self) :
        # Bio.Nexus cannot read from a StringIO object.
        f0 = testdata_stream("nexus/protein.nex")
        f = StringIO(f0.read() )
        n = nexus_io.read(f)  
          
    def test_parse_fasta_fail(self) :
        f = testdata_stream("globin.fa")
        self.failUnlessRaises(ValueError, 
            nexus_io.read, f  )
            

    def test_parse_clustal_fail(self) :
        # should fail with parse error
        f = StringIO(clustal_io.example)
        self.failUnlessRaises(ValueError, 
            nexus_io.read, f , protein_alphabet )
   
    def test_parse_plain_fail(self) :
        # should fail with parse error
        f = StringIO(plain_io.example)
        self.failUnlessRaises(ValueError, 
            nexus_io.read, f  )
   
             
if __name__ == '__main__':
    unittest.main()
