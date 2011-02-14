#!/usr/bin/env python
 
#  Copyright (c) 2006, The Regents of the University of California, through 
#  Lawrence Berkeley National Laboratory (subject to receipt of any required
#  approvals from the U.S. Dept. of Energy).  All rights reserved.

#  This software is distributed under the new BSD Open Source License.
#  <http://www.opensource.org/licenses/bsd-license.html>
#
#  Redistribution and use in source and binary forms, with or without 
#  modification, are permitted provided that the following conditions are met: 
#
#  (1) Redistributions of source code must retain the above copyright notice, 
#  this list of conditions and the following disclaimer. 
#
#  (2) Redistributions in binary form must reproduce the above copyright 
#  notice, this list of conditions and the following disclaimer in the 
#  documentation and or other materials provided with the distribution. 
#
#  (3) Neither the name of the University of California, Lawrence Berkeley 
#  National Laboratory, U.S. Dept. of Energy nor the names of its contributors 
#  may be used to endorse or promote products derived from this software 
#  without specific prior written permission. 
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
#  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
#  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
#  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
#  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
#  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
#  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
#  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
#  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
#  POSSIBILITY OF SUCH DAMAGE. 


from corebio import *
from corebio.seq import *
from corebio.seq_io import *
from corebio.seq_io import stockholm_io
from test_corebio import *

from StringIO import StringIO

import unittest

class test_stockholm_io(unittest.TestCase) :

    def test_parse1(self) :
        f = testdata_stream("pfam.txt")
        seqs = stockholm_io.read(f)
        #self.assertEquals(len(seqs), 7)
        self.assertEquals(seqs[1].name, "O61132/1-232")
        self.assertEquals(len(seqs[1]), 265)
  
    def test_parse2(self) :
        f = StringIO(stockholm_io.example)
        seqs = stockholm_io.read(f)
        self.assertEquals(len(seqs), 5)
        self.assertEquals(seqs[1].name, "O83071/259-312")
        self.assertEquals(len(seqs[1]), 43)
  
#12345678901234567890123456789012345678901234567890
#12*50 +6 = 606
#QYVTVFYGVPAWRNATIPLFCATKNR.......DTWGTTQCLPDNDDYSE
  
  
    def test_parse3(self) :
        f = testdata_stream("pfam_example.txt")
        seqs = stockholm_io.read(f)
        self.assertEquals(len(seqs), 24)
        self.assertEquals(seqs[5].name, "ENV_HV2BE/24-510")
        self.assertEquals(len(seqs[1]), 606)
        self.assertEquals( str(seqs[0][-6:]) , 'TSRNKR')
  
  
    def test_parse_error(self) :
        """ Wrong alphabet should throw a parsing error """
        f = StringIO(stockholm_io.example)
        self.failUnlessRaises(ValueError, 
            clustal_io.read, f, nucleic_alphabet )

    def test_parse_fasta_fail(self) :
        # should fail with parse error
        f = StringIO(fasta_io.example)
        self.failUnlessRaises(ValueError, 
            stockholm_io.read, f , protein_alphabet )
        
    def test_parse_fasta_fail2(self) :
        # should fail with parse error
        f = testdata_stream("globin.fa")
        self.failUnlessRaises(ValueError, 
            stockholm_io.read, f )


    def test_parse_fail(self) :
        # should fail with parse error
        examples = (
            StringIO(clustal_io.example),
        )
        
        for f in examples :
            self.failUnlessRaises(ValueError, 
                stockholm_io.read, f )

             
if __name__ == '__main__':
    unittest.main()
