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
from corebio.seq import *
from corebio.secstruc.stride import *
from test_corebio import *



import unittest

class test_stride_io(unittest.TestCase) :

    def test_1(self) :
        g = StrideRecord(testdata_stream("stride/stride_test_1.txt"))
        self.assertEquals(g.pdbid, "1A65")
        self.assertEquals(g.residues[0].chainid, "A")
        self.assertEquals(g.residues[0].secstruc, "C")
        self.assertEquals(g.residues[0].solvent_acc_area, float(95.4))
        self.assertEquals(g.residues[0].phi, float(360.00))
        self.assertEquals(g.residues[0].psi, float(157.13))  
        self.assertEquals(g.residues[0].resid, "1")
        self.assertEquals(g.primary(), Seq("QIVNSVDTMT", protein_alphabet))
        self.assertEquals(g.secondary(), Seq("CEETTEEEEE",stride_alphabet))
        self.assertEquals(g.total_area(), float(483.6))
        assert(g.get_residue('A','1') is g.residues[0] ) 
        
    def test_2(self) :
        g = StrideRecord(testdata_stream("stride/stride_test_2.txt"))
        self.assertEquals(g.pdbid, "1A59")
        self.assertEquals(g.residues[1].chainid, " ")
        self.assertEquals(g.residues[1].secstruc, "C")
        self.assertEquals(g.residues[1].solvent_acc_area, float(85.3))
        self.assertEquals(g.residues[1].phi, float(-78.20))
        self.assertEquals(g.residues[1].psi, float(165.31))  
        self.assertEquals(g.residues[1].resid, "3")
        self.assertEquals(g.primary(), Seq("EPTIH", protein_alphabet))
        self.assertEquals(g.secondary(), Seq("CCCCC",stride_alphabet))
        self.assertEquals(g.total_area(), float(610.9))
        assert(g.get_residue(' ','3') is g.residues[1] ) 
    
    def test_3(self):
        g = StrideRecord(testdata_stream("stride/stride_test_3.txt"))
        self.assertEquals(g.pdbid, "1A59")
        self.assertEquals(g.residues[0].chainid, " ")
        self.assertEquals(g.residues[0].secstruc, "T")
        self.assertEquals(g.residues[3].solvent_acc_area, float(60.4))
        self.assertEquals(g.residues[4].phi, float(-103.84))
        self.assertEquals(g.residues[2].psi, float(-27.45))  
        self.assertEquals(g.residues[0].resid, "12")
        self.assertEquals(g.primary(), Seq("VTADV", protein_alphabet))
        self.assertEquals(g.secondary(), Seq("TCCCC",stride_alphabet))
        self.assertEquals(g.total_area(), float(404))
        assert(g.get_residue(' ','12') is g.residues[0] ) 
   


             
             
if __name__ == '__main__':
    
    print "Running additional tests of RunStride. These require that the STRIDE program is installed locally"
    try :
        stride = RunStride()
        fn = testdata_filename('1CGP.pdb')
        #print fn
        data = stride.process_pdb(fn)
        #print data
        record = stride.record(fn)
        #print record
        # Now run standard unittests
    except e:
        print e
    
    unittest.main()
    
    