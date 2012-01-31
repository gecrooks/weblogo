#!/usr/bin/env python
 
#  Copyright (c) 2006 John Gilman
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


from corebio.transform import *
from corebio.seq import *
import unittest

      
class test_mask_low_complexity(unittest.TestCase):
    def test_segging(self):
        before = "mgnrafkshhghflsaegeavkthhghhdhhthfhvenhggkvalkthcgkylsigdhkqvylshhlhgdhslfhlehhggkvsikghhhhyisadhhghvstkehhdhdttfeeiii".upper()
        after = "MGNRAFKSHHGHFLSAEGEAVxxxxxxxxxxxxxxxENHGGKVALKTHCGKYLSIGDHKQVYLSHHLHGDHSLFHLEHHGGKVSIKGHHHHYISADHHGHVSTKEHHDHDTTFEEIII".upper()
        
        
        bseq = Seq(before, protein_alphabet)
        aseq = Seq(after, protein_alphabet)
        xseq = Seq( 'X' * len(bseq), protein_alphabet)
        
        sseq = mask_low_complexity(bseq)
        self.assertEquals(aseq, sseq)

        # Nothing should be segged
        sseq = mask_low_complexity(bseq, 12, 0,0)
        self.assertEquals(bseq, sseq)
     
     
        # Everthing should be segged
        sseq = mask_low_complexity(bseq, 12, 4.3,4.3)
        self.assertEquals(sseq, xseq)



    def test_seg_invalid(self):
        seq = Seq("KTHCGKYLSIGDHKQVYLSHH", protein_alphabet)
        self.failUnlessRaises(ValueError, mask_low_complexity,seq, 12, -1, 0 )
        self.failUnlessRaises(ValueError, mask_low_complexity, seq, 12, 1, 10 )
        self.failUnlessRaises(ValueError, mask_low_complexity, seq, 6, 12, 13 )
        self.failUnlessRaises(ValueError, mask_low_complexity, seq, 6, 2.0, 1.9 )


class test_transform(unittest.TestCase):
    def test_transform(self):
        trans = Transform( Seq("ACGTURYSWKMBDHVN",nucleic_alphabet),                    
                           Seq("ACGTTNNNNNNNNNNN", dna_alphabet) )
        s0 = Seq("AAAAAR",nucleic_alphabet)
        s1 = trans(s0)              # Callable ob
        self.assertEquals(s1.alphabet, dna_alphabet)
        self.assertEquals(s1 ,Seq("AAAAAN",  dna_alphabet))

      
   # def test_translations(self) :
   #     s  = Seq("ACGTURYSWKMBDHVNACGTURYSWKMBDHVN", nucleic_alphabet )
   #     s2 = dna_ext_to_std(s)
   #     s3 = Seq("ACGTTNNNNNNNNNNNACGTTNNNNNNNNNNN", dna_alphabet )
   #     self.assertEquals(s2, s3)
    
    def test_reduced_protein_alphabets(self):
        seq = Seq("ENHGGKVALKTHCGKYLSIGDHKQVYLSHHLHGDHSLFHLEHHGGKVSIKGHHHHYISADHHGHVSTKEHHDHDTTFEEIII", reduced_protein_alphabet)
        
        for t in reduced_protein_alphabets.values():
            s = t(seq)


    


class test_geneticcode(unittest.TestCase) :    
    def test_repr(self) :
        for t in GeneticCode.std_list() :
            r = repr(t)
            gc = eval(r)
            self.assertEquals( r, repr(gc) )
            self.assertEquals( str(gc), str(t) )
            #print r
            #print t
            #print gc
            
    def test_translate_std(self) :
        dna = 'GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA'
        t = GeneticCode.std()     
        s = t.translate(dna)             
        self.assertEquals( str(s), 'AIVMGR*KGAR')
        

        for t in GeneticCode.std_list() :
            p = t.translate(dna)

    def test_translate(self) :
        # Ref: http://lists.open-bio.org/pipermail/biopython/2006-March/002960.html
        
        cft = (
            ( "Vertebrate Mitochondrial",
                    'GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA', 'AIVMGRWKGAR'),
            ( 11,   'CAAGGCGTCGAAYAGCTTCAGGAACAGGAC',   'QGVE?LQEQD'),
            ( 1 ,   'GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA', 'AIVMGR*KGAR'),
            )
            
        for code, dna, protein in cft :                                                                    
            c = GeneticCode.by_name(code)
            trans = c.translate(dna)
            self.assertEquals( str(trans), protein)    

        
    def test_back_translate(self) :
        prot = 'ACDEFGHIKLMNPQRSTVWY*'
        t = GeneticCode.std()
        t.table['CGA']
        s = t.back_translate(prot)
        self.assertEquals( prot, str( t.translate(s) ) )



if __name__ == '__main__':
    unittest.main()

