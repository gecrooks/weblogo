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

from corebio import seq_io
from corebio.seq import *
import unittest
from StringIO import StringIO
from test_corebio import *


class test_alphabet(unittest.TestCase) :

    def test_create_alphabet(self) :
        # Alphabet contains repeated character
        self.failUnlessRaises(ValueError, Alphabet, "alphabet")

        # Alphbaet contains null character
        self.failUnlessRaises(ValueError, Alphabet, "alph\x00")
        
        a1 = Alphabet("alphbet")

    def test_alphabet_alphabetic(self):
        a = Alphabet("alphbet")
        self.failUnless(a.alphabetic("alphbet"))
        self.failUnless(not a.alphabetic("alphbetX"))

    def test_alphabet_ord(self):
        a = generic_alphabet
        for i,c in enumerate(a) :
            self.failUnlessEqual( a.ord(c), i)

        a = Alphabet("alph")
        self.failUnlessEqual(2, a.ord("p"))

    def test_alphabet_chr(self):
        a = generic_alphabet
        for i,c in enumerate(a) :
            self.failUnlessEqual( ord(a.chr(i)), i+32 )            
    
        a = Alphabet("alph")
        self.failUnlessEqual("h", a.chr(3))


    def test_alphabet_ords(self):
        a = Alphabet("alph")
        self.failUnlessEqual(0, a.ords("alphalph")[4])

        a = generic_alphabet
        o = a.ords(a)
        for i,c in enumerate(o) :
            self.failUnlessEqual( c, i)

    def test_alphabet_chrs(self):
        a = Alphabet("alph")
        self.failUnlessEqual(Seq("ppla",a), a.chrs((2,2,1,0)) )
                
    def test_none(self) :
        a1 = Alphabet(None)
        self.failUnlessEqual(a1, generic_alphabet)
    
    def test_create_from_alphabet(self) :
        """ If we pass an alphabet to the constuctor, it's passed 
        right back """
        a1 = Alphabet("kjdahf")
        a2 = Alphabet(a1)
        self.failUnless(a1 == a2)
        
    def test_repr(self):
        a = Alphabet("kjdahf")
        r = repr(a)
        s = str(a)
   
    def test_str(self):
        self.failUnlessEqual( str(Alphabet("kjdahf")), "kjdahf")
        
    def test_letters(self) :        
        self.failUnlessEqual( Alphabet("kjdahf").letters(), "kjdahf")
    
    def test_normalize(self) :
        a = Alphabet("ABCDE")
        s = 'aBbc'
        n = a.normalize(s)
        self.failUnlessEqual( str(n), 'ABBC')

    def test_alt(self) :
        alt = zip('12345', 'ABCED')
        alpha = Alphabet('ABCDE',  alt)
        assert alpha.ord('A')== 0
        for a, c in alt :
            assert alpha.ord(a) == alpha.ord(c)
    
    def test_which_alphabet(self):
        a = Alphabet.which(Seq("ARNDCQEGHILKMFPSTWYVX"))
        assert  a == unambiguous_protein_alphabet
        
        tests = ( 
             (seq_io.read(testdata_stream('cap.fa')), unambiguous_dna_alphabet),
             (seq_io.read(testdata_stream('cox2.msf')), unambiguous_protein_alphabet),
             (seq_io.read(testdata_stream('Rv3829c.fasta')), unambiguous_protein_alphabet),
             )
        for t in tests :
            self.failUnlessEqual(Alphabet.which(t[0]), t[1])


      

class test_seq(unittest.TestCase):
    
    def test_create_seq(self) :
        self.failUnless(Seq("alphabet", "alphbet"))
        self.assertRaises(ValueError, Seq, "not alphabetic", "alphabet" )
        
        a = "Any printable Ascii character `1234567890-=~!@#$%^&*()_+{}|[]\:;'<>?,./QWERTYUIOPASDFGHJKLZXCVBNMqwertyuiopasdfghjklzxcvbnm "
        
        for l in a :
            self.failUnless( l  in generic_alphabet )        
        self.failUnless(Seq(a, generic_alphabet))
        self.failUnlessRaises(ValueError, Seq,
                              "Not zero. \x00", generic_alphabet)

   
    def test_std_alphabets(self) :
        s = Seq("dskjjfskdjbfKJJSJKSKJDjk123u768erbm", generic_alphabet)
        s = Seq("ARNDCQEGHILKMFPSTWYVX", protein_alphabet)
        self.assertRaises(ValueError, Seq, "1234", protein_alphabet)
        s = Seq("AGTCAGCTACGACGCGC", dna_alphabet)
        self.assertEquals( str(s[1]), 'G')
  
    
    def test_ords(self) :
        a = Alphabet("ABC")
        s = Seq("ABCCBA", a)
        self.assertEquals( list(s.ords()), [0,1,2,2,1,0] )
        
    def test_tally(self) :
        s = Seq("AGTCAGCTACGACGCGC", unambiguous_dna_alphabet)
        c = s.tally()
        self.assertEquals(len(unambiguous_dna_alphabet), len(c))
        self.assertEquals(list(c), [4,6,5,2])
        
    def test_tally_nonalphabetic(self) :
        s = Seq("AGTCAGCTACGACGCGC", dna_alphabet)
        c = s.tally( Alphabet("AC"))
        self.assertEquals(2, len(c))
        self.assertEquals(list(c), [4,6])


    def test_words(self) :
        s = Seq("AGTCAGCTACGACGcgcx", dna_alphabet)
        w = list(s.words(2,unambiguous_dna_alphabet ))
        self.assertEquals(len(w),len(s)-2)
        self.assertEquals(w, ['AG', 'GT', 'TC', 'CA', 'AG', 'GC', 'CT', 'TA', 'AC', 'CG', 'GA', 'AC', 'CG', 'GC', 'CG', 'GC'])
        
        self.assertEquals(list(s.words(len(s), unambiguous_dna_alphabet)), [])
        self.assertEquals(list(s.words(len(s)-1,unambiguous_dna_alphabet)), ["AGTCAGCTACGACGCGC",])

    def test_words2(self) :
        s = Seq("AGTCAGCTACGACGCGC", unambiguous_dna_alphabet)
        wc = s.word_count(2)
        count = zip(*wc)[1]
        self.assertEquals(count, (2,2,1,3,1,1,3,1,1,1) )
                
    def test_getslice(self):
        s = Seq("AGTCAGCTACGACGCGC", dna_alphabet)
        slice = s[2:4]
        self.assertEquals( s.alphabet, slice.alphabet)

    def test_create_annotated(self):
        s= "ACGTURYSWKMBDHVNACGTURYSWKMBDHVNAAAAA"
        a  = Seq(s, nucleic_alphabet,
            name = "ID", description = "DESCRIPTION" )
        self.assertEquals(a.name, "ID")
        self.assertEquals(a.description, "DESCRIPTION")
        self.assertEquals(s, str(a))


    def test_add(self):
        s1 = Seq("AAAA", dna_alphabet)
        s2 = Seq("TTTT", dna_alphabet)
        
        s3 = s1 +s2
        self.assertEquals( s3.alphabet, dna_alphabet)
        self.assertEquals( s3, Seq("AAAATTTT", dna_alphabet))
    
        s4 = "AA"
        s5 = s4 + s1
        s6 = s1 + s4
        self.assertEquals( s5.alphabet, s6.alphabet)
        self.assertEquals( s5, s6)
    
    
    def test_join(self):
        s1 = Seq("AAAA", dna_alphabet)
        s2 = Seq("TTTT", dna_alphabet)
        s3 = "GGGG"
        s0 = Seq("", dna_alphabet)
        
        j = s0.join ([s1,s2,s3])
        self.assertEquals( j, Seq("AAAATTTTGGGG", dna_alphabet))
    

    def test_repr(self):
        s1 = Seq("AAAA", dna_alphabet)
        r = repr(s1)
        
    def test_str(self) :
        s1 = Seq("AGCTA", dna_alphabet)
        self.assertEquals(str(s1), "AGCTA")
        # Uncased alpahebt
        self.assertEquals(str(Seq("AgcTA", dna_alphabet)), "AgcTA")
    
    def test_tostring(self) :        
        self.assertEquals(Seq("AgcTAAAA", dna_alphabet).tostring(), "AgcTAAAA")
      

    def test_reverse(self) :
        s = Seq("ACGT", dna_alphabet)
        self.assertEquals(s, s.reverse().reverse())
        self.assertEquals(s.reverse(), Seq("TGCA", dna_alphabet))
   
    def test_translate(self) :
        s = dna('GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA')
        p = s.translate()             
        self.assertEquals( str(p), 'AIVMGR*KGAR')
   
    def test_reverse_complement(self) :
        s = dna('GGGGaaaaaaaatttatatat')
        rc = s.reverse_complement()
        self.assertEquals(rc, dna('atatataaattttttttCCCC'))
   
        self.assertEquals( dna('ACGTRYSWKMBDHVN').reverse_complement(),
                    dna('NBDHVKMWSRYACGT') )
   
   
    def test_ungap(self):
        s = Seq("T-T", dna_alphabet).ungap()
        self.assertEquals(str(s), 'TT')
        s = Seq("T-~---T...~~~--", dna_alphabet).ungap()
        self.assertEquals(str(s), 'TT')
   
    def test_mask(self):
        s = dna('AAaaaaAAA').mask()
        self.assertEquals(str(s), 'AAXXXXAAA')
        
    
    def test_shortcuts(self) :
        sprotein = protein('gGGGGG-PPPPP') 
        d1 = dna('ACGTRYSWKMBDHVN')
        d2 = dna('t')
        d3 = dna('actgta')
        d4 = dna('acgtrysw-kmb-dhvn')
        
        srna = rna('ACUUUUU')

    def test_casechange(self) :
        s1 = dna('ACGTRYSWKMBDHVN')
        s2 = s1.lower().upper()
        self.assertEquals( s1, s2)
        
        
   
class test_seqlist(unittest.TestCase):
    
    
    def test_create(self):
                 # 1234567890123456789012345678
        s0 = Seq("ACGTURYBDHVNACGTURYSWKMBDHVN", nucleic_alphabet )    
        s1 = Seq("ACGTURYSWKMBDHVNACGTURYSWKMBDHVN", nucleic_alphabet )    
        s2 = Seq("ACGTURSWKMBDHVNACGTURKMBDHVN", nucleic_alphabet )    
        seqs = SeqList( [ s0,s1,s2])
        
        self.assertEquals( len(seqs), 3)

    def test_create_annotated(self):
        s0 = Seq("ACGTURYBDHVNACGTURYSWKMBDHVN", nucleic_alphabet )    
        s1 = Seq("ACGTURYSWKMBDHVNACGTURYSWKMBDHVN", nucleic_alphabet )    
        s2 = Seq("ACGTURSWKMBDHVNACGTURKMBDHVN", nucleic_alphabet )    
   
        seqs = SeqList([ s0,s1,s2], alphabet = nucleic_alphabet,
            name = "alsdf", description='a')
        self.assertEquals( seqs.name, 'alsdf')
        self.assertEquals( seqs.description, 'a')
        self.assertEquals( seqs.alphabet, nucleic_alphabet)


    def test_ords(self) :
        s0 = Seq("ACGTURYBDHVNACGTURYSWKMBDHVN", nucleic_alphabet )    
        s1 = Seq("ACGTURYSDHVNACGTURYSWKMBDHVN", nucleic_alphabet )    
        s2 = Seq("ACGTURSWKMBDHVNACGTURKMBDHVN", nucleic_alphabet )    
        seqs = SeqList( [ s0,s1,s2],  nucleic_alphabet)
        a = seqs.ords()
        #self.assertEquals( a.shape, (3, 28) )

        # Fails if seqs are of different lengths
        # FIXME?
        #s3 = Seq("ACGTUR", nucleic_alphabet )   
        #seqs2 = SeqList( [ s0,s1,s3,s2],  nucleic_alphabet)
        #self.failUnlessRaises(ValueError, seqs2.ords )
        
        # Use a different alphabet
        a2 = seqs.ords(nucleic_alphabet)

        # No alphabet
        seqs3 = SeqList( [ s0,s1,s2])
        a3 = seqs3.ords(alphabet = Alphabet("ABC"))
        
        # Fail if no alphabet
        self.failUnlessRaises(ValueError, seqs3.ords )
    
    def test_profile(self) :
        a = Alphabet("ABCD")
        
        s0 = Seq("ABCDD", a )    
        s1 = Seq("AAAAD", a )
        s2 = Seq("AAABD", a )
        s3 = Seq("AAACD", a )
            
        seqs = SeqList( [ s0,s1,s2,s3],  a)

        tally = seqs.profile()

        self.assertEquals( list(tally[0]), [4,0,0,0] )
        self.assertEquals( list(tally[1]), [3,1,0,0] )
        self.assertEquals( list(tally[2]), [3,0,1,0] )
        self.assertEquals( list(tally[3]), [1,1,1,1] )                
        self.assertEquals( list(tally[4]), [0,0,0,4] )

        self.assertEquals( tally[4,'D'], 4 )

    def test_tally(self):
        # 1234567890123456789012345678
        s0 = Seq("ACTTT", nucleic_alphabet )    
        s1 = Seq("ACCCC", nucleic_alphabet )    
        s2 = Seq("GGGG", nucleic_alphabet )    
        seqs = SeqList( [ s0,s1,s2],nucleic_alphabet)

        counts = seqs.tally()
        assert counts == [2, 5, 4, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        
       
             
    def test_create_empty(self):
        s0 = Seq("ACGTURYBDHVNACGTURYSWKMBDHVN", nucleic_alphabet )    
        s1 = Seq("ACGTURYSWKMBDHVNACGTURYSWKMBDHVN", nucleic_alphabet )    
        s2 = Seq("ACGTURSWKMBDHVNACGTURKMBDHVN", nucleic_alphabet )    
        
        seqs = SeqList()
        seqs.append(s0)
        seqs.extend( (s1,s2))
        
        self.assertEquals( len(seqs), 3)
        self.assertEquals( type(seqs), SeqList)
    
    
    def test_repr(self):
        s0 = Seq("ACGTURYBDHVNACGTURYSWKMBDHVN", nucleic_alphabet )    
        s1 = Seq("ACGTURYSWKMBDHVNACGTURYSWKMBDHVN", nucleic_alphabet )    
        s2 = Seq("ACGTURSWKMBDHVNACGTURKMBDHVN", nucleic_alphabet )    
        seqs = SeqList( [ s0,s1,s2])
    
        r = repr(seqs)
        
    def test_str(self) :
        s0 = Seq("ACGTURYBDHVNACGTURYSWKMBDHVN", nucleic_alphabet )    
        s1 = Seq("ACGTURYSWKMBDHVNACGTURYSWKMBDHVN", nucleic_alphabet )    
        s2 = Seq("ACGTURSWKMBDHVNACGTURKMBDHVN", nucleic_alphabet )    
        seqs = SeqList( [ s0,s1,s2])
        s = str(seqs)
        
      
if __name__ == '__main__':
    unittest.main()
