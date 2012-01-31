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


from corebio.ssearch_io import *
from corebio.ssearch_io import fasta, blastxml
from corebio.utils import *
from test_corebio import testdata_stream
import unittest
from StringIO import StringIO

class test_search_io(unittest.TestCase) :    
   # def test_report(self) :
   #     rep = Report()
        
    def test_Annotation(self) :
        meta = Annotation()        

    def test_Hit(self) :
        h = Hit()        

    def test_Alignment(self) :
        a = Alignment()        
        s = str(a)   

class test_read(unittest.TestCase) :
    def test_scan(self) :
        files = [
            testdata_stream("ssearch/xbt001.xml"),        
            testdata_stream("ssearch/xbt002.xml"), 
            testdata_stream("ssearch/xbt003.xml"), 
            testdata_stream("ssearch/xbt004.xml"), 
            testdata_stream("ssearch/xbt005.xml"), 
            testdata_stream("ssearch/ssearch_out.txt"),        
            testdata_stream("ssearch/ssearch_out_compact.txt"),
            testdata_stream("ssearch/fasta_out_compact.txt"),
            testdata_stream("ssearch/ssearch_out_compact2_dbvdb.txt"),
            testdata_stream("ssearch/ssearch_out_compact_dbvdb.txt"),    
        ]
        results = [read(f) for f in files]
        self.assertEquals( len(files), len(results) )
      


class test_blastxml_read(unittest.TestCase) : 
    def examples(self) :
        files = [
            testdata_stream("ssearch/xbt001.xml"),        
            testdata_stream("ssearch/xbt002.xml"), 
            testdata_stream("ssearch/xbt003.xml"), 
            testdata_stream("ssearch/xbt004.xml"), 
            testdata_stream("ssearch/xbt005.xml"), 
            testdata_stream("ssearch/megablast.xml"), 
        ]
        return files        
        
    def test_read_empty(self) :
        f = StringIO()
        self.assertRaises(ValueError, blastxml.read, f)
                  
    def test_scan(self) :
        files = self.examples()
        reports = [blastxml.read(f) for f in files]
        self.assertEquals( len(files), len(reports) )
        
        # First report
        r = reports[0]
        self.assertEquals( r.algorithm, "blastp" )           
        self.assertEquals( r.algorithm_version, "2.2.12" )  
        self.assertEquals( len(r.results[0].hits), 212 )
        hit = r.results[0].hits[6]
        self.assertEquals( hit.target.name, "gi|49609685|emb|CAG73118.1|" )
        self.assertEquals( hit.target.length, 86 )
        self.assertAlmostEquals( hit.alignments[0].raw_score, 219 )
        self.assertAlmostEquals( hit.alignments[0].bit_score, 88.9669 )
        self.assertAlmostEquals( hit.alignments[0].significance, 4.10205e-17 )
        
        self.assertEquals("MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEP--KQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV", hit.alignments[0].query_seq)

        self.assertEquals( "MGGISIWQLLI+AVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDD+ KQDKTSQDADFTAK+IADKQ      +AK EDAK  DKEQV", r.results[0].hits[5].alignments[0].mid_seq)

        
        # Test megablast parseing with two queroes.
        r = reports[5] 
        self.assertEquals( r.algorithm, "blastn" )           
        self.assertEquals( r.algorithm_version, "2.2.14" ) 
        self.assertEquals( len(r.results), 2)
        

class test_fasta_read(unittest.TestCase) :  
    def examples(self) :
        files = [
            testdata_stream("ssearch/ssearch_out.txt"),        
            testdata_stream("ssearch/ssearch_out_compact.txt"),
            testdata_stream("ssearch/fasta_out_compact.txt"),
            testdata_stream("ssearch/ssearch_out_compact2_dbvdb.txt"),
            testdata_stream("ssearch/ssearch_out_compact_dbvdb.txt"),    
        ]
        return files

    def test_read_empty(self) :
        f = StringIO()
        self.assertRaises(ValueError, fasta.read, f)
                 
    def test_read(self) :
        files = self.examples()

        reports = [fasta.read(f) for f in files]
        self.assertEquals( len(files), len(reports) )
 
        # First report
        r = reports[0]
        #print repr(report)
        
        
        self.assertEquals( r.algorithm, "SSEARCH" )           
        self.assertEquals( r.algorithm_version, "3.4t24" )  
        self.assertEquals( len(r.results[0].hits), 16 )
        # d8rxna_ 7.34.3.1.1 Rubredoxin [Desulfovibrio vulg  (  52)   34  15.7     6.9
        hit = r.results[0].hits[10]
        self.assertEquals( hit.target.name, "d8rxna_" )
        self.assertEquals( hit.target.length, 52 )
        self.assertAlmostEquals( hit.raw_score, 34 )
        self.assertAlmostEquals( hit.bit_score, 15.7 )
        self.assertAlmostEquals( hit.significance, 6.9 )
        #print repr(hit)
        self.assertEquals( hit.alignments[0].target_seq, 'EGFLHLEDKPHPLQCQFFVESVIPAGSYQVPYRINVNNG-RPELAFDFKAMKRA..............')

        # Compact report
        r = reports[1]
        self.assertEquals( r.algorithm, "SSEARCH" )           
        self.assertEquals( r.algorithm_version, "3.4t24" )  
        self.assertEquals( len(r.results[0].hits), 16 )
        # d8rxna_ 7.34.3.1.1 Rubredoxin [Desulfovibrio vulg  (  52)   34  15.7     6.9
        hit = r.results[0].hits[10]
        self.assertEquals( hit.target.name, "d8rxna_" )
        self.assertEquals( hit.target.length, 52 )
        self.assertAlmostEquals( hit.raw_score, 34 )
        self.assertAlmostEquals( hit.bit_score, 15.7 )
        self.assertAlmostEquals( hit.significance, 6.9 )



        # DBvDB report
        #result = reports[4][1]
        r = reports[4]
        self.assertEquals( r.algorithm, "SSEARCH" )           
        self.assertEquals( r.algorithm_version, "3.4t24" )  
        self.assertEquals( len(r.results), 23 )
        self.assertEquals( len(r.results[1].hits), 21 )
        # d2sn3__ 7.3.6.1.1 scorpion toxin [Centruroides sc  (  65)   58  22.0   0.058 ...
        hit = r.results[1].hits[1]
        self.assertEquals( hit.target.name, "d2sn3__" )
        self.assertEquals( hit.target.length, 65 )
        self.assertAlmostEquals( hit.raw_score, 58 )
        self.assertAlmostEquals( hit.bit_score, 22 )
        self.assertAlmostEquals( hit.significance, 0.058  )


   # def test_this(self) :
    #    files = self.examples()
     #   for f in files :
     #       print f
      #      reports = fasta.read(f)
       #     print repr(reports)

if __name__ == '__main__':
    unittest.main()

