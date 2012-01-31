#!/usr/bin/env python

"""Unit test for Raf"""

import unittest

from corebio.db.astral  import *
from test_corebio import *

class Astraltest(unittest.TestCase):
   
   
    def testParseDomain(self) :
        s=">d1tpt_1 a.46.2.1 (1-70) Thymidine phosphorylase {Escherichia coli}"
        dom = parse_domain(s)

        assert dom.sid == 'd1tpt_1'
        assert dom.sccs == 'a.46.2.1'
        assert dom.residues.pdbid == '1tpt'
        assert dom.description == 'Thymidine phosphorylase {Escherichia coli}'

        s2="d1tpt_1 a.46.2.1 (1tpt 1-70) Thymidine phosphorylase {E. coli}"
        assert s2 == str(parse_domain(s2)), str(parse_domain(s2))



        #Genetic domains (See Astral release notes)
        s3="g1cph.1 g.1.1.1 (1cph B:,A:) Insulin {Cow (Bos taurus)}"
        assert s3 == str(parse_domain(s3)), str(parse_domain(s3))

        s4="e1cph.1a g.1.1.1 (1cph A:) Insulin {Cow (Bos taurus)}"
        assert s4 == str(parse_domain(s4))

        #Raw Astral header
        s5=">e1cph.1a g.1.1.1 (A:) Insulin {Cow (Bos taurus)}"
        assert s4 ==  str(parse_domain(s5))

        try:
            dom = parse_domain("Totally wrong")
            assert 0, "Should never get here"
        except ValueError, e :
            pass


class RafTest(unittest.TestCase):
    rafLine = "101m_ 0.01 38 010301 111011    0  153    0 mm   1 vv   2 ll   3 ss   4 ee   5 gg   6 ee   7 ww   8 qq   9 ll  10 vv  11 ll  12 hh  13 vv  14 ww  15 aa  16 kk  17 vv  18 ee  19 aa  20 dd  21 vv  22 aa  23 gg  24 hh  25 gg  26 qq  27 dd  28 ii  29 ll  30 ii  31 rr  32 ll  33 ff  34 kk  35 ss  36 hh  37 pp  38 ee  39 tt  40 ll  41 ee  42 kk  43 ff  44 dd  45 rr  46 vv  47 kk  48 hh  49 ll  50 kk  51 tt  52 ee  53 aa  54 ee  55 mm  56 kk  57 aa  58 ss  59 ee  60 dd  61 ll  62 kk  63 kk  64 hh  65 gg  66 vv  67 tt  68 vv  69 ll  70 tt  71 aa  72 ll  73 gg  74 aa  75 ii  76 ll  77 kk  78 kk  79 kk  80 gg  81 hh  82 hh  83 ee  84 aa  85 ee  86 ll  87 kk  88 pp  89 ll  90 aa  91 qq  92 ss  93 hh  94 aa  95 tt  96 kk  97 hh  98 kk  99 ii 100 pp 101 ii 102 kk 103 yy 104 ll 105 ee 106 ff 107 ii 108 ss 109 ee 110 aa 111 ii 112 ii 113 hh 114 vv 115 ll 116 hh 117 ss 118 rr 119 hh 120 pp 121 gg 122 nn 123 ff 124 gg 125 aa 126 dd 127 aa 128 qq 129 gg 130 aa 131 mm 132 nn 133 kk 134 aa 135 ll 136 ee 137 ll 138 ff 139 rr 140 kk 141 dd 142 ii 143 aa 144 aa 145 kk 146 yy 147 kk 148 ee 149 ll 150 gg 151 yy 152 qq 153 gg"

    rafLine2 = "101mA 0.01 38 010301 111011    0  153    0 mm   1 vv   2 ll   3 ss   4 ee   5 gg   6Aee   7Aww   8Aqq"

    rafLine3 = "101mB 0.01 38 010301 111011    0  153   90 mm  91 vv  92 ll  939ss  94 ee  95 gg"


    def test_Parse(self):
        """Can we parse a RAF record?"""
        r = RafSeqMap(self.rafLine)

        assert r.pdbid == "101m"
        assert r.pdb_datestamp =="010301"  
        assert r.flags =="111011"          
      
        i = r.index("143")
        res = r.res[i]
        assert res.chainid =="_"        
        assert res.resid =="143"
        assert res.seqres =="A"
        assert res.atom =="A"

        r = RafSeqMap(self.rafLine2)   
        res = r.res[r.index("6A", chainid="A")]
        assert res.resid =="6A"
        assert res.atom=="E"

    def test_SeqMapAdd(self) :
        r2 = RafSeqMap(self.rafLine2)
        r3 = RafSeqMap(self.rafLine3)

        l = len(r2.res) + len(r3.res)
        r2 += r3
        assert len(r2.res) == l

        r2.extend(r2)
        assert len(r2.res) == l*2

        r4 = r2 + r2
        assert len(r4.res) == l*4

        r4.append(Res())
        assert len(r4.res) == (l*4)+1
        

    def test_SeqMapSlice(self) :
        r = RafSeqMap(self.rafLine)
        r = r[ r.index("124"): r.index("135")+1]
        assert len(r.res) ==12

    def test_SeqMapIndex(self) :
        filename = testdata_stream("scop/raftest.txt").name
        f = open(filename)
        index = Raf(f)
        r = index.get_seqmap("103m")
        assert r.pdbid == "103m", r.pdbid
        assert len(r.res) ==154, len(r.res)
        assert r.pdb_datestamp =="010301"  
        assert r.flags =="111011"

        r = index.get_seqmap("103m 1-10")
        assert r.pdbid == "103m", r.pdbid
        assert len(r.res) ==10, len(r.res)
        assert r.pdb_datestamp =="010301"  
        assert r.flags =="111011"        

        r = index.get_seqmap("104l A:")
        assert r.pdbid == "104l", r.pdbid

        r = index.get_seqmap("104l A:112-113")
        assert r.pdbid == "104l", r.pdbid        
        assert len(r.res)== 2

        r = index.get_seqmap("104l A:112-113,B:146-148")
        assert r.pdbid == "104l", r.pdbid        
        assert len(r.res)== 5        


        assert "103m_" in index
        r = index["103m_"]
        assert r.pdbid == "103m", r.pdbid


    def test_Parse_error(self):
        self.assertRaises(ValueError, RafSeqMap, "tooshort")

    def test_records(self) :
        f = testdata_stream("scop/raftest.txt")
        i=0
        for rsm in RafSeqMap.records(f) : i +=1
        assert i == 16

if __name__ == '__main__':
    unittest.main()










