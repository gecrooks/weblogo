#!/usr/bin/env python


# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
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
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
#  IN THE SOFTWARE.

"""Unit test for Scop"""

import os.path
import unittest

from corebio.db.scop import *
from corebio._py3k import StringIO

from test_corebio import *


class test_scop(unittest.TestCase):

    def testParse(self):
        f = testdata_stream('scop/dir.cla.scop.txt_test')
        try:
            cla = f.read()
            f.close()
            f = testdata_stream('scop/dir.des.scop.txt_test')
            des = f.read()
            f.close()
            f = testdata_stream('scop/dir.hie.scop.txt_test')
            hie = f.read()
        finally:
            f.close()

        scop = Scop.parse_files(StringIO(cla), StringIO(des), StringIO(hie))

        cla_out = StringIO()
        scop.write_cla(cla_out)
        assert cla_out.getvalue() == cla, cla_out.getvalue()
        
        des_out = StringIO()
        scop.write_des(des_out)
        assert des_out.getvalue() == des, des_out.getvalue()

        hie_out = StringIO()
        scop.write_hie(hie_out)
        assert hie_out.getvalue() == hie, hie_out.getvalue()

        domain = scop.domains_by_sid["d1hbia_"]
        self.assertEqual(domain.sunid, 14996)

        domains = scop.domains
        self.assertEqual(len(domains), 14)
        self.assertEqual(domains[4].sunid, 14988)

        self.assertFalse(-111 in scop.nodes_by_sunid)
        self.assertFalse("no such domain" in scop.domains_by_sid )


    def testSccsOrder(self) :
        self.assertEqual(cmp_sccs("a.1.1.1", "a.1.1.1"), 0)
        self.assertEqual(cmp_sccs("a.1.1.2", "a.1.1.1"), 1)
        self.assertEqual(cmp_sccs("a.1.1.2", "a.1.1.11"), -1)
        self.assertEqual(cmp_sccs("a.1.2.2", "a.1.1.11"), 1)
        self.assertEqual(cmp_sccs("a.1.2.2", "a.5.1.11"), -1)
        self.assertEqual(cmp_sccs("b.1.2.2", "a.5.1.11"), 1)
        self.assertEqual(cmp_sccs("b.1.2.2", "b.1.2"), 1)

    def test_sccs_relations(self):
        self.assertEqual(sccs_relation("a.1.1.1", "a.1.1.1"), 1)
        self.assertEqual(sccs_relation("a.1.1.2", "a.1.1.1"), 1)
        self.assertEqual(sccs_relation("a.1.1.2", "a.1.1.11"), 1)
        self.assertEqual(sccs_relation("a.1.2.2", "a.1.1.11"), 0)
        self.assertEqual(sccs_relation("a.1.2.2", "a.5.1.11"), -1)
        self.assertEqual(sccs_relation("b.1.2.2", "a.5.1.11"), -1)
        self.assertEqual(sccs_relation("b.1.2.2", "b.1.2"), 1)



    def testConstructFromDirectory(self):
        dir_path = os.path.join( os.path.split(__file__)[0], 'data/scop')
        scop = Scop.parse(dir_path= dir_path, version="test")
        self.assertTrue(isinstance(scop, Scop))
        domain = scop.domains_by_sid["d1hbia_"]
        self.assertEqual(domain.sunid, 14996)

    def testGetAscendent(self):
        dir_path = os.path.join( os.path.split(__file__)[0], 'data/scop')
        scop = Scop.parse(dir_path=dir_path, version="test")
        domain = scop.domains_by_sid["d1hbia_"]
        # get the fold
        fold = domain.ascendent('cf')
        self.assertEqual(fold.sunid, 46457)
        #get the superfamily
        sf = domain.ascendent('superfamily')
        self.assertEqual(sf.sunid, 46458)
        # px has no px ascendent
        px = domain.ascendent('px')
        self.assertTrue(px is None)
        # an sf has no px ascendent
        px2 = sf.ascendent('px')
        self.assertTrue(px2 is None)


    def test_get_descendents(self):
        """Test getDescendents method"""
        dir_path = os.path.join( os.path.split(__file__)[0], 'data/scop')
        scop = Scop.parse(dir_path=dir_path, version="test")
        fold = scop.nodes_by_sunid[46457]
        # get px descendents
        domains = fold.descendents('px')
        self.assertEqual(len(domains), 14)
        for d in domains:
            self.assertEqual(d.type, 'px')
        sfs = fold.descendents('superfamily')
        self.assertEqual(len(sfs), 1)
        for d in sfs:
            self.assertEqual(d.type, 'sf')
        # cl has no cl descendent
        cl = fold.descendents('cl')
        self.assertEqual(cl, [])


class DesTests(unittest.TestCase):

    def setUp(self) :
        file = testdata_stream("scop/dir.des.scop.txt_test")
        self.filename = file.name
        file.close()

    def test_parse(self):
        with open(self.filename) as f:
            count = 0
            for rec in DesRecord.records(f):
                count +=1
            self.assertEqual(count, 20)

    def testStr(self):
        with open(self.filename) as f:
            for line in f :
                rec = DesRecord(line)
                #End of line is platform dependent. Strip it off
                self.assertEqual(str(rec).rstrip(), line.rstrip())

    def testError(self) :
        corruptRec = "49268\tsp\tb.1.2.1\t-\n"
        self.assertRaises(ValueError, DesRecord, corruptRec)

    def testRecord(self) :
        recLine = '49268\tsp\tb.1.2.1\t-\tHuman (Homo sapiens)    \n'
        recFields = (49268, 'sp', 'b.1.2.1', '', 'Human (Homo sapiens)')
        rec = DesRecord(recLine)
        self.assertEqual(rec.sunid, recFields[0])
        self.assertEqual(rec.nodetype, recFields[1])
        self.assertEqual(rec.sccs, recFields[2])
        self.assertEqual(rec.name, recFields[3])
        self.assertEqual(rec.description, recFields[4])


class test_scop_cla(unittest.TestCase):

    def setUp(self) :
        file = testdata_stream("scop/dir.cla.scop.txt_test")
        self.filename = file.name
        file.close()

    def testParse(self):
        """Can we parse a CLA file?"""
        with open(self.filename) as f:
            count = 0
            for rec in ClaRecord.records(f):
                count +=1
            self.assertEqual(count, 14)

    def testStr(self):
        with open(self.filename) as f:
            for line in f :
                rec = ClaRecord(line)
                #End of line is platform dependent. Strip it off
                self.assertEqual(str(rec).rstrip(), line.rstrip())

    def testError(self) :
        corruptRec = "49268\tsp\tb.1.2.1\t-\n"
        self.assertRaises(ValueError, ClaRecord, corruptRec)

    def testRecord(self) :
        recLine = 'd1dan.1\t1dan\tT:,U:91-106\tb.1.2.1\t21953\tcl=48724,cf=48725,sf=49265,fa=49266,dm=49267,sp=49268,px=21953'

        rec = ClaRecord(recLine)
        self.assertEqual(rec.sid, 'd1dan.1')
        self.assertEqual(rec.residues.pdbid, '1dan')
        self.assertEqual(rec.residues.fragments, (('T','',''),('U','91','106')))
        self.assertEqual(rec.sccs, 'b.1.2.1')
        self.assertEqual(rec.sunid, 21953)
        self.assertEqual(rec.hierarchy, [
            ['cl',48724], ['cf',48725], ['sf',49265], ['fa',49266],
            ['dm',49267], ['sp',49268], ['px',21953]])


class DomTests(unittest.TestCase):
    def setUp(self) :
        file = testdata_stream('scop/domtest.txt')
        self.filename = file.name
        file.close()
    
    def testParse(self):
        with open(self.filename) as f:
            count = 0
            for rec in DomRecord.records(f):
                count +=1
            self.assertEqual(count, 10)

    def testStr(self):
        with open(self.filename) as f:
            for line in f:
                if line:
                    rec = DomRecord(line)
                    self.assertEqual(str(rec).rstrip(), line.rstrip())

    def testError(self) :
        corruptDom = "49xxx268\tsp\tb.1.2.1\t-\n"
        self.assertRaises(ValueError, DomRecord, corruptDom)

    def testRecord(self) :
        recLine = 'd7hbib_\t7hbi\tb:\t1.001.001.001.001.001'
        rec = DomRecord(recLine)
        self.assertEqual(rec.sid, 'd7hbib_')
        self.assertEqual(rec.residues.pdbid,'7hbi')
        self.assertEqual(rec.residues.fragments, (('b','',''),))
        self.assertEqual(rec.hierarchy,'1.001.001.001.001.001')


class ResiduesTests(unittest.TestCase):
    res = (
        ( "-",           () ),
        ( "A:",          (("A", "", ""),) ),
        ( "1:",          (("1", "", ""),) ),
        ( "1-100",       (("", "1", "100"),)  ),
        ( "B:1-101",     (("B",   "1" ,"101"),) ),
        ( "1:1a-100a",   (("1", "1a", "100a"),) ),
        ( "a:-100a--1a", (("a", "-100a", "-1a"),) ),
        ( "-1-100",      (("", "-1", "100"),) ),
        ( "-1-100",      (("", "-1", "100"),) ),
        ( "A:12-19,A:23-25", (("A","12","19"),("A","23","25")) ),
        ( "12-19,1:23-25", (("","12","19"),("1","23","25")) ),
        ( "0-1,1:-1a-25a,T:", (("","0","1"),("1","-1a","25a"),("T","","")) ),
        )


    def testParse(self):
        for loc in self.res :
            r = Residues(loc[0])
            assert r.fragments == loc[1], str(r.locations)

    def testStr(self):
        for loc in self.res :
            r = Residues(loc[0])
            self.assertEqual(str(r), loc[0])

    def testAstralParse(self) :
        """Astral encloses residue subsets in brackets. Lets make sure we
        can parse those too.
        """
        for loc in self.res :
            r = Residues("("+loc[0]+")")
            assert r.fragments == loc[1], str(r.locations)

    def testPdbId(self):
        pdbid ="1ddf"
        for loc in self.res :
            r = Residues("\t 1ddf \t"+loc[0]+"\t\n\n\n")
            self.assertEqual(r.pdbid, pdbid)
            self.assertEqual(str(r), pdbid+" "+loc[0])

            r = Residues(pdbid+" "+loc[0])
            self.assertEqual(r.pdbid, pdbid)
            self.assertEqual(str(r), pdbid+" "+loc[0])

            r = Residues("104l A:112-113")
            self.assertEqual(r.pdbid, "104l")
            self.assertEqual(r.fragments, (('A', '112', '113'),))

    def testJustPdbId(self) :
        r = Residues("1sds")
        self.assertEqual(r.pdbid, "1sds")
        assert not r.fragments


    def testParseError(self) :
        self.assertRaises(ValueError, Residues, "09324923423hh./;,.389")


class HieTests(unittest.TestCase):

    def setUp(self) :
        file = testdata_stream("scop/dir.hie.scop.txt_test")
        self.filename = file.name
        file.close()

    def testParse(self):
        with open(self.filename) as f:
            count = 0
            for rec in HieRecord.records(f):
                count += 1
            self.assertEqual(count, 21, "Wrong number of records?!")

    def testStr(self):
        with open(self.filename) as f:
            for line in f :
                rec = HieRecord(line)
                #End of line is platform dependent. Strip it off
                self.assertEqual(str(rec).rstrip(), line.rstrip())

    def testError(self):
        corruptRec = "4926sdfhjhfgyjdfyg"
        self.assertRaises(ValueError, HieRecord, corruptRec)



if __name__ == '__main__':
    unittest.main()
