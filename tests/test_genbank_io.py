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


import unittest

from weblogo import seq_io
from weblogo.seq_io import genbank_io
from . import data_stream


def examples():
    return (
        data_stream('genbank/NT_019265.gb'),
        data_stream('genbank/cox2.gb'),
        data_stream('genbank/iro.gb'),
        data_stream('genbank/pri1.gb'),
        data_stream('genbank/dbsource_wrap.gb'),
        data_stream('genbank/noref.gb'),
        data_stream('genbank/protein_refseq.gb'),
        data_stream('genbank/cor6_6.gb'),
        data_stream('genbank/origin_line.gb'),

        # These files are too large to include in the distribution
        # data_stream('genbank/arab1.gb'),
        # data_stream('genbank/NC_005213.gbk'),
        # data_stream('genbank/NC_003888.gbk'),
    )


class test_genbank_io(unittest.TestCase):
    # Useful for debugging
    # def test_scan(self) :
    #    for f in examples():
    #        for t in genbank_io._scan(f):
    #            print t
    #        print
    #        print

    def test_parse(self):
        for f in examples():
            # print f.name
            genbank_io.read(f)
            f.close()
            # print seqs

    def test_read(self):
        f = data_stream("genbank/cox2.gb")
        seqs = genbank_io.read(f)

        self.assertEqual(len(seqs), 5)
        self.assertEqual(len(seqs[1]), 210)

        f.seek(0)
        seqs = seq_io.read(f)
        self.assertEqual(len(seqs), 5)
        self.assertEqual(len(seqs[1]), 210)
        f.close()

        f = data_stream('genbank/NT_019265.gb')
        seqs = genbank_io.read(f)
        self.assertEqual(len(seqs), 0)
        f.close()


if __name__ == '__main__':
    unittest.main()
