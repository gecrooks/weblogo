

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

from io import StringIO

import pytest

import weblogo.seq_io.intelligenetics_io as ig_io
from weblogo.seq import rna_alphabet
from . import data_stream


def test_read():
    f = StringIO(ig_io.example)
    seqs = ig_io.read(f)

    assert len(seqs) == 2
    assert seqs[0].description == "H.sapiens fau mRNA, 518 bases"
    assert seqs[1].name == "HSFAU1"
    assert len(seqs[1]) == 299


def test_read2():
    f = data_stream("intelligenetics.txt")
    seqs = ig_io.read(f)
    assert len(seqs[0]) == 518
    assert len(seqs[1]) == 2016
    f.close()


def test_write_seq():
    f = StringIO(ig_io.example)
    seqs = ig_io.read(f)
    fout = StringIO()
    ig_io.write(fout, seqs)
    fout.seek(0)
    seqs2 = ig_io.read(fout)
    assert seqs == seqs2


def test_read_fail():
    with pytest.raises(ValueError):
        f = data_stream("intelligenetics.txt")
        ig_io.read(f, alphabet=rna_alphabet)


def test_write_seq_fail():
    f = StringIO(ig_io.example)
    seqs = ig_io.read(f)
    seqs[0].name = ''
    fout = StringIO()
    with pytest.raises(ValueError):
        ig_io.write(fout, seqs)
