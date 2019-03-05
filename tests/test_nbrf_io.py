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

from io import StringIO

import unittest

from weblogo.seq import protein_alphabet, dna_alphabet
from weblogo.seq_io import nbrf_io, clustal_io, plain_io
from . import data_stream


class test_nbrf_io(unittest.TestCase):
    def test_parse_cox2(self):
        f = data_stream('cox2.nbrf')
        seqs = nbrf_io.read(f)
        self.assertEqual(len(seqs), 5)
        self.assertEqual(len(seqs[1]), 210)
        self.assertEqual(str(seqs[0]),
                         "MAFILSFWMIFLLDSVIVLLSFVCFVCVWICALLFSTVLLVSKLNNIYCTWDFTASKFIDVYWFTIGGMFSLG"
                         "LLLRLCLLLYFGHLNFVSFDLCKVVGFQWYWVYFIFGETTIFSNLILESDYMIGDLRLLQCNHVLTLLSLVIY"
                         "KLWLSAVDVIHSFAISSLGVKVENLVAVMK")
        self.assertEqual(seqs[0].alphabet, protein_alphabet)
        f.close()

    def test_parse_crab(self):
        f = data_stream('crab.nbrf')
        seqs = nbrf_io.read(f)
        self.assertEqual(seqs[0].alphabet, protein_alphabet)
        self.assertEqual(len(seqs), 9)
        self.assertEqual(seqs[2].name, "CRAB_CHICK")
        self.assertEqual(seqs[2].description,
                         "ALPHA CRYSTALLIN B CHAIN (ALPHA(B)-CRYSTALLIN).")
        f.close()

    def test_parse_dna(self):
        f = data_stream('dna.pir')
        seqs = nbrf_io.read(f)
        self.assertEqual(seqs[0].alphabet, dna_alphabet)
        self.assertEqual(len(seqs), 10)
        f.close()

    def test_parse_examples(self):
        f = data_stream('rhod.pir')
        seqs = nbrf_io.read(f)
        self.assertEqual(seqs[0].alphabet, protein_alphabet)
        self.assertEqual(len(seqs), 3)
        f.close()

    def test_parse_protein(self):
        f = data_stream('protein.pir')
        seqs = nbrf_io.read(f)
        self.assertEqual(seqs[0].alphabet, protein_alphabet)
        self.assertEqual(len(seqs), 10)
        f.close()

    def test_parse_clustal_fail(self):
        # should fail with parse error
        f = StringIO(clustal_io.example)
        self.assertRaises(ValueError,
                          nbrf_io.read, f, protein_alphabet)

    def test_parse_plain_fail(self):
        # should fail with parse error
        f = StringIO(plain_io.example)
        self.assertRaises(ValueError,
                          nbrf_io.read, f)

    def test_pir_file_from_clustal(self):
        f = data_stream('clustalw.pir')
        seqs = nbrf_io.read(f)
        self.assertEqual(len(seqs), 2)
        self.assertEqual(seqs[1].endswith(
                'C-AATC-G-CAATG-G--CTTGAACCGGGTAAAAGTCGT-A----------------------------------------'
                '-----------------------------------------'),
                True)
        f.close()

    def test_parse_examples_alphabet(self):
        f = data_stream('rhod.pir')
        seqs = nbrf_io.read(f, alphabet=protein_alphabet)
        self.assertEqual(seqs[0].alphabet, protein_alphabet)
        self.assertEqual(len(seqs), 3)
        f.close()


if __name__ == '__main__':
    unittest.main()
