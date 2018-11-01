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

from weblogo.seq import nucleic_alphabet, protein_alphabet, rna_alphabet
from weblogo.seq_io import stockholm_io, clustal_io, fasta_io
from . import data_stream


class test_stockholm_io(unittest.TestCase):
    def test_parse1(self):
        with data_stream("pfam.txt") as f:
            seqs = stockholm_io.read(f)
        # self.assertEqual(len(seqs), 7)
        self.assertEqual(seqs[1].name, "O61132/1-232")
        self.assertEqual(len(seqs[1]), 265)

    def test_parse2(self):
        f = StringIO(stockholm_io.example)
        seqs = stockholm_io.read(f)
        self.assertEqual(len(seqs), 5)
        self.assertEqual(seqs[1].name, "O83071/259-312")
        self.assertEqual(len(seqs[1]), 43)

    def test_iterseq(self):
        with data_stream("pfam.txt") as f:
            for seq in stockholm_io.iterseq(f):
                pass

    # 12345678901234567890123456789012345678901234567890
    # 12*50 +6 = 606
    # QYVTVFYGVPAWRNATIPLFCATKNR.......DTWGTTQCLPDNDDYSE

    def test_parse3(self):
        with data_stream("pfam_example.txt") as f:
            seqs = stockholm_io.read(f)
        self.assertEqual(len(seqs), 24)
        self.assertEqual(seqs[5].name, "ENV_HV2BE/24-510")
        self.assertEqual(len(seqs[1]), 606)
        self.assertEqual(str(seqs[0][-6:]), 'TSRNKR')

    def test_parse_error(self):
        """ Wrong alphabet should throw a parsing error """
        f = StringIO(stockholm_io.example)
        self.assertRaises(ValueError,
                          clustal_io.read, f, nucleic_alphabet)

    def test_parse_fasta_fail(self):
        # should fail with parse error
        f = StringIO(stockholm_io.example)
        self.assertRaises(ValueError,
                          stockholm_io.read, f, rna_alphabet)

    def test_parse_alphabet_fail(self):
        # should fail with parse error
        f = StringIO(fasta_io.example)
        self.assertRaises(ValueError,
                          stockholm_io.read, f, protein_alphabet)

    def test_parse_fasta_fail2(self):
        # should fail with parse error
        with data_stream("globin.fa") as f:
            self.assertRaises(ValueError, stockholm_io.read, f)

    def test_parse_fail(self):
        # should fail with parse error
        examples = (
            StringIO(clustal_io.example),
        )

        for f in examples:
            self.assertRaises(ValueError,
                              stockholm_io.read, f)


if __name__ == '__main__':
    unittest.main()
