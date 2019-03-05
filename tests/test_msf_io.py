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

from weblogo.seq import nucleic_alphabet
from weblogo.seq_io import msf_io, plain_io
from . import data_stream


class test_msf_io(unittest.TestCase):
    def test_parse_msf(self):
        f = data_stream("dna.msf")
        seqs = msf_io.read(f)
        self.assertEqual(len(seqs), 10)
        self.assertEqual(seqs[1].name, "Carp")
        self.assertEqual(len(seqs[1]), 705)
        self.assertEqual(str(seqs[2][0:10]), 'ATGGCCAACC')
        f.close()

    def test_parse_msf2(self):
        f = data_stream("cox2.msf")
        seqs = msf_io.read(f)
        self.assertEqual(len(seqs), 5)
        self.assertEqual(seqs[1].name, "cox2_crifa")
        self.assertEqual(len(seqs[1]), 166)
        self.assertEqual(str(seqs[2][0:10]), 'MSFILTFWMI')
        f.close()

    def test_parse_1beo(self):
        f = data_stream("1beo.msf")
        msf_io.read(f)
        f.close()

    """ Wrong alphabet should throw a parsing error """

    def test_parse_error(self):
        f = data_stream("cox2.msf")
        self.assertRaises(ValueError,
                          msf_io.read, f, nucleic_alphabet)
        f.close()

    def test_parse_fasta_fail2(self):
        # should fail with parse error
        f = data_stream("globin.fa")
        self.assertRaises(ValueError,
                          msf_io.read, f)
        f.close()

    def test_parse_plain_fail(self):
        # should fail with parse error
        f = StringIO(plain_io.example)
        self.assertRaises(ValueError,
                          msf_io.read, f)
        f.close()

    def test_parse_phylip_fail(self):
        # should fail with parse error
        f = data_stream("phylip_test_2.phy")
        self.assertRaises(ValueError,
                          msf_io.read, f)
        f.close()

    def test_iter(self):
        f = data_stream("1beo.msf")
        for seq in msf_io.iterseq(f):
            pass
        f.close()


if __name__ == '__main__':
    unittest.main()
