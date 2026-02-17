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

import pytest
from io import StringIO

from weblogo.seq import nucleic_alphabet
from weblogo.seq_io import msf_io, plain_io

from . import data_stream


def test_msf_io_parse_msf() -> None:
    f = data_stream("dna.msf")
    seqs = msf_io.read(f)
    assert len(seqs) == 10
    assert seqs[1].name == "Carp"
    assert len(seqs[1]) == 705
    assert str(seqs[2][0:10]) == "ATGGCCAACC"
    f.close()


def test_msf_io_parse_msf2() -> None:
    f = data_stream("cox2.msf")
    seqs = msf_io.read(f)
    assert len(seqs) == 5
    assert seqs[1].name == "cox2_crifa"
    assert len(seqs[1]) == 166
    assert str(seqs[2][0:10]) == "MSFILTFWMI"
    f.close()


def test_msf_io_parse_1beo() -> None:
    f = data_stream("1beo.msf")
    msf_io.read(f)
    f.close()


# Wrong alphabet should throw a parsing error


def test_msf_io_parse_error() -> None:
    f = data_stream("cox2.msf")
    with pytest.raises(ValueError):
        msf_io.read(f, nucleic_alphabet)
    f.close()


def test_msf_io_parse_fasta_fail2() -> None:
    # should fail with parse error
    f = data_stream("globin.fa")
    with pytest.raises(ValueError):
        msf_io.read(f)
    f.close()


def test_msf_io_parse_plain_fail() -> None:
    # should fail with parse error
    f = StringIO(plain_io.example)
    with pytest.raises(ValueError):
        msf_io.read(f)
    f.close()


def test_msf_io_parse_phylip_fail() -> None:
    # should fail with parse error
    f = data_stream("phylip_test_2.phy")
    with pytest.raises(ValueError):
        msf_io.read(f)
    f.close()


def test_msf_io_iter() -> None:
    f = data_stream("1beo.msf")
    for seq in msf_io.iterseq(f):
        pass
    f.close()
