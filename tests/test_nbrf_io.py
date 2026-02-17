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

from weblogo.seq import dna_alphabet, protein_alphabet
from weblogo.seq_io import clustal_io, nbrf_io, plain_io

from . import data_stream


def test_parse_cox2() -> None:
    f = data_stream("cox2.nbrf")
    seqs = nbrf_io.read(f)
    assert len(seqs) == 5
    assert len(seqs[1]) == 210
    assert (
        str(seqs[0])
        == "MAFILSFWMIFLLDSVIVLLSFVCFVCVWICALLFSTVLLVSKLNNIYCTWDFTASKFIDVYWFTIGGMFSLG"
        "LLLRLCLLLYFGHLNFVSFDLCKVVGFQWYWVYFIFGETTIFSNLILESDYMIGDLRLLQCNHVLTLLSLVIY"
        "KLWLSAVDVIHSFAISSLGVKVENLVAVMK"
    )
    assert seqs[0].alphabet == protein_alphabet
    f.close()


def test_parse_crab() -> None:
    f = data_stream("crab.nbrf")
    seqs = nbrf_io.read(f)
    assert seqs[0].alphabet == protein_alphabet
    assert len(seqs) == 9
    assert seqs[2].name == "CRAB_CHICK"
    assert seqs[2].description == "ALPHA CRYSTALLIN B CHAIN (ALPHA(B)-CRYSTALLIN)."
    f.close()


def test_parse_dna() -> None:
    f = data_stream("dna.pir")
    seqs = nbrf_io.read(f)
    assert seqs[0].alphabet == dna_alphabet
    assert len(seqs) == 10
    f.close()


def test_parse_examples() -> None:
    f = data_stream("rhod.pir")
    seqs = nbrf_io.read(f)
    assert seqs[0].alphabet == protein_alphabet
    assert len(seqs) == 3
    f.close()


def test_parse_protein() -> None:
    f = data_stream("protein.pir")
    seqs = nbrf_io.read(f)
    assert seqs[0].alphabet == protein_alphabet
    assert len(seqs) == 10
    f.close()


def test_parse_clustal_fail() -> None:
    # should fail with parse error
    f = StringIO(clustal_io.example)
    with pytest.raises(ValueError):
        nbrf_io.read(f, protein_alphabet)


def test_parse_plain_fail() -> None:
    # should fail with parse error
    f = StringIO(plain_io.example)
    with pytest.raises(ValueError):
        nbrf_io.read(f)


def test_pir_file_from_clustal() -> None:
    f = data_stream("clustalw.pir")
    seqs = nbrf_io.read(f)
    assert len(seqs) == 2
    assert (
        seqs[1].endswith(
            "C-AATC-G-CAATG-G--CTTGAACCGGGTAAAAGTCGT-A----------------------------------------"
            "-----------------------------------------"
        )
        == True
    )
    f.close()


def test_parse_examples_alphabet() -> None:
    f = data_stream("rhod.pir")
    seqs = nbrf_io.read(f, alphabet=protein_alphabet)
    assert seqs[0].alphabet == protein_alphabet
    assert len(seqs) == 3
    f.close()
