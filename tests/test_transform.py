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


import pytest

from weblogo.seq import (
    Seq,
    dna_alphabet,
    nucleic_alphabet,
    protein_alphabet,
    reduced_protein_alphabet,
)
from weblogo.transform import (
    GeneticCode,
    Transform,
    mask_low_complexity,
    reduced_protein_alphabets,
)


def test_mask_low_complexity_segging() -> None:
    before = (
        "mgnrafkshhghflsaegeavkthhghhdhhthfhvenhggkvalkthcgkylsigdhkqvylshhlhgdhslfhlehhg"
        "gkvsikghhhhyisadhhghvstkehhdhdttfeeiii".upper()
    )
    after = (
        "MGNRAFKSHHGHFLSAEGEAVxxxxxxxxxxxxxxxENHGGKVALKTHCGKYLSIGDHKQVYLSHHLHGDHSLFHLEHHGG"
        "KVSIKGHHHHYISADHHGHVSTKEHHDHDTTFEEIII".upper()
    )

    bseq = Seq(before, protein_alphabet)
    aseq = Seq(after, protein_alphabet)
    xseq = Seq("X" * len(bseq), protein_alphabet)

    sseq = mask_low_complexity(bseq)
    assert aseq == sseq

    # Nothing should be segged
    sseq = mask_low_complexity(bseq, 12, 0, 0)
    assert bseq == sseq

    # Everthing should be segged
    sseq = mask_low_complexity(bseq, 12, 4.3, 4.3)
    assert sseq == xseq

    mask_low_complexity(bseq, 100000, 4.3, 4.3)


def test_mask_low_complexity_seg_invalid() -> None:
    seq = Seq("KTHCGKYLSIGDHKQVYLSHH", protein_alphabet)
    with pytest.raises(ValueError):
        mask_low_complexity(seq, 12, -1, 0)
    with pytest.raises(ValueError):
        mask_low_complexity(seq, -1, 0, 0)
    with pytest.raises(ValueError):
        mask_low_complexity(seq, 12, 1, 10)
    with pytest.raises(ValueError):
        mask_low_complexity(seq, 6, 12, 13)
    with pytest.raises(ValueError):
        mask_low_complexity(seq, 6, 2.0, 1.9)


def test_transform_transform() -> None:
    trans = Transform(
        Seq("ACGTURYSWKMBDHVN", nucleic_alphabet),
        Seq("ACGTTNNNNNNNNNNN", dna_alphabet),
    )
    s0 = Seq("AAAAAR", nucleic_alphabet)
    s1 = trans(s0)  # Callable ob
    assert s1.alphabet == dna_alphabet
    assert s1 == Seq("AAAAAN", dna_alphabet)

    s2 = Seq(str(protein_alphabet), protein_alphabet)
    with pytest.raises(ValueError):
        trans(s2)

    # def test_translations(self) -> None:

    #     s = Seq("ACGTURYSWKMBDHVNACGTURYSWKMBDHVN", nucleic_alphabet)
    #     s2 = dna_ext_to_std(s)
    #     s3 = Seq("ACGTTNNNNNNNNNNNACGTTNNNNNNNNNNN", dna_alphabet)
    #     assert s2 == s3


def test_transform_reduced_protein_alphabets() -> None:
    seq = Seq(
        "ENHGGKVALKTHCGKYLSIGDHKQVYLSHHLHGDHSLFHLEHHGGKVSIKGHHHHYISADHHGHVSTKEHHDHDT"
        "TFEEIII",
        reduced_protein_alphabet,
    )

    for t in reduced_protein_alphabets.values():
        t(seq)


def test_geneticcode_repr() -> None:
    for t in GeneticCode.std_list():
        r = repr(t)
        gc = eval(r)
        assert r == repr(gc)
        assert str(gc) == str(t)
        # print r
        # print t
        # print gc


def test_geneticcode_translate_std() -> None:
    dna = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
    t = GeneticCode.std()
    s = t.translate(dna)
    assert str(s) == "AIVMGR*KGAR"

    for t in GeneticCode.std_list():
        t.translate(dna)


def test_geneticcode_translate() -> None:
    # Ref: http://lists.open-bio.org/pipermail/biopython/2006-March/002960.html

    cft = (
        (
            "Vertebrate Mitochondrial",
            "GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA",
            "AIVMGRWKGAR",
        ),
        (11, "CAAGGCGTCGAAYAGCTTCAGGAACAGGAC", "QGVE?LQEQD"),
        (1, "GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA", "AIVMGR*KGAR"),
    )

    for code, dna, protein in cft:
        c = GeneticCode.by_name(code)  # type: ignore
        trans = c.translate(Seq(dna))
        assert str(trans) == protein

    with pytest.raises(ValueError):
        GeneticCode.by_name("not_a_name")


def test_geneticcode_back_translate() -> None:
    prot = Seq("ACDEFGHIKLMNPQRSTVWY*")
    t = GeneticCode.std()
    t.table["CGA"]  # type: ignore
    s = t.back_translate(prot)
    assert str(prot) == str(t.translate(s))

    GeneticCode.std().back_table
