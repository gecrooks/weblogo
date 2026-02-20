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

import pytest

from weblogo import seq_io
from weblogo.seq import (
    Alphabet,
    Seq,
    SeqList,
    dna,
    dna_alphabet,
    generic_alphabet,
    nucleic_alphabet,
    protein,
    protein_alphabet,
    rna,
    unambiguous_dna_alphabet,
    unambiguous_protein_alphabet,
)

from . import data_ref


# --- Tests from test_alphabet ---


def test_create_alphabet() -> None:
    # Alphabet contains repeated character
    with pytest.raises(ValueError):
        Alphabet("alphabet")

    # Alphabet contains null character
    with pytest.raises(ValueError):
        Alphabet("alph\x00")

    Alphabet("alphbet")


def test_alphabet_alphabetic() -> None:
    a = Alphabet("alphbet")
    assert a.alphabetic("alphbet")
    assert not a.alphabetic("alphbetX")


def test_alphabet_ord() -> None:
    a = generic_alphabet
    for i, c in enumerate(a):
        assert a.ord(c) == i

    a = Alphabet("alph")
    assert 2 == a.ord("p")


def test_alphabet_chr() -> None:
    a = generic_alphabet
    for i, c in enumerate(a):
        assert ord(a.chr(i)) == i + 32

    a = Alphabet("alph")
    assert "h" == a.chr(3)


def test_alphabet_ords() -> None:
    a = Alphabet("alph")
    assert 0 == a.ords("alphalph")[4]

    a = generic_alphabet
    o = a.ords(str(a))
    for i, c in enumerate(o):
        assert c == i


def test_alphabet_chrs() -> None:
    a = Alphabet("alph")
    assert Seq("ppla", a) == a.chrs((2, 2, 1, 0))


def test_alphabet_none() -> None:
    a1 = Alphabet(None)
    assert a1 == generic_alphabet


def test_alphabet_create_from_alphabet() -> None:
    """If we pass an alphabet to the constuctor, it's passed
    right back"""
    a1 = Alphabet("kjdahf")
    a2 = Alphabet(str(a1))
    assert a1 == a2

    assert not (a1 == "not an alphabet")


def test_alphabet_repr() -> None:
    a = Alphabet("kjdahf")
    repr(a)
    str(a)


def test_alphabet_str() -> None:
    assert str(Alphabet("kjdahf")) == "kjdahf"


def test_alphabet_letters() -> None:
    assert Alphabet("kjdahf").letters() == "kjdahf"


def test_alphabet_normalize() -> None:
    a = Alphabet("ABCDE")
    s = "aBbc"
    n = a.normalize(s)
    assert str(n) == "ABBC"

    with pytest.raises(ValueError):
        a.normalize("aslkfdnnr33")


def test_alphabet_alt() -> None:
    alt = tuple(zip("12345", "ABCED"))
    alpha = Alphabet("ABCDE", alt)
    assert alpha.ord("A") == 0
    for a, c in alt:
        assert alpha.ord(a) == alpha.ord(c)


def test_alphabet_which_alphabet() -> None:
    a = Alphabet.which(Seq("ARNDCQEGHILKMFPSTWYVX"))
    assert a == unambiguous_protein_alphabet

    test_cases = (
        ("cap.fa", unambiguous_dna_alphabet),
        ("cox2.msf", unambiguous_protein_alphabet),
        ("Rv3829c.fasta", unambiguous_protein_alphabet),
        ("chain_B.fasta", unambiguous_protein_alphabet),
    )
    for filename, expected_alphabet in test_cases:
        with data_ref(filename).open() as f:
            seqs = seq_io.read(f)
            assert Alphabet.which(seqs) == expected_alphabet


# --- Tests from test_seq ---


def test_create_seq() -> None:
    assert Seq("alphabet", Alphabet("alphbet"))
    with pytest.raises(ValueError):
        Seq("not alphabetic", "alphabet")

    a = (
        "Any printable Ascii character `1234567890-=~!@#$%^&*()_+{}|[]\\:;'<>?,./QWERTYUIOPASD"
        "FGHJKLZXCVBNMqwertyuiopasdfghjklzxcvbnm "
    )

    for x in a:
        assert x in generic_alphabet
    assert Seq(a, generic_alphabet)
    with pytest.raises(ValueError):
        Seq("Not zero. \x00", generic_alphabet)


def test_std_alphabets() -> None:
    s = Seq("dskjjfskdjbfKJJSJKSKJDjk123u768erbm", generic_alphabet)
    s = Seq("ARNDCQEGHILKMFPSTWYVX", protein_alphabet)
    with pytest.raises(ValueError):
        Seq("1234", protein_alphabet)
    s = Seq("AGTCAGCTACGACGCGC", dna_alphabet)
    assert str(s[1]) == "G"


def test_seq_ords() -> None:
    a = Alphabet("ABC")
    s = Seq("ABCCBA", a)
    assert list(s.ords()) == [0, 1, 2, 2, 1, 0]


def test_seq_tally() -> None:
    s = Seq("AGTCAGCTACGACGCGC", unambiguous_dna_alphabet)
    c = s.tally()
    assert len(unambiguous_dna_alphabet) == len(c)
    assert list(c) == [4, 6, 5, 2]


def test_seq_tally_nonalphabetic() -> None:
    s = Seq("AGTCAGCTACGACGCGC", dna_alphabet)
    c = s.tally(Alphabet("AC"))
    assert 2 == len(c)
    assert list(c) == [4, 6]


def test_seq_words() -> None:
    s = Seq("AGTCAGCTACGACGcgcx", dna_alphabet)
    w = list(s.words(2, unambiguous_dna_alphabet))
    assert len(w) == len(s) - 2
    assert w == [
        "AG",
        "GT",
        "TC",
        "CA",
        "AG",
        "GC",
        "CT",
        "TA",
        "AC",
        "CG",
        "GA",
        "AC",
        "CG",
        "GC",
        "CG",
        "GC",
    ]

    assert list(s.words(len(s), unambiguous_dna_alphabet)) == []
    assert list(s.words(len(s) - 1, unambiguous_dna_alphabet)) == [
        "AGTCAGCTACGACGCGC",
    ]

    w = list(s.words(200, unambiguous_dna_alphabet))


def test_seq_words2() -> None:
    s = Seq("AGTCAGCTACGACGCGC", unambiguous_dna_alphabet)
    wc = s.word_count(2)
    count = list(zip(*wc))[1]
    assert count == (2, 2, 1, 3, 1, 1, 3, 1, 1, 1)


def test_seq_getslice() -> None:
    s = Seq("AGTCAGCTACGACGCGC", dna_alphabet)
    slice = s[2:4]
    assert s.alphabet == slice.alphabet


def test_seq_create_annotated() -> None:
    s = "ACGTURYSWKMBDHVNACGTURYSWKMBDHVNAAAAA"
    a = Seq(s, nucleic_alphabet, name="ID", description="DESCRIPTION")
    assert a.name == "ID"
    assert a.description == "DESCRIPTION"
    assert s == str(a)


def test_seq_add() -> None:
    s1 = Seq("AAAA", dna_alphabet)
    s2 = Seq("TTTT", dna_alphabet)

    s3 = s1 + s2
    assert s3.alphabet == dna_alphabet
    assert s3 == Seq("AAAATTTT", dna_alphabet)

    assert s3 == Seq("AAAATTTT", dna_alphabet)
    assert s3 != Seq("AAAATTTT", protein_alphabet)
    assert s3 != "not a seq"

    s4 = "AA"
    s5 = s4 + s1
    s6 = s1 + s4
    assert s5.alphabet == s6.alphabet
    assert s5 == s6

    assert s5 == s6
    assert not (s5 != s6)


def test_seq_join() -> None:
    s1 = Seq("AAAA", dna_alphabet)
    s2 = Seq("TTTT", dna_alphabet)
    s3 = Seq("GGGG")
    s0 = Seq("", dna_alphabet)

    j = s0.join([s1, s2, s3])
    assert j == Seq("AAAATTTTGGGG", dna_alphabet)


def test_seq_repr() -> None:
    s1 = Seq("AAAA", dna_alphabet)
    repr(s1)


def test_seq_str() -> None:
    s1 = Seq("AGCTA", dna_alphabet)
    assert str(s1) == "AGCTA"
    # Uncased alpahebt
    assert str(Seq("AgcTA", dna_alphabet)) == "AgcTA"


def test_seq_tostring() -> None:
    assert Seq("AgcTAAAA", dna_alphabet).tostring() == "AgcTAAAA"


def test_seq_reverse() -> None:
    s = Seq("ACGT", dna_alphabet)
    assert s == s.reverse().reverse()
    assert s.reverse() == Seq("TGCA", dna_alphabet)


def test_seq_translate() -> None:
    s = dna("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
    p = s.translate()
    assert str(p) == "AIVMGR*KGAR"
    p.back_translate()


def test_seq_reverse_complement() -> None:
    s = dna("GGGGaaaaaaaatttatatat")
    rc = s.reverse_complement()
    assert rc == dna("atatataaattttttttCCCC")

    assert dna("ACGTRYSWKMBDHVN").reverse_complement() == dna("NBDHVKMWSRYACGT")

    with pytest.raises(ValueError):
        protein("G").reverse_complement()


def test_seq_ungap() -> None:
    s = Seq("T-T", dna_alphabet).ungap()
    assert str(s) == "TT"
    s = Seq("T-~---T...~~~--", dna_alphabet).ungap()
    assert str(s) == "TT"


def test_seq_mask() -> None:
    s = dna("AAaaaaAAA").mask()
    assert str(s) == "AAXXXXAAA"


def test_seq_shortcuts() -> None:
    protein("gGGGGG-PPPPP")
    dna("ACGTRYSWKMBDHVN")
    dna("t")
    dna("actgta")
    dna("acgtrysw-kmb-dhvn")

    rna("ACUUUUU")


def test_seq_casechange() -> None:
    s1 = dna("ACGTRYSWKMBDHVN")
    s2 = s1.lower().upper()
    assert s1 == s2


def test_seq_slice() -> None:
    s1 = dna("ACGTRYSWKMBDHVN")
    s2 = s1[4:6]
    assert s2 == dna("RY")


# --- Tests from test_seqlist ---


def test_seqlist_create() -> None:
    # 1234567890123456789012345678
    s0 = Seq("ACGTURYBDHVNACGTURYSWKMBDHVN", nucleic_alphabet)
    s1 = Seq("ACGTURYSWKMBDHVNACGTURYSWKMBDHVN", nucleic_alphabet)
    s2 = Seq("ACGTURSWKMBDHVNACGTURKMBDHVN", nucleic_alphabet)
    seqs = SeqList([s0, s1, s2])

    assert len(seqs) == 3


def test_seqlist_create_annotated() -> None:
    s0 = Seq("ACGTURYBDHVNACGTURYSWKMBDHVN", nucleic_alphabet)
    s1 = Seq("ACGTURYSWKMBDHVNACGTURYSWKMBDHVN", nucleic_alphabet)
    s2 = Seq("ACGTURSWKMBDHVNACGTURKMBDHVN", nucleic_alphabet)

    seqs = SeqList(
        [s0, s1, s2], alphabet=nucleic_alphabet, name="alsdf", description="a"
    )
    assert seqs.name == "alsdf"
    assert seqs.description == "a"
    assert seqs.alphabet == nucleic_alphabet


def test_seqlist_ords() -> None:
    s0 = Seq("ACGTURYBDHVNACGTURYSWKMBDHVN", nucleic_alphabet)
    s1 = Seq("ACGTURYSDHVNACGTURYSWKMBDHVN", nucleic_alphabet)
    s2 = Seq("ACGTURSWKMBDHVNACGTURKMBDHVN", nucleic_alphabet)
    seqs = SeqList([s0, s1, s2], nucleic_alphabet)
    seqs.ords()

    # Use a different alphabet
    seqs.ords(nucleic_alphabet)

    # No alphabet
    seqs3 = SeqList([s0, s1, s2])
    seqs3.ords(alphabet=Alphabet("ABC"))

    # Fail if no alphabet
    with pytest.raises(ValueError):
        seqs3.ords()


def test_seqlist_isaligned() -> None:
    a = Alphabet("ABCD")

    s0 = Seq("ABCDD", a)
    s1 = Seq("AAAAD", a)
    s2 = Seq("AAABD", a)
    s3 = Seq("AAACD", a)
    seqs = SeqList([s0, s1, s2, s3], a)
    assert seqs.isaligned()

    seqs = SeqList([s0, s1, s2, s3], Alphabet("ABCDE"))
    assert not seqs.isaligned()


def test_seqlist_profile() -> None:
    a = Alphabet("ABCD")

    s0 = Seq("ABCDD", a)
    s1 = Seq("AAAAD", a)
    s2 = Seq("AAABD", a)
    s3 = Seq("AAACD", a)

    seqs = SeqList([s0, s1, s2, s3], a)

    tally = seqs.profile()

    assert list(tally[0]) == [4, 0, 0, 0]
    assert list(tally[1]) == [3, 1, 0, 0]
    assert list(tally[2]) == [3, 0, 1, 0]
    assert list(tally[3]) == [1, 1, 1, 1]
    assert list(tally[4]) == [0, 0, 0, 4]

    assert tally[4, "D"] == 4

    seqs = SeqList([Seq("AAACD", a), Seq("AAACDA", a)], a)
    with pytest.raises(ValueError):
        seqs.profile()

    seqs = SeqList([Seq("AAACD", a), Seq("AAACD", a)])
    with pytest.raises(ValueError):
        seqs.profile()


def test_seqlist_tally() -> None:
    # 1234567890123456789012345678
    s0 = Seq("ACTTT", nucleic_alphabet)
    s1 = Seq("ACCCC", nucleic_alphabet)
    s2 = Seq("GGGG", nucleic_alphabet)
    seqs = SeqList([s0, s1, s2], nucleic_alphabet)

    counts = seqs.tally()
    assert counts == [2, 5, 4, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    seqs = SeqList([Seq("AAACD", nucleic_alphabet), Seq("AAACD", nucleic_alphabet)])
    with pytest.raises(ValueError):
        seqs.tally()


def test_seqlist_create_empty() -> None:
    s0 = Seq("ACGTURYBDHVNACGTURYSWKMBDHVN", nucleic_alphabet)
    s1 = Seq("ACGTURYSWKMBDHVNACGTURYSWKMBDHVN", nucleic_alphabet)
    s2 = Seq("ACGTURSWKMBDHVNACGTURKMBDHVN", nucleic_alphabet)

    seqs = SeqList()
    seqs.append(s0)
    seqs.extend((s1, s2))

    assert len(seqs) == 3
    assert type(seqs) is SeqList


def test_seqlist_repr() -> None:
    s0 = Seq("ACGTURYBDHVNACGTURYSWKMBDHVN", nucleic_alphabet)
    s1 = Seq("ACGTURYSWKMBDHVNACGTURYSWKMBDHVN", nucleic_alphabet)
    s2 = Seq("ACGTURSWKMBDHVNACGTURKMBDHVN", nucleic_alphabet)
    seqs = SeqList([s0, s1, s2])

    repr(seqs)


def test_seqlist_str() -> None:
    s0 = Seq("ACGTURYBDHVNACGTURYSWKMBDHVN", nucleic_alphabet)
    s1 = Seq("ACGTURYSWKMBDHVNACGTURYSWKMBDHVN", nucleic_alphabet)
    s2 = Seq("ACGTURSWKMBDHVNACGTURKMBDHVN", nucleic_alphabet)
    seqs = SeqList([s0, s1, s2])
    str(seqs)


# --- Standalone test ---


def test_alphabet_which_custom_alphabets() -> None:
    """Alphabet.which with an explicit alphabets list."""
    a = Alphabet("AC")
    b = Alphabet("GT")
    s = SeqList([Seq("AACCAA", generic_alphabet)], generic_alphabet)
    result = Alphabet.which(s, alphabets=[a, b])
    assert result == a


def test_seqlist_profile_explicit_alphabet() -> None:
    """SeqList.profile with an explicit alphabet argument."""
    a = Alphabet("ABCD")
    s0 = Seq("ABCDD", a)
    s1 = Seq("AAAAD", a)
    seqs = SeqList([s0, s1], a)
    tally = seqs.profile(a)
    assert list(tally[0]) == [2, 0, 0, 0]


def test_seqlist_profile_skip_non_alphabet_chars() -> None:
    """Profile with a narrower alphabet skips out-of-range chars."""
    wide = Alphabet("ABCD")
    narrow = Alphabet("AB")
    s0 = Seq("ABCD", wide)
    s1 = Seq("AABB", wide)
    seqs = SeqList([s0, s1], wide)
    tally = seqs.profile(narrow)
    # Column 0: A,A -> A=2, B=0. Column 1: B,A -> A=1, B=1
    # Column 2: C,B -> C not in narrow (skipped), B=1 -> A=0, B=1
    # Column 3: D,B -> D not in narrow (skipped), B=1 -> A=0, B=1
    assert list(tally[0]) == [2, 0]
    assert list(tally[1]) == [1, 1]
    assert list(tally[2]) == [0, 1]
    assert list(tally[3]) == [0, 1]


def test_bad_mask() -> None:
    with pytest.raises(ValueError):
        dna("AAaaaaAAA").mask(mask="ABC")


def test_seq_contains() -> None:
    s = dna("ACGT")
    assert "A" in s
    assert "X" not in s


def test_seq_hash() -> None:
    s1 = dna("ACGT")
    s2 = dna("ACGT")
    assert hash(s1) == hash(s2)
    d = {s1: "test"}
    assert d[s2] == "test"
