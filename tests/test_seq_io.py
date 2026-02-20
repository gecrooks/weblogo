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
from io import StringIO

from weblogo import seq_io
from weblogo.seq import nucleic_alphabet, protein_alphabet
from weblogo.seq_io import (
    array_io,
    clustal_io,
    fasta_io,
    genbank_io,
    msf_io,
    nbrf_io,
    phylip_io,
    plain_io,
    stockholm_io,
    table_io,
)

from . import data_ref, data_stream, test_genbank_io


def test_attr() -> None:
    f = seq_io.formats
    for i in f:
        i.names[0]  # Ensure every format has a name.

    seq_io.format_names()
    seq_io.format_extensions()


def test_parse_clustal() -> None:
    with data_ref("clustal.aln").open() as f:
        seqs = seq_io.read(f)
    assert len(seqs) == 7
    assert seqs[1].name == "CATH_HUMAN"
    assert len(seqs[1]) == 395


def test_parse_error() -> None:
    """Wrong alphabet should throw a parsing error"""
    with data_ref("clustal.aln").open() as f:
        with pytest.raises(ValueError):
            seq_io.read(f, nucleic_alphabet)


def test_parse_clustal181() -> None:
    with data_ref("clustal181.aln").open() as f:
        seq_io.read(f, protein_alphabet)


def test_parse_clustal_glualign() -> None:
    with data_ref("clustal_glualign.aln").open() as f:
        seq_io.read(f, nucleic_alphabet)


def test_parse_clustalw182() -> None:
    with data_ref("clustalw182.aln").open() as f:
        seq_io.read(f, protein_alphabet)


def test_read_example_array() -> None:
    f = StringIO(array_io.example)
    seqs = seq_io.read(f)
    # print seqs
    assert len(seqs) == 8
    assert seqs[0].name == ""
    assert len(seqs[1]) == 60


def test_read_fasta() -> None:
    f = StringIO(fasta_io.example)
    seqs = seq_io.read(f)
    # print seqs
    assert len(seqs) == 3
    assert seqs[0].description == "Lamprey GLOBIN V - SEA LAMPREY"
    assert len(seqs[1]) == 231


def test_parse_globin_fasta() -> None:
    with data_ref("globin.fa").open() as f:
        seqs = seq_io.read(f)
    assert len(seqs) == 56


def test_get_parsers_unknown_extension() -> None:
    """Extension not in fnames or fext: best_guess stays as default."""
    fin = StringIO()
    fin.name = "data.xyz"
    parsers = seq_io._get_parsers(fin)
    assert parsers[0] is seq_io._parsers[0]


def test_get_parsers_non_parser_format() -> None:
    """Extension maps to a format not in _parsers (e.g. intelligenetics_io)."""
    from weblogo.seq_io import intelligenetics_io

    fin = StringIO()
    fin.name = "data.ig"
    parsers = seq_io._get_parsers(fin)
    assert parsers[0] is intelligenetics_io


def test_parser_extensions() -> None:
    # Test that the list of extension is a list.
    # Very easy with one extension list to write ('txt') rather than ('txt',)
    for p in seq_io._parsers:
        assert type(p.extensions) is tuple


def test_parser_names() -> None:
    # Same for names
    for p in seq_io._parsers:
        assert type(p.names) is tuple


def test_parsers() -> None:
    # seq_io._parsers is an ordered  list of sequence parsers that are
    # tried, in turn, on files of unknown format. Each parser must raise
    # an exception when fed a format further down the list.

    # We may include examples here for parsers that are not currently in
    # seq_io._parsers

    # TODO: Refactor these examples as with test_genbank_io.examples()
    # TODO: Then test that each example works with read() and iterseq()
    # TODO: Also autotest Write and writeseq, where available.

    fasta_examples = (StringIO(fasta_io.example), data_stream("globin.fa"))

    clustal_examples = (
        StringIO(clustal_io.example),
        data_stream("clustal.aln"),
        data_stream("clustal181.aln"),
        data_stream("clustal_glualign.aln"),
        data_stream("clustalw182.aln"),
    )
    plain_examples = (StringIO(plain_io.example),)
    phylip_examples = (
        data_stream("phylip_test_1.phy"),
        data_stream("phylip_test_2.phy"),
        data_stream("phylip_test_3.phy"),
        data_stream("phylip_test_4.phy"),
        data_stream("phylip_test_5.phy"),
        data_stream("dna.phy"),
    )
    msf_examples = (
        data_stream("dna.msf"),
        data_stream("cox2.msf"),
        data_stream("1beo.msf"),
    )
    nbrf_examples = (
        data_stream("cox2.nbrf"),
        data_stream("crab.nbrf"),
        data_stream("dna.pir"),
        data_stream("rhod.pir"),
        data_stream("protein.pir"),
    )
    stockholm_examples = (
        StringIO(stockholm_io.example),
        data_stream("pfam_example.txt"),
        data_stream("pfam.txt"),
    )
    table_examples = (StringIO(table_io.example),)
    array_examples = (StringIO(array_io.example),)

    examples = {
        fasta_io: fasta_examples,
        clustal_io: clustal_examples,
        plain_io: plain_examples,
        phylip_io: phylip_examples,
        msf_io: msf_examples,
        nbrf_io: nbrf_examples,
        stockholm_io: stockholm_examples,
        table_io: table_examples,
        array_io: array_examples,
        genbank_io: test_genbank_io.examples(),
    }

    parsers = seq_io._parsers

    for i in range(0, len(parsers)):
        for j in range(i + 1, len(parsers)):
            # Check that parser[i] cannot read files intended for parser[j] (where j>i)
            for f in examples[parsers[j]]:
                # print parsers[i].names[0], parsers[j].names[0]
                f.seek(0)
                with pytest.raises(ValueError):
                    parsers[i].read(f)

    # When fed an empty file, the parser should either raise a ValueError
    # or return an empty SeqList
    e = StringIO()
    for p in seq_io._parsers:
        try:
            s = p.read(e)
            assert len(s) == 0
        except ValueError:
            pass

    for e in examples.values():
        for f in e:
            f.close()
