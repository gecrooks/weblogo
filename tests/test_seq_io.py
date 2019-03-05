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

from io import StringIO
import unittest

from weblogo.seq import protein_alphabet, nucleic_alphabet
from weblogo.seq_io import (fasta_io, clustal_io, plain_io, phylip_io, msf_io,
                            nbrf_io, nexus_io, stockholm_io, table_io, array_io,
                            genbank_io)
from weblogo import seq_io
from . import data_stream
from . import test_genbank_io


class test_seq_io(unittest.TestCase):
    def test_attr(self):
        f = seq_io.formats
        for i in f:
            i.names[0]  # Ensure every format has a name.

        seq_io.format_names()
        seq_io.format_extensions()

    def test_parse_clustal(self):
        with data_stream("clustal.aln") as f:
            seqs = seq_io.read(f)
        self.assertEqual(len(seqs), 7)
        self.assertEqual(seqs[1].name, "CATH_HUMAN")
        self.assertEqual(len(seqs[1]), 395)

    def test_parse_error(self):
        """ Wrong alphabet should throw a parsing error """
        with data_stream("clustal.aln") as f:
            self.assertRaises(ValueError,
                              seq_io.read, f, nucleic_alphabet)

    def test_parse_clustal181(self):
        with data_stream("clustal181.aln") as f:
            seq_io.read(f, protein_alphabet)

    def test_parse_clustal_glualign(self):
        with data_stream("clustal_glualign.aln") as f:
            seq_io.read(f, nucleic_alphabet)

    def test_parse_clustalw182(self):
        with data_stream("clustalw182.aln") as f:
            seq_io.read(f, protein_alphabet)

    def test_read_example_array(self):
        f = StringIO(array_io.example)
        seqs = seq_io.read(f)
        # print seqs
        self.assertEqual(len(seqs), 8)
        self.assertEqual(seqs[0].name, None)
        self.assertEqual(len(seqs[1]), 60)

    def test_read_fasta(self):
        f = StringIO(fasta_io.example)
        seqs = seq_io.read(f)
        # print seqs
        self.assertEqual(len(seqs), 3)
        self.assertEqual(seqs[0].description, "Lamprey GLOBIN V - SEA LAMPREY")
        self.assertEqual(len(seqs[1]), 231)

    def test_parse_globin_fasta(self):
        with data_stream("globin.fa") as f:
            seqs = seq_io.read(f)
        self.assertEqual(len(seqs), 56)

    def test_parser_extensions(self):
        # Test that the list of extension is a list.
        # Very easy with one extension list to write ('txt') rather than ('txt',)
        for p in seq_io._parsers:
            self.assertTrue(type(p.extensions) == tuple)

    def test_parser_names(self):
        # Same for names
        for p in seq_io._parsers:
            self.assertTrue(type(p.names) == tuple)

    def test_parsers(self):
        # seq_io._parsers is an ordered  list of sequence parsers that are
        # tried, in turn, on files of unknown format. Each parser must raise
        # an exception when fed a format further down the list.

        # We may include examples here for parsers that are not currently in
        # seq_io._parsers

        # TODO: Refactor these examples as with test_genbank_io.examples()
        # TODO: Then test that each example works with read() and iterseq()
        # TODO: Also autotest Write and writeseq, where available.

        fasta_examples = (StringIO(fasta_io.example),
                          data_stream("globin.fa"))

        clustal_examples = (StringIO(clustal_io.example),
                            data_stream("clustal.aln"),
                            data_stream("clustal181.aln"),
                            data_stream("clustal_glualign.aln"),
                            data_stream("clustalw182.aln"),
                            )
        plain_examples = (StringIO(plain_io.example),)
        phylip_examples = (
            data_stream("phylip_test_1.phy"),
            data_stream('phylip_test_2.phy'),
            data_stream('phylip_test_3.phy'),
            data_stream('phylip_test_4.phy'),
            data_stream('phylip_test_5.phy'),
            data_stream('dna.phy'),
        )
        msf_examples = (
            data_stream("dna.msf"),
            data_stream("cox2.msf"),
            data_stream("1beo.msf"),
        )
        nbrf_examples = (
            data_stream('cox2.nbrf'),
            data_stream('crab.nbrf'),
            data_stream('dna.pir'),
            data_stream('rhod.pir'),
            data_stream('protein.pir'),
        )
        nexus_examples = (
            data_stream("nexus/protein.nex"),
            data_stream("nexus/dna.nex"),
        )
        stockholm_examples = (
            StringIO(stockholm_io.example),
            data_stream("pfam_example.txt"),
            data_stream("pfam.txt"),
        )
        table_examples = (
            StringIO(table_io.example),
        )
        array_examples = (StringIO(array_io.example),)

        examples = {
            fasta_io: fasta_examples,
            clustal_io: clustal_examples,
            plain_io: plain_examples,
            phylip_io: phylip_examples,
            msf_io: msf_examples,
            nbrf_io: nbrf_examples,
            nexus_io: nexus_examples,
            stockholm_io: stockholm_examples,
            table_io: table_examples,
            array_io: array_examples,
            genbank_io: test_genbank_io.examples()
        }

        parsers = seq_io._parsers

        for i in range(0, len(parsers)):
            for j in range(i + 1, len(parsers)):
                # Check that parser[i] cannot read files intended for parser[j] (where j>i)
                for f in examples[parsers[j]]:
                    # print parsers[i].names[0], parsers[j].names[0]
                    f.seek(0)
                    self.assertRaises(ValueError, parsers[i].read, f)

        # When fed an empty file, the parser should either raise a ValueError
        # or return an empty SeqList
        e = StringIO()
        for p in seq_io._parsers:
            try:
                s = p.read(e)
                self.assertEqual(len(s), 0)
            except ValueError:
                pass

        for e in examples.values():
            for f in e:
                f.close()


if __name__ == '__main__':
    unittest.main()
