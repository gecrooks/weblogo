from io import StringIO

import pytest

from weblogo.seq import protein_alphabet
from weblogo.seq_io import nexus_io, plain_io, clustal_io
from . import data_stream


def test_read():
    f = data_stream("nexus/protein.nex")
    seqs = nexus_io.read(f)
    # print seqs
    assert len(seqs) == 10
    assert seqs[0].name == "Cow"
    assert len(seqs[1]) == 234
    assert str(seqs[0][0:10]) == 'MAYPMQLGFQ'
    f.close()


def test_parse_StringIO():
    # Bio.Nexus cannot read from a StringIO object.
    f0 = data_stream("nexus/protein.nex")
    f = StringIO(f0.read())
    nexus_io.read(f)
    f0.close()


def test_parse_fasta_fail():
    with pytest.raises(ValueError):
        with data_stream("globin.fa") as f:
            nexus_io.read(f)


def test_parse_clustal_fail():
    # should fail with parse error
    with pytest.raises(ValueError):
        f = StringIO(clustal_io.example)
        nexus_io.read(f, protein_alphabet)


def test_parse_plain_fail():
    # should fail with parse error
    f = StringIO(plain_io.example)
    with pytest.raises(ValueError):
        nexus_io.read(f)


def test_iterseq():
    f = data_stream("nexus/protein.nex")
    for seq in nexus_io.iterseq(f):
        pass


def test_read_alphabet():
    f = data_stream("nexus/protein.nex")
    seqs = nexus_io.read(f, alphabet=protein_alphabet)
    # print seqs
    assert len(seqs) == 10
