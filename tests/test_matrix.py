#!/usr/bin/env python

import numpy as np

from weblogo.matrix import AlphabeticArray, Motif
from weblogo.seq import protein_alphabet

from . import data_stream


def test_alphabeticarray_create() -> None:
    matrix = AlphabeticArray((protein_alphabet, protein_alphabet))
    matrix["A", "C"] = 10
    assert matrix[0, 1] == 10.0


def test_motif_read_transfac_alphabet_superset() -> None:
    with data_stream("transfac_matrix.txt") as f:
        Motif.read_transfac(f, alphabet="TCGA")

    # Supplied alphabet can be superset of defacto alphabet.
    # Reverts to defacto alphabet
    with data_stream("transfac_matrix.txt") as f:
        Motif.read_transfac(f, alphabet="TCGAXYZ")


def test_motif_read_transfac() -> None:
    f = data_stream("transfac_matrix.txt")
    m = Motif.read_transfac(f)
    f.close()
    assert m[3, "A"] == 0.0
    assert m[0, "G"] == 2.0
    assert np.shape(m.array) == (12, 4)
    f.close()

    f = data_stream("transfac_matrix2.txt")
    m = Motif.read_transfac(f)
    f.close()
    assert m[3, "A"] == 3.0
    assert m[0, "G"] == 152.0
    assert np.shape(m.array) == (15, 4)

    # this one has extra Ps on start of each line
    f = data_stream("transfac_matrix3.txt")
    m = Motif.read_transfac(f)
    f.close()


def test_motif_reindex() -> None:
    f = data_stream("transfac_matrix.txt")
    m = Motif.read_transfac(f)
    f.close()
    m2 = m.reindex("TCGA")

    assert str(m2.alphabet) == "TCGA"

    for k in range(0, 12):
        for i, a in enumerate("AGCT"):
            assert m[k, a] == m2[k, a]


def test_motif_read_transfac_nonstandard_alphabet() -> None:
    """Header with a digit column and non-standard alphabet.

    Covers branch 374->368 (isint(h) True, skip setting position_header)
    and branch 414->418 (no standard alphabet matches, loop exhausts).
    """
    from io import StringIO

    data = "PO\t1\tX\n01\t10\t20\n02\t30\t40\n03\t50\t60\nXX\n"
    m = Motif.read_transfac(StringIO(data))
    assert str(m.alphabet) == "1X"
    assert m.array.shape == (3, 2)
    assert m[0, "1"] == 10.0
    assert m[2, "X"] == 60.0


def test_motif_reverse() -> None:
    f = data_stream("transfac_matrix.txt")
    m = Motif.read_transfac(f)
    f2 = data_stream("transfac_matrix.txt")
    m2 = Motif.read_transfac(f2)
    m2.reverse()

    (K, N) = np.shape(m2)
    for k in range(0, K):
        for n in range(0, N):
            assert m[k, n] == m2[K - k - 1, n]

    f.close()
    f2.close()


def test_motif_complement() -> None:
    f = data_stream("transfac_matrix.txt")
    m = Motif.read_transfac(f)
    f2 = data_stream("transfac_matrix.txt")
    m2 = Motif.read_transfac(f2)
    m2.complement()

    (K, N) = np.shape(m2)
    for k in range(0, K):
        assert m[k, "A"] == m2[k, "T"]
        assert m[k, "G"] == m2[k, "C"]
        assert m[k, "C"] == m2[k, "G"]
        assert m[k, "T"] == m2[k, "A"]
    f.close()
    f2.close()


def test_motif_reverse_complement() -> None:
    f = data_stream("transfac_matrix.txt")
    m = Motif.read_transfac(f)

    f2 = data_stream("transfac_matrix.txt")
    m2 = Motif.read_transfac(f2)

    m.complement()
    m.reverse()

    m2.reverse_complement()

    assert (m.array == m2.array).all()
    f.close()
    f2.close()
