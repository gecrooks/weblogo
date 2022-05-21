#!/usr/bin/env python

import unittest

import numpy as np

from weblogo.matrix import AlphabeticArray, Motif
from weblogo.seq import protein_alphabet

from . import data_stream


class test_AlphabeticArray(unittest.TestCase):
    def test_create(self) -> None:
        matrix = AlphabeticArray((protein_alphabet, protein_alphabet))
        matrix["A", "C"] = 10
        assert matrix[0, 1] == 10.0


class test_Motif(unittest.TestCase):
    def test_read_transfac_alphabet_superset(self) -> None:
        with data_stream("transfac_matrix.txt") as f:
            Motif.read_transfac(f, alphabet="TCGA")

        # Supplied alphabet can be superset of defacto alphabet.
        # Reverts to defacto alphabet
        with data_stream("transfac_matrix.txt") as f:
            Motif.read_transfac(f, alphabet="TCGAXYZ")

    def test_read_transfac(self) -> None:
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

    def test_reindex(self) -> None:
        f = data_stream("transfac_matrix.txt")
        m = Motif.read_transfac(f)
        f.close()
        m2 = m.reindex("TCGA")

        assert str(m2.alphabet) == "TCGA"

        for k in range(0, 12):
            for i, a in enumerate("AGCT"):
                assert m[k, a] == m2[k, a]

    def test_reverse(self) -> None:
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

    def test_complement(self) -> None:
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

    def test_reverse_complement(self) -> None:
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


if __name__ == "__main__":
    unittest.main()
