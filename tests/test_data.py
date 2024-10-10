#!/usr/bin/env python

import unittest

from weblogo.data import amino_acid_composition


class test_data(unittest.TestCase):
    def test_amino_acid_composition(self) -> None:
        cl = [amino_acid_composition[k] for k in "ARNDCQEGHILKMFPSTWYV"]
        self.assertAlmostEqual(sum(cl), 1)


if __name__ == "__main__":
    unittest.main()
