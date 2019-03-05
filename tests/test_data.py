#!/usr/bin/env python

import unittest

from weblogo.data import (amino_acid_composition, resource_names, data_stream, data_string,
                          data_filename)


class test_data(unittest.TestCase):
    def test_amino_acid_composition(self):
        cl = [amino_acid_composition[k] for k in "ARNDCQEGHILKMFPSTWYV"]
        self.assertAlmostEqual(sum(cl), 1)

    def test_resources(self):
        for n in resource_names:
            data_string(n)
            f = data_stream(n)
            f.close()
            data_filename(n)


if __name__ == '__main__':
    unittest.main()
