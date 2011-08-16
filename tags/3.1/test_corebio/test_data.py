#!/usr/bin/env python

import unittest
import time

from corebio.data import *
from corebio.seq import *

class test_data(unittest.TestCase) :

    def test_amino_acid_composition(self) :
        cl = [ amino_acid_composition[k] for k in "ARNDCQEGHILKMFPSTWYV"]
        self.assertAlmostEquals( sum(cl), 1)
    
   
    def test_resources(self) :
        for n in resource_names :
            s = data_string(n)
            f = data_stream(n)
            fn = data_filename(n)
                    
if __name__ == '__main__':
    unittest.main()
