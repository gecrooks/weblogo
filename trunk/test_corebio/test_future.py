#!/usr/bin/env python




import unittest
from corebio import _future
from corebio._future import resource_filename
from corebio._future import resource_stream
from corebio._future import resource_string
from corebio._future import Template
from corebio._future import subprocess





class test_future(unittest.TestCase) :
    assertTrue  = unittest.TestCase.failUnless
    assertFalse = unittest.TestCase.failIf
 
 
 
    def test_resource(self) :
        fn = resource_filename(__name__, 'data/cap.fa', __file__)
        assert fn.endswith('test_corebio/data/cap.fa')
        f = resource_stream(__name__, 'data/cap.fa', __file__)
        s = resource_string(__name__, 'data/cap.fa', __file__)
        assert s.startswith('>aldB')

 
 
 
if __name__ == '__main__':
    unittest.main()

