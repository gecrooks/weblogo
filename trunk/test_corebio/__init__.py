


import unittest
from corebio.utils import resource_string, resource_stream, resource_filename

def testdata_string( name ): 
    return resource_string(__name__, 'data/'+name, __file__)    

def testdata_stream( name ): 
    return resource_stream(__name__, 'data/'+name, __file__)    

def testdata_filename( name ): 
    return resource_filename(__name__, 'data/'+name, __file__)            

# 2.3 Python compatibility
assertTrue  = unittest.TestCase.failUnless
assertFalse = unittest.TestCase.failIf
