

from corebio import *
from corebio.seq import *
from corebio.seq_io import *
from corebio.seq import *
from corebio.secstruc import *
from test_corebio import *



import unittest

class test_secstruc(unittest.TestCase) :

    def test_1(self) :
        record = dssp.DsspRecord( testdata_stream('1crn.dssp') )
        reduced = fa_reduce_secstruc_to_ehl(record.secondary())
        assert str(reduced) == 'LEELLLHHHHHHHHHHHLLLLLHHHHHHHHLLEELLLLLLLLLLLL'
             
if __name__ == '__main__':
      
    unittest.main()
    
    