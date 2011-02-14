#!/usr/bin/env python



from corebio import *
from corebio.seq import *
from corebio.seq_io import *
from corebio.seq import *
from corebio.secstruc.dssp import *
from test_corebio import *



import unittest

class test_dssp_io(unittest.TestCase) :

    def test_1(self) :
        r = DsspRecord(testdata_stream('1crn.dssp'))
        self.assertEquals(r.pdbid, "1crn")
        self.assertEquals(len(r.residues), 46 )

        res = r.residues[11]
        self.assertEquals(res.num, 12)
        self.assertEquals(res.resid, '12')
        self.assertEquals(res.chainid, ' ')
        self.assertEquals(res.aa, 'N')
        self.assertEquals(res.secstruc, "H")
        self.assertEquals(res.solvent_acc_area, float(82))
        self.assertEquals(res.phi, float(-64.9))
        self.assertEquals(res.psi, float(-39.5))  
        self.assertEquals(res.coord, (float( 3.5 ), float(3.9),float(6.8) )  )  


        self.assertEquals(str(r.primary()),
                            "TTabPSIVARSNFNVcRLPGTPEAIcATYTGbIIIPGATaPGDYAN")
        self.assertEquals(str(r.secondary()), 
                            " EE SSHHHHHHHHHHHTTT  HHHHHHHHS EE SSS   GGG  ")
        self.assertEquals(r.total_area(), float(3010.0))
        

  
             
if __name__ == '__main__':
    
    print "Running additional tests of RunDSSP. These require that the DSSP program is installed locally"
    try :
        dssp = RunDssp()
        print dssp.version()
        fin = testdata_stream('1CGP.pdb')
        data = dssp.process_pdb(fin)
        #print data
        fin = testdata_stream('1CGP.pdb')
        record = dssp.record(fin)
        #print record
    except e:
        print e
    
    # Now run standard unittests
    unittest.main()
    

