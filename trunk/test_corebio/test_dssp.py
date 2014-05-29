#!/usr/bin/env python

from __future__ import print_function

import unittest

from corebio import *
from corebio.seq import *
from corebio.seq_io import *
from corebio.secstruc.dssp import *
from test_corebio import *


class test_dssp_io(unittest.TestCase) :

    def test_1(self) :
        f = testdata_stream('1crn.dssp')
        r = DsspRecord(f)
        self.assertEqual(r.pdbid, "1crn")
        self.assertEqual(len(r.residues), 46 )

        res = r.residues[11]
        self.assertEqual(res.num, 12)
        self.assertEqual(res.resid, '12')
        self.assertEqual(res.chainid, ' ')
        self.assertEqual(res.aa, 'N')
        self.assertEqual(res.secstruc, "H")
        self.assertEqual(res.solvent_acc_area, float(82))
        self.assertEqual(res.phi, float(-64.9))
        self.assertEqual(res.psi, float(-39.5))
        self.assertEqual(res.coord, (float(3.5), float(3.9), float(6.8)))

        self.assertEqual(str(r.primary()),
                          "TTabPSIVARSNFNVcRLPGTPEAIcATYTGbIIIPGATaPGDYAN")
        self.assertEqual(str(r.secondary()),
                          " EE SSHHHHHHHHHHHTTT  HHHHHHHHS EE SSS   GGG  ")
        self.assertEqual(r.total_area(), float(3010.0))

        f.close()

if __name__ == '__main__':
    print("Running additional tests of RunDSSP. These require that the DSSP program is installed locally")
    try:
        dssp = RunDssp()
        print(dssp.version())
        fin = testdata_stream('1CGP.pdb')
        data = dssp.process_pdb(fin)
        #print data
        fin = testdata_stream('1CGP.pdb')
        record = dssp.record(fin)
        #print record
    except Exception as exc:
        print(exc)
    # Now run standard unittests
    unittest.main()
