#!/usr/bin/env python

import unittest

def suite():
    modules_to_test = (
        'test_corebio.test_array_io',
        'test_corebio.test_astral',
        'test_corebio.test_clustal_io',
        'test_corebio.test_data',
        'test_corebio.test_db',
        'test_corebio.test_dssp',
        'test_corebio.test_fasta_io',
        'test_corebio.test_genbank_io',
        'test_corebio.test_intelligenetics_io',
        'test_corebio.test_matrix',
        'test_corebio.test_moremath',
        'test_corebio.test_msf_io',
        'test_corebio.test_nbrf_io',
        'test_corebio.test_nexus',
        'test_corebio.test_nexus_io',
        'test_corebio.test_null_io',
        'test_corebio.test_phylip_io',
        'test_corebio.test_plain_io',
        'test_corebio.test_ssearch_io',        
        'test_corebio.test_scop',
        'test_corebio.test_secstruc',
        'test_corebio.test_seq',
        'test_corebio.test_seq_io',
        'test_corebio.test_stockholm_io',
        'test_corebio.test_stride',
        'test_corebio.test_table_io',
        'test_corebio.test_transform',
        'test_corebio.test_utils',
        
    ) 

    alltests = unittest.TestSuite()
    for module in modules_to_test : 
        alltests.addTest(unittest.defaultTestLoader.loadTestsFromName(module))
    return alltests

if __name__ == '__main__':
    unittest.main(defaultTest='suite')