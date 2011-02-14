#!/usr/bin/env python



import unittest

from corebio.db  import *
from test_corebio import *

class DBTest(unittest.TestCase):
   
   
    def test_dbxref(self):
        ref = Dbxref('/db_xref="GDB:39999"')
        assert ref.database == 'GDB'
        assert ref.identifier == '39999'

        ref = Dbxref("GDB:39999")
        assert ref.database == 'GDB'
        assert ref.identifier == '39999'
        
        ref = Dbxref('GDB','39999')
        assert ref.database == 'GDB'
        assert ref.identifier == '39999'
        
        s = repr(ref)
        s2 = str(ref)
       
    def test_datasource(self) :
        ds = DataSource(
            abbrev='abbrev', 
            alt_abbrev=('a','b') ,
            url = 'abc',
            resource_url = 'abc%s',
            parser = None,
            description = 'blah',
            name = 'somename'
            )
            
    def test_databases(self):
        db = default_registry['pdb']

        
if __name__ == '__main__':
    
    print "## Know databases."
    print
    print default_registry
    print
    print
    
    print "## Running Non-unit tests."
    print
    
    tests = [
        ('embl', 'ab050095'),
        ('pdb', '1hlb'),
        ('swissprot', 'p50105'),        
    ]
    
    for t in tests:
        ref = Dbxref(t[0], t[1])
        print ref, ref.data_url()

        data = ref.data_stream()
        print data.readline(),
        print data.readline(),
        print data.readline(),
        print data.readline(),
        print

    print
    print
    print '# Running Unittests'
    unittest.main()



