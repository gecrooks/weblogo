#!/usr/bin/env python
 
#  Copyright (c) 2006, The Regents of the University of California, through 
#  Lawrence Berkeley National Laboratory (subject to receipt of any required
#  approvals from the U.S. Dept. of Energy).  All rights reserved.

#  This software is distributed under the new BSD Open Source License.
#  <http://www.opensource.org/licenses/bsd-license.html>
#
#  Redistribution and use in source and binary forms, with or without 
#  modification, are permitted provided that the following conditions are met: 
#
#  (1) Redistributions of source code must retain the above copyright notice, 
#  this list of conditions and the following disclaimer. 
#
#  (2) Redistributions in binary form must reproduce the above copyright 
#  notice, this list of conditions and the following disclaimer in the 
#  documentation and or other materials provided with the distribution. 
#
#  (3) Neither the name of the University of California, Lawrence Berkeley 
#  National Laboratory, U.S. Dept. of Energy nor the names of its contributors 
#  may be used to endorse or promote products derived from this software 
#  without specific prior written permission. 
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
#  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
#  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
#  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
#  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
#  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
#  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
#  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
#  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
#  POSSIBILITY OF SUCH DAMAGE. 
 

from corebio.utils import *
import unittest
from StringIO import StringIO
import re




class test_utils(unittest.TestCase) :
    assertTrue  = unittest.TestCase.failUnless
    assertFalse = unittest.TestCase.failIf
    
    def test_isfloat(self) :
        self.failUnless( isfloat('0.5'))
        self.failUnless( isfloat( ' 0') )
        self.failUnless( isfloat( '+1000000000  ') ) 
        self.failUnless( isfloat( '2') )
        self.failUnless( isfloat( '0000.2323') )
        self.failUnless( isfloat( '0.1e-23') )
        self.failUnless( isfloat( ' -0.5e+23') )

        self.assertFalse( isfloat( None) ) 
        self.assertFalse( isfloat( '') )
        self.assertFalse( isfloat( 'asdad') )
        self.assertFalse( isfloat( 'q34sd') )
        self.assertFalse( isfloat( '92384.kjdfghiksw')) 
        self.assertFalse( isfloat( 'adf!@#nn') )
                
    def test_isint(self) :
        self.failUnless( isint('0'))
        self.failUnless( isint('-1')) 
        self.failUnless( isint('10'))
        self.failUnless( isint('100101012234')) 
        self.failUnless( isint('000'))      
        
        self.assertFalse( isint( None) ) 
        self.assertFalse( isint( '') )
        self.assertFalse( isint( 'asdad') )
        self.assertFalse( isint( 'q34sd') )
        self.assertFalse( isint( '0.23')) 
        self.assertFalse( isint( 'adf!@#nn') )
               
    def test_remove_whitespace(self) :
        self.assertEquals( remove_whitespace("  kjashd askjdh askjdh\tasdf"),
        "kjashdaskjdhaskjdhasdf"),
   
    def test_isblank(self) :
        blank = ('', ' ', '\n', '\t \n\n')
        not_blank = (' a',)
        
        for s in blank :
            self.assertTrue( isblank(s))
        for s in not_blank:
            self.assertFalse( isblank(s))
 
    def test_group_count(self):
        test = 'aaabbbbcccddea'
        out = group_count(test)     
        self.assertTrue( tuple(out) == (('a',3),('b',4),('c',3),('d',2),('e',1),('a',1)) )    



        
    def test_reiterate(self) :
        i = Reiterate( iter("123456") )
        for item in i :
            pass
        self.assertRaises(StopIteration, i.next )
        self.assertFalse( i.has_item() )
        self.assertTrue( i.peek() is None )
        
        # pushback
        i = Reiterate( iter("123456") )
        i.next()
        i.push("0")
        self.assertEquals( "0", i.next())
        p = i.peek()
        n = i.next()
        self.assertEquals(p,n)
        self.assertEquals( i.index() ,2)
        self.assertTrue( i.has_item() )
        
        # Repeated application of Reiterate should return same iterator.
        assert i is iter(i)
        assert i is Reiterate(i)
        
        
        
    def test_token(self) :
        t = Token( 'kind', 'some data', 4, 3)
        s = str(t)
        r = repr(t)
        t2 = eval(r)
        assert t2.typeof == 'kind'

        
    def test_struct(self) :
        s = Struct(a=3,b=4)
        s2 = eval(repr(s))
        assert s2.a == 3
    
    def test_invert_dict(self) :
        d = dict( a=3, b=4)
        invd =invert_dict(d)
        assert 3 in invd
        assert invd[3] == 'a'
    
    def test_crc64(self) :     
        assert crc64("IHATEMATH") == "E3DCADD69B01ADD1"

    def test_crc32(self) :
        assert crc32("Test the CRC-32 of this string.") == "%08X"%1571220330

    def test_find_command(self) :
        p = find_command('more')
        p = find_command('python')
        
        self.assertRaises( EnvironmentError, find_command, 'NOSUCH')
        #print p

    def test_ArgumentValueError(self):
        message = "Some message"
        component = "whatsit"

        try:
            raise ArgumentError(message, component)
        except ArgumentError, err:
            assert err.msg == message
            assert err.key == component
            
        try:
            raise ArgumentError(message, component, 10)
        except ValueError, err:
            assert err.msg == message
            assert err.key == component
            assert err.value ==10
  
    def test_frozendict(self) :
        d = frozendict( a='b', c='c')
        assert d['a'] == 'b'

        self.assertRaises(AttributeError, lambda D: D.pop(), d)

    def test_file_index(self):
        stream = StringIO(tfile)
        
        idx = FileIndex(stream)
        assert idx[0].startswith('line 0')
        assert idx[4].startswith('line 4')

        def parser(f) :
            return int(f.readline().split()[1])

        idx = FileIndex(stream, parser=parser)
        assert len(idx) == 5
        assert idx[0] == 0
        assert idx[4] == 4


        key = re.compile(r"(line \d)")
        def linekey(line) :
            k = key.search(line)
            if k is None: return None
            return k.group(1)
            
        idx = FileIndex(stream, linekey=linekey, parser=parser)
        self.assertEquals(len(idx), 4)
        assert idx[0] == 0
        assert idx[3] == 4
        self.assertRaises(IndexError, idx.__getitem__, 5)  
        
        #print idx._key_dict
        assert idx['line 1'] == 1
        assert idx['line 4'] == 4  
        assert 'line 1' in idx
        assert not ('Blah' in idx)
        assert not (20 in idx)

        # Test iteration over values
        t = 0
        for v in idx : t+=v
        assert t == 8



        def test_resource(self) :
            fn = resource_filename(__name__, 'data/cap.fa', __file__)
            assert fn.endswith('test_corebio/data/cap.fa')
            f = resource_stream(__name__, 'data/cap.fa', __file__)
            s = resource_string(__name__, 'data/cap.fa', __file__)
            assert s.startswith('>aldB')



tfile = """line 0   
line 1
Blah
line 3
line 4 
"""




if __name__ == '__main__':
    unittest.main()

