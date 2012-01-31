#!/usr/bin/env python
 


from StringIO import StringIO

from corebio import *
from corebio.seq import *
from corebio.utils import *
from corebio import data

import unittest
from test_corebio import *
from numpy import *

from corebio.matrix import *

class test_AlphabeticArray(unittest.TestCase) :
    def test_create(self): 

        matrix = AlphabeticArray( (protein_alphabet, protein_alphabet) )
        matrix['A', 'C'] = 10 
        assert matrix[0,1]==10.0
        
 
class test_Motif(unittest.TestCase) :
    def test_read_transfac(self):
        f = testdata_stream("transfac_matrix.txt")
        m = Motif.read_transfac(f)

        assert m[3, 'A']==0.0
        assert m[0, 'G']==2.0         
        assert shape(m.array) == (12,4)

        f = testdata_stream("transfac_matrix2.txt")
        m = Motif.read_transfac(f)

        assert m[3, 'A']==3.0
        assert m[0, 'G']==152.0         
        assert shape(m.array) == (15,4)



class test_SubMatrix(unittest.TestCase) :

    def test_create(self):
        ab = 'ABCD'
        ar = asarray( [[1,2,3,4], [5,6,7,8],[9,10,11,12],[13,14,15,16]])
        s = SubMatrix(ab, ar)
        
        assert s[0,0] == 1
        assert s['A', 'A'] ==1
        assert s['B', 'C'] ==7
        s['B', 'C'] =-1
        assert s['B', 'C'] ==-1

    def test_get(self):
        ab = Alphabet('ABCD')
        ar = asarray( [[1,2,3,4], [5,6,7,8],[9,10,11,12],[13,14,15,16]])
        s = SubMatrix(ab, ar)
        s1 = 'DCCBBBAAA'
        s2 = 'BA'
        v = s.index( (s1,s2) )
        #print v
        for m,i in enumerate(s1) :
            for n,j in enumerate(s2) :
                assert s[i,j]== v[m,n]

    def test_get_subMatrix(self):
        ab = Alphabet('ABCD')
        ar = asarray( [[1,2,3,4], [5,6,7,8],[9,10,11,12],[13,14,15,16]])
        mat = SubMatrix(ab, ar)
        
        mat2 = mat.reindex( 'ABC')
        assert all( mat2.array == asarray([[1,2,3], [5,6,7],[9,10,11]]) )

        mat2 = mat.reindex( 'BA')
        assert all(mat2.array == asarray([[6,5],[2,1]]) )

        mat2 = mat.reindex( Alphabet('BA'))
        assert all(mat2.array == asarray([[6,5],[2,1]]) )

                
    def test_fail_get(self) :
        ab = Alphabet('ABCD')
        ar = asarray( [[1,2,3,4], [5,6,7,8],[9,10,11,12],[13,14,15,16]])
        s = SubMatrix(ab, ar)
        
        self.assertRaises( IndexError, s.__getitem__ , ('E','A') )
        self.assertRaises( IndexError, s.__getitem__ , ('5','6') )
        
        #FIXME
        self.assertRaises( IndexError, s.index , ('E','A') )
 
    def test_repr(self) :
        ab = Alphabet('ABCD')
        ar = asarray( [[1,2,3,4], [5,6,7,8],[9,10,11,12],[13,14,15,16]])
        s = SubMatrix(ab, ar)
        
        string = repr(s)
        #print string        


    def test_read(self) :
        f = StringIO(test_matrix1)
        mat = SubMatrix.read(f)
        assert mat['a','a'] ==4

    def test_read_asymmetric_fail(self) :
        f = StringIO(test_matrix4)
        self.failUnlessRaises(ValueError, SubMatrix.read, f )

    def test_read_alphabets(self) :

        # incompatable alphabets
        f = StringIO(test_matrix3)
        self.failUnlessRaises(ValueError, 
            SubMatrix.read, f )

        f = StringIO(test_matrix3)
        mat = SubMatrix.read(f, alphabet = Alphabet('ARNDCQEGHILKMFPSTWYV'))
        
        f2 = StringIO(test_matrix1)
        self.failUnlessRaises(ValueError, 
            SubMatrix.read, f2 , unambiguous_protein_alphabet)

    def test_read_corrupt(self) :
        f = StringIO(test_matrix2)
        self.failUnlessRaises(ValueError, 
            SubMatrix.read, f )
 
    def test_read_pam(self):
        mat = SubMatrix.read( data.data_stream("pam250") )
        self.assertEquals(mat[0,0], 2.0)

        mat = SubMatrix.read( data.data_stream("pam120") )
        self.assertEquals(mat[4,5], -7)

    def test_read_blosum(self):
        mat = SubMatrix.read( data.data_stream("blosum80") )
        self.assertEquals(mat[0,10], -3)

        mat = SubMatrix.read( data.data_stream("blosum62") )
        self.assertEquals(mat[4,5], -4)

    def test_read_blast(self):
         # New style blast matrices have letters at beginning of lines and a '*'
        mat = SubMatrix.read(testdata_stream("blosum35.blast.new") )
        self.assertEquals(mat[4,5], -3)
        
        # Matrices formatted for old blast have a '*' (stop)
        # column and no letters at the beggining of lines 
        mat = SubMatrix.read(testdata_stream("blosum35.blast") )
        self.assertEquals(mat[0,10], -2)
        self.assertEquals(mat.array.shape, (23,23))

        # For comparison, we'll also parse a matrix without '*'
        mat = SubMatrix.read(testdata_stream("pam250.mat") )
        self.assertEquals(mat[4,5], -5)
        

    
             
test_matrix1 = """# A Test Matrix
# More comments

# And blank line should be ignored

A    4  -2  -2  -2   0  -1  -1  -1  -2  -2  -2  -2  -1  -2  -1   0  -1  -3  -2  -1  -2  -1  -1
R   -2   6  -1  -1  -4   1   0  -3   0  -3  -3   2  -2  -3  -2  -1  -1  -2  -2  -3   3  -1  -1
N   -2  -1   7   1  -3   0   0  -1   0  -5  -4   0  -3  -4  -2   0  -1  -3  -2  -4   3  -1  -1
D   -2  -1   1   7  -4   0   1  -1  -1  -6  -5   0  -4  -5  -1   0  -1  -4  -3  -5   0  -2  -2
C    0  -4  -3  -4  12  -3  -4  -3  -3  -1  -2  -4  -1  -2  -3  -2  -1  -2  -2   0  -3   5  -2
Q   -1   1   0   0  -3   6   1  -2   0  -3  -3   1  -2  -3  -1   0  -1  -2  -2  -3   0   1  -1
E   -1   0   0   1  -4   1   5  -2  -1  -4  -4   1  -3  -4  -1  -1  -1  -3  -3  -4   0  -1  -1
G   -1  -3  -1  -1  -3  -2  -2   7  -2  -6  -5  -2  -4  -5  -2  -1  -2  -4  -4  -5  -2  -2  -2
H   -2   0   0  -1  -3   0  -1  -2   9  -3  -3  -1  -2  -1  -2  -1  -1   0   0  -3   0  -1  -1
I   -2  -3  -5  -6  -1  -3  -4  -6  -3   5   2  -4   1   0  -4  -4  -2  -1  -1   3  -4  -2  -2
L   -2  -3  -4  -5  -2  -3  -4  -5  -3   2   5  -3   2   1  -3  -3  -2  -1  -1   1  -4  -2  -2
K   -2   2   0   0  -4   1   1  -2  -1  -4  -3   5  -2  -4  -1  -1  -1  -3  -3  -3   1  -1  -1
M   -1  -2  -3  -4  -1  -2  -3  -4  -2   1   2  -2   7   1  -3  -2  -1   0   0   1  -3  -2  -1
F   -2  -3  -4  -5  -2  -3  -4  -5  -1   0   1  -4   1   7  -3  -3  -2   3   3   0  -3  -2  -1
P   -1  -2  -2  -1  -3  -1  -1  -2  -2  -4  -3  -1  -3  -3   8  -1  -2  -3  -3  -3  -2  -2  -2
S    0  -1   0   0  -2   0  -1  -1  -1  -4  -3  -1  -2  -3  -1   4   1  -3  -2  -3   0  -1  -1
T   -1  -1  -1  -1  -1  -1  -1  -2  -1  -2  -2  -1  -1  -2  -2   1   5  -2  -2  -1  -1  -1  -1
W   -3  -2  -3  -4  -2  -2  -3  -4   0  -1  -1  -3   0   3  -3  -3  -2  12   3  -2  -3  -2  -1
Y   -2  -2  -2  -3  -2  -2  -3  -4   0  -1  -1  -3   0   3  -3  -2  -2   3   8  -2  -2  -2  -1
V   -1  -3  -4  -5   0  -3  -4  -5  -3   3   1  -3   1   0  -3  -3  -1  -2  -2   5  -4  -2  -2
B   -2   3   3   0  -3   0   0  -2   0  -4  -4   1  -3  -3  -2   0  -1  -3  -2  -4   3  -1  -1
Z   -1  -1  -1  -2   5   1  -1  -2  -1  -2  -2  -1  -2  -2  -2  -1  -1  -2  -2  -2  -1   3  -1
X   -1  -1  -1  -2  -2  -1  -1  -2  -1  -2  -2  -1  -1  -1  -2  -1  -1  -1  -1  -2  -1  -1  -1
"""

test_matrix2 = """# An invalid Test Matrix
# Its got a non-numerical value in it. Is the correct exception raised?

# And blank line should be ignored

A    4  -2  -2  -2   0  -1  -1  -1  -2  -2  -2  -2  -1  -2  -1   0  -1  -3  -2  -1  -2  -1  -1
R   -2   6  -1  -1  -4   1   0  -3   0  -3  -3   2  -2  -3  -2  -1  -1  -2  -2  -3   3  -1  -1
N   -2  -1   7   1  -3   0   0  -1   0  -5  -4   0  -3  -4  -2   0  -1  -3  -2  -4   3  -1  -1
D   -2  -1   1   7  -4   0   1  -1  -1  -6  -5   0  -4  -5  -1   0  -1  -4  -3  -5   0  -2  -2
C    0  -4  -3  -4  12  -3  -4  -3  -3  -1  -2  -4  -1  -2  -3  -2  -1  -2  -2   0  -3   5  -2
Q   -1   1   0   0  -3   6   1  -2   0  -3  -3   1  -2  -3  -1   0  -1  -2  -2  -3   0   1  -1
E   -1   0   0   1  -4   x   5  -2  -1  -4  -4   1  -3  -4  -1  -1  -1  -3  -3  -4   0  -1  -1
G   -1  -3  -1  -1  -3  -2  -2   7  -2  -6  -5  -2  -4  -5  -2  -1  -2  -4  -4  -5  -2  -2  -2
H   -2   0   0  -1  -3   0  -1  -2   9  -3  -3  -1  -2  -1  -2  -1  -1   0   0  -3   0  -1  -1
I   -2  -3  -5  -6  -1  -3  -4  -6  -3   5   2  -4   1   0  -4  -4  -2  -1  -1   3  -4  -2  -2
L   -2  -3  -4  -5  -2  -3  -4  -5  -3   2   5  -3   2   1  -3  -3  -2  -1  -1   1  -4  -2  -2
K   -2   2   0   0  -4   1   1  -2  -1  -4  -3   5  -2  -4  -1  -1  -1  -3  -3  -3   1  -1  -1
M   -1  -2  -3  -4  -1  -2  -3  -4  -2   1   2  -2   7   1  -3  -2  -1   0   0   1  -3  -2  -1
F   -2  -3  -4  -5  -2  -3  -4  -5  -1   0   1  -4   1   7  -3  -3  -2   3   3   0  -3  -2  -1
P   -1  -2  -2  -1  -3  -1  -1  -2  -2  -4  -3  -1  -3  -3   8  -1  -2  -3  -3  -3  -2  -2  -2
S    0  -1   0   0  -2   0  -1  -1  -1  -4  -3  -1  -2  -3  -1   4   1  -3  -2  -3   0  -1  -1
T   -1  -1  -1  -1  -1  -1  -1  -2  -1  -2  -2  -1  -1  -2  -2   1   5  -2  -2  -1  -1  -1  -1
W   -3  -2  -3  -4  -2  -2  -3  -4   0  -1  -1  -3   0   3  -3  -3  -2  12   3  -2  -3  -2  -1
Y   -2  -2  -2  -3  -2  -2  -3  -4   0  -1  -1  -3   0   3  -3  -2  -2   3   8  -2  -2  -2  -1
V   -1  -3  -4  -5   0  -3  -4  -5  -3   3   1  -3   1   0  -3  -3  -1  -2  -2   5  -4  -2  -2
B   -2   3   3   0  -3   0   0  -2   0  -4  -4   1  -3  -3  -2   0  -1  -3  -2  -4   3  -1  -1
Z   -1  -1  -1  -2   5   1  -1  -2  -1  -2  -2  -1  -2  -2  -2  -1  -1  -2  -2  -2  -1   3  -1
X   -1  -1  -1  -2  -2  -1  -1  -2  -1  -2  -2  -1  -1  -1  -2  -1  -1  -1  -1  -2  -1  -1  -1
"""

test_matrix3 = """#
# This test matrix has a smaller alphabet
A    4  -2  -2  -2   0  -1  -1  -1  -2  -2  -2  -2  -1  -2  -1   0  -1  -3  -2  -1  
R   -2   6  -1  -1  -4   1   0  -3   0  -3  -3   2  -2  -3  -2  -1  -1  -2  -2  -3  
N   -2  -1   7   1  -3   0   0  -1   0  -5  -4   0  -3  -4  -2   0  -1  -3  -2  -4  
D   -2  -1   1   7  -4   0   1  -1  -1  -6  -5   0  -4  -5  -1   0  -1  -4  -3  -5  
C    0  -4  -3  -4  12  -3  -4  -3  -3  -1  -2  -4  -1  -2  -3  -2  -1  -2  -2   0  
Q   -1   1   0   0  -3   6   4  -2   0  -3  -3   1  -2  -3  -1   0  -1  -2  -2  -3  
E   -1   0   0   1  -4   4   5  -2  -1  -4  -4   1  -3  -4  -1  -1  -1  -3  -3  -4  
G   -1  -3  -1  -1  -3  -2  -2   7  -2  -6  -5  -2  -4  -5  -2  -1  -2  -4  -4  -5  
H   -2   0   0  -1  -3   0  -1  -2   9  -3  -3  -1  -2  -1  -2  -1  -1   0   0  -3  
I   -2  -3  -5  -6  -1  -3  -4  -6  -3   5   2  -4   1   0  -4  -4  -2  -1  -1   3  
L   -2  -3  -4  -5  -2  -3  -4  -5  -3   2   5  -3   2   1  -3  -3  -2  -1  -1   1  
K   -2   2   0   0  -4   1   1  -2  -1  -4  -3   5  -2  -4  -1  -1  -1  -3  -3  -3  
M   -1  -2  -3  -4  -1  -2  -3  -4  -2   1   2  -2   7   1  -3  -2  -1   0   0   1  
F   -2  -3  -4  -5  -2  -3  -4  -5  -1   0   1  -4   1   7  -3  -3  -2   3   3   0  
P   -1  -2  -2  -1  -3  -1  -1  -2  -2  -4  -3  -1  -3  -3   8  -1  -2  -3  -3  -3  
S    0  -1   0   0  -2   0  -1  -1  -1  -4  -3  -1  -2  -3  -1   4   1  -3  -2  -3  
T   -1  -1  -1  -1  -1  -1  -1  -2  -1  -2  -2  -1  -1  -2  -2   1   5  -2  -2  -1  
W   -3  -2  -3  -4  -2  -2  -3  -4   0  -1  -1  -3   0   3  -3  -3  -2  12   3  -2  
Y   -2  -2  -2  -3  -2  -2  -3  -4   0  -1  -1  -3   0   3  -3  -2  -2   3   8  -2  
V   -1  -3  -4  -5   0  -3  -4  -5  -3   3   1  -3   1   0  -3  -3  -1  -2  -2   5  
"""

test_matrix4 = """# This matrix is invalid because it is asymetric! (AR, RA)

A    4   2  -2  -2   0  -1  -1  -1  -2  -2  -2  -2  -1  -2  -1   0  -1  -3  -2  -1  -2  -1  -1
R   -2   6  -1  -1  -4   1   0  -3   0  -3  -3   2  -2  -3  -2  -1  -1  -2  -2  -3   3  -1  -1
N   -2  -1   7   1  -3   0   0  -1   0  -5  -4   0  -3  -4  -2   0  -1  -3  -2  -4   3  -1  -1
D   -2  -1   1   7  -4   0   1  -1  -1  -6  -5   0  -4  -5  -1   0  -1  -4  -3  -5   0  -2  -2
C    0  -4  -3  -4  12  -3  -4  -3  -3  -1  -2  -4  -1  -2  -3  -2  -1  -2  -2   0  -3   5  -2
Q   -1   1   0   0  -3   6   1  -2   0  -3  -3   1  -2  -3  -1   0  -1  -2  -2  -3   0   1  -1
E   -1   0   0   1  -4   1   5  -2  -1  -4  -4   1  -3  -4  -1  -1  -1  -3  -3  -4   0  -1  -1
G   -1  -3  -1  -1  -3  -2  -2   7  -2  -6  -5  -2  -4  -5  -2  -1  -2  -4  -4  -5  -2  -2  -2
H   -2   0   0  -1  -3   0  -1  -2   9  -3  -3  -1  -2  -1  -2  -1  -1   0   0  -3   0  -1  -1
I   -2  -3  -5  -6  -1  -3  -4  -6  -3   5   2  -4   1   0  -4  -4  -2  -1  -1   3  -4  -2  -2
L   -2  -3  -4  -5  -2  -3  -4  -5  -3   2   5  -3   2   1  -3  -3  -2  -1  -1   1  -4  -2  -2
K   -2   2   0   0  -4   1   1  -2  -1  -4  -3   5  -2  -4  -1  -1  -1  -3  -3  -3   1  -1  -1
M   -1  -2  -3  -4  -1  -2  -3  -4  -2   1   2  -2   7   1  -3  -2  -1   0   0   1  -3  -2  -1
F   -2  -3  -4  -5  -2  -3  -4  -5  -1   0   1  -4   1   7  -3  -3  -2   3   3   0  -3  -2  -1
P   -1  -2  -2  -1  -3  -1  -1  -2  -2  -4  -3  -1  -3  -3   8  -1  -2  -3  -3  -3  -2  -2  -2
S    0  -1   0   0  -2   0  -1  -1  -1  -4  -3  -1  -2  -3  -1   4   1  -3  -2  -3   0  -1  -1
T   -1  -1  -1  -1  -1  -1  -1  -2  -1  -2  -2  -1  -1  -2  -2   1   5  -2  -2  -1  -1  -1  -1
W   -3  -2  -3  -4  -2  -2  -3  -4   0  -1  -1  -3   0   3  -3  -3  -2  12   3  -2  -3  -2  -1
Y   -2  -2  -2  -3  -2  -2  -3  -4   0  -1  -1  -3   0   3  -3  -2  -2   3   8  -2  -2  -2  -1
V   -1  -3  -4  -5   0  -3  -4  -5  -3   3   1  -3   1   0  -3  -3  -1  -2  -2   5  -4  -2  -2
B   -2   3   3   0  -3   0   0  -2   0  -4  -4   1  -3  -3  -2   0  -1  -3  -2  -4   3  -1  -1
Z   -1  -1  -1  -2   5   1  -1  -2  -1  -2  -2  -1  -2  -2  -2  -1  -1  -2  -2  -2  -1   3  -1
X   -1  -1  -1  -2  -2  -1  -1  -2  -1  -2  -2  -1  -1  -1  -2  -1  -1  -1  -1  -2  -1  -1  -1
"""
             
             
if __name__ == '__main__':
    unittest.main()

