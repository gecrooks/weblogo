
#  Copyright (c) 2005 Gavin E. Crooks
#  Copyright (c) 2006 John Gilman

#  This software is distributed under the MIT Open Source License.
#  <http://www.opensource.org/licenses/mit-license.html>
#
#  Permission is hereby granted, free of charge, to any person obtaining a 
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
#  IN THE SOFTWARE.


""" 
Arrays indexed by alphabetic strings.
"""
from __future__ import absolute_import, print_function


import numpy as na
from corebio._py3k import zip

from .seq import Alphabet
from .seq import (unambiguous_dna_alphabet, unambiguous_rna_alphabet,
                  unambiguous_protein_alphabet)

from .utils import isint

__all__= 'AlphabeticArray', 'submatrix_alphabet', 'SubMatrix', 'Motif'

class AlphabeticArray(object) :
    """An alphabetic array. Wraps a numpy array so that each dimension
    can be associated with an alphabet and indexed with characters or strings.
    
    Attributes :
    - alphabets -- A sequence of alphabets used to index the array
    - array     -- The underlying array object that is indexed.
    
    Examples : 
    
    >>> from corebio.seq import *
    >>> from corebio.matrix import AlphabeticArray
    >>>
    >>> str(protein_alphabet)
    'ACDEFGHIKLMNOPQRSTUVWYBJZX*-'
    >>> matrix = AlphabeticArray( (protein_alphabet, protein_alphabet) )
    >>>
    >>> # Index by character or integer:
    >>> matrix['A', 'C'] = 10 
    >>> matrix[0,1]
    10
    >>>
    >>> # Different alphabets on each dimension:
    >>> import numpy as na    
    >>> a234 = na.zeros( shape = (2,3,4) )
    >>> alpha = ( "AB", "ABC", "ABCD")
    >>> aa = AlphabeticArray(alpha,a234)
    >>> aa['A', 'B', 'C'] = 22
    >>>
    >>> # String indices are converted to integer index arrays:
    ...
    >>> aa['A', 'B', 'ABCD']
    array([ 0,  0, 22,  0])
    
    
    Authors: 
    o GEC 2005, JXG 2006
    
    """
    # Design note: Subclassing numpy arrays is hard, so instead we
    # build this proxy wrapper.
    
    __slots__ = ['alphabets', 'array']
    
    def __init__(self, alphabets, values=None, dtype=None) :
        """
        Args:
        - alphabets -- a list of alphabets (as string or Alphabet objects) to
                    be used to convert strings into indices. The lengths of 
                    the alphabets match the shape of the indexed array. 
                    Alternatively, an integer or None in the list indicate a 
                    non-alphabetic dimension. If None the dimension length is 
                    taken from values argument.
        - values -- An array of values to be indexed. If None a new  
                 array is created. If this argument is not a numpy array
                 then the alphabet list must be explicit (cannot contain 
                 None.)
        - dtype -- An optional numpy type code.
        """
    
        # A dummy object to be used in place of None in the alphabets list
        # so that we get meaningful error messages if we try to index a 
        # nonalphabetic dimension with a string.
        class NullAlphabet(object) :
            def ord(self, key) :
                raise IndexError('This dimension does not have an alphabet')
            def ords(self, key) :
                raise IndexError('This dimension does not have an alphabet')
        
        alpha = []
        shape = []
        for a in alphabets :
            if isinstance(a, str) : a =  Alphabet(a)
           
            if a is None :
                shape.append(None)
                alpha.append(NullAlphabet())
            elif isinstance(a, Alphabet) : 
                shape.append(len(a) )
                alpha.append(a)
            else :
                shape.append(int(a) )
                alpha.append(None)
        
        shape = tuple(shape)
        if values is None :
            values = na.zeros( shape=shape, dtype=dtype)
        else :
            values = na.asarray(values, dtype=dtype)
            vshape = values.shape
            if len(shape) != len(vshape) :
                raise ValueError("The values array is the wrong shape.")
            for s1, s2 in zip(shape, vshape):
                if s1 is not None and s1 != s2:
                    raise ValueError("The values array is the wrong shape.")
        self.array = values
        self.alphabets = tuple(alpha)
             
             
    def __getitem__(self, key) :
        return self.array.__getitem__( self._ordkey(key) )
        
    def __setitem__(self, key, value) :
        self.array.__setitem__( self._ordkey(key) , value ) 

    def _ordkey(self,  key) :
        """Convert string indices into integers. Handles characters, strings
        slices with strings, and tuples of the same. Anything else is 
        unchanged.
        """
        
        def norm(key, alpha) :
            if key is None : 
                return None
            elif isinstance(key, str) or isinstance(key, Alphabet) : 
                key = str(key)
                if len(key) ==1 : return alpha.ord(key)
                if len(key) ==0 : return None
                return na.asarray(alpha.ords(key))
            elif isinstance(key, slice) :
                start = norm(key.start, alpha)
                stop = norm(key.stop, alpha)
                step = key.step
                return slice(start, stop, step)
            else :
                return key            
            
        if isinstance(key, tuple) :
            return tuple([norm(k,a) for k, a in zip(key, self.alphabets)])
        else :
            return norm(key, self.alphabets[0])
 
    def index(self, keys) :
        """ Return an array of shape (len(key1), len(key2), ...) whose values       
        are indexed by keys.

        a.outerindex( (I,J,K) )[i,j,k] == a.array[I[i],J[j],K[k]]  

        """        
        # TODO: Above docstring is not very clear. 
        # Deep voodoo using numpy indexing     
        keys =  self._ordkey( keys)
        outerkeys = []
        for i, k in enumerate(keys) :
            if k is None: 
                k = range(0,self.array.shape[i])
            k = na.asarray(k)
            for j in range(len(keys)-i-1) : 
                k = k[...,na.newaxis]
            outerkeys.append(k)
        return self.array.__getitem__(  tuple(outerkeys) ) 
 
 
    def reindex(self, new_alphabets) :
        """Create a new AlphabeticArray with the given alphabets. The new 
        alphabet must be a subset of the current alphabet. Useful for
        extracting a submatrix or for permuting the alphabet.
        """
        new_array = self.index(new_alphabets)
        return AlphabeticArray(new_alphabets,new_array)

 
    # The following code is designed to proxy all attributes
    # of the wrapped array. But I'm not entirely sure that this will work as
    # intended.       
    def __getattr__(self, name) :
        try :
            return object.__getattr__(self, name) 
        except AttributeError:
            return getattr(self.array, name)

    def __setattr__(self, name, value) :
        try :
            return object.__setattr__(self, name, value) 
        except AttributeError:
            return setattr(self.array, name, value)

# End class AlphabeticArray



#TODO: move to seq?
submatrix_alphabet = Alphabet("ARNDCQEGHILKMFPSTWYVBZX")

class SubMatrix(AlphabeticArray) :
    """A two dimensional array indexed by an Alphabet. Used to hold substitution
    matrices and similar information. 
    
    Various standard substitution matrices are available from the data package
    >>> from corebio import data
    >>> mat = SubMatrix.read(data.data_stream('blosum100'))   
    
    Attr:
    - alphabet     -- An Alphabet
    - array        -- A numpy array
    - name         -- The name of this matrix (if any) as a string.
    - description  -- The description, if any.
    - scale        -- The scale constant of a log-odds matrix, if known.     
    
    Authors: 
    o GEC 2005, JXG 2006
    
    """
    # TODO: __str__
    # TODO: __repr__
    # TODO: normalize
    # TODO: freq->log_odds (With additional ambiguity characters?)
    # TODO: from_seqs   

 
    __slots__ = ['alphabet', 'array', 'name', 'description', 'scale' ] 
        
    def __init__(self, alphabet, array=None, typeof=None, name=None,
            description = None, scale=None) :
        AlphabeticArray.__init__(self, (alphabet, alphabet), array, typeof)
        self.alphabet = Alphabet(alphabet)
        self.name = name
        self.description = description
        self.scale = scale    
        
        
    def reindex(self, alphabet) :
        return  AlphabeticArray.reindex(self, (alphabet, alphabet))
        
        
        
    @staticmethod
    def read(fin, alphabet=None, typeof=na.float64) :    
        """ Parse and return a substitution matrix
        
        Arguments:
        - fin       --  matrix file 
        - alphabet  -- The set of substitution characters. Default: ''
        -  typeof    -- A numpy type or typecode.
        Returns:
        -  A numpy matrix of substitution scores
        Raises:
        -  ValueError on unreadable input
        """
        # TODO: Parse name, description, scale, where avaliable.
        # TODO: Include '*' in submatrix_alphabet
        # TODO: Read DNA substitution matrixes
        if alphabet is None : alphabet =submatrix_alphabet
        L = len(alphabet)
        matrix = na.zeros( (L,L), typeof)
        
        i = 0 
        
        #print(">", alphabet)
        
        for linenum,line in enumerate(fin) :
            #print(">>", linenum, i, line)
            if line.isspace() or line[0] =='#' or line[0]=='*': 
                continue
                
            cells = line.split()
            
            # Header line? "A  R  N  D  C  Q  E..."
            if cells[1] == alphabet[1] : continue
                
            # Lines look like this:
    # A  5 -1 -1 -1 -2  0 -1  0 -2 -1 -2  0  0 -2 -2  1  0 -2 -1  0 -1 -1  0 -5 
            # The initial character and final number (corresponds to '*' stop)
            # are optional.
            
            if cells[0].isalpha() and cells[0] != alphabet[i] :
                raise ValueError("Incompatible alphabet: line %d : %s %s: " 
                    % (linenum, line[0], alphabet[i]) )
    
            if cells[0].isalpha() : cells = cells[1:]
            if len(cells) == 24: cells = cells[:23]   # Chop off '*' state
            if len(cells) != L:
                    raise ValueError(
                        "SubMatrix matrix parse error: line %d"% linenum )
            
            for j in range(0, L) :
                matrix[i,j] = float(cells[j])
                # FIXME Should catch and rethrow parsing error here?
    
            i +=1
            if i == L : break
            
        if i != L :
                raise ValueError("Premature EOF")
    
        for i in range(0,L) :
            for j in range(0,L) :
                if matrix[i,j] != matrix[j,i] :
                    raise ValueError(
                        "Substitution matrix is asymmetric! (%d,%d)" %(i,j) )
            
        return SubMatrix(alphabet,matrix)    
  
# End class SubMatrix




#TODO
# Separate PWM (Position weight matrix. (Log odds?) 
#, ICM (Information content matrix
#, PFM (Position frequency matrix)


class Motif(AlphabeticArray) :
    """A two dimensional array where the second dimension is indexed by an 
    Alphabet. Used to represent sequence motifs and similar information.

    
    Attr:
    - alphabet     -- An Alphabet
    - array        -- A numpy array
    - name         -- The name of this motif (if any) as a string.
    - description  -- The description, if any.
    
    """
            
    def __init__(self, alphabet, array=None, dtype=None, name=None,
            description = None, scale=None) :
        AlphabeticArray.__init__(self, (None, alphabet), array, dtype)
        self.name = name
        self.description = description
        self.scale = scale    

    @property
    def alphabet(self):
        return self.alphabets[1]
        
    def reindex(self, alphabet) :
        return  Motif(alphabet, AlphabeticArray.reindex(self, (None, alphabet)))
  
    # These methods alter self, and therefore do not return a value.
    # (Compare to Seq objects, where the data is immutable and therefore methods return a new Seq.)
    # TODO: Should reindex (above) also act on self?
    
    def reverse(self):
        """Reverse sequence data"""
#        self.array = na.array(self.array[::-1]) # This is a view into the origional numpy array.
        self.array = self.array[::-1] # This is a view into the origional numpy array.
    
    def complement(self) :
        """Complement nucleic acid sequence."""
        from corebio.seq import Seq, Alphabet
        alphabet = self.alphabet
        complement_alphabet = Alphabet(Seq(alphabet, alphabet).complement())
        self.alphabets = (None, complement_alphabet) 
        
        m = self.reindex(alphabet)
        self.alphabets = (None, alphabet)
        self.array = m.array
        

    def reverse_complement(self):
         """Complements and reverses nucleic acid sequence (i.e. the other strand
         of a DNA sequence.) 
         """
         self.reverse()
         self.complement()
     
      
      
      
    @staticmethod #TODO: should be classmethod?
    def read_transfac( fin, alphabet = None) :
        """ Parse a sequence matrix from a file. 
        Returns a tuple of (alphabet, matrix)
        """
   
        items = []

        start=True
        for line in fin :
            if line.isspace() or line[0] =='#' : continue
            stuff = line.split()
            if start and stuff[0] != 'PO' and stuff[0] != 'P0': continue
            if stuff[0]=='XX' or stuff[0]=='//': break
            start = False
            items.append(stuff)
        if len(items) < 2  :
            raise ValueError("Vacuous file.")

        # Is the first line a header line?
        header = items.pop(0)
        hcols = len(header)
        rows = len(items)
        cols = len(items[0])
        if not( header[0] == 'PO' or header[0] =='P0' or hcols == cols-1 or hcols == cols-2) :
            raise ValueError("Missing header line!")

        # Do all lines (except the first) contain the same number of items?
        cols = len(items[0])
        for i in range(1, len(items)) :
            if cols != len(items[i]) :
                raise ValueError("Inconsistant length, row %d: " % i)

        # Vertical or horizontal arrangement?
        if header[0] == 'PO' or header[0] == 'P0': header.pop(0)

        position_header = True    
        alphabet_header = True    
        for h in header :
            if not isint(h) : position_header = False
            if not str.isalpha(h) : alphabet_header = False

        if not position_header and not alphabet_header :
            raise ValueError("Can't parse header: %s" % str(header))

        if position_header and alphabet_header :
            raise ValueError("Can't parse header")


        # Check row headers
        if alphabet_header :
            for i,r in enumerate(items) :
                if not isint(r[0]) and r[0][0]!='P' : 
                    raise ValueError(
                        "Expected position as first item on line %d" % i)
                r.pop(0)
                defacto_alphabet = ''.join(header)
        else :
            a = []
            for i,r in enumerate(items) :
                if not ischar(r[0]) and r[0][0]!='P' : 
                    raise ValueError(
                        "Expected position as first item on line %d" % i)
                a.append(r.pop(0))
            defacto_alphabet = ''.join(a)                

        # Check defacto_alphabet
        defacto_alphabet = Alphabet(defacto_alphabet)

        if alphabet :
            if not defacto_alphabet.alphabetic(alphabet) :
                raise ValueError("Incompatible alphabets: %s , %s (defacto)"
                                 % (alphabet, defacto_alphabet))
        else :            
            alphabets = (unambiguous_rna_alphabet,
                        unambiguous_dna_alphabet,                      
                        unambiguous_protein_alphabet,
                      )
            for a in alphabets :
                if defacto_alphabet.alphabetic(a) :
                    alphabet = a
                    break
            if not alphabet :
                alphabet = defacto_alphabet
   

        # The last item of each row may be extra cruft. Remove
        if len(items[0]) == len(header) +1 :
            for r in items :
                r.pop()

        # items should now be a list of lists of numbers (as strings) 
        rows = len(items)
        cols = len(items[0])
        matrix = na.zeros( (rows,cols) , dtype=na.float64) 
        for r in range( rows) :
            for c in range(cols):
                matrix[r,c] = float( items[r][c]) 

        if position_header :
            matrix.transpose() 

        return Motif(defacto_alphabet, matrix).reindex(alphabet)

        

   
        
