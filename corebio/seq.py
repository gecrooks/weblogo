
#  Copyright (c) 2005 Gavin E. Crooks <gec@compbio.berkeley.edu>
#
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
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
#  THE SOFTWARE.
#



""" Alphabetic sequences and associated tools and data.

Seq is a subclass of a python string with additional annotation and an alphabet.
The characters in string must be contained in the alphabet. Various standard
alphabets are provided.


Classes :
    Alphabet    -- A subset of non-null ascii characters
    Seq         -- An alphabetic string
    SeqList     -- A collection of Seq's
  
Alphabets :    
    o generic_alphabet  -- A generic alphabet. Any printable ASCII character.
    o protein_alphabet -- IUCAP/IUB Amino Acid one letter codes. 
    o nucleic_alphabet -- IUPAC/IUB Nucleic Acid codes 'ACGTURYSWKMBDHVN-'
    o dna_alphabet -- Same as nucleic_alphabet, with 'U' (Uracil) an 
        alternative for 'T' (Thymidine).
    o rna_alphabet -- Same as nucleic_alphabet, with 'T' (Thymidine) an
        alternative for 'U' (Uracil).
    o reduced_nucleic_alphabet -- All ambiguous codes in 'nucleic_alphabet' are
        alternative to 'N' (aNy)
    o reduced_protein_alphabet -- All ambiguous ('BZJ') and non-canonical amino 
        acids codes ( 'U', Selenocysteine and 'O', Pyrrolysine)  in 
        'protein_alphabet' are alternative to 'X'.
    o unambiguous_dna_alphabet -- 'ACGT'
    o unambiguous_rna_alphabet -- 'ACGU'
    o unambiguous_protein_alphabet -- The twenty canonical amino acid one letter
        codes, in alphabetic order, 'ACDEFGHIKLMNPQRSTVWY'

Amino Acid Codes:
    Code  Alt.  Meaning
    -----------------
    A           Alanine
    B           Aspartic acid or Asparagine
    C           Cysteine
    D           Aspartate
    E           Glutamate
    F           Phenylalanine
    G           Glycine
    H           Histidine
    I           Isoleucine
    J           Leucine or Isoleucine    
    K           Lysine
    L           Leucine
    M           Methionine
    N           Asparagine
    O           Pyrrolysine    
    P           Proline
    Q           Glutamine
    R           Arginine
    S           Serine
    T           Threonine
    U           Selenocysteine
    V           Valine
    W           Tryptophan
    Y           Tyrosine
    Z           Glutamate or Glutamine
    X    ?      any
    *           translation stop
    -    .~     gap 

Nucleotide Codes:
    Code  Alt.  Meaning
    ------------------------------
    A           Adenosine
    C           Cytidine
    G           Guanine
    T           Thymidine
    U           Uracil
    R           G A (puRine)
    Y           T C (pYrimidine)
    K           G T (Ketone)
    M           A C (aMino group)
    S           G C (Strong interaction)
    W           A T (Weak interaction)
    B           G T C (not A) (B comes after A)
    D           G A T (not C) (D comes after C)
    H           A C T (not G) (H comes after G)
    V           G C A (not T, not U) (V comes after U)
    N   X?      A G C T (aNy)
    -   .~      A gap 
    



Refs:
    http://www.chem.qmw.ac.uk/iupac/AminoAcid/A2021.html
    http://www.chem.qmw.ac.uk/iubmb/misc/naseq.html    
Status:
    Beta    
Authors:
    GEC 2004,2005
"""

# TODO: Add this to docstring somewhere.
# To replace all ambiguous nucleic code by 'N', replace alphabet and then n 
# normalize.
# 
# >>> Seq( 'ACGT-RYKM', reduced_nucleic_alphabet).normalized()
# 'ACGT-NNNN'
    
from array import array
from string import maketrans
from corebio.moremath import argmax, sqrt

__all__ = [
    'Alphabet', 
    'Seq', 
    'rna', 'dna', 'protein',
    'SeqList', 
    'generic_alphabet', 
    'protein_alphabet', 
    'nucleic_alphabet',
    'dna_alphabet',
    'rna_alphabet', 
    'reduced_nucleic_alphabet',
    'reduced_protein_alphabet',
    'unambiguous_dna_alphabet',
    'unambiguous_dna_alphabet', 
    'unambiguous_rna_alphabet', 
    'unambiguous_protein_alphabet',
    'generic_alphabet',
    ]



class Alphabet(object) :
    """An ordered subset of printable ascii characters.

    Status:
        Beta
    Authors: 
        - GEC 2005
    """
    __slots__ = ['_letters', '_alternatives','_ord_table', '_chr_table']
 
    # We're immutable, so use __new__ not __init__
    def __new__(cls, letters, alternatives= None) :
        """Create a new, immutable Alphabet.
        
        arguments:
        - letters -- the letters in the alphabet. The ordering determines
            the ordinal position of each character in this alphabet.
        - alt -- A list of (alternative, canonical) letters. The alternatives
            are given the same ordinal position as the canonical characters. 
            e.g. (('?','X'),('x', 'X')) states that '?' and 'x' are synonomous 
            with 'X'.  Values that are not in 'letters' are ignored. Alternatives
            that are already in 'letters' are also ignored. If the same
            alternative character is used twice then the alternative is assigned
            to the canonical character that occurs first in 'letters'. The 
            default is to assume that upper and lower case characters are
            equivalent, unless both cases are included in 'letters'.                   
        raises:
            ValueError : Repetitive or otherwise illegal set of letters.        
        """
        self = object.__new__(cls)

        # Printable Ascii characters
        ascii_letters = "".join([chr(__i) for __i in range(32,128)])

        if letters is None : letters = ascii_letters
        self._letters = letters

        equivalent_by_case = zip( 'abcdefghijklmnopqrstuvwxyz',
                                  'ABCDEFGHIJKLMNOPQRSTUVWXYZ')

        if alternatives is None : alternatives = equivalent_by_case

        
        # The ord_table maps between the ordinal position of a character in ascii
        # and the ordinal position in this alphabet. Characters not in the
        # alphabet are given a position of 255. The ord_table is stored as a 
        # string. 
        ord_table = ["\xff",] * 256
        for i,a in enumerate(letters) :
            n = ord(a)
            if n == 0 :
                raise ValueError("Alphabet cannot contain null character \\0")
            if ord_table[ n ] != "\xff":
                raise ValueError("Repetitive alphabet")
            ord_table[ n ] = chr(i)

        # Add alternatives
        _from = []
        _to = []
        for e, c in alternatives :
            if c in letters :
                n = ord(e)
                if ord_table[ n ] == "\xff" : # empty  
                    ord_table[ n ] = ord_table[ ord(c)]     
                    _from.append(e)
                    _to.append(c)
        self._alternatives = (''.join(_from), ''.join(_to))
                                
        ord_table = "".join(ord_table)
        assert( ord_table[0] == "\xff")
        self._ord_table = ord_table

        # The chr_table maps between ordinal position in the alphabet letters
        # and the ordinal position in ascii. This map is not the inverse of
        # ord_table if there are alternatives.
        chr_table = ["\x00"]*256
        for i,a in enumerate(letters) :
            chr_table[ i ] = a
        chr_table = "".join(chr_table)
        self._chr_table = chr_table

        return self


    def alphabetic(self, string) :
        """True if all characters of the string are in this alphabet."""
        table = self._ord_table
        for s in str(string):
            if table[ord(s)] == "\xff" :
                return False
        return True
        
    def chr(self, n) :
        """ The n'th character in the alphabet (zero indexed) or \\0 """
        return self._chr_table[n]

    def ord(self, c) :
        """The ordinal position of the character c in this alphabet,
        or 255 if no such character.
        """
        return ord(self._ord_table[ord(c)])
       
    def chrs(self, sequence_of_ints) :
        """Convert a sequence of ordinals into an alphabetic string."""
        if not isinstance(sequence_of_ints, array) :
            sequence_of_ints = array('B', sequence_of_ints)
        s = sequence_of_ints.tostring().translate(self._chr_table)
        return Seq(s, self)

    def ords(self, string) :
        """Convert an alphabetic string into a byte array of ordinals."""
        string = str(string)
        s = string.translate(self._ord_table)
        a = array('B',s)
        return a

    
    def normalize(self, string) :
        """Normalize an alphabetic string by converting all alternative symbols 
        to the canonical equivalent in 'letters'.
        """
        if not self.alphabetic(string) :
            raise ValueError("Not an alphabetic string.")
        return self.chrs(self.ords(string))
    
    def letters(self) :
        """ Letters of the alphabet as a string."""
        return str(self)

    def _all_letters(self) :
        """ All allowed letters, including alternatives."""
        let = []
        let.append(self._letters)
        for key, value in self._alternatives :
            let.append(value)
        return ''.join(let)

    def __repr__(self) :
        return "Alphabet( '" + self._letters +"', zip"+ repr(self._alternatives)+" )" 
    
    def __str__(self) :
        return str(self._letters)

    def __len__(self) :
        return len(self._letters)

    def __eq__(self, other) :
        if not hasattr(other, "_ord_table") : return False
        return self._ord_table == other._ord_table

    def __ne__(self, other) :
        return not self.__eq__(other)

    def __iter__(self) :
        return iter(self._letters)

    def __getitem__(self, key) :
        return self._letters[key]

    @staticmethod 
    def which(seqs, alphabets=None) :
        """ Returns the most appropriate unambiguous protein, RNA or DNA alphabet
        for a Seq or SeqList. If a list of alphabets is supplied, then the best alphabet
        is selected from that list.

        The heuristic is to count the occurrences of letters for each alphabet and 
        downweight longer alphabets by the square root of the alphabet length. Ties
        go to the first alphabet in the list.

        """
        if alphabets is None :
            alphabets = (unambiguous_dna_alphabet,
                    unambiguous_rna_alphabet,
                    unambiguous_protein_alphabet,
                    )
        score= [sum(seqs.tally(a))/sqrt(len(a)) for a in alphabets]
        best = score.index(max(score))
        a = alphabets[best]
        return a

# End class Alphabet
        
# ------------------- Standard ALPHABETS -------------------
# Standard alphabets are defined here, after Alphabet class.

generic_alphabet = Alphabet(None, None)


protein_alphabet = Alphabet('ACDEFGHIKLMNOPQRSTUVWYBJZX*-', 
                        zip('acdefghiklmnopqrstuvwybjzx?.~',
                            'ACDEFGHIKLMNOPQRSTUVWYBJZXX--') )


nucleic_alphabet     =  Alphabet("ACGTURYSWKMBDHVN-", 
                            zip("acgturyswkmbdhvnXx?.~",
                                "ACGTURYSWKMBDHVNNNN--") )

dna_alphabet  =    Alphabet("ACGTRYSWKMBDHVN-", 
                        zip('acgtryswkmbdhvnXx?.~Uu', 
                            'ACGTRYSWKMBDHVNNNN--TT') )

rna_alphabet  =    Alphabet("ACGURYSWKMBDHVN-", 
                        zip('acguryswkmbdhvnXx?.~Tt', 
                            'ACGURYSWKMBDHVNNNN--UU') )

reduced_nucleic_alphabet  =  Alphabet("ACGTN-", 
                            zip('acgtryswkmbdhvnXx?.~TtRYSWKMBDHV', 
                                'ACGTNNNNNNNNNNNNNN--TTNNNNNNNNNN') )

reduced_protein_alphabet = Alphabet('ACDEFGHIKLMNPQRSTVWYX*-', 
                                zip('acdefghiklmnpqrstvwyx?.~BbZzUu',
                                    'ACDEFGHIKLMNPQRSTVWYXX--XXXXCC') )

unambiguous_dna_alphabet    =  Alphabet("ACGT", zip('acgt','ACGT') )

unambiguous_rna_alphabet    =  Alphabet("ACGU", zip('acgu','ACGU') )

unambiguous_protein_alphabet =  Alphabet("ACDEFGHIKLMNPQRSTVWY",
                        zip('acdefghiklmnopqrstuvwy',
                            'ACDEFGHIKLMNOPQRSTUVWY') )

   
_complement_table = maketrans("ACGTRYSWKMBDHVN-acgtUuryswkmbdhvnXx?.~",
                              "TGCAYRSWMKVHDBN-tgcaAayrswmkvhdbnXx?.~")



class Seq(str):
    """ An alphabetic string. A subclass of "str" consisting solely of
    letters from the same alphabet.

    Attributes:
        alphabet    -- A string or Alphabet of allowed characters.
        name        -- A short string used to identify the sequence.
        description -- A string describing the sequence   
        
    Authors :
        GEC 2005
    """
    # TODO: need a method to return a copy of the string with a new alphabet,
    # preserving the sequence, name and alphabet?
    
    def __new__(cls, obj, 
            alphabet= generic_alphabet, 
            name =None,  description=None,
            ):
        self = str.__new__(cls, obj)
        if alphabet is None: 
            alphabet = generic_alphabet
        if  not isinstance(alphabet, Alphabet): 
            alphabet = Alphabet(alphabet)
        if not alphabet.alphabetic(self) :
            raise ValueError("Sequence not alphabetic %s, '%s'" %(alphabet, self))
        
        self._alphabet=alphabet
        self.name = name
        self.description = description
                           
        return self
 
    # BEGIN PROPERTIES
            
    # Make alphabet constant 
    def _get_alphabet(self):
        return self._alphabet
    alphabet = property(_get_alphabet)     

    # END PROPERTIES        


    def ords(self) :
        """ Convert sequence to an array of integers 
        in the range [0, len(alphabet) ) 
        """
        return self.alphabet.ords(self) 
        
    def tally(self, alphabet = None):
        """Counts the occurrences of alphabetic characters.
                
        Arguments:
        - alphabet -- an optional alternative alphabet

        Returns :
            A list of character counts in alphabetic order.
        """
        # Renamed from count() since this conflicts with str.count().
        if not alphabet : alphabet = self.alphabet 
        L = len(alphabet)
        counts = [0,] * L
        
        ords = alphabet.ords(self) 
        
        for n in ords:
            if n<L : counts[n] +=1
        
        return counts
        
        

    def __getslice__(self, i, j):    
        cls = self.__class__
        return cls( str.__getslice__(self,i,j), self.alphabet)    
     
    def __getitem__(self, key) :
        cls = self.__class__
        return cls( str.__getitem__(self,key), self.alphabet)
        
    def __add__(self, other) :
        # called for "self + other"
        cls = self.__class__
        return cls( str.__add__(self, other), self.alphabet) 
    
    def __radd__(self, other) :
        # Called when "other + self" and other is superclass of self
        cls = self.__class__
        return cls( str.__add__(self, other), self.alphabet) 
    
    def join(self, str_list) :
        cls = self.__class__
        return cls( super(Seq, self).join(str_list), self.alphabet)        
    
    def __eq__(self, other) :
        if not hasattr(other, "alphabet") : return False
        if self.alphabet != other.alphabet :
            return False
        return str.__eq__(self, other)

    def __ne__(self, other) :
        return not self.__eq__(other)

    def tostring(self) :
        """ Converts Seq to a raw string. 
        """
        # Compatibility with biopython
        return str(self)

    # ---- Transformations of Seq ----
    def reverse(self) :
        """Return the reversed sequence. 
        
        Note that this method returns a new object, in contrast to
        the in-place reverse() method of list objects.
        """
        cls = self.__class__
        return cls( self[::-1], self.alphabet) 

    def ungap(self) :
        # FIXME: Gap symbols should be specified by the Alphabet?
        return self.remove( '-.~')

    def remove(self, delchars) :
        """Return a new alphabetic sequence with all characters in 'delchars'
         removed.
        """
        cls = self.__class__
        return cls( str(self).translate(maketrans('',''), delchars), self.alphabet) 

    def lower(self) :
        """Return a lower case copy of the sequence. """
        cls = self.__class__
        trans = maketrans('ABCDEFGHIJKLMNOPQRSTUVWXYZ','abcdefghijklmnopqrstuvwxyz')
        return cls(str(self).translate(trans), self.alphabet)
        
    def upper(self) :
        """Return a lower case copy of the sequence. """
        cls = self.__class__
        trans = maketrans('abcdefghijklmnopqrstuvwxyz','ABCDEFGHIJKLMNOPQRSTUVWXYZ')
        return cls(str(self).translate(trans), self.alphabet)
        
    def mask(self, letters= 'abcdefghijklmnopqrstuvwxyz', mask='X') :
        """Replace all occurrences of letters with the mask character.
        The default is to replace all lower case letters with 'X'.
        """
        LL = len(letters)
        if len(mask) !=1 : 
            raise ValueError("Mask should be single character") 
        to = mask * LL
        trans = maketrans( letters, to)
        cls = self.__class__
        return cls(str(self).translate(trans), self.alphabet)
 
    def translate(self) :
        """Translate a nucleotide sequence to a polypeptide using full
        IUPAC ambiguities in DNA/RNA and amino acid codes, using the
        standard genetic code. See corebio.transform.GeneticCode for
        details and more options.
        """
        # Note: masks str.translate
        from transform import GeneticCode
        return GeneticCode.std().translate(self)

    def back_translate(self) :
        """Translate a protein sequence back into coding DNA, using the
        standard genetic code. See corebio.transform.GeneticCode for
        details and more options.
        """
        from transform import GeneticCode
        return GeneticCode.std().back_translate(self)

       
    def reverse_complement(self) :
        """Returns reversed complementary nucleic acid sequence (i.e. the other
        strand of a DNA sequence.) 
        """
        return self.reverse().complement()
       
    def complement(self) :
        """Returns complementary nucleic acid sequence."""
        if not nucleic_alphabet.alphabetic(self.alphabet):
            raise ValueError("Incompatable alphabets")
        s = str.translate(self, _complement_table)
        cls = self.__class__
        return cls(s, self.alphabet, self.name, self.description) 
 
    def words(self, k, alphabet=None) :
        """Return an iteration over all subwords of length k in the sequence. If an optional
        alphabet is provided, only words from that alphabet are returned.
        
        >>> list(Seq("abcabc").words(3))
        ['abc', 'bca', 'cab', 'abc']
        """
        
        if len(self) < k : return
        
        # An optimization. Chopping up strings is faster.
        seq = self.alphabet.normalize(self).tostring()
        #seq = self.tostring() 
    
        for i in range(0,len(seq)-k+1) :
            word = seq[i:i+k]
            if alphabet is None or alphabet.alphabetic(word) :
                yield word

    def word_count(self, k, alphabet=None):
        """Return a count of all subwords in the sequence.
        
        >>> from corebio.seq import *
        >>> Seq("abcabc").word_count(3)
        [('abc', 2), ('bca', 1), ('cab', 1)]
        """
        from corebio.utils import group_count
        words = sorted(self.words(k,alphabet))
        return group_count(words)
 
        
# end class Seq


class SeqList(list):
    """ A list of sequences. 
    """

    __slots__ =["alphabet", "name", "description"]

    def __init__(self, alist=[], alphabet=None, name=None, description=None):
        list.__init__(self, alist)
        self.alphabet = alphabet
        self.name = name
        self.description = description

    # TOOWTDI. Replicates seq_io.read()
    #@classmethod
    #def read(cls, afile, alphabet = None):
    #    return corebio.seq_io.read(afile, alphabet)
    #read = classmethod(read)    
     
    def isaligned(self) :
        """Are all sequences of the same length and alphabet?"""
        if len(self)==0: return True
        A = self.alphabet
        if A is None : A = self[0].alphabet
        L = len(self[0])
        
        for s in self:
            if len(s)!=L : return False
            if s.alphabet != A : return False
        return True
        
        
         
    def ords(self, alphabet=None) :
        """ Convert sequence list into a 2D array of ordinals.
        """
        if not alphabet : alphabet = self.alphabet
        if not alphabet : raise ValueError("No alphabet")
        k = []
        for s in self:
            k.append( alphabet.ords(s) )
        return k
 
    def tally(self, alphabet = None):
        """Counts the occurrences of alphabetic characters.

        Arguments:
            - alphabet -- an optional alternative alphabet

        Returns :
        A list of character counts in alphabetic order.
        """
        if not alphabet : alphabet = self.alphabet
        if not alphabet : raise ValueError("No alphabet")

        counts = [sum(c) for c in zip(* [ s.tally(alphabet) for s in self])]
        return counts
 
        
    def profile(self, alphabet = None):
        """Counts the occurrences of characters in each column.

        Returns: Motif(counts, alphabet)
        """
        if not alphabet : alphabet = self.alphabet
        if not alphabet : raise ValueError("No alphabet")        
        
        N = len(alphabet) 
        ords = self.ords(alphabet)
        L = len(ords[0])
        counts = [ [0,]*N for l in range(0,L)]
        
        for o in ords :
            if len(o)!=L : raise ValueError("Sequences are of incommensurate lengths. Cannot tally.")
            for j,n in enumerate(o) :
                if n<N : counts[ j][n] +=1
 
		from corebio.matrix import Motif

        return Motif(alphabet, counts)
# end class SeqList


def dna(string) :
    """Create an alphabetic sequence representing a stretch of DNA.    
    """
    return Seq(string, alphabet = dna_alphabet)
    
def rna(string) :
    """Create an alphabetic sequence representing a stretch of RNA.    
    """
    return Seq(string, alphabet = rna_alphabet)

def protein(string) :
    """Create an alphabetic sequence representing a stretch of polypeptide.    
    """
    return Seq(string, alphabet = protein_alphabet)



 
