#  Copyright (c) 2006 John Gilman
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

""" Transformations of Seqs (alphabetic sequences).



Classes :
- Transform   -- Simple transforms of alphabetic strings.
- GeneticCode -- The genetic mapping of dna to protein.

Functions :
-  mask_low_complexity -- Implementation of Seg algorithm to remove low complexity  
        regions from protein sequences.
   
    
"""


from corebio.data import dna_extended_letters, dna_ambiguity
from corebio.seq import Seq, protein_alphabet, nucleic_alphabet, dna_alphabet
from string import maketrans
from corebio.moremath import log2 , entropy

__all__ = [
    'Transform',
    'mask_low_complexity',
    'GeneticCode'
    ]

class Transform(object) :
    """A translation between alphabetic strings.
    (This class is not called 'Translation' to avoid confusion with the
    biological translation of rna to protein.)
    
    Example:
    trans = Transform( 
        Seq("ACGTRYSWKMBDHVN-acgtUuryswkmbdhvnXx?.~'", dna_alphabet),                    
        Seq("ACGTRYSWKMNNNNN-acgtUuryswkmbnnnnXx?.~", reduced_nucleic_alphabet)         
        )
    s0 = Seq("AAAAAV", nucleic_alphabet)
    s1 = trans(s0)              
    assert(s1.alphabet == reduced_nucleic_alphabet)
    assert(s2 == Seq("AAAAAN",  reduced_nucleic_alphabet)
        
    Status : Beta 
    """
    
    __slots__ = ["table", "source", "target"]
    def __init__(self, source, target) :

        self.table = maketrans(source, target)
        self.source = source
        self.target = target
     
     
    def __call__(self, seq) :
        """Translate sequence."""
        if not self.source.alphabet.alphabetic(seq) :
            raise ValueError("Incompatable alphabets")
        s = str.translate(seq, self.table)
        cls = self.target.__class__
        return cls(s, self.target.alphabet, seq.name, seq.description)
# End class Translation

# FIXME: Test, document, add to seq.
dna_complement = Transform(
        Seq("ACGTRYSWKMBDHVN-acgtUuryswkmbdhvnXx?.~", dna_alphabet),  
        Seq("TGCAYRSWMKVHDBN-tgcaAayrswmkvhdbnXx?.~", dna_alphabet),  
    )

     

def mask_low_complexity(seq, width =12, trigger=1.8, extension=2.0, mask='X') :
    """ Mask low complexity regions in protein sequences.
    
    Uses the method of Seg [1] by Wootton & Federhen [2] to divide a sequence   
    into regions of high and low complexity. The sequence is divided into
    overlapping windows. Low complexity windows either have a sequence entropy
    less that the trigger complexity, or have an entropy less than the extension    
    complexity and neighbor other low-complexity windows. The sequence within   
    low complexity regions are replaced with the mask character (default 'X'), 
    and the masked alphabetic sequence is returned.
    
    The default parameters, width=12, trigger=1.8, extension=2.0, mask='X' are
    suitable for masking protein sequences before a database search. The 
    standard default seg parameters are width=12, trigger=2.2, extension=2.5
    
    Arguments:
        Seq seq         -- An alphabetic sequence
        int width       -- Window width
        float trigger   -- Entropy in bits between 0 and 4.3.. ( =log_2(20) )
        float extension -- Entropy in bits between 0 and 4.3.. ( =log_2(20) )
        char mask       -- The mask character (default: 'X') 
    Returns :
        Seq         -- A masked alphabetic sequence
    Raises :
        ValueError  -- On invalid arguments
    Refs:
        [1] seg man page: 
            http://bioportal.weizmann.ac.il/education/materials/gcg/seg.html
        [2] Wootton & Federhen (Computers and Chemistry 17; 149-163, (1993)) 
    Authors:
        GEC 2005
    Future :
        - Optional mask character.
        - Option to lower case masked symbols.
        - Remove arbitary restriction to protein.
    """
    
    lg20 = log2(20)
    if trigger<0 or trigger>lg20 :
        raise ValueError("Invalid trigger complexity: %f"% trigger) 
    if extension<0 or extension>lg20 or extension<trigger:
        raise ValueError("Invalid extension complexity: %f"% extension)
    if width<0 :
        raise ValueError("Invalid width: %d"% width)

    if width > len(seq) : return seq
    
    s = seq.ords()

    X = seq.alphabet.ord(mask)

    
    nwindows = len(seq)- width +1
    ent = [ 0 for x in range(0, nwindows)]
    count = [ 0 for x in range(0, len(seq.alphabet) )]
    
    for c in s[0:width] : count[c] +=1
    ent[0] = entropy(count,2)
    
    for i in range(1, nwindows) :
        count[ s[i-1] ] -= 1
        count[ s[i+width-1] ] +=1
        ent[i] = entropy(count,2)
    
    prev_segged = False 
    for i in range(0, nwindows) :
        if ((prev_segged and ent[i]< extension) or 
            ent[i]< trigger) :
            for j in range(0, width) : s[i+j]=X
            prev_segged=True
        else :
            prev_segged = False


    # Redo, only backwards
    prev_segged = False 
    for i in range(nwindows-1, -1, -1) :
        if ((prev_segged and ent[i]< extension) or 
            ent[i]< trigger) :
            for j in range(0, width) : s[i+j]=X
            prev_segged=True
        else :
            prev_segged = False

    
    return  seq.alphabet.chrs(s)
# end mask_low_complexity()    
     

class GeneticCode(object):
    """An encoding of amino acids by DNA triplets.
 
    Example : 
    
    Genetic Code [1]: Standard   
          T         C         A         G      
       +---------+---------+---------+---------+
     T | TTT F   | TCT S   | TAT Y   | TGT C   | T
     T | TTC F   | TCC S   | TAC Y   | TGC C   | C
     T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
     T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
       +---------+---------+---------+---------+
     C | CTT L   | CCT P   | CAT H   | CGT R   | T
     C | CTC L   | CCC P   | CAC H   | CGC R   | C
     C | CTA L   | CCA P   | CAA Q   | CGA R   | A
     C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
       +---------+---------+---------+---------+
     A | ATT I   | ACT T   | AAT N   | AGT S   | T
     A | ATC I   | ACC T   | AAC N   | AGC S   | C
     A | ATA I   | ACA T   | AAA K   | AGA R   | A
     A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
       +---------+---------+---------+---------+
     G | GTT V   | GCT A   | GAT D   | GGT G   | T
     G | GTC V   | GCC A   | GAC D   | GGC G   | C
     G | GTA V   | GCA A   | GAA E   | GGA G   | A
     G | GTG V   | GCG A   | GAG E   | GGG G   | G
       +---------+---------+---------+---------+

    
    See Also :
    -- http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c
    -- http://www.ncbi.nlm.nih.gov/projects/collab/FT/index.html#7.5
    Authors:
        JXG, GEC
    """
    # TODO: Explain use of '?' in translated sequence.
    # TODO: Does translate fails with aproriate execption when fed gaps?
    # TODO: Can back_translate handle gaps?
    
    def __init__(self, ident, description,
        amino_acid, start, base1, base2, base3):
        """Create a new GeneticCode.

        Args:
        -- ident - Standarad identifier (Or zero). An integer
        -- description 
        -- amino acid - A sequecne of amino acids and stop codons. e.g.
            "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
        -- start - A sequence indicating start codons, e.g.,
            "---M---------------M---------------M----------------------------"
        -- base1 - The first base of each codon. e.g., 
            "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
        -- base2 - The second base of each codon. e.g.,
            "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
        -- base3 - The last base of each codon. e.g., 
            "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"            
        """
        self.ident = ident
        self.description = description
   
        self.amino_acid = amino_acid
        self.start = start
        self.base1 = base1
        self.base2 = base2
        self.base3 = base3
    
        stop_codons = []
        start_codons = []
        for i, a in enumerate(amino_acid) :
            codon = base1[i] + base2[i] + base3[i]
            if a=='*' :  stop_codons.append(codon)
            if start[i] == 'M': start_codons.append(codon)
            
        self.stop_codons = tuple(stop_codons)
        self.start_codons = tuple(start_codons)
        
        # Building the full translation table is expensive,
        # so we avoid doing so until necessary.
        self._table = None
        self._back_table = None

    #@staticmethod
    def std_list():
        "Return a list of standard genetic codes."
        return _codon_tables
    std_list = staticmethod(std_list)

    #@staticmethod
    def std():
        "The standard 'universal' genetic code."
        return _codon_tables[0]
    std = staticmethod(std)

    
    #@staticmethod
    def by_name(name) :
        """Find a genetic code in the code list by name or identifier.
        """
        for t in _codon_tables :
            if t.ident == name or t.description == name :
                return t
        raise ValueError("No such translation table: %s" % str(name) )
    by_name = staticmethod(by_name)        
    
  
    def _get_table(self) :
        if self._table is None : self._create_table() 
        return self._table  
    table = property(_get_table, None, "A map between codons and amino acids")

    def _get_back_table(self) :
        if self._back_table is None : 
            self._create_table() 
        return self._back_table  
    back_table = property(_get_back_table, None, "A map between amino acids and codons")


    def _create_table(self) :
        aa = self.amino_acid
        base1 = self.base1
        base2 = self.base2
        base3 = self.base3
        
        # Construct a table of unambiguous codon translations
        table = {}
        for i, a in enumerate(aa) :
            codon = base1[i] + base2[i] + base3[i]
            table[codon] = a
        
        # Build the back table.
        back_table = {}
        items = table.items()
        items.sort()
        for codon, aa in items[::-1] :
            back_table[aa] = codon   # Use first codon, alphabetically.
        back_table['X'] = 'NNN'
        back_table['B'] = 'NNN'         
        back_table['Z'] = 'NNN'         
        back_table['J'] = 'NNN'         
        self._back_table = back_table
                
        ltable = {}
        letters = dna_extended_letters+'U' # include RNA in table

        # Create a list of all possble codons
        codons = []
        for c1 in letters:
            for c2 in letters:
                for c3 in letters :
                    codons.append( c1+c2+c3)        

        # For each ambiguous codon, construct all compatible unambiguous codons.
        # Translate and collect a set of all possible translated amino acids.
        # If more than one translation look for possible amino acid ambiguity       
        # codes. 
        for C in codons :
            translated = dict() # Use dict, because no set in py2.3
            c = C.replace('U', 'T') # Convert rna codon to dna
            for c1 in dna_ambiguity[c[0]]:
                for c2 in dna_ambiguity[c[1]]:
                    for c3 in dna_ambiguity[c[2]]:
                        aa = table[ c1+c2+c3 ]
                        translated[aa] = ''
            translated = list(translated.keys())
            translated.sort()
            if len(translated) ==1 :
                trans = list(translated)[0]
            elif translated == ['D','N'] :
                trans = 'B'      
            elif translated == ['E','Q'] :
                trans = 'Z' 
            elif translated == ['I','L'] :
                trans = 'J'             
            elif '*' in translated:
                trans = '?'
            else :
                trans = 'X'
            ltable[C] = trans

        self._table = ltable
    # End create tables

    def translate(self, seq, frame=0) :
        """Translate a DNA sequence to a polypeptide using full
        IUPAC ambiguities in DNA/RNA and amino acid codes.
        
        Returns : 
        -- Seq - A polypeptide sequence 
        """
        # TODO: Optimize.
        # TODO: Insanity check alphabet.
        seq = str(seq)
        table = self.table
        trans = []
        L = len(seq)
        for i in range(frame, L-2, 3) :   
            codon = seq[i:i+3].upper()
            trans.append( table[codon])
        return Seq(''.join(trans), protein_alphabet)

             
    def back_translate(self, seq) :
        """Convert protein back into coding DNA.
        
        Args:
        -- seq - A polypeptide sequence.
        
        Returns :
        -- Seq - A dna sequence
        """
        # TODO: Optimzie
        # TODO: Insanity check alphabet.
        table = self.back_table
        seq = str(seq)
        trans = [ table[a] for a in seq]
        return Seq(''.join(trans), dna_alphabet)
        
 #TODO: translate_orf(self, seq, start) ?
 #TODO: translate_to_stop(self, seq, frame) ?       
 #TODO: translate_all_frames(self,seq) -> 6 translations.
 
    def __repr__(self) :
        string = []
        string += 'GeneticCode( %d, "' % self.ident
        string += self.description
        string += '", \n'
        string += '    amino_acid = "'
        string += self.amino_acid
        string += '",\n'
        string += '    start =      "'
        string += self.start
        string += '",\n'
        string += '    base1 =      "'
        string += self.base1
        string += '",\n'
        string += '    base2 =      "'
        string += self.base2
        string += '",\n'
        string += '    base3 =      "'
        string += self.base3
        string += '" )'
        return ''.join(string)
        
        
    def __str__(self) :
        """Returns a text representation of this genetic code."""
        # Inspired by http://bugzilla.open-bio.org/show_bug.cgi?id=1963
        letters = "TCAG" # Convectional ordering for codon tables.
        string = []
        
        if self.ident :
            string += 'Genetic Code [%d]: ' % self.ident
        else :
            string += 'Genetic Code: '
        string +=  self.description or ''        

        string += "\n    "
        string += " ".join( ["  %s      " % c2 for c2 in letters] ) 
        
        string += "\n   +" 
        string +=  "+".join(["---------" for c2 in letters]) + "+  "
        
        table = self.table
        
        for c1 in letters :
            for c3 in letters :
                string += '\n '
                string  += c1 
                string  += " |"
                for c2 in letters :
                    codon = c1+c2+c3
                    string += " " + codon
                    if codon in self.stop_codons :
                        string  += " Stop|"
                    else :
                        amino = table.get(codon, '?')
                        if codon in self.start_codons :
                            string += " %s(s)|" % amino
                        else :
                            string += " %s   |" % amino
                string += " " + c3
                
            string += "\n   +"
            string += "+".join(["---------" for c2 in letters]) 
            string += "+  "
        string += '\n'
        return ''.join(string)
# end class GeneticCode

 
# Data from http://www.ncbi.nlm.nih.gov/projects/collab/FT/index.html#7.5
# Aug. 2006
# Genetic Code Tables
# 
# Authority      International Sequence Databank Collaboration
# Contact        NCBI
# Scope          /transl_table qualifier
# URL            http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c   
_codon_tables = ( 
    GeneticCode(1, "Standard",
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "---M---------------M---------------M----------------------------",
        "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
        "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
        "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"),
        
    GeneticCode(2, "Vertebrate Mitochondrial",
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
        "--------------------------------MMMM---------------M------------",
        "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
        "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
        "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"),
 
    GeneticCode(3, "Yeast Mitochondrial",
        "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "----------------------------------MM----------------------------",
        "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
        "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
        "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"), 

    GeneticCode(4, "Mold, Protozoan, Coelenterate Mitochondrial & Mycoplasma/Spiroplasma",
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "--MM---------------M------------MMMM---------------M------------",
        "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
        "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
        "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"),
        
    GeneticCode(5, "Invertebrate Mitochondrial", 
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
        "---M----------------------------MMMM---------------M------------",
        "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
        "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
        "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"),

    GeneticCode(6, "Ciliate, Dasycladacean and Hexamita Nuclear",    
        "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "-----------------------------------M----------------------------",
        "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
        "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
        "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"),
 
    GeneticCode(9, "Echinoderm and Flatworm Mitochondrial", 
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        "-----------------------------------M---------------M------------",
        "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
        "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
        "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"),

    GeneticCode(10, "Euplotid Nuclear",
        "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "-----------------------------------M----------------------------",
        "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
        "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
        "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"),
  
    GeneticCode(11, "Bacterial and Plant Plastid",
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "---M---------------M------------MMMM---------------M------------",
        "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
        "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
        "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"),

    GeneticCode(12, "Alternative Yeast Nuclear",
        "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "-------------------M---------------M----------------------------",
        "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
        "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
        "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"),
                                                                   
    GeneticCode(13,"Ascidian Mitochondrial", 
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
        "-----------------------------------M----------------------------",
        "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
        "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
        "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"),

    GeneticCode(14, "Alternative Flatworm Mitochondrial",
        "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        "-----------------------------------M----------------------------",
        "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
        "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
        "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"),

    GeneticCode(15, "Blepharisma Nuclear",
        "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "-----------------------------------M----------------------------",
        "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
        "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
        "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"),

    GeneticCode(16, "Chlorophycean Mitochondrial",
        "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "-----------------------------------M----------------------------",
        "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
        "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
        "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"),

    GeneticCode(21, "Trematode Mitochondrial",
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        "-----------------------------------M---------------M------------",
        "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
        "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
        "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"),

    GeneticCode(22, "Scenedesmus obliquus Mitochondrial",
        "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "-----------------------------------M----------------------------",
        "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
        "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
        "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"),

    GeneticCode(23,"Thraustochytrium Mitochondrial",
        "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "--------------------------------M--M---------------M------------",
        "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
        "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
        "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG",),
    )
    
          
        
           