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
- GeneticCode -- The genetic mapping of DNA to protein.

Functions :
-  mask_low_complexity -- Implementation of Seg algorithm to remove low complexity  
        regions from protein sequences.

Other:
-   reduced_protein_alphabets -- A dictionary of transforms that reduce the size of the protein alphabet,
        merging various amino acids into classes.
        
        "LiBn" where n is 2 to 19 are from Li et al (2003), table I, 2 to 19 groups.
        "LiBn" where n is 2 to 19 are from Li et al (2003), table II (no interlacing), 2 to 19 groups.
        
        Ref: Li et al Reduction of protein sequence complexity by residue grouping, 
            Prot. Eng. 16 323-330 (2003)  

    
"""


from corebio.data import dna_extended_letters, dna_ambiguity
from corebio.seq import Seq, protein_alphabet, nucleic_alphabet, dna_alphabet, Alphabet

from corebio.seq import reduced_protein_alphabet as std_protein_alphabet 
from string import maketrans
from corebio.moremath import log2 , entropy

__all__ = [
    'Transform',
    'mask_low_complexity',
    'GeneticCode',
    'reduced_protein_alphabets'
    ]

class Transform(object) :
    """A translation between alphabetic strings.
    (This class is not called 'Translation' to avoid confusion with the
    biological translation of RNA to protein.)
    
    Example:
    trans = Transform( 
        Seq("ACGTRYSWKMBDHVN-acgtUuryswkmbdhvnXx?.~", dna_alphabet),                    
        Seq("ACGTRYSWKMNNNNN-acgtUuryswkmbnnnnXx?.~", reduced_nucleic_alphabet)         
        )
    s0 = Seq("AAAAAV", nucleic_alphabet)
    s1 = trans(s0)              
    assert(s1.alphabet == reduced_nucleic_alphabet)
    assert(s2 == Seq("AAAAAN",  reduced_nucleic_alphabet)
        
    Status : Beta 
    """
    
    __slots__ = ["table", "source", "target","name","description"]
    def __init__(self, source, target, name=None,description=None) :

        self.table = maketrans(source, target)
        self.source = source
        self.target = target
        self.name = name
        self.description = description
     
    def __call__(self, seq) :
        """Translate sequence."""
        if not self.source.alphabet.alphabetic(seq) :
            raise ValueError("Incompatible alphabets")
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
    less than the trigger complexity, or have an entropy less than the extension    
    complexity and neighbor other low-complexity windows. The sequence within   
    a low complexity region is replaced with the mask character (default 'X'), 
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

    segged = seq.alphabet.chrs(s)
    segged.name =seq.name
    segged.description = seq.description
    return segged
    
    
    
    
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
        -- ident - Standard identifier (or zero). An integer
        -- description 
        -- amino acid - A sequence of amino acids and stop codons. e.g.
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
            c = C.replace('U', 'T') # Convert RNA codon to DNA
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
        -- Seq - A DNA sequence
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
        letters = "TCAG" # Conventional ordering for codon tables.
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
    


reduced_protein_alphabets = {
#
    "LiB2": Transform(
    	Seq("CFYWMLIV-GPATSNHQEDRKX*-", std_protein_alphabet  ),
    	Seq("IIIIIIII-SSSSSSSSSSSSX*-", Alphabet("ISX*-") ),
    	"Li et al (2003), table II, group 2"),
#
    "LiB3": Transform(
    	Seq("CFYWMLIV-GPATS-NHQEDRKX*-", std_protein_alphabet  ),
    	Seq("IIIIIIII-SSSSS-EEEEEEEX*-", Alphabet("ISEX*-") ),
    	"Li et al (2003), table II, group 3"),
#
    "LiB4": Transform(
    	Seq("CFYW-MLIV-GPATS-NHQEDRKX*-", std_protein_alphabet  ),
    	Seq("YYYY-IIII-SSSSS-EEEEEEEX*-", Alphabet("YISEX*-") ),
    	"Li et al (2003), table II, group 4"),	
#
    "LiB5": Transform(
    	Seq("CFYW-MLIV-G-PATS-NHQEDRKX*-", std_protein_alphabet ) ,
    	Seq("YYYY-IIII-G-SSSS-EEEEEEEX*-", Alphabet("YIGSEX*-") ),
    	"Li et al (2003), table II, group 5"),	
#
    "LiB6": Transform(
    	Seq("CFYW-MLIV-G-P-ATS-NHQEDRKX*-", std_protein_alphabet ) ,
    	Seq("YYYY-IIII-G-P-SSS-EEEEEEEX*-", Alphabet("YIGPSEX*-") ),
    	"Li et al (2003), table II, group 6"),	
#
    "LiB7": Transform(
    	Seq("CFYW-MLIV-G-P-ATS-NHQED-RKX*-", std_protein_alphabet  ),
    	Seq("YYYY-IIII-G-P-SSS-EEEEE-KKX*-", Alphabet("YIGPSEKX*-") ),
    	"Li et al (2003), table II, group 7"),	
#
    "LiB8": Transform(
    	Seq("CFYW-MLIV-G-P-ATS-NH-QED-RKX*-", std_protein_alphabet  ),
    	Seq("YYYY-IIII-G-P-SSS-NN-EEE-KKX*-", Alphabet("YIGPSNEKX*-") ),
    	"Li et al (2003), table II, group 8"),	
#
    "LiB9": Transform(
    	Seq("CFYW-ML-IV-G-P-ATS-NH-QED-RKX*-", std_protein_alphabet  ),
    	Seq("YYYY-LL-II-G-P-SSS-NN-EEE-KKX*-", Alphabet("YLIGPSNEKX*-") ),
    	"Li et al (2003), table II, group 9"),	
#
    "LiB10": Transform(
    	Seq("C-FYW-ML-IV-G-P-ATS-NH-QED-RKX*-", std_protein_alphabet  ),
    	Seq("C-YYY-LL-II-G-P-SSS-NN-EEE-KKX*-", Alphabet("CYLIGPSNEKX*-") ),
    	"Li et al (2003), table II, group 10"),	
#
    "LiB11": Transform(
    	Seq("C-FYW-ML-IV-G-P-A-TS-NH-QED-RKX*-", std_protein_alphabet  ),
    	Seq("C-YYY-LL-II-G-P-A-SS-NN-EEE-KKX*-", Alphabet("CYLIGPASNEKX*-") ),
    	"Li et al (2003), table II, group 11"),	
#
    "LiB12": Transform(
    	Seq("C-FYW-ML-IV-G-P-A-TS-NH-QE-D-RKX*-", std_protein_alphabet  ),
    	Seq("C-YYY-LL-II-G-P-A-SS-NN-EE-D-KKX*-", Alphabet("CYLIGPASNEDKX*-") ),
    	"Li et al (2003), table II, group 12"),	
#
    "LiB13": Transform(
    	Seq("C-FYW-ML-IV-G-P-A-T-S-NH-QE-D-RKX*-", std_protein_alphabet ),
    	Seq("C-YYY-LL-II-G-P-A-T-S-NN-EE-D-KKX*-", Alphabet("CYLIGPATSNEDKX*-") ),
    	"Li et al (2003), table II, group 13"),	
#
    "LiB14": Transform(
    	Seq("C-FYW-ML-IV-G-P-A-T-S-N-H-QE-D-RKX*-", std_protein_alphabet  ),
    	Seq("C-YYY-LL-II-G-P-A-T-S-N-H-EE-D-KKX*-", Alphabet("CYLIGPATSNHEDKX*-") ),
    	"Li et al (2003), table II, group 14"),	
#
    "LiB15": Transform(
    	Seq("C-FYW-ML-IV-G-P-A-T-S-N-H-QE-D-R-KX*-", std_protein_alphabet  ),
    	Seq("C-YYY-LL-II-G-P-A-T-S-N-H-EE-D-R-KX*-", Alphabet("CYLIGPATSNHEDRKX*-") ),
    	"Li et al (2003), table II, group 15"),	
#
    "LiB16": Transform(
    	Seq("C-FY-W-ML-IV-G-P-A-T-S-N-H-QE-D-R-KX*-", std_protein_alphabet  ),
    	Seq("C-YY-W-LL-II-G-P-A-T-S-N-H-EE-D-R-KX*-", Alphabet("CYWLIGPATSNHEDRKX*-") ),
    	"Li et al (2003), table II, group 16"),	
#
    "LiB17": Transform(
    	Seq("C-FY-W-ML-IV-G-P-A-T-S-N-H-Q-E-D-R-KX*-", std_protein_alphabet  ),
    	Seq("C-YY-W-LL-II-G-P-A-T-S-N-H-Q-E-D-R-KX*-", Alphabet("CYWLIGPATSNHQEDRKX*-") ),
    	"Li et al (2003), table II, group 17"),	
#
    "LiB18": Transform(
    	Seq("C-FY-W-M-L-IV-G-P-A-T-S-N-H-Q-E-D-R-KX*-", std_protein_alphabet  ),
    	Seq("C-YY-W-M-L-II-G-P-A-T-S-N-H-Q-E-D-R-KX*-", Alphabet("CYWMLIGPATSNHQEDRKX*-") ),
    	"Li et al (2003), table II, group 18"),	
#
    "LiB19": Transform(
    	Seq("C-F-Y-W-M-L-IV-G-P-A-T-S-N-H-Q-E-D-R-KX*-", std_protein_alphabet  ),
    	Seq("C-F-Y-W-M-L-II-G-P-A-T-S-N-H-Q-E-D-R-KX*-", Alphabet("CFYWMLIGPATSNHQEDRKX*-") ),
    	"Li et al (2003), table II, group 19"),	
#
    "LiB19": Transform(
        Seq("C-F-Y-W-M-L-I-V-G-P-A-T-S-N-H-Q-E-D-R-KX*-", std_protein_alphabet  ),
        Seq("C-F-Y-W-M-L-I-V-G-P-A-T-S-N-H-Q-E-D-R-KX*-", Alphabet("CFYWMLIVGPATSNHQEDRKX*-") ),
            "Li et al (2003), table II, group 20"),	
    	
#
    "LiA2": Transform(
    	Seq("CMFILVWY-AGTSNQDEHRKPX*-", std_protein_alphabet  ),
    	Seq("IIIIIIII-SSSSSSSSSSSSX*-", Alphabet("ISX*-") ),
    	"Li et al (2003), table I, group 2"),
#
    "LiA3": Transform(
    	Seq("CMFILVWY-AGTSP-NQDEHRKX*-", std_protein_alphabet  ),
    	Seq("IIIIIIII-SSSSS-EEEEEEEX*-", Alphabet("ISEX*-") ),
    	"Li et al (2003), table I, group 3"),
#
    "LiA4": Transform(
    	Seq("CMFWY-ILV-AGTS-NQDEHRKPX*-", std_protein_alphabet  ),
    	Seq("YYYYY-III-SSSS-EEEEEEEEX*-", Alphabet("YISEX*-") ),
    	"Li et al (2003), table I, group 4"),	
#
    "LiA5": Transform(
    	Seq("FWYH-MILV-CATSP-G-NQDERKX*-", std_protein_alphabet ) ,
    	Seq("YYYY-IIII-SSSSS-G-EEEEEEX*-", Alphabet("YISGEX*-") ),
    	"Li et al (2003), table I, group 5"),	
#
    "LiA6": Transform(
    	Seq("FWYH-MILV-CATS-P-G-NQDERKX*-", std_protein_alphabet ) ,
    	Seq("YYYY-IIII-SSSS-P-G-EEEEEEX*-", Alphabet("YISPGEX*-") ),
    	"Li et al (2003), table I, group 6"),	
#
    "LiA7": Transform(
    	Seq("FWYH-MILV-CATS-P-G-NQDE-RKX*-", std_protein_alphabet  ),
    	Seq("YYYY-IIII-SSSS-P-G-EEEE-KKX*-", Alphabet("YISPGEKX*-") ),
    	"Li et al (2003), table I, group 7"),	
#
    "LiA8": Transform(
    	Seq("FWYH-MILV-CA-NTS-P-G-DE-QRKX*-", std_protein_alphabet  ),
    	Seq("YYYY-IIII-AA-SSS-P-G-NN-KKKX*-", Alphabet("YIASPGNKX*-") ),
    	"Li et al (2003), table I, group 8"),	
#
    "LiA9": Transform(
    	Seq("FWYH-ML-IV-CA-NTS-P-G-DE-QRKX*-", std_protein_alphabet  ),
    	Seq("YYYY-LL-VV-AA-SSS-P-G-NN-KKKX*-", Alphabet("YLVASPGNKX*-") ),
    	"Li et al (2003), table I, group 9"),	
#
    "LiA10": Transform(
    	Seq("FWY-ML-IV-CA-TS-NH-P-G-DE-QRKX*-", std_protein_alphabet  ),
    	Seq("YYY-LL-VV-AA-TT-NN-P-G-DD-KKKX*-", Alphabet("YLVATNPGDKX*-") ),
    	"Li et al (2003), table I, group 10"),	
#
    "LiA11": Transform(
    	Seq("FWY-ML-IV-CA-TS-NH-P-G-D-QE-RKX*-", std_protein_alphabet  ),
    	Seq("YYY-LL-VV-AA-TT-NN-P-G-D-EE-KKX*-", Alphabet("YLVATNPGDEKX*-") ),
    	"Li et al (2003), table I, group 11"),	
#
    "LiA12": Transform(
    	Seq("FWY-ML-IV-C-A-TS-NH-P-G-D-QE-RKX*-", std_protein_alphabet  ),
    	Seq("YYY-LL-VV-C-A-TT-NN-P-G-D-EE-KKX*-", Alphabet("YLVCATNPGDEKX*-") ),
    	"Li et al (2003), table I, group 12"),	
#
    "LiA13": Transform(
    	Seq("FWY-ML-IV-C-A-T-S-NH-P-G-D-QE-RKX*-", std_protein_alphabet ),
    	Seq("YYY-LL-VV-C-A-T-S-NN-P-G-D-EE-KKX*-", Alphabet("YLVCATSNPGDEKX*-") ),
    	"Li et al (2003), table I, group 13"),	
#
    "LiA14": Transform(
    	Seq("FWY-ML-IV-C-A-T-S-NH-P-G-D-QE-R-KX*-", std_protein_alphabet  ),
    	Seq("YYY-LL-VV-C-A-T-S-NN-P-G-D-EE-R-KX*-", Alphabet("YLVCATSNPGDERKX*-") ),
    	"Li et al (2003), table I, group 14"),	
#
    "LiA15": Transform(
    	Seq("FWY-ML-IV-C-A-T-S-N-H-P-G-D-QE-R-KX*-", std_protein_alphabet  ),
    	Seq("YYY-LL-VV-C-A-T-S-N-H-P-G-D-EE-R-KX*-", Alphabet("YLVCATSNHPGDERKX*-") ),
    	"Li et al (2003), table I, group 15"),	
#
    "LiA16": Transform(
    	Seq("W-FY-ML-IV-C-A-T-S-N-H-P-G-D-QE-R-KX*-", std_protein_alphabet  ),
    	Seq("W-YY-LL-VV-C-A-T-S-N-H-P-G-D-EE-R-KX*-", Alphabet("WYLVCATSNHPGDERKX*-") ),
    	"Li et al (2003), table I, group 16"),	
#
    "LiA17": Transform(
    	Seq("W-FY-ML-IV-C-A-T-S-N-H-P-G-D-Q-E-R-KX*-", std_protein_alphabet  ),
    	Seq("W-YY-LL-VV-C-A-T-S-N-H-P-G-D-Q-E-R-KX*-", Alphabet("WYLVCATSNHPGDQERKX*-") ),
    	"Li et al (2003), table I, group 17"),	
#
    "LiA18": Transform(
    	Seq("W-FY-M-L-IV-C-A-T-S-N-H-P-G-D-Q-E-R-KX*-", std_protein_alphabet  ),
    	Seq("W-YY-M-L-VV-C-A-T-S-N-H-P-G-D-Q-E-R-KX*-", Alphabet("WYMLVCATSNHPGDQERKX*-") ),
    	"Li et al (2003), table I, group 18"),	
#
    "LiA19": Transform(
    	Seq("W-F-Y-M-L-IV-C-A-T-S-N-H-P-G-D-Q-E-R-KX*-", std_protein_alphabet  ),
    	Seq("W-F-Y-M-L-VV-C-A-T-S-N-H-P-G-D-Q-E-R-KX*-", Alphabet("WFYMLVCATSNHPGDQERKX*-") ),
    	"Li et al (2003), table I, group 19"),
#
    "LiA20": Transform(
        Seq("W-F-Y-M-L-I-V-C-A-T-S-N-H-P-G-D-Q-E-R-KX*-", std_protein_alphabet  ),
        Seq("W-F-Y-M-L-I-V-C-A-T-S-N-H-P-G-D-Q-E-R-KX*-", Alphabet("WFYMLIVCATSNHPGDQERKX*-") ),
        "Li et al (2003), table I, group 20"),

}
 
        
           
