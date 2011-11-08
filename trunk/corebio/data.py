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

"""
Standard data used in computational biology.


To convert a property dictionary to a list :
>>> comp = [ amino_acid_composition[k] for k in amino_acid_letters]



Resources: 
    Various standard data files are included in the corebio distribution. These
    may be loaded with the data_string, data_stream or data_filename methods.
    A complete set of names is stored in 'resource_names'
 
BLOSUM Scoring Matrices
    Source: ftp://ftp.ncbi.nih.gov/repository/blocks/unix/blosum
    These are all new blast style with 1/3 bit scaling
    - blosum35
    - blosum45    
    - blosum62    
    - blosum40    
    - blosum50    
    - blosum80    
    - blosum100   

Other substitution scoring matrices:
    - dist20_comp 
    - pam250
    - pam120
    - vtml160
    
Description of database cross references :
    - dbxref.txt (http://www.expasy.org/cgi-bin/lists?dbxref.txt)

    
Attributes:
    - amino_acid_letters
        -- Standard codes for the 20 canonical amino acids, in alphabetic
        order.
        
    - amino_acid_alternative_letters
        -- Amino acid one letter codes, alphabetic by three letter codes.

    - amino_acid_extended_letters

    - dna_letters

    - dna_extended_letters

    - rna_letters
    
    - rna_extended_letters

    - dna_ambiguity 

    - rna_ambiguity
    
    - amino_acid_ambiguity
    
    - amino_acid_mass
        -- Monomer isotopically averaged molecular mass 
    
    - dna_mass
    
    - rna_mass
        
    - one_to_three      
        -- Map from standard 1 letter amino acid codes to standard three
        letter codes. 
        Ref: http://www.ebi.ac.uk/RESID/faq.html
      
    - standard_three_to_one
        -- Map from standard 3 letter amino acid codes to standard 1
        letter codes.
         
    - extended_three_to_one
        -- Map between three letter amino acid codes (first letter capitalized) 
        and standard one letter codes. This map contains many nonstandard three
        letter codes, used, for example, to specify chemically modified amino
        acids in PDB files.
        Ref: http://astral.berkeley.edu/ 
        Ref: http://www.ebi.ac.uk/RESID/faq.html

    - amino_acid_names

    - amino_acid_composition
        -- Average amino acid composition of proteins.
        Ref: McCaldon P., Argos P. Proteins 4:99-122 (1988).

    - kyte_doolittle_hydrophobicity 
        -- Kyte-Doolittle hydrophobicity scale.
        Ref: Kyte J., Doolittle R.F. J. Mol. Biol. 157:105-132 (1982)
        
    - nucleotide_names
    
    - amino_acid_accesible_surface_area
        -- Nominal maximum solvent accessoble area for unmodified amino acids,
        in square Angstroms.
        Ref: Sander & Rost, (1994), Proteins, 20:216-226


Status: Beta (Data needs to be proof checked.)    
"""

# FIXME: Proof check data
# FIXME: Add __all__

# The ExPasy ProtScale tool is a great source of amino acid properties.
# http://au.expasy.org/cgi-bin/protscale.pl       

from StringIO import StringIO
from corebio.utils import resource_string, resource_stream,resource_filename
import utils

# Explicitly list set of available data resources. We want to be able to access
# these resources in, for example, a webapp, without inadvertently allowing
# unrestricted read access to the local file system.

resource_names = [
    'blosum35',
    'blosum45',    
    'blosum62',    
    'blosum40',    
    'blosum50',    
    'blosum80',    
    'blosum100',   
    'dist20_comp', 
    'pam250',
    'pam120', 
    'dbxref.txt',
    'vtml160',
    ]

_resource_filenames = {
    'blosum35':    'data/blosum35.mat',
    'blosum45':    'data/blosum45.mat',    
    'blosum62':    'data/blosum62.mat',    
    'blosum40':    'data/blosum40.mat',    
    'blosum50':    'data/blosum50.mat',    
    'blosum80':    'data/blosum80.mat',    
    'blosum100':   'data/blosum100.mat',   
    'dist20_comp': 'data/dist20_comp.mat', 
    'pam250':      'data/pam250.mat',
    'pam120':      'data/pam120.mat', 
    'dbxref.txt' : 'data/dbxref.txt',
    'vtml160' :    'data/vtml160',
    }

# TODO: Substitution matrix parser, SeqMatrix.read
# _resource_parsers = {}

def data_string( name ): 
    """Load the specified resource as a string."""
    fn = _resource_filenames[name]
    return resource_string(__name__, fn , __file__)    

def data_stream( name ):
    """Provide an open file handle to the specified resource."""
    fn = _resource_filenames[name]
    return resource_stream(__name__, fn , __file__)    

def data_filename( name ): 
    """Provide a filename for the given resource in the local filesystem."""
    fn = _resource_filenames[name]
    return resource_filename(__name__, fn, __file__)            

#def data_object( name, parser = None) :
#    if parser is None : 
#        if name in _resource_parsers :
#            parser = _resource_parsers[name]
#        else :
#            parser = str    
#    return parser( data_stream(name) )


amino_acid_letters = "ACDEFGHIKLMNPQRSTVWY"

amino_acid_alternative_letters = "ARNDCQEGHILKMFPSTWYV"

amino_acid_extended_letters = "ACDEFGHIKLMNOPQRSTUVWYBJZX*-"


dna_letters = "GATC"
dna_extended_letters = "GATCRYWSMKHBVDN"

rna_letters = "GAUC"  
rna_extended_letters = "GAUCRYWSMKHBVDN"


dna_ambiguity = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "M": "AC",
    "R": "AG",
    "W": "AT",
    "S": "CG",
    "Y": "CT",
    "K": "GT",
    "V": "ACG",
    "H": "ACT",
    "D": "AGT",
    "B": "CGT",
    "X": "GATC",
    "N": "GATC",
}

rna_ambiguity = {
    "A": "A",
    "C": "C",
    "G": "G",
    "U": "U",
    "M": "AC",
    "R": "AG",
    "W": "AU",
    "S": "CG",
    "Y": "CU",
    "K": "GU",
    "V": "ACG",
    "H": "ACU",
    "D": "AGU",
    "B": "CGU",
    "X": "GAUC",
    "N": "GAUC",
}

amino_acid_ambiguity = {
    "A": "A",
    "B": "ND",
    "C": "C",
    "D": "D",
    "E": "E",
    "F": "F",
    "G": "G",
    "H": "H",
    "I": "I",
    "K": "K",
    "L": "L",
    "M": "M",
    "N": "N",
    "P": "P",
    "Q": "Q",
    "R": "R",
    "S": "S",
    "T": "T",
    "V": "V",
    "W": "W",
    "X": "ACDEFGHIKLMNPQRSTVWY",
    "Y": "Y",
    "Z": "QE",
    "J": "IL",
    'U': 'U',
    'O': 'O',
}


# Monomer isotopically averaged molecular mass 
# Data Checked GEC Nov 2006
amino_acid_mass = {
    "A": 89.09,
    "B" : 132.66,  # Averaged proportional to amino_acid_composition
    "C": 121.16,
    "D": 133.10,
    "E": 147.13,
    "F": 165.19,
    "G": 75.07,
    "H": 155.16, 
    "I": 131.18,
    "J": 131.18,
    "K": 146.19,
    "L": 131.18,
    "M": 149.21,
    "N": 132.12,
    # "O" : ???, # TODO
    "P": 115.13,
    "Q": 146.15,
    "R": 174.20,
    "S": 105.09,
    "T": 119.12,
    "U" : 168.05,
    "V": 117.15,
    "W": 204.23,
    "X" : 129.15, # Averaged proportional to amino_acid_composition  
    "Y": 181.19,
    "Z" : 146.76, # Averaged proportional to amino_acid_composition    
    }
    
dna_mass = {
    "A": 347.,
    "C": 323.,
    "G": 363.,
    "T": 322.,
    }

rna_mass = {
    "A": 363.,
    "C": 319.,
    "G": 379.,
    "U": 340.,
}

one_to_three = {
    'A':'Ala', 'B':'Asx', 'C':'Cys', 'D':'Asp',
    'E':'Glu', 'F':'Phe', 'G':'Gly', 'H':'His',
    'I':'Ile', 'K':'Lys', 'L':'Leu', 'M':'Met',
    'N':'Asn', 'P':'Pro', 'Q':'Gln', 'R':'Arg',
    'S':'Ser', 'T':'Thr', 'V':'Val', 'W':'Trp',
    'Y':'Tyr', 'Z':'Glx', 'X':'Xaa', 
    'U':'Sec', 'J':'Xle', 'O':'Pyl'
    }


standard_three_to_one = utils.invert_dict(one_to_three)

extended_three_to_one= {
'2as':'D', '3ah':'H', '5hp':'E', 'Acl':'R', 'Agm':'R', 'Aib':'A', 'Ala':'A', 'Alm':'A', 'Alo':'T', 'Aly':'K', 'Arg':'R', 'Arm':'R', 'Asa':'D', 'Asb':'D', 'Ask':'D', 'Asl':'D', 'Asn':'N', 'Asp':'D', 'Asq':'D', 'Asx':'B', 'Aya':'A', 'Bcs':'C', 'Bhd':'D', 'Bmt':'T', 'Bnn':'A', 'Buc':'C', 'Bug':'L', 'C5c':'C', 'C6c':'C', 'Ccs':'C', 'Cea':'C', 'Cgu':'E', 'Chg':'A', 'Cle':'L', 'Cme':'C', 'Csd':'A', 'Cso':'C', 'Csp':'C', 'Css':'C', 'Csw':'C', 'Csx':'C', 'Cxm':'M', 'Cy1':'C', 'Cy3':'C', 'Cyg':'C', 'Cym':'C', 'Cyq':'C', 'Cys':'C', 'Dah':'F', 'Dal':'A', 'Dar':'R', 'Das':'D', 'Dcy':'C', 'Dgl':'E', 'Dgn':'Q', 'Dha':'A', 'Dhi':'H', 'Dil':'I', 'Div':'V', 'Dle':'L', 'Dly':'K', 'Dnp':'A', 'Dpn':'F', 'Dpr':'P', 'Dsn':'S', 'Dsp':'D', 'Dth':'T', 'Dtr':'W', 'Dty':'Y', 'Dva':'V', 'Efc':'C', 'Fla':'A', 'Fme':'M', 'Ggl':'E', 'Gl3':'G', 'Gln':'Q', 'Glu':'E', 'Glx':'Z', 'Gly':'G', 'Glz':'G', 'Gma':'E', 'Gsc':'G', 'Hac':'A', 'Har':'R', 'Hic':'H', 'Hip':'H', 'His':'H', 'Hmr':'R', 'Hpq':'F', 'Htr':'W', 'Hyp':'P', 'Iil':'I', 'Ile':'I', 'Iyr':'Y', 'Kcx':'K', 'Leu':'L', 'Llp':'K', 'Lly':'K', 'Ltr':'W', 'Lym':'K', 'Lys':'K', 'Lyz':'K', 'Maa':'A', 'Men':'N', 'Met':'M', 'Mhs':'H', 'Mis':'S', 'Mle':'L', 'Mpq':'G', 'Msa':'G', 'Mse':'M', 'Mva':'V', 'Nem':'H', 'Nep':'H', 'Nle':'L', 'Nln':'L', 'Nlp':'L', 'Nmc':'G', 'Oas':'S', 'Ocs':'C', 'Omt':'M', 'Paq':'Y', 'Pca':'E', 'Pec':'C', 'Phe':'F', 'Phi':'F', 'Phl':'F', 'Pr3':'C', 'Pro':'P', 'Prr':'A', 'Ptr':'Y', 'Pyl':'O', 'Sac':'S', 'Sar':'G', 'Sch':'C', 'Scs':'C', 'Scy':'C', 'Sec':'U', 'Sel':'U', 'Sep':'S', 'Ser':'S', 'Set':'S', 'Shc':'C', 'Shr':'K', 'Smc':'C', 'Soc':'C', 'Sty':'Y', 'Sva':'S', 'Ter':'*', 'Thr':'T', 'Tih':'A', 'Tpl':'W', 'Tpo':'T', 'Tpq':'A', 'Trg':'K', 'Tro':'W', 'Trp':'W', 'Tyb':'Y', 'Tyq':'Y', 'Tyr':'Y', 'Tys':'Y', 'Tyy':'Y', 'Unk':'X', 'Val':'V', 'Xaa':'X', 'Xer':'X', 'Xle':'J'}
# Initial table is from the ASTRAL RAF release notes.
# added UNK
# Extra IUPAC: Xle, Xaa, Sec, Pyl
# The following have been seen in biopython code.
# Ter : '*'     Termination
# Sel : 'U'     A typo for Sec, selenocysteine? 
# Xer : 'X'     Another alternative for unknown?


amino_acid_names = {
    'A'	: 'alanine',	
    'M'	: 'methionine',  
    'C'	: 'cysteine',
    'N'	: 'asparagine',
    'D'	: 'aspartic acid',
    'P'	: 'proline',
    'E'	: 'glutamic acid',
    'Q'	: 'glutamine',
    'F'	: 'phenylalanine',
    'R'	: 'arginine',
    'G'	: 'glycine',	
    'S'	: 'serine',
    'H'	: 'histidine',	
    'T' : 'threonine',
    'I'	: 'isoleucine',	
    'V'	: 'valine',
    'K'	: 'lysine',
    'W'	: 'tryptophan', 
    'L'	: 'leucine',	
    'Y'	: 'tyrosine', 
    'B' : 'aspartic acid or asparagine',
    'J' : 'leucine or isoleucine',
    'X' : 'unknown',
    'Z' : 'glutamic acid or glutamine',
    'U' : 'selenocysteine',
    'O' : 'pyrrolysine',
    '*' : 'translation stop',
    '-' : 'gap'
    }

amino_acid_composition = dict(
    A = .082, R = .057, N = .044, D = .053, C = .017, 
    Q = .040, E = .062, G = .072, H = .022, I = .052,  
    L = .090, K = .057, M = .024, F =.039, P = .051, 
    S = .069, T = .058, W = .013, Y= .032, V =.066 )
      

kyte_doolittle_hydrophobicity = dict(
    A=1.8, R=-4.5, N=-3.5, D=-3.5,  C=2.5, 
    Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5,
    L=3.8, K=-3.9,  M=1.9,  F=2.8, P=-1.6,
    S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2 )


nucleotide_names = { 
    'A' : 'Adenosine',
    'C'	: 'Cytidine',
    'G'	: 'Guanine',
    'T'	: 'Thymidine',
    'U'	: 'Uracil',
    'R'	: 'G A (puRine)',
    'Y'	: 'T C (pYrimidine)',
    'K'	: 'G T (Ketone)',
    'M'	: 'A C (aMino group)',
    'S'	: 'G C (Strong interaction)',
    'W'	: 'A T (Weak interaction)',
    'B'	: 'G T C (not A) (B comes after A)',
    'D'	: 'G A T (not C) (D comes after C)',
    'H'	: 'A C T (not G) (H comes after G)',
    'V'	: 'G C A (not T, not U) (V comes after U)',
    'N' : 'A G C T (aNy)',
    '-' : 'gap', 
    }
    

# TODO: CHECK VALUES, UNITS    
amino_acid_accesible_surface_area = {
    'A' : 106.0,
    'C' : 135.0,
    'D' : 163.0,
    'E' : 194.0,
    'F' : 197.0,
    'G' : 84.0,
    'H' : 184.0,
    'I' : 169.0,
    'K' : 205.0,
    'L' : 164.0,
    'M' : 188.0,
    'N' : 157.0,
    'P' : 136.0,
    'Q' : 198.0,
    'R' : 248.0,
    'S' : 130.0,
    'T' : 142.0,
    'V' : 142.0,
    'W' : 227.0,
    'Y' : 222.0
    }
    
    
    
    
    
    
