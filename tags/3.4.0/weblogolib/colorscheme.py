
#  Copyright (c) 2003-2005 The Regents of the University of California.
#  Copyright (c) 2005 Gavin E. Crooks

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

""" Popular color codings for nucleic and amino acids. 

Classes:
    ColorScheme -- A color scheme
    ColorGroup  
    
    
Generic
    monochrome

Nucleotides
    nucleotide
    base pairing

Amino Acid
    hydrophobicity
    chemistry
    charge
    taylor

Status : Beta - Needs documentation.

"""
# Good online references include bioruby and the JalView alignment editor.
# Clamp, M., Cuff, J., Searle, S. M. and Barton, G. J. (2004), 
# "The Jalview Java Alignment Editor," Bioinformatics, 12, 426-7
# http://www.jalview.org


from corebio import seq
from .color import Color

class ColorScheme(object):
    """ A coloring of an alphabet.
    
    title : string            -- A human readable description
    defualt_color : Color           --
    groups : list of color groups 
    alphabet : string               -- The set of colored symbols
    color -- A map between a symbol and a Coloring
    

    """
    
    def __init__(self, 
                groups = [], 
                title = "", 
                description = "",
                default_color = "black", 
                alphabet = seq.generic_alphabet) :
        """  """
        self.title= title
        self.description = description
        self.default_color = Color.from_string(default_color)
        self.groups = groups
        self.alphabet = alphabet
            
        color = {}
        for cg in groups :
            for s in cg.symbols :
                color[s] = cg.color
                if s not in alphabet :
                    raise KeyError("Colored symbol does not exist in alphabet.")
        self._color = color

    def color(self, symbol) :
        if symbol in self._color :
            return self._color[symbol]
        return self.default_color
        
class ColorGroup(object) :
    """Associate a group of symbols with a color"""
    def __init__(self, symbols, color, description=None) :
        self.symbols = symbols              
        self.color =  Color.from_string(color)
        self.description = description


         
monochrome = ColorScheme([]) # This list intentionally left blank
               
# From makelogo
nucleotide = ColorScheme([
    ColorGroup("G", "orange"),
    ColorGroup("TU", "red"),
    ColorGroup("C",  "blue"),
    ColorGroup("A",  "green") 
    ]) 

base_pairing = ColorScheme([
    ColorGroup("TAU",  "darkorange", "Weak (2 Watson-Crick hydrogen bonds)"),
    ColorGroup("GC",    "blue", "Strong (3 Watson-Crick hydrogen bonds)")],
    )

# From Crooks2004c-Proteins-SeqStr.pdf
hydrophobicity = ColorScheme([
    ColorGroup( "RKDENQ",   "blue", "hydrophilic"),
    ColorGroup( "SGHTAP",   "green", "neutral"  ),
    ColorGroup( "YVMCLFIW", "black",  "hydrophobic") ],
    alphabet = seq.unambiguous_protein_alphabet
    )

# from makelogo
chemistry = ColorScheme([
  ColorGroup( "GSTYC",  "green",   "polar"),
  ColorGroup( "NQ",      "purple", "neutral"), 
  ColorGroup( "KRH",     "blue",   "basic"),
  ColorGroup( "DE",      "red",    "acidic"),
  ColorGroup("PAWFLIMV", "black",  "hydrophobic") ],
  alphabet = seq.unambiguous_protein_alphabet
  )   

charge = ColorScheme([
    ColorGroup("KRH", "blue", "Positive" ),
    ColorGroup( "DE", "red", "Negative") ],
    alphabet = seq.unambiguous_protein_alphabet
    )


taylor = ColorScheme([
    ColorGroup( 'A', '#CCFF00' ),
    ColorGroup( 'C', '#FFFF00' ),
    ColorGroup( 'D', '#FF0000'),
    ColorGroup( 'E', '#FF0066' ),
    ColorGroup( 'F', '#00FF66'),
    ColorGroup( 'G', '#FF9900'),
    ColorGroup( 'H', '#0066FF'),
    ColorGroup( 'I', '#66FF00'),
    ColorGroup( 'K', '#6600FF'),
    ColorGroup( 'L', '#33FF00'),
    ColorGroup( 'M', '#00FF00'),
    ColorGroup( 'N', '#CC00FF'),
    ColorGroup( 'P', '#FFCC00'),
    ColorGroup( 'Q', '#FF00CC'),
    ColorGroup( 'R', '#0000FF'),
    ColorGroup( 'S', '#FF3300'),
    ColorGroup( 'T', '#FF6600'),
    ColorGroup( 'V', '#99FF00'),
    ColorGroup( 'W', '#00CCFF'),
    ColorGroup( 'Y', '#00FFCC')],
    title = "Taylor",
    description = "W. Taylor, Protein Engineering, Vol 10 , 743-746 (1997)",
    alphabet = seq.unambiguous_protein_alphabet
    )
    


