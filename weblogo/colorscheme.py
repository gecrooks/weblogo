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
    SymbolColor
    IndexColor
    RefSeqColor

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


from typing import Sequence, List, Optional
from . import seq
from .color import Color
from .seq import Alphabet


# TODO: Make as abstract
class ColorRule(object):
    """
    Define an interface for coloring individual symbols based on their position
    and identity.  Subclasses should reimplement the symbol_color() method to
    return a Color object based on the given parameters.
    """

    def symbol_color(self, seq_index: int, symbol: str, rank: int) -> Optional[Color]:
        raise NotImplementedError   # pragma: no cover


class ColorScheme(ColorRule):
    """
    Specify which color each symbol in a sequence logo should be.

    A color scheme is primarily a container of color rules.  These rules would
    be along the lines of "hydrophobic residues are blue" or "indices 5-10 are
    red" or "the wildtype sequence is black".  When a color is requested for a
    particular symbol, each rule is consulted in turn until one provides a
    color.  If no rule provides a color, the given default color will be used.
    """

    def __init__(self,
                 rules: List[ColorRule] = [],
                 title: str = "",
                 description: str = "",
                 default_color: str = "black",
                 alphabet: Alphabet = seq.generic_alphabet) -> None:

        self.rules = rules
        self.title = title
        self.description = description
        self.default_color = Color.from_string(default_color)
        self.alphabet = alphabet

    def symbol_color(self, seq_index: int, symbol: str, rank: int) -> Color:
        if symbol not in self.alphabet:
            raise KeyError("Colored symbol '%s' does not exist in alphabet." % symbol)

        for rule in self.rules:
            color = rule.symbol_color(seq_index, symbol, rank)
            if color is not None:
                return color

        return self.default_color


class SymbolColor(ColorRule):
    """
    Represent the given set of symbols (e.g. "DEHKR" for charged residues) with
    a single color.
    """

    def __init__(self, symbols: str, color: str, description: str = None) -> None:
        self.symbols = symbols.upper()
        self.color = Color.from_string(color)
        self.description = description

    def symbol_color(self, seq_index: int, symbol: str, rank: int) -> Optional[Color]:
        if symbol.upper() in self.symbols:
            return self.color
        return None


class IndexColor(ColorRule):
    """
    Represent the given set of indices (e.g. range(10) for the first ten
    residues) with a single color.
    """

    def __init__(self, indices: Sequence[list], color: str, description: str = None) -> None:
        self.indices = indices
        self.color = Color.from_string(color)
        self.description = description

    def symbol_color(self, seq_index: int, symbol: str, rank: int) -> Optional[Color]:
        if seq_index in self.indices:
            return self.color
        return None


class RefSeqColor(ColorRule):
    """
    Color the given reference sequence in its own color, so you can easily see
    which positions match that sequence and which don't.
    """

    def __init__(self, ref_seq: str, color: str, description: str = None) -> None:
        self.ref_seq = ref_seq.upper()
        self.color = Color.from_string(color)
        self.description = description

    def symbol_color(self, seq_index: int, symbol: str, rank: int) -> Optional[Color]:
        if symbol.upper() == self.ref_seq[seq_index]:
            return self.color
        return None


monochrome = ColorScheme([])  # This list intentionally left blank

# From makelogo
nucleotide = ColorScheme(
        [
            SymbolColor("G", "orange"),
            SymbolColor("TU", "red"),
            SymbolColor("C", "blue"),
            SymbolColor("A", "green")
        ],
)

base_pairing = ColorScheme(
        [
            SymbolColor("TAU", "darkorange", "Weak (2 Watson-Crick hydrogen bonds)"),
            SymbolColor("GC", "blue", "Strong (3 Watson-Crick hydrogen bonds)")
        ],
)

# From Crooks2004c-Proteins-SeqStr.pdf
hydrophobicity = ColorScheme(
        [
            SymbolColor("RKDENQ", "blue", "hydrophilic"),
            SymbolColor("SGHTAP", "green", "neutral"),
            SymbolColor("YVMCLFIW", "black", "hydrophobic")
        ],
        alphabet=seq.unambiguous_protein_alphabet
)

# from makelogo
chemistry = ColorScheme(
        [
            SymbolColor("GSTYC", "green", "polar"),
            SymbolColor("NQ", "purple", "neutral"),
            SymbolColor("KRH", "blue", "basic"),
            SymbolColor("DE", "red", "acidic"),
            SymbolColor("PAWFLIMV", "black", "hydrophobic")
        ],
        alphabet=seq.unambiguous_protein_alphabet
)

charge = ColorScheme(
        [
            SymbolColor("KRH", "blue", "Positive"),
            SymbolColor("DE", "red", "Negative")
        ],
        alphabet=seq.unambiguous_protein_alphabet
)

taylor = ColorScheme(
        [
            SymbolColor('A', '#CCFF00'),
            SymbolColor('C', '#FFFF00'),
            SymbolColor('D', '#FF0000'),
            SymbolColor('E', '#FF0066'),
            SymbolColor('F', '#00FF66'),
            SymbolColor('G', '#FF9900'),
            SymbolColor('H', '#0066FF'),
            SymbolColor('I', '#66FF00'),
            SymbolColor('K', '#6600FF'),
            SymbolColor('L', '#33FF00'),
            SymbolColor('M', '#00FF00'),
            SymbolColor('N', '#CC00FF'),
            SymbolColor('P', '#FFCC00'),
            SymbolColor('Q', '#FF00CC'),
            SymbolColor('R', '#0000FF'),
            SymbolColor('S', '#FF3300'),
            SymbolColor('T', '#FF6600'),
            SymbolColor('V', '#99FF00'),
            SymbolColor('W', '#00CCFF'),
            SymbolColor('Y', '#00FFCC')
        ],
        title="Taylor",
        description="W. Taylor, Protein Engineering, Vol 10 , 743-746 (1997)",
        alphabet=seq.unambiguous_protein_alphabet
)
