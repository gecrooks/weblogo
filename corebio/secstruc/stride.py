
#  Copyright (c) 2003 Gavin E. Crooks
#  Copyright (c) 2005 David D. Ding <dding@berkeley.edu>
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

"""Command and control of the program STRIDE: Protein secondary structure assignment from atomic coordinates.

This module provides an interface to STRIDE, a program used to recognize
secondary structural elements in proteins from their atomic coordinates.

Refs:
    - http://wolf.bi.umist.ac.uk/unix/stride.html
"""
from __future__ import absolute_import

from subprocess import * 

from ..db.astral import to_one_letter_code
from ..seq import Seq, protein_alphabet, Alphabet
from ..utils import stdrepr, find_command
from .._py3k import StringIO

# alphabet for stride secondary structure
stride_alphabet = Alphabet("HGIEBC12345678@&T")

# Dictionary for conversion between names and alphabet
stride_alphabet_names  = { 
    "H": "AlphaHelix",
    "G": "310Helix",
    "I": "PiHelix",
    "E": "Strand",
    "b": "Bridge",
    "B": "Bridge",
    "C": "Coil",
    "1": "TurnI",
    "2": "TurnI'",
    "3": "TurnII",
    "4": "TurnII'",
    "5": "TurnVIa",
    "6": "TurnVIb",
    "7": "TurnVIII",
    "8": "TurnIV",
    "@": "GammaClassic",
    "&": "GammaInv",
    "T": "Turn"
    }

class RunStride(object) :
    """Run a locally installed STRIDE executable."""
    def __init__(self, path=None) :
        command = find_command('stride', path=path) 
        self.command = command
    
    
    def process_pdb(self, pdb_filename) :
        """Process a protein structure file in PDB format through the STRIDE 
        program and return the raw output.
        """
        args = [self.command, pdb_filename]
        try :
            p = Popen(args, stdout=PIPE)
            (out,err) = p.communicate() 
        except OSError :
            raise RuntimeError("Cannot communicate with STRIDE.")  
        return out
    
    def record(self, pdb_filename) :
        """Process a protein structure file in PDB format through the STRIDE
        program and return the parsed output as a StrideRecord object.
        """
        data = self.process_pdb(pdb_filename) 
        return StrideRecord(StringIO(data))

#end class        

class StrideRecord(object) :
    """Representation of stride output."""
    __slots__ = ['pdbid', 'residues', '_res_dict']
    
    def __init__(self, stride_file) :
        """ Read and parse a STRIDE output file.
        
        args:
            - stride_file   : An open file handle
        attributes :
            - pdbid     : The PDB id.
            - res       : A list of Res objects, one per PDB resiude
        """     
        res =[]
        f=stride_file
        self.pdbid = f.readline()[75:79]
        for l in f:
            if l[0:3] =="ASG":
                res.append(StrideResidue(l)) 
                
        self.residues = res # A list of Res objects
        
        self._res_dict = None

    def total_area(self) :
        """ Return the solvent accessible area """
        area = 0
        for i in self.residues :
            area += i.solvent_acc_area
        return area
    
    def primary(self):
        """ Return the protein primary sequence as a Seq object."""
        return Seq(''.join([r.aa for r in self.residues]), protein_alphabet)
        
    def secondary(self):
        """Return the secondary structure of the protein as a Seq object"""
        return Seq(''.join([r.secstruc for r in self.residues]), stride_alphabet)
        
        
    def get_residue(self, chainid, resid) :
        """ Return the given residue """
        if not self._res_dict :
            d = {}
            for r in self.residues :
                d[ (r.chainid, r.resid)] = r
            self._res_dict =d
        
        return self._res_dict[(chainid, resid)]
        
    def __repr__(self): 
        return stdrepr(self)
# end class
    
    
class StrideResidue(object):                 
    """ Structural information of a single residue. An ASG line from a stride
        output file.
        
        Attributes :
         - chainid 
         - resid   
         - aa 
         - secstruc 
         - solvent_acc_area 
         - phi 
         - psi
    """
    __slots__ = ['chainid', 'resid', 'aa', 'secstruc', 'solvent_acc_area', 
                    'phi', 'psi', ]
    def __init__(self, res_line) :
        """ Eats a single 'ASG' line from a stride file, splits it up  
        into parts and return a Res object."""
            
        if (len(res_line)<70): 
            raise ValueError("Line not long enough")
        try: 
            self.chainid = res_line[9:10]
            # STRIDE converts blank chain ids into dashes. Undo.
            if self.chainid=="-" : self.chainid = " "
                
            # In rare cases STRIDE columns can be misaligned. Grab extra 
            # white space to compensate.
            self.resid = res_line[10:15].strip() 
            self.aa = to_one_letter_code[res_line[5:8].capitalize()]
            self.secstruc = res_line[24:25]
            self.solvent_acc_area = float(res_line[64:71]) 
            self.phi = float(res_line[42:49].strip())
            self.psi = float(res_line[52:59].strip())
        except FloatingPointError:
            raise FloatingPointError("Can't float phi, psi, or area")
        except KeyError:
            raise KeyError("Can't find three letter code in dictionary")
        except LookupError:
            raise LookupError("One of the values is out of index of res_line")
            
                
    # End __init__

    def __repr__(self): 
        return stdrepr(self)         
# End class










        
