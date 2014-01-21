
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
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
#  THE SOFTWARE.

"""Command and control of the program DSSP: Dictionary of protein secondary structure. The program DSSP defines secondary structure, geometrical features and solvent exposure of proteins, given atomic coordinates in Protein Data Bank (PDB) format.


See also :
- http://swift.cmbi.ru.nl/gv/dssp/

"""
from __future__ import absolute_import

from subprocess import * 

from ..seq import Seq, protein_alphabet, Alphabet
from ..utils import stdrepr, find_command
from .._py3k import StringIO

# alphabet for stride secondary structure
dssp_alphabet = Alphabet("HBEGITS ")

# Dictionary for conversion between alphabet and secondary structure names
dssp_alphabet_names  = {
    'H' : 'alpha helix',
    'B' : 'residue in isolated beta-bridge',
    'E' : 'extended strand, participates in beta ladder',
    'G' : '3-helix (3/10 helix)',
    'I' : '5 helix (pi helix)',
    'T' : 'hydrogen bonded turn',
    'S' : 'bend',
    ' ' : 'loop or irregular', 
    }

_dssp_header = "==== Secondary Structure Definition by the program DSSP, updated CMBI version by ElmK / April 1,2000"


class RunDssp(object) :
    """Run a locally installed DSSP executable."""
    def __init__(self, path=None) :
        command = find_command('dsspcmbi', path=path) 
        self.command = command
    
    def version(self) :
        args = [self.command, '-V']
        try :
            p = Popen(args, stdout=PIPE)
            (out,err) = p.communicate() 
        except OSError :
            raise RuntimeError("Cannot communicate with DSSP.")  
        return out.strip()
    
    def process_pdb(self, pdb_file) :
        """Process a protein structure fle in PDB format through the DSSP 
        program and return the raw output.
        """
        args = [self.command, '--']
        try :
            p = Popen(args, stdin= pdb_file, stdout=PIPE)
            (out,err) = p.communicate() 
        except OSError :
            raise RuntimeError("Cannot communicate with DSSP.")  
        return out
    
    def record(self, pdb_file) :
        """Process a protein structure fle in PDB format through the DSSP
        program and return the parsed output as a DsspRecord object.
        """
        data = self.process_pdb(pdb_file) 
        return DsspRecord(StringIO(data))

#end class        


class DsspRecord(object) :
    """Representation of DSSP output."""
    __slots__ = ['pdbid', 'residues', '_res_dict']
    
    def __init__(self, dssp_file) :
        """Construct a DsspRecord from a DSSP output file.
        
        args:
            - dssp_file  : An open file handle
        attributes :
            - pdbid     : The PDB id.
            - residues  : A list of DsspResidue objects, one per PDB resiude
        """     
        res =[]
        lines = iter(dssp_file)

        header = next(lines)
        if not header.startswith(_dssp_header) :
            raise ValueError("Unrecognized file type: "+ header)
        line = next(lines)
        line = next(lines)
        self.pdbid = line[62:66].lower()

        for line in lines:
            if line.startswith('  #  RESIDUE') : break
        
        for line in lines :
            if line and line[13] != '!':
                res.append( DsspResidue(line) )
        self.residues = res 
                     
    def total_area(self) :
        """Return the total solvent accessible area """
        area = 0
        for r in self.residues :
            area += r.solvent_acc_area
        return area
    
    def primary(self):
        """ Return the protein primary sequence as a Seq object."""
        return Seq(''.join([r.aa for r in self.residues]), protein_alphabet)
        
    def secondary(self):
        """Return the secondary structure of the protein as a Seq object"""
        return Seq(''.join([r.secstruc for r in self.residues]), dssp_alphabet)
        
        
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



class DsspResidue(object):
    """Structural information for a single protein residue derived from a
    DSSP output file.

    Attributes :
    - num
        Sequential residue number
    - chainid 
    - resid   
    - aa 
        one letter amino acid code,
    - secstruc
        compromise summary of secondary structure, intended to approximate
        crystallographers' intuition 
    - secstruc_details
        Seven characters, 3-turns/helix, 4-turns/helix, 5-turns/helix, 
        geometrical bend, chirality, beta bridge label, beta bridge label
    - BP1, BP2
        residue number of first and second bridge partner followed by one letter
        sheet label        
    - solvent_acc_area
        residue water exposed surface in Angstrom**2 
    - hbond_NHtoO
    - hbond_OtoNH
    - tco 
        cosine of angle between C=O of residue I and C=O of residue I-1. For
        alpha-helices, TCO is near +1, for beta-sheets TCO is near -1. Not used
        for structure definition.
    - kappa
        virtual bond angle (bend angle) defined by the three C-alpha atoms of
        residues I-2,I,I+2. Used to define bend
    - alpha
        virtual torsion angle (dihedral angle) defined by the four C-alpha atoms
        of residues I-1,I,I+1,I+2. Used to define chirality 
    - phi, psi    
        IUPAC peptide backbone torsion angles
    - coord
         C-alpha atom coordinates    
    """
    
    __slots__ = ['num', 'chainid', 'resid', 'aa', 'secstruc',
                'secstruc_details',
                'BP1', 'BP2', 'solvent_acc_area', 'hbond_NHtoO', 'hbond_OtoNH',
                'tco', 'kappa', 'alpha', 'phi', 'psi', 'coord']
                
    def __init__(self, line) :  
        """Parse a single line of a DSSP file."""            
                
        try:
            self.num      = int(line[0:5].strip())
            self.chainid  = line[11]
            self.resid    = line[6:11].strip()
            self.aa         = line[13]
            self.secstruc  = line[16]
            self.secstruc_details = line[18:24]
            self.BP1= line[26:30].strip()
            self.BP2 = line[30:34].strip()
            self.solvent_acc_area = float(line[34:38].strip() )

            self.hbond_NHtoO = ((int(line[39:45]),float(line[46:50])), 
                                    ( int(line[61:67]),float(line[68:72]) ) )
            self.hbond_OtoNH = ((int(line[50:56]),float(line[57:61])), 
                                    ( int(line[72:78]),float(line[79:83]) ) )
            self.tco  = float(line[85:91].strip() )
            self.kappa= float(line[91:97 ].strip() )
            self.alpha= float(line[97:103].strip() )
            self.phi= float(line[103:109].strip() )
            self.psi= float(line[109:115].strip() )
            x = float(line[115:122].strip() )
            y= float(line[122:129].strip() )
            z= float(line[129:236].strip() )
            self.coord = (x,y,z)              
        except FloatingPointError:
            raise ValueError("Cannot parse line")

        # Additional validation            
        if self.secstruc not in dssp_alphabet : 
            raise ValueError("Unknown secondary structure code")
                
    # End __init__

    def __repr__(self): 
        return stdrepr(self)
# End class
                
                
                
