
# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# Modifications Copyright 2004/2005 James Casbon. 
# Copyright 2005 by Regents of the University of California. All rights reserved
#   (Major rewrite for conformance to corebio. Gavin Crooks)
#
# This code is derived from the Biopython distribution and is governed by it's
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""ASTRAL:  Compendium for Sequence and Structure Analysis. 

The ASTRAL compendium provides databases and tools useful for analyzing protein structures and their sequences.  It is partially derived from, and augments the SCOP: Structural Classification of Proteins database. Most of the resources depend upon the coordinate files maintained and distributed by the Protein Data Bank.

Ref:
    http://astral.berkeley.edu/ 

* Classes :
    - Raf       -- A file of ASTRAL RAF (Rapid Access Format) Sequence Maps.
    - RafSeqMap -- A sequence map, a RAF record.
    - Res       -- A single residue mapping from a RAF record.
    
* Functions :
    - parse_domain  -- Convert an ASTRAL fasta header string into a Scop domain.
    - normalize_letters -- Normalize RAF amino acid codes.

"""

import re
from copy import copy

from corebio.db.scop import Domain, Residues
from corebio.data import extended_three_to_one as to_one_letter_code
from corebio.utils import FileIndex

__all__ = ('astral_evalues', 'astral_percent_identities', 
            'astral_evalues_filenames', 'normalize_letters', 'parse_domain', 
            'Raf', 'RafSeqMap', 'Res')

# Percentage identity filtered ASTRAL SCOP genetic domain sequence subset
astral_percent_identities = [10,20,25,30,35,40,50,70,90,95,100]

# E-value filtered ASTRAL SCOP genetic domain sequence subsets, based on PDB SEQRES records.
astral_evalues = [10, 5, 1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 1e-4, 1e-5, 1e-10, 1e-15,1e-20, 1e-25, 1e-50]

# A map between evalues and astral filename suffixes.
astral_evalues_filenames = { 
    10: 'e+1', 5: 'e+0,7', 1: 'e+0', 0.5: 'e-0,3', 0.1: 'e-1',
    0.05: 'e-1,3', 0.01: 'e-2', 0.005: 'e-2,3', 0.001: 'e-3',
    1e-4: 'e-4',  1e-5: 'e-5', 1e-10: 'e-10', 1e-15: 'e-15',
    1e-20: 'e-20', 1e-25: 'e-25', 1e-50: 'e-50' }



def normalize_letters(one_letter_code) :
    """Convert RAF one-letter amino acid codes into IUPAC standard codes.
    Letters are uppercased, and "." ("Unknown") is converted to "X".
    """
    if one_letter_code == '.' :
        return 'X'
    else :
        return one_letter_code.upper()
     
_domain_re = re.compile(r">?([\w_\.]*)\s+([\w\.]*)\s+\(([^)]*)\) (.*)")
def parse_domain(str) :
    """Convert an ASTRAL fasta header string into a SCOP domain.

    An ASTRAL (http://astral.stanford.edu/) header contains a concise
    description of a SCOP domain. A very similar format is used when a
    Domain object is converted into a string.  The Domain returned by this
    method contains most of the SCOP information, but it will not be located
    within the SCOP hierarchy (i.e. the parent node will be None). The
    description is composed of the SCOP protein and species descriptions.

    A typical ASTRAL header looks like --
    >d1tpt_1 a.46.2.1 (1-70) Thymidine phosphorylase {Escherichia coli}
    """

    m = _domain_re.match(str)
    if (not m) : raise ValueError("Domain: "+ str)

    dom = Domain()
    dom.sid = m.group(1)
    dom.sccs = m.group(2)
    dom.residues = Residues(m.group(3))
    if not dom.residues.pdbid :
        dom.residues.pdbid= dom.sid[1:5]
    dom.description = m.group(4).strip()

    return dom

    
class Raf(FileIndex) :
    """ASTRAL RAF (Rapid Access Format) Sequence Maps.

    The ASTRAL RAF Sequence Maps record the relationship between the PDB SEQRES
    records (representing the sequence of the molecule used in an experiment) 
    and the ATOM records (representing the atoms experimentally observed). 

    This data is derived from the Protein Data Bank CIF files. Known errors in 
    the CIF files are corrected manually, with the original PDB file serving as
    the final arbiter in case of discrepancies. 

    Residues are referenced by residue ID. This consists of the PDB residue
    sequence number (up to 4 digits) and an optional PDB insertion code (an
    ascii alphabetic character, a-z, A-Z). e.g. "1", "10A", "1010b", "-1"

    See "ASTRAL RAF Sequence Maps":http://astral.stanford.edu/raf.html

    The RAF file itself is about 50 MB. Each line consists of a sequence map of
    a different protein chain. This index provides rapid, random
    access of RAF records without having to load the entire file into memory.

    This class does not load the entire RAF file into memory. Instead, it
    reads the file once, noting the location and content of each RafSeqMap.
    The index key is a concatenation of the  PDB ID and chain ID. e.g
    "2drcA", "155c_". RAF uses an underscore to indicate blank
    chain IDs. Custom maps of subsequences or spanning multiple chains can
    be constructed with the get_seqmap method. 
    
    """
    def __init__(self, raf_file) :
        def linekey(line) :
            if not line or len(line)<5 or line.isspace() or line[0]=='#':
                return None
            return line[0:5]
        def parser( f) : return RafSeqMap(f.readline())
        
        FileIndex.__init__(self, raf_file, linekey, parser)
        

    def get_seqmap(self, residues) :
        """Get the sequence map for a collection of residues.

        residues -- A SCOP style description of a collection of residues from a
                    PDB strucure, (e.g. '(1bba A:10-20,B:)'), as a string or a
                    scop.Residues instance.
        """
        if type(residues)== str :
            residues = Residues(residues)

        pdbid  = residues.pdbid
        frags = residues.fragments
        if not frags: frags =(('_','',''),) # All residues of unnamed chain

        seqMap = None
        for frag in frags :
            chainid = frag[0]
            if chainid=='' or chainid=='-' or chainid==' ' or chainid=='_':
                chainid = '_'
            sid = pdbid + chainid
            
            sm = self[sid]
            
            # Cut out fragment of interest
            start = 0
            end = len(sm.res)
            if frag[1] : start = int(sm.index(frag[1], chainid))
            if frag[2] : end = int(sm.index(frag[2], chainid)+1)
            
            sm = sm[start:end]

            if seqMap is None :
                seqMap = sm
            else :
                seqMap += sm
                            
        return seqMap
 # End Raf

class RafSeqMap(object) :
    """ASTRAL RAF (Rapid Access Format) Sequence Maps.

    RafSeqMap is a list like object; you can find the location of particular
    residues with index(), slice this RafSeqMap into fragments, and glue 
    fragments back together with extend().

    - pdbid -- The PDB 4 character ID
    - pdb_datestamp -- From the PDB file
    - version -- The RAF format version. e.g. 0.01
    - flags -- RAF flags. (See release notes for more information.)
    - res -- A list of Res objects, one for each residue in this sequence map
    """

    def __init__(self, raf_record=None) :
        """Parses a RAF record into a RafSeqMap object."""
        
        self.pdbid = ''
        self.pdb_datestamp = ''
        self.version = ''
        self.flags = ''
        self.res = []

        if not raf_record : return
        
        header_len = 38
        line = raf_record.rstrip()  # no trailing whitespace        

        if len(line)<header_len: 
            raise ValueError("Incomplete header: "+line)

        self.pdbid = line[0:4]
        chainid = line[4:5]
        
        self.version = line[6:10]

        # Raf format versions 0.01 and 0.02 are identical for practical purposes
        if(self.version != "0.01" and  self.version !="0.02") :
            raise ValueError("Incompatible RAF version: "+self.version) 

        self.pdb_datestamp = line[14:20]
        self.flags = line[21:27]

        for i in range(header_len, len(line), 7) :
            f = line[i : i+7]
            if len(f)!=7:
                raise ValueError("Corrupt Field: ("+f+")" )
            r = Res()
            r.chainid = chainid
            r.resid =  f[0:5].strip()
            r.atom = normalize_letters(f[5:6])
            r.seqres = normalize_letters(f[6:7])

            self.res.append(r)
    # end __init__

    #@staticmethod
    def records(raf_file) :
        """Iterates over a Raf file, generating RafSeqMaps """
        for line in raf_file:
            if line[0] =='#':  continue  # A comment 
            if line.isspace() : continue
            yield RafSeqMap(line)        
    records = staticmethod(records)      
        
    def index(self, resid, chainid="_") :
        for i in range(0, len(self.res)) :
            if self.res[i].resid == resid and self.res[i].chainid == chainid :
                return i
        raise KeyError("No such residue "+chainid+resid)

    def __getslice__(self, i, j) :
        s = copy(self)
        s.res = s.res[i:j]
        return s

    def append(self, res) :
        """Append another Res object onto the list of residue mappings."""
        self.res.append(res)

    def extend(self, other) :
        """Append another RafSeqMap onto the end of self.

        Both RafSeqMaps must have the same PDB ID, PDB datestamp and
        RAF version.  The RAF flags are erased if they are inconsistent. This
        may happen when fragments are taken from different chains.
        """
        if not isinstance(other, RafSeqMap):
            raise TypeError("Can only extend a RafSeqMap with a RafSeqMap.")
        if self.pdbid != other.pdbid :
            raise TypeError("Cannot add fragments from different proteins.")
        if self.version != other.version :
            raise TypeError("Incompatible rafs.")
        if self.pdb_datestamp != other.pdb_datestamp :
            raise TypeError("Different pdb dates!")
        if self.flags != other.flags :
            self.flags = ''
        self.res += other.res

    def __iadd__(self, other) :
        self.extend(other)
        return self

    def __add__(self, other) :
        s = copy(self)
        s.extend(other)
        return s

    def extract_atoms(self, pdb_handle, out_handle) :
        """Extract all relevant ATOM and HETATOM records from a PDB file.

        The PDB file is scanned for ATOM and HETATOM records. If the
        chain ID, residue ID (seqNum and iCode), and residue type match
        a residue in this sequence map, then the record is echoed to the
        output handle.

        This is typically used to find the coordinates of a domain, or other
        residue subset.

        pdb_file -- A handle to the relevant PDB file.
        out_file -- All output is written to this stream.
        """
        resSet = {}
        for r in self.res :
            if r.atom=='X' : # Unknown residue type
                continue
            chainid = r.chainid
            if chainid == '_':
                chainid = ' '
            resid = r.resid
            resSet[(chainid,resid)] = r

        resFound = {}
        for line in pdb_handle :
            if line.startswith("ATOM  ") or line.startswith("HETATM") :
                chainid = line[21:22]
                resid = line[22:27].strip()
                key = (chainid, resid)
                if key in resSet:
                    res = resSet[key]
                    atom_aa = res.atom
                    resName = line[17:20].capitilize()
                    if resName in to_one_letter_code :
                        if to_one_letter_code[resName] == atom_aa :
                            out_handle.write(line)
                            resFound[key] = res

        if len(resSet) != len(resFound) :
            raise RuntimeError('I could not find at least one ATOM or ' 
               'HETATM record for each and every residue in this sequence map.')

        
class Res(object) :
    """ A single residue mapping from a RAF record.

    - chainid -- A single character chain ID.
    - resid   -- The residue ID. 
    - atom    -- amino acid one-letter code from ATOM records. 
    - seqres  -- amino acid one-letter code from SEQRES records.
    """
    def __init__(self) :
        self.chainid = ''
        self.resid = ''
        self.atom = ''
        self.seqres = ''

   
