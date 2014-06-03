
# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# Modifications Copyright 2004/2005 James Casbon. 
# Copyright 2005 by Regents of the University of California. All rights Reserved.
#   (Major rewrite for conformance to corebio. Gavin Crooks)
#
# This code is derived from the Biopython distribution and is governed by it's
# license.  Please see the LICENSE file that should have been included
# as part of this package.


""" SCOP: Structural Classification of Proteins.

The SCOP database aims to provide a manually constructed classification of
all know protein structures into a hierarchy, the main levels of which
are family, superfamily and fold.

* SCOP: http://scop.mrc-lmb.cam.ac.uk/scop/
* Introduction: http://scop.mrc-lmb.cam.ac.uk/scop/intro.html
* SCOP parsable files: http://scop.mrc-lmb.cam.ac.uk/scop/parse/

The Scop object in this module represents the entire SCOP classification. It
can be built from the three SCOP parsable files (see DesRecord, HieRecord and
ClaRecord), modified if so desired, and converted back to the same file formats.
A single SCOP domain (represented by the Domain class) can be obtained from
SCOP using the domain's SCOP identifier (sid).

Classes:
 - Scop     -- The entire SCOP hierarchy.
 - Node     -- A node in the SCOP hierarchy.
 - Domain   -- A SCOP domain.
 - Residues -- A collection of residues from a PDB structure.
 - HieRecord -- Handle the SCOP HIErarchy files.
 - DesRecord -- Handle the SCOP DEScription file.
 - ClaRecord -- Handle the SCOP CLAssification file.


nodeCodeDict  -- A mapping between known 2 letter node codes and a longer
                  description. The known node types are 'cl' (class), 'cf'
                  (fold), 'sf' (superfamily), 'fa' (family), 'dm' (domain), 
                  'sp' (species), 'px' (domain). Additional node types may
                  be added in the future.
"""
from __future__ import print_function

import os, re

from .._py3k import iteritems, cmp

nodeCodeDict = { 'cl':'class', 'cf':'fold', 'sf':'superfamily',
                 'fa':'family', 'dm':'protein', 'sp':'species', 'px':'domain'}


_nodetype_to_code = dict((v, k) for k, v in iteritems(nodeCodeDict))


nodeCodeOrder = [ 'ro', 'cl', 'cf', 'sf', 'fa', 'dm', 'sp', 'px' ] 


def cmp_sccs(sccs1, sccs2) :
    """Order SCOP concise classification strings (sccs).

    a.4.5.1 < a.4.5.11 < b.1.1.1 

    A sccs (e.g. a.4.5.11) compactly represents a domain's classification.
    The letter represents the class, and the numbers are the fold,
    superfamily, and family, respectively.

    """
    s1 = sccs1.split(".")
    s2 = sccs2.split(".")
    if s1[0] != s2[0]:
        return cmp(s1[0], s2[0])

    s1 = tuple(map(int, s1[1:]))
    s2 = tuple(map(int, s2[1:]))
    return cmp(s1, s2)


def sccs_relation(sccs1, sccs2) :
    """Are two SCOP domains related? Returns +1 if homologous , -1 if not homologous, 
    and 0 if ambiguous.

    Protein domains are classified into a hierarchy of class, fold, superfamily
    and family. Homologous domains are placed into the  same superfamily, 
    whereas domains belonging to different classes or folds are considered 
    unrelated. The evolutionary relationship of domains classified in the same
    fold but different superfamilies is ambiguous. (Green and Brenner, 2002).           

     - sccs1, sccs2 : Two SCOP concise classification string (e.g. a.4.5.11) which represent 
    a domain's classification. The letter represents the class, and the numbers are the fold,
    superfamily, and family, respectively.
    """

    s1 = sccs1.split('.')
    s2 = sccs2.split('.')

    if s1[0] != s2[0] :
        #  Different classes: not homologs
        return -1
    elif s1[1] != s2[1] :
        #  Different folds: not homologs
        return -1
    elif s1[2] == s2[2] :
        # Same class, fold and superfamily: homologs
        return +1
    else :
        # Same class and fold, but different superfamilies: Ambiguous
        return 0



def _open_scop_file(scop_dir_path, version, filetype) :
    filename = "dir.%s.scop.txt_%s" % (filetype,version)
    afile = open(os.path.join( scop_dir_path, filename))
    return afile


class Scop(object):
    """The entire SCOP hierarchy.

    root             -- The root node of the hierarchy 
    domains          -- A list of all domains
    nodes_by_sid     -- A dictionary of nodes indexed by SCOP identifier 
                        (e.g. 'd1hbia_')
    domains_by_sunid -- A dictionary of domains indexed by SCOP unique
                        identifiers (e.g. 14996)
    """
    def __init__(self):
        """ An empty SCOP object.
        
        See also Scop.parse() and Scop.parse_files()
        """
        self.root = None
        self.domains = []
        self.nodes_by_sunid = dict()
        self.domains_by_sid = dict()

    @classmethod
    def parse(cls, dir_path, version) :
        """Build the SCOP hierarchy from the SCOP parsable files.
            
         - dir_path -- A directory that contains the SCOP files
         - version  -- The SCOP version (as a string)
         
        The SCOP files are named dir.XXX.scop.txt_VERSION, where XXX
        is 'cla', 'des' or 'hie'.
        """
        cla_file = None
        des_file = None
        hie_file = None
        try :
            cla_file = _open_scop_file( dir_path, version, 'cla')
            des_file = _open_scop_file( dir_path, version, 'des')
            hie_file = _open_scop_file( dir_path, version, 'hie')
            scop = cls.parse_files(cla_file, des_file, hie_file) 
        finally :
            # If we opened the files, we close the files
            if cla_file : cla_file.close()
            if des_file : des_file.close()
            if hie_file : hie_file.close()
        
        return scop

    @classmethod
    def parse_files(cls, cla_file, des_file, hie_file):
        """Build the SCOP hierarchy from the SCOP parsable files.
        
            - cla_file -- the CLA clasification file
            - des_file -- the DES description file
            - hie_file -- the HIE hierarchy file
        """

        self = cls()

        sunidDict = {}

        root = Node()
        domains = []
        root.sunid=0
        root.type='ro'
        sunidDict[root.sunid] = root

        root.description = 'SCOP Root'

        # Build the rest of the nodes using the DES file
        for rec in DesRecord.records(des_file):
            if rec.nodetype =='px' :
                n = Domain()
                n.sid = rec.name
                domains.append(n)
            else : 
                n = Node()
            n.sunid = rec.sunid
            n.type = rec.nodetype
            n.sccs = rec.sccs
            n.description = rec.description
                
            sunidDict[n.sunid] = n
 
        # Glue all of the Nodes together using the HIE file
        for rec in HieRecord.records(hie_file):
            if not rec.sunid in sunidDict :
                # print(rec.sunid)
                raise ValueError("Incomplete data?")

                    
            n = sunidDict[rec.sunid]
            if rec.parent !='': # Not root node
                if not rec.parent in sunidDict :
                    raise ValueError("Incomplete data?")
                n.parent = sunidDict[rec.parent]
                
            for c in rec.children:
                if not c in sunidDict :
                    raise ValueError("Incomplete data?")
                n.children.append(sunidDict[c])

                        
        # Fill in the gaps with information from the CLA file
        sidDict = {}
        for rec in ClaRecord.records(cla_file):
            n = sunidDict[rec.sunid]
            assert n.sccs == rec.sccs
            assert n.sid == rec.sid
            n.residues = rec.residues
            sidDict[n.sid] = n

        # Clean up
        self.root = root
        self.nodes_by_sunid = sunidDict
        self.domains_by_sid = sidDict
        self.domains = tuple(domains)

        return self


    def write_hie(self, stream) :
        """Build an HIE SCOP parsable file from this object"""
        # We order nodes to ease comparison with original file
        nodes = sorted(self.nodes_by_sunid.values(), key=lambda n: n.sunid)
        for n in nodes :
            stream.write(str(n.to_hie_record()))


    def write_des(self, stream) :
        """Build a DES SCOP parsable file from this object""" 
        # Origional SCOP file is not ordered?
        nodes = sorted(self.nodes_by_sunid.values(), key=lambda n: n.sunid)
        for n in nodes :
            if n != self.root :
                stream.write(str(n.to_des_record()))


    def write_cla(self, stream) :
        """Build a CLA SCOP parsable file from this object"""                
        # We order nodes to ease comparison with original file
        nodes = sorted(self.domains_by_sid.values(), key=lambda n: n.sunid)
        for n in nodes :
            stream.write(str(n.to_cla_record()))
# End Scop
     

  
class Node(object) :
    """ A node in the SCOP hierarchy

    sunid  -- SCOP unique identifiers. e.g. '14986'
    parent -- The parent node
    children -- A list of child nodes
    sccs     -- SCOP concise classification string. e.g. 'a.1.1.2'
    type     -- A 2 letter node type code. e.g. 'px' for domains
    description -- 
        
    """
    def __init__(self) :
        """A new, uninitialized SCOP node."""
        self.sunid=''    
        self.parent = None
        self.children=[]
        self.sccs = ''   
        self.type =''    
        self.description =''

    def __str__(self) :
        s = []
        s.append(str(self.sunid))
        s.append(self.sccs)
        s.append(self.type)
        s.append(self.description)

        return " ".join(s)

    def to_hie_record(self):
        """Return an Hie.Record"""
        rec = HieRecord()
        rec.sunid = str(self.sunid)
        if self.parent : # Not root node
            rec.parent = str(self.parent.sunid)
        else:
            rec.parent = '-'
        for c in self.children :
            rec.children.append(str(c.sunid))
        return rec
    
    def to_des_record(self):
        """Return a Des.Record"""        
        rec = DesRecord()
        rec.sunid = str(self.sunid)
        rec.nodetype = self.type
        rec.sccs = self.sccs
        rec.description = self.description
        return rec
        
    def descendents( self, node_type) :
        """ Return a list of all descendant nodes of the given type. Node type
        can be a two letter code or longer description. e.g. 'fa' or 'family'
        """
        if node_type in _nodetype_to_code:
            node_type = _nodetype_to_code[node_type]
            
        nodes = [self]

        while nodes[0].type != node_type:
            if nodes[0].type == 'px' : 
                return [] # Fell of the bottom of the hierarchy
            child_list = []
            for n in nodes:
                for child in n.children:
                    child_list.append( child )
                nodes = child_list
                
        return nodes
                    

    def ascendent( self, node_type) :
        """ Return the ancestor node of the given type, or None. Node type can 
        be a two letter code or longer description. e.g. 'fa' or 'family'
        """
        if node_type in _nodetype_to_code :
            node_type = _nodetype_to_code[node_type]

        n = self
        if n.type == node_type: return None
        while n.type != node_type:
            if n.type == 'ro': 
                return None # Fell of the top of the hierarchy
            n = n.parent            
                
        return n
# End Node                                                            
    

class Domain(Node) :
    """ A SCOP domain. A leaf node in the Scop hierarchy.

    - sid      -- The SCOP domain identifier. e.g. 'd5hbib_'
    - residues -- A Residue object. It defines the collection
                  of PDB atoms that make up this domain.
    """
    def __init__(self) :
        Node.__init__(self)
        self.sid = ''         
        self.residues = None

    def __str__(self) :
        s = []
        s.append(self.sid)
        s.append(self.sccs)
        s.append("("+str(self.residues)+")")

        if not self.parent :
            s.append(self.description)
        else :
            sp = self.parent
            dm = sp.parent
            s.append(dm.description)
            s.append("{"+sp.description+"}")

        return " ".join(s)


    def relation(self, other) :
        """BETA"""
        if self.sid == other.sid : return 0
        return sccs_relation(self.sccs, other.sccs)
        

    def to_des_record(self):
        """Return a des.Record"""
        rec = Node.to_des_record(self)
        rec.name = self.sid
        return rec

    def to_cla_record(self) :
        """Return a cla.Record"""        
        rec = ClaRecord()
        rec.sid = self.sid
        rec.residues = self.residues
        rec.sccs = self.sccs
        rec.sunid = self.sunid
        
        n = self
        while n.sunid != 0: # Not root node
            rec.hierarchy.append( (n.type, str(n.sunid)) )
            n = n.parent

        rec.hierarchy.reverse()
       
        return rec
# End Domain



class DesRecord(object):
    """ Handle the SCOP DEScription file.

    The file format is described in the SCOP
    "release notes.":http://scop.berkeley.edu/release-notes-1.55.html 
    The latest DES file can be found
    "elsewhere at SCOP.":http://scop.mrc-lmb.cam.ac.uk/scop/parse/
  
    The DES file consists of one DES record per line. Each record 
    holds information for one node in the SCOP hierarchy, and consists
    of 5 tab-delimited fields,
    sunid, node type, sccs, node name, node description.

    For example ::
    
    21953   px      b.1.2.1 d1dan.1 1dan T:,U:91-106
    48724   cl      b       -       All beta proteins
    48725   cf      b.1     -       Immunoglobulin-like beta-sandwich
    49265   sf      b.1.2   -       Fibronectin type III
    49266   fa      b.1.2.1 -       Fibronectin type III

   
    - sunid       -- SCOP unique identifiers
    - nodetype    -- One of 'cl' (class), 'cf' (fold), 'sf' (superfamily),
                   'fa' (family), 'dm' (protein), 'sp' (species),
                   'px' (domain). Additional node types may be added.
    - sccs        -- SCOP concise classification strings. e.g. b.1.2.1
    - name        -- The SCOP ID (sid) for domains (e.g. d1anu1),
                   currently empty for other node types
    - description --  e.g. "All beta proteins","Fibronectin type III", 
    """
    def __init__(self, record=None):
        
        if not record :
            self.sunid = ''
            self.nodetype = ''
            self.sccs = ''
            self.name = ''
            self.description =''
        else :
            entry = record.rstrip()  # no trailing whitespace
            columns = entry.split("\t")  # separate the tab-delimited cols
            if len(columns) != 5:
                raise ValueError("I don't understand the format of %s" % entry)
        
            self.sunid, self.nodetype, self.sccs, self.name, self.description \
                = columns
            if self.name == '-' : self.name =''
            self.sunid = int(self.sunid)
        
    def __str__(self):
        s = []
        s.append(self.sunid)
        s.append(self.nodetype)        
        s.append(self.sccs)        
        if self.name :
            s.append(self.name)
        else :
            s.append("-")
        s.append(self.description)        
        return "\t".join(map(str,s)) + "\n"

    @staticmethod
    def records(des_file):
        """Iterates over a DES file, generating DesRecords """
        for line in des_file:
            if line[0] =='#':  continue  # A comment 
            if line.isspace() : continue
            yield DesRecord(line)

# End DesRecord        


class HieRecord(object):
    """Handle the SCOP HIErarchy files, which describe the SCOP hierarchy in
    terms of SCOP unique identifiers (sunid).

    The file format is described in the SCOP
    "release notes.":http://scop.berkeley.edu/release-notes-1.55.html 
    The latest HIE file can be found
    "elsewhere at SCOP.":http://scop.mrc-lmb.cam.ac.uk/scop/parse/
  
    "Release 1.55":http://scop.berkeley.edu/parse/dir.hie.scop.txt_1.55     
    Records consist of 3 tab-delimited fields: node's sunid,
    parent's sunid, and a list of children's sunids. For example ::
    
    0       -       46456,48724,51349,53931,56572,56835,56992,57942
    21953   49268   -
    49267   49266   49268,49269
    
    Each record holds information for one node in the SCOP hierarchy.

    sunid      -- SCOP unique identifiers of this node
    parent     -- Parents sunid
    children   -- Sequence of childrens sunids
    """
    def __init__(self, record = None):
        self.sunid = None
        self.parent = None
        self.children = []
        
        if not record : return
        
        # Parses HIE records.
        entry = record.rstrip()        # no trailing whitespace
        columns = entry.split('\t')   # separate the tab-delimited cols
        if len(columns) != 3:
            raise ValueError("I don't understand the format of %s" % entry)
        
        self.sunid, self.parent, children = columns

        if self.sunid =='-' : self.sunid = ''
        if self.parent =='-' : self.parent = ''
        else : self.parent = int( self.parent )

        if children =='-' :
            self.children = ()
        else :
            self.children = children.split(',')
            self.children = map ( int, self.children )

        self.sunid = int(self.sunid)
        
    def __str__(self):
        s = []
        s.append(str(self.sunid))

        if self.parent:
            s.append(str(self.parent))
        else:
            if self.sunid != 0:
                s.append('0')
            else:
                s.append('-')
                
        if self.children :
            child_str = map(str, self.children)
            s.append(",".join(child_str))
        else:
            s.append('-')

        return "\t".join(s) + "\n"


    @staticmethod
    def records(hie_file):
        """Iterates over a DOM file, generating DomRecords """
        for line in hie_file:
            if line[0] =='#':  continue  # A comment 
            if line.isspace() : continue
            yield HieRecord(line)

# End HieRecord



class ClaRecord(object):
    """Handle the SCOP CLAssification file, which describes SCOP domains.

    The file format is described in the SCOP
    "release notes.":http://scop.berkeley.edu/release-notes-1.55.html 
    The latest CLA file can be found
    "elsewhere at SCOP.":http://scop.mrc-lmb.cam.ac.uk/scop/parse/

    sid         --  SCOP identifier. e.g. d1danl2
    residues    --  The domain definition as a Residues object
    sccs        --  SCOP concise classification strings.  e.g. b.1.2.1
    sunid       --  SCOP unique identifier for this domain
    hierarchy   --  A sequence of tuples (nodetype, sunid) describing the
                    location of this domain in the SCOP hierarchy.
                    See the Scop module for a description of nodetypes.
    """
    def __init__(self, record=None):
        self.sid = ''
        self.residues = None 
        self.sccs = ''
        self.sunid =''
        self.hierarchy = []
        
        if not record: return
        
        # Parse a tab-delimited CLA record.
        entry = record.rstrip()        # no trailing whitespace
        columns = entry.split('\t')   # separate the tab-delimited cols
        if len(columns) != 6:
            raise ValueError("I don't understand the format of %s" % entry)

        self.sid, pdbid, residues, self.sccs, self.sunid, hierarchy = columns
        self.residues = Residues(residues)
        self.residues.pdbid = pdbid
        self.sunid = int(self.sunid)
        
        h = []
        for ht in hierarchy.split(",") :
            h.append( ht.split('='))        
        for ht in h:
            ht[1] = int(ht[1])
        self.hierarchy = h
    
    def __str__(self):
        s = []
        s.append(self.sid)
        s += str(self.residues).split(" ")
        s.append(self.sccs)
        s.append(self.sunid)

        h=[]
        for ht in self.hierarchy:
             h.append("=".join(map(str,ht))) 
        s.append(",".join(h))
       
        return "\t".join(map(str,s)) + "\n"

    @staticmethod
    def records(cla_file):
        """Iterates over a DOM file, generating DomRecords """
        for line in cla_file:
            if line[0] =='#':  continue  # A comment 
            if line.isspace() : continue
            yield ClaRecord(line)

# End ClaRecord
   
   

class DomRecord(object):
    """Handle the SCOP DOMain file.

    The DOM file has been officially deprecated. For more information see
    the SCOP"release notes.":http://scop.berkeley.edu/release-notes-1.55.html 
    The DOM files for older releases can be found 
    "elsewhere at SCOP.":http://scop.mrc-lmb.cam.ac.uk/scop/parse/

    DOM records consist of 4 tab-delimited fields:
    sid, pdbid, residues, hierarchy
    For example ::
    
    d1sctg_ 1sct    g:      1.001.001.001.001.001
    d1scth_ 1sct    h:      1.001.001.001.001.001
    d1flp__ 1flp    -       1.001.001.001.001.002
    d1moh__ 1moh    -       1.001.001.001.001.002

    sid -- The SCOP ID of the entry, e.g. d1anu1
    residues -- The domain definition as a Residues object
    hierarchy -- A string specifying where this domain is in the hierarchy.
    """
    def __init__(self, record= None):
        self.sid = ''
        self.residues = []
        self.hierarchy = ''
        
        if record:
            entry = record.rstrip()  # no trailing whitespace
            columns = entry.split("\t")  # separate the tab-delimited cols
            if len(columns) != 4:
                raise ValueError("I don't understand the format of %s" % entry)
            self.sid, pdbid, res, self.hierarchy = columns
            self.residues = Residues(res)
            self.residues.pdbid = pdbid
        
    def __str__(self):
        s = []
        s.append(self.sid)
        s.append(str(self.residues).replace(" ","\t") )
        s.append(self.hierarchy)
        return "\t".join(s) + "\n"

    @staticmethod
    def records(dom_file):
        """Iterates over a DOM file, generating DomRecords """
        for line in dom_file:
            if line[0] =='#':  continue  # A comment 
            if line.isspace() : continue
            yield DomRecord(line)

# End DomRecord
    

    

_pdbid_re = re.compile(r"^(\w\w\w\w)(?:$|\s+|_)(.*)")
_fragment_re = re.compile(r"\(?(\w:)?(-?\w*)-?(-?\w*)\)?(.*)")

class Residues(object) :
    """A collection of residues from a PDB structure.

    This class provides code to work with SCOP domain definitions. These
    are concisely expressed as one or more chain fragments. For example,
    "(1bba A:10-20,B:)" indicates residue 10 through 20 (inclusive) of
    chain A, and every residue of chain B in the pdb structure 1bba. The pdb
    id and brackets are optional. In addition "-" indicates every residue of
    a pbd structure with one unnamed chain.

    Start and end residue ids consist of the residue sequence number and an
    optional single letter insertion code. e.g. "12", "-1", "1a", "1000"


    pdbid -- An optional PDB id, e.g. "1bba"
    fragments -- A sequence of tuples (chainID, startResID, endResID)
    """


    def __init__(self, str=None) :
        self.pdbid = ''
        self.fragments = ()
        if str is not None : self._parse(str)


    def _parse(self, string):
        string = string.strip()

        #Is there a pdbid at the front? e.g. 1bba A:1-100
        m = _pdbid_re.match(string)
        if m is not None :
            self.pdbid = m.group(1)
            string = m.group(2) # Everything else

        if string=='' or string == '-' or string=='(-)':  # no fragments, whole sequence
            return
    
        fragments = []
        for l in string.split(",") :
            m = _fragment_re.match(l)
            if m is None:
                raise ValueError("I don't understand the format of %s" % l)
            chain, start, end, postfix = m.groups()

            if postfix != "" :
                 raise ValueError("I don't understand the format of %s" % l )

            if chain:
                if chain[-1] != ':':
                    raise ValueError("I don't understand the chain in %s" % l)
                chain = chain[:-1]   # chop off the ':'
            else :
                chain ="" 
            
            fragments.append((chain, start, end))
        self.fragments = tuple(fragments)
            
    def __str__(self):
        prefix =""
        if self.pdbid :
            prefix =self.pdbid +' '
            
        if not self.fragments: return prefix+'-'
        strs = []
        for chain, start, end in self.fragments:
            s = []
            if chain: s.append("%s:" % chain)
            if start: s.append("%s-%s" % (start, end))
            strs.append("".join(s))
        return prefix+ ",".join(strs)
# End Residues
    



















